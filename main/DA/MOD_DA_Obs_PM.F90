#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Obs_PM
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Passive-microwave observation handling for data assimilation.
!
!    This module defines one passive-microwave observation stream (PM_type)
!    and a collection of streams (PM_set_type). A PM_type instance represents
!    one source/sensor/variable configuration and owns the observation reader,
!    patch-space forward calculation, and mapping to observation-grid cells.
!
!    The observation-space workflow is:
!      1. clear time-step data and read current observations through
!         grid_obs_type;
!      2. map observation incidence angle from grid cells to patch-grid
!         intersection parts;
!      3. compute ensemble H(x) on valid parts with the passive-microwave
!         radiative-transfer operator, then area-weight valid parts back to
!         patch-space H(x);
!      4. map patch-space H(x) to the observation grid;
!      5. crop mapped H(x) to valid observation locations;
!      6. gather IO-task results on the master and broadcast this%hx_obs so all
!         MPI tasks share the same observation-space H(x).
!
!    PM_set_type applies this workflow to all configured PM streams and
!    concatenates their observation values, errors, metadata, and H(x) into
!    set-level arrays for assimilation.
!
!    PM observation operator supports both satellite brightness temperature
!    assimilation through microwave radiative transfer modeling and the direct
!    assimilation of retrieved soil moisture products.
!
! HISTORY:
!    Lu Li, 05/2026: First implementation
!    Lu Li, Zhilong Fan, 11/2026: DA v1.0
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Const_Physical
   USE MOD_DataType
   USE MOD_Vars_TimeInvariants
   USE MOD_Vars_1DForcing
   USE MOD_Vars_Global
   USE MOD_LandPatch
   USE MOD_Block
   USE MOD_Namelist
   USE MOD_NetCDFSerial
   USE MOD_SPMD_Task
   USE MOD_Pixelset
   USE MOD_TimeManager
   USE MOD_DA_Obs_GridType
   USE MOD_DA_Operator_RTM
   USE MOD_DA_Vars_TimeVars
   IMPLICIT NONE

!-----------------------------------------------------------------------------
! One passive-microwave observation stream
!-----------------------------------------------------------------------------
   type :: PM_type

      ! Metadata for this observation stream.
      real(r8)              :: fghz

      ! Model-space H(x) on local patches, retained for diagnostic output.
      real(r8), allocatable :: hx_model(:,:)               ! [0:ens, patch]

      ! Observations and metadata at the locations read for this time step.
      type(grid_obs_type)   :: y                           ! Observation reader and grid mapping

      ! Model equivalents at observation locations; unmapped values are spval.
      real(r8), allocatable :: hx_obs(:,:)                 ! [nobs, ens]

   CONTAINS

      procedure, PUBLIC :: init       => PM_init           ! Initialize the observation reader and grid mapping
      procedure, PUBLIC :: calc_hx    => PM_calc_hx        ! Compute H(x) on patches and map it to observation locations
      procedure, PUBLIC :: clear      => PM_clear          ! Clear time-step data while retaining the mapping

   END type PM_type


   type :: PM_config_type

      character(len=16)  :: source_name
      character(len=32)  :: sensor_name
      character(len=16)  :: var_name
      real(r8)           :: fghz

   END type PM_config_type

!-----------------------------------------------------------------------------
! Collection of passive-microwave observation streams and their combined data
!-----------------------------------------------------------------------------
   type :: PM_set_type

      ! Each stream represents one source/sensor/variable configuration.
      integer :: nsensor = 0                               ! Number of configured observation streams
      type(PM_type), allocatable :: sensors(:)             ! Stream objects [nsensor]
      type(PM_config_type), allocatable :: cfgs(:)         ! Corresponding stream configurations [nsensor]

      ! Current time-step data concatenated in stream order for assimilation.
      integer :: nobs = 0                                  ! Total number of observations across streams in the combined arrays
      real(r8), allocatable :: lat(:)                      ! Observation latitude [nobs], degrees
      real(r8), allocatable :: lon(:)                      ! Observation longitude [nobs], degrees
      real(r8), allocatable :: y(:)                        ! Actual observation values [nobs]
      real(r8), allocatable :: r(:)                        ! Observation error standard deviations [nobs]
      real(r8), allocatable :: hx_obs(:,:)                 ! Simulated observations at y locations [nobs, ens]

   CONTAINS

      procedure, PUBLIC :: init       => PM_set_init       ! Initialize all streams from their configurations
      procedure, PUBLIC :: concat     => PM_set_concat     ! Combine stream data into the set-level arrays
      procedure, PUBLIC :: calc_hx    => PM_set_calc_hx    ! Compute both H(x) representations and combine streams
      procedure, PUBLIC :: clear      => PM_set_clear      ! Clear time-step data while retaining stream configurations

   END type PM_set_type


!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE PM_init (this, pixelset, source_name, sensor_name, var_name, fghz)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(PM_type),      intent(inout) :: this
      type(pixelset_type), intent(in)    :: pixelset       ! target pixel/patch set
      character(len=*),    intent(in)    :: source_name    ! observation source class, e.g., "PM"
      character(len=*),    intent(in)    :: sensor_name    ! sensor/subdirectory name
      character(len=*),    intent(in)    :: var_name       ! observation variable/file prefix
      real(r8),            intent(in)    :: fghz           ! frequency [GHz]

!-----------------------------------------------------------------------------
      CALL this%y%init(pixelset, source_name, sensor_name, var_name)
      this%fghz = fghz

   END SUBROUTINE PM_init

!-----------------------------------------------------------------------------

   SUBROUTINE PM_calc_hx (this, idate, deltim)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(PM_type),      intent(inout) :: this
      integer,            intent(in)    :: idate(3)
      real(r8),           intent(in)    :: deltim

!------------------------ Local Variables ------------------------------------
      integer :: iobs                                      ! observation index in this%y
      integer :: ip                                        ! IO task loop index on master
      integer :: ib, jb                                    ! block indices of an observation grid cell
      integer :: il, jl                                    ! local indices inside the block

      integer :: ndata                                     ! valid obs count handled by one IO task
      integer :: isrc                                      ! MPI source rank for received obs data
      integer :: smesg(2)                                  ! IO-to-master message header: rank, ndata
      integer :: rmesg(2)                                  ! master receive message header: rank, ndata

      real(r8) :: area_obs                                 ! mapped patch area in one observation grid cell
      type(block_data_real8_2d) :: area_wgrid              ! mapped patch area on obs world grid (this%y%grid)
      type(block_data_real8_3d) :: hx_wgrid                ! ensemble H(x) mapped to obs world grid (this%y%grid)

      integer,  allocatable :: iloc(:)                     ! obs indices handled by one IO task
      integer,  allocatable :: tmp_idx(:)                  ! obs indices received by master
      real(r8), allocatable :: tmp_data(:,:)               ! H(x) values sent/received by MPI
      logical,  allocatable :: filter(:)                   ! local patch mask used for mapping
      integer :: iens                                      ! Ensemble member index
      integer :: np                                        ! Patch index
      integer :: ipart                                     ! Patch-grid intersection part index
      integer :: ipol                                      ! Polarization index
      real(r8) :: hx_tmp(2)                                ! Brightness temperatures for both polarizations
      real(r8) :: sumarea(0:DEF_DA_ENS_NUM)                ! Valid part area for each member
      type(block_data_real8_2d) :: theta_wgrid             ! Incidence angles on the observation grid
      type(pointer_real8_1d), allocatable :: theta_part(:) ! Incidence angles on patch-grid parts

!-----------------------------------------------------------------------------
      ! Reuse observations within the time step; clear releases them after analysis.
      IF (.not. allocated(this%y%y)) CALL this%y%read(idate, deltim)
      IF (allocated(this%hx_obs)) deallocate(this%hx_obs)

      IF (allocated(this%hx_model)) deallocate(this%hx_model)
      allocate(this%hx_model(0:DEF_DA_ENS_NUM, numpatch))
      this%hx_model = spval

      IF (this%y%nobs == 0) RETURN

      ! Map incidence angle from obs grid cells to patch-grid intersection parts.
      IF (p_is_io) THEN
         CALL allocate_block_data(this%y%grid, theta_wgrid)
         CALL flush_block_data(theta_wgrid, spval)

         DO iobs = 1, this%y%nobs
            ib = this%y%grid%xblk(this%y%jj(iobs))
            jb = this%y%grid%yblk(this%y%ii(iobs))
            il = this%y%grid%xloc(this%y%jj(iobs))
            jl = this%y%grid%yloc(this%y%ii(iobs))
            IF (ib /= 0 .and. jb /= 0) THEN
               IF (gblock%pio(ib,jb) == p_iam_glb) THEN
                  theta_wgrid%blk(ib,jb)%val(il,jl) = this%y%theta(iobs)
               ENDIF
            ENDIF
         ENDDO
      ENDIF
      CALL this%y%mg2p%allocate_part(theta_part)
      CALL this%y%mg2p%grid2part(theta_wgrid, theta_part)

      ! Compute part-level H(x), then area-average valid parts to patch space.
      IF (p_is_worker) THEN
         allocate(filter(numpatch))
         filter(:) = patchtype(:) <= 2

         DO np = 1, numpatch
            IF (this%y%mg2p%npart(np) == 0) CYCLE
            IF (.not. filter(np)) CYCLE

            ! the H(x) operator.
            SELECT CASE (trim(this%y%var_name))
!-----------------------------------------------------------------------------
! Level-1 brightness temperature observations from passive microwave radiometers.
! Forward modeling uses the passive microwave radiative transfer operator
!-----------------------------------------------------------------------------
            CASE ('TB_H', 'TB_V')
               ipol = merge(2, 1, this%y%var_name == 'TB_V')
               sumarea(:) = 0.0_r8

               ! TB is evaluated only where the obs file supplies theta,
               ! so a domain-wide theta field gives domain-wide TB.
               DO ipart = 1, this%y%mg2p%npart(np)
                  IF (theta_part(np)%val(ipart) == spval) CYCLE

                  DO iens = 0, merge(DEF_DA_ENS_NUM, 0, DEF_DA)
                     CALL PM_RTM( &
                        patchtype(np), patchclass(np), dz_sno_ens(:,iens,np), &
                        forc_topo(np), htop(np), &
                        tref_ens(iens,np), t_soisno_ens(:,iens,np), tleaf_ens(iens,np), &
                        wliq_soisno_ens(:,iens,np), wice_soisno_ens(:,iens,np), &
                        snowdp_ens(iens,np), &
                        lai_ens(iens,np), sai_ens(iens,np), &
                        wf_clay(:,np), wf_sand(:,np), wf_silt(:,np), &
                        BD_all(:,np), porsl(:,np), &
                        theta_part(np)%val(ipart), this%fghz, &
                        hx_tmp(1), hx_tmp(2))

                     IF (hx_tmp(ipol) == spval) CYCLE
                     IF (sumarea(iens) == 0) this%hx_model(iens,np) = 0
                     sumarea(iens) = sumarea(iens) + this%y%mg2p%areapart(np)%val(ipart)
                     this%hx_model(iens,np) = this%hx_model(iens,np) + &
                        hx_tmp(ipol) * this%y%mg2p%areapart(np)%val(ipart)
                  ENDDO
               ENDDO

               WHERE (sumarea > 0.0)
                  this%hx_model(:,np) = this%hx_model(:,np) / sumarea
               ENDWHERE

!-----------------------------------------------------------------------------
! Retrieved soil moisture products from passive microwave observations.
! H(x) is the mean liquid volumetric soil moisture over the upper 5 cm [m3/m3].
!-----------------------------------------------------------------------------
            CASE ('SM')
               DO iens = 0, merge(DEF_DA_ENS_NUM, 0, DEF_DA)
                  IF (wliq_soisno_ens(1,iens,np) == spval) CYCLE
                  this%hx_model(iens,np) = &
                     (wliq_soisno_ens(1,iens,np) + &
                      wliq_soisno_ens(2,iens,np) + &
                      wliq_soisno_ens(3,iens,np)*(0.05-0.0451)/(0.0906-0.0451)) / &
                     (0.05*denh2o)
               ENDDO
            CASE DEFAULT
               CYCLE
            END SELECT
         ENDDO

         deallocate(filter)
      ENDIF

      CALL this%y%mg2p%deallocate_part(theta_part)
      IF (p_is_io) THEN
         IF (allocated(theta_wgrid%blk)) deallocate(theta_wgrid%blk)
      ENDIF

      IF (.not. DEF_DA) RETURN

      ! Keep the mapped numerator and area denominator on the same valid patches.
      allocate (filter(numpatch))
      IF (p_is_worker) filter(:) = patchtype(:) <= 2 .and. &
         all(this%hx_model(1:DEF_DA_ENS_NUM,:) /= spval, dim=1)

      ! Map model-space H(x) to the observation grid.
      IF (p_is_io) THEN
         CALL allocate_block_data(this%y%grid, area_wgrid)
         CALL allocate_block_data(this%y%grid, hx_wgrid, DEF_DA_ENS_NUM)
      ENDIF
      CALL this%y%mg2p%get_sumarea(area_wgrid, filter)
      CALL this%y%mg2p%pset2grid(this%hx_model(1:DEF_DA_ENS_NUM,:), hx_wgrid, spv=spval, msk=filter)
      deallocate (filter)

      ! Extract H(x) at the current observation locations.
      allocate (this%hx_obs(this%y%nobs, DEF_DA_ENS_NUM))

      IF (p_is_io) THEN
         allocate (iloc(this%y%nobs))
         this%hx_obs = spval

         ndata = 0
         DO iobs = 1, this%y%nobs
            ib = this%y%grid%xblk(this%y%jj(iobs))
            jb = this%y%grid%yblk(this%y%ii(iobs))
            il = this%y%grid%xloc(this%y%jj(iobs))
            jl = this%y%grid%yloc(this%y%ii(iobs))
            IF (ib /= 0 .and. jb /= 0) THEN
               IF (gblock%pio(ib,jb) == p_iam_glb) THEN
                  area_obs = area_wgrid%blk(ib,jb)%val(il,jl)
                  IF (area_obs > 0 .and. &
                     all(hx_wgrid%blk(ib,jb)%val(:,il,jl) /= spval)) THEN
#ifdef USEMPI
                     ndata = ndata + 1
                     iloc(ndata) = iobs
                     this%hx_obs(ndata,:) = hx_wgrid%blk(ib,jb)%val(:,il,jl) / area_obs
#endif
                  ENDIF
               ENDIF
            ENDIF
         ENDDO

#ifdef USEMPI
         smesg = (/p_iam_glb, ndata/)
         CALL mpi_send(smesg, 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)

         IF (ndata > 0) THEN
            CALL mpi_send(iloc(1:ndata), ndata, MPI_INTEGER, p_address_master, mpi_tag_data, p_comm_glb, p_err)
            allocate (tmp_data(ndata, DEF_DA_ENS_NUM))
            tmp_data = this%hx_obs(1:ndata,:)
            CALL mpi_send(tmp_data, ndata*DEF_DA_ENS_NUM, MPI_REAL8, p_address_master, mpi_tag_data+1, p_comm_glb, p_err)
            deallocate (tmp_data)
         ENDIF
#endif
         deallocate (iloc)
         IF (allocated(hx_wgrid%blk))   deallocate (hx_wgrid%blk)
         IF (allocated(area_wgrid%blk)) deallocate (area_wgrid%blk)
      ENDIF

#ifdef USEMPI
      IF (p_is_master) THEN
         this%hx_obs = spval

         DO ip = 0, p_np_io - 1
            CALL mpi_recv(rmesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            ndata = rmesg(2)
            IF (ndata > 0) THEN
               allocate (tmp_idx(ndata))
               allocate (tmp_data(ndata, DEF_DA_ENS_NUM))

               isrc = rmesg(1)
               CALL mpi_recv(tmp_idx, ndata, MPI_INTEGER, isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL mpi_recv(tmp_data, ndata*DEF_DA_ENS_NUM, MPI_REAL8, isrc, mpi_tag_data+1, p_comm_glb, p_stat, p_err)
               this%hx_obs(tmp_idx,:) = tmp_data

               deallocate (tmp_idx)
               deallocate (tmp_data)
            ENDIF
         ENDDO
      ENDIF
      CALL mpi_bcast(this%hx_obs, this%y%nobs*DEF_DA_ENS_NUM, MPI_REAL8, p_address_master, p_comm_glb, p_err)
#endif

   END SUBROUTINE PM_calc_hx

!-----------------------------------------------------------------------------

   SUBROUTINE PM_clear (this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(PM_type), intent(inout) :: this

!-----------------------------------------------------------------------------
      CALL this%y%clear()

      IF (allocated(this%hx_model)) deallocate(this%hx_model)
      IF (allocated(this%hx_obs))      deallocate(this%hx_obs)

   END SUBROUTINE PM_clear

!-----------------------------------------------------------------------------

   SUBROUTINE PM_set_init (this, pixelset, configs)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(PM_set_type),   intent(inout) :: this
      type(pixelset_type),  intent(in)    :: pixelset
      type(PM_config_type), intent(in)    :: configs(:)

      integer :: i

      CALL this%clear()
      IF (allocated(this%sensors)) THEN
         DO i = 1, size(this%sensors)
            CALL this%sensors(i)%clear()
         ENDDO
         deallocate(this%sensors)
      ENDIF
      IF (allocated(this%cfgs)) deallocate(this%cfgs)

      this%nsensor = 0
      IF (size(configs) == 0) return

      this%nsensor = size(configs)

      allocate(this%cfgs(size(configs)))
      allocate(this%sensors(size(configs)))
      DO i = 1, size(configs)
         this%cfgs(i) = configs(i)
         CALL this%sensors(i)%init( &
            pixelset, &
            this%cfgs(i)%source_name, &
            this%cfgs(i)%sensor_name, &
            this%cfgs(i)%var_name, &
            this%cfgs(i)%fghz)
      ENDDO

   END SUBROUTINE PM_set_init

!-----------------------------------------------------------------------------

   SUBROUTINE PM_set_clear (this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(PM_set_type), intent(inout) :: this
      integer :: is

!-----------------------------------------------------------------------------
      this%nobs = 0
      IF (allocated(this%lat))       deallocate(this%lat)
      IF (allocated(this%lon))       deallocate(this%lon)
      IF (allocated(this%y))         deallocate(this%y)
      IF (allocated(this%r))         deallocate(this%r)
      IF (allocated(this%hx_obs))        deallocate(this%hx_obs)

      IF (.not. allocated(this%sensors)) return

      DO is = 1, this%nsensor
         CALL this%sensors(is)%clear()
      ENDDO

   END SUBROUTINE PM_set_clear

!-----------------------------------------------------------------------------

   SUBROUTINE PM_set_concat (this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(PM_set_type), intent(inout) :: this
      integer :: is, n, offset

!-----------------------------------------------------------------------------
      this%nobs = 0
      IF (allocated(this%lat)) deallocate(this%lat)
      IF (allocated(this%lon)) deallocate(this%lon)
      IF (allocated(this%y)) deallocate(this%y)
      IF (allocated(this%r)) deallocate(this%r)
      IF (allocated(this%hx_obs)) deallocate(this%hx_obs)
      IF (.not. allocated(this%sensors)) RETURN

      ! Count observations with model equivalents and allocate combined arrays.
      DO is = 1, this%nsensor
         IF (this%sensors(is)%y%nobs == 0) CYCLE
         IF (.not. allocated(this%sensors(is)%hx_obs)) CYCLE
         this%nobs = this%nobs + this%sensors(is)%y%nobs
      ENDDO
      IF (this%nobs == 0) return

      allocate(this%lat(this%nobs))
      allocate(this%lon(this%nobs))
      allocate(this%y(this%nobs))
      allocate(this%r(this%nobs))
      allocate(this%hx_obs(this%nobs, DEF_DA_ENS_NUM))

      ! Combine observations and H(x) in stream order.
      offset = 0
      DO is = 1, this%nsensor
         IF (.not. allocated(this%sensors(is)%hx_obs)) CYCLE

         n = this%sensors(is)%y%nobs
         IF (n == 0) CYCLE

         this%lat(offset+1:offset+n) = this%sensors(is)%y%lat
         this%lon(offset+1:offset+n) = this%sensors(is)%y%lon
         this%y(offset+1:offset+n) = this%sensors(is)%y%y
         this%r(offset+1:offset+n) = this%sensors(is)%y%r
         this%hx_obs(offset+1:offset+n,:) = this%sensors(is)%hx_obs

         offset = offset + n
      ENDDO

   END SUBROUTINE PM_set_concat

!-----------------------------------------------------------------------------

   SUBROUTINE PM_set_calc_hx(this, idate, deltim)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(PM_set_type),  intent(inout) :: this
      integer,            intent(in)    :: idate(3)
      real(r8),           intent(in)    :: deltim

      integer :: is

!-----------------------------------------------------------------------------
      IF (.not. allocated(this%sensors)) RETURN

      DO is = 1, this%nsensor
         CALL this%sensors(is)%calc_hx(idate, deltim)
      ENDDO

      IF (DEF_DA) CALL this%concat()

   END SUBROUTINE PM_set_calc_hx

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Obs_PM
#endif
