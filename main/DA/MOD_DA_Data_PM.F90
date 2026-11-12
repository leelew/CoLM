#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Data_PM
!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Passive-microwave observation handling for data assimilation.
!
!   This module defines one passive-microwave observation stream (PM_type)
!   and a collection of streams (PM_set_type). A PM_type instance represents
!   one source/sensor/variable configuration and owns the observation reader,
!   patch-space forward calculation, and mapping to observation-grid cells.
!
! WORKFLOW:
!   1. Init: Store the stream configuration and microwave frequency.
!      Initialize GridType to read the reference grid and build the mapping.
!   2. Forward: Read current observations through GridType when needed.
!      Apply the selected operator to compute H(x) on model patches.
!      For assimilation, map patch H(x) to the observation grid and extract
!      Values at the current observation locations.
!      Gather and broadcast observation-space H(x) across all processes.
!      For multiple streams, concatenate observations and matching H(x).
!   3. Clear: Release current observations and H(x), retaining the stream
!      Configuration, reference grid, and spatial mapping.
!
! SUPPORTED:
!   TB_H / TB_V:
!      Horizontally / vertically polarized brightness temperatures.
!      Forward modeling uses MOD_DA_Operator_RTM with the configured
!      Frequency and observation-specific incidence angles.
!   SM:
!      Retrieved soil moisture.
!      H(x) is the mean liquid volumetric soil moisture in the upper 5 cm.
!
! NOTES:
!   Forward-only runs compute model-space H(x) for diagnostic output.
!   Assimilation runs additionally prepare observation-space ensemble H(x).
!   Both modes use the configured observation files.
!
! HISTORY:
!   Lu Li, 05/2026: First implementation
!   Lu Li, Zhilong Fan, 11/2026: DA v1.0
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
   USE MOD_DA_Data_GridType
   USE MOD_DA_Operator_RTM
   USE MOD_DA_Vars_TimeVars
   IMPLICIT NONE

!-----------------------------------------------------------------------------
! Passive-microwave observation derived type
!-----------------------------------------------------------------------------
   type :: PM_type

      real(r8)                          :: fghz                          ! Microwave frequency [GHz]
      real(r8), allocatable             :: hx_model(:,:)                 ! Model-space H(x) on local patches
      type(grid_obs_type)               :: y                             ! Observations, reference grid, and patch mapping
      real(r8), allocatable             :: hx_obs(:,:)                   ! Observation-space ensemble H(x)

   CONTAINS

      procedure, PUBLIC                 :: init       => PM_init         ! Initialize the observation reader and grid mapping
      procedure, PUBLIC                 :: forward    => PM_forward      ! Compute patch H(x) and map it to observation locations
      procedure, PUBLIC                 :: clear      => PM_clear        ! Clear current observations and H(x), retaining the mapping

   END type PM_type

   type :: PM_config_type

      character(len=16)                 :: source_name
      character(len=32)                 :: sensor_name
      character(len=16)                 :: var_name
      real(r8)                          :: fghz

   END type PM_config_type

!-----------------------------------------------------------------------------
! Collection of passive-microwave observation streams and their combined data
!-----------------------------------------------------------------------------
   type :: PM_set_type

      integer                           :: nsensor = 0                   ! Number of configured observation streams
      type(PM_type), allocatable        :: sensors(:)                    ! Individual observation streams
      type(PM_config_type), allocatable :: cfgs(:)                       ! Configuration of each observation stream

      integer                           :: nobs = 0                      ! Total number of observations across streams
      real(r8), allocatable             :: lat(:)                        ! Observation latitudes
      real(r8), allocatable             :: lon(:)                        ! Observation longitudes
      real(r8), allocatable             :: y(:)                          ! Observed values
      real(r8), allocatable             :: r(:)                          ! Observation error standard deviations
      real(r8), allocatable             :: hx_obs(:,:)                   ! Observation-space ensemble H(x)

   CONTAINS

      procedure, PUBLIC                 :: init       => PM_set_init     ! Initialize all observation streams
      procedure, PUBLIC                 :: concat     => PM_set_concat   ! Combine observations and matching H(x) in stream order
      procedure, PUBLIC                 :: forward    => PM_set_forward  ! Run forward calculations for all streams and combine their records
      procedure, PUBLIC                 :: clear      => PM_set_clear    ! Clear current observations and H(x), retaining stream configurations

   END type PM_set_type

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE PM_init(this, pixelset, source_name, sensor_name, var_name, fghz)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(PM_type),      intent(inout) :: this
      type(pixelset_type), intent(in)    :: pixelset
      character(len=*),    intent(in)    :: source_name
      character(len=*),    intent(in)    :: sensor_name
      character(len=*),    intent(in)    :: var_name
      real(r8),            intent(in)    :: fghz

      CALL this%y%init(pixelset, source_name, sensor_name, var_name)
      this%fghz = fghz

   END SUBROUTINE PM_init

!-----------------------------------------------------------------------------

   SUBROUTINE PM_forward(this, idate, deltim)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------------ Dummy Arguments ------------------------------
      class(PM_type),     intent(inout)   :: this
      integer,            intent(in)      :: idate(3)
      real(r8),           intent(in)      :: deltim

!------------------------------ Local Variables ------------------------------
      type(block_data_real8_2d)           :: theta_wgrid                ! Incidence angles on the observation grid
      type(pointer_real8_1d), allocatable :: theta_part(:)              ! Incidence angles on patch-grid parts
      logical,  allocatable               :: filter(:)                  ! Local patch mask used for mapping
      real(r8)                            :: tb_part(2)                 ! Brightness temperatures of the current patch-grid intersection
      real(r8)                            :: sumarea(0:DEF_DA_ENS_NUM)  ! Valid part area for each member
      real(r8)                            :: area_obs                   ! Mapped patch area in one observation grid cell
      type(block_data_real8_2d)           :: area_wgrid                 ! Mapped patch area on obs world grid
      type(block_data_real8_3d)           :: hx_wgrid                   ! Ensemble H(x) mapped to obs world grid
      integer,  allocatable               :: iloc(:)
      integer,  allocatable               :: tmp_idx(:)
      real(r8), allocatable               :: tmp_data(:,:)
      integer                             :: iobs, iens, np, ipart, ipol
      integer                             :: ib, jb, il, jl
      integer                             :: ip, ndata, isrc, smesg(2), rmesg(2)

!-----------------------------------------------------------------------------
      CALL this%clear()

      ! Read observations
      CALL this%y%read(idate, deltim)
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
      allocate(this%hx_model(0:DEF_DA_ENS_NUM, numpatch))
      this%hx_model = spval

      IF (p_is_worker) THEN
         DO np = 1, numpatch
            IF (patchtype(np) /= 1) CYCLE

            SELECT CASE (trim(this%y%var_name))
!-----------------------------------------------------------------------------
! Level-1 brightness temperature observations from passive microwave radiometers.
!-----------------------------------------------------------------------------
            CASE ('TB_H', 'TB_V')
               ipol = merge(2, 1, this%y%var_name == 'TB_V')
               sumarea(:) = 0.0_r8

               DO ipart = 1, this%y%mg2p%npart(np)
                  IF (theta_part(np)%val(ipart) == spval) CYCLE

                  DO iens = 0, DEF_DA_ENS_NUM
                     CALL PM_RTM(                                                                     &
                        patchtype(np), patchclass(np), dz_sno_ens(:,iens,np),                         &
                        forc_topo(np), htop(np),                                                      &
                        tref_ens(iens,np), t_soisno_ens(:,iens,np), tleaf_ens(iens,np),               &
                        wliq_soisno_ens(:,iens,np), wice_soisno_ens(:,iens,np),                       &
                        snowdp_ens(iens,np),                                                          &
                        lai_ens(iens,np), sai_ens(iens,np),                                           &
                        wf_clay(:,np), wf_sand(:,np), wf_silt(:,np),                                  &
                        BD_all(:,np), porsl(:,np),                                                    &
                        theta_part(np)%val(ipart), this%fghz,                                         &
                        tb_part(1), tb_part(2))

                     IF (tb_part(ipol) == spval) CYCLE
                     IF (sumarea(iens) == 0) this%hx_model(iens,np) = 0
                     sumarea(iens) = sumarea(iens) + this%y%mg2p%areapart(np)%val(ipart)
                     this%hx_model(iens,np) = this%hx_model(iens,np) +                                &
                        tb_part(ipol) * this%y%mg2p%areapart(np)%val(ipart)
                  ENDDO
               ENDDO

               WHERE (sumarea > 0.0)
                  this%hx_model(:,np) = this%hx_model(:,np) / sumarea
               ENDWHERE

!-----------------------------------------------------------------------------
! Retrieved soil moisture products from passive microwave observations.
!-----------------------------------------------------------------------------
            CASE ('SM')
               DO iens = 0, DEF_DA_ENS_NUM
                  IF (wliq_soisno_ens(1,iens,np) == spval) CYCLE
                  this%hx_model(iens,np) = (wliq_soisno_ens(1,iens,np) + wliq_soisno_ens(2,iens,np) + &
                     wliq_soisno_ens(3,iens,np)*(0.05-0.0451)/(0.0906-0.0451)) / (0.05*denh2o)
               ENDDO
            END SELECT
         ENDDO
      ENDIF

      IF (.not. DEF_DA) RETURN

      ! Keep the mapped numerator and area denominator on the same valid patches.
      allocate (filter(numpatch))
      IF (p_is_worker) filter(:) = all(this%hx_model(1:DEF_DA_ENS_NUM,:) /= spval, dim=1)

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
                  IF (area_obs > 0 .and. all(hx_wgrid%blk(ib,jb)%val(:,il,jl) /= spval)) THEN
                     ndata = ndata + 1
                     iloc(ndata) = iobs
                     this%hx_obs(ndata,:) = hx_wgrid%blk(ib,jb)%val(:,il,jl) / area_obs
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
         deallocate (iloc)
         IF (allocated(hx_wgrid%blk  )) deallocate(hx_wgrid%blk  )
         IF (allocated(area_wgrid%blk)) deallocate(area_wgrid%blk)
      ENDIF

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

   END SUBROUTINE PM_forward

!-----------------------------------------------------------------------------

   SUBROUTINE PM_clear(this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(PM_type), intent(inout) :: this

      CALL this%y%clear()
      IF (allocated(this%hx_model)) deallocate(this%hx_model)
      IF (allocated(this%hx_obs  )) deallocate(this%hx_obs  )

   END SUBROUTINE PM_clear

!-----------------------------------------------------------------------------

   SUBROUTINE PM_set_init(this, pixelset, configs)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(PM_set_type),   intent(inout) :: this
      type(pixelset_type),  intent(in)    :: pixelset
      type(PM_config_type), intent(in)    :: configs(:)
      integer                             :: i

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
         CALL this%sensors(i)%init(pixelset, this%cfgs(i)%source_name, &
            this%cfgs(i)%sensor_name, this%cfgs(i)%var_name, this%cfgs(i)%fghz)
      ENDDO

   END SUBROUTINE PM_set_init

!-----------------------------------------------------------------------------

   SUBROUTINE PM_set_clear(this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(PM_set_type), intent(inout) :: this
      integer                           :: is

      this%nobs = 0
      IF (allocated(this%lat   )) deallocate(this%lat   )
      IF (allocated(this%lon   )) deallocate(this%lon   )
      IF (allocated(this%y     )) deallocate(this%y     )
      IF (allocated(this%r     )) deallocate(this%r     )
      IF (allocated(this%hx_obs)) deallocate(this%hx_obs)

      IF (.not. allocated(this%sensors)) return

      DO is = 1, this%nsensor
         CALL this%sensors(is)%clear()
      ENDDO

   END SUBROUTINE PM_set_clear

!-----------------------------------------------------------------------------

   SUBROUTINE PM_set_concat(this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(PM_set_type), intent(inout) :: this
      integer                           :: is, n, offset

      this%nobs = 0
      IF (allocated(this%lat   )) deallocate(this%lat   )
      IF (allocated(this%lon   )) deallocate(this%lon   )
      IF (allocated(this%y     )) deallocate(this%y     )
      IF (allocated(this%r     )) deallocate(this%r     )
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

   SUBROUTINE PM_set_forward(this, idate, deltim)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(PM_set_type),  intent(inout) :: this
      integer,            intent(in)     :: idate(3)
      real(r8),           intent(in)     :: deltim
      integer                            :: is

      IF (.not. allocated(this%sensors)) RETURN

      DO is = 1, this%nsensor
         CALL this%sensors(is)%forward(idate, deltim)
      ENDDO

      IF (DEF_DA) CALL this%concat()

   END SUBROUTINE PM_set_forward

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Data_PM
#endif
