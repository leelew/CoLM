#include <define.h>
#ifdef DataAssimilation
MODULE MOD_DA_Data_MR
!-----------------------------------------------------------------------------
! DESCRIPTION:
!   MODIS snow-cover observation source for snow-depth data assimilation.
!
! WORKFLOW:
!   1. Init: Store stream configurations and initialize the reference grid
!      and spatial mapping through GridType.
!   2. Forward: Read observations, compute model-space H(x) using the snow
!      depletion curve, map ensemble H(x) to observation locations, and
!      gather and broadcast the results. Concatenate configured streams.
!   3. Clear: Release current observations and H(x), retaining configurations
!      and spatial mappings.
!
! HISTORY:
!   Zhuo Wang, 11/2026: Snow observation source refactoring
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_Pixelset
   USE MOD_Vars_Global
   USE MOD_Vars_TimeInvariants
   USE MOD_Vars_TimeVariables, ONLY: z0m
   USE MOD_DA_Data_GridType
   USE MOD_SnowFraction
   USE MOD_Block
   USE MOD_DA_Vars_TimeVars
   IMPLICIT NONE

!-----------------------------------------------------------------------------
! MODIS observation derived type.
!-----------------------------------------------------------------------------
   type :: MR_type

      real(r8), allocatable               :: hx_model(:,:)                 ! Model-space H(x) on local patches
      type(grid_obs_type)                 :: y                             ! Observations, reference grid, and patch mapping
      real(r8), allocatable               :: hx_obs(:,:)                   ! Observation-space ensemble H(x)

   CONTAINS

      procedure, PUBLIC                   :: init       => MR_init         ! Initialize the observation reader and grid mapping
      procedure, PUBLIC                   :: forward    => MR_forward      ! Compute model- and observation-space H(x)
      procedure, PUBLIC                   :: clear      => MR_clear        ! Clear current observations and H(x), retaining the mapping

   END type MR_type

   type :: MR_config_type

      character(len=16)                   :: source_name
      character(len=32)                   :: sensor_name
      character(len=16)                   :: var_name

   END type MR_config_type

!-----------------------------------------------------------------------------
! Collection of MODIS observation streams and their combined data.
!-----------------------------------------------------------------------------
   type :: MR_set_type

      integer                             :: nsensor = 0                   ! Number of configured observation streams
      type(MR_type), allocatable          :: sensors(:)                    ! Individual observation streams
      type(MR_config_type), allocatable   :: cfgs(:)                       ! Configuration of each observation stream

      integer                             :: nobs = 0                      ! Total number of observations across streams
      real(r8), allocatable               :: lat(:)                        ! Observation latitudes
      real(r8), allocatable               :: lon(:)                        ! Observation longitudes
      real(r8), allocatable               :: y(:)                          ! Observed values
      real(r8), allocatable               :: sigma(:)                      ! Observation error standard deviations
      real(r8), allocatable               :: hx_obs(:,:)                   ! Observation-space ensemble H(x)

   CONTAINS

      procedure, PUBLIC                   :: init       => MR_set_init
      procedure, PUBLIC                   :: concat     => MR_set_concat
      procedure, PUBLIC                   :: forward    => MR_set_forward
      procedure, PUBLIC                   :: clear      => MR_set_clear

   END type MR_set_type

!-----------------------------------------------------------------------------
CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE MR_init(this, pixelset, source_name, sensor_name, var_name)

!-----------------------------------------------------------------------------
      class(MR_type), intent(inout)       :: this
      type(pixelset_type), intent(in)     :: pixelset ! target pixel/patch set
      character(len=*), intent(in)        :: source_name ! observation source class
      character(len=*), intent(in)        :: sensor_name ! sensor/subdirectory name
      character(len=*), intent(in)        :: var_name    ! observation variable/file prefix
!-----------------------------------------------------------------------------
      CALL this%y%init(pixelset, source_name, sensor_name, var_name)

   END SUBROUTINE MR_init

!-----------------------------------------------------------------------------

   SUBROUTINE MR_forward(this, idate, deltim)

!-----------------------------------------------------------------------------
      class(MR_type), intent(inout)       :: this   ! MODIS observation source object
      integer,        intent(in)          :: idate(3) ! Current model date
      real(r8),       intent(in)          :: deltim   ! Model time step (s)

      integer                             :: iens, np, nens
      real(r8)                            :: fwet_snow_local, sigf_local
      logical, allocatable                :: filter(:)       ! Valid patch mask for grid mapping
      integer                             :: i, ip, ndata, isrc           ! Loop indices, valid count, and MPI source
      integer                             :: ilat, ilon, ib, jb, il, jl   ! Observation grid and block indices
      integer                             :: smesg(2), rmesg(2)           ! MPI message headers (rank and count)
      integer, allocatable                :: iloc(:)         ! Local observation indices on IO process
      integer, allocatable                :: itemp(:)        ! Received observation indices on master
      real(r8), allocatable               :: dtemp(:,:)     ! MPI buffer for ensemble H(x)
      real(r8)                            :: dummy_pdata(1,1)            ! Placeholder for non-worker processes
      real(r8)                            :: area_obs                   ! Total mapped patch area in grid cell
      type(block_data_real8_2d)           :: area_wgrid ! Area-weight sums on observation grid
      type(block_data_real8_3d)           :: hx_wgrid   ! Gridded ensemble H(x) before normalization

!-----------------------------------------------------------------------------
      CALL this%clear()

      ! Read observations.
      CALL this%y%read(idate, deltim)
      IF (this%y%nobs == 0) RETURN

      IF (p_is_worker) THEN
         IF (allocated(this%hx_model)) deallocate(this%hx_model)
         allocate(this%hx_model(0:DEF_DA_ENS_NUM,numpatch))
         this%hx_model = spval
         ! Compute model-space H(x) without changing the ensemble states.
         nens = 0
         IF (DEF_DA) nens = DEF_DA_ENS_NUM
         DO iens = 0, nens
            DO np = 1, numpatch
               IF (patchtype(np) >= 2) CYCLE
               IF (snowdp_ens(iens,np) == spval .or. scv_ens(iens,np) == spval) CYCLE
               fwet_snow_local = fwet_snow_ens(iens,np)
               sigf_local = sigf_ens(iens,np)
               CALL snowfraction( &
                  lai_ens(iens,np), &
                  sai_ens(iens,np), &
                  z0m(np), &
                  zlnd, &
                  scv_ens(iens,np), &
                  snowdp_ens(iens,np), &
                  fwet_snow_local, &
                  sigf_local, &
                  this%hx_model(iens,np))
            ENDDO
         ENDDO
      ENDIF

      IF (.not. DEF_DA) RETURN
#ifndef SinglePoint
      IF (p_is_io) THEN
         CALL allocate_block_data(this%y%grid, area_wgrid)
      ENDIF
      ! Calculate sum area
      IF (p_is_worker) THEN
         allocate(filter(numpatch))
         filter = patchtype < 2 .and. &
                  all(this%hx_model(1:DEF_DA_ENS_NUM,:) /= spval, DIM=1)
         CALL this%y%mg2p%get_sumarea(area_wgrid, filter)
      ELSE
         CALL this%y%mg2p%get_sumarea(area_wgrid)
      ENDIF

      ! Map patch data to grid
      IF (p_is_io) THEN
         CALL allocate_block_data(this%y%grid, hx_wgrid, DEF_DA_ENS_NUM)
      ENDIF
      IF (p_is_worker) THEN
         CALL this%y%mg2p%pset2grid( &
            this%hx_model(1:DEF_DA_ENS_NUM,:), &
            hx_wgrid, spv=spval, msk=filter)
      ELSE
         dummy_pdata = spval
         CALL this%y%mg2p%pset2grid( &
            dummy_pdata, hx_wgrid, spv=spval)
      ENDIF
      IF (allocated(filter)) deallocate(filter)
      IF (.not. allocated(this%hx_obs)) THEN
         allocate(this%hx_obs(this%y%nobs,DEF_DA_ENS_NUM))
      ENDIF
      this%hx_obs = spval

      IF (p_is_io) THEN
         allocate(iloc(this%y%nobs))
         ndata = 0

         DO i = 1, this%y%nobs
            ! Read one-based indices in the observation reference grid
            ilat = this%y%ii(i)
            ilon = this%y%jj(i)
            IF (ilat < 1 .or. ilat > this%y%grid%nlat) CYCLE
            IF (ilon < 1 .or. ilon > this%y%grid%nlon) CYCLE

            ! Convert reference-grid indices to block indices
            ib = this%y%grid%xblk(ilon)
            jb = this%y%grid%yblk(ilat)
            IF (ib < 1 .or. ib > size(gblock%pio,1)) CYCLE
            IF (jb < 1 .or. jb > size(gblock%pio,2)) CYCLE

            ! Extract only blocks owned by this IO process
            IF (gblock%pio(ib,jb) /= p_iam_glb) CYCLE
            il = this%y%grid%xloc(ilon)
            jl = this%y%grid%yloc(ilat)
            IF (.not. allocated(area_wgrid%blk(ib,jb)%val)) CYCLE
            IF (.not. allocated(hx_wgrid%blk(ib,jb)%val)) CYCLE
            IF (il < 1 .or. il > size(area_wgrid%blk(ib,jb)%val,1)) CYCLE
            IF (jl < 1 .or. jl > size(area_wgrid%blk(ib,jb)%val,2)) CYCLE

            ! Normalize the mapped H(x) by the valid patch area
            area_obs = area_wgrid%blk(ib,jb)%val(il,jl)
            IF (area_obs <= 0.0_r8) CYCLE
            IF (any(hx_wgrid%blk(ib,jb)%val(:,il,jl) == spval)) CYCLE
            ndata = ndata + 1
            iloc(ndata) = i
#ifdef USEMPI
            this%hx_obs(ndata,:) = hx_wgrid%blk(ib,jb)%val(:,il,jl) / area_obs
#else
            this%hx_obs(i,:) = hx_wgrid%blk(ib,jb)%val(:,il,jl) / area_obs
#endif
         ENDDO

         ! Send IO-process results to master
#ifdef USEMPI
         smesg = (/p_iam_glb, ndata/)
         CALL mpi_send(smesg, 2, MPI_INTEGER, p_address_master, &
                       mpi_tag_mesg, p_comm_glb, p_err)
         IF (ndata > 0) THEN
            CALL mpi_send(iloc(1:ndata), ndata, MPI_INTEGER, &
                          p_address_master, mpi_tag_data, p_comm_glb, p_err)
            allocate(dtemp(ndata,DEF_DA_ENS_NUM))
            dtemp = this%hx_obs(1:ndata,:)
            CALL mpi_send(dtemp, ndata*DEF_DA_ENS_NUM, MPI_REAL8, &
                          p_address_master, mpi_tag_data+1, p_comm_glb, p_err)
            deallocate(dtemp)
         ENDIF
#endif
         deallocate(iloc)

         IF (allocated(hx_wgrid%blk)) deallocate(hx_wgrid%blk)
         IF (allocated(area_wgrid%blk)) deallocate(area_wgrid%blk)
      ENDIF

      ! Master: receive results from IO processes
#ifdef USEMPI
      IF (p_is_master) THEN
         this%hx_obs = spval
         DO ip = 0, p_np_io-1
            CALL mpi_recv(rmesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, &
                          mpi_tag_mesg, p_comm_glb, p_stat, p_err)
            ndata = rmesg(2)
            isrc = rmesg(1)
            IF (ndata > 0) THEN
               allocate(itemp(ndata))
               allocate(dtemp(ndata,DEF_DA_ENS_NUM))
               CALL mpi_recv(itemp, ndata, MPI_INTEGER, isrc, &
                             mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL mpi_recv(dtemp, ndata*DEF_DA_ENS_NUM, MPI_REAL8, isrc, &
                             mpi_tag_data+1, p_comm_glb, p_stat, p_err)
               this%hx_obs(itemp,:) = dtemp
               deallocate(itemp)
               deallocate(dtemp)
            ENDIF
         ENDDO
      ENDIF
      CALL mpi_bcast(this%hx_obs, this%y%nobs*DEF_DA_ENS_NUM, MPI_REAL8, &
                     p_address_master, p_comm_glb, p_err)
#endif
#endif

   END SUBROUTINE MR_forward

!-----------------------------------------------------------------------------

   SUBROUTINE MR_clear(this)

!-----------------------------------------------------------------------------
      class(MR_type), intent(inout)       :: this
!-----------------------------------------------------------------------------
      CALL this%y%clear()
      IF (allocated(this%hx_model)) deallocate(this%hx_model)
      IF (allocated(this%hx_obs  )) deallocate(this%hx_obs  )

   END SUBROUTINE MR_clear

!-----------------------------------------------------------------------------

   SUBROUTINE MR_set_init(this, pixelset, configs)

!-----------------------------------------------------------------------------
      class(MR_set_type), intent(inout)   :: this
      type(pixelset_type), intent(in)     :: pixelset
      type(MR_config_type), intent(in)    :: configs(:)
      integer                             :: is
!-----------------------------------------------------------------------------
      IF (allocated(this%sensors)) THEN
         DO is = 1, size(this%sensors)
            CALL this%sensors(is)%clear()
         ENDDO
         deallocate(this%sensors)
      ENDIF
      IF (allocated(this%cfgs)) deallocate(this%cfgs)
      this%nsensor = size(configs)
      this%nobs = 0
      IF (this%nsensor == 0) RETURN
      allocate(this%cfgs(this%nsensor), this%sensors(this%nsensor))
      this%cfgs = configs
      DO is = 1, this%nsensor
         CALL this%sensors(is)%init(pixelset, &
            this%cfgs(is)%source_name, this%cfgs(is)%sensor_name, &
            this%cfgs(is)%var_name)
      ENDDO

   END SUBROUTINE MR_set_init

!-----------------------------------------------------------------------------

   SUBROUTINE MR_set_clear(this)

!-----------------------------------------------------------------------------
      class(MR_set_type), intent(inout)   :: this
      integer                             :: is
!-----------------------------------------------------------------------------
      this%nobs = 0
      IF (allocated(this%lat   )) deallocate(this%lat   )
      IF (allocated(this%lon   )) deallocate(this%lon   )
      IF (allocated(this%y     )) deallocate(this%y     )
      IF (allocated(this%sigma )) deallocate(this%sigma )
      IF (allocated(this%hx_obs)) deallocate(this%hx_obs)
      IF (.not. allocated(this%sensors)) RETURN
      DO is = 1, this%nsensor
         CALL this%sensors(is)%clear()
      ENDDO

   END SUBROUTINE MR_set_clear

!-----------------------------------------------------------------------------

   SUBROUTINE MR_set_concat(this)

!-----------------------------------------------------------------------------
      class(MR_set_type), intent(inout)   :: this
      integer                             :: is, n, offset
!-----------------------------------------------------------------------------
      this%nobs = 0
      IF (allocated(this%lat   )) deallocate(this%lat   )
      IF (allocated(this%lon   )) deallocate(this%lon   )
      IF (allocated(this%y     )) deallocate(this%y     )
      IF (allocated(this%sigma )) deallocate(this%sigma )
      IF (allocated(this%hx_obs)) deallocate(this%hx_obs)
      DO is = 1, this%nsensor
         IF (.not. allocated(this%sensors(is)%hx_obs)) CYCLE
         this%nobs = this%nobs + this%sensors(is)%y%nobs
      ENDDO
      IF (this%nobs == 0) RETURN
      allocate(this%lat(this%nobs))
      allocate(this%lon(this%nobs))
      allocate(this%y(this%nobs))
      allocate(this%hx_obs(this%nobs, DEF_DA_ENS_NUM))
      allocate(this%sigma(this%nobs))

      ! Concatenate observations and matching H(x) in stream order.
      offset = 0
      DO is = 1, this%nsensor
         IF (.not. allocated(this%sensors(is)%hx_obs)) CYCLE
         n = this%sensors(is)%y%nobs
         IF (n == 0) CYCLE
         this%lat(offset+1:offset+n) = this%sensors(is)%y%lat
         this%lon(offset+1:offset+n) = this%sensors(is)%y%lon
         this%y(offset+1:offset+n) = this%sensors(is)%y%y
         this%sigma(offset+1:offset+n) = this%sensors(is)%y%sigma
         this%hx_obs(offset+1:offset+n,:) = this%sensors(is)%hx_obs
         offset = offset + n
      ENDDO

   END SUBROUTINE MR_set_concat

!-----------------------------------------------------------------------------

   SUBROUTINE MR_set_forward(this, idate, deltim)

!-----------------------------------------------------------------------------
      class(MR_set_type), intent(inout)   :: this
      integer,            intent(in)      :: idate(3)
      real(r8),           intent(in)      :: deltim

      integer                             :: is
!-----------------------------------------------------------------------------
      DO is = 1, this%nsensor
         CALL this%sensors(is)%forward(idate, deltim)
      ENDDO
      CALL this%concat()

   END SUBROUTINE MR_set_forward

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Data_MR
#endif
