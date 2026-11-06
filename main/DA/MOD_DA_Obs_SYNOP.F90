#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Obs_SYNOP
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    SYNOP observation handling for 2 m temperature and humidity assimilation.
!
!    SYNOP_type represents one source/sensor/variable stream. SYNOP_set_type
!    applies the observation workflow to all configured streams and combines
!    their observations, metadata, errors, and H(x) for assimilation.
!
!    The observation-space workflow is:
!      1. clear time-step data and read observations through site_obs_type;
!      2. compute model-space H(x) on local patches;
!      3. extract H(x) at the mapped station patches;
!      4. share observation-space H(x) across MPI ranks;
!      5. concatenate SYNOP streams in configuration order.
!
!    Temperature H(x) is tref_ens; humidity H(x) is relative humidity derived
!    from qref_ens, tref_ens, and surface pressure. The assimilation driver
!    analyzes temperature before humidity, so the humidity operator uses the
!    updated temperature.
!
!    Station coordinates are read from <DEF_DA_OBS_DIR>/<sensor>/sites.nc.
!    Hourly files <var_name>_YYYY_MM_DD_HH.nc contain 'latitude', 'longitude',
!    'time', 'id', <var_name>, and 'r'. Variables '2MT' and '2MQ' denote
!    temperature [K] and relative humidity [0-1]. Observation errors 'r' are
!    standard deviations; the upper observation collection handles localization.
!
! HISTORY:
!    Lu Li, 05/2026: Initial implementation
!    Lu Li, Zhilong Fan, 11/2026: DA v1.0
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_LandPatch
   USE MOD_Namelist
   USE MOD_Pixelset
   USE MOD_Vars_Global
   USE MOD_Vars_1DForcing
   USE MOD_SPMD_Task
   USE MOD_DA_Obs_SiteType
   USE MOD_DA_Vars_TimeVars
   IMPLICIT NONE

!-----------------------------------------------------------------------------
! One SYNOP observation stream
!-----------------------------------------------------------------------------
   type :: SYNOP_type

      ! Model-space H(x) on local patches, retained for diagnostic output.
      real(r8), allocatable :: hx_model(:,:)                  ! [0:ens, patch]

      ! Observations and metadata at the locations read for this time step.
      type(site_obs_type)   :: y

      ! Model equivalents at observation locations; unmapped values are spval.
      real(r8), allocatable :: hx_obs(:,:)                    ! [nobs, ens]

   CONTAINS

      procedure, PUBLIC :: init       => SYNOP_init           ! Initialize the observation reader and patch mapping
      procedure, PUBLIC :: calc_model => SYNOP_calc_model     ! Compute model-space H(x) on local patches
      procedure, PUBLIC :: calc_obs   => SYNOP_calc_obs       ! Read y, compute H(x), and map to observation space
      procedure, PUBLIC :: clear      => SYNOP_clear          ! Clear time-step data while retaining the mapping

   END type SYNOP_type

!-----------------------------------------------------------------------------
! SYNOP observation stream configuration
!-----------------------------------------------------------------------------
   type :: SYNOP_config_type

      character(len=16) :: source_name                        ! Observation source class
      character(len=32) :: sensor_name                        ! Sensor or station network name
      character(len=16) :: var_name                           ! Observation variable and file prefix

   END type SYNOP_config_type

!-----------------------------------------------------------------------------
! Collection of SYNOP observation streams and their combined data
!-----------------------------------------------------------------------------
   type :: SYNOP_set_type

      integer :: nsensor = 0                                  ! Number of configured observation streams
      type(SYNOP_type), allocatable :: sensors(:)             ! Stream objects [nsensor]
      type(SYNOP_config_type), allocatable :: cfgs(:)         ! Stream configurations [nsensor]

      ! Current time-step data concatenated in stream order for assimilation.
      integer :: nobs = 0                                     ! Total number of observations
      real(r8), allocatable :: lat(:)                         ! Observation latitude [nobs], degrees
      real(r8), allocatable :: lon(:)                         ! Observation longitude [nobs], degrees
      real(r8), allocatable :: y(:)                           ! Observation values [nobs]
      real(r8), allocatable :: r(:)                           ! Observation error standard deviations [nobs]
      real(r8), allocatable :: hx_obs(:,:)                    ! Model equivalents at y locations [nobs, ens]

   CONTAINS

      procedure, PUBLIC :: init       => SYNOP_set_init       ! Initialize all streams from their configurations
      procedure, PUBLIC :: calc_model => SYNOP_set_calc_model ! Recompute model-space H(x) for all streams
      procedure, PUBLIC :: calc_obs   => SYNOP_set_calc_obs   ! Prepare observation-space y and H(x), then combine streams
      procedure, PUBLIC :: concat     => SYNOP_set_concat     ! Combine stream data into the set-level arrays
      procedure, PUBLIC :: clear      => SYNOP_set_clear      ! Clear time-step data while retaining configurations

   END type SYNOP_set_type

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE SYNOP_init (this, pixelset, source_name, sensor_name, var_name)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(SYNOP_type), intent(inout) :: this
      type(pixelset_type), intent(in) :: pixelset
      character(len=*), intent(in) :: source_name
      character(len=*), intent(in) :: sensor_name
      character(len=*), intent(in) :: var_name

!-----------------------------------------------------------------------------
      CALL this%y%init(pixelset, source_name, sensor_name, var_name)

   END SUBROUTINE SYNOP_init

!-----------------------------------------------------------------------------

   SUBROUTINE SYNOP_calc_model (this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(SYNOP_type), intent(inout) :: this

!------------------------ Local Variables ------------------------------------
      integer :: np                                           ! Local patch index
      integer :: iens                                         ! Ensemble member index
      real(r8) :: tc                                          ! 2 m air temperature [degrees Celsius]
      real(r8) :: es                                          ! Saturation vapor pressure [hPa]
      real(r8) :: e                                           ! Water vapor pressure [hPa]

!-----------------------------------------------------------------------------
      IF (.not. p_is_worker) RETURN
      IF (allocated(this%hx_model)) deallocate(this%hx_model)
      allocate(this%hx_model(0:DEF_DA_ENS_NUM, numpatch))
      this%hx_model = spval
      IF (numpatch == 0) RETURN

      SELECT CASE (trim(this%y%var_name))
!-----------------------------------------------------------------------------
! 2 m air temperature observations [K].
! H(x) is the model 2 m air temperature.
!-----------------------------------------------------------------------------
      CASE ('2MT')
         this%hx_model = tref_ens

!-----------------------------------------------------------------------------
! 2 m relative humidity observations [0-1].
! H(x) is derived from specific humidity, temperature, and surface pressure.
!-----------------------------------------------------------------------------
      CASE ('2MQ')
         ! Convert specific humidity to relative humidity [0-1].
         DO np = 1, numpatch
            DO iens = 0, merge(DEF_DA_ENS_NUM, 0, DEF_DA)
               rhref_ens(iens,np) = spval

               IF (tref_ens(iens,np) == spval .or. &
                   qref_ens(iens,np) == spval .or. &
                   forc_psrf(np)    == spval .or. &
                   forc_psrf(np)    <= 0.0_r8) CYCLE

               tc = tref_ens(iens,np) - 273.15_r8
               es = 6.112_r8 * exp((17.67_r8*tc)/(tc + 243.5_r8))
               e  = (qref_ens(iens,np)*forc_psrf(np)/100.0_r8) / &
                    (0.378_r8*qref_ens(iens,np) + 0.622_r8)
               IF (es > 0.0_r8) THEN
                  rhref_ens(iens,np) = max(min(1.0_r8, e/es), 0.0_r8)
               ENDIF
            ENDDO
         ENDDO
         this%hx_model = rhref_ens

!-----------------------------------------------------------------------------
! snow depth.
!-----------------------------------------------------------------------------
      CASE ('SD')
         RETURN
      END SELECT

   END SUBROUTINE SYNOP_calc_model

!-----------------------------------------------------------------------------

   SUBROUTINE SYNOP_calc_obs (this, idate, deltim)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(SYNOP_type), intent(inout) :: this
      integer, intent(in) :: idate(3)
      real(r8), intent(in) :: deltim

!------------------------ Local Variables ------------------------------------
      integer :: iobs                                         ! Observation index
      integer :: isite                                        ! Site-table index
      integer :: np                                           ! Local patch index
      integer, allocatable :: owner(:)                        ! Number of workers contributing to each observation

!-----------------------------------------------------------------------------
      CALL this%clear()

      ! Read current observations.
      CALL this%y%read(idate, deltim)

      ! Compute model-space H(x).
      CALL this%calc_model()
      IF (.not. DEF_DA) RETURN
      allocate (this%hx_obs(this%y%nobs, 1:DEF_DA_ENS_NUM))
      this%hx_obs(:,:) = spval
      IF (this%y%nobs == 0) RETURN

      ! Extract H(x) at the mapped station patches.
      allocate (owner(this%y%nobs))
      owner(:) = 0
      this%hx_obs(:,:) = 0.0_r8

      IF (p_is_worker .and. numpatch > 0 .and. allocated(this%hx_model)) THEN
         IF (size(this%hx_model, 2) >= numpatch) THEN
            DO iobs = 1, this%y%nobs
               isite = this%y%isite(iobs)
               IF (this%y%site_iworker(isite) == p_iam_glb) THEN
                  np = this%y%site_iworker_ipatch(isite)
                  IF (np >= 1 .and. np <= numpatch) THEN
                     this%hx_obs(iobs, 1:DEF_DA_ENS_NUM) = this%hx_model(1:DEF_DA_ENS_NUM, np)
                     owner(iobs) = 1
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
      ENDIF

      ! Share observation-space H(x) across all ranks.
#ifdef USEMPI
      CALL mpi_allreduce (MPI_IN_PLACE, this%hx_obs, this%y%nobs*DEF_DA_ENS_NUM, MPI_REAL8, &
         MPI_SUM, p_comm_glb, p_err)
      CALL mpi_allreduce (MPI_IN_PLACE, owner, this%y%nobs, MPI_INTEGER, &
         MPI_SUM, p_comm_glb, p_err)
#endif

      ! Mark observations without exactly one contributing worker as missing.
      DO iobs = 1, this%y%nobs
         IF (owner(iobs) /= 1) this%hx_obs(iobs, 1:DEF_DA_ENS_NUM) = spval
      ENDDO

      deallocate (owner)

   END SUBROUTINE SYNOP_calc_obs

!-----------------------------------------------------------------------------

   SUBROUTINE SYNOP_clear (this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(SYNOP_type), intent(inout) :: this

!-----------------------------------------------------------------------------
      CALL this%y%clear()
      IF (allocated(this%hx_model)) deallocate(this%hx_model)
      IF (allocated(this%hx_obs)) deallocate(this%hx_obs)

   END SUBROUTINE SYNOP_clear

!-----------------------------------------------------------------------------

   SUBROUTINE SYNOP_set_init (this, pixelset, configs)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(SYNOP_set_type),    intent(inout) :: this
      type(pixelset_type),      intent(in)    :: pixelset
      type(SYNOP_config_type),  intent(in)    :: configs(:)

!------------------------ Local Variables ------------------------------------
      integer :: i

!-----------------------------------------------------------------------------
      CALL this%clear()
      IF (allocated(this%cfgs)) deallocate(this%cfgs)
      IF (allocated(this%sensors)) deallocate(this%sensors)
      this%nsensor = 0

      IF (size(configs) == 0) RETURN

      allocate(this%cfgs(size(configs)))
      this%cfgs = configs

      ! Initialize one observation reader and operator per stream.
      this%nsensor = size(configs)
      allocate(this%sensors(this%nsensor))
      DO i = 1, this%nsensor
         CALL this%sensors(i)%init(pixelset, configs(i)%source_name, &
            configs(i)%sensor_name, configs(i)%var_name)
      ENDDO

   END SUBROUTINE SYNOP_set_init

!-----------------------------------------------------------------------------

   SUBROUTINE SYNOP_set_clear (this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(SYNOP_set_type), intent(inout) :: this

!------------------------ Local Variables ------------------------------------
      integer :: i

!-----------------------------------------------------------------------------
      this%nobs = 0
      IF (allocated(this%lat)) deallocate(this%lat)
      IF (allocated(this%lon)) deallocate(this%lon)
      IF (allocated(this%y))   deallocate(this%y)
      IF (allocated(this%hx_obs))  deallocate(this%hx_obs)
      IF (allocated(this%r))   deallocate(this%r)

      IF (allocated(this%sensors)) THEN
         DO i = 1, size(this%sensors)
            CALL this%sensors(i)%clear()
         ENDDO
      ENDIF

   END SUBROUTINE SYNOP_set_clear

!-----------------------------------------------------------------------------

   SUBROUTINE SYNOP_set_calc_model (this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(SYNOP_set_type), intent(inout) :: this

!------------------------ Local Variables ------------------------------------
      integer :: is

!-----------------------------------------------------------------------------
      IF (.not. allocated(this%sensors)) RETURN
      DO is = 1, this%nsensor
         CALL this%sensors(is)%calc_model()
      ENDDO

   END SUBROUTINE SYNOP_set_calc_model

!-----------------------------------------------------------------------------

   SUBROUTINE SYNOP_set_calc_obs (this, idate, deltim)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(SYNOP_set_type), intent(inout) :: this
      integer, intent(in) :: idate(3)
      real(r8), intent(in) :: deltim

!------------------------ Local Variables ------------------------------------
      integer :: is

!-----------------------------------------------------------------------------
      CALL this%clear()
      IF (.not. allocated(this%sensors)) RETURN
      DO is = 1, this%nsensor
         CALL this%sensors(is)%calc_obs(idate, deltim)
      ENDDO
      IF (DEF_DA) CALL this%concat()

   END SUBROUTINE SYNOP_set_calc_obs

!-----------------------------------------------------------------------------

   SUBROUTINE SYNOP_set_concat (this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(SYNOP_set_type), intent(inout) :: this
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

   END SUBROUTINE SYNOP_set_concat

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Obs_SYNOP
#endif
