#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Obs_INSITU
!-----------------------------------------------------------------------------
! DESCRIPTION:
!   In-situ observation handling for data assimilation.
!
!   This module defines one site observation stream (INSITU_type) and a
!   collection of streams (INSITU_set_type). Each stream represents one
!   source/network/variable configuration and owns the observation reader,
!   patch-space forward calculation, and extraction at observation locations.
!   Station networks such as SYNOP, ISMN, and SNOTEL share this structure;
!   the network name selects the dataset, while the variable selects H(x).
!
! WORKFLOW:
!   1. Init: Store the stream configuration. Initialize SiteType to read
!      the site table and build the site-to-patch mapping.
!   2. Forward: Read current observations through SiteType when needed.
!      Apply the selected operator to compute H(x) on local patches.
!      For assimilation, extract H(x) at the mapped site patches and share
!      observation-space H(x) across all processes.
!      For multiple streams, concatenate observations and matching H(x).
!   3. Clear: Release current observations and H(x), retaining the stream
!      configuration, site table, and site-to-patch mapping.
!
! SUPPORTED:
!   2MT:
!      2 m air temperature [K]. H(x) is the model 2 m air temperature.
!   2MQ:
!      2 m relative humidity [0-1]. H(x) is derived from model specific
!      humidity, temperature, and surface pressure.
!
! HISTORY:
!   Lu Li, 05/2026: Initial implementation
!   Lu Li, Zhilong Fan, 11/2026: DA v1.0
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
! In-situ observation derived type
!-----------------------------------------------------------------------------
   type :: INSITU_type

      real(r8), allocatable                 :: hx_model(:,:)                     ! Model-space H(x) on local patches
      type(site_obs_type)                   :: y                                 ! Observations, site table, and patch mapping
      real(r8), allocatable                 :: hx_obs(:,:)                       ! Observation-space ensemble H(x)

   CONTAINS

      procedure, PUBLIC                     :: init       => INSITU_init         ! Initialize the observation reader and patch mapping
      procedure, PUBLIC                     :: forward    => INSITU_forward      ! Compute patch H(x) and map it to observation locations
      procedure, PUBLIC                     :: clear      => INSITU_clear        ! Clear current observations and H(x), retaining the mapping

   END type INSITU_type

!-----------------------------------------------------------------------------
! In-situ observation stream configuration
!-----------------------------------------------------------------------------
   type :: INSITU_config_type

      character(len=16)                     :: source_name                       ! Observation source class
      character(len=32)                     :: sensor_name                       ! Sensor or station network name
      character(len=16)                     :: var_name                          ! Observation variable and file prefix

   END type INSITU_config_type

!-----------------------------------------------------------------------------
! Collection of in-situ observation streams and their combined data
!-----------------------------------------------------------------------------
   type :: INSITU_set_type

      integer                               :: nsensor = 0                       ! Number of configured observation streams
      type(INSITU_type), allocatable        :: sensors(:)                        ! Individual observation streams
      type(INSITU_config_type), allocatable :: cfgs(:)                           ! Configuration of each observation stream

      integer                               :: nobs = 0                          ! Total number of observations across streams
      real(r8), allocatable                 :: lat(:)                            ! Observation latitudes
      real(r8), allocatable                 :: lon(:)                            ! Observation longitudes
      real(r8), allocatable                 :: y(:)                              ! Observed values
      real(r8), allocatable                 :: r(:)                              ! Observation error standard deviations
      real(r8), allocatable                 :: hx_obs(:,:)                       ! Observation-space ensemble H(x)

   CONTAINS

      procedure, PUBLIC                     :: init       => INSITU_set_init     ! Initialize all observation streams
      procedure, PUBLIC                     :: concat     => INSITU_set_concat   ! Combine observations and matching H(x) in stream order
      procedure, PUBLIC                     :: forward    => INSITU_set_forward  ! Run forward calculations for all streams and combine their records
      procedure, PUBLIC                     :: clear      => INSITU_set_clear    ! Clear current observations and H(x), retaining stream configurations

   END type INSITU_set_type

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE INSITU_init(this, pixelset, source_name, sensor_name, var_name)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(INSITU_type), intent(inout) :: this
      type(pixelset_type), intent(in)   :: pixelset
      character(len=*), intent(in)      :: source_name
      character(len=*), intent(in)      :: sensor_name
      character(len=*), intent(in)      :: var_name

      CALL this%y%init(pixelset, source_name, sensor_name, var_name)

   END SUBROUTINE INSITU_init

!-----------------------------------------------------------------------------

   SUBROUTINE INSITU_forward(this, idate, deltim)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------------ Dummy Arguments ------------------------------
      class(INSITU_type), intent(inout) :: this
      integer, intent(in)               :: idate(3)
      real(r8), intent(in)              :: deltim

!------------------------------ Local Variables ------------------------------
      integer                           :: iobs      ! Observation index
      integer                           :: isite     ! Site-table index
      integer                           :: np        ! Local patch index
      integer                           :: iens      ! Ensemble member index
      real(r8)                          :: tc        ! 2 m air temperature [degrees Celsius]
      real(r8)                          :: es        ! Saturation vapor pressure [hPa]
      real(r8)                          :: e         ! Water vapor pressure [hPa]
      integer, allocatable              :: owner(:)  ! Number of workers contributing to each observation

!-----------------------------------------------------------------------------
      CALL this%clear()

      ! Reuse observations until clear is called at the end of the time step.
      CALL this%y%read(idate, deltim)
      IF (this%y%nobs == 0) RETURN

      ! Compute local patch H(x); all ranks then participate in observation mapping.
      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            allocate(this%hx_model(0:DEF_DA_ENS_NUM, numpatch))
            this%hx_model = spval

            SELECT CASE (trim(this%y%var_name))
!-----------------------------------------------------------------------------
! 2 m air temperature observations [K].
!-----------------------------------------------------------------------------
            CASE ('2MT')
               this%hx_model = tref_ens

!-----------------------------------------------------------------------------
! 2 m relative humidity observations [0-1].
!-----------------------------------------------------------------------------
            CASE ('2MQ')
               ! Convert specific humidity to relative humidity [0-1].
               DO np = 1, numpatch
                  DO iens = 0, DEF_DA_ENS_NUM
                     rhref_ens(iens,np) = spval

                     IF (tref_ens(iens,np) == spval .or.                                   &
                        qref_ens(iens,np) == spval .or.                                    &
                        forc_psrf(np)    == spval .or.                                     &
                        forc_psrf(np)    <= 0.0_r8) CYCLE

                     tc = tref_ens(iens,np) - 273.15_r8
                     es = 6.112_r8 * exp((17.67_r8*tc)/(tc + 243.5_r8))
                     e  = (qref_ens(iens,np)*forc_psrf(np)/100.0_r8) /                     &
                        (0.378_r8*qref_ens(iens,np) + 0.622_r8)
                     IF (es > 0.0_r8) THEN
                        rhref_ens(iens,np) = max(min(1.0_r8, e/es), 0.0_r8)
                     ENDIF
                  ENDDO
               ENDDO
               this%hx_model = rhref_ens

            END SELECT
         ENDIF
      ENDIF

      IF (.not. DEF_DA) RETURN

      ! Extract H(x) at the mapped station patches.
      allocate (this%hx_obs(this%y%nobs, 1:DEF_DA_ENS_NUM))
      this%hx_obs(:,:) = 0.0_r8
      allocate (owner(this%y%nobs))
      owner(:) = 0

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
      CALL mpi_allreduce(MPI_IN_PLACE, this%hx_obs, this%y%nobs*DEF_DA_ENS_NUM, MPI_REAL8, &
         MPI_SUM, p_comm_glb, p_err)
      CALL mpi_allreduce(MPI_IN_PLACE, owner, this%y%nobs, MPI_INTEGER,                    &
         MPI_SUM, p_comm_glb, p_err)
#endif

      ! Mark observations without exactly one contributing worker as missing.
      DO iobs = 1, this%y%nobs
         IF (owner(iobs) /= 1) this%hx_obs(iobs, 1:DEF_DA_ENS_NUM) = spval
      ENDDO

      deallocate (owner)

   END SUBROUTINE INSITU_forward

!-----------------------------------------------------------------------------

   SUBROUTINE INSITU_clear(this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(INSITU_type), intent(inout) :: this

      CALL this%y%clear()
      IF (allocated(this%hx_model)) deallocate(this%hx_model)
      IF (allocated(this%hx_obs  )) deallocate(this%hx_obs  )

   END SUBROUTINE INSITU_clear

!-----------------------------------------------------------------------------

   SUBROUTINE INSITU_set_init(this, pixelset, configs)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(INSITU_set_type),    intent(inout) :: this
      type(pixelset_type),      intent(in)     :: pixelset
      type(INSITU_config_type),  intent(in)    :: configs(:)
      integer                                  :: i

      CALL this%clear()
      IF (allocated(this%cfgs   )) deallocate(this%cfgs   )
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

   END SUBROUTINE INSITU_set_init

!-----------------------------------------------------------------------------

   SUBROUTINE INSITU_set_clear(this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(INSITU_set_type), intent(inout) :: this
      integer                               :: i

      this%nobs = 0
      IF (allocated(this%lat   )) deallocate(this%lat   )
      IF (allocated(this%lon   )) deallocate(this%lon   )
      IF (allocated(this%y     )) deallocate(this%y     )
      IF (allocated(this%hx_obs)) deallocate(this%hx_obs)
      IF (allocated(this%r     )) deallocate(this%r     )

      IF (allocated(this%sensors)) THEN
         DO i = 1, size(this%sensors)
            CALL this%sensors(i)%clear()
         ENDDO
      ENDIF

   END SUBROUTINE INSITU_set_clear

!-----------------------------------------------------------------------------

   SUBROUTINE INSITU_set_forward(this, idate, deltim)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(INSITU_set_type), intent(inout) :: this
      integer, intent(in)                   :: idate(3)
      real(r8), intent(in)                  :: deltim
      integer                               :: is

      IF (.not. allocated(this%sensors)) RETURN
      DO is = 1, this%nsensor
         CALL this%sensors(is)%forward(idate, deltim)
      ENDDO
      IF (DEF_DA) CALL this%concat()

   END SUBROUTINE INSITU_set_forward

!-----------------------------------------------------------------------------

   SUBROUTINE INSITU_set_concat(this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(INSITU_set_type), intent(inout) :: this
      integer                               :: is, n, offset

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

   END SUBROUTINE INSITU_set_concat

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Obs_INSITU
#endif
