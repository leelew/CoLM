#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Data_TQ
!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Observation handling for 2 m temperature and humidity data assimilation.
!
! WORKFLOW:
!   1. Init: Configure INSITU observation streams and initialize the source
!      module. Store observation search and localization parameters.
!   2. Forward: Call INSITU to read observations and compute matching H(x).
!      For assimilation, INSITU combines records from all configured streams.
!   3. Prepare: Select nearby valid observations for each target patch and
!      convert error standard deviations to localized variances for LETKF.
!   4. Save: Copy INSITU model-space H(x) into the shared diagnostic arrays.
!   5. Clear: Release current observations and H(x), retaining source
!      configurations and site-to-patch mappings for subsequent forward calls.
!
! SUPPORTED:
!   INSITU:
!      Station observations of 2 m temperature (2MT) and relative humidity
!      (2MRH), such as SYNOP records.
!
! NOTES:
!   Separate TQ_type instances serve temperature and humidity analysis
!   targets. MOD_DA_Assim_TQ controls their analysis order and state updates.
!
! HISTORY:
!   Lu Li, 11/2026: DA v1.0
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_LandPatch
   USE MOD_Namelist
   USE MOD_Pixelset
   USE MOD_SPMD_Task
   USE MOD_Vars_Global
   USE MOD_Vars_TimeInvariants
   USE MOD_DA_Data_INSITU
   USE MOD_DA_Vars_TimeVars
   IMPLICIT NONE

!-----------------------------------------------------------------------------
! Observation collection for one 2 m analysis target
!-----------------------------------------------------------------------------
   type :: TQ_type

      type(TQ_config_type), allocatable     :: cfgs(:)                        ! Observation stream configurations
      logical                               :: use_insitu = .false.           ! Include site observations

      ! Observation selection and LETKF parameters.
      real(r8)                              :: dres                           ! Observation search half-width [degrees]
      real(r8)                              :: locr                           ! Gaussian localization length scale [km]
      real(r8)                              :: infl                           ! Ensemble covariance inflation factor [-]

      ! Source-specific observation containers.
      type(INSITU_set_type)                 :: insitu                         ! Station observations and ensemble H(x)

   CONTAINS

      procedure                             :: init       => TQ_type_init     ! Configure sources and assimilation parameters
      procedure                             :: forward    => TQ_type_forward  ! Compute model- and observation-space H(x) from enabled sources
      procedure                             :: save       => TQ_type_save     ! Copy model-space H(x) into shared diagnostics
      procedure                             :: clear      => TQ_type_clear    ! Clear time-step data while retaining configurations
      procedure                             :: prepare    => TQ_type_prepare  ! Select nearby observations and localize error variances

   END type TQ_type

   type :: TQ_config_type

      character(len=16)                     :: source_name
      character(len=32)                     :: sensor_name
      character(len=16)                     :: var_name
      integer                               :: istream

   END type TQ_config_type

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE TQ_type_init(this, pixelset, configs)

!-----------------------------------------------------------------------------
      class(TQ_type), intent(inout)         :: this
      type(pixelset_type), intent(in)       :: pixelset
      type(TQ_config_type), intent(in)      :: configs(:)

      integer                               :: i
      integer                               :: iinsitu
      integer                               :: ninsitu
      type(INSITU_config_type), allocatable :: insitu_configs(:)
!-----------------------------------------------------------------------------
      CALL this%clear()
      this%use_insitu = .false.

      allocate(this%cfgs(size(configs)))
      this%cfgs = configs
      this%dres = DEF_DA_ASSIM_TQ_DRES
      this%locr = DEF_DA_ASSIM_TQ_LOCR
      this%infl = DEF_DA_ASSIM_TQ_INFL

      ! Configure the station source for site streams.
      ninsitu = count(configs%source_name == 'site')
      allocate(insitu_configs(ninsitu))
      iinsitu = 0
      DO i = 1, size(configs)
         IF (trim(configs(i)%source_name) /= 'site') CYCLE
         iinsitu = iinsitu + 1
         insitu_configs(iinsitu)%source_name = configs(i)%source_name
         insitu_configs(iinsitu)%sensor_name = configs(i)%sensor_name
         insitu_configs(iinsitu)%var_name = configs(i)%var_name
      ENDDO
      CALL this%insitu%init(pixelset, insitu_configs)
      this%use_insitu = ninsitu > 0

   END SUBROUTINE TQ_type_init

!-----------------------------------------------------------------------------

   SUBROUTINE TQ_type_forward(this, idate, deltim)

!-----------------------------------------------------------------------------
      class(TQ_type), intent(inout)         :: this
      integer, intent(in)                   :: idate(3)
      real(r8), intent(in)                  :: deltim
!-----------------------------------------------------------------------------
      IF (this%use_insitu) CALL this%insitu%forward(idate, deltim)

   END SUBROUTINE TQ_type_forward

!-----------------------------------------------------------------------------

   SUBROUTINE TQ_type_clear(this)

!-----------------------------------------------------------------------------
      class(TQ_type), intent(inout)         :: this
!-----------------------------------------------------------------------------
      CALL this%insitu%clear()

   END SUBROUTINE TQ_type_clear

!-----------------------------------------------------------------------------

   SUBROUTINE TQ_type_save(this, is_analysis)

!-----------------------------------------------------------------------------
      class(TQ_type), intent(in)            :: this
      logical, intent(in)                   :: is_analysis

      integer                               :: i
      integer                               :: is
      integer                               :: iinsitu
      integer                               :: np
      logical                               :: valid_hx
!-----------------------------------------------------------------------------
      iinsitu = 0
      DO i = 1, size(this%cfgs)
         IF (trim(this%cfgs(i)%source_name) /= 'site') CYCLE
         iinsitu = iinsitu + 1
         is = this%cfgs(i)%istream
         IF (is_analysis) THEN
            hx_a_ens(is,:,:) = spval
            hx_a(is,:) = spval
         ELSE
            hx_ol(is,:) = spval
            hx_f_ens(is,:,:) = spval
            hx_f(is,:) = spval
         ENDIF

         DO np = 1, numpatch
            IF (.not. DEF_DA) THEN
               hx_ol(is,np) = this%insitu%sensors(iinsitu)%hx_model(0,np)
               CYCLE
            ENDIF
            valid_hx = all(this%insitu%sensors(iinsitu)%hx_model(1:DEF_DA_ENS_NUM,np) /= spval)
            IF (is_analysis) THEN
               hx_a_ens(is,1:DEF_DA_ENS_NUM,np) = this%insitu%sensors(iinsitu)%hx_model(1:DEF_DA_ENS_NUM,np)
               IF (valid_hx) hx_a(is,np) = &
                  sum(this%insitu%sensors(iinsitu)%hx_model(1:DEF_DA_ENS_NUM,np)) / DEF_DA_ENS_NUM
            ELSE
               hx_ol(is,np) = this%insitu%sensors(iinsitu)%hx_model(0,np)
               hx_f_ens(is,1:DEF_DA_ENS_NUM,np) = this%insitu%sensors(iinsitu)%hx_model(1:DEF_DA_ENS_NUM,np)
               IF (valid_hx) hx_f(is,np) = &
                  sum(this%insitu%sensors(iinsitu)%hx_model(1:DEF_DA_ENS_NUM,np)) / DEF_DA_ENS_NUM
            ENDIF
         ENDDO
      ENDDO

   END SUBROUTINE TQ_type_save

!-----------------------------------------------------------------------------

   SUBROUTINE TQ_type_prepare(this, np, nobs_p, y_p, hx_p, r_p)

!-----------------------------------------------------------------------------
      class(TQ_type),         intent(in)    :: this
      integer,                intent(in)    :: np
      integer,                intent(out)   :: nobs_p
      real(r8), allocatable,  intent(out)   :: y_p(:)
      real(r8), allocatable,  intent(out)   :: hx_p(:,:)
      real(r8), allocatable,  intent(out)   :: r_p(:)

      integer                               :: iens, iobs
      real(r8)                              :: lat_p, lon_p, dlat, dlon
      logical,  allocatable                 :: valid(:)
      real(r8), allocatable                 :: lat(:), lon(:), dist(:), hav(:)
!-----------------------------------------------------------------------------
      nobs_p = 0
      lat_p = patchlatr(np)*180.0_r8/pi
      lon_p = patchlonr(np)*180.0_r8/pi

      ! All source records are globally available before this local search.
      allocate(valid(this%insitu%nobs))
      valid = .false.
      DO iobs = 1, this%insitu%nobs
         dlat = abs(this%insitu%lat(iobs) - lat_p)
         dlon = modulo(abs(this%insitu%lon(iobs) - lon_p), 360.0_r8)
         IF (dlon > 180.0_r8) dlon = 360.0_r8 - dlon
         valid(iobs) = dlat < this%dres .and. dlon < this%dres .and. &
            this%insitu%y(iobs) /= spval .and. all(this%insitu%hx_obs(iobs,:) /= spval)
      ENDDO
      nobs_p = count(valid)
      IF (nobs_p == 0) RETURN

      allocate(y_p(nobs_p), hx_p(nobs_p,DEF_DA_ENS_NUM), r_p(nobs_p))
      allocate(lat(nobs_p), lon(nobs_p), dist(nobs_p), hav(nobs_p))
      y_p = pack(this%insitu%y, valid)
      DO iens = 1, DEF_DA_ENS_NUM
         hx_p(:,iens) = pack(this%insitu%hx_obs(:,iens), valid)
      ENDDO
      r_p = pack(this%insitu%r, valid)
      lat = pack(this%insitu%lat, valid)
      lon = pack(this%insitu%lon, valid)

      ! Convert the error standard deviation to variance, then inflate R.
      hav = sin((lat*pi/180.0_r8 - patchlatr(np))/2.0_r8)**2 + &
         cos(lat*pi/180.0_r8)*cos(patchlatr(np))* &
         sin((lon*pi/180.0_r8 - patchlonr(np))/2.0_r8)**2
      dist = 2.0_r8*6378.1_r8*asin(sqrt(max(0.0_r8,min(1.0_r8,hav))))
      r_p = r_p**2 / max(exp(-(dist**2)/(2.0_r8*this%locr**2)), 1.0e-6_r8)

   END SUBROUTINE TQ_type_prepare

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Data_TQ
#endif
