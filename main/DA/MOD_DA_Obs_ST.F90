#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Obs_ST
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Observation aggregation module for soil-temperature data assimilation.
!
!    This module provides ST_type, a soil-temperature-specific observation
!    container used by MOD_DA_Assim_ST. It sits above source-specific
!    observation modules, currently PM.
!
!    ST_type owns the soil-temperature observation workflow:
!      1. receive a unified list of observation configs from the assimilation
!         driver;
!      2. split configs by source and initialize the corresponding source
!         modules, e.g., PM_set_type;
!      3. call each enabled source to read observations and calculate ensemble
!         H(x);
!      4. concatenate all enabled source/sensor observations into one vector;
!      5. prepare localized patch-level y, R, and H(x) arrays for LETKF.
!
!    Source-level parameters, such as PM sensor, variable, and frequency,
!    are passed to the corresponding source module config type.
!    Soil-temperature assimilation parameters, such as observation search
!    radius and localization radius, are kept in ST_type.
!
! HISTORY:
!    Lu Li, 05/2026: First implementation
!    Lu Li, Zhilong Fan, 11/2026: DA v1.0
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_DA_Obs_PM
   USE MOD_DA_Vars_TimeVars
   USE MOD_Namelist
   USE MOD_Pixelset
   USE MOD_Vars_Global
   USE MOD_Vars_TimeInvariants
   IMPLICIT NONE

!-----------------------------------------------------------------------------
! Configuration for one observation stream used in soil temperature assimilation.
!-----------------------------------------------------------------------------
   type :: ST_config_type

      character(len=16) :: source_name                  ! Observation source, currently 'PM'
      character(len=32) :: sensor_name                  ! Sensor or dataset identifier used by the source module
      character(len=16) :: var_name                     ! Observation variable passed to the source operator
      real(r8) :: fghz                                  ! Microwave frequency for the PM brightness temperature operator [GHz]
      integer  :: istream                               ! Stream index in the shared hx_ol, hx_f* and hx_a* diagnostic arrays

   END type ST_config_type

!-----------------------------------------------------------------------------
! Observation collection for soil temperature assimilation.
! Source modules prepare observations and ensemble H(x); this type combines
! their records and selects nearby observations for each target patch.
!-----------------------------------------------------------------------------
   type :: ST_type

      type(ST_config_type), allocatable :: cfgs(:)      ! Configuration for each observation stream

      logical :: use_pm = .false.                       ! Include passive microwave observations

      ! Observation selection and LETKF parameters.
      real(r8) :: dres                                  ! Observation search half-width in latitude and longitude [degrees]
      real(r8) :: locr                                  ! Gaussian localization length scale for distance-based R inflation [km]
      real(r8) :: infl                                  ! Ensemble background covariance inflation factor

      ! Source-specific observation containers.
      type(PM_set_type) :: pm                           ! Passive microwave observations and ensemble H(x)

      ! Combined records on all ranks; units of y, hx_obs and r depend on the variable.
      integer :: nobs = 0                               ! Total number of records from all enabled sources
      real(r8), allocatable :: lat(:)                   ! Observation latitudes [nobs], degrees north
      real(r8), allocatable :: lon(:)                   ! Observation longitudes [nobs], degrees east
      real(r8), allocatable :: y(:)                     ! Observed values [nobs]
      real(r8), allocatable :: hx_obs(:,:)              ! Forecast observation equivalents [nobs, DEF_DA_ENS_NUM]
      real(r8), allocatable :: r(:)                     ! Observation error standard deviations [nobs]

   CONTAINS

      procedure :: init       => ST_type_init           ! Configure enabled sources and assimilation parameters
      procedure :: calc_obs   => ST_type_calc_obs       ! Prepare observation-space y and H(x) from enabled sources
      procedure :: calc_model => ST_type_calc_model     ! Recompute PM model-space H(x) after state updates
      procedure :: save       => ST_type_save           ! Copy PM patch H(x) into forecast or analysis diagnostics
      procedure :: clear      => ST_type_clear          ! Deallocate combined records and clear source observations
      procedure :: concat     => ST_type_concat         ! Combine enabled sources in matching y, hx_obs and r order
      procedure :: prepare    => ST_type_prepare        ! Select nearby valid records and localize their error variances


   END type ST_type

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE ST_type_init(this, pixelset, configs)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(ST_type),       intent(inout) :: this
      type(pixelset_type),  intent(in)    :: pixelset
      type(ST_config_type), intent(in)    :: configs(:)

!------------------------ Local Variables ------------------------------------
      integer :: i, ipm, npm
      type(PM_config_type), allocatable  :: pm_configs(:)

!-----------------------------------------------------------------------------
      CALL this%clear()
      this%use_pm = .false.

      IF (allocated(this%cfgs)) deallocate(this%cfgs)
      IF (size(configs) == 0) RETURN

      allocate(this%cfgs(size(configs)))
      this%cfgs = configs
      this%dres = DEF_DA_ASSIM_ST_DRES
      this%locr = DEF_DA_ASSIM_ST_LOCR
      this%infl = DEF_DA_ASSIM_ST_INFL

      ! the number of each source type
      npm = 0
      DO i = 1, size(configs)
         IF (trim(configs(i)%source_name) == 'PM') npm = npm + 1
      ENDDO

      ! initialize each source type module
      IF (npm > 0) THEN
         this%use_pm = .true.
         allocate(pm_configs(npm))

         ipm = 0
         DO i = 1, size(configs)
            IF (trim(configs(i)%source_name) == 'PM') THEN
               ipm = ipm + 1
               pm_configs(ipm)%source_name = configs(i)%source_name
               pm_configs(ipm)%sensor_name = configs(i)%sensor_name
               pm_configs(ipm)%var_name    = configs(i)%var_name
               pm_configs(ipm)%fghz        = configs(i)%fghz
            ENDIF
         ENDDO

         CALL this%pm%init(pixelset, pm_configs)

         deallocate(pm_configs)
      ENDIF

   END SUBROUTINE ST_type_init

!-----------------------------------------------------------------------------

   SUBROUTINE ST_type_clear(this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(ST_type), intent(inout) :: this

!-----------------------------------------------------------------------------
      this%nobs = 0
      IF (allocated(this%lat   )) deallocate(this%lat   )
      IF (allocated(this%lon   )) deallocate(this%lon   )
      IF (allocated(this%y     )) deallocate(this%y     )
      IF (allocated(this%hx_obs)) deallocate(this%hx_obs)
      IF (allocated(this%r     )) deallocate(this%r     )

      CALL this%pm%clear()

   END SUBROUTINE ST_type_clear

!-----------------------------------------------------------------------------

   SUBROUTINE ST_type_concat(this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(ST_type), intent(inout) :: this

!------------------------ Local Variables ------------------------------------
      integer :: n
      integer :: offset

!-----------------------------------------------------------------------------
      this%nobs = 0
      IF (allocated(this%lat   )) deallocate(this%lat   )
      IF (allocated(this%lon   )) deallocate(this%lon   )
      IF (allocated(this%y     )) deallocate(this%y     )
      IF (allocated(this%hx_obs)) deallocate(this%hx_obs)
      IF (allocated(this%r     )) deallocate(this%r     )

      IF (this%use_pm) this%nobs = this%nobs + this%pm%nobs
      IF (this%nobs == 0) RETURN

      allocate(this%lat(this%nobs))
      allocate(this%lon(this%nobs))
      allocate(this%y(this%nobs))
      allocate(this%hx_obs(this%nobs, DEF_DA_ENS_NUM))
      allocate(this%r(this%nobs))

      offset = 0

      IF (this%use_pm) THEN
         n = this%pm%nobs
         IF (n > 0) THEN
            this%lat(offset+1:offset+n) = this%pm%lat
            this%lon(offset+1:offset+n) = this%pm%lon
            this%y(offset+1:offset+n) = this%pm%y
            this%hx_obs(offset+1:offset+n,:) = this%pm%hx_obs
            this%r(offset+1:offset+n) = this%pm%r
            offset = offset + n
         ENDIF
      ENDIF

   END SUBROUTINE ST_type_concat

!-----------------------------------------------------------------------------

   SUBROUTINE ST_type_calc_obs(this, idate, deltim)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(ST_type), intent(inout) :: this
      integer,        intent(in)    :: idate(3)
      real(r8),       intent(in)    :: deltim

!-----------------------------------------------------------------------------
      CALL this%clear()

      IF (this%use_pm) CALL this%pm%calc_obs(idate, deltim)

      IF (DEF_DA) CALL this%concat()

   END SUBROUTINE ST_type_calc_obs

!-----------------------------------------------------------------------------

   SUBROUTINE ST_type_calc_model(this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(ST_type), intent(inout) :: this

!-----------------------------------------------------------------------------
      IF (this%use_pm) CALL this%pm%calc_model()

   END SUBROUTINE ST_type_calc_model

!-----------------------------------------------------------------------------

   SUBROUTINE ST_type_save (this, is_analysis)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(ST_type), intent(in) :: this
      logical,        intent(in) :: is_analysis

!------------------------ Local Variables ------------------------------------
      integer :: istream
      integer :: is
      integer :: ipm
      integer :: np
      integer :: iens
      logical :: valid_hx

!-----------------------------------------------------------------------------
      IF (nsource <= 0) RETURN
      IF (numpatch <= 0) RETURN
      IF (.not. allocated(hx_f_ens)) RETURN

      DO istream = 1, size(this%cfgs)
         is = this%cfgs(istream)%istream
         IF (is_analysis) THEN
            hx_a_ens(is,:,:) = spval
            hx_a    (is,:)   = spval
         ELSE
            hx_ol   (is,:)   = spval
            hx_f_ens(is,:,:) = spval
            hx_f    (is,:)   = spval
         ENDIF
      ENDDO

      ipm = 0
      DO istream = 1, size(this%cfgs)
         IF (trim(this%cfgs(istream)%source_name) == 'PM') THEN
            ipm = ipm + 1
            is = this%cfgs(istream)%istream

            IF (is_analysis) THEN
               DO np = 1, numpatch
                  valid_hx = all(this%pm%sensors(ipm)%hx_model(1:DEF_DA_ENS_NUM,np) /= spval)
                  DO iens = 1, DEF_DA_ENS_NUM
                     hx_a_ens(is,iens,np) = this%pm%sensors(ipm)%hx_model(iens,np)
                  ENDDO
                  IF (valid_hx) hx_a(is,np) = &
                     sum(this%pm%sensors(ipm)%hx_model(1:DEF_DA_ENS_NUM,np)) / DEF_DA_ENS_NUM
               ENDDO
            ELSE
               DO np = 1, numpatch
                  hx_ol(is,np) = this%pm%sensors(ipm)%hx_model(0,np)
                  IF (.not. DEF_DA) CYCLE
                  valid_hx = all(this%pm%sensors(ipm)%hx_model(1:DEF_DA_ENS_NUM,np) /= spval)
                  DO iens = 1, DEF_DA_ENS_NUM
                     hx_f_ens(is,iens,np) = this%pm%sensors(ipm)%hx_model(iens,np)
                  ENDDO
                  IF (valid_hx) hx_f(is,np) = &
                     sum(this%pm%sensors(ipm)%hx_model(1:DEF_DA_ENS_NUM,np)) / DEF_DA_ENS_NUM
               ENDDO
            ENDIF
         ENDIF
      ENDDO

   END SUBROUTINE ST_type_save

!-----------------------------------------------------------------------------

   SUBROUTINE ST_type_prepare (this, np, nobs_p, y_p, hx_p, r_p)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(ST_type),         intent(in)  :: this
      integer,                intent(in)  :: np
      integer,                intent(out) :: nobs_p
      real(r8), allocatable,  intent(out) :: y_p(:)
      real(r8), allocatable,  intent(out) :: hx_p(:,:)
      real(r8), allocatable,  intent(out) :: r_p(:)

!------------------------ Local Variables ------------------------------------
      integer :: iens, iobs
      real(r8) :: lat_p, lon_p, dlat, dlon
      logical,  allocatable :: valid(:)
      real(r8), allocatable :: lat(:), lon(:), dist(:), hav(:)

!-----------------------------------------------------------------------------
      nobs_p = 0
      lat_p = patchlatr(np)*180.0_r8/pi
      lon_p = patchlonr(np)*180.0_r8/pi

      ! All source records are globally available before this local search.
      allocate(valid(this%nobs))
      valid = .false.
      DO iobs = 1, this%nobs
         dlat = abs(this%lat(iobs) - lat_p)
         dlon = modulo(abs(this%lon(iobs) - lon_p), 360.0_r8)
         IF (dlon > 180.0_r8) dlon = 360.0_r8 - dlon
         valid(iobs) = dlat < this%dres .and. dlon < this%dres .and. &
                       this%y(iobs) /= spval .and. all(this%hx_obs(iobs,:) /= spval)
      ENDDO
      nobs_p = count(valid)
      IF (nobs_p == 0) RETURN

      allocate(y_p(nobs_p), hx_p(nobs_p,DEF_DA_ENS_NUM), r_p(nobs_p))
      allocate(lat(nobs_p), lon(nobs_p), dist(nobs_p), hav(nobs_p))
      y_p = pack(this%y, valid)
      DO iens = 1, DEF_DA_ENS_NUM
         hx_p(:,iens) = pack(this%hx_obs(:,iens), valid)
      ENDDO
      r_p = pack(this%r, valid)
      lat = pack(this%lat, valid)
      lon = pack(this%lon, valid)

      ! Convert the error standard deviation to variance, then inflate R.
      hav = sin((lat*pi/180.0_r8 - patchlatr(np))/2.0_r8)**2 + &
            cos(lat*pi/180.0_r8)*cos(patchlatr(np))* &
            sin((lon*pi/180.0_r8 - patchlonr(np))/2.0_r8)**2
      dist = 2.0_r8*6378.1_r8*asin(sqrt(max(0.0_r8,min(1.0_r8,hav))))
      r_p = r_p**2 / max(exp(-(dist**2)/(2.0_r8*this%locr**2)), 1.0e-6_r8)

   END SUBROUTINE ST_type_prepare

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Obs_ST
#endif
