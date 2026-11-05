#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Obs_TQ
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Observation collection for 2 m temperature and humidity assimilation.
!
!    This module manages SYNOP observations for one 2 m analysis target.
!    SYNOP stores the combined observations and H(x); this module selects
!    nearby records and localizes their observation error variances without
!    copying the source arrays.
!
!    The assimilation driver uses separate collections for temperature and
!    humidity and controls their analysis order. Source modules compute H(x);
!    this collection saves their model-space H(x) in the shared diagnostics.
!
! HISTORY:
!    Lu Li, 11/2026: DA v1.0
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_LandPatch
   USE MOD_Namelist
   USE MOD_Pixelset
   USE MOD_SPMD_Task
   USE MOD_Vars_Global
   USE MOD_Vars_TimeInvariants
   USE MOD_DA_Obs_SYNOP
   USE MOD_DA_Vars_TimeVars
   IMPLICIT NONE

!-----------------------------------------------------------------------------
! Configuration for one stream used in 2 m temperature or humidity assimilation
!-----------------------------------------------------------------------------
   type :: TQ_config_type

      character(len=16) :: source_name                         ! Observation source class, e.g., 'site'
      character(len=32) :: sensor_name                         ! Station network or dataset identifier
      character(len=16) :: var_name                            ! Observation variable: '2MT' or '2MQ'
      integer :: istream                                       ! Stream index in the shared H(x) diagnostics

   END type TQ_config_type

!-----------------------------------------------------------------------------
! Observation collection for one 2 m analysis target
!-----------------------------------------------------------------------------
   type :: TQ_type

      type(TQ_config_type), allocatable :: cfgs(:)             ! Observation stream configurations
      logical :: use_synop = .false.                           ! Include SYNOP station observations

      ! Observation selection and LETKF parameters.
      real(r8) :: dres                                         ! Observation search half-width [degrees]
      real(r8) :: locr                                         ! Gaussian localization length scale [km]
      real(r8) :: infl                                         ! Ensemble covariance inflation factor [-]

      ! Source-specific observation containers.
      type(SYNOP_set_type) :: synop                            ! Station observations and ensemble H(x)

   CONTAINS

      procedure :: init       => TQ_type_init                  ! Configure sources and assimilation parameters
      procedure :: calc_obs   => TQ_type_calc_obs              ! Prepare observation-space y and H(x) from enabled sources
      procedure :: calc_model => TQ_type_calc_model            ! Recompute source model-space H(x) after state updates
      procedure :: save       => TQ_type_save                  ! Copy model-space H(x) into shared diagnostics
      procedure :: clear      => TQ_type_clear                 ! Clear time-step data while retaining configurations
      procedure :: prepare    => TQ_type_prepare               ! Select nearby observations and localize error variances

   END type TQ_type

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE TQ_type_init (this, pixelset, configs)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(TQ_type), intent(inout)    :: this
      type(pixelset_type), intent(in)   :: pixelset
      type(TQ_config_type), intent(in) :: configs(:)

!------------------------ Local Variables ------------------------------------
      integer :: i
      integer :: isynop
      integer :: nsynop
      type(SYNOP_config_type), allocatable :: synop_configs(:)

!-----------------------------------------------------------------------------
      CALL this%clear()
      this%use_synop = .false.
      IF (allocated(this%cfgs)) deallocate(this%cfgs)
      allocate(this%cfgs(size(configs)))
      this%cfgs = configs
      this%dres = DEF_DA_ASSIM_TQ_DRES
      this%locr = DEF_DA_ASSIM_TQ_LOCR
      this%infl = DEF_DA_ASSIM_TQ_INFL

      ! Configure the station source for site streams.
      nsynop = count(configs%source_name == 'site')
      allocate(synop_configs(nsynop))
      isynop = 0
      DO i = 1, size(configs)
         IF (trim(configs(i)%source_name) /= 'site') CYCLE
         isynop = isynop + 1
         synop_configs(isynop)%source_name = configs(i)%source_name
         synop_configs(isynop)%sensor_name = configs(i)%sensor_name
         synop_configs(isynop)%var_name = configs(i)%var_name
      ENDDO
      CALL this%synop%init(pixelset, synop_configs)
      this%use_synop = nsynop > 0

   END SUBROUTINE TQ_type_init

!-----------------------------------------------------------------------------

   SUBROUTINE TQ_type_calc_obs (this, idate, deltim)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(TQ_type), intent(inout) :: this
      integer, intent(in)            :: idate(3)
      real(r8), intent(in)           :: deltim

!-----------------------------------------------------------------------------
      CALL this%clear()
      IF (this%use_synop) CALL this%synop%calc_obs(idate, deltim)

   END SUBROUTINE TQ_type_calc_obs

!-----------------------------------------------------------------------------

   SUBROUTINE TQ_type_calc_model (this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(TQ_type), intent(inout) :: this

!-----------------------------------------------------------------------------
      IF (this%use_synop) CALL this%synop%calc_model()

   END SUBROUTINE TQ_type_calc_model

!-----------------------------------------------------------------------------

   SUBROUTINE TQ_type_clear (this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(TQ_type), intent(inout) :: this

!-----------------------------------------------------------------------------
      CALL this%synop%clear()

   END SUBROUTINE TQ_type_clear

!-----------------------------------------------------------------------------

   SUBROUTINE TQ_type_save (this, is_analysis)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(TQ_type), intent(in) :: this
      logical, intent(in)        :: is_analysis

!------------------------ Local Variables ------------------------------------
      integer :: i
      integer :: is
      integer :: isynop                                        
      integer :: np
      logical :: valid_hx

!-----------------------------------------------------------------------------
      IF (.not. p_is_worker) RETURN
      IF (nsource <= 0 .or. numpatch <= 0) RETURN
      IF (.not. allocated(hx_f_ens)) RETURN
      IF (.not. this%use_synop) RETURN
      IF (.not. allocated(this%synop%sensors)) RETURN

      isynop = 0
      DO i = 1, size(this%cfgs)
         IF (trim(this%cfgs(i)%source_name) /= 'site') CYCLE
         isynop = isynop + 1
         is = this%cfgs(i)%istream
         IF (is_analysis) THEN
            hx_a_ens(is,:,:) = spval
            hx_a(is,:) = spval
         ELSE
            hx_ol(is,:) = spval
            hx_f_ens(is,:,:) = spval
            hx_f(is,:) = spval
         ENDIF
         IF (.not. allocated(this%synop%sensors(isynop)%hx_model)) CYCLE

         DO np = 1, numpatch
            valid_hx = all(this%synop%sensors(isynop)%hx_model(1:DEF_DA_ENS_NUM,np) /= spval)
            IF (is_analysis) THEN
               hx_a_ens(is,1:DEF_DA_ENS_NUM,np) = this%synop%sensors(isynop)%hx_model(1:DEF_DA_ENS_NUM,np)
               IF (valid_hx) hx_a(is,np) = &
                  sum(this%synop%sensors(isynop)%hx_model(1:DEF_DA_ENS_NUM,np)) / DEF_DA_ENS_NUM
            ELSE
               hx_ol(is,np) = this%synop%sensors(isynop)%hx_model(0,np)
               hx_f_ens(is,1:DEF_DA_ENS_NUM,np) = this%synop%sensors(isynop)%hx_model(1:DEF_DA_ENS_NUM,np)
               IF (valid_hx) hx_f(is,np) = &
                  sum(this%synop%sensors(isynop)%hx_model(1:DEF_DA_ENS_NUM,np)) / DEF_DA_ENS_NUM
            ENDIF
         ENDDO
      ENDDO

   END SUBROUTINE TQ_type_save

!-----------------------------------------------------------------------------

   SUBROUTINE TQ_type_prepare (this, np, nobs_p, y_p, hx_p, r_p)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(TQ_type),         intent(in)  :: this
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
      allocate(valid(this%synop%nobs))
      valid = .false.
      DO iobs = 1, this%synop%nobs
         dlat = abs(this%synop%lat(iobs) - lat_p)
         dlon = modulo(abs(this%synop%lon(iobs) - lon_p), 360.0_r8)
         IF (dlon > 180.0_r8) dlon = 360.0_r8 - dlon
         valid(iobs) = dlat < this%dres .and. dlon < this%dres .and. &
                       this%synop%y(iobs) /= spval .and. all(this%synop%hx_obs(iobs,:) /= spval)
      ENDDO
      nobs_p = count(valid)
      IF (nobs_p == 0) RETURN

      allocate(y_p(nobs_p), hx_p(nobs_p,DEF_DA_ENS_NUM), r_p(nobs_p))
      allocate(lat(nobs_p), lon(nobs_p), dist(nobs_p), hav(nobs_p))
      y_p = pack(this%synop%y, valid)
      DO iens = 1, DEF_DA_ENS_NUM
         hx_p(:,iens) = pack(this%synop%hx_obs(:,iens), valid)
      ENDDO
      r_p = pack(this%synop%r, valid)
      lat = pack(this%synop%lat, valid)
      lon = pack(this%synop%lon, valid)

      ! Convert the error standard deviation to variance, then inflate R.
      hav = sin((lat*pi/180.0_r8 - patchlatr(np))/2.0_r8)**2 + &
            cos(lat*pi/180.0_r8)*cos(patchlatr(np))* &
            sin((lon*pi/180.0_r8 - patchlonr(np))/2.0_r8)**2
      dist = 2.0_r8*6378.1_r8*asin(sqrt(max(0.0_r8,min(1.0_r8,hav))))
      r_p = r_p**2 / max(exp(-(dist**2)/(2.0_r8*this%locr**2)), 1.0e-6_r8)

   END SUBROUTINE TQ_type_prepare

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Obs_TQ
#endif
