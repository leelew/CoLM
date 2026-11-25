#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Data_SNOW
!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Observation handling for snow-depth data assimilation.
!
! WORKFLOW:
!   1. Init: Group stream configurations by source and initialize the source
!      module. Store observation search and localization parameters.
!   2. Forward: Prepare MODIS observations and H(x), then concatenate the
!      observations, error standard deviations, and matching ensemble H(x).
!   3. Prepare: Select nearby valid observations for each target patch and
!      convert error standard deviations to localized variances for LETKF.
!   4. Save: Copy MR model-space H(x) into the shared diagnostic arrays.
!   5. Clear: Release current observations and H(x), retaining configurations
!      and spatial mappings for subsequent forward calls.
!
! SUPPORTED:
!   MR / MODIS:
!      Fractional snow cover (FSC) observations.
!
! HISTORY:
!   Zhuo Wang, 11/2026: DA v1.0
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_Pixelset
   USE MOD_Vars_Global
   USE MOD_Vars_TimeInvariants
   USE MOD_DA_Data_MR
   USE MOD_DA_Vars_TimeVars
   IMPLICIT NONE

!-----------------------------------------------------------------------------
! Observation collection for snow-depth assimilation.
!-----------------------------------------------------------------------------
   type :: Snow_type

      type(Snow_config_type), allocatable :: cfgs(:)                        ! Configuration for each observation stream
      logical                             :: use_mr = .false.               ! Include MODIS observations

      real(r8)                            :: dres                           ! Observation search radius in latitude and longitude [deg]
      real(r8)                            :: locr                           ! Gaussian localization length scale [km]
      real(r8)                            :: infl                           ! Background error covariance inflation factor

      type(MR_set_type)                   :: mr                             ! MODIS observations and ensemble H(x)

      integer                             :: nobs = 0                       ! Total number of observations across enabled sources
      real(r8), allocatable               :: lat(:)                         ! Observation latitudes [deg]
      real(r8), allocatable               :: lon(:)                         ! Observation longitudes [deg]
      real(r8), allocatable               :: y(:)                           ! Observed fractional snow cover
      real(r8), allocatable               :: hx_obs(:,:)                        ! Observation-space ensemble H(x)
      real(r8), allocatable               :: sigma(:)                       ! Observation error standard deviations

   CONTAINS

      procedure                           :: init       => Snow_type_init    ! Configure enabled sources and assimilation parameters
      procedure                           :: forward    => Snow_type_forward ! Compute model- and observation-space H(x)
      procedure                           :: save       => Snow_type_save    ! Copy MR patch H(x) into forecast or analysis diagnostics
      procedure                           :: clear      => Snow_type_clear   ! Clear current observations and H(x)
      procedure                           :: concat     => Snow_type_concat  ! Combine observations and matching ensemble H(x)
      procedure                           :: prepare    => Snow_type_prepare ! Prepare the local observation package for LETKF

   END type Snow_type

   type :: Snow_config_type

      character(len=16)                   :: source_name
      character(len=32)                   :: sensor_name
      character(len=16)                   :: var_name
      integer                             :: istream

   END type Snow_config_type

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE Snow_type_init(this, pixelset, configs)

!-----------------------------------------------------------------------------
      class(Snow_type),      intent(inout):: this
      type(pixelset_type),  intent(in)    :: pixelset
      type(Snow_config_type), intent(in)  :: configs(:)

      integer                             :: i, imr, nmr
      type(MR_config_type), allocatable   :: mr_configs(:)
!-----------------------------------------------------------------------------
      CALL this%clear()
      this%use_mr = .false.
      IF (allocated(this%cfgs)) deallocate(this%cfgs)
      allocate(this%cfgs(size(configs)))
      this%cfgs = configs
      this%dres = DEF_DA_ASSIM_SNOW_DRES
      this%locr = DEF_DA_ASSIM_SNOW_LOCR
      this%infl = DEF_DA_ASSIM_SNOW_INFL

      ! Initialize MODIS observation streams.
      nmr = count(configs%source_name == 'MR' .or. configs%source_name == 'MODIS')
      IF (nmr > 0) THEN
         this%use_mr = .true.
         allocate(mr_configs(nmr))
         imr = 0
         DO i = 1, size(configs)
            IF (trim(configs(i)%source_name) == 'MR' .or. trim(configs(i)%source_name) == 'MODIS') THEN
               imr = imr + 1
               mr_configs(imr)%source_name = configs(i)%source_name
               mr_configs(imr)%sensor_name = configs(i)%sensor_name
               mr_configs(imr)%var_name    = configs(i)%var_name
            ENDIF
         ENDDO

         CALL this%mr%init(pixelset, mr_configs)

         deallocate(mr_configs)
      ENDIF

   END SUBROUTINE Snow_type_init

!-----------------------------------------------------------------------------

   SUBROUTINE Snow_type_clear(this)

!-----------------------------------------------------------------------------
      class(Snow_type), intent(inout)     :: this
!-----------------------------------------------------------------------------
      this%nobs = 0
      IF (allocated(this%lat   )) deallocate(this%lat   )
      IF (allocated(this%lon   )) deallocate(this%lon   )
      IF (allocated(this%y     )) deallocate(this%y     )
      IF (allocated(this%hx_obs)) deallocate(this%hx_obs)
      IF (allocated(this%sigma )) deallocate(this%sigma )

      CALL this%mr%clear()

   END SUBROUTINE Snow_type_clear

!-----------------------------------------------------------------------------

   SUBROUTINE Snow_type_concat(this)

!-----------------------------------------------------------------------------
      class(Snow_type), intent(inout)     :: this
!-----------------------------------------------------------------------------
      this%nobs = 0
      IF (allocated(this%lat   )) deallocate(this%lat   )
      IF (allocated(this%lon   )) deallocate(this%lon   )
      IF (allocated(this%y     )) deallocate(this%y     )
      IF (allocated(this%hx_obs)) deallocate(this%hx_obs)
      IF (allocated(this%sigma )) deallocate(this%sigma )

      IF (this%use_mr) this%nobs = this%mr%nobs
      IF (this%nobs == 0) RETURN

      allocate(this%lat(this%nobs))
      allocate(this%lon(this%nobs))
      allocate(this%y(this%nobs))
      allocate(this%hx_obs(this%nobs, DEF_DA_ENS_NUM))
      allocate(this%sigma(this%nobs))

      this%lat   = this%mr%lat
      this%lon   = this%mr%lon
      this%y     = this%mr%y
      this%hx_obs    = this%mr%hx_obs
      this%sigma = this%mr%sigma

   END SUBROUTINE Snow_type_concat

!-----------------------------------------------------------------------------

   SUBROUTINE Snow_type_forward(this, idate, deltim)

!-----------------------------------------------------------------------------
      class(Snow_type), intent(inout)     :: this
      integer,          intent(in)        :: idate(3)
      real(r8),         intent(in)        :: deltim
!-----------------------------------------------------------------------------
      IF (this%use_mr) CALL this%mr%forward(idate, deltim)
      CALL this%concat()

   END SUBROUTINE Snow_type_forward

!-----------------------------------------------------------------------------

   SUBROUTINE Snow_type_save(this, is_analysis)

!-----------------------------------------------------------------------------
      class(Snow_type), intent(in)        :: this
      logical,        intent(in)          :: is_analysis

      integer                             :: istream, is, imr, np, iens
      logical                             :: valid_hx
!-----------------------------------------------------------------------------
      IF (numpatch <= 0) RETURN
      IF (is_analysis .and. .not. DEF_DA) RETURN
      IF (.not. allocated(this%cfgs)) RETURN
      IF (.not. allocated(hx_ol)) RETURN
      IF (DEF_DA) THEN
         IF (.not. allocated(hx_f_ens) .or. .not. allocated(hx_a_ens)) RETURN
         IF (.not. allocated(hx_f) .or. .not. allocated(hx_a)) RETURN
      ENDIF

      DO istream = 1, size(this%cfgs)
         is = this%cfgs(istream)%istream
         IF (is < 1 .or. is > size(hx_ol,1)) CYCLE
         IF (is_analysis) THEN
            IF (.not. DEF_DA) CYCLE
            hx_a_ens(is,:,:) = spval
            hx_a    (is,:)   = spval
         ELSE
            hx_ol   (is,:)   = spval
            IF (DEF_DA) THEN
               hx_f_ens(is,:,:) = spval
               hx_f    (is,:)   = spval
            ENDIF
         ENDIF
      ENDDO

      imr = 0
      DO istream = 1, size(this%cfgs)
         IF (trim(this%cfgs(istream)%source_name) == 'MR' .or. &
            trim(this%cfgs(istream)%source_name) == 'MODIS') THEN
            imr = imr + 1
            is = this%cfgs(istream)%istream
            IF (is < 1 .or. is > size(hx_ol,1)) CYCLE
            IF (this%mr%sensors(imr)%y%nobs == 0) CYCLE
            IF (.not. allocated(this%mr%sensors(imr)%hx_model)) CYCLE

            IF (is_analysis) THEN
               DO np = 1, numpatch
                  valid_hx = all(this%mr%sensors(imr)%hx_model(1:DEF_DA_ENS_NUM,np) /= spval)
                  DO iens = 1, DEF_DA_ENS_NUM
                     hx_a_ens(is,iens,np) = this%mr%sensors(imr)%hx_model(iens,np)
                  ENDDO
                  IF (valid_hx) hx_a(is,np) = &
                     sum(this%mr%sensors(imr)%hx_model(1:DEF_DA_ENS_NUM,np)) / DEF_DA_ENS_NUM
               ENDDO
            ELSE
               DO np = 1, numpatch
                  hx_ol(is,np) = this%mr%sensors(imr)%hx_model(0,np)
                  IF (.not. DEF_DA) CYCLE
                  valid_hx = all(this%mr%sensors(imr)%hx_model(1:DEF_DA_ENS_NUM,np) /= spval)
                  DO iens = 1, DEF_DA_ENS_NUM
                     hx_f_ens(is,iens,np) = this%mr%sensors(imr)%hx_model(iens,np)
                  ENDDO
                  IF (valid_hx) hx_f(is,np) = &
                     sum(this%mr%sensors(imr)%hx_model(1:DEF_DA_ENS_NUM,np)) / DEF_DA_ENS_NUM
               ENDDO
            ENDIF

         ENDIF
      ENDDO

   END SUBROUTINE Snow_type_save

!-----------------------------------------------------------------------------

   SUBROUTINE Snow_type_prepare(this, np, nobs_p, y_p, hx_p, r_p)

!-----------------------------------------------------------------------------
      class(Snow_type),         intent(in):: this
      integer,                intent(in)  :: np
      integer,                intent(out) :: nobs_p
      real(r8), allocatable,  intent(out) :: y_p(:)
      real(r8), allocatable,  intent(out) :: hx_p(:,:)
      real(r8), allocatable,  intent(out) :: r_p(:)

      integer                             :: iens, iobs
      real(r8)                            :: lat_p, lon_p, dlat, dlon
      logical,  allocatable               :: valid(:)
      real(r8), allocatable               :: lat(:), lon(:), dist(:), hav(:)
!-----------------------------------------------------------------------------
      nobs_p = 0
      IF (this%nobs == 0) RETURN
      lat_p = patchlatr(np)*180.0_r8/pi
      lon_p = patchlonr(np)*180.0_r8/pi

      ! All source records are globally available before this local search.
      allocate(valid(this%nobs))
      valid = .false.
      DO iobs = 1, this%nobs
         dlat = abs(this%lat(iobs) - lat_p)
         dlon = modulo(abs(this%lon(iobs) - lon_p), 360.0_r8)
         IF (dlon > 180.0_r8) dlon = 360.0_r8 - dlon
         valid(iobs) = dlat <= this%dres .and. dlon <= this%dres .and. &
            this%y(iobs) >= 0.0_r8 .and. this%y(iobs) <= 1.0_r8 .and. &
            this%sigma(iobs) /= spval .and. this%sigma(iobs) > 0.0_r8 .and. &
            all(this%hx_obs(iobs,:) /= spval)
      ENDDO
      nobs_p = count(valid)
      IF (nobs_p == 0) RETURN

      allocate(y_p(nobs_p), hx_p(nobs_p,DEF_DA_ENS_NUM), r_p(nobs_p))
      allocate(lat(nobs_p), lon(nobs_p), dist(nobs_p), hav(nobs_p))
      y_p = pack(this%y, valid)
      DO iens = 1, DEF_DA_ENS_NUM
         hx_p(:,iens) = pack(this%hx_obs(:,iens), valid)
      ENDDO
      r_p = pack(this%sigma, valid)
      lat = pack(this%lat, valid)
      lon = pack(this%lon, valid)

      ! Convert the error standard deviation to variance, then inflate R.
      hav = sin((lat*pi/180.0_r8 - patchlatr(np))/2.0_r8)**2 + &
         cos(lat*pi/180.0_r8)*cos(patchlatr(np))* &
         sin((lon*pi/180.0_r8 - patchlonr(np))/2.0_r8)**2
      dist = 2.0_r8*6378.1_r8*asin(sqrt(max(0.0_r8,min(1.0_r8,hav))))
      r_p = r_p**2 / max(exp(-(dist**2)/(2.0_r8*this%locr**2)), 1.0e-6_r8)

   END SUBROUTINE Snow_type_prepare

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Data_SNOW
#endif
