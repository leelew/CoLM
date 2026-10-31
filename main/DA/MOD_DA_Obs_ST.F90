#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Obs_ST
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Observation aggregation module for soil-temperature data assimilation.
!
!    This module provides ST_type, a soil-temperature-specific observation
!    container used by MOD_DA_Assim_ST. It mirrors MOD_DA_Obs_SM and sits
!    above source-specific observation modules such as PM, MR, and in-situ
!    observations.
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
! HISTORY:
!    Lu Li, 05/2026: First implementation
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_DA_Obs_PM
   USE MOD_DA_Vars_TimeVariables
   USE MOD_Namelist
   USE MOD_Pixelset
   USE MOD_Vars_Global
   USE MOD_Vars_TimeInvariants
   IMPLICIT NONE

   type :: ST_config_type

      character(len=16) :: source_name ! 'PM', 'MR', 'situ'
      character(len=32) :: sensor_name ! 'SMAP', 'MODIS', 'FY3D', 'ISMN', etc.
      character(len=16) :: var_name    ! 'TB_H', 'TB_V', 'LST'
      real(r8) :: fghz                 ! frequency in GHz
      integer  :: istream              ! index of this stream in the shared hx_* diagnostics

   END type ST_config_type


   type :: ST_type

      type(ST_config_type), allocatable :: cfgs(:)

      logical :: use_pm = .false.

      ! assimilation parameters
      real(r8) :: dres      ! search radius for observations around a target patch [deg]
      real(r8) :: locr      ! localization radius [km]
      real(r8) :: infl      ! inflation factor for error covariance

      ! module-level data
      type(PM_set_type) :: pm

      ! outputs used for assimilation
      integer :: nobs = 0
      real(r8), allocatable :: lat(:)
      real(r8), allocatable :: lon(:)
      real(r8), allocatable :: y(:)
      real(r8), allocatable :: hx(:,:)
      real(r8), allocatable :: r(:)

   CONTAINS

      procedure :: init       => ST_type_init
      procedure :: calcg      => ST_type_calc_on_grid
      procedure :: calcp      => ST_type_calc_on_pset
      procedure :: save       => ST_type_save
      procedure :: clear      => ST_type_clear
      procedure :: concat     => ST_type_concat
      procedure :: prepare    => ST_type_prepare_on_pset

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

      !//TODO: initialize other source types (MR, situ)

   END SUBROUTINE ST_type_init

!-----------------------------------------------------------------------------

   SUBROUTINE ST_type_clear(this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(ST_type), intent(inout) :: this

!-----------------------------------------------------------------------------
      this%nobs = 0
      IF (allocated(this%lat))       deallocate(this%lat)
      IF (allocated(this%lon))       deallocate(this%lon)
      IF (allocated(this%y))         deallocate(this%y)
      IF (allocated(this%hx))        deallocate(this%hx)
      IF (allocated(this%r))         deallocate(this%r)

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
      IF (allocated(this%lat)) deallocate(this%lat)
      IF (allocated(this%lon)) deallocate(this%lon)
      IF (allocated(this%y))   deallocate(this%y)
      IF (allocated(this%hx))  deallocate(this%hx)
      IF (allocated(this%r))   deallocate(this%r)

      IF (this%use_pm) this%nobs = this%nobs + this%pm%nobs
      IF (this%nobs == 0) RETURN

      allocate(this%lat(this%nobs))
      allocate(this%lon(this%nobs))
      allocate(this%y(this%nobs))
      allocate(this%hx(this%nobs, DEF_DA_ENS_NUM))
      allocate(this%r(this%nobs))

      offset = 0

      IF (this%use_pm) THEN
         n = this%pm%nobs
         IF (n > 0) THEN
            this%lat(offset+1:offset+n) = this%pm%lat
            this%lon(offset+1:offset+n) = this%pm%lon
            this%y(offset+1:offset+n) = this%pm%y
            this%hx(offset+1:offset+n,:) = this%pm%hx
            this%r(offset+1:offset+n) = this%pm%r
            offset = offset + n
         ENDIF
      ENDIF

   END SUBROUTINE ST_type_concat

!-----------------------------------------------------------------------------

   SUBROUTINE ST_type_calc_on_grid(this, idate, deltim)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(ST_type), intent(inout) :: this
      integer,        intent(in)    :: idate(3)
      real(r8),       intent(in)    :: deltim

!-----------------------------------------------------------------------------
      CALL this%clear()

      IF (this%use_pm) CALL this%pm%calcg(idate, deltim)

      CALL this%concat()

   END SUBROUTINE ST_type_calc_on_grid

!-----------------------------------------------------------------------------

   SUBROUTINE ST_type_calc_on_pset(this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(ST_type), intent(inout) :: this

!-----------------------------------------------------------------------------
      IF (this%use_pm) CALL this%pm%calcp()

   END SUBROUTINE ST_type_calc_on_pset

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
                  valid_hx = all(this%pm%sensors(ipm)%hx_pset(1:DEF_DA_ENS_NUM,np) /= spval)
                  DO iens = 1, DEF_DA_ENS_NUM
                     hx_a_ens(is,iens,np) = this%pm%sensors(ipm)%hx_pset(iens,np)
                  ENDDO
                  IF (valid_hx) hx_a(is,np) = &
                     sum(this%pm%sensors(ipm)%hx_pset(1:DEF_DA_ENS_NUM,np)) / DEF_DA_ENS_NUM
               ENDDO
            ELSE
               DO np = 1, numpatch
                  hx_ol(is,np) = this%pm%sensors(ipm)%hx_pset(0,np)
                  valid_hx = all(this%pm%sensors(ipm)%hx_pset(1:DEF_DA_ENS_NUM,np) /= spval)
                  DO iens = 1, DEF_DA_ENS_NUM
                     hx_f_ens(is,iens,np) = this%pm%sensors(ipm)%hx_pset(iens,np)
                  ENDDO
                  IF (valid_hx) hx_f(is,np) = &
                     sum(this%pm%sensors(ipm)%hx_pset(1:DEF_DA_ENS_NUM,np)) / DEF_DA_ENS_NUM
               ENDDO
            ENDIF
         ENDIF
      ENDDO

   END SUBROUTINE ST_type_save

!-----------------------------------------------------------------------------

   SUBROUTINE ST_type_prepare_on_pset(this, np, nobs_p, y_p, hx_p, r_p)

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
      integer :: iens, iobs, n
      real(r8) :: lat_p, lon_p, dlat, dlon
      logical,  allocatable :: valid(:)
      real(r8), allocatable :: lat(:)
      real(r8), allocatable :: lon(:)
      real(r8), allocatable :: dist(:)

!-----------------------------------------------------------------------------
      nobs_p = 0

      ! find the nearest y values around target patch
      allocate(valid(this%nobs))
      valid = .false.
      IF (this%nobs > 0) THEN
         lat_p = patchlatr(np)*180/pi
         lon_p = patchlonr(np)*180/pi
         DO iobs = 1, this%nobs
            dlat = abs(this%lat(iobs) - lat_p)
            dlon = abs(this%lon(iobs) - lon_p)
            IF (dlon > 180) dlon = 360 - dlon

            valid(iobs) = dlat < this%dres .and. dlon < this%dres .and. &
                          this%y(iobs) /= spval .and. all(this%hx(iobs,:) /= spval)
         ENDDO
      ENDIF
      n = count(valid)
      IF (n == 0) RETURN
      nobs_p = n

      ! pack the valid y and H(x)
      allocate(y_p(n))
      allocate(hx_p(n, DEF_DA_ENS_NUM))
      y_p = pack(this%y, valid)
      DO iens = 1, DEF_DA_ENS_NUM
         hx_p(:,iens) = pack(this%hx(:,iens), valid)
      ENDDO

      ! Convert the error standard deviation to variance and localize R.
      allocate(r_p(n))
      allocate(lat(n))
      allocate(lon(n))
      allocate(dist(n))
      lat = pack(this%lat, valid)
      lon = pack(this%lon, valid)
      dist = 2.0*6378.1*asin(sqrt( &
         sin((lat*pi/180.0 - patchlatr(np))/2.0)**2 + &
         cos(lat*pi/180.0)*cos(patchlatr(np))* &
         sin((lon*pi/180.0 - patchlonr(np))/2.0)**2))
      r_p = pack(this%r, valid)**2 / max(exp(-(dist**2)/(2.0*this%locr**2)), 1.0e-6)

   END SUBROUTINE ST_type_prepare_on_pset

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Obs_ST
#endif
