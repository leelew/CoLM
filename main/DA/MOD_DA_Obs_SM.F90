#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Obs_SM
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Observation aggregation module for soil-moisture data assimilation.
!
!    This module provides SM_type, a soil-moisture-specific observation
!    container used by MOD_DA_Assim_SM. It sits above source-specific
!    observation modules such as PM, GNOS, AM, and in-situ observations.
!
!    SM_type owns the soil-moisture observation workflow:
!      1. receive a unified list of observation configs from the assimilation
!         driver;
!      2. split configs by source and initialize the corresponding source
!         modules, e.g., PM_set_type;
!      3. call each enabled source to read observations and calculate ensemble
!         H(x);
!      4. concatenate all enabled source/sensor observations into one vector;
!      5. prepare localized patch-level y, R, and H(x) arrays for LETKF.
!    ANA2MTQ uses patch-mean SYNOP analysis products held in TimeVariables.
!    Its nearby-patch rows are appended in prepare(), not concat(), so
!    nobs below counts only observations selected by the spatial search.
!
!    Source-level parameters, such as PM sensor, variable, and frequency,
!    are passed to the corresponding source module config type.
!    Soil-moisture assimilation parameters, such as observation search radius
!    and localization radius, are kept in SM_config_type.
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

   type :: SM_config_type

      character(len=16) :: source_name ! 'PM', 'GNOS', 'AM', 'situ', 'ANA2MTQ'
      character(len=32) :: sensor_name ! 'SMAP', 'SMOS', 'CYGNSS', 'FY3D', 'ISMN', etc.
      character(len=16) :: var_name    ! 'TB_H', 'TB_V', 'SM', 'DDM', '2MTQ'
      real(r8) :: fghz                 ! frequency in GHz
      integer  :: istream              ! index of this stream in the shared hx_* diagnostics

   END type SM_config_type


   type :: SM_type

      type(SM_config_type), allocatable :: cfgs(:)

      logical :: use_pm = .false.
      logical :: use_gnos = .false.
      logical :: use_ana2mtq = .false.

      ! assimilation parameters
      real(r8) :: dres      ! search radius for observations around a target patch [deg]
      real(r8) :: locr      ! localization radius [km]
      real(r8) :: infl      ! inflation factor for error covariance

      ! module-level data
      type(PM_set_type) :: pm
      !type(GNOS_set_type) :: gnos

      ! outputs used for assimilation
      integer :: nobs = 0
      real(r8), allocatable :: lat(:)
      real(r8), allocatable :: lon(:)
      real(r8), allocatable :: y(:)
      real(r8), allocatable :: hx(:,:)
      real(r8), allocatable :: r(:)

   CONTAINS

      procedure :: init       => SM_type_init
      procedure :: calcg      => SM_type_calc_on_grid
      procedure :: calcp      => SM_type_calc_on_pset
      procedure :: save       => SM_type_save
      procedure :: clear      => SM_type_clear
      procedure :: concat     => SM_type_concat
      procedure :: prepare    => SM_type_prepare_on_pset

   END type SM_type


!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE SM_type_init(this, pixelset, configs)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(SM_type),       intent(inout) :: this
      type(pixelset_type),  intent(in)    :: pixelset
      type(SM_config_type), intent(in)    :: configs(:)

!------------------------ Local Variables ------------------------------------
      integer :: i, ipm, npm
      type(PM_config_type), allocatable  :: pm_configs(:)

!-----------------------------------------------------------------------------
      this%use_pm = .false.
      this%use_gnos = .false.
      this%use_ana2mtq = .false.

      IF (allocated(this%cfgs)) deallocate(this%cfgs)

      allocate(this%cfgs(size(configs)))
      this%cfgs = configs
      this%dres = DEF_DA_ASSIM_SM_DRES
      this%locr = DEF_DA_ASSIM_SM_LOCR
      this%infl = DEF_DA_ASSIM_SM_INFL

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

      this%use_ana2mtq = any(configs%source_name == 'ANA2MTQ')
      IF (this%use_ana2mtq) THEN
         IF (.not. (any(DEF_DA_OBS_SOURCE == 'SYNOP' .and. DEF_DA_OBS_TARGET == '2MT') .and. &
                    any(DEF_DA_OBS_SOURCE == 'SYNOP' .and. DEF_DA_OBS_TARGET == '2MQ'))) &
            ERROR STOP 'ANA2MTQ requires SYNOP 2MT and 2MQ screen level analyses'
      ENDIF

   END SUBROUTINE SM_type_init

!-----------------------------------------------------------------------------

   SUBROUTINE SM_type_clear(this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(SM_type), intent(inout) :: this

!-----------------------------------------------------------------------------
      this%nobs = 0
      IF (allocated(this%lat))       deallocate(this%lat)
      IF (allocated(this%lon))       deallocate(this%lon)
      IF (allocated(this%y))         deallocate(this%y)
      IF (allocated(this%hx))        deallocate(this%hx)
      IF (allocated(this%r))         deallocate(this%r)

      CALL this%pm%clear()

   END SUBROUTINE SM_type_clear

!-----------------------------------------------------------------------------

   SUBROUTINE SM_type_concat(this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(SM_type), intent(inout) :: this

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

   END SUBROUTINE SM_type_concat

!-----------------------------------------------------------------------------

   SUBROUTINE SM_type_calc_on_grid(this, idate, deltim)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(SM_type), intent(inout) :: this
      integer,        intent(in)    :: idate(3)
      real(r8),       intent(in)    :: deltim

!-----------------------------------------------------------------------------
      CALL this%clear()

      IF (this%use_pm) CALL this%pm%calcg(idate, deltim)

      CALL this%concat()

   END SUBROUTINE SM_type_calc_on_grid

!-----------------------------------------------------------------------------

   SUBROUTINE SM_type_calc_on_pset(this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(SM_type), intent(inout) :: this

!-----------------------------------------------------------------------------
      IF (this%use_pm) CALL this%pm%calcp()

   END SUBROUTINE SM_type_calc_on_pset

!-----------------------------------------------------------------------------

   SUBROUTINE SM_type_save (this, is_analysis)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(SM_type), intent(in) :: this
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

   END SUBROUTINE SM_type_save

!-----------------------------------------------------------------------------

   SUBROUTINE SM_type_prepare_on_pset(this, np, nobs_p, y_p, hx_p, r_p)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(SM_type),         intent(in)  :: this
      integer,                intent(in)  :: np
      integer,                intent(out) :: nobs_p
      real(r8), allocatable,  intent(out) :: y_p(:)
      real(r8), allocatable,  intent(out) :: hx_p(:,:)
      real(r8), allocatable,  intent(out) :: r_p(:)

!------------------------ Local Variables ------------------------------------
      integer :: iens, iobs, ip, n, k
      real(r8) :: lat_p, lon_p, dlat, dlon
      real(r8) :: delta_t, delta_q
      logical,  allocatable :: valid(:)
      logical,  allocatable :: valid_ana2mtq(:)
      real(r8), allocatable :: lat(:)
      real(r8), allocatable :: lon(:)
      real(r8), allocatable :: dist(:)

!-----------------------------------------------------------------------------
      nobs_p = 0
      lat_p = patchlatr(np)*180/pi
      lon_p = patchlonr(np)*180/pi

      ! Select the existing PM/source observations around the target patch.
      allocate(valid(this%nobs))
      valid = .false.
      DO iobs = 1, this%nobs
         dlat = abs(this%lat(iobs) - lat_p)
         dlon = abs(this%lon(iobs) - lon_p)
         IF (dlon > 180) dlon = 360 - dlon

         valid(iobs) = dlat < this%dres .and. dlon < this%dres .and. &
                       this%y(iobs) /= spval .and. all(this%hx(iobs,:) /= spval)
      ENDDO
      n = count(valid)
      nobs_p = n

      ! ANA2MTQ observations are located at patch centres. Select T/q as a
      ! pair; every selected nearby patch contributes two observation rows.
      allocate(valid_ana2mtq(size(patchlatr)))
      valid_ana2mtq = .false.
      IF (this%use_ana2mtq) THEN
         DO ip = 1, size(patchlatr)
            dlat = abs(patchlatr(ip)*180/pi - lat_p)
            dlon = abs(patchlonr(ip)*180/pi - lon_p)
            IF (dlon > 180) dlon = 360 - dlon

            valid_ana2mtq(ip) = dlat < this%dres .and. dlon < this%dres .and. &
                                tref_f(ip) /= tref_a(ip) .and. qref_f(ip) /= qref_a(ip) .and. &
                                tref_a(ip) /= spval .and. qref_a(ip) /= spval .and. &
                                tref_f(ip) /= spval .and. qref_f(ip) /= spval .and. &
                                all(tref_ens(1:DEF_DA_ENS_NUM,ip) /= spval) .and. &
                                all(qref_ens(1:DEF_DA_ENS_NUM,ip) /= spval)

         ENDDO

         nobs_p = nobs_p + 2*count(valid_ana2mtq)
      ENDIF
      IF (nobs_p == 0) RETURN

      allocate(y_p(nobs_p))
      allocate(hx_p(nobs_p, DEF_DA_ENS_NUM))
      allocate(r_p(nobs_p))
      allocate(lat(nobs_p))
      allocate(lon(nobs_p))
      allocate(dist(nobs_p))

      ! Keep the original source rows first.
      IF (n > 0) THEN
         y_p(1:n) = pack(this%y, valid)
         DO iens = 1, DEF_DA_ENS_NUM
            hx_p(1:n,iens) = pack(this%hx(:,iens), valid)
         ENDDO
         r_p(1:n) = pack(this%r, valid)
         lat(1:n) = pack(this%lat, valid)
         lon(1:n) = pack(this%lon, valid)
      ENDIF

      ! Append paired T/q increments and ensemble values for nearby patches.
      IF (this%use_ana2mtq) THEN
         k = n
         DO ip = 1, size(valid_ana2mtq)
            IF (.not. valid_ana2mtq(ip)) CYCLE
            delta_t = tref_a(ip) - tref_f(ip)
            delta_q = qref_a(ip) - qref_f(ip)
            y_p(k+1) = delta_t
            y_p(k+2) = delta_q
            hx_p(k+1,:) = tref_ens(1:DEF_DA_ENS_NUM,ip) - delta_t
            hx_p(k+2,:) = qref_ens(1:DEF_DA_ENS_NUM,ip) - delta_q
            r_p(k+1:k+2) = 0.01_r8
            lat(k+1:k+2) = patchlatr(ip)*180/pi
            lon(k+1:k+2) = patchlonr(ip)*180/pi
            k = k + 2
         ENDDO
      ENDIF

      ! Apply the same distance-dependent R inflation to all observations.
      dist = 2.0*6378.1*asin(sqrt( &
         sin((lat*pi/180.0 - patchlatr(np))/2.0)**2 + &
         cos(lat*pi/180.0)*cos(patchlatr(np))* &
         sin((lon*pi/180.0 - patchlonr(np))/2.0)**2))
      r_p = r_p / max(exp(-(dist**2)/(2.0*this%locr**2)), 1.0e-6)

   END SUBROUTINE SM_type_prepare_on_pset

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Obs_SM
#endif
