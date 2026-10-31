#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Obs_SM
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Observation aggregation module for soil-moisture data assimilation.
!
!    This module provides SM_type, a soil-moisture-specific observation
!    container used by MOD_DA_Assim_SM. It sits above source-specific
!    observation modules such as PM, GNOS, AM, ANA2MTQ.
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
!
!    ANA2MTQ selects patches with nonzero T and q analysis increments,
!    gathers their paired data on master, and broadcasts to all workers.
!    Its nobs counts patches; the SM arrays contain two rows per selected patch.
!
!    Source-level parameters, such as PM sensor, variable, and frequency,
!    are passed to the corresponding source module config type.
!    Soil-moisture assimilation parameters, such as observation search radius
!    and localization radius, are kept in SM_type.
!
! HISTORY:
!    Lu Li, 05/2026: First implementation
!    Zhilong Fan, Lu Li, 10/2026: First operational version
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_DA_Obs_PM
   USE MOD_DA_Obs_ANA2MTQ, ONLY: ANA2MTQ_type
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

      ! source-specific observation containers
      type(PM_set_type) :: pm
      type(ANA2MTQ_type) :: ana2mtq
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
      CALL this%clear()
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
         CALL this%ana2mtq%init()
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
      CALL this%ana2mtq%clear()

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
      IF (this%use_ana2mtq) this%nobs = this%nobs + 2*this%ana2mtq%nobs
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

      ! Each selected ANA2MTQ patch supplies one T row and one q row.
      IF (this%use_ana2mtq) THEN
         n = this%ana2mtq%nobs
         IF (n > 0) THEN
            this%lat(offset+1:offset+2*n:2) = this%ana2mtq%lat
            this%lat(offset+2:offset+2*n:2) = this%ana2mtq%lat
            this%lon(offset+1:offset+2*n:2) = this%ana2mtq%lon
            this%lon(offset+2:offset+2*n:2) = this%ana2mtq%lon
            this%y(offset+1:offset+2*n:2) = this%ana2mtq%t%y
            this%y(offset+2:offset+2*n:2) = this%ana2mtq%q%y
            this%hx(offset+1:offset+2*n:2,:) = this%ana2mtq%t%hx
            this%hx(offset+2:offset+2*n:2,:) = this%ana2mtq%q%hx
            this%r(offset+1:offset+2*n:2) = this%ana2mtq%t%r
            this%r(offset+2:offset+2*n:2) = this%ana2mtq%q%r
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
      IF (this%use_ana2mtq) CALL this%ana2mtq%calcp()

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
                       this%y(iobs) /= spval .and. all(this%hx(iobs,:) /= spval)
      ENDDO
      nobs_p = count(valid)
      IF (nobs_p == 0) RETURN

      allocate(y_p(nobs_p), hx_p(nobs_p,DEF_DA_ENS_NUM), r_p(nobs_p))
      allocate(lat(nobs_p), lon(nobs_p), dist(nobs_p), hav(nobs_p))
      y_p = pack(this%y, valid)
      DO iens = 1, DEF_DA_ENS_NUM
         hx_p(:,iens) = pack(this%hx(:,iens), valid)
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

   END SUBROUTINE SM_type_prepare_on_pset

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Obs_SM
#endif
