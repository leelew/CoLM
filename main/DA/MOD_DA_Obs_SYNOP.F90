#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Obs_SYNOP
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    SYNOP in-situ observation handling for 2 m temperature and humidity
!    data assimilation. Temperature and humidity are assimilated by two
!    independent streams, configured with
!       DEF_DA_OBS_SOURCE = 'situ'
!       DEF_DA_OBS_SENSOR = 'SYNOP'
!       DEF_DA_OBS_VAR    = '2MT'
!       DEF_DA_OBS_TARGET = '2MT'
!    and
!       DEF_DA_OBS_SOURCE = 'situ'
!       DEF_DA_OBS_SENSOR = 'SYNOP'
!       DEF_DA_OBS_VAR    = '2MQ'
!       DEF_DA_OBS_TARGET = '2MQ'
!
!    DEF_DA_OBS_VAR is the observed variable and the hourly file variable
!    name; DEF_DA_OBS_TARGET selects the analysis variable.
!
!    The hourly files are
!       <DEF_DA_OBS_DIR>/<sensor>/2MT_YYYY_MM_DD_HH.nc
!       <DEF_DA_OBS_DIR>/<sensor>/2MQ_YYYY_MM_DD_HH.nc
!    containing the station observations '2MT' (2 m air temperature [K])
!    and '2MQ' (2 m relative humidity [0-1]) respectively, together with
!    'latitude', 'longitude', 'time', and 'id'. The station table is
!       <DEF_DA_OBS_DIR>/<sensor>/sites.nc
!
!    Observation operators:
!       H_T   = tref_ens
!       H_RH  = RH(qref_ens, tref_ens, forc_psrf)
!    where RH is calculated from the model 2 m specific humidity and
!    temperature following the legacy SYNOP implementation. The analysis
!    variables are the model 2 m air temperature and specific humidity,
!    tref_ens and qref_ens. A '2MT' stream reads only the 2MT file and
!    updates tref_ens; a '2MQ' stream reads only the 2MQ file and updates
!    qref_ens. The 2MT analysis runs first, so H_RH for a '2MQ' stream is
!    evaluated with the analyzed tref_ens.
!
!    SYNOP_set_type owns the observation workflow:
!      1. read the 2 m temperature and/or relative humidity observations;
!      2. derive rhref_ens from qref_ens, tref_ens, and surface pressure;
!      3. gather tref_ens and/or rhref_ens to station locations through
!         site_obs_type;
!      4. concatenate the enabled variables into one observation vector;
!      5. prepare localized patch-level y, R, and H(x) arrays for LETKF.
!
!    Observation errors 'r' are standard deviations: K for 2 m temperature
!    and RH fraction for relative humidity.
!
! HISTORY:
!    Lu Li, 05/2026: Initial implementation
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_DA_Obs_SiteType
   USE MOD_DA_Vars_TimeVariables
   USE MOD_Namelist
   USE MOD_Pixelset
   USE MOD_Vars_Global
   USE MOD_Vars_TimeInvariants
   USE MOD_Vars_1DForcing, ONLY: forc_psrf
   USE MOD_SPMD_Task
   IMPLICIT NONE

   PUBLIC :: SYNOP_calc_rh

   type :: SYNOP_config_type

      character(len=16) :: source_name ! 'situ'
      character(len=32) :: sensor_name ! 'SYNOP'
      character(len=16) :: var_name    ! observed variable: '2MT' / '2MQ'
      integer  :: istream              ! namelist index of this stream; slot in the
                                       ! shared hx_* diagnostics
      logical  :: use_t  = .false.     ! assimilate the 2 m temperature variable
      logical  :: use_rh = .false.     ! assimilate the 2 m relative humidity variable

   END type SYNOP_config_type


   type :: SYNOP_sensor_type

      character(len=32) :: sensor_name = ''
      character(len=16) :: var_name    = ''

      ! station observations of 2 m temperature and relative humidity
      type(site_obs_type) :: t
      type(site_obs_type) :: rh

      ! ensemble H(x) gathered at the station locations [nobs, 0:ens]
      real(r8), allocatable :: hx_t(:,:)
      real(r8), allocatable :: hx_rh(:,:)

   END type SYNOP_sensor_type


   type :: SYNOP_set_type

      type(SYNOP_config_type), allocatable :: cfgs(:)
      type(SYNOP_sensor_type), allocatable :: sensors(:)

      ! assimilation parameters
      real(r8) :: dres      ! search radius for observations around a target patch [deg]
      real(r8) :: locr      ! localization radius [km]
      real(r8) :: infl      ! inflation factor for error covariance

      ! enabled variables
      logical  :: use_t  = .false.     ! any configured stream assimilates 2 m temperature
      logical  :: use_rh = .false.     ! any configured stream assimilates 2 m relative humidity

      ! outputs used for assimilation
      integer :: nobs = 0
      real(r8), allocatable :: lat(:)
      real(r8), allocatable :: lon(:)
      real(r8), allocatable :: y(:)
      real(r8), allocatable :: hx(:,:)
      real(r8), allocatable :: r(:)

   CONTAINS

      procedure :: init    => SYNOP_type_init
      procedure :: calcg   => SYNOP_type_calc_on_grid
      procedure :: save    => SYNOP_type_save
      procedure :: clear   => SYNOP_type_clear
      procedure :: concat  => SYNOP_type_concat
      procedure :: prepare => SYNOP_type_prepare_on_pset

   END type SYNOP_set_type


!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE SYNOP_calc_rh ()

!-----------------------------------------------------------------------------
! Derive 2 m relative humidity from the model specific humidity and
! temperature, using the legacy SYNOP formula.
!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Local Variables ------------------------------------
      integer :: iens, np
      real(r8) :: tc, es, e

!-----------------------------------------------------------------------------
      IF (.not. p_is_worker) RETURN
      IF (numpatch <= 0) RETURN
      IF (.not. allocated(tref_ens))  RETURN
      IF (.not. allocated(qref_ens))  RETURN
      IF (.not. allocated(rhref_ens)) RETURN
      IF (.not. allocated(forc_psrf)) RETURN

      DO np = 1, numpatch
         DO iens = 0, DEF_DA_ENS_NUM
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

   END SUBROUTINE SYNOP_calc_rh

!-----------------------------------------------------------------------------

   SUBROUTINE SYNOP_type_init (this, pixelset, configs)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(SYNOP_set_type),    intent(inout) :: this
      type(pixelset_type),      intent(in)    :: pixelset
      type(SYNOP_config_type),  intent(in)    :: configs(:)

!------------------------ Local Variables ------------------------------------
      integer :: i

!-----------------------------------------------------------------------------
      IF (allocated(this%cfgs))    deallocate(this%cfgs)
      IF (allocated(this%sensors)) deallocate(this%sensors)
      this%nobs = 0

      IF (size(configs) == 0) RETURN

      allocate(this%cfgs(size(configs)))
      this%cfgs = configs

      this%use_t  = any(configs%use_t)
      this%use_rh = any(configs%use_rh)

      this%dres = DEF_DA_ASSIM_2MTQ_DRES
      this%locr = DEF_DA_ASSIM_2MTQ_LOCR
      this%infl = DEF_DA_ASSIM_2MTQ_INFL

      ! One sensor owns the temperature and relative humidity pair. Each
      ! variable is stored in its own hourly file, named after the file
      ! variable, following the same convention as the gridded observation
      ! type. Only the variables enabled by the stream are initialized.
      ! The namelist VAR ('2MT' / '2MQ') is also the hourly file variable.
      allocate(this%sensors(size(configs)))
      DO i = 1, size(configs)
         this%sensors(i)%sensor_name = configs(i)%sensor_name
         this%sensors(i)%var_name    = configs(i)%var_name

         IF (configs(i)%use_t) THEN
            CALL this%sensors(i)%t%init(pixelset, configs(i)%source_name, &
               configs(i)%sensor_name, configs(i)%var_name)
         ENDIF
         IF (configs(i)%use_rh) THEN
            CALL this%sensors(i)%rh%init(pixelset, configs(i)%source_name, &
               configs(i)%sensor_name, configs(i)%var_name)
         ENDIF
      ENDDO

   END SUBROUTINE SYNOP_type_init

!-----------------------------------------------------------------------------

   SUBROUTINE SYNOP_type_clear (this)

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
      IF (allocated(this%hx))  deallocate(this%hx)
      IF (allocated(this%r))   deallocate(this%r)

      IF (allocated(this%sensors)) THEN
         DO i = 1, size(this%sensors)
            CALL this%sensors(i)%t%clear()
            CALL this%sensors(i)%rh%clear()
            IF (allocated(this%sensors(i)%hx_t))  deallocate(this%sensors(i)%hx_t)
            IF (allocated(this%sensors(i)%hx_rh)) deallocate(this%sensors(i)%hx_rh)
         ENDDO
      ENDIF

   END SUBROUTINE SYNOP_type_clear

!-----------------------------------------------------------------------------

   SUBROUTINE SYNOP_type_calc_on_grid (this, idate, deltim)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(SYNOP_set_type), intent(inout) :: this
      integer,               intent(in)    :: idate(3)
      real(r8),              intent(in)    :: deltim

!------------------------ Local Variables ------------------------------------
      integer :: i

!-----------------------------------------------------------------------------
      CALL this%clear()

      IF (.not. allocated(this%sensors)) RETURN

      ! H(x) for humidity is relative humidity, derived from qref_ens;
      ! only needed when a relative humidity stream is enabled
      IF (this%use_rh) CALL SYNOP_calc_rh()

      DO i = 1, size(this%sensors)
         IF (this%cfgs(i)%use_t) THEN
            CALL this%sensors(i)%t%read(idate, deltim)
            CALL this%sensors(i)%t%gather(tref_ens, this%sensors(i)%hx_t)
         ENDIF
         IF (this%cfgs(i)%use_rh) THEN
            CALL this%sensors(i)%rh%read(idate, deltim)
            CALL this%sensors(i)%rh%gather(rhref_ens, this%sensors(i)%hx_rh)
         ENDIF
      ENDDO

      CALL this%concat()

   END SUBROUTINE SYNOP_type_calc_on_grid

!-----------------------------------------------------------------------------

   SUBROUTINE SYNOP_type_concat (this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(SYNOP_set_type), intent(inout) :: this

!------------------------ Local Variables ------------------------------------
      integer :: i, n, nt, nrh
      integer :: offset

!-----------------------------------------------------------------------------
      this%nobs = 0
      IF (allocated(this%lat)) deallocate(this%lat)
      IF (allocated(this%lon)) deallocate(this%lon)
      IF (allocated(this%y))   deallocate(this%y)
      IF (allocated(this%hx))  deallocate(this%hx)
      IF (allocated(this%r))   deallocate(this%r)

      IF (.not. allocated(this%sensors)) RETURN

      DO i = 1, size(this%sensors)
         IF (this%cfgs(i)%use_t)  this%nobs = this%nobs + this%sensors(i)%t%nobs
         IF (this%cfgs(i)%use_rh) this%nobs = this%nobs + this%sensors(i)%rh%nobs
      ENDDO
      IF (this%nobs == 0) RETURN

      allocate(this%lat(this%nobs))
      allocate(this%lon(this%nobs))
      allocate(this%y(this%nobs))
      allocate(this%hx(this%nobs, DEF_DA_ENS_NUM))
      allocate(this%r(this%nobs))

      offset = 0

      DO i = 1, size(this%sensors)
         ! 2 m temperature observations
         IF (this%cfgs(i)%use_t) THEN
            nt = this%sensors(i)%t%nobs
            IF (nt > 0) THEN
               n = offset + 1
               this%lat(n:offset+nt) = this%sensors(i)%t%lat
               this%lon(n:offset+nt) = this%sensors(i)%t%lon
               this%y  (n:offset+nt) = this%sensors(i)%t%y
               this%hx (n:offset+nt,:) = this%sensors(i)%hx_t(:,1:DEF_DA_ENS_NUM)
               this%r  (n:offset+nt) = this%sensors(i)%t%r
               offset = offset + nt
            ENDIF
         ENDIF

         ! 2 m relative humidity observations
         IF (this%cfgs(i)%use_rh) THEN
            nrh = this%sensors(i)%rh%nobs
            IF (nrh > 0) THEN
               n = offset + 1
               this%lat(n:offset+nrh) = this%sensors(i)%rh%lat
               this%lon(n:offset+nrh) = this%sensors(i)%rh%lon
               this%y  (n:offset+nrh) = this%sensors(i)%rh%y
               this%hx (n:offset+nrh,:) = this%sensors(i)%hx_rh(:,1:DEF_DA_ENS_NUM)
               this%r  (n:offset+nrh) = this%sensors(i)%rh%r
               offset = offset + nrh
            ENDIF
         ENDIF
      ENDDO

   END SUBROUTINE SYNOP_type_concat

!-----------------------------------------------------------------------------

   SUBROUTINE SYNOP_type_save (this, is_analysis)

!-----------------------------------------------------------------------------
! Save patch-space H(x) diagnostics. The shared hx_* arrays have one slot per
! observation stream; a '2MT' stream stores 2 m temperature in its slot and a
! '2MQ' stream stores 2 m relative humidity. Since tref_ens and rhref_ens are
! derived from the analysis variables, the analysis H(x) is available directly
! after the update.
!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(SYNOP_set_type), intent(in) :: this
      logical,               intent(in) :: is_analysis

!------------------------ Local Variables ------------------------------------
      integer :: i
      integer :: is
      integer :: np
      logical :: valid_hx

!-----------------------------------------------------------------------------
      IF (nsource <= 0) RETURN
      IF (numpatch <= 0) RETURN
      IF (.not. allocated(hx_f_ens)) RETURN
      IF (this%use_t  .and. .not. allocated(tref_ens))  RETURN
      IF (this%use_rh .and. .not. allocated(rhref_ens)) RETURN

      DO i = 1, size(this%cfgs)
         is = this%cfgs(i)%istream
         IF (is_analysis) THEN
            hx_a_ens(is,:,:) = spval
            hx_a    (is,:)   = spval
         ELSE
            hx_ol   (is,:)   = spval
            hx_f_ens(is,:,:) = spval
            hx_f    (is,:)   = spval
         ENDIF
      ENDDO

      DO i = 1, size(this%cfgs)
         is = this%cfgs(i)%istream
         DO np = 1, numpatch
            IF (this%cfgs(i)%use_t) THEN
               valid_hx = all(tref_ens(1:DEF_DA_ENS_NUM,np) /= spval)

               IF (is_analysis) THEN
                  hx_a_ens(is,1:DEF_DA_ENS_NUM,np) = tref_ens(1:DEF_DA_ENS_NUM,np)
                  IF (valid_hx) hx_a(is,np) = &
                     sum(tref_ens(1:DEF_DA_ENS_NUM,np)) / DEF_DA_ENS_NUM
               ELSE
                  hx_ol(is,np) = tref_ens(0,np)
                  hx_f_ens(is,1:DEF_DA_ENS_NUM,np) = tref_ens(1:DEF_DA_ENS_NUM,np)
                  IF (valid_hx) hx_f(is,np) = &
                     sum(tref_ens(1:DEF_DA_ENS_NUM,np)) / DEF_DA_ENS_NUM
               ENDIF
            ELSE IF (this%cfgs(i)%use_rh) THEN
               valid_hx = all(rhref_ens(1:DEF_DA_ENS_NUM,np) /= spval)

               IF (is_analysis) THEN
                  hx_a_ens(is,1:DEF_DA_ENS_NUM,np) = rhref_ens(1:DEF_DA_ENS_NUM,np)
                  IF (valid_hx) hx_a(is,np) = &
                     sum(rhref_ens(1:DEF_DA_ENS_NUM,np)) / DEF_DA_ENS_NUM
               ELSE
                  hx_ol(is,np) = rhref_ens(0,np)
                  hx_f_ens(is,1:DEF_DA_ENS_NUM,np) = rhref_ens(1:DEF_DA_ENS_NUM,np)
                  IF (valid_hx) hx_f(is,np) = &
                     sum(rhref_ens(1:DEF_DA_ENS_NUM,np)) / DEF_DA_ENS_NUM
               ENDIF
            ENDIF
         ENDDO
      ENDDO

   END SUBROUTINE SYNOP_type_save

!-----------------------------------------------------------------------------

   SUBROUTINE SYNOP_type_prepare_on_pset (this, np, nobs_p, y_p, hx_p, r_p)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(SYNOP_set_type),  intent(in)  :: this
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

      ! find the observations around the target patch
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

   END SUBROUTINE SYNOP_type_prepare_on_pset

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Obs_SYNOP
#endif
