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
!    ANA2MTQ builds fixed patch-neighbor and worker communication indices at
!    initialization. Each time step exchanges the required T/q data before
!    prepare() appends valid neighbor rows. this%nobs excludes these rows.
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
   USE MOD_DA_Vars_TimeVariables
   USE MOD_Namelist
   USE MOD_Pixelset
   USE MOD_Vars_Global
   USE MOD_Vars_TimeInvariants
   USE MOD_WorkerPushData, ONLY: worker_pushdata_type, build_worker_pushdata, worker_push_data
   USE MOD_Utils, ONLY: quicksort
   USE MOD_SPMD_Task
   USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
   IMPLICIT NONE
   PRIVATE :: SM_patch_lat_bounds

   type :: SM_config_type

      character(len=16) :: source_name ! 'PM', 'GNOS', 'AM', 'situ', 'ANA2MTQ'
      character(len=32) :: sensor_name ! 'SMAP', 'SMOS', 'CYGNSS', 'FY3D', 'ISMN', etc.
      character(len=16) :: var_name    ! 'TB_H', 'TB_V', 'SM', 'DDM', '2MTQ'
      real(r8) :: fghz                 ! frequency in GHz
      integer  :: istream              ! index of this stream in the shared hx_* diagnostics

   END type SM_config_type

   type :: SM_patch_neighbors_type

      integer, allocatable :: idx(:)        ! Indices in the requested-patch cache
      real(r8), allocatable :: dist(:)      ! Distance to the target patch [km]

   END type SM_patch_neighbors_type

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
      !type(GNOS_set_type) :: gnos

      ! Fixed ANA2MTQ indices; retained by clear() between time steps.
      type(SM_patch_neighbors_type), allocatable :: ana_neighbors(:)
      type(worker_pushdata_type), allocatable :: ana_push
      integer :: nana_patch = 0             ! Unique patches requested by this worker
      real(r8) :: ana_r_t = 0.01_r8         ! Temperature error variance [K^2]
      real(r8) :: ana_r_q = 0.01_r8         ! Specific-humidity error variance [(kg/kg)^2]

      ! Current-step data at requested patches; member 0 is excluded.
      logical, allocatable :: ana_valid(:)
      real(r8), allocatable :: ana_y(:,:)    ! [T/q, requested patch]
      real(r8), allocatable :: ana_hx(:,:,:) ! [ensemble member, T/q, requested patch]

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
      procedure :: free_neighbors => SM_type_free_patch_neighbors
      procedure, PRIVATE :: build_neighbors => SM_type_build_patch_neighbors
      procedure, PRIVATE :: exchange_neighbors => SM_type_exchange_patch_data

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
      CALL this%free_neighbors()
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

   SUBROUTINE SM_type_build_patch_neighbors(this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(SM_type), intent(inout) :: this

      integer :: np, i, nglobal, istt, iend
      integer, allocatable :: counts(:), displs(:), ids_me(:), ids_req(:), order(:), cache_index(:)
      real(r8) :: lat_p, lon_p, radius
      real(r8), allocatable :: lat_me(:), lon_me(:), lat_all(:), lon_all(:), dlon(:), hav(:)
      logical, allocatable :: needed(:), keep(:)

!-----------------------------------------------------------------------------
      IF (.not. p_is_worker) RETURN
      IF (.not. ieee_is_finite(this%dres) .or. this%dres <= 0.0_r8) &
         ERROR STOP 'ANA2MTQ requires a positive patch search half-width'
      IF (.not. ieee_is_finite(this%locr) .or. this%locr <= 0.0_r8) &
         ERROR STOP 'ANA2MTQ requires a positive localization length'

      ! Global IDs are worker offsets plus local patch indices.
      allocate(counts(0:p_np_worker-1), displs(0:p_np_worker-1))
#ifdef USEMPI
      CALL mpi_allgather(numpatch, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, p_comm_worker, p_err)
#else
      counts(0) = numpatch
#endif
      displs(0) = 0
      DO i = 1, p_np_worker-1
         displs(i) = displs(i-1) + counts(i-1)
      ENDDO
      nglobal = sum(counts)
      allocate(ids_me(numpatch), lat_me(numpatch), lon_me(numpatch))
      lat_me = spval
      lon_me = spval
      DO np = 1, numpatch
         ids_me(np) = displs(p_iam_worker) + np
         IF (patchtype(np) >= 3) CYCLE
         IF (.not. ieee_is_finite(patchlatr(np)) .or. .not. ieee_is_finite(patchlonr(np))) CYCLE
         IF (patchlatr(np) == spval .or. patchlonr(np) == spval) CYCLE
         IF (abs(patchlatr(np)) > pi/2.0_r8) CYCLE
         lat_me(np) = patchlatr(np)
         lon_me(np) = modulo(patchlonr(np)+pi, 2.0_r8*pi) - pi
      ENDDO

      ! Coordinates are shared only while building the fixed neighbor lists.
      allocate(lat_all(nglobal), lon_all(nglobal), order(nglobal))
#ifdef USEMPI
      CALL mpi_allgatherv(lat_me, numpatch, MPI_REAL8, lat_all, counts, displs, &
         MPI_REAL8, p_comm_worker, p_err)
      CALL mpi_allgatherv(lon_me, numpatch, MPI_REAL8, lon_all, counts, displs, &
         MPI_REAL8, p_comm_worker, p_err)
#else
      lat_all = lat_me
      lon_all = lon_me
#endif
      order = [(i, i=1,nglobal)]
      CALL quicksort(nglobal, lat_all, order)
      lon_all = lon_all(order)

      allocate(this%ana_neighbors(numpatch), needed(nglobal))
      needed = .false.
      radius = this%dres*pi/180.0_r8
      DO np = 1, numpatch
         allocate(this%ana_neighbors(np)%idx(0), this%ana_neighbors(np)%dist(0))
         IF (lat_me(np) == spval) CYCLE
         lat_p = lat_me(np)
         lon_p = lon_me(np)
         CALL SM_patch_lat_bounds(lat_all, lat_p-radius, lat_p+radius, istt, iend)
         dlon = abs(lon_all(istt:iend) - lon_p)
         dlon = min(dlon, 2.0_r8*pi-dlon)
         keep = dlon < radius
         this%ana_neighbors(np)%idx = pack(order(istt:iend), keep)
         needed(this%ana_neighbors(np)%idx) = .true.

         hav = sin((lat_all(istt:iend)-lat_p)/2.0_r8)**2 + &
               cos(lat_all(istt:iend))*cos(lat_p)*sin(dlon/2.0_r8)**2
         this%ana_neighbors(np)%dist = &
            pack(2.0_r8*6378.1_r8*asin(sqrt(max(0.0_r8,min(1.0_r8,hav)))), keep)
      ENDDO

      ! Deduplicate requests across all local target patches.
      this%nana_patch = count(needed)
      ids_req = pack([(i, i=1,nglobal)], needed)
      allocate(cache_index(nglobal))
      cache_index = 0
      DO i = 1, this%nana_patch
         cache_index(ids_req(i)) = i
      ENDDO
      DO np = 1, numpatch
         this%ana_neighbors(np)%idx = cache_index(this%ana_neighbors(np)%idx)
      ENDDO
      allocate(this%ana_push)
      CALL build_worker_pushdata(numpatch, ids_me, this%nana_patch, ids_req, this%ana_push)

   END SUBROUTINE SM_type_build_patch_neighbors

!-----------------------------------------------------------------------------

   SUBROUTINE SM_patch_lat_bounds(lat, south, north, istt, iend)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      real(r8), intent(in) :: lat(:), south, north
      integer, intent(out) :: istt, iend
      integer :: left, right, mid

!-----------------------------------------------------------------------------
      left = 1
      right = size(lat) + 1
      DO WHILE (left < right)
         mid = left + (right-left)/2
         IF (lat(mid) <= south) THEN
            left = mid + 1
         ELSE
            right = mid
         ENDIF
      ENDDO
      istt = left
      right = size(lat) + 1
      DO WHILE (left < right)
         mid = left + (right-left)/2
         IF (lat(mid) < north) THEN
            left = mid + 1
         ELSE
            right = mid
         ENDIF
      ENDDO
      iend = left - 1

   END SUBROUTINE SM_patch_lat_bounds

!-----------------------------------------------------------------------------

   SUBROUTINE SM_type_exchange_patch_data(this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(SM_type), intent(inout) :: this

      integer :: np, iens, ivar
      real(r8) :: tref_mean_f, tref_mean_a, qref_mean_f, qref_mean_a
      logical :: ready
      logical, allocatable :: valid(:)
      real(r8), allocatable :: y_local(:,:), values(:)

!-----------------------------------------------------------------------------
      IF (.not. p_is_worker) RETURN
      IF (.not. allocated(this%ana_push)) ERROR STOP 'ANA2MTQ patch indices are not initialized'
      IF (DEF_DA_ENS_NUM < 1) ERROR STOP 'ANA2MTQ requires at least one ensemble member'

      ready = .true.
      IF (numpatch > 0) THEN
         ready = allocated(tref_ens_f) .and. allocated(tref_ens_a) .and. &
                 allocated(qref_ens_f) .and. allocated(qref_ens_a)
         IF (ready) THEN
            ready = lbound(tref_ens_f,1) == 0 .and. ubound(tref_ens_f,1) == DEF_DA_ENS_NUM .and. &
                    lbound(tref_ens_a,1) == 0 .and. ubound(tref_ens_a,1) == DEF_DA_ENS_NUM .and. &
                    lbound(qref_ens_f,1) == 0 .and. ubound(qref_ens_f,1) == DEF_DA_ENS_NUM .and. &
                    lbound(qref_ens_a,1) == 0 .and. ubound(qref_ens_a,1) == DEF_DA_ENS_NUM .and. &
                    size(tref_ens_f,2) == numpatch .and. size(tref_ens_a,2) == numpatch .and. &
                    size(qref_ens_f,2) == numpatch .and. size(qref_ens_a,2) == numpatch
         ENDIF
      ENDIF
#ifdef USEMPI
      CALL mpi_allreduce(MPI_IN_PLACE, ready, 1, MPI_LOGICAL, MPI_LAND, p_comm_worker, p_err)
#endif
      IF (.not. ready) ERROR STOP 'ANA2MTQ requires forecast/analysis ensembles [0:nens, local patch]'

      allocate(valid(numpatch), y_local(2,numpatch), values(numpatch))
      valid = .false.
      y_local = spval
      DO np = 1, numpatch
         IF (patchtype(np) >= 3) CYCLE
         IF (.not. all(ieee_is_finite(tref_ens_f(1:DEF_DA_ENS_NUM,np)))) CYCLE
         IF (.not. all(ieee_is_finite(tref_ens_a(1:DEF_DA_ENS_NUM,np)))) CYCLE
         IF (.not. all(ieee_is_finite(qref_ens_f(1:DEF_DA_ENS_NUM,np)))) CYCLE
         IF (.not. all(ieee_is_finite(qref_ens_a(1:DEF_DA_ENS_NUM,np)))) CYCLE
         IF (.not. all(tref_ens_f(1:DEF_DA_ENS_NUM,np) /= spval)) CYCLE
         IF (.not. all(tref_ens_a(1:DEF_DA_ENS_NUM,np) /= spval)) CYCLE
         IF (.not. all(qref_ens_f(1:DEF_DA_ENS_NUM,np) /= spval)) CYCLE
         IF (.not. all(qref_ens_a(1:DEF_DA_ENS_NUM,np) /= spval)) CYCLE
         tref_mean_f = sum(tref_ens_f(1:DEF_DA_ENS_NUM,np)) / DEF_DA_ENS_NUM
         tref_mean_a = sum(tref_ens_a(1:DEF_DA_ENS_NUM,np)) / DEF_DA_ENS_NUM
         qref_mean_f = sum(qref_ens_f(1:DEF_DA_ENS_NUM,np)) / DEF_DA_ENS_NUM
         qref_mean_a = sum(qref_ens_a(1:DEF_DA_ENS_NUM,np)) / DEF_DA_ENS_NUM
         ! Retain paired selection with nonzero T and q mean analysis increments.
         IF (tref_mean_f == tref_mean_a .or. qref_mean_f == qref_mean_a) CYCLE
         valid(np) = .true.
         y_local(:,np) = [tref_mean_a, qref_mean_a]
      ENDDO

      allocate(this%ana_y(2,this%nana_patch), this%ana_hx(DEF_DA_ENS_NUM,2,this%nana_patch))
      allocate(this%ana_valid(this%nana_patch))
      DO ivar = 1, 2
         CALL worker_push_data(this%ana_push, y_local(ivar,:), this%ana_y(ivar,:), spval)
      ENDDO
      DO iens = 1, DEF_DA_ENS_NUM
         values = spval
         IF (numpatch > 0) THEN
            WHERE (valid) values = tref_ens_f(iens,:)
         ENDIF
         CALL worker_push_data(this%ana_push, values, this%ana_hx(iens,1,:), spval)
         values = spval
         IF (numpatch > 0) THEN
            WHERE (valid) values = qref_ens_f(iens,:)
         ENDIF
         CALL worker_push_data(this%ana_push, values, this%ana_hx(iens,2,:), spval)
      ENDDO
      this%ana_valid = all(this%ana_y /= spval, dim=1) .and. &
                      all(all(this%ana_hx /= spval, dim=1), dim=1)

   END SUBROUTINE SM_type_exchange_patch_data

!-----------------------------------------------------------------------------

   SUBROUTINE SM_type_free_patch_neighbors(this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(SM_type), intent(inout) :: this

!-----------------------------------------------------------------------------
      IF (allocated(this%ana_push))      deallocate(this%ana_push)
      IF (allocated(this%ana_neighbors)) deallocate(this%ana_neighbors)
      IF (allocated(this%ana_valid))     deallocate(this%ana_valid)
      IF (allocated(this%ana_y))         deallocate(this%ana_y)
      IF (allocated(this%ana_hx))        deallocate(this%ana_hx)
      this%nana_patch = 0

   END SUBROUTINE SM_type_free_patch_neighbors

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
      IF (allocated(this%ana_valid)) deallocate(this%ana_valid)
      IF (allocated(this%ana_y))     deallocate(this%ana_y)
      IF (allocated(this%ana_hx))    deallocate(this%ana_hx)

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
      IF (this%use_ana2mtq) CALL this%exchange_neighbors()

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
      integer :: iens, iobs, ip, icache, n, k
      real(r8) :: lat_p, lon_p, dlat, dlon
      logical,  allocatable :: valid(:)
      logical,  allocatable :: valid_ana2mtq(:)
      real(r8), allocatable :: lat(:)
      real(r8), allocatable :: lon(:)
      real(r8), allocatable :: dist(:)

!-----------------------------------------------------------------------------
      nobs_p = 0
      lat_p = patchlatr(np)*180/pi
      lon_p = patchlonr(np)*180/pi

      ! Select nearby source observations.
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

      ! Geometry is fixed; only current-step data validity is checked here.
      IF (this%use_ana2mtq) THEN
         IF (.not. allocated(this%ana_valid)) ERROR STOP 'ANA2MTQ data must be exchanged before prepare'
         valid_ana2mtq = this%ana_valid(this%ana_neighbors(np)%idx)
         nobs_p = nobs_p + 2*count(valid_ana2mtq)
      ELSE
         allocate(valid_ana2mtq(0))
      ENDIF
      IF (nobs_p == 0) RETURN

      allocate(y_p(nobs_p))
      allocate(hx_p(nobs_p, DEF_DA_ENS_NUM))
      allocate(r_p(nobs_p))
      allocate(lat(n))
      allocate(lon(n))
      allocate(dist(nobs_p))

      IF (n > 0) THEN
         y_p(1:n) = pack(this%y, valid)
         DO iens = 1, DEF_DA_ENS_NUM
            hx_p(1:n,iens) = pack(this%hx(:,iens), valid)
         ENDDO
         r_p(1:n) = pack(this%r, valid)
         lat(1:n) = pack(this%lat, valid)
         lon(1:n) = pack(this%lon, valid)
         dist(1:n) = 2.0_r8*6378.1_r8*asin(sqrt( &
            sin((lat*pi/180.0_r8 - patchlatr(np))/2.0_r8)**2 + &
            cos(lat*pi/180.0_r8)*cos(patchlatr(np))* &
            sin((lon*pi/180.0_r8 - patchlonr(np))/2.0_r8)**2))
      ENDIF

      k = n
      DO ip = 1, size(valid_ana2mtq)
         IF (.not. valid_ana2mtq(ip)) CYCLE
         icache = this%ana_neighbors(np)%idx(ip)
         y_p(k+1:k+2) = this%ana_y(:,icache)
         hx_p(k+1,:) = this%ana_hx(:,1,icache)
         hx_p(k+2,:) = this%ana_hx(:,2,icache)
         r_p(k+1) = this%ana_r_t
         r_p(k+2) = this%ana_r_q
         dist(k+1:k+2) = this%ana_neighbors(np)%dist(ip)
         k = k + 2
      ENDDO

      ! Apply the same distance-dependent R inflation to all observations.
      r_p = r_p / max(exp(-(dist**2)/(2.0*this%locr**2)), 1.0e-6)

   END SUBROUTINE SM_type_prepare_on_pset

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Obs_SM
#endif
