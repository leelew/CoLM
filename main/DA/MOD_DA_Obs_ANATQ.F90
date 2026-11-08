#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Obs_ANATQ
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    2 m temperature and relative humidity analysis means as pseudo-observations
!
!    Each ANATQ_type handles one variable and stores its pseudo-observations
!    and observation-space H(x). ANATQ_set_type combines the fixed temperature
!    and relative humidity collections for land-surface assimilation.
!
! HISTORY:
!    Lu Li, 11/2026: DA v1.0
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_LandPatch
   USE MOD_SPMD_Task
   USE MOD_Vars_Global
   USE MOD_Vars_TimeInvariants
   USE MOD_DA_Obs_PseudoType
   USE MOD_DA_Assim_TQ
   IMPLICIT NONE

!-----------------------------------------------------------------------------
! Patch-based pseudo-observations from SYNOP 2 m temperature and humidity analyses.
! Temperature and humidity select changed patches independently.
! Combined arrays store temperature records first, followed by humidity records.
! Patches are both model and observation locations; no spatial mapping is needed.
!-----------------------------------------------------------------------------
   type :: ANATQ_type

      type(pseudo_obs_type) :: y                    ! Pseudo-observations and metadata for one variable
      real(r8), allocatable :: hx_obs(:,:)          ! Forecast values at y locations [nobs, ens]

   CONTAINS

      procedure :: init     => ANATQ_init           ! Initialize an empty observation collection
      procedure :: calc_hx => ANATQ_calc_hx         ! Construct y and matching observation-space H(x)
      procedure :: clear    => ANATQ_clear          ! Clear observations and H(x) for this variable

   END type ANATQ_type

!-----------------------------------------------------------------------------
! Fixed temperature and humidity collections and their combined records
!-----------------------------------------------------------------------------
   type :: ANATQ_set_type

      type(ANATQ_type) :: tref                      ! 2 m temperature pseudo-observations [K]
      type(ANATQ_type) :: rhref                     ! 2 m relative humidity pseudo-observations [0-1]

      integer :: nobs = 0                           ! Total record count: tref%y%nobs + rhref%y%nobs
      real(r8), allocatable :: lat(:)               ! Patch latitudes [nobs], degrees
      real(r8), allocatable :: lon(:)               ! Patch longitudes [nobs], degrees
      real(r8), allocatable :: y(:)                 ! Analysis ensemble means [nobs]
      real(r8), allocatable :: r(:)                 ! Observation error standard deviations [nobs]
      real(r8), allocatable :: hx_obs(:,:)          ! Forecast values at y locations [nobs, ens]

   CONTAINS

      procedure :: init     => ANATQ_set_init       ! Initialize the temperature and humidity collections
      procedure :: calc_hx => ANATQ_set_calc_hx     ! Prepare both variables and combine their records
      procedure :: concat   => ANATQ_set_concat     ! Combine temperature records followed by humidity
      procedure :: clear    => ANATQ_set_clear      ! Clear both variables and their combined arrays

   END type ANATQ_set_type

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE ANATQ_init (this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(ANATQ_type), intent(inout) :: this

!-----------------------------------------------------------------------------
      CALL this%clear()
      CALL this%y%init()

   END SUBROUTINE ANATQ_init

!-----------------------------------------------------------------------------

   SUBROUTINE ANATQ_calc_hx (this, var_name)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(ANATQ_type), intent(inout) :: this
      character(len=*), intent(in) :: var_name

!------------------------ Local Variables ------------------------------------
      real(r8) :: r                                 ! Observation error standard deviation
      integer :: iens                               ! Ensemble member index
      integer :: iproc                              ! Global process rank
      integer :: nlocal                             ! Selected-patch count on this rank
      integer :: nobs                               ! Global observation count for this variable
      logical, allocatable :: mask(:)               ! Local patches with changed analysis means
      real(r8), allocatable :: y(:)                 ! Analysis means at patch observation locations
      real(r8), allocatable :: hx_obs(:,:)          ! Local observation-space H(x) on each worker [ens, nlocal]
      real(r8), allocatable :: hx_recv(:,:)         ! Gathered H(x) receive buffer [ens, nobs]
      integer, allocatable :: counts(:)             ! Selected-patch counts by rank
      integer, allocatable :: displs(:)             ! Receive offsets in patch units

!-----------------------------------------------------------------------------
      CALL this%clear()
      IF (p_is_worker) THEN
         allocate(mask(numpatch), y(numpatch))
      ELSE
         allocate(mask(0), y(0))
      ENDIF
      allocate(counts(0:p_np_glb-1), displs(0:p_np_glb-1))
      mask = .false.
      y = spval

      ! Select changed patches and calculate their analysis means.
      SELECT CASE (trim(var_name))
      CASE ('2MT')
         r = 1.0_r8
         IF (p_is_worker) THEN
            IF (allocated(tref_ens_f) .and. allocated(tref_ens_a)) THEN
               IF (.not. any(tref_ens_f == spval) .and. .not. any(tref_ens_a == spval)) THEN
                  y = sum(tref_ens_a(1:DEF_DA_ENS_NUM,:), dim=1) / DEF_DA_ENS_NUM
                  mask = patchtype < 3 .and. &
                     y /= sum(tref_ens_f(1:DEF_DA_ENS_NUM,:), dim=1) / DEF_DA_ENS_NUM
               ENDIF
            ENDIF
         ENDIF
      CASE ('2MRH')
         r = 0.04_r8
         IF (p_is_worker) THEN
            IF (allocated(rhref_ens_f) .and. allocated(rhref_ens_a)) THEN
               IF (.not. any(rhref_ens_f == spval) .and. .not. any(rhref_ens_a == spval)) THEN
                  y = sum(rhref_ens_a(1:DEF_DA_ENS_NUM,:), dim=1) / DEF_DA_ENS_NUM
                  mask = patchtype < 3 .and. &
                     y /= sum(rhref_ens_f(1:DEF_DA_ENS_NUM,:), dim=1) / DEF_DA_ENS_NUM
               ENDIF
            ENDIF
         ENDIF
      CASE DEFAULT
         RETURN
      END SELECT

      ! PseudoType handles observations only; all ranks call read.
      CALL this%y%read(mask, y, r)
      nobs = this%y%nobs
      IF (nobs == 0) RETURN

      ! Pack H(x) with the same mask and patch order used for y.
      nlocal = count(mask)
      allocate(hx_obs(DEF_DA_ENS_NUM,nlocal))
      IF (nlocal > 0) THEN
         SELECT CASE (trim(var_name))
         CASE ('2MT')
            DO iens = 1, DEF_DA_ENS_NUM
               hx_obs(iens,:) = pack(tref_ens_f(iens,:), mask)
            ENDDO
         CASE ('2MRH')
            DO iens = 1, DEF_DA_ENS_NUM
               hx_obs(iens,:) = pack(rhref_ens_f(iens,:), mask)
            ENDDO
         END SELECT
      ENDIF

      ! Gather H(x) in rank order, matching the observation records.
      counts = 0
      displs = 0
      CALL mpi_gather(nlocal, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, &
         p_address_master, p_comm_glb, p_err)
      IF (p_is_master) THEN
         DO iproc = 1, p_np_glb-1
            displs(iproc) = displs(iproc-1) + counts(iproc-1)
         ENDDO
         allocate(hx_recv(DEF_DA_ENS_NUM,nobs))
      ELSE
         allocate(hx_recv(DEF_DA_ENS_NUM,1))
      ENDIF
      CALL mpi_gatherv(hx_obs, DEF_DA_ENS_NUM*nlocal, MPI_REAL8, hx_recv, &
         DEF_DA_ENS_NUM*counts, DEF_DA_ENS_NUM*displs, MPI_REAL8, &
         p_address_master, p_comm_glb, p_err)

      ! Store the complete observation-space H(x) on every rank.
      allocate(this%hx_obs(nobs,DEF_DA_ENS_NUM))
      IF (p_is_master) this%hx_obs = transpose(hx_recv)
      CALL mpi_bcast(this%hx_obs, nobs*DEF_DA_ENS_NUM, MPI_REAL8, &
         p_address_master, p_comm_glb, p_err)
      deallocate(hx_obs, hx_recv)

   END SUBROUTINE ANATQ_calc_hx

!-----------------------------------------------------------------------------

   SUBROUTINE ANATQ_clear (this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(ANATQ_type), intent(inout) :: this

!-----------------------------------------------------------------------------
      CALL this%y%clear()
      IF (allocated(this%hx_obs)) deallocate(this%hx_obs)

   END SUBROUTINE ANATQ_clear

!-----------------------------------------------------------------------------

   SUBROUTINE ANATQ_set_init (this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(ANATQ_set_type), intent(inout) :: this

!-----------------------------------------------------------------------------
      CALL this%clear()
      CALL this%tref%init()
      CALL this%rhref%init()

   END SUBROUTINE ANATQ_set_init

!-----------------------------------------------------------------------------

   SUBROUTINE ANATQ_set_calc_hx (this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(ANATQ_set_type), intent(inout) :: this

!-----------------------------------------------------------------------------
      CALL this%clear()
      CALL this%tref%calc_hx('2MT')
      CALL this%rhref%calc_hx('2MRH')
      CALL this%concat()

   END SUBROUTINE ANATQ_set_calc_hx

!-----------------------------------------------------------------------------

   SUBROUTINE ANATQ_set_concat (this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(ANATQ_set_type), intent(inout) :: this
      integer :: nt, nq

!-----------------------------------------------------------------------------
      nt = this%tref%y%nobs
      nq = this%rhref%y%nobs
      this%nobs = nt + nq
      IF (this%nobs == 0) RETURN

      allocate(this%lat(this%nobs), this%lon(this%nobs), this%y(this%nobs), this%r(this%nobs))
      allocate(this%hx_obs(this%nobs,DEF_DA_ENS_NUM))
      IF (nt > 0) THEN
         this%lat(1:nt) = this%tref%y%lat
         this%lon(1:nt) = this%tref%y%lon
         this%y(1:nt) = this%tref%y%y
         this%r(1:nt) = this%tref%y%r
         this%hx_obs(1:nt,:) = this%tref%hx_obs
      ENDIF
      IF (nq > 0) THEN
         this%lat(nt+1:nt+nq) = this%rhref%y%lat
         this%lon(nt+1:nt+nq) = this%rhref%y%lon
         this%y(nt+1:nt+nq) = this%rhref%y%y
         this%r(nt+1:nt+nq) = this%rhref%y%r
         this%hx_obs(nt+1:nt+nq,:) = this%rhref%hx_obs
      ENDIF

   END SUBROUTINE ANATQ_set_concat

!-----------------------------------------------------------------------------

   SUBROUTINE ANATQ_set_clear (this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(ANATQ_set_type), intent(inout) :: this

!-----------------------------------------------------------------------------
      CALL this%tref%clear()
      CALL this%rhref%clear()
      this%nobs = 0
      IF (allocated(this%lat)) deallocate(this%lat)
      IF (allocated(this%lon)) deallocate(this%lon)
      IF (allocated(this%y))   deallocate(this%y)
      IF (allocated(this%r))   deallocate(this%r)
      IF (allocated(this%hx_obs)) deallocate(this%hx_obs)

   END SUBROUTINE ANATQ_set_clear

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Obs_ANATQ
#endif
