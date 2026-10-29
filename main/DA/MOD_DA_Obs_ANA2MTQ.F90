#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Obs_ANA2MTQ
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    SYNOP analysis means as pseudo-observations for soil-moisture assimilation.
!    Select changed patches, gather on master, and broadcast paired T/q data.
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_DA_Vars_TimeVariables, ONLY: tref_ens_f, tref_ens_a, qref_ens_f, qref_ens_a
   USE MOD_LandPatch, ONLY: numpatch
   USE MOD_Namelist, ONLY: DEF_DA_ENS_NUM
   USE MOD_SPMD_Task
   USE MOD_Vars_Global, ONLY: pi, spval
   USE MOD_Vars_TimeInvariants, ONLY: patchlatr, patchlonr, patchtype
   USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: ANA2MTQ_var_type, ANA2MTQ_type

   type :: ANA2MTQ_var_type

      real(r8) :: r = 0.01_r8               ! Error variance in squared variable units
      real(r8), allocatable :: y(:)          ! Analysis means [nobs]
      real(r8), allocatable :: hx(:,:)       ! Forecast H(x) [nobs, DEF_DA_ENS_NUM]

   END type ANA2MTQ_var_type

!-----------------------------------------------------------------------------
! Screen-level pseudo-observation collection
!-----------------------------------------------------------------------------
   type :: ANA2MTQ_type

      integer :: nobs = 0                 ! Number of patches with nonzero T and q analysis increments
      real(r8), allocatable :: lat(:)        ! Patch latitude [nobs], degrees
      real(r8), allocatable :: lon(:)        ! Patch longitude [nobs], degrees
      type(ANA2MTQ_var_type) :: t            ! Temperature [K], r in K^2
      type(ANA2MTQ_var_type) :: q            ! Specific humidity [kg/kg], r in (kg/kg)^2

   CONTAINS

      procedure, PUBLIC :: init  => ANA2MTQ_init
      procedure, PUBLIC :: calcp => ANA2MTQ_calc_on_pset
      procedure, PUBLIC :: clear => ANA2MTQ_clear

   END type ANA2MTQ_type

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE ANA2MTQ_init (this, r_t, r_q)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(ANA2MTQ_type), intent(inout) :: this
      real(r8), optional, intent(in) :: r_t, r_q

!-----------------------------------------------------------------------------
      CALL this%clear()
      IF (present(r_t)) this%t%r = r_t
      IF (present(r_q)) this%q%r = r_q
      IF (.not. ieee_is_finite(this%t%r) .or. this%t%r <= 0.0_r8) &
         ERROR STOP 'ANA2MTQ requires a positive temperature error variance'
      IF (.not. ieee_is_finite(this%q%r) .or. this%q%r <= 0.0_r8) &
         ERROR STOP 'ANA2MTQ requires a positive humidity error variance'

   END SUBROUTINE ANA2MTQ_init

!-----------------------------------------------------------------------------

   SUBROUTINE ANA2MTQ_calc_on_pset (this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(ANA2MTQ_type), intent(inout) :: this

!------------------------ Local Variables ------------------------------------
      integer :: np, nlocal_patch, nlocal, ipatch
      integer :: nfield
      real(r8) :: tref_mean_f, tref_mean_a, qref_mean_f, qref_mean_a
      logical :: ready
      logical, allocatable :: valid(:)
      real(r8), allocatable :: local_data(:,:), global_data(:,:)
#ifdef USEMPI
      integer :: iproc
      integer, allocatable :: counts(:), displs(:)
#endif

!-----------------------------------------------------------------------------
      CALL this%clear()
      IF (DEF_DA_ENS_NUM < 1) ERROR STOP 'ANA2MTQ requires at least one ensemble member'
      nfield = 4 + 2*DEF_DA_ENS_NUM
      nlocal_patch = 0
      IF (p_is_worker) nlocal_patch = numpatch

      ! Shared ensembles retain member 0; assimilation uses members 1:nens.
      ready = .true.
      IF (nlocal_patch > 0) THEN
         ready = allocated(tref_ens_f) .and. allocated(tref_ens_a) .and. &
                 allocated(qref_ens_f) .and. allocated(qref_ens_a)
         IF (ready) THEN
            ready = lbound(tref_ens_f,1) == 0 .and. ubound(tref_ens_f,1) == DEF_DA_ENS_NUM .and. &
                    lbound(tref_ens_a,1) == 0 .and. ubound(tref_ens_a,1) == DEF_DA_ENS_NUM .and. &
                    lbound(qref_ens_f,1) == 0 .and. ubound(qref_ens_f,1) == DEF_DA_ENS_NUM .and. &
                    lbound(qref_ens_a,1) == 0 .and. ubound(qref_ens_a,1) == DEF_DA_ENS_NUM .and. &
                    size(tref_ens_f,2) == nlocal_patch .and. size(tref_ens_a,2) == nlocal_patch .and. &
                    size(qref_ens_f,2) == nlocal_patch .and. size(qref_ens_a,2) == nlocal_patch
         ENDIF
      ENDIF
#ifdef USEMPI
      CALL mpi_allreduce(MPI_IN_PLACE, ready, 1, MPI_LOGICAL, MPI_LAND, p_comm_glb, p_err)
#endif
      IF (.not. ready) ERROR STOP 'ANA2MTQ requires forecast/analysis ensembles [0:nens, local patch]'

      allocate(valid(nlocal_patch))
      valid = .false.
      DO np = 1, nlocal_patch
         IF (patchtype(np) >= 3) CYCLE
         IF (.not. ieee_is_finite(patchlatr(np)) .or. .not. ieee_is_finite(patchlonr(np))) CYCLE
         IF (patchlatr(np) == spval .or. patchlonr(np) == spval) CYCLE
         IF (abs(patchlatr(np)) > pi/2.0_r8) CYCLE
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
         valid(np) = tref_mean_f /= tref_mean_a .and. qref_mean_f /= qref_mean_a
      ENDDO

      nlocal = count(valid)
      allocate(local_data(nfield,max(1,nlocal)))
      local_data = spval
      ipatch = 0
      DO np = 1, nlocal_patch
         IF (.not. valid(np)) CYCLE
         ipatch = ipatch + 1
         local_data(1,ipatch) = patchlatr(np)*180.0_r8/pi
         local_data(2,ipatch) = (modulo(patchlonr(np)+pi,2.0_r8*pi)-pi)*180.0_r8/pi
         local_data(3,ipatch) = sum(tref_ens_a(1:DEF_DA_ENS_NUM,np)) / DEF_DA_ENS_NUM
         local_data(4,ipatch) = sum(qref_ens_a(1:DEF_DA_ENS_NUM,np)) / DEF_DA_ENS_NUM
         local_data(5:4+DEF_DA_ENS_NUM,ipatch) = tref_ens_f(1:DEF_DA_ENS_NUM,np)
         local_data(5+DEF_DA_ENS_NUM:nfield,ipatch) = qref_ens_f(1:DEF_DA_ENS_NUM,np)
      ENDDO

      ! Keep both variables and their forecast members together for each patch.
#ifdef USEMPI
      allocate(counts(0:p_np_glb-1), displs(0:p_np_glb-1))
      counts = 0
      displs = 0
      CALL mpi_gather(nlocal, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
      IF (p_is_master) this%nobs = sum(counts)
      CALL mpi_bcast(this%nobs, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
      IF (this%nobs == 0) RETURN

      IF (p_is_master) THEN
         counts = counts*nfield
         DO iproc = 1, p_np_glb-1
            displs(iproc) = displs(iproc-1) + counts(iproc-1)
         ENDDO
         allocate(global_data(nfield,this%nobs))
      ELSE
         allocate(global_data(nfield,1))
      ENDIF
      CALL mpi_gatherv(local_data, nfield*nlocal, MPI_REAL8, global_data, counts, displs, &
         MPI_REAL8, p_address_master, p_comm_glb, p_err)
#else
      this%nobs = nlocal
      IF (this%nobs == 0) RETURN
      allocate(global_data(nfield,this%nobs))
      global_data = local_data
#endif

      allocate(this%lat(this%nobs), this%lon(this%nobs))
      allocate(this%t%y(this%nobs), this%q%y(this%nobs))
      allocate(this%t%hx(this%nobs,DEF_DA_ENS_NUM), this%q%hx(this%nobs,DEF_DA_ENS_NUM))
      IF (p_is_master) THEN
         this%lat = global_data(1,:)
         this%lon = global_data(2,:)
         this%t%y = global_data(3,:)
         this%q%y = global_data(4,:)
         this%t%hx = transpose(global_data(5:4+DEF_DA_ENS_NUM,:))
         this%q%hx = transpose(global_data(5+DEF_DA_ENS_NUM:nfield,:))
      ENDIF
#ifdef USEMPI
      CALL mpi_bcast(this%lat, this%nobs, MPI_REAL8, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast(this%lon, this%nobs, MPI_REAL8, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast(this%t%y, this%nobs, MPI_REAL8, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast(this%q%y, this%nobs, MPI_REAL8, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast(this%t%hx, this%nobs*DEF_DA_ENS_NUM, MPI_REAL8, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast(this%q%hx, this%nobs*DEF_DA_ENS_NUM, MPI_REAL8, p_address_master, p_comm_glb, p_err)
#endif

   END SUBROUTINE ANA2MTQ_calc_on_pset

!-----------------------------------------------------------------------------

   SUBROUTINE ANA2MTQ_clear (this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(ANA2MTQ_type), intent(inout) :: this

!-----------------------------------------------------------------------------
      this%nobs = 0
      IF (allocated(this%lat)) deallocate(this%lat)
      IF (allocated(this%lon)) deallocate(this%lon)
      IF (allocated(this%t%y))  deallocate(this%t%y)
      IF (allocated(this%q%y))  deallocate(this%q%y)
      IF (allocated(this%t%hx)) deallocate(this%t%hx)
      IF (allocated(this%q%hx)) deallocate(this%q%hx)

   END SUBROUTINE ANA2MTQ_clear

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Obs_ANA2MTQ
#endif
