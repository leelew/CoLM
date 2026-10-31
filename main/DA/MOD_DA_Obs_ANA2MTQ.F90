#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Obs_ANA2MTQ
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    2m t/q analysis means as pseudo-observations
!
!    Use for land-atmosphere coupling assimilation
!
! HISTORY:
!    Lu Li, 10/2026: First operational version
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_DA_Obs_PseudoType
   USE MOD_Namelist
   USE MOD_DA_Vars_TimeVariables
   IMPLICIT NONE

!-----------------------------------------------------------------------------
! Screen-level pseudo-observation collection
!-----------------------------------------------------------------------------
   type :: ANA2MTQ_type

      type(pseudo_obs_type) :: t         ! 2 m temperature [K]
      type(pseudo_obs_type) :: q         ! 2 m specific humidity [kg/kg]

      ! Concatenated temperature and humidity records for assimilation.
      integer :: nobs = 0                ! t%nobs + q%nobs
      real(r8), allocatable :: lat(:)    ! [nobs], degrees
      real(r8), allocatable :: lon(:)    ! [nobs], degrees
      real(r8), allocatable :: y(:)      ! [nobs]
      real(r8), allocatable :: r(:)      ! Error standard deviations [nobs]
      real(r8), allocatable :: hx(:,:)   ! Forecast H(x) [nobs, ens]

   CONTAINS

      procedure, PUBLIC  :: init   => ANA2MTQ_init
      procedure, PUBLIC  :: calcp  => ANA2MTQ_calc_on_pset
      procedure, PUBLIC  :: clear  => ANA2MTQ_clear
      procedure, PRIVATE :: concat => ANA2MTQ_concat

   END type ANA2MTQ_type

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE ANA2MTQ_init (this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(ANA2MTQ_type), intent(inout) :: this

!-----------------------------------------------------------------------------
      CALL this%clear()
      CALL this%t%init()
      CALL this%q%init()

   END SUBROUTINE ANA2MTQ_init

!-----------------------------------------------------------------------------

   SUBROUTINE ANA2MTQ_calc_on_pset (this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(ANA2MTQ_type), intent(inout) :: this

!-----------------------------------------------------------------------------
      CALL this%clear()
      ! Observation errors are standard deviations.
      CALL this%t%read(tref_ens_f, tref_ens_a, 1.0_r8)
      CALL this%q%read(qref_ens_f, qref_ens_a, 0.2_r8)
      CALL this%concat()

   END SUBROUTINE ANA2MTQ_calc_on_pset

!-----------------------------------------------------------------------------

   SUBROUTINE ANA2MTQ_concat (this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(ANA2MTQ_type), intent(inout) :: this
      integer :: nt, nq

!-----------------------------------------------------------------------------
      nt = this%t%nobs
      nq = this%q%nobs
      this%nobs = nt + nq
      IF (this%nobs == 0) RETURN

      allocate(this%lat(this%nobs), this%lon(this%nobs))
      allocate(this%y(this%nobs), this%r(this%nobs), this%hx(this%nobs,DEF_DA_ENS_NUM))
      IF (nt > 0) THEN
         this%lat(1:nt) = this%t%lat
         this%lon(1:nt) = this%t%lon
         this%y(1:nt) = this%t%y
         this%r(1:nt) = this%t%r
         this%hx(1:nt,:) = this%t%hx
      ENDIF
      IF (nq > 0) THEN
         this%lat(nt+1:nt+nq) = this%q%lat
         this%lon(nt+1:nt+nq) = this%q%lon
         this%y(nt+1:nt+nq) = this%q%y
         this%r(nt+1:nt+nq) = this%q%r
         this%hx(nt+1:nt+nq,:) = this%q%hx
      ENDIF

   END SUBROUTINE ANA2MTQ_concat

!-----------------------------------------------------------------------------

   SUBROUTINE ANA2MTQ_clear (this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(ANA2MTQ_type), intent(inout) :: this

!-----------------------------------------------------------------------------
      CALL this%t%clear()
      CALL this%q%clear()
      this%nobs = 0
      IF (allocated(this%lat)) deallocate(this%lat)
      IF (allocated(this%lon)) deallocate(this%lon)
      IF (allocated(this%y))   deallocate(this%y)
      IF (allocated(this%r))   deallocate(this%r)
      IF (allocated(this%hx)) deallocate(this%hx)

   END SUBROUTINE ANA2MTQ_clear

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Obs_ANA2MTQ
#endif
