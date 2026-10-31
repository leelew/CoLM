#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Obs_ANA2MTQ
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    2 m temperature and relative humidity analysis means as pseudo-observations
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
! Patch-based pseudo-observations from SYNOP 2 m temperature and humidity analyses.
! Temperature and humidity select changed patches independently.
! Combined arrays store temperature records first, followed by humidity records.
!-----------------------------------------------------------------------------
   type :: ANA2MTQ_type

      type(pseudo_obs_type) :: tref      ! 2 m temperature pseudo-observations [K]
      type(pseudo_obs_type) :: rhref     ! 2 m relative humidity pseudo-observations [0-1]

      ! Combined records available on all ranks for land-surface assimilation.
      integer :: nobs = 0                ! Total record count: tref%nobs + rhref%nobs
      real(r8), allocatable :: lat(:)    ! Patch latitudes [nobs], degrees north
      real(r8), allocatable :: lon(:)    ! Patch longitudes [nobs], degrees east in [-180, 180)
      real(r8), allocatable :: y(:)      ! Analysis ensemble means [nobs]; K for T, fraction for RH
      real(r8), allocatable :: r(:)      ! Error variances [nobs]; K^2 for T, fraction squared for RH
      real(r8), allocatable :: hx(:,:)   ! Forecast values [nobs, DEF_DA_ENS_NUM]; same units as y

   CONTAINS

      procedure, PUBLIC  :: init   => ANA2MTQ_init         ! Initialize empty temperature and humidity collections
      procedure, PUBLIC  :: calcp  => ANA2MTQ_calc_on_pset ! Construct and communicate both collections, then combine
      procedure, PUBLIC  :: clear  => ANA2MTQ_clear        ! Clear both collections and deallocate combined arrays
      procedure, PRIVATE :: concat => ANA2MTQ_concat       ! Combine temperature and humidity records in this order

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
      CALL this%tref%init()
      CALL this%rhref%init()

   END SUBROUTINE ANA2MTQ_init

!-----------------------------------------------------------------------------

   SUBROUTINE ANA2MTQ_calc_on_pset (this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(ANA2MTQ_type), intent(inout) :: this

!-----------------------------------------------------------------------------
      CALL this%clear()
      CALL this%tref%read(tref_ens_f, tref_ens_a, 1.0_r8)
      CALL this%rhref%read(rhref_ens_f, rhref_ens_a, 0.04_r8**2)
      CALL this%concat()

   END SUBROUTINE ANA2MTQ_calc_on_pset

!-----------------------------------------------------------------------------

   SUBROUTINE ANA2MTQ_concat (this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(ANA2MTQ_type), intent(inout) :: this
      integer :: nt, nq

!-----------------------------------------------------------------------------
      nt = this%tref%nobs
      nq = this%rhref%nobs
      this%nobs = nt + nq
      IF (this%nobs == 0) RETURN

      allocate(this%lat(this%nobs), this%lon(this%nobs), this%y(this%nobs), this%r(this%nobs))
      allocate(this%hx(this%nobs,DEF_DA_ENS_NUM))
      IF (nt > 0) THEN
         this%lat(1:nt) = this%tref%lat
         this%lon(1:nt) = this%tref%lon
         this%y(1:nt) = this%tref%y
         this%r(1:nt) = this%tref%r
         this%hx(1:nt,:) = this%tref%hx
      ENDIF
      IF (nq > 0) THEN
         this%lat(nt+1:nt+nq) = this%rhref%lat
         this%lon(nt+1:nt+nq) = this%rhref%lon
         this%y(nt+1:nt+nq) = this%rhref%y
         this%r(nt+1:nt+nq) = this%rhref%r
         this%hx(nt+1:nt+nq,:) = this%rhref%hx
      ENDIF

   END SUBROUTINE ANA2MTQ_concat

!-----------------------------------------------------------------------------

   SUBROUTINE ANA2MTQ_clear (this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(ANA2MTQ_type), intent(inout) :: this

!-----------------------------------------------------------------------------
      CALL this%tref%clear()
      CALL this%rhref%clear()
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
