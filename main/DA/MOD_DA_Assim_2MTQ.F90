#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Assim_2MTQ
!-----------------------------------------------------------------------------
! Screen-analysis driver interface. SYNOP analysis is not yet implemented.
! Products remain missing until valid, time-matched analyses are produced.
!-----------------------------------------------------------------------------
   USE MOD_Precision, ONLY: r8
   USE MOD_SPMD_Task, ONLY: p_is_master, p_is_worker
   USE MOD_DA_Vars_TimeVariables, ONLY: tref_f, qref_f, tref_a, qref_a
   USE MOD_Vars_Global, ONLY: spval
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: init_assim_2MTQ, run_assim_2MTQ, end_assim_2MTQ

CONTAINS

   SUBROUTINE init_assim_2MTQ ()
      IF (p_is_master) WRITE (*,'(A)') &
         '[CoLM-DA] 2MTQ analysis is not implemented; ANA2MTQ products remain missing.'
      ! TODO: initialize SYNOP observation streams and the screen analysis.
   END SUBROUTINE init_assim_2MTQ

   SUBROUTINE run_assim_2MTQ (idate, deltim)
      integer, intent(in) :: idate(3)
      real(r8), intent(in) :: deltim

      ! Invalidate the previous cycle even when no new observations exist.
      IF (p_is_worker) THEN
         IF (allocated(tref_f)) tref_f = spval
         IF (allocated(qref_f)) qref_f = spval
         IF (allocated(tref_a)) tref_a = spval
         IF (allocated(qref_a)) qref_a = spval
      ENDIF

      ! TODO: read/QC SYNOP at idate, analyse and fill valid patch products.
      ! Keep background tref_ens/qref_ens unchanged for downstream SM H(x).
      ! Only publish products matching the downstream analysis time.
   END SUBROUTINE run_assim_2MTQ

   SUBROUTINE end_assim_2MTQ ()
      ! TODO: finalize SYNOP streams. Patch products are owned and released
      ! by MOD_DA_Vars_TimeVariables, not this driver.
   END SUBROUTINE end_assim_2MTQ

END MODULE MOD_DA_Assim_2MTQ
#endif
