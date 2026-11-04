#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Assim_DISC
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Route ensemble-mean land runoff and accumulate DA discharge diagnostics.
!    Retain the DA river state between calls and restore the ordinary river
!    state after the additional routing calculation.
!
!    The current scheme propagates land-analysis information through runoff;
!    it does not directly assimilate discharge observations.
!-----------------------------------------------------------------------------
   USE MOD_Precision, ONLY: r8
#ifdef GridRiverLakeFlow
   USE MOD_SPMD_Task, ONLY: p_is_worker
   USE MOD_Namelist, ONLY: DEF_DA_ENS_NUM, DEF_Reservoir_Method
   USE MOD_LandPatch, ONLY: numpatch
   USE MOD_Vars_Global, ONLY: spval
   USE MOD_Vars_1DFluxes, ONLY: rnof
   USE MOD_Grid_RiverLakeNetwork, ONLY: numucat
   USE MOD_Grid_RiverLakeFlow, ONLY: grid_riverlake_flow
   USE MOD_Grid_RiverLakeTimeVars, ONLY: wdsrf_ucat, veloc_riv, volresv
   USE MOD_Grid_RiverLakeHistState, ONLY: a_discharge, acctime_ucat
   USE MOD_DA_Vars_TimeVars, ONLY: rnof_ens, wdsrf_ucat_a, veloc_riv_a, volresv_a, &
      wdsrf_ucat_ol, veloc_riv_ol, volresv_ol
   USE MOD_DA_Vars_AccVars, ONLY: &
      a_discharge_ol, acctime_ucat_ol, a_discharge_a, acctime_ucat_a
#endif
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: run_assim_DISC

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE run_assim_DISC (idate, deltim)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      integer, intent(in) :: idate(3)  ! Model time: year, Julian day, seconds of day
      real(r8), intent(in) :: deltim   ! Land-model time step [s]

!-----------------------------------------------------------------------------
#ifdef GridRiverLakeFlow
      IF (p_is_worker) THEN
         ! Back up the ordinary river state and history accumulators.
         IF (allocated(wdsrf_ucat_ol)) THEN
            wdsrf_ucat_ol = wdsrf_ucat
            veloc_riv_ol = veloc_riv
         ENDIF
         IF (DEF_Reservoir_Method > 0) THEN
            IF (allocated(volresv_ol)) volresv_ol = volresv
         ENDIF
         IF (numucat > 0) THEN
            a_discharge_ol  = a_discharge
            acctime_ucat_ol = acctime_ucat
         ENDIF

         ! Ensemble-mean forecast runoff for river routing
         IF (numpatch > 0) &
            rnof = sum(rnof_ens(1:DEF_DA_ENS_NUM,:), dim=1) / DEF_DA_ENS_NUM

         ! Load the DA river state retained from the previous time step.
         IF (allocated(wdsrf_ucat_a)) THEN
            WHERE (wdsrf_ucat_a /= spval) wdsrf_ucat = wdsrf_ucat_a
            WHERE (veloc_riv_a  /= spval) veloc_riv  = veloc_riv_a
         ENDIF
         IF (DEF_Reservoir_Method > 0) THEN
            IF (allocated(volresv_a)) THEN
               WHERE (volresv_a /= spval) volresv = volresv_a
            ENDIF
         ENDIF
      ENDIF

      CALL grid_riverlake_flow (idate(1), deltim)

      IF (p_is_worker) THEN
         ! Accumulate only the discharge and time contributed by this routing call.
         IF (numucat > 0) THEN
            a_discharge_a  = a_discharge_a  + (a_discharge  - a_discharge_ol )
            acctime_ucat_a = acctime_ucat_a + (acctime_ucat - acctime_ucat_ol)
         ENDIF

         ! Retain the updated DA river state for the next time step.
         IF (allocated(wdsrf_ucat_a)) THEN
            wdsrf_ucat_a = wdsrf_ucat
            veloc_riv_a  = veloc_riv
         ENDIF
         IF (DEF_Reservoir_Method > 0) THEN
            IF (allocated(volresv_a)) volresv_a = volresv
         ENDIF
      ENDIF

      ! Restore the ordinary river state and discharge accumulators.
      IF (p_is_worker) THEN
         IF (allocated(wdsrf_ucat_ol)) THEN
            wdsrf_ucat = wdsrf_ucat_ol
            veloc_riv = veloc_riv_ol
         ENDIF
         IF (DEF_Reservoir_Method > 0) THEN
            IF (allocated(volresv_ol)) volresv = volresv_ol
         ENDIF
         IF (numucat > 0) THEN
            a_discharge  = a_discharge_ol
            acctime_ucat = acctime_ucat_ol
         ENDIF
      ENDIF
#endif

   END SUBROUTINE run_assim_DISC

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Assim_DISC
#endif
