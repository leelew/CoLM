#include <define.h>

#if defined(DataAssimilation) && defined(GridRiverLakeFlow)
MODULE MOD_DA_Assim_DISC
!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Discharge diagnostic driver using ensemble-mean forecast runoff.
!
! WORKFLOW:
!   1. Back up the ordinary river state and discharge accumulators.
!   2. Load the retained DA river state and ensemble-mean forecast runoff.
!   3. Advance river routing and accumulate discharge diagnostics.
!   4. Retain the updated DA river state and restore the ordinary river state.
!
! NOTES:
!   No discharge observations are assimilated.
!
! HISTORY:
!   Lu Li, 05/2026: Initial implementation
!   Lu Li, Zhilong Fan, 11/2026: DA v1.0
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_LandPatch
   USE MOD_Vars_Global
   USE MOD_Vars_1DFluxes
   USE MOD_Grid_RiverLakeNetwork
   USE MOD_Grid_RiverLakeFlow
   USE MOD_Grid_RiverLakeTimeVars
   USE MOD_Grid_RiverLakeHistState
   USE MOD_DA_Vars_TimeVars
   USE MOD_DA_Vars_AccVars
   IMPLICIT NONE

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE run_DA_DISC(idate, deltim)

!-----------------------------------------------------------------------------
      integer,  intent(in) :: idate(3)
      real(r8), intent(in) :: deltim
!-----------------------------------------------------------------------------
      IF (p_is_worker) THEN
         ! Back up ordinary catchment states and load the retained DA states.
         IF (numucat > 0) THEN
            wdsrf_ucat_ol   = wdsrf_ucat
            veloc_riv_ol    = veloc_riv
            a_discharge_ol  = a_discharge
            acctime_ucat_ol = acctime_ucat

            WHERE (wdsrf_ucat_a /= spval) wdsrf_ucat = wdsrf_ucat_a
            WHERE (veloc_riv_a  /= spval) veloc_riv  = veloc_riv_a
         ENDIF
         IF (DEF_Reservoir_Method > 0 .and. numresv > 0) THEN
            volresv_ol = volresv
            WHERE (volresv_a /= spval) volresv = volresv_a
         ENDIF

         ! Drive routing with ensemble-mean forecast runoff.
         IF (numpatch > 0) rnof = sum(rnof_ens(1:DEF_DA_ENS_NUM,:), dim=1) / DEF_DA_ENS_NUM
      ENDIF

      CALL grid_riverlake_flow(idate(1), deltim)

      IF (p_is_worker) THEN
         IF (numucat > 0) THEN
            ! Accumulate discharge and routing time from this integration.
            a_discharge_a  = a_discharge_a  + (a_discharge  - a_discharge_ol )
            acctime_ucat_a = acctime_ucat_a + (acctime_ucat - acctime_ucat_ol)

            ! Retain DA states and restore ordinary states and accumulators.
            wdsrf_ucat_a = wdsrf_ucat
            veloc_riv_a  = veloc_riv

            wdsrf_ucat   = wdsrf_ucat_ol
            veloc_riv    = veloc_riv_ol
            a_discharge  = a_discharge_ol
            acctime_ucat = acctime_ucat_ol
         ENDIF
         IF (DEF_Reservoir_Method > 0 .and. numresv > 0) THEN
            volresv_a = volresv
            volresv   = volresv_ol
         ENDIF
      ENDIF

   END SUBROUTINE run_DA_DISC

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Assim_DISC
#endif
