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
!   No discharge observations are assimilated. The additional routing uses
!   forecast runoff; land-state analysis does not recompute runoff here.
!   The suffix _a denotes the separate DA routing trajectory, rather than
!   a direct discharge-observation analysis. Vars modules manage its storage.
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

   SUBROUTINE run_DA_DISC (idate, deltim)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      integer, intent(in) :: idate(3)
      real(r8), intent(in) :: deltim

!-----------------------------------------------------------------------------
      IF (.not. DEF_DA) RETURN

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

   END SUBROUTINE run_DA_DISC

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Assim_DISC
#endif
