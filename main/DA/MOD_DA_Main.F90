#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Main
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Coordinate ensemble forecasts, data assimilation, and model-state updates.
!
! NOTES:
!    Enabled analyses follow a fixed order, independent of the order of
!    observation streams in the namelist:
!       2 m temperature -> 2 m relative humidity -> 
!       snow -> soil moisture -> soil & snow temperature.
!
!    After analysis, online mode applies the ensemble mean to the model,
!    while offline mode restores the unperturbed trajectory.
!
! AUTHOR:
!     Lu Li, 12/2024: Initial implemention
!     Zhilong Fan, Lu Li, 03/2024: Clean codes
!     Lu Li, Zhilong Fan, 11/2026: DA v1.0
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Spmd_task
   USE MOD_Namelist
   USE MOD_LandPatch
#ifdef GridRiverLakeFlow
   USE MOD_Vars_1DFluxes
   USE MOD_Vars_Global
   USE MOD_Grid_RiverLakeFlow
   USE MOD_Grid_RiverLakeNetwork
   USE MOD_Grid_RiverLakeTimeVars
   USE MOD_DA_Vars_AccVars
   USE MOD_Grid_RiverLakeHist
   USE MOD_DA_Assim_DISC
#endif
   USE MOD_DA_Assim_TQ
   USE MOD_DA_Assim_SM
   USE MOD_DA_Assim_ST
   USE MOD_DA_Ens
   IMPLICIT NONE
   SAVE

   logical :: use_assim_2mt  = .false.
   logical :: use_assim_2mq  = .false.
   logical :: use_assim_disc = .false.
   logical :: use_assim_sm   = .false.
   logical :: use_assim_st   = .false.

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE init_DA ()

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!-----------------------------------------------------------------------------

      use_assim_2mt  = any(DEF_DA_OBS_TARGET == '2MT')
      use_assim_2mq  = any(DEF_DA_OBS_TARGET == '2MQ')
#ifdef GridRiverLakeFlow
      use_assim_disc = .true.
#endif
      use_assim_sm   = any(DEF_DA_OBS_TARGET == 'SM' )
      use_assim_st   = any(DEF_DA_OBS_TARGET == 'ST' )

      IF (use_assim_2mt .or. use_assim_2mq) CALL init_DA_TQ ()
      IF (use_assim_sm) CALL init_DA_SM ()
      IF (use_assim_st) CALL init_DA_ST ()

   END SUBROUTINE init_DA

!-----------------------------------------------------------------------------

   SUBROUTINE run_DA (idate, deltim, dolai, doalb, dosst, oro, is_spinup)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      integer,  intent(in)    :: idate(3)
      real(r8), intent(in)    :: deltim
      logical,  intent(in)    :: dolai
      logical,  intent(in)    :: doalb
      logical,  intent(in)    :: dosst
      real(r8), intent(inout) :: oro(numpatch)
      logical,  intent(in)    :: is_spinup

      integer :: i

!-----------------------------------------------------------------------------
      ! Save the open-loop trajectory
      CALL save_member(0)

      ! Generate ensemble members
      IF (DEF_DA_ENS_NUM > 1) THEN
         CALL ensemble (deltim)
         DO i = 1, DEF_DA_ENS_NUM
            CALL load_member(i)
            IF (p_is_worker) CALL CoLMDRIVER (idate, deltim, dolai, doalb, dosst, oro)
            CALL save_member(i)
         ENDDO
      ENDIF

      ! Run data assimilation
      IF (.not. is_spinup .and. DEF_DA_ENS_NUM > 1) THEN
         CALL mpi_barrier (p_comm_glb, p_err)
         IF (use_assim_2mt .or. use_assim_2mq) CALL run_DA_TQ (idate, deltim)
         IF (use_assim_sm) CALL run_DA_SM (idate, deltim)
         IF (use_assim_st) CALL run_DA_ST (idate, deltim)
#ifdef GridRiverLakeFlow
         IF (use_assim_disc) CALL run_DA_DISC (idate, deltim)
#endif
      ENDIF

      ! Restore open-loop (if offline) or analysis (if online) to model states
      CALL postprocess()

   END SUBROUTINE run_DA

!-----------------------------------------------------------------------------

   SUBROUTINE end_DA ()

!-----------------------------------------------------------------------------
   IMPLICIT NONE

!-----------------------------------------------------------------------------

      IF (use_assim_2mt .or. use_assim_2mq) CALL end_DA_TQ ()
      use_assim_2mt  = .false.
      use_assim_2mq  = .false.
      use_assim_disc = .false.
      IF (use_assim_sm) CALL end_DA_SM ()
      use_assim_sm   = .false.
      IF (use_assim_st) CALL end_DA_ST ()
      use_assim_st   = .false.

   END SUBROUTINE end_DA

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Main
#endif