#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Main
!-----------------------------------------------------------------------------
! DESCRIPTION:
!     Main procedures for data assimilation
!
! AUTHOR:
!     Lu Li, 12/2024
!     Zhilong Fan, Lu Li, 03/2024
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Spmd_task
   USE MOD_Namelist
   USE MOD_LandPatch
   USE MOD_DA_Assim_SM
   USE MOD_DA_Assim_ST
   USE MOD_DA_Ens
#ifdef GridRiverLakeFlow
   USE MOD_Grid_RiverLakeFlow
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   USE MOD_Grid_Reservoir, only: numresv
   USE MOD_Grid_RiverLakeHist, only: a_discharge, acctime_ucat, a_wdsrf_ucat, a_veloc_riv, a_floodarea, &
      a_volresv, a_qresv_in, a_qresv_out, acctime_resv
   USE MOD_DA_Vars_TimeVariables, only: a_discharge_ol, acctime_ucat_ol, a_wdsrf_ucat_ol, a_veloc_riv_ol, &
      a_floodarea_ol, a_discharge_mean, acctime_ucat_mean, &
      a_wdsrf_ucat_mean, a_veloc_riv_mean, a_floodarea_mean, &
      a_volresv_ol, a_qresv_in_ol, a_qresv_out_ol, acctime_resv_ol, &
      a_volresv_mean, a_qresv_in_mean, a_qresv_out_mean, acctime_resv_mean
#endif
   IMPLICIT NONE
   SAVE

   logical :: use_assim_sm = .false.
   logical :: use_assim_st = .false.

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE init_DA ()

!-----------------------------------------------------------------------------
   IMPLICIT NONE

   integer :: i

!-----------------------------------------------------------------------------

      use_assim_sm = .false.
      DO i = 1, size(DEF_DA_OBS_TARGET)
         IF (trim(DEF_DA_OBS_TARGET(i)) == 'SM') THEN
            use_assim_sm = .true.
            EXIT
         ENDIF
      ENDDO
      use_assim_st = .false.
      DO i = 1, size(DEF_DA_OBS_TARGET)
         IF (trim(DEF_DA_OBS_TARGET(i)) == 'ST') THEN
            use_assim_st = .true.
            EXIT
         ENDIF
      ENDDO

      IF (use_assim_sm) CALL init_Assim_SM ()
      IF (use_assim_st) CALL init_Assim_ST ()

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

   integer :: i, n

!-----------------------------------------------------------------------------
      ! generate ensemble members
      CALL model_to_ens_member(0)
      IF (DEF_DA_ENS_NUM > 1) THEN
         CALL ensemble (deltim)
         DO i = 1, DEF_DA_ENS_NUM
             CALL ens_member_to_model(i)
             IF (p_is_worker) THEN
               CALL CoLMDRIVER (idate, deltim, dolai, doalb, dosst, oro)
             ENDIF
             CALL model_to_ens_member(i)
         ENDDO
         CALL ens_mean_to_model()
#if (defined GridRiverLakeFlow)
         IF (.not. is_spinup) THEN
             IF (p_is_worker) THEN
                 IF (numucat > 0) THEN
                     a_discharge_ol  = a_discharge
                     acctime_ucat_ol = acctime_ucat
                     a_wdsrf_ucat_ol = a_wdsrf_ucat
                     a_veloc_riv_ol  = a_veloc_riv
                     a_floodarea_ol  = a_floodarea
                 ENDIF
                 IF (numresv > 0) THEN
                     acctime_resv_ol = acctime_resv
                     a_volresv_ol    = a_volresv
                     a_qresv_in_ol   = a_qresv_in
                     a_qresv_out_ol  = a_qresv_out
                 ENDIF
             ENDIF

             CALL grid_riverlake_flow (idate(1), deltim)

             IF (p_is_worker) THEN
                 IF (numucat > 0) THEN
                     a_discharge_mean  = a_discharge_mean  + (a_discharge  - a_discharge_ol)
                     acctime_ucat_mean = acctime_ucat_mean + (acctime_ucat - acctime_ucat_ol)
                     a_wdsrf_ucat_mean = a_wdsrf_ucat_mean + (a_wdsrf_ucat - a_wdsrf_ucat_ol)
                     a_veloc_riv_mean  = a_veloc_riv_mean  + (a_veloc_riv  - a_veloc_riv_ol)
                     a_floodarea_mean  = a_floodarea_mean  + (a_floodarea  - a_floodarea_ol)

                     a_discharge  = a_discharge_ol
                     acctime_ucat = acctime_ucat_ol
                     a_wdsrf_ucat = a_wdsrf_ucat_ol
                     a_veloc_riv  = a_veloc_riv_ol
                     a_floodarea  = a_floodarea_ol
                 ENDIF
                 IF (numresv > 0) THEN
                     acctime_resv_mean = acctime_resv_mean + (acctime_resv - acctime_resv_ol)
                     a_volresv_mean    = a_volresv_mean    + (a_volresv    - a_volresv_ol)
                     a_qresv_in_mean   = a_qresv_in_mean   + (a_qresv_in   - a_qresv_in_ol)
                     a_qresv_out_mean  = a_qresv_out_mean  + (a_qresv_out  - a_qresv_out_ol)

                     acctime_resv = acctime_resv_ol
                     a_volresv    = a_volresv_ol
                     a_qresv_in   = a_qresv_in_ol
                     a_qresv_out  = a_qresv_out_ol
                 ENDIF
             ENDIF
         ENDIF
#endif
         CALL model_to_ens_mean()
      ELSE
         ! open-loop only: fill the single ensemble slot with the model
         ! state so 1:N diagnostics and means stay valid (= open-loop).
         CALL model_to_ens_member(1)
      ENDIF

      CALL ens_member_to_model(0)

      IF (is_spinup) RETURN

      ! data assimilation
      IF (p_is_master) THEN
         IF (use_assim_sm) write(*,'(A)') '[CoLM-DA] Start surface soil moisture (SM) data assimilation.'
         IF (use_assim_st) write(*,'(A)') '[CoLM-DA] Start soil temperature (ST) data assimilation.'

         n = 0
         write(*,'(A)') '[CoLM-DA] Assimilation observation sources:'
         write(*,'(A)') '[CoLM-DA] +------+----------+----------+------------------+----------+------------+'
         write(*,'(A)') '[CoLM-DA] | No.  | Target   | Source   | Sensor           | Var      | Freq(GHz)  |'
         write(*,'(A)') '[CoLM-DA] +------+----------+----------+------------------+----------+------------+'

         DO i = 1, size(DEF_DA_OBS_TARGET)
            IF (len_trim(DEF_DA_OBS_TARGET(i)) > 0) THEN
               n = n + 1
               write(*,'(A,I4,A,A8,A,A8,A,A16,A,A8,A,F10.3,A)') &
                  '[CoLM-DA] | ', n, ' | ', &
                  DEF_DA_OBS_TARGET(i), ' | ', &
                  DEF_DA_OBS_SOURCE(i), ' | ', &
                  DEF_DA_OBS_SENSOR(i), ' | ', &
                  DEF_DA_OBS_VAR(i),    ' | ', &
                  DEF_DA_OBS_FGHZ(i),         ' |'
            ENDIF
         ENDDO

         IF (n == 0) THEN
            write(*,'(A)') '[CoLM-DA] | none |          |          |                  |          |            |'
         ENDIF

         write(*,'(A)') '[CoLM-DA] +------+----------+----------+------------------+----------+------------+'
      ENDIF


      CALL mpi_barrier (p_comm_glb, p_err)
      IF (use_assim_sm) THEN
         CALL run_Assim_SM (idate, deltim)
      ENDIF
      IF (use_assim_st) THEN
         CALL run_Assim_ST (idate, deltim)
      ENDIF

   END SUBROUTINE run_DA

!-----------------------------------------------------------------------------

   SUBROUTINE end_DA ()

!-----------------------------------------------------------------------------
   IMPLICIT NONE

!-----------------------------------------------------------------------------

      IF (use_assim_sm) CALL end_Assim_SM ()
      use_assim_sm = .false.
      IF (use_assim_st) CALL end_Assim_ST ()
      use_assim_st = .false.

   END SUBROUTINE end_DA

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Main
#endif
