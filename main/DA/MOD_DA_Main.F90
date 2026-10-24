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
   USE MOD_Vars_1DFluxes, only: rnof
   USE MOD_DA_Vars_1DFluxes
   USE MOD_DA_Log
#ifdef GridRiverLakeFlow
   USE MOD_Grid_RiverLakeFlow
   USE MOD_Grid_RiverLakeNetwork
   USE MOD_Grid_Reservoir
   USE MOD_Grid_RiverLakeHist
   USE MOD_DA_Vars_TimeVariables
   USE MOD_DA_Vars_Hist
#endif
   IMPLICIT NONE
   SAVE

   logical :: use_assim_sm = .false.
   logical :: use_assim_st = .false.

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE init_DA (idate)

!-----------------------------------------------------------------------------
   IMPLICIT NONE

   integer, intent(in) :: idate(3)

   integer :: i, n
   character(len=256) :: msg

!-----------------------------------------------------------------------------

      CALL DA_log_set_date (idate)
      CALL DA_log_key ('Initialize data assimilation.')
      CALL DA_log_key ('DA log file: '//trim(DEF_DA_LOG_FILE))

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

      IF (use_assim_sm) CALL DA_log_key ('Start surface soil moisture (SM) data assimilation.')
      IF (use_assim_st) CALL DA_log_key ('Start soil temperature (ST) data assimilation.')

      n = 0
      CALL DA_log_key ('Assimilation observation sources:')
      CALL DA_log_key ('+------+----------+----------+------------------+----------+------------+')
      CALL DA_log_key ('| No.  | Target   | Source   | Sensor           | Var      | Freq(GHz)  |')
      CALL DA_log_key ('+------+----------+----------+------------------+----------+------------+')

      DO i = 1, size(DEF_DA_OBS_TARGET)
         IF (len_trim(DEF_DA_OBS_TARGET(i)) > 0) THEN
            n = n + 1
            write (msg,'(A,I4,A,A8,A,A8,A,A16,A,A8,A,F10.3,A)') '| ', n, ' | ', &
               DEF_DA_OBS_TARGET(i), ' | ', &
               DEF_DA_OBS_SOURCE(i), ' | ', &
               DEF_DA_OBS_SENSOR(i), ' | ', &
               DEF_DA_OBS_VAR(i),    ' | ', &
               DEF_DA_OBS_FGHZ(i),         ' |'
            CALL DA_log_key (msg)
         ENDIF
      ENDDO

      IF (n == 0) THEN
         CALL DA_log_key ('| none |          |          |                  |          |            |')
      ENDIF

      CALL DA_log_key ('+------+----------+----------+------------------+----------+------------+')

      IF (use_assim_sm) CALL init_assim_SM ()
      IF (use_assim_st) CALL init_assim_ST ()

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
      CALL DA_log_set_date (idate)

      ! generate ensemble members
      CALL model_to_ens(0)
      IF (DEF_DA_ENS_NUM > 1) THEN
         CALL ensemble (deltim)
         DO i = 1, DEF_DA_ENS_NUM
             CALL ens_to_model(i)
             IF (p_is_worker) THEN
               CALL CoLMDRIVER (idate, deltim, dolai, doalb, dosst, oro)
             ENDIF
             CALL model_to_ens(i)
         ENDDO
#if (defined GridRiverLakeFlow)
         CALL run_DA_riverlake (idate, deltim, is_spinup)
#endif
      ENDIF

      CALL ens_to_model(0)

      IF (is_spinup) RETURN

      ! data assimilation
      CALL mpi_barrier (p_comm_glb, p_err)
      IF (use_assim_sm) THEN
         CALL run_assim_SM (idate, deltim)
      ENDIF
      IF (use_assim_st) THEN
         CALL run_assim_ST (idate, deltim)
      ENDIF

   END SUBROUTINE run_DA

!-----------------------------------------------------------------------------

#if (defined GridRiverLakeFlow)
   SUBROUTINE run_DA_riverlake (idate, deltim, is_spinup)

!-----------------------------------------------------------------------------
   IMPLICIT NONE

   integer,  intent(in) :: idate(3)
   real(r8), intent(in) :: deltim
   logical,  intent(in) :: is_spinup

!-----------------------------------------------------------------------------
      IF (is_spinup) RETURN

      IF (p_is_worker) THEN
         rnof = sum(rnof_ens(1:DEF_DA_ENS_NUM,:), dim=1) / real(DEF_DA_ENS_NUM, r8)
         IF (allocated(wdsrf_ucat_a)) THEN
            WHERE (wdsrf_ucat_a /= spval) wdsrf_ucat = wdsrf_ucat_a
            WHERE (veloc_riv_a  /= spval) veloc_riv  = veloc_riv_a
         ENDIF
         IF (DEF_Reservoir_Method > 0) THEN
            IF (allocated(volresv_a)) THEN
               WHERE (volresv_a /= spval) volresv = volresv_a
            ENDIF
         ENDIF
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
            a_discharge_a  = a_discharge_a  + (a_discharge  - a_discharge_ol )
            acctime_ucat_a = acctime_ucat_a + (acctime_ucat - acctime_ucat_ol)
            a_wdsrf_ucat_a = a_wdsrf_ucat_a + (a_wdsrf_ucat - a_wdsrf_ucat_ol)
            a_veloc_riv_a  = a_veloc_riv_a  + (a_veloc_riv  - a_veloc_riv_ol )
            a_floodarea_a  = a_floodarea_a  + (a_floodarea  - a_floodarea_ol )

            a_discharge  = a_discharge_ol
            acctime_ucat = acctime_ucat_ol
            a_wdsrf_ucat = a_wdsrf_ucat_ol
            a_veloc_riv  = a_veloc_riv_ol
            a_floodarea  = a_floodarea_ol
         ENDIF
         IF (numresv > 0) THEN
            acctime_resv_a = acctime_resv_a + (acctime_resv - acctime_resv_ol)
            a_volresv_a    = a_volresv_a    + (a_volresv    - a_volresv_ol   )
            a_qresv_in_a   = a_qresv_in_a   + (a_qresv_in   - a_qresv_in_ol  )
            a_qresv_out_a  = a_qresv_out_a  + (a_qresv_out  - a_qresv_out_ol )

            acctime_resv = acctime_resv_ol
            a_volresv    = a_volresv_ol
            a_qresv_in   = a_qresv_in_ol
            a_qresv_out  = a_qresv_out_ol
         ENDIF
         IF (allocated(wdsrf_ucat_a)) THEN
            wdsrf_ucat_a = wdsrf_ucat
            veloc_riv_a  = veloc_riv
         ENDIF
         IF (DEF_Reservoir_Method > 0) THEN
            IF (allocated(volresv_a)) volresv_a = volresv
         ENDIF
      ENDIF

   END SUBROUTINE run_DA_riverlake
#endif

!-----------------------------------------------------------------------------

   SUBROUTINE end_DA ()

!-----------------------------------------------------------------------------
   IMPLICIT NONE

!-----------------------------------------------------------------------------

      IF (use_assim_sm) CALL end_assim_SM ()
      use_assim_sm = .false.
      IF (use_assim_st) CALL end_assim_ST ()
      use_assim_st = .false.

      CALL end_DA_log ()

   END SUBROUTINE end_DA

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Main
#endif
