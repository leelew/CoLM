#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Vars_Hist
!-----------------------------------------------------------------------------
! DESCRIPTION:
!   History output of DA-related quantities.
!
!   Write history diagnostics of model states and observation equivalents
!   computed through observation operators. Observation equivalents are
!   output only where observations are available, ensuring consistent
!   spatial coverage for comparison with observations.
!
! AUTHOR:
!   Zhilong Fan, 09/2026
!   Lu Li, Zhilong Fan, 11/2026: DA v1.0
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_Vars_Global
   USE MOD_NetCDFSerial
   USE MOD_LandPatch
   USE MOD_DataType
   USE MOD_HistGridded
   USE MOD_HistInterface
   USE MOD_Vars_TimeInvariants
   USE MOD_Forcing
#ifdef GridRiverLakeFlow
   USE MOD_Vector_ReadWrite
   USE MOD_Grid_RiverLakeNetwork
#endif
   USE MOD_DA_Vars_TimeVars
   USE MOD_DA_Vars_AccVars
   IMPLICIT NONE
   SAVE

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE save_DA_hist (file_hist, itime_in_file, idate)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      character(len=*),         intent(in) :: file_hist
      integer,                  intent(in) :: itime_in_file
      integer,                  intent(in) :: idate(3)

      logical, allocatable                 :: filter_model(:)        ! Valid patches for model-state diagnostics
      logical, allocatable                 :: filter_obs(:)          ! Valid patches with H(x) for the current observation stream
      type(block_data_real8_2d)            :: sumarea_model          ! Mapped area of valid model patches
      type(block_data_real8_2d)            :: sumarea_obs            ! Mapped area of valid H(x) patches for the current observation stream
      character(len=96)                    :: stream_suffix
      character(len=96)                    :: stream_label
      character(len=256)                   :: file_hist_ucat
      integer                              :: itime_in_file_ucat
      integer                              :: i, s

!-----------------------------------------------------------------------------
      IF (HistForm == 'Gridded' .and. p_is_io) CALL allocate_block_data (ghist, sumarea_model)
      IF (HistForm == 'Gridded' .and. p_is_io) CALL allocate_block_data (ghist, sumarea_obs  )
      IF (p_is_worker .and. numpatch > 0) THEN
         allocate (filter_model(numpatch))
         allocate (filter_obs  (numpatch))
      ENDIF

      ! Write model-space diagnostics for all valid patches.
      IF (p_is_worker .and. numpatch > 0) THEN
         filter_model(:) = patchtype <= 3 .and. patchmask
         IF (DEF_forcing%has_missing_value) filter_model = filter_model .and. forcmask_pch
      ENDIF
      IF (HistForm == 'Gridded') CALL mp2g_hist%get_sumarea (sumarea_model, filter_model)

      CALL write_history_variable_3d ( DEF_hist_vars%DA_wliq_soisno_ol,                                   &
         a_wliq_soisno_ol, file_hist, 'f_wliq_soisno_ol', itime_in_file, 'soilsnow',                      &
         maxsnl + 1, nl_soil - maxsnl, sumarea_model, filter_model,                                       &
         'Open-loop liquid water in soil and snow layers','kg/m2')
      CALL write_history_variable_3d ( DEF_hist_vars%DA_wliq_soisno_f,                                    &
         a_wliq_soisno_f, file_hist, 'f_wliq_soisno_f', itime_in_file, 'soilsnow',                        &
         maxsnl + 1, nl_soil - maxsnl, sumarea_model, filter_model,                                       &
         'Forecast ensemble-mean liquid water in soil and snow layers','kg/m2')
      CALL write_history_variable_3d ( DEF_hist_vars%DA_wliq_soisno_a,                                    &
         a_wliq_soisno_a, file_hist, 'f_wliq_soisno_a', itime_in_file, 'soilsnow',                        &
         maxsnl + 1, nl_soil - maxsnl, sumarea_model, filter_model,                                       &
         'Analysis ensemble-mean liquid water in soil and snow layers','kg/m2')
      CALL write_history_variable_3d ( DEF_hist_vars%DA_t_soisno_ol,                                      &
         a_t_soisno_ol, file_hist, 'f_t_soisno_ol', itime_in_file, 'soilsnow',                            &
         maxsnl + 1, nl_soil - maxsnl, sumarea_model, filter_model,                                       &
         'Open-loop temperature in soil and snow layers','K')
      CALL write_history_variable_3d ( DEF_hist_vars%DA_t_soisno_f,                                       &
         a_t_soisno_f, file_hist, 'f_t_soisno_f', itime_in_file, 'soilsnow',                              &
         maxsnl + 1, nl_soil - maxsnl, sumarea_model, filter_model,                                       &
         'Forecast ensemble-mean temperature in soil and snow layers','K')
      CALL write_history_variable_3d ( DEF_hist_vars%DA_t_soisno_a,                                       &
         a_t_soisno_a, file_hist, 'f_t_soisno_a', itime_in_file, 'soilsnow',                              &
         maxsnl + 1, nl_soil - maxsnl, sumarea_model, filter_model,                                       &
         'Analysis ensemble-mean temperature in soil and snow layers','K')
      CALL write_history_variable_4d ( DEF_hist_vars%DA_wliq_soisno_ens,                                  &
         a_wliq_soisno_ens, file_hist, 'f_wliq_soisno_ens', itime_in_file,                                &
         'soilsnow', maxsnl + 1, nl_soil - maxsnl, 'ens', 1, DEF_DA_ENS_NUM, sumarea_model, filter_model, &
         'Ensemble liquid water in soil layers','kg/m2')

      ! Write H(x) diagnostics for each observation stream, excluding ANATQ.
      s = 0
      DO i = 1, size(DEF_DA_OBS_TARGET)
         s = s + 1
         IF (trim(DEF_DA_OBS_SOURCE(i)) == 'ANATQ') CYCLE

         stream_suffix = '_'//trim(DEF_DA_OBS_SENSOR(i))//'_'//                                           &
            trim(DEF_DA_OBS_VAR(i))//'_'//trim(DEF_DA_OBS_TARGET(i))
         stream_label = trim(DEF_DA_OBS_SENSOR(i))//' '//                                                 &
            trim(DEF_DA_OBS_VAR(i))//' for '//trim(DEF_DA_OBS_TARGET(i))

         IF (p_is_worker .and. numpatch > 0) THEN
            filter_obs = filter_model .and. a_hx_ol(s,:) /= spval
         ENDIF
         IF (HistForm == 'Gridded') CALL mp2g_hist%get_sumarea (sumarea_obs, filter_obs)

         CALL write_history_variable_2d ( DEF_hist_vars%DA_hx_ol,                                         &
            a_hx_ol(s,:), file_hist, 'f_hx_ol'//trim(stream_suffix), itime_in_file,                       &
            sumarea_obs, filter_obs,                                                                      &
            'Open-loop H(x) of '//trim(stream_label), '', nac_obs(s,:))
         CALL write_history_variable_2d ( DEF_hist_vars%DA_hx_f,                                          &
            a_hx_f(s,:), file_hist, 'f_hx_f'//trim(stream_suffix), itime_in_file,                         &
            sumarea_obs, filter_obs,                                                                      &
            'Prior ensemble-mean H(x) of '//trim(stream_label), '', nac_obs(s,:))
         CALL write_history_variable_2d ( DEF_hist_vars%DA_hx_a,                                          &
            a_hx_a(s,:), file_hist, 'f_hx_a'//trim(stream_suffix), itime_in_file,                         &
            sumarea_obs, filter_obs,                                                                      &
            'Analysis ensemble-mean H(x) of '//trim(stream_label), '', nac_obs(s,:))
         CALL write_history_variable_3d ( DEF_hist_vars%DA_hx_f_ens,                                      &
            a_hx_f_ens(s,:,:), file_hist, 'f_hx_f_ens'//trim(stream_suffix), itime_in_file,               &
            'ens', 1, DEF_DA_ENS_NUM, sumarea_obs, filter_obs,                                            &
            'Prior ensemble H(x) of '//trim(stream_label), '', nac_obs(s,:))
         CALL write_history_variable_3d ( DEF_hist_vars%DA_hx_a_ens,                                      &
            a_hx_a_ens(s,:,:), file_hist, 'f_hx_a_ens'//trim(stream_suffix), itime_in_file,               &
            'ens', 1, DEF_DA_ENS_NUM, sumarea_obs, filter_obs,                                            &
            'Analysis ensemble H(x) of '//trim(stream_label), '', nac_obs(s,:))
      ENDDO

      ! Write discharge diagnostics.
#ifdef GridRiverLakeFlow
      IF (DEF_hist_vars%DA_discharge) THEN
         IF (p_is_master) THEN
            i = index(trim(file_hist), '_', back=.true.)
            file_hist_ucat = file_hist(1:i) // 'unitcat_' // file_hist(i+1:)
            CALL ncio_write_time (file_hist_ucat, 'time', idate, itime_in_file_ucat, DEF_HIST_FREQ)
         ENDIF

         IF (p_is_worker .and. numucat > 0) THEN
            WHERE (acctime_ucat_a > 0)
               a_discharge_a = a_discharge_a / acctime_ucat_a
            ELSE WHERE
               a_discharge_a = spval
            END WHERE
         ENDIF
         CALL vector_gather_map2grid_and_write ( a_discharge_a, numucat,                                  &
            totalnumucat, ucat_data_address, griducat%nlon, x_ucat, griducat%nlat, y_ucat,                &
            file_hist_ucat, 'f_discharge_a', 'lon_ucat', 'lat_ucat', itime_in_file_ucat,                  &
            'Discharge in river and flood plain from the ensemble-mean state', 'm^3/s')
      ENDIF
#endif

      IF (allocated(filter_model  )) deallocate (filter_model)
      IF (allocated(filter_obs    )) deallocate (filter_obs  )
      IF (allocated(a_discharge_a )) a_discharge_a  = 0.0_r8
      IF (allocated(acctime_ucat_a)) acctime_ucat_a = 0.0_r8

   END SUBROUTINE save_DA_hist

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Vars_Hist
#endif