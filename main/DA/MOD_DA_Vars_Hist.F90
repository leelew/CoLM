#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Vars_Hist
!-----------------------------------------------------------------------------
! DESCRIPTION:
!   History output of DA-related quantities.
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

   type(block_data_real8_2d) :: sumarea_hx  ! mapped area of each H(x) stream

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE save_DA_hist (file_hist, itime_in_file, idate, sumarea, filter)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      character(len=*),         intent(in)    :: file_hist
      integer,                  intent(in)    :: itime_in_file
      integer,                  intent(in)    :: idate(3)
      type(block_data_real8_2d), intent(inout) :: sumarea
      logical,                  intent(inout) :: filter(:)

!------------------------ Local Variables ------------------------------------
      logical, allocatable :: filter_s(:)
      character(len=96) :: stream_suffix, stream_label
#ifdef GridRiverLakeFlow
      character(len=256) :: file_hist_ucat
      integer :: itime_in_file_ucat
#endif
      integer :: i, s

!-----------------------------------------------------------------------------

      ! Build the valid-patch mask and mapped area.
      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            filter(:) = patchtype <= 3
            IF (DEF_forcing%has_missing_value) filter = filter .and. forcmask_pch
            filter = filter .and. patchmask
         ENDIF
      ENDIF
      IF (HistForm == 'Gridded') CALL mp2g_hist%get_sumarea (sumarea, filter)

      IF (DEF_DA) THEN
         ! Write soil and snow diagnostics.
         CALL write_history_variable_3d ( DEF_hist_vars%DA_wliq_soisno_ol, &
            a_wliq_soisno_ol, file_hist, 'f_wliq_soisno_ol', itime_in_file, 'soilsnow', &
            maxsnl + 1, nl_soil - maxsnl, sumarea, filter, 'Open-loop liquid water in soil and snow layers','kg/m2')
         CALL write_history_variable_3d ( DEF_hist_vars%DA_wliq_soisno_f, &
            a_wliq_soisno_f, file_hist, 'f_wliq_soisno_f', itime_in_file, 'soilsnow', &
            maxsnl + 1, nl_soil - maxsnl, sumarea, filter, 'Forecast ensemble-mean liquid water in soil and snow layers','kg/m2')
         CALL write_history_variable_3d ( DEF_hist_vars%DA_wliq_soisno_a, &
            a_wliq_soisno_a, file_hist, 'f_wliq_soisno_a', itime_in_file, 'soilsnow', &
            maxsnl + 1, nl_soil - maxsnl, sumarea, filter, 'Analysis ensemble-mean liquid water in soil and snow layers','kg/m2')
         CALL write_history_variable_3d ( DEF_hist_vars%DA_t_soisno_ol, &
            a_t_soisno_ol, file_hist, 'f_t_soisno_ol', itime_in_file, 'soilsnow', &
            maxsnl + 1, nl_soil - maxsnl, sumarea, filter, 'Open-loop temperature in soil and snow layers','K')
         CALL write_history_variable_3d ( DEF_hist_vars%DA_t_soisno_f, &
            a_t_soisno_f, file_hist, 'f_t_soisno_f', itime_in_file, 'soilsnow', &
            maxsnl + 1, nl_soil - maxsnl, sumarea, filter, 'Forecast ensemble-mean temperature in soil and snow layers','K')
         CALL write_history_variable_3d ( DEF_hist_vars%DA_t_soisno_a, &
            a_t_soisno_a, file_hist, 'f_t_soisno_a', itime_in_file, 'soilsnow', &
            maxsnl + 1, nl_soil - maxsnl, sumarea, filter, 'Analysis ensemble-mean temperature in soil and snow layers','K')
         CALL write_history_variable_4d ( DEF_hist_vars%DA_wliq_soisno_ens, &
            a_wliq_soisno_ens, file_hist, 'f_wliq_soisno_ens', itime_in_file, &
            'soilsnow', maxsnl + 1, nl_soil - maxsnl, 'ens', 1, DEF_DA_ENS_NUM, sumarea, filter, 'ensemble liquid water in soil layers','kg/m2')
      ENDIF

      ! Write H(x) diagnostics, one variable per observation stream.
      ! ANATQ slots are deliberately left as spval, so they are skipped.
      ! Streams are assumed to be contiguous namelist entries from 1, which
      ! is also how the source modules fill the hx slots.
      IF (nsource > 0) THEN
         IF (p_is_worker) THEN
            allocate (filter_s(numpatch))
         ELSE
            allocate (filter_s(0))
         ENDIF
         IF (HistForm == 'Gridded' .and. p_is_io) THEN
            IF (.not. allocated(sumarea_hx%blk)) CALL allocate_block_data (ghist, sumarea_hx)
         ENDIF

         s = 0
         DO i = 1, size(DEF_DA_OBS_TARGET)
            IF (len_trim(DEF_DA_OBS_TARGET(i)) == 0) CYCLE
            s = s + 1
            IF (trim(DEF_DA_OBS_SOURCE(i)) == 'ANATQ') CYCLE

            stream_suffix = '_'//trim(DEF_DA_OBS_SENSOR(i))//'_'// &
               trim(DEF_DA_OBS_VAR(i))//'_'//trim(DEF_DA_OBS_TARGET(i))
            stream_label = trim(DEF_DA_OBS_SENSOR(i))//' '// &
               trim(DEF_DA_OBS_VAR(i))//' for '//trim(DEF_DA_OBS_TARGET(i))

            ! The three H(x) means of a stream share the open-loop validity mask.
            IF (p_is_worker) filter_s = filter .and. a_hx_ol(s,:) /= spval
            IF (HistForm == 'Gridded') CALL mp2g_hist%get_sumarea (sumarea_hx, filter_s)

            CALL write_history_variable_2d ( DEF_hist_vars%DA_hx_ol, &
               a_hx_ol(s,:), file_hist, 'f_hx_ol'//trim(stream_suffix), itime_in_file, &
               sumarea_hx, filter_s, 'Open-loop H(x) of '//trim(stream_label), 'source-dependent', nac_hx_ol(s,:))
            IF (.not. DEF_DA) CYCLE
            CALL write_history_variable_2d ( DEF_hist_vars%DA_hx_f, &
               a_hx_f(s,:), file_hist, 'f_hx_f'//trim(stream_suffix), itime_in_file, &
               sumarea_hx, filter_s, 'Prior H(x) of '//trim(stream_label), 'source-dependent')
            CALL write_history_variable_2d ( DEF_hist_vars%DA_hx_a, &
               a_hx_a(s,:), file_hist, 'f_hx_a'//trim(stream_suffix), itime_in_file, &
               sumarea_hx, filter_s, 'Analysis H(x) of '//trim(stream_label), 'source-dependent')
            CALL write_history_variable_3d ( DEF_hist_vars%DA_hx_f_ens, &
               a_hx_f_ens(s,:,:), file_hist, 'f_hx_f_ens'//trim(stream_suffix), itime_in_file, &
               'ens', 1, DEF_DA_ENS_NUM, sumarea, filter, 'Prior ensemble H(x) of '//trim(stream_label), 'source-dependent')
            CALL write_history_variable_3d ( DEF_hist_vars%DA_hx_a_ens, &
               a_hx_a_ens(s,:,:), file_hist, 'f_hx_a_ens'//trim(stream_suffix), itime_in_file, &
               'ens', 1, DEF_DA_ENS_NUM, sumarea, filter, 'Analysis ensemble H(x) of '//trim(stream_label), 'source-dependent')
         ENDDO

         deallocate (filter_s)
      ENDIF

      IF (.not. DEF_DA) RETURN

      ! Normalize and write the analysis discharge.
#ifdef GridRiverLakeFlow
      itime_in_file_ucat = 1
      IF (p_is_master) THEN
         i = len_trim (file_hist)
         DO WHILE (file_hist(i:i) /= '_')
            i = i - 1
         ENDDO
         file_hist_ucat = file_hist(1:i) // 'unitcat_' // file_hist(i+1:)
         CALL ncio_write_time (file_hist_ucat, 'time', idate, itime_in_file_ucat, DEF_HIST_FREQ)
      ENDIF

      IF (DEF_hist_vars%DA_discharge) THEN
         IF (p_is_worker) THEN
            IF (numucat > 0) THEN
               WHERE (acctime_ucat_a > 0)
                  a_discharge_a = a_discharge_a / acctime_ucat_a
               ELSE WHERE
                  a_discharge_a = spval
               END WHERE
            ENDIF
         ENDIF
         CALL vector_gather_map2grid_and_write ( a_discharge_a, numucat,       &
            totalnumucat, ucat_data_address, griducat%nlon, x_ucat, griducat%nlat, y_ucat, &
            file_hist_ucat, 'f_discharge_a', 'lon_ucat', 'lat_ucat', itime_in_file_ucat, &
            'discharge in river and flood plain from the ensemble-mean state', 'm^3/s')
      ENDIF

      IF (p_is_worker) THEN
         IF (numucat > 0) THEN
            IF (allocated(a_discharge_a )) a_discharge_a  (:) = 0.0_r8
            IF (allocated(acctime_ucat_a)) acctime_ucat_a (:) = 0.0_r8
         ENDIF
      ENDIF
#endif

   END SUBROUTINE save_DA_hist

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Vars_Hist
#endif
