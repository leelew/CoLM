#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Vars_Hist
!-----------------------------------------------------------------------------
! DESCRIPTION:
!   History output of DA-related quantities.
!
! AUTHOR:
!   Zhilong Fan, 09/2026
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_Vars_Global
   USE MOD_NetCDFSerial
   USE MOD_LandPatch
   USE MOD_DataType
   USE MOD_HistGridded
   USE MOD_Vars_TimeInvariants
   USE MOD_Forcing
#ifdef GridRiverLakeFlow
   USE MOD_Vector_ReadWrite
   USE MOD_Grid_RiverLakeNetwork
#endif
   USE MOD_DA_Vars_TimeVars
   USE MOD_DA_Vars_AccVars
   USE MOD_DA_Vars_Utils
   IMPLICIT NONE
   SAVE

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE save_DA_hist (file_hist, HistForm, itime_in_file, idate, sumarea, filter)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      character(len=*),             intent(in) :: file_hist
      character(len=*),             intent(in) :: HistForm
      integer,                      intent(in) :: itime_in_file
      integer,                      intent(in) :: idate(3)
      type(block_data_real8_2d), intent(inout) :: sumarea
      logical,                   intent(inout) :: filter(:)

      logical, allocatable :: filter_hx(:,:)
#ifdef GridRiverLakeFlow
      character(len=256) :: file_hist_ucat
      integer :: itime_in_file_ucat
#endif
      integer :: i

!-----------------------------------------------------------------------------

      ! Set filter
      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            filter(:) = patchtype <= 3
            IF (DEF_forcing%has_missing_value) filter = filter .and. forcmask_pch
            filter = filter .and. patchmask
         ENDIF
      ENDIF
      IF (HistForm == 'Gridded') CALL mp2g_hist%get_sumarea (sumarea, filter)

      ! Save physics variables
      CALL write_DA_hist (HistForm, DEF_hist_vars%DA_wliq_soisno_ol, &
         a_wliq_soisno_ol, file_hist, 'f_wliq_soisno_ol', itime_in_file, &
         'soilsnow', maxsnl + 1, nl_soil - maxsnl, sumarea, filter, &
         'Open-loop liquid water in soil and snow layers', 'kg/m2')
      CALL write_DA_hist (HistForm, DEF_hist_vars%DA_wliq_soisno_f, &
         a_wliq_soisno_f, file_hist, 'f_wliq_soisno_f', itime_in_file, &
         'soilsnow', maxsnl + 1, nl_soil - maxsnl, sumarea, filter, &
         'Forecast ensemble-mean liquid water in soil and snow layers', 'kg/m2')
      CALL write_DA_hist (HistForm, DEF_hist_vars%DA_wliq_soisno_a, &
         a_wliq_soisno_a, file_hist, 'f_wliq_soisno_a', itime_in_file, &
         'soilsnow', maxsnl + 1, nl_soil - maxsnl, sumarea, filter, &
         'Analysis ensemble-mean liquid water in soil and snow layers', 'kg/m2')
      CALL write_DA_hist (HistForm, DEF_hist_vars%DA_t_soisno_ol, &
         a_t_soisno_ol, file_hist, 'f_t_soisno_ol', itime_in_file, &
         'soilsnow', maxsnl + 1, nl_soil - maxsnl, sumarea, filter, &
         'Open-loop temperature in soil and snow layers', 'K')
      CALL write_DA_hist (HistForm, DEF_hist_vars%DA_t_soisno_f, &
         a_t_soisno_f, file_hist, 'f_t_soisno_f', itime_in_file, &
         'soilsnow', maxsnl + 1, nl_soil - maxsnl, sumarea, filter, &
         'Forecast ensemble-mean temperature in soil and snow layers', 'K')
      CALL write_DA_hist (HistForm, DEF_hist_vars%DA_t_soisno_a, &
         a_t_soisno_a, file_hist, 'f_t_soisno_a', itime_in_file, &
         'soilsnow', maxsnl + 1, nl_soil - maxsnl, sumarea, filter, &
         'Analysis ensemble-mean temperature in soil and snow layers', 'K')
      CALL write_DA_hist (HistForm, DEF_hist_vars%DA_wliq_soisno_ens, &
         a_wliq_soisno_ens, file_hist, 'f_wliq_soisno_ens', itime_in_file, &
         'soilsnow', maxsnl + 1, nl_soil - maxsnl, 'ens', 1, DEF_DA_ENS_NUM, &
         sumarea, filter, 'ensemble liquid water in soil layers', 'kg/m2')

      ! Save h(x)
      IF (nsource > 0) THEN
         allocate (filter_hx(nsource,numpatch))

         IF (p_is_worker) THEN
            DO i = 1, nsource
               filter_hx(i,:) = filter .and. a_hx_ol(i,:) /= spval
            ENDDO
         ENDIF
         CALL write_DA_hist (HistForm, DEF_hist_vars%DA_hx_ol, &
            a_hx_ol, file_hist, 'f_hx_ol', itime_in_file, &
            'source', 1, nsource, sumarea, filter, &
            'Open-loop H(x) by DA observation source', 'source-dependent', &
            filter_3d = filter_hx)
         CALL write_DA_hist (HistForm, DEF_hist_vars%DA_hx_f, &
            a_hx_f, file_hist, 'f_hx_f', itime_in_file, &
            'source', 1, nsource, sumarea, filter, &
            'Prior H(x) by DA observation source', 'source-dependent', &
            filter_3d = filter_hx)
         CALL write_DA_hist (HistForm, DEF_hist_vars%DA_hx_a, &
            a_hx_a, file_hist, 'f_hx_a', itime_in_file, &
            'source', 1, nsource, sumarea, filter, &
            'Analysis H(x) by DA observation source', 'source-dependent', &
            filter_3d = filter_hx)
         CALL write_DA_hist (HistForm, DEF_hist_vars%DA_hx_f_ens, &
            a_hx_f_ens, file_hist, 'f_hx_f_ens', itime_in_file, &
            'source', 1, nsource, 'ens', 1, DEF_DA_ENS_NUM, &
            sumarea, filter, 'Prior ensemble H(x) by DA observation source', &
            'source-dependent')
         CALL write_DA_hist (HistForm, DEF_hist_vars%DA_hx_a_ens, &
            a_hx_a_ens, file_hist, 'f_hx_a_ens', itime_in_file, &
            'source', 1, nsource, 'ens', 1, DEF_DA_ENS_NUM, &
            sumarea, filter, 'Analysis ensemble H(x) by DA observation source', &
            'source-dependent')

         deallocate (filter_hx)
      ENDIF

      ! Save discharge
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
            a_discharge_a  (:) = 0.
            acctime_ucat_a (:) = 0.
         ENDIF
      ENDIF
#endif

   END SUBROUTINE save_DA_hist

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Vars_Hist
#endif
