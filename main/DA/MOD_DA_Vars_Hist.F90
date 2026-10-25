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
   USE MOD_LandPatch, only: numpatch
   USE MOD_DataType
   USE MOD_HistGridded
   USE MOD_Vars_TimeInvariants
   USE MOD_Forcing
   USE MOD_DA_Vars_TimeVariables
   USE MOD_DA_Vars_1DAccFluxes
#ifdef GridRiverLakeFlow
   USE MOD_Vector_ReadWrite
   USE MOD_Grid_RiverLakeNetwork
#endif
   IMPLICIT NONE
   SAVE

   PUBLIC :: hist_DA_history_out

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE hist_DA_history_out (file_hist, HistForm, itime_in_file, idate, &
         sumarea, filter)

!-----------------------------------------------------------------------------
   IMPLICIT NONE

   character(len=*), intent(in) :: file_hist
   character(len=*), intent(in) :: HistForm
   integer,          intent(in) :: itime_in_file
   integer,          intent(in) :: idate(3)

   type(block_data_real8_2d), intent(inout) :: sumarea
   logical,                   intent(inout) :: filter(:)

   logical, allocatable :: filter_hx(:,:)
   integer :: i
#ifdef GridRiverLakeFlow
   character(len=256) :: file_hist_ucat
   integer :: itime_in_file_ucat
#endif

!-----------------------------------------------------------------------------
      ! --------------------------------------------------------------------
      ! Soil and snow state diagnostics, same filter as the non-DA
      ! temperature and water variables.
      ! [soil => 0; urban and built-up => 1; wetland => 2; land ice => 3;
      !  land water bodies => 4; ocean => 99]
      ! --------------------------------------------------------------------

      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            filter(:) = patchtype <= 3

            IF (DEF_forcing%has_missing_value) THEN
               filter = filter .and. forcmask_pch
            ENDIF

            filter = filter .and. patchmask
         ENDIF
      ENDIF

      IF (HistForm == 'Gridded') THEN
         CALL mp2g_hist%get_sumarea (sumarea, filter)
      ENDIF

      CALL write_DA_history_variable_3d (HistForm, DEF_hist_vars%DA_wliq_soisno_ol, &
         a_wliq_soisno_ol, file_hist, 'f_wliq_soisno_ol', itime_in_file, &
         'soilsnow', maxsnl + 1, nl_soil - maxsnl, sumarea, filter, &
         'Open-loop liquid water in soil and snow layers', 'kg/m2')
      CALL write_DA_history_variable_3d (HistForm, DEF_hist_vars%DA_wliq_soisno_f, &
         a_wliq_soisno_f, file_hist, 'f_wliq_soisno_f', itime_in_file, &
         'soilsnow', maxsnl + 1, nl_soil - maxsnl, sumarea, filter, &
         'Forecast ensemble-mean liquid water in soil and snow layers', 'kg/m2')
      CALL write_DA_history_variable_3d (HistForm, DEF_hist_vars%DA_wliq_soisno_a, &
         a_wliq_soisno_a, file_hist, 'f_wliq_soisno_a', itime_in_file, &
         'soilsnow', maxsnl + 1, nl_soil - maxsnl, sumarea, filter, &
         'Analysis ensemble-mean liquid water in soil and snow layers', 'kg/m2')

      CALL write_DA_history_variable_3d (HistForm, DEF_hist_vars%DA_t_soisno_ol, &
         a_t_soisno_ol, file_hist, 'f_t_soisno_ol', itime_in_file, &
         'soilsnow', maxsnl + 1, nl_soil - maxsnl, sumarea, filter, &
         'Open-loop temperature in soil and snow layers', 'K')
      CALL write_DA_history_variable_3d (HistForm, DEF_hist_vars%DA_t_soisno_f, &
         a_t_soisno_f, file_hist, 'f_t_soisno_f', itime_in_file, &
         'soilsnow', maxsnl + 1, nl_soil - maxsnl, sumarea, filter, &
         'Forecast ensemble-mean temperature in soil and snow layers', 'K')
      CALL write_DA_history_variable_3d (HistForm, DEF_hist_vars%DA_t_soisno_a, &
         a_t_soisno_a, file_hist, 'f_t_soisno_a', itime_in_file, &
         'soilsnow', maxsnl + 1, nl_soil - maxsnl, sumarea, filter, &
         'Analysis ensemble-mean temperature in soil and snow layers', 'K')

      ! ensemble soil moisture & temperature in soil layers [kg/m2]
      CALL write_DA_history_variable_4d (HistForm, DEF_hist_vars%DA_wliq_soisno_ens, &
         a_wliq_soisno_ens, file_hist, 'f_wliq_soisno_ens', itime_in_file, &
         'soilsnow', maxsnl + 1, nl_soil - maxsnl, 'ens', 1, DEF_DA_ENS_NUM, &
         sumarea, filter, 'ensemble liquid water in soil layers', 'kg/m2')

      ! --------------------------------------------------------------------
      ! brightness temperature (excluding land ice, land water bodies and
      ! ocean patches)
      ! [soil => 0; urban and built-up => 1; wetland => 2; land ice => 3;
      !  land water bodies => 4; ocean => 99]
      ! --------------------------------------------------------------------

      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN

            filter(:) = patchtype <= 2

            IF (DEF_forcing%has_missing_value) THEN
               filter = filter .and. forcmask_pch
            END IF

            filter = filter .and. patchmask
         END IF
      END IF

      IF (HistForm == 'Gridded') THEN
         CALL mp2g_hist%get_sumarea (sumarea, filter)
      END IF

      ! Generic DA observation-space diagnostics. The source dimension follows
      ! the order of the non-empty DEF_DA_OBS_* entries.
      IF (nsource > 0) THEN
         allocate (filter_hx(nsource,numpatch))

         IF (p_is_worker) THEN
            DO i = 1, nsource
               filter_hx(i,:) = filter .and. a_hx_ol(i,:) /= spval
            ENDDO
         ENDIF
         CALL write_DA_history_variable_3d (HistForm, DEF_hist_vars%DA_hx_ol, &
            a_hx_ol, file_hist, 'f_hx_ol', itime_in_file, &
            'source', 1, nsource, sumarea, filter, &
            'Open-loop H(x) by DA observation source', 'source-dependent', &
            filter_3d = filter_hx)

         IF (p_is_worker) THEN
            DO i = 1, nsource
               filter_hx(i,:) = filter .and. a_hx_f(i,:) /= spval
            ENDDO
         ENDIF
         CALL write_DA_history_variable_3d (HistForm, DEF_hist_vars%DA_hx_f, &
            a_hx_f, file_hist, 'f_hx_f', itime_in_file, &
            'source', 1, nsource, sumarea, filter, &
            'Prior H(x) by DA observation source', 'source-dependent', &
            filter_3d = filter_hx)

         IF (p_is_worker) THEN
            DO i = 1, nsource
               filter_hx(i,:) = filter .and. a_hx_a(i,:) /= spval
            ENDDO
         ENDIF
         CALL write_DA_history_variable_3d (HistForm, DEF_hist_vars%DA_hx_a, &
            a_hx_a, file_hist, 'f_hx_a', itime_in_file, &
            'source', 1, nsource, sumarea, filter, &
            'Analysis H(x) by DA observation source', 'source-dependent', &
            filter_3d = filter_hx)

         deallocate (filter_hx)

         IF (DEF_DA_ENS_NUM > 1) THEN
            CALL write_DA_history_variable_4d (HistForm, DEF_hist_vars%DA_hx_f_ens, &
               a_hx_f_ens, file_hist, 'f_hx_f_ens', itime_in_file, &
               'source', 1, nsource, 'ens', 1, DEF_DA_ENS_NUM, &
               sumarea, filter, 'Prior ensemble H(x) by DA observation source', &
               'source-dependent')
            CALL write_DA_history_variable_4d (HistForm, DEF_hist_vars%DA_hx_a_ens, &
               a_hx_a_ens, file_hist, 'f_hx_a_ens', itime_in_file, &
               'source', 1, nsource, 'ens', 1, DEF_DA_ENS_NUM, &
               sumarea, filter, 'Analysis ensemble H(x) by DA observation source', &
               'source-dependent')
         ENDIF
      ENDIF

#ifdef GridRiverLakeFlow
      ! --------------------------------------------------------------------
      ! River discharge on unit catchments
      ! --------------------------------------------------------------------
      IF (DEF_DA_ENS_NUM > 1) THEN
         ! Same unitcat history file as hist_grid_riverlake_out (which creates it).
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
      ENDIF
#endif

   END SUBROUTINE hist_DA_history_out

!-----------------------------------------------------------------------------

   SUBROUTINE write_DA_history_variable_3d ( HistForm, is_hist, &
         acc_vec, file_hist, varname, itime_in_file, dim1name, lb1, ndim1, &
         sumarea, filter, longname, units, filter_3d)

   IMPLICIT NONE

   character(len=*), intent(in) :: HistForm
   logical, intent(in) :: is_hist

   real(r8), intent(inout) :: acc_vec(:,:)
   character(len=*), intent(in) :: file_hist
   character(len=*), intent(in) :: varname
   integer, intent(in) :: itime_in_file
   character(len=*), intent(in) :: dim1name
   integer, intent(in) :: lb1, ndim1

   type(block_data_real8_2d), intent(in) :: sumarea
   logical, intent(in) :: filter(:)
   logical, intent(in), optional :: filter_3d(:,:)
   character (len=*), intent(in) :: longname
   character (len=*), intent(in) :: units

      IF (.not. is_hist) RETURN

      select CASE (HistForm)
      CASE ('Gridded')
         IF (present(filter_3d)) THEN
            CALL flux_map_and_write_3d ( &
               acc_vec, file_hist, varname, itime_in_file, dim1name, lb1, ndim1, &
               sumarea, filter, longname, units, filter_3d)
         ELSE
            CALL flux_map_and_write_3d ( &
               acc_vec, file_hist, varname, itime_in_file, dim1name, lb1, ndim1, &
               sumarea, filter, longname, units)
         ENDIF
      END select

   END SUBROUTINE write_DA_history_variable_3d

!-----------------------------------------------------------------------------

   SUBROUTINE write_DA_history_variable_4d ( HistForm, is_hist,   &
         acc_vec, file_hist, varname, itime_in_file, &
         dim1name, lb1, ndim1, dim2name, lb2, ndim2, &
         sumarea, filter, longname, units, acc_num)

   USE MOD_Vars_1DAccFluxes, only: nac

   IMPLICIT NONE

   character(len=*), intent(in) :: HistForm
   logical, intent(in) :: is_hist

   real(r8), intent(inout) :: acc_vec(:,:,:)
   character(len=*), intent(in) :: file_hist
   character(len=*), intent(in) :: varname
   integer, intent(in) :: itime_in_file
   character(len=*), intent(in) :: dim1name, dim2name
   integer, intent(in) :: lb1, ndim1, lb2, ndim2

   type(block_data_real8_2d), intent(in) :: sumarea
   logical, intent(in) :: filter(:)
   character (len=*), intent(in) :: longname
   character (len=*), intent(in) :: units
   real(r8), intent(in), optional:: acc_num(:)

   integer :: i1, i2

      IF (.not. is_hist) RETURN

      IF ( .not. present(acc_num) ) THEN
         IF (p_is_worker) &
            WHERE (acc_vec /= spval) acc_vec = acc_vec / nac
      ELSE
         IF (p_is_worker) THEN
            DO i1 = lbound(acc_vec,1), ubound(acc_vec,1)
               DO i2 = lbound(acc_vec,2), ubound(acc_vec,2)
                  WHERE (acc_vec(i1,i2,:)/=spval .and. acc_num>0) &
                        acc_vec(i1,i2,:) = acc_vec(i1,i2,:) / acc_num
               ENDDO
            ENDDO
         ENDIF
      ENDIF

      select CASE (HistForm)
      CASE ('Gridded')
         CALL flux_map_and_write_4d ( &
            acc_vec, file_hist, varname, itime_in_file, dim1name, lb1, ndim1, &
            dim2name, lb2, ndim2, sumarea, filter, longname, units)
      END select

   END SUBROUTINE write_DA_history_variable_4d

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Vars_Hist
#endif
