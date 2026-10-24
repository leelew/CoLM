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
   USE MOD_WorkerPushData
   USE MOD_Vector_ReadWrite
   USE MOD_Grid_RiverLakeNetwork
   USE MOD_Grid_RiverLakeHist
   USE MOD_Grid_Reservoir
#endif
   IMPLICIT NONE
   SAVE

   PUBLIC :: hist_DA_history_out

#ifdef GridRiverLakeFlow
   PUBLIC :: hist_DA_riverlake_init
   PUBLIC :: hist_DA_riverlake_out
   PUBLIC :: hist_DA_riverlake_final

!=============================================================================
! Routing output accumulators: window sums shared by the OL and a-track calls,
! and the patch-level values remapped from them at history output.
!=============================================================================
   real(r8), allocatable :: a_discharge_ol     (:) ! open-loop discharge sum                            [m^3/s*s]
   real(r8), allocatable :: acctime_ucat_ol    (:) ! open-loop averaging time                           [s]
   real(r8), allocatable :: a_wdsrf_ucat_ol    (:) ! open-loop water depth sum                          [m*s]
   real(r8), allocatable :: a_veloc_riv_ol     (:) ! open-loop water velocity sum                       [m/s*s]
   real(r8), allocatable :: a_floodarea_ol     (:) ! open-loop flooded area sum                         [m^2*s]
   real(r8), allocatable :: a_discharge_a      (:) ! prior & posterior discharge sum                    [m^3/s*s]
   real(r8), allocatable :: acctime_ucat_a     (:) ! prior & posterior averaging time                   [s]
   real(r8), allocatable :: a_wdsrf_ucat_a     (:) ! prior & posterior water depth sum                  [m*s]
   real(r8), allocatable :: a_veloc_riv_a      (:) ! prior & posterior water velocity sum               [m/s*s]
   real(r8), allocatable :: a_floodarea_a      (:) ! prior & posterior flooded area sum                 [m^2*s]

   real(r8), allocatable :: a_discharge_pch_a  (:) ! prior & posterior discharge on patches             [m^3/s]
   real(r8), allocatable :: a_dis_rmth_pch_a   (:) ! prior & posterior river mouth discharge on patches [m^3/s]
   real(r8), allocatable :: a_wdsrf_ucat_pch_a (:) ! prior & posterior water depth on patches           [m]
   real(r8), allocatable :: a_veloc_riv_pch_a  (:) ! prior & posterior water velocity on patches        [m/s]
   real(r8), allocatable :: a_floodfrc_pch_a   (:) ! prior & posterior flooded fraction on patches      [-]

   real(r8), allocatable :: acctime_resv_ol    (:) ! open-loop averaging time                           [s]
   real(r8), allocatable :: a_volresv_ol       (:) ! open-loop reservoir volume sum                     [m^3*s]
   real(r8), allocatable :: a_qresv_in_ol      (:) ! open-loop reservoir inflow sum                     [m^3/s*s]
   real(r8), allocatable :: a_qresv_out_ol     (:) ! open-loop reservoir outflow sum                    [m^3/s*s]
   real(r8), allocatable :: acctime_resv_a     (:) ! prior & posterior averaging time                   [s]
   real(r8), allocatable :: a_volresv_a        (:) ! prior & posterior reservoir volume sum             [m^3*s]
   real(r8), allocatable :: a_qresv_in_a       (:) ! prior & posterior reservoir inflow sum             [m^3/s*s]
   real(r8), allocatable :: a_qresv_out_a      (:) ! prior & posterior reservoir outflow sum            [m^3/s*s]
#endif

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

#ifdef GridRiverLakeFlow
   SUBROUTINE hist_DA_riverlake_init ()

!-----------------------------------------------------------------------------
   IMPLICIT NONE

!-----------------------------------------------------------------------------
      IF (p_is_worker) THEN
         IF (numucat > 0) THEN
            allocate (a_discharge_ol      (numucat)); a_discharge_ol      (:) = spval
            allocate (acctime_ucat_ol     (numucat)); acctime_ucat_ol     (:) = spval
            allocate (a_wdsrf_ucat_ol     (numucat)); a_wdsrf_ucat_ol     (:) = spval
            allocate (a_veloc_riv_ol      (numucat)); a_veloc_riv_ol      (:) = spval
            allocate (a_floodarea_ol      (numucat)); a_floodarea_ol      (:) = spval
            allocate (a_discharge_a       (numucat)); a_discharge_a       (:) = 0.
            allocate (acctime_ucat_a      (numucat)); acctime_ucat_a      (:) = 0.
            allocate (a_wdsrf_ucat_a      (numucat)); a_wdsrf_ucat_a      (:) = 0.
            allocate (a_veloc_riv_a       (numucat)); a_veloc_riv_a       (:) = 0.
            allocate (a_floodarea_a       (numucat)); a_floodarea_a       (:) = 0.
         ENDIF
         IF (numpatch > 0) THEN
            allocate (a_discharge_pch_a     (numpatch)); a_discharge_pch_a     (:) = spval
            allocate (a_dis_rmth_pch_a      (numpatch)); a_dis_rmth_pch_a      (:) = spval
            allocate (a_wdsrf_ucat_pch_a    (numpatch)); a_wdsrf_ucat_pch_a    (:) = spval
            allocate (a_veloc_riv_pch_a     (numpatch)); a_veloc_riv_pch_a     (:) = spval
            allocate (a_floodfrc_pch_a      (numpatch)); a_floodfrc_pch_a      (:) = spval
         ENDIF
         IF (numresv > 0) THEN
            allocate (acctime_resv_ol     (numresv)); acctime_resv_ol     (:) = spval
            allocate (a_volresv_ol        (numresv)); a_volresv_ol        (:) = spval
            allocate (a_qresv_in_ol       (numresv)); a_qresv_in_ol       (:) = spval
            allocate (a_qresv_out_ol      (numresv)); a_qresv_out_ol      (:) = spval
            allocate (acctime_resv_a      (numresv)); acctime_resv_a      (:) = 0.
            allocate (a_volresv_a         (numresv)); a_volresv_a         (:) = 0.
            allocate (a_qresv_in_a        (numresv)); a_qresv_in_a        (:) = 0.
            allocate (a_qresv_out_a       (numresv)); a_qresv_out_a       (:) = 0.
         ENDIF
      ENDIF

   END SUBROUTINE hist_DA_riverlake_init

!-----------------------------------------------------------------------------

   SUBROUTINE hist_DA_riverlake_final ()

!-----------------------------------------------------------------------------
   IMPLICIT NONE

!-----------------------------------------------------------------------------
      IF (allocated(a_discharge_ol      )) deallocate (a_discharge_ol      )
      IF (allocated(acctime_ucat_ol     )) deallocate (acctime_ucat_ol     )
      IF (allocated(a_wdsrf_ucat_ol     )) deallocate (a_wdsrf_ucat_ol     )
      IF (allocated(a_veloc_riv_ol      )) deallocate (a_veloc_riv_ol      )
      IF (allocated(a_floodarea_ol      )) deallocate (a_floodarea_ol      )
      IF (allocated(a_discharge_a       )) deallocate (a_discharge_a       )
      IF (allocated(acctime_ucat_a      )) deallocate (acctime_ucat_a      )
      IF (allocated(a_wdsrf_ucat_a      )) deallocate (a_wdsrf_ucat_a      )
      IF (allocated(a_veloc_riv_a       )) deallocate (a_veloc_riv_a       )
      IF (allocated(a_floodarea_a       )) deallocate (a_floodarea_a       )
      IF (allocated(a_discharge_pch_a   )) deallocate (a_discharge_pch_a   )
      IF (allocated(a_dis_rmth_pch_a    )) deallocate (a_dis_rmth_pch_a    )
      IF (allocated(a_wdsrf_ucat_pch_a  )) deallocate (a_wdsrf_ucat_pch_a  )
      IF (allocated(a_veloc_riv_pch_a   )) deallocate (a_veloc_riv_pch_a   )
      IF (allocated(a_floodfrc_pch_a    )) deallocate (a_floodfrc_pch_a    )
      IF (allocated(acctime_resv_ol     )) deallocate (acctime_resv_ol     )
      IF (allocated(a_volresv_ol        )) deallocate (a_volresv_ol        )
      IF (allocated(a_qresv_in_ol       )) deallocate (a_qresv_in_ol       )
      IF (allocated(a_qresv_out_ol      )) deallocate (a_qresv_out_ol      )
      IF (allocated(acctime_resv_a      )) deallocate (acctime_resv_a      )
      IF (allocated(a_volresv_a         )) deallocate (a_volresv_a         )
      IF (allocated(a_qresv_in_a        )) deallocate (a_qresv_in_a        )
      IF (allocated(a_qresv_out_a       )) deallocate (a_qresv_out_a       )

   END SUBROUTINE hist_DA_riverlake_final

!-----------------------------------------------------------------------------

   SUBROUTINE hist_DA_riverlake_out (file_hist, idate)

!-----------------------------------------------------------------------------
   IMPLICIT NONE

   character(len=*), intent(in) :: file_hist
   integer, intent(in) :: idate(3)

   character(len=256) :: file_hist_ucat
   real(r8), allocatable :: acc_vec_grid (:)
   real(r8), allocatable :: a_floodfrc_ucat_a    (:)
   real(r8), allocatable :: a_floodfrc_inpm_a    (:)
   integer :: itime_in_file_ucat, i

!-----------------------------------------------------------------------------
      IF (DEF_DA_ENS_NUM <= 1) RETURN

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

      IF (p_is_worker) THEN
         IF (numinpm  > 0) allocate (acc_vec_grid (numinpm ))
      ENDIF

      IF (DEF_hist_vars%riv_height) THEN
         IF (p_is_worker) THEN
            IF (numucat > 0) THEN
               WHERE (acctime_ucat_a > 0)
                  a_wdsrf_ucat_a = a_wdsrf_ucat_a / acctime_ucat_a
               ELSE WHERE
                  a_wdsrf_ucat_a = spval
               END WHERE
            ENDIF
            CALL worker_push_data (push_ucat2grid, a_wdsrf_ucat_a, acc_vec_grid, fillvalue = spval)
            CALL worker_remap_data_grid2pset (remap_patch2inpm, acc_vec_grid, a_wdsrf_ucat_pch_a, &
               fillvalue = spval, mode = 'average')
         ENDIF

         CALL vector_gather_map2grid_and_write ( a_wdsrf_ucat_a, numucat,          &
            totalnumucat, ucat_data_address, griducat%nlon, x_ucat, griducat%nlat, y_ucat, &
            file_hist_ucat, 'f_wdpth_ucat_a', 'lon_ucat', 'lat_ucat', itime_in_file_ucat, &
            'deepest water depth in river and flood plain from the ensemble-mean state', 'm')
      ENDIF

      IF (DEF_hist_vars%riv_veloct) THEN
         IF (p_is_worker) THEN
            IF (numucat > 0) THEN
               WHERE (acctime_ucat_a > 0)
                  a_veloc_riv_a = a_veloc_riv_a / acctime_ucat_a
               ELSE WHERE
                  a_veloc_riv_a = spval
               END WHERE
            ENDIF
            CALL worker_push_data (push_ucat2grid, a_veloc_riv_a, acc_vec_grid, fillvalue = spval)
            CALL worker_remap_data_grid2pset (remap_patch2inpm, acc_vec_grid, a_veloc_riv_pch_a, &
               fillvalue = spval, mode = 'average')
         ENDIF

         CALL vector_gather_map2grid_and_write ( a_veloc_riv_a, numucat,           &
            totalnumucat, ucat_data_address, griducat%nlon, x_ucat, griducat%nlat, y_ucat, &
            file_hist_ucat, 'f_veloc_riv_a', 'lon_ucat', 'lat_ucat', itime_in_file_ucat, &
            'water velocity in river from the ensemble-mean state', 'm/s')
      ENDIF

      IF (DEF_hist_vars%discharge) THEN
         IF (p_is_worker) THEN
            IF (numucat > 0) THEN
               WHERE (acctime_ucat_a > 0)
                  a_discharge_a = a_discharge_a / acctime_ucat_a
               ELSE WHERE
                  a_discharge_a = spval
               END WHERE
            ENDIF

            CALL worker_push_data (push_ucat2grid, a_discharge_a, acc_vec_grid, fillvalue = spval)

            IF (numinpm > 0)  THEN
               WHERE ((sum_grid_area /= spval) .and. (acc_vec_grid /= spval))
                  acc_vec_grid = acc_vec_grid / sum_grid_area
               ELSE WHERE
                  acc_vec_grid = spval
               END WHERE
            ENDIF

            CALL worker_remap_data_grid2pset (remap_patch2inpm, acc_vec_grid, a_discharge_pch_a, &
               fillvalue = spval, mode = 'sum')

         ENDIF

         CALL vector_gather_map2grid_and_write ( a_discharge_a, numucat,       &
            totalnumucat, ucat_data_address, griducat%nlon, x_ucat, griducat%nlat, y_ucat, &
            file_hist_ucat, 'f_discharge_a', 'lon_ucat', 'lat_ucat', itime_in_file_ucat, &
            'discharge in river and flood plain from the ensemble-mean state', 'm^3/s')

         ! ----- river mouth discharge: same accumulator, ocean outlets only -----
         ! The patch-level regrid uses the river mouth area, not sum_grid_area.
         IF (p_is_worker) THEN
            IF (numucat > 0) THEN
               WHERE (ucat_next /= -9) a_discharge_a = spval
            ENDIF

            CALL worker_push_data (push_ucat2grid, a_discharge_a, acc_vec_grid, fillvalue = spval)

            IF (numinpm > 0)  THEN
               WHERE ((sum_rmth_area /= spval) .and. (acc_vec_grid /= spval))
                  acc_vec_grid = acc_vec_grid / sum_rmth_area
               ELSE WHERE
                  acc_vec_grid = spval
               END WHERE
            ENDIF

            CALL worker_remap_data_grid2pset (remap_patch2inpm, acc_vec_grid, a_dis_rmth_pch_a, &
               fillvalue = spval, mode = 'sum')

         ENDIF

         CALL vector_gather_map2grid_and_write ( a_discharge_a, numucat,          &
            totalnumucat, ucat_data_address, griducat%nlon, x_ucat, griducat%nlat, y_ucat, &
            file_hist_ucat, 'f_discharge_rivermouth_a', 'lon_ucat', 'lat_ucat', itime_in_file_ucat, &
            'river mouth discharge into ocean from the ensemble-mean state', 'm^3/s')

      ENDIF

      IF (DEF_hist_vars%floodfrc .or. DEF_hist_vars%floodarea) THEN

         IF (p_is_worker) THEN

            IF (numucat > 0) THEN
               WHERE (acctime_ucat_a > 0)
                  a_floodarea_a = a_floodarea_a / acctime_ucat_a
               ELSE WHERE
                  a_floodarea_a = spval
               END WHERE

               allocate (a_floodfrc_ucat_a (numucat))
               WHERE (topo_area > 0)
                  a_floodfrc_ucat_a = a_floodarea_a / topo_area
               ELSE WHERE
                  a_floodfrc_ucat_a = spval
               END WHERE
            ENDIF

            IF (numinpm > 0) allocate (a_floodfrc_inpm_a (numinpm))

            CALL worker_push_data (push_ucat2inpm, a_floodfrc_ucat_a, a_floodfrc_inpm_a, &
               fillvalue = spval, mode = 'average')

            CALL worker_remap_data_grid2pset (remap_patch2inpm, a_floodfrc_inpm_a, &
               a_floodfrc_pch_a, fillvalue = spval, mode = 'average')

            IF (allocated (a_floodfrc_ucat_a)) deallocate (a_floodfrc_ucat_a)
            IF (allocated (a_floodfrc_inpm_a)) deallocate (a_floodfrc_inpm_a)
         ENDIF

      ENDIF

      IF (allocated (acc_vec_grid   )) deallocate (acc_vec_grid   )

      ! ----- reservoir variables, same layout as hist_grid_riverlake_out -----
      IF (DEF_Reservoir_Method > 0 .and. totalnumresv > 0) THEN

         IF (DEF_hist_vars%volresv) THEN
            IF (p_is_worker) THEN
               IF (numresv > 0) THEN
                  WHERE (acctime_resv_a > 0)
                     a_volresv_a = a_volresv_a / acctime_resv_a
                  ELSEWHERE
                     a_volresv_a = spval
                  END WHERE
               ENDIF
            ENDIF

            CALL vector_gather_and_write ( a_volresv_a, numresv, totalnumresv, resv_data_address, &
               file_hist_ucat, 'volresv_a', 'reservoir', itime_in_file_ucat, &
               'reservoir water volume from the ensemble-mean state', 'm^3')
         ENDIF

         IF (DEF_hist_vars%qresv_in) THEN
            IF (p_is_worker) THEN
               IF (numresv > 0) THEN
                  WHERE (acctime_resv_a > 0)
                     a_qresv_in_a = a_qresv_in_a / acctime_resv_a
                  ELSEWHERE
                     a_qresv_in_a = spval
                  END WHERE
               ENDIF
            ENDIF

            CALL vector_gather_and_write ( a_qresv_in_a, numresv, totalnumresv, resv_data_address, &
               file_hist_ucat, 'qresv_in_a', 'reservoir', itime_in_file_ucat, &
               'reservoir inflow from the ensemble-mean state', 'm^3/s')
         ENDIF

         IF (DEF_hist_vars%qresv_out) THEN
            IF (p_is_worker) THEN
               IF (numresv > 0) THEN
                  WHERE (acctime_resv_a > 0)
                     a_qresv_out_a = a_qresv_out_a / acctime_resv_a
                  ELSEWHERE
                     a_qresv_out_a = spval
                  END WHERE
               ENDIF
            ENDIF

            CALL vector_gather_and_write ( a_qresv_out_a, numresv, totalnumresv, resv_data_address, &
               file_hist_ucat, 'qresv_out_a', 'reservoir', itime_in_file_ucat, &
               'reservoir outflow from the ensemble-mean state', 'm^3/s')
         ENDIF

      ENDIF

      IF (p_is_worker) THEN
         IF (numucat > 0) THEN
            a_discharge_a     (:) = 0.
            acctime_ucat_a    (:) = 0.
            a_wdsrf_ucat_a    (:) = 0.
            a_veloc_riv_a     (:) = 0.
            a_floodarea_a     (:) = 0.
         ENDIF
         IF (numresv > 0) THEN
            acctime_resv_a    (:) = 0.
            a_volresv_a       (:) = 0.
            a_qresv_in_a      (:) = 0.
            a_qresv_out_a     (:) = 0.
         ENDIF
      ENDIF

   END SUBROUTINE hist_DA_riverlake_out

!-----------------------------------------------------------------------------

   SUBROUTINE hist_DA_riverlake_write_patch (file_hist, HistForm, itime_in_file, sumarea_one, nac_one)

!-----------------------------------------------------------------------------
   IMPLICIT NONE

   character(len=*), intent(in) :: file_hist
   character(len=*), intent(in) :: HistForm
   integer,          intent(in) :: itime_in_file
   type(block_data_real8_2d), intent(in) :: sumarea_one
   real(r8),         intent(in) :: nac_one(:)

!-----------------------------------------------------------------------------
      IF (DEF_DA_ENS_NUM <= 1) RETURN

      CALL write_DA_history_variable_2d ( HistForm, DEF_hist_vars%riv_height, &
         a_wdsrf_ucat_pch_a, file_hist, 'f_wdpth_ucat_regrid_a', itime_in_file, &
         sumarea_ucat, filter_ucat, &
         'regridded deepest water depth in river and flood plain from the ensemble-mean state', &
         'm', nac_one)

      CALL write_DA_history_variable_2d ( HistForm, DEF_hist_vars%riv_veloct, &
         a_veloc_riv_pch_a, file_hist, 'f_veloc_riv_regrid_a', itime_in_file, &
         sumarea_ucat, filter_ucat, &
         'regridded water velocity in river from the ensemble-mean state', &
         'm/s', nac_one)

      CALL write_DA_history_variable_2d ( HistForm, DEF_hist_vars%discharge, &
         a_discharge_pch_a, file_hist, 'f_discharge_a', itime_in_file, &
         sumarea_one, filter_ucat, &
         'regridded discharge in river and flood plain from the ensemble-mean state', &
         'm^3/s', nac_one, input_mode = 'total')

      CALL write_DA_history_variable_2d ( HistForm, DEF_hist_vars%discharge, &
         a_dis_rmth_pch_a, file_hist, 'f_discharge_rivermouth_regrid_a', itime_in_file, &
         sumarea_one, filter_ucat, &
         'regridded river mouth discharge into ocean from the ensemble-mean state', &
         'm^3/s', nac_one, input_mode = 'total')

      CALL write_DA_history_variable_2d ( HistForm, DEF_hist_vars%floodfrc, &
         a_floodfrc_pch_a, file_hist, 'f_floodfrc_a', itime_in_file, &
         sumarea_inpm, filter_inpm, &
         'flooded area fraction from the ensemble-mean state', '100%', nac_one)

      IF (trim(HistForm) == 'Gridded') THEN
         CALL write_DA_history_variable_2d ( HistForm, DEF_hist_vars%floodarea, &
            a_floodfrc_pch_a, file_hist, 'f_floodarea_a', itime_in_file, &
            sumarea_one, filter_inpm, &
            'flooded area from the ensemble-mean state', 'km^2', nac_one)
      ENDIF

   END SUBROUTINE hist_DA_riverlake_write_patch

#endif

!-----------------------------------------------------------------------------

   SUBROUTINE hist_DA_history_out (file_hist, HistForm, itime_in_file, &
         sumarea, filter, sumarea_one, nac_one)

!-----------------------------------------------------------------------------
   IMPLICIT NONE

   character(len=*), intent(in) :: file_hist
   character(len=*), intent(in) :: HistForm
   integer,          intent(in) :: itime_in_file

   type(block_data_real8_2d), intent(inout) :: sumarea
   logical,                   intent(inout) :: filter(:)
   type(block_data_real8_2d), intent(in), optional :: sumarea_one
   real(r8), allocatable,     intent(in), optional :: nac_one(:)

   logical, allocatable :: filter_hx(:,:)
   integer :: i

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
      ! the order of DEF_DA_OBS_* entries whose TARGET is 'SM' or 'ST'.
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
      ! River and lake diagnostics on patches
      ! --------------------------------------------------------------------

      IF (DEF_DA_ENS_NUM > 1) THEN
         CALL hist_DA_riverlake_write_patch (file_hist, HistForm, itime_in_file, sumarea_one, nac_one)
      ENDIF
#endif

   END SUBROUTINE hist_DA_history_out

!-----------------------------------------------------------------------------

   SUBROUTINE write_DA_history_variable_2d ( HistForm, is_hist, &
         acc_vec, file_hist, varname, itime_in_file, sumarea, filter, &
         longname, units, acc_num, input_mode)

   USE MOD_Vars_1DAccFluxes, only: nac

   IMPLICIT NONE

   character(len=*), intent(in) :: HistForm
   logical, intent(in) :: is_hist

   real(r8),         intent(inout) :: acc_vec(:)
   character(len=*), intent(in)    :: file_hist
   character(len=*), intent(in)    :: varname
   integer,          intent(in)    :: itime_in_file
   character(len=*), intent(in)    :: longname
   character(len=*), intent(in)    :: units

   type(block_data_real8_2d), intent(in) :: sumarea
   logical, intent(in) :: filter(:)
   real(r8), intent(in), optional  :: acc_num(:)

   character(len=*), intent(in), optional :: input_mode

      IF (.not. is_hist) RETURN

      IF ( .not. present(acc_num) ) THEN
         IF (p_is_worker) &
            WHERE (acc_vec /= spval) acc_vec = acc_vec / nac
      ELSE
         IF (p_is_worker)  &
            WHERE (acc_vec/=spval .and. acc_num>0) acc_vec = acc_vec / acc_num
      ENDIF

      select CASE (HistForm)
      CASE ('Gridded')
         IF (present(input_mode)) THEN
            CALL flux_map_and_write_2d ( &
               acc_vec, file_hist, varname, itime_in_file, sumarea, filter, longname, units, input_mode)
         ELSE
            CALL flux_map_and_write_2d ( &
               acc_vec, file_hist, varname, itime_in_file, sumarea, filter, longname, units)
         ENDIF
      END select

   END SUBROUTINE write_DA_history_variable_2d

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
