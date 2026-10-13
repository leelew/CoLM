#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Hist
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
   USE MOD_Vars_Global, only: spval
   USE MOD_NetCDFSerial, only: ncio_write_time
#ifdef GridRiverLakeFlow
   USE MOD_WorkerPushData, only: worker_push_data, worker_remap_data_grid2pset
   USE MOD_Vector_ReadWrite, only: vector_gather_map2grid_and_write, vector_gather_and_write
   USE MOD_Grid_RiverLakeNetwork, only: numucat, numinpm, totalnumucat, ucat_data_address, &
      griducat, x_ucat, y_ucat, topo_area, ucat_next, push_ucat2grid, push_ucat2inpm, remap_patch2inpm
   USE MOD_Grid_RiverLakeHist, only: sum_grid_area, sum_rmth_area
   USE MOD_Grid_Reservoir, only: numresv, totalnumresv, resv_data_address
   USE MOD_DA_Vars_TimeVariables, only: a_discharge_mean, acctime_ucat_mean, a_discharge_mean_pch, &
      a_dis_rmth_mean_pch, &
      a_wdsrf_ucat_mean, a_veloc_riv_mean, a_floodarea_mean, &
      a_wdsrf_ucat_mean_pch, a_veloc_riv_mean_pch, a_floodfrc_mean_pch, &
      a_volresv_mean, a_qresv_in_mean, a_qresv_out_mean, acctime_resv_mean
#endif
   IMPLICIT NONE
   SAVE

#ifdef GridRiverLakeFlow
   PUBLIC :: hist_DA_riverlake_out
#endif

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

#ifdef GridRiverLakeFlow
   SUBROUTINE hist_DA_riverlake_out (file_hist, idate)

!-----------------------------------------------------------------------------
   IMPLICIT NONE

   character(len=*), intent(in) :: file_hist
   integer, intent(in) :: idate(3)

   character(len=256) :: file_hist_ucat
   real(r8), allocatable :: acc_vec_grid (:)
   real(r8), allocatable :: a_floodfrc_mean_ucat (:)
   real(r8), allocatable :: a_floodfrc_mean_inpm (:)
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
               WHERE (acctime_ucat_mean > 0)
                  a_wdsrf_ucat_mean = a_wdsrf_ucat_mean / acctime_ucat_mean
               ELSE WHERE
                  a_wdsrf_ucat_mean = spval
               END WHERE
            ENDIF
            CALL worker_push_data (push_ucat2grid, a_wdsrf_ucat_mean, acc_vec_grid, fillvalue = spval)
            CALL worker_remap_data_grid2pset (remap_patch2inpm, acc_vec_grid, a_wdsrf_ucat_mean_pch, &
               fillvalue = spval, mode = 'average')
         ENDIF

         CALL vector_gather_map2grid_and_write ( a_wdsrf_ucat_mean, numucat,          &
            totalnumucat, ucat_data_address, griducat%nlon, x_ucat, griducat%nlat, y_ucat, &
            file_hist_ucat, 'f_wdpth_ucat_mean', 'lon_ucat', 'lat_ucat', itime_in_file_ucat, &
            'deepest water depth in river and flood plain from the ensemble-mean state', 'm')
      ENDIF

      IF (DEF_hist_vars%riv_veloct) THEN
         IF (p_is_worker) THEN
            IF (numucat > 0) THEN
               WHERE (acctime_ucat_mean > 0)
                  a_veloc_riv_mean = a_veloc_riv_mean / acctime_ucat_mean
               ELSE WHERE
                  a_veloc_riv_mean = spval
               END WHERE
            ENDIF
            CALL worker_push_data (push_ucat2grid, a_veloc_riv_mean, acc_vec_grid, fillvalue = spval)
            CALL worker_remap_data_grid2pset (remap_patch2inpm, acc_vec_grid, a_veloc_riv_mean_pch, &
               fillvalue = spval, mode = 'average')
         ENDIF

         CALL vector_gather_map2grid_and_write ( a_veloc_riv_mean, numucat,           &
            totalnumucat, ucat_data_address, griducat%nlon, x_ucat, griducat%nlat, y_ucat, &
            file_hist_ucat, 'f_veloc_riv_mean', 'lon_ucat', 'lat_ucat', itime_in_file_ucat, &
            'water velocity in river from the ensemble-mean state', 'm/s')
      ENDIF

      IF (DEF_hist_vars%discharge) THEN
         IF (p_is_worker) THEN
            IF (numucat > 0) THEN
               WHERE (acctime_ucat_mean > 0)
                  a_discharge_mean = a_discharge_mean / acctime_ucat_mean
               ELSE WHERE
                  a_discharge_mean = spval
               END WHERE
            ENDIF

            CALL worker_push_data (push_ucat2grid, a_discharge_mean, acc_vec_grid, fillvalue = spval)

            IF (numinpm > 0)  THEN
               WHERE ((sum_grid_area /= spval) .and. (acc_vec_grid /= spval))
                  acc_vec_grid = acc_vec_grid / sum_grid_area
               ELSE WHERE
                  acc_vec_grid = spval
               END WHERE
            ENDIF

            CALL worker_remap_data_grid2pset (remap_patch2inpm, acc_vec_grid, a_discharge_mean_pch, &
               fillvalue = spval, mode = 'sum')

         ENDIF

         CALL vector_gather_map2grid_and_write ( a_discharge_mean, numucat,       &
            totalnumucat, ucat_data_address, griducat%nlon, x_ucat, griducat%nlat, y_ucat, &
            file_hist_ucat, 'f_discharge_mean', 'lon_ucat', 'lat_ucat', itime_in_file_ucat, &
            'discharge in river and flood plain from the ensemble-mean state', 'm^3/s')

         ! ----- river mouth discharge: same accumulator, ocean outlets only -----
         ! The patch-level regrid uses the river mouth area, not sum_grid_area.
         IF (p_is_worker) THEN
            IF (numucat > 0) THEN
               WHERE (ucat_next /= -9) a_discharge_mean = spval
            ENDIF

            CALL worker_push_data (push_ucat2grid, a_discharge_mean, acc_vec_grid, fillvalue = spval)

            IF (numinpm > 0)  THEN
               WHERE ((sum_rmth_area /= spval) .and. (acc_vec_grid /= spval))
                  acc_vec_grid = acc_vec_grid / sum_rmth_area
               ELSE WHERE
                  acc_vec_grid = spval
               END WHERE
            ENDIF

            CALL worker_remap_data_grid2pset (remap_patch2inpm, acc_vec_grid, a_dis_rmth_mean_pch, &
               fillvalue = spval, mode = 'sum')

         ENDIF

         CALL vector_gather_map2grid_and_write ( a_discharge_mean, numucat,          &
            totalnumucat, ucat_data_address, griducat%nlon, x_ucat, griducat%nlat, y_ucat, &
            file_hist_ucat, 'f_discharge_rivermouth_mean', 'lon_ucat', 'lat_ucat', itime_in_file_ucat, &
            'river mouth discharge into ocean from the ensemble-mean state', 'm^3/s')

      ENDIF

      IF (DEF_hist_vars%floodfrc .or. DEF_hist_vars%floodarea) THEN

         IF (p_is_worker) THEN

            IF (numucat > 0) THEN
               WHERE (acctime_ucat_mean > 0)
                  a_floodarea_mean = a_floodarea_mean / acctime_ucat_mean
               ELSE WHERE
                  a_floodarea_mean = spval
               END WHERE

               allocate (a_floodfrc_mean_ucat (numucat))
               WHERE (topo_area > 0)
                  a_floodfrc_mean_ucat = a_floodarea_mean / topo_area
               ELSE WHERE
                  a_floodfrc_mean_ucat = spval
               END WHERE
            ENDIF

            IF (numinpm > 0) allocate (a_floodfrc_mean_inpm (numinpm))

            CALL worker_push_data (push_ucat2inpm, a_floodfrc_mean_ucat, a_floodfrc_mean_inpm, &
               fillvalue = spval, mode = 'average')

            CALL worker_remap_data_grid2pset (remap_patch2inpm, a_floodfrc_mean_inpm, &
               a_floodfrc_mean_pch, fillvalue = spval, mode = 'average')

            IF (allocated (a_floodfrc_mean_ucat)) deallocate (a_floodfrc_mean_ucat)
            IF (allocated (a_floodfrc_mean_inpm)) deallocate (a_floodfrc_mean_inpm)
         ENDIF

      ENDIF

      IF (allocated (acc_vec_grid   )) deallocate (acc_vec_grid   )

      ! ----- reservoir variables, same layout as hist_grid_riverlake_out -----
      IF (DEF_Reservoir_Method > 0 .and. totalnumresv > 0) THEN

         IF (DEF_hist_vars%volresv) THEN
            IF (p_is_worker) THEN
               IF (numresv > 0) THEN
                  WHERE (acctime_resv_mean > 0)
                     a_volresv_mean = a_volresv_mean / acctime_resv_mean
                  ELSEWHERE
                     a_volresv_mean = spval
                  END WHERE
               ENDIF
            ENDIF

            CALL vector_gather_and_write ( a_volresv_mean, numresv, totalnumresv, resv_data_address, &
               file_hist_ucat, 'volresv_mean', 'reservoir', itime_in_file_ucat, &
               'reservoir water volume from the ensemble-mean state', 'm^3')
         ENDIF

         IF (DEF_hist_vars%qresv_in) THEN
            IF (p_is_worker) THEN
               IF (numresv > 0) THEN
                  WHERE (acctime_resv_mean > 0)
                     a_qresv_in_mean = a_qresv_in_mean / acctime_resv_mean
                  ELSEWHERE
                     a_qresv_in_mean = spval
                  END WHERE
               ENDIF
            ENDIF

            CALL vector_gather_and_write ( a_qresv_in_mean, numresv, totalnumresv, resv_data_address, &
               file_hist_ucat, 'qresv_in_mean', 'reservoir', itime_in_file_ucat, &
               'reservoir inflow from the ensemble-mean state', 'm^3/s')
         ENDIF

         IF (DEF_hist_vars%qresv_out) THEN
            IF (p_is_worker) THEN
               IF (numresv > 0) THEN
                  WHERE (acctime_resv_mean > 0)
                     a_qresv_out_mean = a_qresv_out_mean / acctime_resv_mean
                  ELSEWHERE
                     a_qresv_out_mean = spval
                  END WHERE
               ENDIF
            ENDIF

            CALL vector_gather_and_write ( a_qresv_out_mean, numresv, totalnumresv, resv_data_address, &
               file_hist_ucat, 'qresv_out_mean', 'reservoir', itime_in_file_ucat, &
               'reservoir outflow from the ensemble-mean state', 'm^3/s')
         ENDIF

      ENDIF

      IF (p_is_worker) THEN
         IF (numucat > 0) THEN
            a_discharge_mean  (:) = 0.
            acctime_ucat_mean (:) = 0.
            a_wdsrf_ucat_mean (:) = 0.
            a_veloc_riv_mean  (:) = 0.
            a_floodarea_mean  (:) = 0.
         ENDIF
         IF (numresv > 0) THEN
            acctime_resv_mean (:) = 0.
            a_volresv_mean    (:) = 0.
            a_qresv_in_mean   (:) = 0.
            a_qresv_out_mean  (:) = 0.
         ENDIF
      ENDIF

   END SUBROUTINE hist_DA_riverlake_out
#endif

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Hist
#endif
