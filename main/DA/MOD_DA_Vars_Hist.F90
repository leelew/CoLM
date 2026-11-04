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
   USE MOD_Block
   USE MOD_HistGridded
   USE MOD_Vars_TimeInvariants
   USE MOD_Forcing
   USE MOD_Vars_1DAccFluxes, only: nac
#ifdef GridRiverLakeFlow
   USE MOD_Vector_ReadWrite
   USE MOD_Grid_RiverLakeNetwork
#endif
   USE MOD_DA_Vars_TimeVars
   USE MOD_DA_Vars_AccVars
   IMPLICIT NONE
   SAVE

   INTERFACE write_DA_hist
      MODULE PROCEDURE write_DA_hist_2d
      MODULE PROCEDURE write_DA_hist_3d
   END INTERFACE write_DA_hist

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE save_DA_hist (file_hist, HistForm, itime_in_file, idate, sumarea, filter)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      character(len=*),         intent(in)    :: file_hist
      character(len=*),         intent(in)    :: HistForm
      integer,                  intent(in)    :: itime_in_file
      integer,                  intent(in)    :: idate(3)
      type(block_data_real8_2d), intent(inout) :: sumarea
      logical,                  intent(inout) :: filter(:)

!------------------------ Local Variables ------------------------------------
      logical, allocatable :: filter_hx(:,:)
#ifdef GridRiverLakeFlow
      character(len=256) :: file_hist_ucat
      integer :: itime_in_file_ucat
#endif
      integer :: i

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

      ! Write soil and snow diagnostics.
      CALL write_DA_hist ( &
         HistForm, &
         DEF_hist_vars%DA_wliq_soisno_ol, &
         a_wliq_soisno_ol, &
         file_hist, &
         'f_wliq_soisno_ol', &
         itime_in_file, &
         'soilsnow', &
         maxsnl + 1, &
         nl_soil - maxsnl, &
         sumarea, &
         filter, &
         'Open-loop liquid water in soil and snow layers', &
         'kg/m2')
      CALL write_DA_hist ( &
         HistForm, &
         DEF_hist_vars%DA_wliq_soisno_f, &
         a_wliq_soisno_f, &
         file_hist, &
         'f_wliq_soisno_f', &
         itime_in_file, &
         'soilsnow', &
         maxsnl + 1, &
         nl_soil - maxsnl, &
         sumarea, &
         filter, &
         'Forecast ensemble-mean liquid water in soil and snow layers', &
         'kg/m2')
      CALL write_DA_hist ( &
         HistForm, &
         DEF_hist_vars%DA_wliq_soisno_a, &
         a_wliq_soisno_a, &
         file_hist, &
         'f_wliq_soisno_a', &
         itime_in_file, &
         'soilsnow', &
         maxsnl + 1, &
         nl_soil - maxsnl, &
         sumarea, &
         filter, &
         'Analysis ensemble-mean liquid water in soil and snow layers', &
         'kg/m2')
      CALL write_DA_hist ( &
         HistForm, &
         DEF_hist_vars%DA_t_soisno_ol, &
         a_t_soisno_ol, &
         file_hist, &
         'f_t_soisno_ol', &
         itime_in_file, &
         'soilsnow', &
         maxsnl + 1, &
         nl_soil - maxsnl, &
         sumarea, &
         filter, &
         'Open-loop temperature in soil and snow layers', &
         'K')
      CALL write_DA_hist ( &
         HistForm, &
         DEF_hist_vars%DA_t_soisno_f, &
         a_t_soisno_f, &
         file_hist, &
         'f_t_soisno_f', &
         itime_in_file, &
         'soilsnow', &
         maxsnl + 1, &
         nl_soil - maxsnl, &
         sumarea, &
         filter, &
         'Forecast ensemble-mean temperature in soil and snow layers', &
         'K')
      CALL write_DA_hist ( &
         HistForm, &
         DEF_hist_vars%DA_t_soisno_a, &
         a_t_soisno_a, &
         file_hist, &
         'f_t_soisno_a', &
         itime_in_file, &
         'soilsnow', &
         maxsnl + 1, &
         nl_soil - maxsnl, &
         sumarea, &
         filter, &
         'Analysis ensemble-mean temperature in soil and snow layers', &
         'K')
      CALL write_DA_hist ( &
         HistForm, &
         DEF_hist_vars%DA_wliq_soisno_ens, &
         a_wliq_soisno_ens, &
         file_hist, &
         'f_wliq_soisno_ens', &
         itime_in_file, &
         'soilsnow', &
         maxsnl + 1, &
         nl_soil - maxsnl, &
         'ens', &
         1, &
         DEF_DA_ENS_NUM, &
         sumarea, &
         filter, &
         'ensemble liquid water in soil layers', &
         'kg/m2')

      ! Write H(x) diagnostics by observation source.
      IF (nsource > 0) THEN
         allocate (filter_hx(nsource,numpatch))

         IF (p_is_worker) THEN
            DO i = 1, nsource
               filter_hx(i,:) = filter .and. a_hx_ol(i,:) /= spval
            ENDDO
         ENDIF
         CALL write_DA_hist ( &
            HistForm, &
            DEF_hist_vars%DA_hx_ol, &
            a_hx_ol, &
            file_hist, &
            'f_hx_ol', &
            itime_in_file, &
            'source', &
            1, &
            nsource, &
            sumarea, &
            filter, &
            'Open-loop H(x) by DA observation source', &
            'source-dependent', &
            filter_3d = filter_hx)
         CALL write_DA_hist ( &
            HistForm, &
            DEF_hist_vars%DA_hx_f, &
            a_hx_f, &
            file_hist, &
            'f_hx_f', &
            itime_in_file, &
            'source', &
            1, &
            nsource, &
            sumarea, &
            filter, &
            'Prior H(x) by DA observation source', &
            'source-dependent', &
            filter_3d = filter_hx)
         CALL write_DA_hist ( &
            HistForm, &
            DEF_hist_vars%DA_hx_a, &
            a_hx_a, &
            file_hist, &
            'f_hx_a', &
            itime_in_file, &
            'source', &
            1, &
            nsource, &
            sumarea, &
            filter, &
            'Analysis H(x) by DA observation source', &
            'source-dependent', &
            filter_3d = filter_hx)
         CALL write_DA_hist ( &
            HistForm, &
            DEF_hist_vars%DA_hx_f_ens, &
            a_hx_f_ens, &
            file_hist, &
            'f_hx_f_ens', &
            itime_in_file, &
            'source', &
            1, &
            nsource, &
            'ens', &
            1, &
            DEF_DA_ENS_NUM, &
            sumarea, &
            filter, &
            'Prior ensemble H(x) by DA observation source', &
            'source-dependent')
         CALL write_DA_hist ( &
            HistForm, &
            DEF_hist_vars%DA_hx_a_ens, &
            a_hx_a_ens, &
            file_hist, &
            'f_hx_a_ens', &
            itime_in_file, &
            'source', &
            1, &
            nsource, &
            'ens', &
            1, &
            DEF_DA_ENS_NUM, &
            sumarea, &
            filter, &
            'Analysis ensemble H(x) by DA observation source', &
            'source-dependent')

         deallocate (filter_hx)
      ENDIF

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
! WARNING (duplicated codes):
!   Keep DA history definitions, accumulation, and output within DA modules
!   without modifying the core model. Necessary utilities are included here
!   to avoid circular module dependencies.
!
!   Remove these copies once the core utilities are moved to an independent
!   module that DA can use directly.
!-----------------------------------------------------------------------------

   SUBROUTINE write_DA_hist_2d (HistForm, is_hist, &
         acc_vec, file_hist, varname, itime_in_file, dim1name, lb1, ndim1, &
         sumarea, filter, longname, units, filter_3d)

!-----------------------------------------------------------------------------

      IMPLICIT NONE

      character(len=*), intent(in)             :: HistForm
      logical, intent(in)                      :: is_hist
      real(r8), intent(inout)                  :: acc_vec(:,:)
      character(len=*), intent(in)             :: file_hist
      character(len=*), intent(in)             :: varname
      integer, intent(in)                      :: itime_in_file
      character(len=*), intent(in)             :: dim1name
      integer, intent(in)                      :: lb1
      integer, intent(in)                      :: ndim1
      type(block_data_real8_2d), intent(in)    :: sumarea
      logical, intent(in)                      :: filter(:)
      character(len=*), intent(in)             :: longname
      character(len=*), intent(in)             :: units
      logical, intent(in), optional            :: filter_3d(:,:)

      type(block_data_real8_3d)                :: flux_xy_3d
      type(block_data_real8_3d)                :: sumarea_3d
      type(block_data_real8_2d)                :: sumarea_one_source
      integer                                  :: iblkme
      integer                                  :: xblk
      integer                                  :: yblk
      integer                                  :: xloc
      integer                                  :: yloc
      integer                                  :: i1
      integer                                  :: compress
      real(r8)                                 :: area_den

!-----------------------------------------------------------------------------
      IF (.not. is_hist) RETURN
      IF (HistForm /= 'Gridded') RETURN

      IF (p_is_worker) THEN
         WHERE (acc_vec /= spval) acc_vec = acc_vec / nac
      ENDIF
      IF (p_is_io) THEN
         CALL allocate_block_data (ghist, flux_xy_3d, ndim1, lb1                        )
         IF (present(filter_3d)) THEN
            CALL allocate_block_data (ghist, sumarea_3d, ndim1, lb1                     )
            CALL allocate_block_data (ghist, sumarea_one_source                         )
         ENDIF
      ENDIF

      IF (present(filter_3d)) THEN
         DO i1 = lb1, lb1 + ndim1 - 1
            CALL mp2g_hist%get_sumarea(sumarea_one_source, filter_3d(i1,:))
            IF (p_is_io) THEN
               DO iblkme = 1, gblock%nblkme
                  xblk = gblock%xblkme(iblkme)
                  yblk = gblock%yblkme(iblkme)
                  sumarea_3d%blk(xblk,yblk)%val(i1,:,:) = &
                     sumarea_one_source%blk(xblk,yblk)%val(:,:)
               ENDDO
            ENDIF
         ENDDO
      ENDIF

      CALL mp2g_hist%pset2grid (acc_vec, flux_xy_3d, spv = spval, msk = filter)

      IF (p_is_io) THEN
         DO iblkme = 1, gblock%nblkme
            xblk = gblock%xblkme(iblkme)
            yblk = gblock%yblkme(iblkme)

            DO yloc = 1, ghist%ycnt(yblk)
               DO xloc = 1, ghist%xcnt(xblk)

                  IF (sumarea%blk(xblk,yblk)%val(xloc,yloc) > 0.00001) THEN
                     DO i1 = flux_xy_3d%lb1, flux_xy_3d%ub1
                        IF (flux_xy_3d%blk(xblk,yblk)%val(i1,xloc,yloc) /= spval) THEN
                           IF (present(filter_3d)) THEN
                              area_den = sumarea_3d%blk(xblk,yblk)%val(i1,xloc,yloc)
                           ELSE
                              area_den = sumarea%blk(xblk,yblk)%val(xloc,yloc)
                           ENDIF

                           IF (area_den > 0.00001_r8) THEN
                              flux_xy_3d%blk(xblk,yblk)%val(i1,xloc,yloc) &
                                 = flux_xy_3d%blk(xblk,yblk)%val(i1,xloc,yloc) &
                                 / area_den
                           ELSE
                              flux_xy_3d%blk(xblk,yblk)%val(i1,xloc,yloc) = spval
                           ENDIF
                        ENDIF
                     ENDDO
                  ELSE
                     flux_xy_3d%blk(xblk,yblk)%val(:,xloc,yloc) = spval
                  ENDIF

               ENDDO
            ENDDO

         ENDDO
      ENDIF

      compress = DEF_HIST_CompressLevel
      CALL hist_write_var_real8_3d (file_hist, varname, dim1name, ghist, &
         itime_in_file, flux_xy_3d, compress, longname, units)

   END SUBROUTINE write_DA_hist_2d

!-----------------------------------------------------------------------------

   SUBROUTINE write_DA_hist_3d (HistForm, is_hist, &
         acc_vec, file_hist, varname, itime_in_file, &
         dim1name, lb1, ndim1, dim2name, lb2, ndim2, &
         sumarea, filter, longname, units, acc_num)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      character(len=*), intent(in)             :: HistForm
      logical, intent(in)                      :: is_hist
      real(r8), intent(inout)                  :: acc_vec(:,:,:)
      character(len=*), intent(in)             :: file_hist
      character(len=*), intent(in)             :: varname
      integer, intent(in)                      :: itime_in_file
      character(len=*), intent(in)             :: dim1name
      integer, intent(in)                      :: lb1
      integer, intent(in)                      :: ndim1
      character(len=*), intent(in)             :: dim2name
      integer, intent(in)                      :: lb2
      integer, intent(in)                      :: ndim2
      type(block_data_real8_2d), intent(in)    :: sumarea
      logical, intent(in)                      :: filter(:)
      character(len=*), intent(in)             :: longname
      character(len=*), intent(in)             :: units
      real(r8), intent(in), optional           :: acc_num(:)

      type(block_data_real8_4d)                :: flux_xy_4d
      integer                                  :: iblkme
      integer                                  :: xblk
      integer                                  :: yblk
      integer                                  :: xloc
      integer                                  :: yloc
      integer                                  :: i1
      integer                                  :: i2
      integer                                  :: compress

!-----------------------------------------------------------------------------
      IF (.not. is_hist) RETURN
      IF (HistForm /= 'Gridded') RETURN

      IF (.not. present(acc_num)) THEN
         IF (p_is_worker) THEN
            WHERE (acc_vec /= spval) acc_vec = acc_vec / nac
         ENDIF
      ELSE
         IF (p_is_worker) THEN
            DO i1 = lbound(acc_vec,1), ubound(acc_vec,1)
               DO i2 = lbound(acc_vec,2), ubound(acc_vec,2)
                  WHERE (acc_vec(i1,i2,:) /= spval .and. acc_num > 0)
                     acc_vec(i1,i2,:) = acc_vec(i1,i2,:) / acc_num
                  END WHERE
               ENDDO
            ENDDO
         ENDIF
      ENDIF

      IF (p_is_io) THEN
         CALL allocate_block_data (ghist, flux_xy_4d, ndim1, ndim2, lb1 = lb1, lb2 = lb2)
      ENDIF

      CALL mp2g_hist%pset2grid (acc_vec, flux_xy_4d, spv = spval, msk = filter)

      IF (p_is_io) THEN
         DO iblkme = 1, gblock%nblkme
            xblk = gblock%xblkme(iblkme)
            yblk = gblock%yblkme(iblkme)

            DO yloc = 1, ghist%ycnt(yblk)
               DO xloc = 1, ghist%xcnt(xblk)

                  IF (sumarea%blk(xblk,yblk)%val(xloc,yloc) > 0.00001) THEN
                     DO i1 = flux_xy_4d%lb1, flux_xy_4d%ub1
                        DO i2 = flux_xy_4d%lb2, flux_xy_4d%ub2
                           IF (flux_xy_4d%blk(xblk,yblk)%val(i1,i2,xloc,yloc) /= spval) THEN
                              flux_xy_4d%blk(xblk,yblk)%val(i1,i2,xloc,yloc) &
                                 = flux_xy_4d%blk(xblk,yblk)%val(i1,i2,xloc,yloc) &
                                 / sumarea%blk(xblk,yblk)%val(xloc,yloc)
                           ENDIF
                        ENDDO
                     ENDDO
                  ELSE
                     flux_xy_4d%blk(xblk,yblk)%val(:,:,xloc,yloc) = spval
                  ENDIF

               ENDDO
            ENDDO

         ENDDO
      ENDIF

      compress = DEF_HIST_CompressLevel
      CALL hist_write_var_real8_4d (file_hist, varname, dim1name, dim2name, &
         ghist, itime_in_file, flux_xy_4d, compress, longname, units)

   END SUBROUTINE write_DA_hist_3d

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Vars_Hist
#endif
