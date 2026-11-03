#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Vars_Utils
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Provides array accumulation and gridded history output utilities for DA.
!    History output includes time normalization, patch-to-grid mapping, area
!    normalization, and writing the resulting fields.
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Vars_Global
   USE MOD_SPMD_Task
   USE MOD_DataType
   USE MOD_Block
   USE MOD_Namelist
   USE MOD_HistGridded
   USE MOD_Vars_1DAccFluxes
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: acc_DA_hist
   PUBLIC :: write_DA_hist

   INTERFACE acc_DA_hist
      MODULE PROCEDURE acc_DA_hist_2d
      MODULE PROCEDURE acc_DA_hist_3d
   END INTERFACE acc_DA_hist

   INTERFACE write_DA_hist
      MODULE PROCEDURE write_DA_hist_2d
      MODULE PROCEDURE write_DA_hist_3d
   END INTERFACE write_DA_hist

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE acc_DA_hist_2d (var, s)

!-----------------------------------------------------------------------------
      IMPLICIT NONE
      real(r8), intent(in)                   :: var(:,:)
      real(r8), intent(inout)                :: s  (:,:)
      integer                                :: i1
      integer                                :: i2

      DO i2 = lbound(var,2), ubound(var,2)
         DO i1 = lbound(var,1), ubound(var,1)
            IF (var(i1,i2) /= spval) THEN
               IF (s(i1,i2) /= spval) THEN
                  s(i1,i2) = s(i1,i2) + var(i1,i2)
               ELSE
                  s(i1,i2) = var(i1,i2)
               ENDIF
            ENDIF
         ENDDO
      ENDDO

   END SUBROUTINE acc_DA_hist_2d

!-----------------------------------------------------------------------------

   SUBROUTINE acc_DA_hist_3d (var, s)

!-----------------------------------------------------------------------------
      IMPLICIT NONE
      real(r8), intent(in)                   :: var(:,:,:)
      real(r8), intent(inout)                :: s  (:,:,:)
      integer                                :: i1
      integer                                :: i2
      integer                                :: i3

      DO i3 = lbound(var,3), ubound(var,3)
         DO i2 = lbound(var,2), ubound(var,2)
            DO i1 = lbound(var,1), ubound(var,1)
               IF (var(i1,i2,i3) /= spval) THEN
                  IF (s(i1,i2,i3) /= spval) THEN
                     s(i1,i2,i3) = s(i1,i2,i3) + var(i1,i2,i3)
                  ELSE
                     s(i1,i2,i3) = var(i1,i2,i3)
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO

   END SUBROUTINE acc_DA_hist_3d

!-----------------------------------------------------------------------------

   SUBROUTINE write_DA_hist_2d (HistForm, is_hist, &
         acc_vec, file_hist, varname, itime_in_file, dim1name, lb1, ndim1, &
         sumarea, filter, longname, units, filter_3d)

!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Normalize diagnostic sums, map patches to the history grid, and write
!    the resulting field. Array rank excludes the history-grid dimensions.
!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      character(len=*), intent(in)           :: HistForm
      logical, intent(in)                    :: is_hist
      real(r8), intent(inout)                :: acc_vec(:,:)
      character(len=*), intent(in)           :: file_hist
      character(len=*), intent(in)           :: varname
      integer, intent(in)                    :: itime_in_file
      character(len=*), intent(in)           :: dim1name
      integer, intent(in)                    :: lb1
      integer, intent(in)                    :: ndim1
      type(block_data_real8_2d), intent(in)  :: sumarea
      logical, intent(in)                    :: filter(:)
      character(len=*), intent(in)           :: longname
      character(len=*), intent(in)           :: units
      logical, intent(in), optional          :: filter_3d(:,:)

!------------------------ Local Variables ------------------------------------
      type(block_data_real8_3d)              :: flux_xy_3d
      type(block_data_real8_3d)              :: sumarea_3d
      type(block_data_real8_2d)              :: sumarea_one_source
      integer                                :: iblkme
      integer                                :: xblk
      integer                                :: yblk
      integer                                :: xloc
      integer                                :: yloc
      integer                                :: i1
      integer                                :: compress
      real(r8)                               :: area_den

!-----------------------------------------------------------------------------
      IF (.not. is_hist) RETURN
      IF (HistForm /= 'Gridded') RETURN

      IF (p_is_worker)  THEN
         WHERE (acc_vec /= spval)  acc_vec = acc_vec / nac
      ENDIF
      IF (p_is_io) THEN
         CALL allocate_block_data (ghist, flux_xy_3d, ndim1, lb1)
         IF (present(filter_3d)) THEN
            CALL allocate_block_data (ghist, sumarea_3d, ndim1, lb1)
            CALL allocate_block_data (ghist, sumarea_one_source)
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
! DESCRIPTION:
!    Normalize diagnostic sums, map patches to the history grid, and write
!    the resulting field. Array rank excludes the history-grid dimensions.
!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      character(len=*), intent(in)           :: HistForm
      logical, intent(in)                    :: is_hist
      real(r8), intent(inout)                :: acc_vec(:,:,:)
      character(len=*), intent(in)           :: file_hist
      character(len=*), intent(in)           :: varname
      integer, intent(in)                    :: itime_in_file
      character(len=*), intent(in)           :: dim1name
      integer, intent(in)                    :: lb1
      integer, intent(in)                    :: ndim1
      character(len=*), intent(in)           :: dim2name
      integer, intent(in)                    :: lb2
      integer, intent(in)                    :: ndim2
      type(block_data_real8_2d), intent(in)  :: sumarea
      logical, intent(in)                    :: filter(:)
      character(len=*), intent(in)           :: longname
      character(len=*), intent(in)           :: units
      real(r8), intent(in), optional         :: acc_num(:)

!------------------------ Local Variables ------------------------------------
      type(block_data_real8_4d)              :: flux_xy_4d
      integer                                :: iblkme
      integer                                :: xblk
      integer                                :: yblk
      integer                                :: xloc
      integer                                :: yloc
      integer                                :: i1
      integer                                :: i2
      integer                                :: compress

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
END MODULE MOD_DA_Vars_Utils
#endif
