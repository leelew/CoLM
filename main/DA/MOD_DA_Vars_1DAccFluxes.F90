#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Vars_1DAccFluxes
!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Accumulators for DA diagnostic variables
!
! AUTHOR:
!   Zhilong Fan, 10/2026
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Vars_Global
   USE MOD_SPMD_Task
   USE MOD_LandPatch
   USE MOD_Namelist, only: DEF_DA_ENS_NUM, DEF_DA_OBS_TARGET
   USE MOD_DA_Vars_TimeVariables
   IMPLICIT NONE
   SAVE

   PUBLIC :: allocate_1D_DAAccFluxes
   PUBLIC :: deallocate_1D_DAAccFluxes
   PUBLIC :: flush_1D_DAAccFluxes
   PUBLIC :: accumulate_1D_DAAccFluxes

   real(r8), allocatable :: a_hx_ol          (:,:  ) ! open-loop H(x)                                 [source, patch]
   real(r8), allocatable :: a_hx_f_ens       (:,:,:) ! prior H(x)                                     [source, ensemble, patch]
   real(r8), allocatable :: a_hx_a_ens       (:,:,:) ! posterior H(x)                                 [source, ensemble, patch]
   real(r8), allocatable :: a_hx_f           (:,:  ) ! ensemble-mean prior H(x)                       [source, patch]
   real(r8), allocatable :: a_hx_a           (:,:  ) ! ensemble-mean posterior H(x)                   [source, patch]
   real(r8), allocatable :: a_wliq_soisno_ol (:,:  ) ! open-loop liquid water in layers               [kg/m2]
   real(r8), allocatable :: a_wliq_soisno_f  (:,:  ) ! ensemble-mean prior liquid water in layers     [kg/m2]
   real(r8), allocatable :: a_wliq_soisno_a  (:,:  ) ! ensemble-mean posterior liquid water in layers [kg/m2]
   real(r8), allocatable :: a_wliq_soisno_ens(:,:,:) ! ensemble liquid water in layers                [kg/m2]
   real(r8), allocatable :: a_t_soisno_ol    (:,:  ) ! open-loop temperature in layers                [K]
   real(r8), allocatable :: a_t_soisno_f     (:,:  ) ! ensemble-mean prior temperature in layers      [K]
   real(r8), allocatable :: a_t_soisno_a     (:,:  ) ! ensemble-mean posterior temperature in layers  [K]

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE allocate_1D_DAAccFluxes ()

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      integer :: i

!-----------------------------------------------------------------------------
      nsource = 0
      DO i = 1, size(DEF_DA_OBS_TARGET)
         IF (trim(DEF_DA_OBS_TARGET(i)) == 'SM' .or. &
             trim(DEF_DA_OBS_TARGET(i)) == 'ST') nsource = nsource + 1
      ENDDO

      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            IF (nsource > 0) THEN
               allocate (a_hx_ol   (nsource,               numpatch))
               allocate (a_hx_f_ens(nsource,DEF_DA_ENS_NUM,numpatch))
               allocate (a_hx_a_ens(nsource,DEF_DA_ENS_NUM,numpatch))
               allocate (a_hx_f    (nsource,               numpatch))
               allocate (a_hx_a    (nsource,               numpatch))
            ENDIF
            allocate (a_wliq_soisno_ol (maxsnl+1:nl_soil,numpatch))
            allocate (a_wliq_soisno_f  (maxsnl+1:nl_soil,numpatch))
            allocate (a_wliq_soisno_a  (maxsnl+1:nl_soil,numpatch))
            allocate (a_wliq_soisno_ens(maxsnl+1:nl_soil,DEF_DA_ENS_NUM,numpatch))
            allocate (a_t_soisno_ol    (maxsnl+1:nl_soil,numpatch))
            allocate (a_t_soisno_f     (maxsnl+1:nl_soil,numpatch))
            allocate (a_t_soisno_a     (maxsnl+1:nl_soil,numpatch))
         ENDIF
      ENDIF

   END SUBROUTINE allocate_1D_DAAccFluxes

!-----------------------------------------------------------------------------

   SUBROUTINE deallocate_1D_DAAccFluxes ()

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!-----------------------------------------------------------------------------
      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            IF (allocated(a_hx_ol))    deallocate (a_hx_ol   )
            IF (allocated(a_hx_f_ens)) deallocate (a_hx_f_ens)
            IF (allocated(a_hx_a_ens)) deallocate (a_hx_a_ens)
            IF (allocated(a_hx_f))     deallocate (a_hx_f    )
            IF (allocated(a_hx_a))     deallocate (a_hx_a    )
            deallocate (a_wliq_soisno_ol )
            deallocate (a_wliq_soisno_f  )
            deallocate (a_wliq_soisno_a  )
            deallocate (a_wliq_soisno_ens)
            deallocate (a_t_soisno_ol    )
            deallocate (a_t_soisno_f     )
            deallocate (a_t_soisno_a     )
         ENDIF
      ENDIF

   END SUBROUTINE deallocate_1D_DAAccFluxes

!-----------------------------------------------------------------------------

   SUBROUTINE flush_1D_DAAccFluxes ()

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!-----------------------------------------------------------------------------
      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            IF (allocated(a_hx_ol))    a_hx_ol   (:,:)   = spval
            IF (allocated(a_hx_f_ens)) a_hx_f_ens(:,:,:) = spval
            IF (allocated(a_hx_a_ens)) a_hx_a_ens(:,:,:) = spval
            IF (allocated(a_hx_f))     a_hx_f    (:,:)   = spval
            IF (allocated(a_hx_a))     a_hx_a    (:,:)   = spval
            a_wliq_soisno_ol (:,:)   = spval
            a_wliq_soisno_f  (:,:)   = spval
            a_wliq_soisno_a  (:,:)   = spval
            a_wliq_soisno_ens(:,:,:) = spval
            a_t_soisno_ol    (:,:)   = spval
            a_t_soisno_f     (:,:)   = spval
            a_t_soisno_a     (:,:)   = spval
         ENDIF
      ENDIF

   END SUBROUTINE flush_1D_DAAccFluxes

!-----------------------------------------------------------------------------

   SUBROUTINE accumulate_1D_DAAccFluxes ()

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!-----------------------------------------------------------------------------
      IF (allocated(hx_f_ens) .and. allocated(a_hx_f_ens)) THEN
         CALL acc2d_DA (hx_ol   , a_hx_ol   )
         CALL acc3d_DA (hx_f_ens, a_hx_f_ens)
         CALL acc3d_DA (hx_a_ens, a_hx_a_ens)
         CALL acc2d_DA (hx_f    , a_hx_f    )
         CALL acc2d_DA (hx_a    , a_hx_a    )
      ENDIF
      CALL acc2d_DA (wliq_soisno_ol, a_wliq_soisno_ol)
      CALL acc2d_DA (wliq_soisno_f , a_wliq_soisno_f )
      CALL acc2d_DA (wliq_soisno_a , a_wliq_soisno_a )
      CALL acc3d_DA (wliq_soisno_ens(:,1:DEF_DA_ENS_NUM,:), a_wliq_soisno_ens)
      CALL acc2d_DA (t_soisno_ol   , a_t_soisno_ol   )
      CALL acc2d_DA (t_soisno_f    , a_t_soisno_f    )
      CALL acc2d_DA (t_soisno_a    , a_t_soisno_a    )

   END SUBROUTINE accumulate_1D_DAAccFluxes

!-----------------------------------------------------------------------------

   SUBROUTINE acc2d_DA (var, s)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      real(r8), intent(in)    :: var(:,:)
      real(r8), intent(inout) :: s  (:,:)

      integer :: i1, i2

!-----------------------------------------------------------------------------
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

   END SUBROUTINE acc2d_DA

!-----------------------------------------------------------------------------

   SUBROUTINE acc3d_DA (var, s)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      real(r8), intent(in)    :: var(:,:,:)
      real(r8), intent(inout) :: s  (:,:,:)

      integer :: i1, i2, i3

!-----------------------------------------------------------------------------
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

   END SUBROUTINE acc3d_DA

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Vars_1DAccFluxes
#endif
