#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Vars_AccVars
!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Stores diagnostic sums for DA history output.
!
!   Diagnostic values are allocated by MOD_DA_Vars_TimeVars, and
!   Valid samples are summed here, the history output procedure converts
!   These sums to means.
!
! AUTHOR:
!   Zhilong Fan, 10/2026
!   Lu Li, Zhilong Fan, 11/2026: DA v1.0
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Vars_Global
   USE MOD_SPMD_Task
   USE MOD_LandPatch
   USE MOD_Namelist
#ifdef GridRiverLakeFlow
   USE MOD_Grid_RiverLakeNetwork
#endif
   USE MOD_DA_Vars_TimeVars
   IMPLICIT NONE
   SAVE

   INTERFACE acc_DA_hist
      MODULE PROCEDURE acc_DA_hist_2d
      MODULE PROCEDURE acc_DA_hist_3d
   END INTERFACE acc_DA_hist

!=============================================================================
! Diagnostic sample sums for history output
!=============================================================================
   real(r8), allocatable :: nac_obs           (  :,:)  ! Valid H(x) sample counts per observation stream and patch [-]
   real(r8), allocatable :: a_hx_ol           (  :,:)  ! Sum of open-loop H(x)
   real(r8), allocatable :: a_hx_f_ens        (:,:,:)  ! Sum of forecast H(x) for each ensemble member
   real(r8), allocatable :: a_hx_a_ens        (:,:,:)  ! Sum of analysis H(x) for each ensemble member
   real(r8), allocatable :: a_hx_f            (  :,:)  ! Sum of ensemble-mean forecast H(x)
   real(r8), allocatable :: a_hx_a            (  :,:)  ! Sum of ensemble-mean analysis H(x)
   real(r8), allocatable :: a_wliq_soisno_ol  (  :,:)  ! Sum of open-loop layer liquid water                      [kg/m2]
   real(r8), allocatable :: a_wliq_soisno_f   (  :,:)  ! Sum of ensemble-mean forecast layer liquid water         [kg/m2]
   real(r8), allocatable :: a_wliq_soisno_a   (  :,:)  ! Sum of ensemble-mean analysis layer liquid water         [kg/m2]
   real(r8), allocatable :: a_wliq_soisno_ens (:,:,:)  ! Sum of layer liquid water for each ensemble member       [kg/m2]
   real(r8), allocatable :: a_t_soisno_ol     (  :,:)  ! Sum of open-loop layer temperature                       [K]
   real(r8), allocatable :: a_t_soisno_f      (  :,:)  ! Sum of ensemble-mean forecast layer temperature          [K]
   real(r8), allocatable :: a_t_soisno_a      (  :,:)  ! Sum of ensemble-mean analysis layer temperature          [K]
#ifdef GridRiverLakeFlow
   real(r8), allocatable :: a_discharge_ol    (    :)  ! Open-loop of the  discharge integral                     [m3]
   real(r8), allocatable :: acctime_ucat_ol   (    :)  ! Open-loop of the  routing duration                       [s]
   real(r8), allocatable :: a_discharge_a     (    :)  ! Analysis discharge integral over the output interval     [m3]
   real(r8), allocatable :: acctime_ucat_a    (    :)  ! Analysis routing duration over the output interval       [s]
#endif

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE allocate_DA_AccVars ()

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            allocate (nac_obs           (nstream,                          numpatch))
            allocate (a_hx_ol           (nstream,                          numpatch))
            allocate (a_hx_f_ens        (nstream,          DEF_DA_ENS_NUM, numpatch))
            allocate (a_hx_a_ens        (nstream,          DEF_DA_ENS_NUM, numpatch))
            allocate (a_hx_f            (nstream,                          numpatch))
            allocate (a_hx_a            (nstream,                          numpatch))
            allocate (a_wliq_soisno_ol  (maxsnl+1:nl_soil,                 numpatch))
            allocate (a_wliq_soisno_f   (maxsnl+1:nl_soil,                 numpatch))
            allocate (a_wliq_soisno_a   (maxsnl+1:nl_soil,                 numpatch))
            allocate (a_wliq_soisno_ens (maxsnl+1:nl_soil, DEF_DA_ENS_NUM, numpatch))
            allocate (a_t_soisno_ol     (maxsnl+1:nl_soil,                 numpatch))
            allocate (a_t_soisno_f      (maxsnl+1:nl_soil,                 numpatch))
            allocate (a_t_soisno_a      (maxsnl+1:nl_soil,                 numpatch))
         ENDIF
#ifdef GridRiverLakeFlow
         IF (numucat > 0) THEN
            allocate (a_discharge_ol    (                                  numucat ))
            allocate (acctime_ucat_ol   (                                  numucat ))
            allocate (a_discharge_a     (                                  numucat ))
            allocate (acctime_ucat_a    (                                  numucat ))
         ENDIF
#endif
      ENDIF

   END SUBROUTINE allocate_DA_AccVars

!-----------------------------------------------------------------------------

   SUBROUTINE deallocate_DA_AccVars ()

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      IF (allocated(nac_obs          )) deallocate (nac_obs          )
      IF (allocated(a_hx_ol          )) deallocate (a_hx_ol          )
      IF (allocated(a_hx_f_ens       )) deallocate (a_hx_f_ens       )
      IF (allocated(a_hx_a_ens       )) deallocate (a_hx_a_ens       )
      IF (allocated(a_hx_f           )) deallocate (a_hx_f           )
      IF (allocated(a_hx_a           )) deallocate (a_hx_a           )
      IF (allocated(a_wliq_soisno_ol )) deallocate (a_wliq_soisno_ol )
      IF (allocated(a_wliq_soisno_f  )) deallocate (a_wliq_soisno_f  )
      IF (allocated(a_wliq_soisno_a  )) deallocate (a_wliq_soisno_a  )
      IF (allocated(a_wliq_soisno_ens)) deallocate (a_wliq_soisno_ens)
      IF (allocated(a_t_soisno_ol    )) deallocate (a_t_soisno_ol    )
      IF (allocated(a_t_soisno_f     )) deallocate (a_t_soisno_f     )
      IF (allocated(a_t_soisno_a     )) deallocate (a_t_soisno_a     )
#ifdef GridRiverLakeFlow
      IF (allocated(a_discharge_ol   )) deallocate (a_discharge_ol   )
      IF (allocated(acctime_ucat_ol  )) deallocate (acctime_ucat_ol  )
      IF (allocated(a_discharge_a    )) deallocate (a_discharge_a    )
      IF (allocated(acctime_ucat_a   )) deallocate (acctime_ucat_a   )
#endif

   END SUBROUTINE deallocate_DA_AccVars

!-----------------------------------------------------------------------------

   SUBROUTINE flush_DA_AccVars ()

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      IF (allocated(nac_obs          )) nac_obs           = 0.0_r8
      IF (allocated(a_hx_ol          )) a_hx_ol           = spval
      IF (allocated(a_hx_f_ens       )) a_hx_f_ens        = spval
      IF (allocated(a_hx_a_ens       )) a_hx_a_ens        = spval
      IF (allocated(a_hx_f           )) a_hx_f            = spval
      IF (allocated(a_hx_a           )) a_hx_a            = spval
      IF (allocated(a_wliq_soisno_ol )) a_wliq_soisno_ol  = spval
      IF (allocated(a_wliq_soisno_f  )) a_wliq_soisno_f   = spval
      IF (allocated(a_wliq_soisno_a  )) a_wliq_soisno_a   = spval
      IF (allocated(a_wliq_soisno_ens)) a_wliq_soisno_ens = spval
      IF (allocated(a_t_soisno_ol    )) a_t_soisno_ol     = spval
      IF (allocated(a_t_soisno_f     )) a_t_soisno_f      = spval
      IF (allocated(a_t_soisno_a     )) a_t_soisno_a      = spval
#ifdef GridRiverLakeFlow
      IF (allocated(a_discharge_ol   )) a_discharge_ol    = spval
      IF (allocated(acctime_ucat_ol  )) acctime_ucat_ol   = spval
      IF (allocated(a_discharge_a    )) a_discharge_a     = spval
      IF (allocated(acctime_ucat_a   )) acctime_ucat_a    = spval
#endif

   END SUBROUTINE flush_DA_AccVars

!-----------------------------------------------------------------------------

   SUBROUTINE accumulate_DA_AccVars ()

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      WHERE (hx_ol /= spval) nac_obs = nac_obs + 1.0_r8
      CALL acc_DA_hist (hx_ol,                   a_hx_ol          )
      CALL acc_DA_hist (hx_f_ens,                a_hx_f_ens       )
      CALL acc_DA_hist (hx_a_ens,                a_hx_a_ens       )
      CALL acc_DA_hist (hx_f,                    a_hx_f           )
      CALL acc_DA_hist (hx_a,                    a_hx_a           )
      CALL acc_DA_hist (wliq_soisno_ol,          a_wliq_soisno_ol )
      CALL acc_DA_hist (wliq_soisno_f,           a_wliq_soisno_f  )
      CALL acc_DA_hist (wliq_soisno_a,           a_wliq_soisno_a  )
      CALL acc_DA_hist (wliq_soisno_ens(:,1:,:), a_wliq_soisno_ens)
      CALL acc_DA_hist (t_soisno_ol,             a_t_soisno_ol    )
      CALL acc_DA_hist (t_soisno_f,              a_t_soisno_f     )
      CALL acc_DA_hist (t_soisno_a,              a_t_soisno_a     )

   END SUBROUTINE accumulate_DA_AccVars

!-----------------------------------------------------------------------------
! WARNING (duplicated codes):
!   Keep DA history definitions, accumulation, and output within DA modules
!   Without modifying the core model. Necessary utilities are included here
!   To avoid circular module dependencies.
!
!   Remove these copies once the core utilities are moved to an independent
!   Module that DA can use directly.
!-----------------------------------------------------------------------------

   SUBROUTINE acc_DA_hist_2d (var, s)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      real(r8), intent(in)    :: var (  :,:)
      real(r8), intent(inout) :: s   (  :,:)
      integer :: i1, i2

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

      real(r8), intent(in)    :: var (:,:,:)
      real(r8), intent(inout) :: s   (:,:,:)
      integer :: i1, i2, i3

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
END MODULE MOD_DA_Vars_AccVars
#endif