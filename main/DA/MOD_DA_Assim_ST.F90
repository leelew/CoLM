#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Assim_ST
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Soil-temperature assimilation driver using observations prepared by Obs_ST.
!
! WORKFLOW:
!    1. Initialize observation streams from the namelist.
!    2. Prepare observations and forecast H(x); save forecast diagnostics.
!    3. Apply a localized LETKF update and temperature constraints per patch.
!    4. Save analysis diagnostics, recompute H(x), and clear time-step data.
!
! NOTES:
!    The LETKF update is applied to the top two soil layers, with temperature
!    increments limited to +/-20 K. Source modules provide the H(x) operators.
!    Model feedback is handled by the DA ensemble postprocessing step.
!
! HISTORY:
!    Lu Li, 05/2026: First implementation
!    Lu Li, Zhilong Fan, 11/2026: DA v1.0
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Vars_TimeInvariants
   USE MOD_Vars_Global
   USE MOD_LandPatch
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_DA_Obs_ST
   USE MOD_DA_EnKF
   USE MOD_DA_Vars_TimeVars
   IMPLICIT NONE

   type(ST_type) :: st

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE init_DA_ST ()

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      integer :: is, i, nsensor
      type(ST_config_type), allocatable :: configs(:)

!-----------------------------------------------------------------------------
      nsensor = 0
      DO i = 1, size(DEF_DA_OBS_TARGET)
         IF (trim(DEF_DA_OBS_TARGET(i)) == 'ST') THEN
            nsensor = nsensor + 1
         ENDIF
      ENDDO
      IF (nsensor == 0) RETURN

      allocate(configs(nsensor))
      is = 0
      DO i = 1, size(DEF_DA_OBS_TARGET)
         IF (trim(DEF_DA_OBS_TARGET(i)) == 'ST') THEN
            is = is + 1
            configs(is)%source_name = DEF_DA_OBS_SOURCE(i)
            configs(is)%sensor_name = DEF_DA_OBS_SENSOR(i)
            configs(is)%var_name    = DEF_DA_OBS_VAR(i)
            configs(is)%fghz        = DEF_DA_OBS_FGHZ(i)
            configs(is)%istream     = i
         ENDIF
      ENDDO

      CALL st%init(landpatch, configs)

      deallocate(configs)

   END SUBROUTINE init_DA_ST

!-----------------------------------------------------------------------------

   SUBROUTINE run_DA_ST (idate, deltim)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      integer,  intent(in) :: idate(3)
      real(r8), intent(in) :: deltim

!------------------------ Local Variables ------------------------------------
      integer :: np
      integer :: iens
      integer :: nobs_p                      ! Number of valid observations selected for the target patch
      real(r8), allocatable :: y_p(:)        ! Selected observation values [nobs_p]
      real(r8), allocatable :: hx_p(:,:)     ! Forecast observation equivalents [nobs_p, ens]
      real(r8), allocatable :: r_p(:)        ! Localized observation error variances [nobs_p]
      real(r8), allocatable :: trans(:,:)    ! LETKF ensemble transform [ens, ens]
      real(r8), allocatable :: state_f(:,:)  ! Forecast soil temperature [soil, ens], K
      real(r8), allocatable :: state_a(:,:)  ! Transformed soil temperature [soil, ens], K

!-----------------------------------------------------------------------------
      ! Prepare observations and forecast H(x).
      CALL st%calc_obs(idate, deltim)
      CALL st%save(.false.)

      IF (p_is_worker .and. numpatch > 0) THEN
         t_soisno_ol(:,:) = t_soisno_ens(:,0,:)
         IF (DEF_DA_ENS_NUM > 1) THEN
            t_soisno_f(:,:) = sum(t_soisno_ens(:,1:DEF_DA_ENS_NUM,:), dim=2) / DEF_DA_ENS_NUM
         ENDIF
      ENDIF

      ! Apply the localized LETKF update to each eligible patch.
      IF (p_is_worker .and. DEF_DA_ENS_NUM > 1) THEN
         DO np = 1, numpatch
            IF (patchtype(np) >= 3) CYCLE

            ! Select observations and localize their error variances.
            CALL st%prepare(np, nobs_p, y_p, hx_p, r_p)
            IF (nobs_p == 0) CYCLE

            ! Compute the ensemble transform.
            allocate(trans(DEF_DA_ENS_NUM, DEF_DA_ENS_NUM))
            CALL letkf(DEF_DA_ENS_NUM, nobs_p, hx_p, y_p, r_p, st%infl, trans)

            allocate(state_f(nl_soil, DEF_DA_ENS_NUM))
            allocate(state_a(nl_soil, DEF_DA_ENS_NUM))
            DO iens = 1, DEF_DA_ENS_NUM
               state_f(:,iens) = t_soisno_ens(1:,iens,np)
            ENDDO
            CALL dgemm('N', 'N', nl_soil, DEF_DA_ENS_NUM, DEF_DA_ENS_NUM, &
                        1.0_r8, state_f, nl_soil, trans, DEF_DA_ENS_NUM, &
                        0.0_r8, state_a, nl_soil)

            ! limit the temperature increment to keep the analysis physical
            DO iens = 1, DEF_DA_ENS_NUM
               t_soisno_ens(1:2,iens,np) = t_soisno_ens(1:2,iens,np) + &
                  max(-20.0_r8, min(20.0_r8, state_a(1:2,iens) - t_soisno_ens(1:2,iens,np)))
            ENDDO

            deallocate(y_p, r_p, hx_p, trans)
            deallocate(state_f, state_a)
         ENDDO
      ENDIF

      ! Save analysis diagnostics and recompute model-space H(x).
      IF (DEF_DA_ENS_NUM > 1) THEN
         IF (p_is_worker .and. numpatch > 0) THEN
            t_soisno_a(:,:) = sum(t_soisno_ens(:,1:DEF_DA_ENS_NUM,:), dim=2) / DEF_DA_ENS_NUM
         ENDIF
         CALL st%calc_model()
         CALL st%save(.true.)
      ENDIF
      CALL st%clear()

   END SUBROUTINE run_DA_ST

!-----------------------------------------------------------------------------

   SUBROUTINE end_DA_ST()

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!-----------------------------------------------------------------------------
      CALL st%clear()
      IF (allocated(st%cfgs)) deallocate(st%cfgs)
      IF (allocated(st%pm%sensors)) deallocate(st%pm%sensors)
      IF (allocated(st%pm%cfgs)) deallocate(st%pm%cfgs)

      st%pm%nsensor = 0
      st%use_pm = .false.

   END SUBROUTINE end_DA_ST

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Assim_ST
#endif
