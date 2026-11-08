#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Assim_SM
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Soil-moisture assimilation driver using observations prepared by Obs_SM,
!    including passive-microwave data and ANATQ pseudo-observations.
!
! WORKFLOW:
!    1. Initialize observation streams from the namelist.
!    2. Prepare observations and forecast H(x); save forecast diagnostics.
!    3. Apply a localized LETKF update and soil-water constraints per patch.
!    4. Save analysis diagnostics, recompute H(x), and clear time-step data.
!
! NOTES:
!    The LETKF update is applied to the top two soil layers. Water-content
!    bounds are enforced and the aquifer compensates liquid-water changes.
!    Model feedback is handled by the DA ensemble postprocessing step.
!
! HISTORY:
!    Lu Li, Zhilong Fan, 11/2026: DA v1.0
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Vars_TimeInvariants
   USE MOD_Vars_Global
   USE MOD_LandPatch
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_Const_Physical
   USE MOD_DA_Obs_SM
   USE MOD_DA_EnKF
   USE MOD_DA_Vars_TimeVars
   IMPLICIT NONE

   type(SM_type) :: sm

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE init_DA_SM ()

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      integer :: is, i, nsensor
      type(SM_config_type), allocatable :: configs(:)

!-----------------------------------------------------------------------------
      nsensor = 0
      DO i = 1, size(DEF_DA_OBS_TARGET)
         IF (trim(DEF_DA_OBS_TARGET(i)) == 'SM') THEN
            nsensor = nsensor + 1
         ENDIF
      ENDDO
      IF (nsensor == 0) RETURN

      allocate(configs(nsensor))
      is = 0
      DO i = 1, size(DEF_DA_OBS_TARGET)
         IF (trim(DEF_DA_OBS_TARGET(i)) == 'SM') THEN
            is = is + 1
            configs(is)%source_name = DEF_DA_OBS_SOURCE(i)
            configs(is)%sensor_name = DEF_DA_OBS_SENSOR(i)
            configs(is)%var_name    = DEF_DA_OBS_VAR(i)
            configs(is)%fghz        = DEF_DA_OBS_FGHZ(i)
            configs(is)%istream     = i
         ENDIF
      ENDDO

      CALL sm%init(landpatch, configs)

      deallocate(configs)

   END SUBROUTINE init_DA_SM

!-----------------------------------------------------------------------------

   SUBROUTINE run_DA_SM (idate, deltim)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      integer,  intent(in) :: idate(3)
      real(r8), intent(in) :: deltim

!------------------------ Local Variables ------------------------------------
      integer :: np
      integer :: iens
      integer :: il
      integer :: nobs_p                      ! Number of valid observations selected for the target patch
      real(r8), allocatable :: y_p(:)        ! Selected observation values [nobs_p]
      real(r8), allocatable :: hx_p(:,:)     ! Forecast observation equivalents [nobs_p, ens]
      real(r8), allocatable :: r_p(:)        ! Localized observation error variances [nobs_p]
      real(r8), allocatable :: trans(:,:)    ! LETKF ensemble transform [ens, ens]
      real(r8), allocatable :: state_f(:,:)  ! Forecast soil liquid water [soil, ens], kg/m2
      real(r8), allocatable :: state_a(:,:)  ! Transformed soil liquid water [soil, ens], kg/m2

!-----------------------------------------------------------------------------
      ! Prepare observations and forecast H(x).
      CALL sm%calc_hx(idate, deltim)
      CALL sm%save(.false.)

      IF (.not. DEF_DA) RETURN

      IF (p_is_worker .and. numpatch > 0) THEN
         wliq_soisno_ol(:,:) = wliq_soisno_ens(:,0,:)
         wliq_soisno_f (:,:) = sum(wliq_soisno_ens(:,1:DEF_DA_ENS_NUM,:), dim=2) / DEF_DA_ENS_NUM
      ENDIF

      ! Apply the localized LETKF update to each eligible patch.
      IF (p_is_worker) THEN
         DO np = 1, numpatch
            IF (patchtype(np) >= 3) CYCLE

            ! Select observations and localize their error variances.
            CALL sm%prepare(np, nobs_p, y_p, hx_p, r_p)
            IF (nobs_p == 0) CYCLE

            ! Compute the ensemble transform.
            allocate(trans(DEF_DA_ENS_NUM, DEF_DA_ENS_NUM))
            CALL letkf(DEF_DA_ENS_NUM, nobs_p, hx_p, y_p, r_p, sm%infl, trans)

            allocate(state_f(nl_soil, DEF_DA_ENS_NUM))
            allocate(state_a(nl_soil, DEF_DA_ENS_NUM))
            DO iens = 1, DEF_DA_ENS_NUM
               state_f(:,iens) = wliq_soisno_ens(1:,iens,np)
            ENDDO
            CALL dgemm('N', 'N', nl_soil, DEF_DA_ENS_NUM, DEF_DA_ENS_NUM, &
                        1.0_r8, state_f, nl_soil, trans, DEF_DA_ENS_NUM, &
                        0.0_r8, state_a, nl_soil)
            state_a = max(0.0_r8, state_a)
            DO iens = 1, DEF_DA_ENS_NUM
               wliq_soisno_ens(1:2,iens,np) = state_a(1:2,iens)
            ENDDO

            ! postprocess soil moisture to ensure physical consistency
            DO il = 1, nl_soil
               DO iens = 1, DEF_DA_ENS_NUM
                  wliq_soisno_ens(il,iens,np) = max(0.0_r8, wliq_soisno_ens(il,iens,np))
                  wice_soisno_ens(il,iens,np) = max(0.0_r8, wice_soisno_ens(il,iens,np))

                  IF (wliq_soisno_ens(il,iens,np) == 0.0_r8 .and. &
                      wice_soisno_ens(il,iens,np) == 0.0_r8) THEN
                     IF (t_soisno_ens(il,iens,np) - tfrz < -5.0_r8) THEN
                        wice_soisno_ens(il,iens,np) = 1.0e-10_r8
                     ELSE
                        wliq_soisno_ens(il,iens,np) = 1.0e-10_r8
                     ENDIF
                  ENDIF

                  wliq_soisno_ens(il,iens,np) = min(porsl(il,np)*(dz_soi(il)*denh2o), wliq_soisno_ens(il,iens,np))
                  wice_soisno_ens(il,iens,np) = min( &
                     max(0.0_r8, porsl(il,np) - wliq_soisno_ens(il,iens,np)/(dz_soi(il)*denh2o)) &
                     * (dz_soi(il)*denice), wice_soisno_ens(il,iens,np))
               ENDDO
            ENDDO

            wa_ens(1:DEF_DA_ENS_NUM,np) = wa_ens(1:DEF_DA_ENS_NUM,np) - &
               sum(wliq_soisno_ens(1:,1:DEF_DA_ENS_NUM,np) - state_f, dim=1)

            deallocate(y_p, r_p, hx_p, trans)
            deallocate(state_f, state_a)
         ENDDO
      ENDIF

      ! Recompute model- and observation-space H(x), then save analysis diagnostics.
      IF (p_is_worker .and. numpatch > 0) THEN
         wliq_soisno_a(:,:) = sum(wliq_soisno_ens(:,1:DEF_DA_ENS_NUM,:), dim=2) / DEF_DA_ENS_NUM
      ENDIF
      CALL sm%calc_hx(idate, deltim)
      CALL sm%save(.true.)
      CALL sm%clear()

   END SUBROUTINE run_DA_SM

!-----------------------------------------------------------------------------

   SUBROUTINE end_DA_SM()

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      CALL sm%clear()
      IF (allocated(sm%cfgs)) deallocate(sm%cfgs)
      IF (allocated(sm%pm%sensors)) deallocate(sm%pm%sensors)
      IF (allocated(sm%pm%cfgs)) deallocate(sm%pm%cfgs)

      sm%pm%nsensor = 0
      sm%use_pm = .false.
      sm%use_anatq = .false.

   END SUBROUTINE end_DA_SM

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Assim_SM
#endif
