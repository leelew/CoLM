#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Assim_SM
!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Soil moisture data assimilation using observations.
!
! WORKFLOW:
!   1. Init: Select soil moisture observation streams from the namelist and
!      initialize Data_SM with their configurations.
!   2. Run: Prepare observations and forecast H(x), then save diagnostics. For
!      assimilation, select nearby observations for each soil patch and apply
!      LETKF to update liquid water in the upper two soil layers. Enforce
!      soil-water constraints and adjust aquifer storage for the change in
!      soil liquid water. Recompute H(x) and save analysis diagnostics.
!   3. End: Release observation data and stream configurations.
!
! HISTORY:
!   Lu Li, Zhilong Fan, 11/2026: DA v1.0
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Vars_TimeInvariants
   USE MOD_Vars_Global
   USE MOD_LandPatch
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_Const_Physical
   USE MOD_DA_Data_SM
   USE MOD_DA_EnKF
   USE MOD_DA_Vars_TimeVars
   IMPLICIT NONE

   type(SM_type)                        :: sm

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE init_DA_SM()

!-----------------------------------------------------------------------------
      integer                           :: is, i, nsensor
      type(SM_config_type), allocatable :: configs(:)
!-----------------------------------------------------------------------------
      nsensor = count(DEF_DA_OBS_TARGET == 'SM')
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

   SUBROUTINE run_DA_SM(idate, deltim)

!-----------------------------------------------------------------------------
      integer,  intent(in)              :: idate(3)
      real(r8), intent(in)              :: deltim

      integer                           :: np, iens, il
      integer                           :: nobs_p                                 ! Number of observations selected for the target patch
      real(r8), allocatable             :: y_p(:)                                 ! Selected observation values
      real(r8), allocatable             :: hx_p(:,:)                              ! Forecast ensemble H(x) at selected observation locations
      real(r8), allocatable             :: r_p(:)                                 ! Localized observation error variances
      real(r8)                          :: trans(DEF_DA_ENS_NUM, DEF_DA_ENS_NUM)  ! LETKF ensemble transform matrix
      real(r8)                          :: state_f(nl_soil, DEF_DA_ENS_NUM)       ! Forecast soil liquid water [kg/m2]
      real(r8)                          :: state_a(nl_soil, DEF_DA_ENS_NUM)       ! Transformed soil liquid water [kg/m2]
!-----------------------------------------------------------------------------
      ! Prepare observations and forecast H(x).
      CALL sm%forward(idate, deltim)
      CALL sm%save(.false.)

      IF (.not. DEF_DA) RETURN

      ! Apply the localized LETKF update to each soil patch.
      IF (p_is_worker) THEN
         DO np = 1, numpatch
            ! Save open-loop and forecast states before analysis.
            wliq_soisno_ol(:,np) = wliq_soisno_ens(:,0,np)
            wliq_soisno_f (:,np) = sum(wliq_soisno_ens(:,1:DEF_DA_ENS_NUM,np), dim=2) / DEF_DA_ENS_NUM

            IF (patchtype(np) >= 1) CYCLE

            ! Select observations and localize their error variances.
            CALL sm%prepare(np, nobs_p, y_p, hx_p, r_p)
            IF (nobs_p == 0) CYCLE

            ! Compute the ensemble transform.
            CALL letkf(DEF_DA_ENS_NUM, nobs_p, hx_p, y_p, r_p, sm%infl, trans)

            ! Transform the soil profile and apply increments to the upper two layers.
            DO iens = 1, DEF_DA_ENS_NUM
               state_f(:,iens) = wliq_soisno_ens(1:,iens,np)
            ENDDO
            CALL dgemm('N', 'N', nl_soil, DEF_DA_ENS_NUM, DEF_DA_ENS_NUM, &
               1.0_r8, state_f, nl_soil, trans, DEF_DA_ENS_NUM,           &
               0.0_r8, state_a, nl_soil)
            state_a = max(0.0_r8, state_a)
            DO iens = 1, DEF_DA_ENS_NUM
               wliq_soisno_ens(1:2,iens,np) = state_a(1:2,iens)
            ENDDO

            ! Enforce nonnegative water content and soil pore-space limits.
            DO il = 1, nl_soil
               DO iens = 1, DEF_DA_ENS_NUM
                  wliq_soisno_ens(il,iens,np) = max(0.0_r8, wliq_soisno_ens(il,iens,np))
                  wice_soisno_ens(il,iens,np) = max(0.0_r8, wice_soisno_ens(il,iens,np))

                  IF (wliq_soisno_ens(il,iens,np) == 0.0_r8 .and. wice_soisno_ens(il,iens,np) == 0.0_r8) THEN
                     IF (t_soisno_ens(il,iens,np) - tfrz < -5.0_r8) THEN
                        wice_soisno_ens(il,iens,np) = 1.0e-10_r8
                     ELSE
                        wliq_soisno_ens(il,iens,np) = 1.0e-10_r8
                     ENDIF
                  ENDIF

                  wliq_soisno_ens(il,iens,np) = min(porsl(il,np)*(dz_soi(il)*denh2o), wliq_soisno_ens(il,iens,np))
                  wice_soisno_ens(il,iens,np) = min(                                             &
                     max(0.0_r8, porsl(il,np) - wliq_soisno_ens(il,iens,np)/(dz_soi(il)*denh2o)) &
                     * (dz_soi(il)*denice), wice_soisno_ens(il,iens,np))
               ENDDO
            ENDDO

            ! Offset the soil liquid-water increment with aquifer storage.
            wa_ens(1:DEF_DA_ENS_NUM,np) = wa_ens(1:DEF_DA_ENS_NUM,np) - &
               sum(wliq_soisno_ens(1:,1:DEF_DA_ENS_NUM,np) - state_f, dim=1)

            deallocate(y_p, r_p, hx_p)
         ENDDO
      ENDIF

      ! Recompute model- and observation-space H(x), then save analysis diagnostics.
      IF (p_is_worker) THEN
         wliq_soisno_a(:,:) = sum(wliq_soisno_ens(:,1:DEF_DA_ENS_NUM,:), dim=2) / DEF_DA_ENS_NUM
      ENDIF
      CALL sm%forward(idate, deltim)
      CALL sm%save(.true.)
      CALL sm%clear()

   END SUBROUTINE run_DA_SM

!-----------------------------------------------------------------------------

   SUBROUTINE end_DA_SM()

!-----------------------------------------------------------------------------
      CALL sm%clear()
      IF (allocated(sm%cfgs      )) deallocate(sm%cfgs      )
      IF (allocated(sm%pm%sensors)) deallocate(sm%pm%sensors)
      IF (allocated(sm%pm%cfgs   )) deallocate(sm%pm%cfgs   )

      sm%pm%nsensor = 0
      sm%use_pm = .false.
      sm%use_anatq = .false.

   END SUBROUTINE end_DA_SM

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Assim_SM
#endif
