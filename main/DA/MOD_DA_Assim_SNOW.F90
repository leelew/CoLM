#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Assim_SNOW
!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Snow-depth data assimilation using MODIS fractional snow cover observations.
!
! WORKFLOW:
!   1. Init: Select snow observation streams from the namelist and initialize
!      Data_SNOW with their configurations.
!   2. Run: Prepare observations and forecast H(x), then save diagnostics. For
!      assimilation, select nearby observations for each supported patch and
!      apply LETKF to update snow depth. Reconstruct snow layers and adjust
!      aquifer storage for the change in snow water equivalent. Recompute H(x)
!      and save analysis diagnostics.
!   3. End: Release observation data and stream configurations.
!
! HISTORY:
!   Lu Li, Zhuo Wang, 11/2026: DA v1.0
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Vars_TimeInvariants
   USE MOD_Vars_Global
   USE MOD_LandPatch
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_Const_Physical, ONLY: tfrz
   USE MOD_Vars_1DForcing, ONLY: forc_t, forc_us, forc_vs
   USE MOD_DA_Data_SNOW
   USE MOD_DA_EnKF
   USE MOD_DA_Vars_TimeVars
   USE MOD_SnowLayersCombineDivide
   IMPLICIT NONE

   type(SNOW_type)                        :: snow
   integer, parameter                     :: lb = maxsnl + 1

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE init_DA_SNOW()

!-----------------------------------------------------------------------------
      integer                             :: is, i, nsensor
      type(Snow_config_type), allocatable :: configs(:)
!-----------------------------------------------------------------------------
      nsensor = count(DEF_DA_OBS_TARGET == 'SNOW')
      IF (nsensor == 0) RETURN

      allocate(configs(nsensor))
      is = 0
      DO i = 1, size(DEF_DA_OBS_TARGET)
         IF (trim(DEF_DA_OBS_TARGET(i)) == 'SNOW') THEN
            is = is + 1
            configs(is)%source_name = DEF_DA_OBS_SOURCE(i)
            configs(is)%sensor_name = DEF_DA_OBS_SENSOR(i)
            configs(is)%var_name    = DEF_DA_OBS_VAR(i)
            configs(is)%istream     = i
         ENDIF
      ENDDO

      CALL snow%init(landpatch, configs)

      deallocate(configs)

   END SUBROUTINE init_DA_SNOW

!-----------------------------------------------------------------------------

   SUBROUTINE run_DA_SNOW(idate, deltim)

!-----------------------------------------------------------------------------
      integer,  intent(in)                :: idate(3)
      real(r8), intent(in)                :: deltim

      integer                             :: np, iens, kk, snl_local
      integer                             :: nobs_p                                 ! Number of observations selected for the target patch
      real(r8), allocatable               :: y_p(:)                                 ! Selected observation values
      real(r8), allocatable               :: hx_p(:,:)                              ! Forecast ensemble H(x) at selected observation locations
      real(r8), allocatable               :: r_p(:)                                 ! Localized observation error variances
      real(r8)                            :: trans(DEF_DA_ENS_NUM, DEF_DA_ENS_NUM)  ! LETKF ensemble transform matrix
      real(r8)                            :: state_f(1, DEF_DA_ENS_NUM)             ! Forecast snow depth [m]
      real(r8)                            :: state_a(1, DEF_DA_ENS_NUM)             ! Transformed snow depth [m]
      real(r8)                            :: old_scv_ens(DEF_DA_ENS_NUM)            ! Forecast snow water equivalent [mm]
      real(r8)                            :: current_rho, bifall                    ! Existing and new snow density [kg/m3]
      real(r8)                            :: patch_wind, patch_tair
      real(r8)                            :: delta_swe, swe_after
      real(r8)                            :: depth_before, mass_before, depth_scale, mass_scale
      real(r8), parameter                 :: thin_snow = 0.01_r8       ! Snow-node initialization threshold [m]
      real(r8), parameter                 :: snow_eps = 1.0e-5_r8
      real(r8)                            :: zi_snow_local(lb-1:0)                  ! Snow-layer interface depths [m]
!-----------------------------------------------------------------------------
      ! Prepare observations and forecast H(x).
      CALL snow%forward(idate, deltim)
      CALL snow%save(.false.)

      IF (.not. DEF_DA) THEN
         CALL snow%clear()
         RETURN
      ENDIF

      ! Apply the localized LETKF update to each supported patch.
      IF (p_is_worker .and. numpatch > 0) THEN
         snowdp_ol(:) = spval
         snowdp_f (:) = spval
         DO np = 1, numpatch
            IF (patchtype(np) >= 2) CYCLE

            ! Save open-loop and forecast states before analysis.
            snowdp_ol(np) = snowdp_ens(0,np)
            IF (any(snowdp_ens(1:DEF_DA_ENS_NUM,np) == spval)) CYCLE
            snowdp_f(np) = sum(snowdp_ens(1:DEF_DA_ENS_NUM,np)) / DEF_DA_ENS_NUM
            IF (DEF_DA_ENS_NUM <= 1) CYCLE
            IF (any(scv_ens(1:DEF_DA_ENS_NUM,np) == spval)) CYCLE

            ! Select observations and localize their error variances.
            CALL snow%prepare(np, nobs_p, y_p, hx_p, r_p)
            IF (nobs_p == 0) CYCLE

            ! Compute the ensemble transform.
            CALL letkf(DEF_DA_ENS_NUM, nobs_p, hx_p, y_p, r_p, snow%infl, trans)

            ! Transform snow depth and enforce nonnegative values.
            state_f(1,:) = snowdp_ens(1:DEF_DA_ENS_NUM,np)
            old_scv_ens = scv_ens(1:DEF_DA_ENS_NUM,np)
            CALL dgemm('N', 'N', 1, DEF_DA_ENS_NUM, DEF_DA_ENS_NUM, &
               1.0_r8, state_f, 1, trans, DEF_DA_ENS_NUM,           &
               0.0_r8, state_a, 1)
            state_a = max(0.0_r8, state_a)

            ! Estimate new snow density from air temperature and wind speed.
            patch_tair = forc_t(np)
            patch_wind = sqrt(forc_us(np)**2 + forc_vs(np)**2)
            bifall = 109.0_r8 + 6.0_r8*(patch_tair - tfrz) + 26.0_r8*sqrt(max(0.0_r8, patch_wind))
            bifall = max(50.0_r8, min(500.0_r8, bifall))

            ! Update snow storage while preserving existing layer states.
            DO iens = 1, DEF_DA_ENS_NUM
               IF (abs(state_a(1,iens) - state_f(1,iens)) <= 1.0e-12_r8) CYCLE
               snowdp_ens(iens,np) = state_a(1,iens)
               IF (snowdp_ens(iens,np) > snow_eps) THEN
                  ! Retain the original bulk density for existing snow.
                  IF (state_f(1,iens) > snow_eps .and. old_scv_ens(iens) > 0.0_r8) THEN
                     current_rho = old_scv_ens(iens) / state_f(1,iens)
                     current_rho = max(50.0_r8, min(900.0_r8, current_rho))
                     scv_ens(iens,np) = snowdp_ens(iens,np) * current_rho
                  ELSE
                     scv_ens(iens,np) = snowdp_ens(iens,np) * bifall
                  ENDIF
               ELSE
                  snowdp_ens(iens,np) = 0.0_r8
                  scv_ens(iens,np) = 0.0_r8
               ENDIF

               ! Diagnose existing snow nodes as in CoLMMAIN.
               snl_local = -count(wice_soisno_ens(lb:0,iens,np) + &
                                  wliq_soisno_ens(lb:0,iens,np) > 0.0_r8)
               IF (snowdp_ens(iens,np) < thin_snow) THEN
                  ! Retain bulk snow storage without explicit snow nodes.
                  snl_local = 0
               ELSE IF (snl_local < 0 .and. state_f(1,iens) > snow_eps .and. &
                        old_scv_ens(iens) > 0.0_r8) THEN
                  ! Scale existing layers without resetting temperature or liquid fraction.
                  depth_before = sum(dz_sno_ens(snl_local+1:0,iens,np))
                  mass_before = sum(wice_soisno_ens(snl_local+1:0,iens,np)) + &
                                sum(wliq_soisno_ens(snl_local+1:0,iens,np))
                  IF (depth_before <= 0.0_r8 .or. mass_before <= 0.0_r8 .or. &
                      any(dz_sno_ens(snl_local+1:0,iens,np) <= 0.0_r8) .or. &
                      any(wice_soisno_ens(snl_local+1:0,iens,np) < 0.0_r8) .or. &
                      any(wliq_soisno_ens(snl_local+1:0,iens,np) < 0.0_r8) .or. &
                      any(wice_soisno_ens(snl_local+1:0,iens,np) + &
                          wliq_soisno_ens(snl_local+1:0,iens,np) <= 0.0_r8)) THEN
                     WRITE(*,*) '[SNOW-INVALID-LAYERS]', np, iens
                     snowdp_ens(iens,np) = state_f(1,iens)
                     scv_ens(iens,np) = old_scv_ens(iens)
                     CYCLE
                  ENDIF
                  depth_scale = snowdp_ens(iens,np) / depth_before
                  mass_scale = scv_ens(iens,np) / mass_before
                  dz_sno_ens(snl_local+1:0,iens,np) = dz_sno_ens(snl_local+1:0,iens,np) * depth_scale
                  wice_soisno_ens(snl_local+1:0,iens,np) = wice_soisno_ens(snl_local+1:0,iens,np) * mass_scale
                  wliq_soisno_ens(snl_local+1:0,iens,np) = wliq_soisno_ens(snl_local+1:0,iens,np) * mass_scale
               ELSE
                  ! Initialize a snow node only when no explicit snow layers exist.
                  snl_local = -1
                  dz_sno_ens(0,iens,np) = snowdp_ens(iens,np)
                  wice_soisno_ens(0,iens,np) = scv_ens(iens,np)
                  wliq_soisno_ens(0,iens,np) = 0.0_r8
                  t_soisno_ens(0,iens,np) = min(tfrz, forc_t(np))
               ENDIF

               ! Clear unused nodes and reconstruct active snow-layer geometry.
               IF (snl_local > maxsnl) THEN
                  dz_sno_ens(lb:snl_local,iens,np) = 0.0_r8
                  wice_soisno_ens(lb:snl_local,iens,np) = 0.0_r8
                  wliq_soisno_ens(lb:snl_local,iens,np) = 0.0_r8
                  t_soisno_ens(lb:snl_local,iens,np) = 0.0_r8
               ENDIF
               zi_snow_local = 0.0_r8
               z_sno_ens(lb:0,iens,np) = 0.0_r8
               IF (snl_local < 0) THEN
                  DO kk = 0, snl_local+1, -1
                     zi_snow_local(kk-1) = zi_snow_local(kk) - dz_sno_ens(kk,iens,np)
                     z_sno_ens(kk,iens,np) = 0.5_r8 * (zi_snow_local(kk-1) + zi_snow_local(kk))
                  ENDDO
                  CALL snowlayersdivide(lb, snl_local, &
                        z_sno_ens(lb:0,iens,np), &
                        dz_sno_ens(lb:0,iens,np), &
                        zi_snow_local, &
                        wliq_soisno_ens(lb:0,iens,np), &
                        wice_soisno_ens(lb:0,iens,np), &
                        t_soisno_ens(lb:0,iens,np))
                  swe_after = sum(wice_soisno_ens(snl_local+1:0,iens,np)) + &
                              sum(wliq_soisno_ens(snl_local+1:0,iens,np))
                  IF (abs(swe_after - scv_ens(iens,np)) > 1.0e-8_r8 * max(1.0_r8, scv_ens(iens,np))) THEN
                     WRITE(*,*) '[SNOW-DIVIDE-ERROR]', np, iens, scv_ens(iens,np), swe_after
                  ENDIF
               ENDIF

               ! Offset the SWE increment with aquifer storage.
               delta_swe = scv_ens(iens,np) - old_scv_ens(iens)
               wa_ens(iens,np) = wa_ens(iens,np) - delta_swe
            ENDDO

            deallocate(y_p, r_p, hx_p)
         ENDDO
      ENDIF

      ! Recompute model- and observation-space H(x), then save analysis diagnostics.
      IF (p_is_worker .and. numpatch > 0) THEN
         snowdp_a(:) = spval
         DO np = 1, numpatch
            IF (patchtype(np) >= 2) CYCLE
            IF (any(snowdp_ens(1:DEF_DA_ENS_NUM,np) == spval)) CYCLE
            snowdp_a(np) = sum(snowdp_ens(1:DEF_DA_ENS_NUM,np)) / DEF_DA_ENS_NUM
         ENDDO
      ENDIF
      CALL snow%forward(idate, deltim)
      CALL snow%save(.true.)
      CALL snow%clear()

   END SUBROUTINE run_DA_SNOW

!-----------------------------------------------------------------------------

   SUBROUTINE end_DA_SNOW()

!-----------------------------------------------------------------------------
      CALL snow%clear()
      IF (allocated(snow%cfgs      )) deallocate(snow%cfgs      )
      IF (allocated(snow%mr%sensors)) deallocate(snow%mr%sensors)
      IF (allocated(snow%mr%cfgs   )) deallocate(snow%mr%cfgs   )

      snow%mr%nsensor = 0
      snow%use_mr = .false.

   END SUBROUTINE end_DA_SNOW

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Assim_SNOW
#endif

