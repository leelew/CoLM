#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Assim_TQ
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Screen-level temperature and humidity assimilation driver using SYNOP
!    observations prepared by Obs_TQ.
!
! WORKFLOW:
!    1. Initialize temperature and humidity observation streams from the namelist.
!    2. Prepare temperature and RH observations from the forecast ensembles.
!    3. Analyze temperature and RH independently, then derive specific humidity.
!    4. Save analysis ensembles and H(x) diagnostics, then clear time-step data.
!
! NOTES:
!    Both analyses use the same forecast background. Temperature increments
!    are limited to +/-20 K and RH is constrained to [0, 1]. Specific humidity
!    is derived from the analyzed temperature and RH. Forecast and analysis
!    ensembles are retained for ANATQ pseudo-observations. Model feedback is
!    handled by the DA ensemble postprocessing step.
!
! HISTORY:
!    Lu Li, 05/2026: Initial implementation
!    Lu Li, Zhilong Fan, 11/2026: DA v1.0
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_LandPatch
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_Vars_Global
   USE MOD_Vars_1DForcing
   USE MOD_DA_Obs_TQ
   USE MOD_DA_EnKF
   USE MOD_DA_Vars_TimeVars
   IMPLICIT NONE

   type(TQ_type) :: tq_tref
   type(TQ_type) :: tq_qref

   ! Products of the current SYNOP assimilation cycle for ANATQ.
   ! Analysis arrays are allocated only after both analysis steps finish.
   real(r8), allocatable :: tref_ens_f  (:,:)
   real(r8), allocatable :: tref_ens_a  (:,:)
   real(r8), allocatable :: rhref_ens_f (:,:)
   real(r8), allocatable :: rhref_ens_a (:,:)

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE init_DA_TQ ()

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Local Variables ------------------------------------
      integer :: is
      integer :: i
      integer :: nsensor_t
      integer :: nsensor_q
      type(TQ_config_type), allocatable :: configs(:)

!-----------------------------------------------------------------------------
      nsensor_t = 0
      nsensor_q = 0
      DO i = 1, size(DEF_DA_OBS_TARGET)
         IF (trim(DEF_DA_OBS_TARGET(i)) == '2MT') nsensor_t = nsensor_t + 1
         IF (trim(DEF_DA_OBS_TARGET(i)) == '2MQ') nsensor_q = nsensor_q + 1
      ENDDO
      IF (nsensor_t == 0 .and. nsensor_q == 0) RETURN

      ! Group streams by analysis target and retain their namelist indices.
      IF (nsensor_t > 0) THEN
         allocate(configs(nsensor_t))
         is = 0
         DO i = 1, size(DEF_DA_OBS_TARGET)
            IF (trim(DEF_DA_OBS_TARGET(i)) == '2MT') THEN
               is = is + 1
               configs(is)%source_name = DEF_DA_OBS_SOURCE(i)
               configs(is)%sensor_name = DEF_DA_OBS_SENSOR(i)
               configs(is)%var_name    = DEF_DA_OBS_VAR(i)
               configs(is)%istream     = i
            ENDIF
         ENDDO
         CALL tq_tref%init(landpatch, configs)
         deallocate(configs)
      ENDIF

      IF (nsensor_q > 0) THEN
         allocate(configs(nsensor_q))
         is = 0
         DO i = 1, size(DEF_DA_OBS_TARGET)
            IF (trim(DEF_DA_OBS_TARGET(i)) == '2MQ') THEN
               is = is + 1
               configs(is)%source_name = DEF_DA_OBS_SOURCE(i)
               configs(is)%sensor_name = DEF_DA_OBS_SENSOR(i)
               configs(is)%var_name    = DEF_DA_OBS_VAR(i)
               configs(is)%istream     = i
            ENDIF
         ENDDO
         CALL tq_qref%init(landpatch, configs)
         deallocate(configs)
      ENDIF

   END SUBROUTINE init_DA_TQ

!-----------------------------------------------------------------------------

   SUBROUTINE run_DA_TQ (idate, deltim)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      integer,  intent(in) :: idate(3)
      real(r8), intent(in) :: deltim

!------------------------ Local Variables ------------------------------------
      integer :: np
      integer :: iens
      integer :: nobs_p                       ! Number of valid observations selected for the target patch
      integer :: ivar                         ! Temperature or RH analysis step
      real(r8) :: infl                        ! Covariance inflation factor for the current target
      real(r8) :: tc                          ! 2 m air temperature [degrees Celsius]
      real(r8) :: es                          ! Saturation vapor pressure [hPa]
      real(r8) :: e                           ! Water vapor pressure [hPa]
      real(r8), allocatable :: y_p(:)         ! Selected observation values [nobs_p]
      real(r8), allocatable :: hx_p(:,:)      ! Forecast observation equivalents [nobs_p, ens]
      real(r8), allocatable :: r_p(:)         ! Localized observation error variances [nobs_p]
      real(r8), allocatable :: trans(:,:)     ! LETKF ensemble transform [ens, ens]
      real(r8), allocatable :: state_f(:,:)   ! Forecast temperature or RH [1, ens], K or 0-1
      real(r8), allocatable :: state_a(:,:)   ! Transformed temperature or RH [1, ens], K or 0-1

!-----------------------------------------------------------------------------
      ! Discard products from the previous assimilation cycle.
      IF (allocated(tref_ens_f )) deallocate (tref_ens_f )
      IF (allocated(tref_ens_a )) deallocate (tref_ens_a )
      IF (allocated(rhref_ens_f)) deallocate (rhref_ens_f)
      IF (allocated(rhref_ens_a)) deallocate (rhref_ens_a)
      IF (.not. allocated(tq_tref%cfgs) .and. .not. allocated(tq_qref%cfgs)) RETURN
      IF (DEF_DA_ENS_NUM <= 1) RETURN

      ! Prepare both observation operators before updating either state.
      IF (allocated(tq_tref%cfgs)) THEN
         CALL tq_tref%calc_obs(idate, deltim)
         CALL tq_tref%save(.false.)
      ENDIF
      IF (allocated(tq_qref%cfgs)) THEN
         CALL tq_qref%calc_obs(idate, deltim)
         CALL tq_qref%save(.false.)
      ENDIF

      IF (p_is_worker .and. numpatch > 0) THEN
         ! Compute relative humidity before saving the prior ensembles.
         DO np = 1, numpatch
            DO iens = 0, DEF_DA_ENS_NUM
               rhref_ens(iens,np) = spval

               IF (tref_ens(iens,np) == spval .or. &
                   qref_ens(iens,np) == spval .or. &
                   forc_psrf(np)    == spval .or. &
                   forc_psrf(np)    <= 0.0_r8) CYCLE

               tc = tref_ens(iens,np) - 273.15_r8
               es = 6.112_r8 * exp((17.67_r8*tc)/(tc + 243.5_r8))
               e  = (qref_ens(iens,np)*forc_psrf(np)/100.0_r8) / &
                    (0.378_r8*qref_ens(iens,np) + 0.622_r8)
               IF (es > 0.0_r8) THEN
                  rhref_ens(iens,np) = max(min(1.0_r8, e/es), 0.0_r8)
               ENDIF
            ENDDO
         ENDDO

         ! Preserve forecasts before either screen-level analysis.
         allocate (tref_ens_f (0:DEF_DA_ENS_NUM, numpatch))
         allocate (rhref_ens_f(0:DEF_DA_ENS_NUM, numpatch))
         tref_ens_f(:,:)  = tref_ens
         rhref_ens_f(:,:) = rhref_ens
      ENDIF

      ! Analyze temperature and RH separately using their forecast backgrounds.
      DO ivar = 1, 2
         SELECT CASE (ivar)
         CASE (1)
            IF (.not. allocated(tq_tref%cfgs)) CYCLE
            infl = tq_tref%infl
         CASE (2)
            IF (.not. allocated(tq_qref%cfgs)) CYCLE
            infl = tq_qref%infl
         END SELECT
         IF (.not. p_is_worker) CYCLE

         DO np = 1, numpatch
            IF (patchtype(np) >= 3) CYCLE

            SELECT CASE (ivar)
            CASE (1)
               IF (any(tref_ens_f(1:DEF_DA_ENS_NUM,np) == spval)) CYCLE
               CALL tq_tref%prepare(np, nobs_p, y_p, hx_p, r_p)
            CASE (2)
               IF (any(rhref_ens_f(1:DEF_DA_ENS_NUM,np) == spval)) CYCLE
               CALL tq_qref%prepare(np, nobs_p, y_p, hx_p, r_p)
            END SELECT
            IF (nobs_p == 0) CYCLE

            allocate(trans(DEF_DA_ENS_NUM, DEF_DA_ENS_NUM))
            CALL letkf(DEF_DA_ENS_NUM, nobs_p, hx_p, y_p, r_p, infl, trans)

            allocate(state_f(1, DEF_DA_ENS_NUM))
            allocate(state_a(1, DEF_DA_ENS_NUM))
            SELECT CASE (ivar)
            CASE (1)
               state_f(1,:) = tref_ens_f(1:DEF_DA_ENS_NUM,np)
            CASE (2)
               state_f(1,:) = rhref_ens_f(1:DEF_DA_ENS_NUM,np)
            END SELECT

            CALL dgemm('N', 'N', 1, DEF_DA_ENS_NUM, DEF_DA_ENS_NUM, &
                        1.0_r8, state_f, 1, trans, DEF_DA_ENS_NUM, &
                        0.0_r8, state_a, 1)

            SELECT CASE (ivar)
            CASE (1)
               ! Limit temperature increments to +/-20 K.
               tref_ens(1:DEF_DA_ENS_NUM,np) = state_f(1,:) + &
                  max(-20.0_r8, min(20.0_r8, state_a(1,:) - state_f(1,:)))
            CASE (2)
               ! Bound the analyzed relative humidity to [0, 1].
               rhref_ens(1:DEF_DA_ENS_NUM,np) = max(0.0_r8, min(1.0_r8, state_a(1,:)))
            END SELECT

            deallocate(y_p, r_p, hx_p, trans, state_f, state_a)
         ENDDO
      ENDDO

      ! Publish the completed analyses and save model-space H(x) diagnostics.
      ! Patches without temperature or RH updates retain their forecast values.
      IF (p_is_worker .and. numpatch > 0) THEN
         ! Convert analyzed RH to specific humidity using analyzed temperature.
         DO np = 1, numpatch
            IF (patchtype(np) >= 3) CYCLE
            IF (forc_psrf(np) == spval .or. forc_psrf(np) <= 0.0_r8) CYCLE
            DO iens = 1, DEF_DA_ENS_NUM
               IF (tref_ens(iens,np) == spval .or. rhref_ens(iens,np) == spval) CYCLE
               IF (tref_ens(iens,np) == tref_ens_f(iens,np) .and. &
                   rhref_ens(iens,np) == rhref_ens_f(iens,np)) CYCLE

               tc = tref_ens(iens,np) - 273.15_r8
               es = 6.112_r8 * exp((17.67_r8*tc)/(tc + 243.5_r8))
               e  = rhref_ens(iens,np)*es
               qref_ens(iens,np) = 0.622_r8*e / (forc_psrf(np)/100.0_r8 - 0.378_r8*e)
            ENDDO
         ENDDO

         ! Preserve completed temperature and RH analyses for ANATQ.
         allocate (tref_ens_a (0:DEF_DA_ENS_NUM, numpatch))
         allocate (rhref_ens_a(0:DEF_DA_ENS_NUM, numpatch))
         tref_ens_a (:,:) = tref_ens
         rhref_ens_a(:,:) = rhref_ens
      ENDIF
      IF (allocated(tq_tref%cfgs)) THEN
         CALL tq_tref%calc_model()
         CALL tq_tref%save(.true.)
      ENDIF
      IF (allocated(tq_qref%cfgs)) THEN
         CALL tq_qref%calc_model()
         CALL tq_qref%save(.true.)
      ENDIF
      IF (allocated(tq_tref%cfgs)) CALL tq_tref%clear()
      IF (allocated(tq_qref%cfgs)) CALL tq_qref%clear()

   END SUBROUTINE run_DA_TQ

!-----------------------------------------------------------------------------

   SUBROUTINE end_DA_TQ ()

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!-----------------------------------------------------------------------------
      IF (allocated(tq_tref%cfgs)) THEN
         CALL tq_tref%clear()
         deallocate(tq_tref%cfgs)
         IF (allocated(tq_tref%synop%cfgs)) deallocate(tq_tref%synop%cfgs)
         IF (allocated(tq_tref%synop%sensors)) deallocate(tq_tref%synop%sensors)
      ENDIF
      IF (allocated(tq_qref%cfgs)) THEN
         CALL tq_qref%clear()
         deallocate(tq_qref%cfgs)
         IF (allocated(tq_qref%synop%cfgs)) deallocate(tq_qref%synop%cfgs)
         IF (allocated(tq_qref%synop%sensors)) deallocate(tq_qref%synop%sensors)
      ENDIF

      IF (allocated(tref_ens_f )) deallocate (tref_ens_f )
      IF (allocated(tref_ens_a )) deallocate (tref_ens_a )
      IF (allocated(rhref_ens_f)) deallocate (rhref_ens_f)
      IF (allocated(rhref_ens_a)) deallocate (rhref_ens_a)

   END SUBROUTINE end_DA_TQ

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Assim_TQ
#endif
