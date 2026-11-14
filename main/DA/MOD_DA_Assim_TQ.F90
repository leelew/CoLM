#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Assim_TQ
!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Screen-level temperature and humidity data assimilation using observations
!
! WORKFLOW:
!   1. Init: Select temperature and humidity observation streams from the
!      namelist and initialize a Data_TQ collection for each analysis target.
!   2. Run: Prepare observations and forecast H(x), then save diagnostics. For
!      assimilation, preserve temperature and relative humidity forecasts.
!   3. End: Release observation data, stream configurations and ANATQ
!      products.
!
! HISTORY:
!   Lu Li, 05/2026: Initial implementation
!   Lu Li, Zhilong Fan, 11/2026: DA v1.0
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_LandPatch
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_Vars_Global
   USE MOD_Vars_1DForcing
   USE MOD_DA_Data_TQ
   USE MOD_DA_EnKF
   USE MOD_DA_Vars_TimeVars
   IMPLICIT NONE

   type(TQ_type)                        :: tq_tref
   type(TQ_type)                        :: tq_qref

   real(r8), allocatable                :: tref_ens_f(:,:)
   real(r8), allocatable                :: tref_ens_a(:,:)
   real(r8), allocatable                :: rhref_ens_f(:,:)
   real(r8), allocatable                :: rhref_ens_a(:,:)

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE init_DA_TQ()

!-----------------------------------------------------------------------------
      integer                           :: is, i, nsensor_t, nsensor_q
      type(TQ_config_type), allocatable :: configs(:)
!-----------------------------------------------------------------------------
      nsensor_t = count(DEF_DA_OBS_TARGET == '2MT')
      nsensor_q = count(DEF_DA_OBS_TARGET == '2MQ')
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

   SUBROUTINE run_DA_TQ(idate, deltim)

!-----------------------------------------------------------------------------
      integer,  intent(in)              :: idate(3)
      real(r8), intent(in)              :: deltim

      integer                           :: np, iens                               ! Loop indices
      integer                           :: nobs_p                                 ! Number of observations selected for the target patch
      logical                           :: has_analysis                           ! Whether temperature or RH analysis was applied to the patch
      real(r8)                          :: tc                                     ! 2 m air temperature [degrees Celsius]
      real(r8)                          :: es                                     ! Saturation vapor pressure [hPa]
      real(r8)                          :: e                                      ! Water vapor pressure [hPa]
      real(r8), allocatable             :: y_p(:)                                 ! Selected observation values
      real(r8), allocatable             :: hx_p(:,:)                              ! Forecast ensemble H(x) at selected observation locations
      real(r8), allocatable             :: r_p(:)                                 ! Localized observation error variances
      real(r8)                          :: trans(DEF_DA_ENS_NUM, DEF_DA_ENS_NUM)  ! LETKF ensemble transform matrix
      real(r8)                          :: state_f(1, DEF_DA_ENS_NUM)             ! Forecast temperature [K] or relative humidity [0-1]
      real(r8)                          :: state_a(1, DEF_DA_ENS_NUM)             ! Transformed temperature [K] or relative humidity [0-1]
!-----------------------------------------------------------------------------
      ! Discard products from the previous assimilation cycle.
      IF (allocated(tref_ens_f )) deallocate(tref_ens_f )
      IF (allocated(tref_ens_a )) deallocate(tref_ens_a )
      IF (allocated(rhref_ens_f)) deallocate(rhref_ens_f)
      IF (allocated(rhref_ens_a)) deallocate(rhref_ens_a)

      ! Prepare both observation operators before updating either state.
      IF (tq_tref%use_insitu) THEN
         CALL tq_tref%forward(idate, deltim)
         CALL tq_tref%save(.false.)
      ENDIF
      IF (tq_qref%use_insitu) THEN
         CALL tq_qref%forward(idate, deltim)
         CALL tq_qref%save(.false.)
      ENDIF

      IF (.not. DEF_DA) THEN
         IF (tq_tref%use_insitu) CALL tq_tref%clear()
         IF (tq_qref%use_insitu) CALL tq_qref%clear()
         RETURN
      ENDIF

      IF (p_is_worker) THEN
         allocate(tref_ens_f (0:DEF_DA_ENS_NUM, numpatch))
         allocate(rhref_ens_f(0:DEF_DA_ENS_NUM, numpatch))
         allocate(tref_ens_a (0:DEF_DA_ENS_NUM, numpatch))
         allocate(rhref_ens_a(0:DEF_DA_ENS_NUM, numpatch))

         DO np = 1, numpatch
            has_analysis = .false.
            ! Compute relative humidity and preserve the forecast background.
            DO iens = 0, DEF_DA_ENS_NUM
               IF (tref_ens(iens,np) == spval .or. &
                   qref_ens(iens,np) == spval .or. &
                   forc_psrf(np)     == spval) CYCLE

               tc = tref_ens(iens,np) - 273.15_r8
               es = 6.112_r8 * exp((17.67_r8*tc)/(tc + 243.5_r8))
               e  = (qref_ens(iens,np)*forc_psrf(np)/100.0_r8) / &
                  (0.378_r8*qref_ens(iens,np) + 0.622_r8)
               IF (es > 0.0_r8) rhref_ens(iens,np) = max(min(1.0_r8, e/es), 0.0_r8)
            ENDDO
            tref_ens_f (:,np) = tref_ens(:,np)
            rhref_ens_f(:,np) = rhref_ens(:,np)

            ! Analyze temperature and RH independently on soil patches.
            IF (patchtype(np) == 0) THEN
               ! Analyze 2 m temperature using its forecast background.
               IF (tq_tref%use_insitu) THEN
                  CALL tq_tref%prepare(np, nobs_p, y_p, hx_p, r_p)
                  IF (nobs_p > 0) THEN
                     state_f(1,:) = tref_ens_f(1:DEF_DA_ENS_NUM,np)
                     CALL letkf(DEF_DA_ENS_NUM, nobs_p, hx_p, y_p, r_p, tq_tref%infl, trans)
                     CALL dgemm('N', 'N', 1, DEF_DA_ENS_NUM, DEF_DA_ENS_NUM, &
                          1.0_r8, state_f, 1, trans, DEF_DA_ENS_NUM,         &
                          0.0_r8, state_a, 1)

                     ! Limit temperature increments to +/-20 K.
                     tref_ens(1:DEF_DA_ENS_NUM,np) = state_f(1,:) + &
                        max(-20.0_r8, min(20.0_r8, state_a(1,:) - state_f(1,:)))
                     has_analysis = .true.
                     deallocate(y_p, r_p, hx_p)
                  ENDIF
               ENDIF

               ! Analyze 2 m relative humidity using its forecast background.
               IF (tq_qref%use_insitu) THEN
                  CALL tq_qref%prepare(np, nobs_p, y_p, hx_p, r_p)
                  IF (nobs_p > 0) THEN
                     state_f(1,:) = rhref_ens_f(1:DEF_DA_ENS_NUM,np)
                     CALL letkf(DEF_DA_ENS_NUM, nobs_p, hx_p, y_p, r_p, tq_qref%infl, trans)
                     CALL dgemm('N', 'N', 1, DEF_DA_ENS_NUM, DEF_DA_ENS_NUM, &
                          1.0_r8, state_f, 1, trans, DEF_DA_ENS_NUM,         &
                          0.0_r8, state_a, 1)

                     ! Bound analyzed relative humidity to [0, 1].
                     rhref_ens(1:DEF_DA_ENS_NUM,np) = max(0.0_r8, min(1.0_r8, state_a(1,:)))
                     has_analysis = .true.
                     deallocate(y_p, r_p, hx_p)
                  ENDIF
               ENDIF

               ! Recompute specific humidity only after a local analysis.
               IF (has_analysis) THEN
                  DO iens = 1, DEF_DA_ENS_NUM
                     tc = tref_ens(iens,np) - 273.15_r8
                     es = 6.112_r8 * exp((17.67_r8*tc)/(tc + 243.5_r8))
                     e  = rhref_ens(iens,np)*es
                     qref_ens(iens,np) = 0.622_r8*e / (forc_psrf(np)/100.0_r8 - 0.378_r8*e)
                  ENDDO
               ENDIF
            ENDIF

            ! Preserve complete patch analyses for ANATQ.
            tref_ens_a (:,np) = tref_ens(:,np)
            rhref_ens_a(:,np) = rhref_ens(:,np)
         ENDDO
      ENDIF

      ! Recompute model- and observation-space H(x), then save analysis diagnostics.
      IF (tq_tref%use_insitu) THEN
         CALL tq_tref%forward(idate, deltim)
         CALL tq_tref%save(.true.)
      ENDIF
      IF (tq_qref%use_insitu) THEN
         CALL tq_qref%forward(idate, deltim)
         CALL tq_qref%save(.true.)
      ENDIF
      IF (tq_tref%use_insitu) CALL tq_tref%clear()
      IF (tq_qref%use_insitu) CALL tq_qref%clear()

   END SUBROUTINE run_DA_TQ

!-----------------------------------------------------------------------------

   SUBROUTINE end_DA_TQ()

!-----------------------------------------------------------------------------
      IF (tq_tref%use_insitu) THEN
         CALL tq_tref%clear()
         deallocate(tq_tref%cfgs)
         IF (allocated(tq_tref%insitu%cfgs   )) deallocate(tq_tref%insitu%cfgs   )
         IF (allocated(tq_tref%insitu%sensors)) deallocate(tq_tref%insitu%sensors)
         tq_tref%insitu%nsensor = 0
         tq_tref%use_insitu = .false.
      ENDIF
      IF (tq_qref%use_insitu) THEN
         CALL tq_qref%clear()
         deallocate(tq_qref%cfgs)
         IF (allocated(tq_qref%insitu%cfgs   )) deallocate(tq_qref%insitu%cfgs   )
         IF (allocated(tq_qref%insitu%sensors)) deallocate(tq_qref%insitu%sensors)
         tq_qref%insitu%nsensor = 0
         tq_qref%use_insitu = .false.
      ENDIF

      IF (allocated(tref_ens_f )) deallocate(tref_ens_f )
      IF (allocated(tref_ens_a )) deallocate(tref_ens_a )
      IF (allocated(rhref_ens_f)) deallocate(rhref_ens_f)
      IF (allocated(rhref_ens_a)) deallocate(rhref_ens_a)

   END SUBROUTINE end_DA_TQ

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Assim_TQ
#endif
