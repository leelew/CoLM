#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Assim_2MTQ
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Main 2 m temperature and humidity assimilation driver for SYNOP
!    observations. The two variables are assimilated independently: a '2MT'
!    stream updates the 2 m air temperature ensemble tref_ens only, and a
!    '2MQ' stream updates the 2 m specific humidity ensemble qref_ens only.
!    The model humidity state is specific humidity, so H(x) maps qref_ens to
!    relative humidity with tref_ens and surface pressure before the LETKF
!    update (see MOD_DA_Obs_SYNOP).
!
!    One LETKF transform is computed per patch and per enabled target and
!    applied directly to the ensemble of the corresponding variable. The 2MT
!    analysis runs before the 2MQ analysis, so the humidity observation
!    operator sees the analyzed temperature. The analyzed humidity is
!    constrained to [0, qsat] and the temperature increment is limited to
!    +/-20 K to keep the analysis physical.
!
!    Because tref_ens and rhref_ens are derived from the analysis variables,
!    the analysis H(x) diagnostics are available immediately after the
!    update.
!
! HISTORY:
!    Lu Li, 05/2026: Initial implementation
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_DA_Obs_2MTQ
   USE MOD_DA_EnKF
   USE MOD_DA_Vars_TimeVars
   USE MOD_LandPatch
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_Vars_Global, ONLY: spval
   USE MOD_Vars_1DForcing, ONLY: forc_psrf
   USE MOD_Qsadv, ONLY: qsadv
   IMPLICIT NONE

   type(TQ_type) :: tq_t                      ! '2MT' streams: 2 m temperature
   type(TQ_type) :: tq_q                      ! '2MQ' target: specific humidity analyzed using 2MRH

   logical :: use_2mt = .false.               ! any '2MT' stream configured
   logical :: use_2mq = .false.               ! any '2MQ' stream configured

   ! Products of the current SYNOP assimilation cycle for ANA2MTQ.
   ! Ensemble arrays use [0:DEF_DA_ENS_NUM, numpatch].
   ! Analysis arrays are allocated only after both analysis steps finish.
   real(r8), allocatable :: tref_ens_f  (:,:) ! Forecast 2 m temperature ensemble [K]
   real(r8), allocatable :: tref_ens_a  (:,:) ! Analysis 2 m temperature ensemble [K]
   real(r8), allocatable :: tref_f      (:)   ! Forecast ensemble-mean 2 m temperature [K]
   real(r8), allocatable :: tref_a      (:)   ! Analysis ensemble-mean 2 m temperature [K]
   real(r8), allocatable :: qref_ens_f  (:,:) ! Forecast 2 m specific humidity ensemble [kg/kg]
   real(r8), allocatable :: qref_ens_a  (:,:) ! Analysis 2 m specific humidity ensemble [kg/kg]
   real(r8), allocatable :: qref_f      (:)   ! Forecast ensemble-mean 2 m specific humidity [kg/kg]
   real(r8), allocatable :: qref_a      (:)   ! Analysis ensemble-mean 2 m specific humidity [kg/kg]
   real(r8), allocatable :: rhref_ens_f (:,:) ! Forecast 2 m relative humidity ensemble [0-1]
   real(r8), allocatable :: rhref_ens_a (:,:) ! Analysis 2 m relative humidity ensemble [0-1]

   PUBLIC :: init_DA_2MTQ
   PUBLIC :: run_DA_2MTQ
   PUBLIC :: end_DA_2MTQ

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE init_DA_2MTQ ()

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      integer :: is, i, nsensor_t, nsensor_q
      type(TQ_config_type), allocatable :: configs(:)

!-----------------------------------------------------------------------------
      nsensor_t = 0
      nsensor_q = 0
      DO i = 1, size(DEF_DA_OBS_TARGET)
         IF (trim(DEF_DA_OBS_TARGET(i)) == '2MT') nsensor_t = nsensor_t + 1
         IF (trim(DEF_DA_OBS_TARGET(i)) == '2MQ') nsensor_q = nsensor_q + 1
      ENDDO
      use_2mt = nsensor_t > 0
      use_2mq = nsensor_q > 0
      IF (.not. use_2mt .and. .not. use_2mq) RETURN

      ! Group SYNOP streams by analysis target; VAR selects the observation operator.
      ! istream is the namelist index of the stream, consistent with the
      ! other assimilation drivers and the shared hx_* diagnostics, so the
      ! DEF_DA_OBS_TARGET entries must be filled contiguously from 1.
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
         CALL tq_t%init(landpatch, configs)
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
         CALL tq_q%init(landpatch, configs)
         deallocate(configs)
      ENDIF

   END SUBROUTINE init_DA_2MTQ

!-----------------------------------------------------------------------------

   SUBROUTINE run_DA_2MTQ (idate, deltim)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      integer,  intent(in) :: idate(3)
      real(r8), intent(in) :: deltim

!------------------------ Local Variables ------------------------------------
      integer :: np
      integer :: iens
      integer :: nobs_p
      real(r8) :: es, esdT, qs, qsdT
      real(r8) :: tc_rh                       ! 2 m air temperature [degrees Celsius]
      real(r8) :: es_rh                       ! Saturation vapor pressure [hPa]
      real(r8) :: e_rh                        ! Water vapor pressure [hPa]
      real(r8), allocatable :: y_p(:)
      real(r8), allocatable :: hx_p(:,:)
      real(r8), allocatable :: r_p(:)
      real(r8), allocatable :: trans(:,:)
      real(r8), allocatable :: tref_work_f(:,:)
      real(r8), allocatable :: tref_work_a(:,:)
      real(r8), allocatable :: qref_work_f(:,:)
      real(r8), allocatable :: qref_work_a(:,:)

!-----------------------------------------------------------------------------
      ! Discard products from the previous assimilation cycle.
      IF (allocated(tref_ens_f )) deallocate (tref_ens_f )
      IF (allocated(tref_ens_a )) deallocate (tref_ens_a )
      IF (allocated(tref_f     )) deallocate (tref_f     )
      IF (allocated(tref_a     )) deallocate (tref_a     )
      IF (allocated(qref_ens_f )) deallocate (qref_ens_f )
      IF (allocated(qref_ens_a )) deallocate (qref_ens_a )
      IF (allocated(qref_f     )) deallocate (qref_f     )
      IF (allocated(qref_a     )) deallocate (qref_a     )
      IF (allocated(rhref_ens_f)) deallocate (rhref_ens_f)
      IF (allocated(rhref_ens_a)) deallocate (rhref_ens_a)
      IF (.not. use_2mt .and. .not. use_2mq) RETURN
      IF (DEF_DA_ENS_NUM <= 1) RETURN

      ! read the SYNOP observations and calculate the prior H(x) at station
      ! sites for each enabled target
      IF (use_2mt) THEN
         CALL tq_t%calc_obs(idate, deltim)
         CALL tq_t%save(.false.)
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

               tc_rh = tref_ens(iens,np) - 273.15_r8
               es_rh = 6.112_r8 * exp((17.67_r8*tc_rh)/(tc_rh + 243.5_r8))
               e_rh  = (qref_ens(iens,np)*forc_psrf(np)/100.0_r8) / &
                    (0.378_r8*qref_ens(iens,np) + 0.622_r8)
               IF (es_rh > 0.0_r8) THEN
                  rhref_ens(iens,np) = max(min(1.0_r8, e_rh/es_rh), 0.0_r8)
               ENDIF
            ENDDO
         ENDDO

         ! Preserve forecasts before either screen-level analysis.
         allocate (tref_ens_f (0:DEF_DA_ENS_NUM, numpatch))
         allocate (tref_f     (numpatch))
         allocate (qref_ens_f (0:DEF_DA_ENS_NUM, numpatch))
         allocate (qref_f     (numpatch))
         allocate (rhref_ens_f(0:DEF_DA_ENS_NUM, numpatch))
         tref_ens_f(:,:) = tref_ens
         qref_ens_f(:,:) = qref_ens
         rhref_ens_f(:,:) = rhref_ens
         tref_f(:) = sum(tref_ens(1:DEF_DA_ENS_NUM,:), dim=1) / DEF_DA_ENS_NUM
         qref_f(:) = sum(qref_ens(1:DEF_DA_ENS_NUM,:), dim=1) / DEF_DA_ENS_NUM
      ENDIF

      ! 2 m temperature analysis: loop over patches and perform LETKF for
      ! each patch with the temperature observations only
      IF (use_2mt .and. p_is_worker .and. DEF_DA_ENS_NUM > 1) THEN
         DO np = 1, numpatch
            IF (patchtype(np) >= 3) CYCLE

            ! prepare localized y, R, and H(x) for the target patch
            CALL tq_t%prepare(np, nobs_p, y_p, hx_p, r_p)
            IF (nobs_p == 0) CYCLE

            ! perform data assimilation
            allocate(trans(DEF_DA_ENS_NUM, DEF_DA_ENS_NUM))
            CALL letkf(DEF_DA_ENS_NUM, nobs_p, hx_p, y_p, r_p, tq_t%infl, trans)

            ! apply the LETKF transform to the 2 m temperature ensemble
            allocate(tref_work_f(1, DEF_DA_ENS_NUM))
            allocate(tref_work_a(1, DEF_DA_ENS_NUM))
            DO iens = 1, DEF_DA_ENS_NUM
               tref_work_f(1,iens) = tref_ens(iens,np)
            ENDDO

            CALL dgemm('N', 'N', 1, DEF_DA_ENS_NUM, DEF_DA_ENS_NUM, &
                        1.0_r8, tref_work_f, 1, trans, DEF_DA_ENS_NUM, &
                        0.0_r8, tref_work_a, 1)

            ! TODO: confirm the post-analysis constraints before finalizing.
            ! The +/-20 K temperature increment is a legacy heuristic from the
            ! soil-temperature DA, not a LETKF or 2 m temperature requirement.
            DO iens = 1, DEF_DA_ENS_NUM
               ! limit the temperature increment to keep the analysis physical
               tref_ens(iens,np) = tref_ens(iens,np) + &
                  max(-20.0_r8, min(20.0_r8, tref_work_a(1,iens) - tref_ens(iens,np)))
            ENDDO

            deallocate(y_p, r_p, hx_p, trans)
            deallocate(tref_work_f, tref_work_a)
         ENDDO
      ENDIF

      ! 2 m humidity analysis: H(x) is relative humidity derived from
      ! qref_ens, so it is calculated after the temperature analysis;
      ! The humidity operator in calc_obs therefore uses the analyzed tref_ens
      IF (use_2mq) THEN
         CALL tq_q%calc_obs(idate, deltim)
         CALL tq_q%save(.false.)
      ENDIF

      IF (use_2mq .and. p_is_worker .and. DEF_DA_ENS_NUM > 1) THEN
         DO np = 1, numpatch
            IF (patchtype(np) >= 3) CYCLE

            ! prepare localized y, R, and H(x) for the target patch
            CALL tq_q%prepare(np, nobs_p, y_p, hx_p, r_p)
            IF (nobs_p == 0) CYCLE

            ! perform data assimilation
            allocate(trans(DEF_DA_ENS_NUM, DEF_DA_ENS_NUM))
            CALL letkf(DEF_DA_ENS_NUM, nobs_p, hx_p, y_p, r_p, tq_q%infl, trans)

            ! apply the LETKF transform to the 2 m specific humidity ensemble
            allocate(qref_work_f(1, DEF_DA_ENS_NUM))
            allocate(qref_work_a(1, DEF_DA_ENS_NUM))
            DO iens = 1, DEF_DA_ENS_NUM
               qref_work_f(1,iens) = qref_ens(iens,np)
            ENDDO

            CALL dgemm('N', 'N', 1, DEF_DA_ENS_NUM, DEF_DA_ENS_NUM, &
                        1.0_r8, qref_work_f, 1, trans, DEF_DA_ENS_NUM, &
                        0.0_r8, qref_work_a, 1)

            ! TODO: confirm the post-analysis constraints before finalizing.
            ! The q in [0, qsat] clipping enforces physical consistency, but
            ! it makes the analysis depart from the direct LETKF solution.
            DO iens = 1, DEF_DA_ENS_NUM
               ! keep specific humidity within [0, qsat]
               qref_ens(iens,np) = max(0.0_r8, qref_work_a(1,iens))
               IF (tref_ens(iens,np) /= spval .and. &
                   forc_psrf(np)    /= spval .and. &
                   forc_psrf(np)    >  0.0_r8) THEN
                  CALL qsadv(tref_ens(iens,np), forc_psrf(np), es, esdT, qs, qsdT)
                  qref_ens(iens,np) = min(qref_ens(iens,np), qs)
               ENDIF
            ENDDO

            deallocate(y_p, r_p, hx_p, trans)
            deallocate(qref_work_f, qref_work_a)
         ENDDO
      ENDIF

      ! Publish the completed analyses and save model-space H(x) diagnostics.
      ! Disabled targets retain their forecast values.
      IF (DEF_DA_ENS_NUM > 1) THEN
         IF (p_is_worker .and. numpatch > 0) THEN
            ! Compute relative humidity from the final analyzed ensembles.
            DO np = 1, numpatch
               DO iens = 0, DEF_DA_ENS_NUM
                  rhref_ens(iens,np) = spval

                  IF (tref_ens(iens,np) == spval .or. &
                      qref_ens(iens,np) == spval .or. &
                      forc_psrf(np)    == spval .or. &
                      forc_psrf(np)    <= 0.0_r8) CYCLE

                  tc_rh = tref_ens(iens,np) - 273.15_r8
                  es_rh = 6.112_r8 * exp((17.67_r8*tc_rh)/(tc_rh + 243.5_r8))
                  e_rh  = (qref_ens(iens,np)*forc_psrf(np)/100.0_r8) / &
                       (0.378_r8*qref_ens(iens,np) + 0.622_r8)
                  IF (es_rh > 0.0_r8) THEN
                     rhref_ens(iens,np) = max(min(1.0_r8, e_rh/es_rh), 0.0_r8)
                  ENDIF
               ENDDO
            ENDDO

            ! Store final ensembles after temperature and saturation constraints.
            allocate (tref_ens_a (0:DEF_DA_ENS_NUM, numpatch))
            allocate (tref_a     (numpatch))
            allocate (qref_ens_a (0:DEF_DA_ENS_NUM, numpatch))
            allocate (qref_a     (numpatch))
            allocate (rhref_ens_a(0:DEF_DA_ENS_NUM, numpatch))
            tref_ens_a (:,:) = tref_ens
            qref_ens_a (:,:) = qref_ens
            rhref_ens_a(:,:) = rhref_ens
            IF (use_2mt) THEN
               tref_a(:) = sum(tref_ens(1:DEF_DA_ENS_NUM,:), dim=1) / DEF_DA_ENS_NUM
            ELSE
               tref_a(:) = tref_f(:)
            ENDIF
            IF (use_2mq) THEN
               qref_a(:) = sum(qref_ens(1:DEF_DA_ENS_NUM,:), dim=1) / DEF_DA_ENS_NUM
            ELSE
               qref_a(:) = qref_f(:)
            ENDIF
         ENDIF
         IF (use_2mt) THEN
            CALL tq_t%calc_model()
            CALL tq_t%save(.true.)
         ENDIF
         IF (use_2mq) THEN
            CALL tq_q%calc_model()
            CALL tq_q%save(.true.)
         ENDIF
      ENDIF
      IF (use_2mt) CALL tq_t%clear()
      IF (use_2mq) CALL tq_q%clear()

   END SUBROUTINE run_DA_2MTQ

!-----------------------------------------------------------------------------

   SUBROUTINE end_DA_2MTQ ()

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!-----------------------------------------------------------------------------
      IF (use_2mt) THEN
         CALL tq_t%clear()
         IF (allocated(tq_t%cfgs))    deallocate(tq_t%cfgs)
         IF (allocated(tq_t%synop%cfgs)) deallocate(tq_t%synop%cfgs)
         IF (allocated(tq_t%synop%sensors)) deallocate(tq_t%synop%sensors)
      ENDIF
      IF (use_2mq) THEN
         CALL tq_q%clear()
         IF (allocated(tq_q%cfgs))    deallocate(tq_q%cfgs)
         IF (allocated(tq_q%synop%cfgs)) deallocate(tq_q%synop%cfgs)
         IF (allocated(tq_q%synop%sensors)) deallocate(tq_q%synop%sensors)
      ENDIF

      IF (allocated(tref_ens_f )) deallocate (tref_ens_f )
      IF (allocated(tref_ens_a )) deallocate (tref_ens_a )
      IF (allocated(tref_f     )) deallocate (tref_f     )
      IF (allocated(tref_a     )) deallocate (tref_a     )
      IF (allocated(qref_ens_f )) deallocate (qref_ens_f )
      IF (allocated(qref_ens_a )) deallocate (qref_ens_a )
      IF (allocated(qref_f     )) deallocate (qref_f     )
      IF (allocated(qref_a     )) deallocate (qref_a     )
      IF (allocated(rhref_ens_f)) deallocate (rhref_ens_f)
      IF (allocated(rhref_ens_a)) deallocate (rhref_ens_a)

      use_2mt = .false.
      use_2mq = .false.

   END SUBROUTINE end_DA_2MTQ

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Assim_2MTQ
#endif
