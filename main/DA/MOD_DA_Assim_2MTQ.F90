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
   USE MOD_DA_Obs_SYNOP
   USE MOD_DA_EnKF
   USE MOD_DA_Vars_TimeVariables
   USE MOD_LandPatch
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_Vars_Global, ONLY: spval
   USE MOD_Vars_1DForcing, ONLY: forc_psrf
   USE MOD_Qsadv, ONLY: qsadv
   IMPLICIT NONE

   type(SYNOP_set_type) :: synop_t   ! '2MT' streams: 2 m temperature
   type(SYNOP_set_type) :: synop_q   ! '2MQ' streams: 2 m relative humidity

   logical :: use_2mt = .false.      ! any '2MT' stream configured
   logical :: use_2mq = .false.      ! any '2MQ' stream configured

   PUBLIC :: init_assim_2MTQ
   PUBLIC :: run_assim_2MTQ
   PUBLIC :: end_assim_2MTQ

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE init_assim_2MTQ ()

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      integer :: is, i, nsensor_t, nsensor_q
      type(SYNOP_config_type), allocatable :: configs(:)

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

      ! Each SYNOP stream is packed with the variable enabled by its target.
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
               configs(is)%use_t       = .true.
            ENDIF
         ENDDO
         CALL synop_t%init(landpatch, configs)
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
               configs(is)%use_rh      = .true.
            ENDIF
         ENDDO
         CALL synop_q%init(landpatch, configs)
         deallocate(configs)
      ENDIF

   END SUBROUTINE init_assim_2MTQ

!-----------------------------------------------------------------------------

   SUBROUTINE run_assim_2MTQ (idate, deltim)

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
      real(r8), allocatable :: y_p(:)
      real(r8), allocatable :: hx_p(:,:)
      real(r8), allocatable :: r_p(:)
      real(r8), allocatable :: trans(:,:)
      real(r8), allocatable :: tref_work_f(:,:)
      real(r8), allocatable :: tref_work_a(:,:)
      real(r8), allocatable :: qref_work_f(:,:)
      real(r8), allocatable :: qref_work_a(:,:)

!-----------------------------------------------------------------------------
      ! read the SYNOP observations and calculate the prior H(x) at station
      ! sites for each enabled target
      IF (use_2mt) THEN
         CALL synop_t%calcg(idate, deltim)
         CALL synop_t%save(.false.)
      ENDIF

      IF (p_is_worker .and. numpatch > 0) THEN
         ! Preserve the complete prior ensemble before either screen-level analysis.
         tref_ens_f(:,:) = tref_ens
         qref_ens_f(:,:) = qref_ens
         tref_ens_a(:,:) = tref_ens_f
         qref_ens_a(:,:) = qref_ens_f
         tref_f(:) = sum(tref_ens(1:DEF_DA_ENS_NUM,:), dim=1) / DEF_DA_ENS_NUM
         qref_f(:) = sum(qref_ens(1:DEF_DA_ENS_NUM,:), dim=1) / DEF_DA_ENS_NUM
      ENDIF

      ! 2 m temperature analysis: loop over patches and perform LETKF for
      ! each patch with the temperature observations only
      IF (use_2mt .and. p_is_worker .and. DEF_DA_ENS_NUM > 1) THEN
         DO np = 1, numpatch
            IF (patchtype(np) >= 3) CYCLE

            ! prepare localized y, R, and H(x) for the target patch
            CALL synop_t%prepare(np, nobs_p, y_p, hx_p, r_p)
            IF (nobs_p == 0) CYCLE

            ! perform data assimilation
            allocate(trans(DEF_DA_ENS_NUM, DEF_DA_ENS_NUM))
            CALL letkf(DEF_DA_ENS_NUM, nobs_p, hx_p, y_p, r_p, synop_t%infl, trans)

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
      ! SYNOP_calc_rh inside calcg therefore sees the analyzed tref_ens
      IF (use_2mq) THEN
         CALL synop_q%calcg(idate, deltim)
         CALL synop_q%save(.false.)
      ENDIF

      IF (use_2mq .and. p_is_worker .and. DEF_DA_ENS_NUM > 1) THEN
         DO np = 1, numpatch
            IF (patchtype(np) >= 3) CYCLE

            ! prepare localized y, R, and H(x) for the target patch
            CALL synop_q%prepare(np, nobs_p, y_p, hx_p, r_p)
            IF (nobs_p == 0) CYCLE

            ! perform data assimilation
            allocate(trans(DEF_DA_ENS_NUM, DEF_DA_ENS_NUM))
            CALL letkf(DEF_DA_ENS_NUM, nobs_p, hx_p, y_p, r_p, synop_q%infl, trans)

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

      ! update H(x) for the analyzed ensembles and save the analysis
      ! diagnostics. The means of a disabled target stay at the prior means,
      ! so the ANA2MTQ products are always defined.
      IF (DEF_DA_ENS_NUM > 1) THEN
         IF (p_is_worker .and. numpatch > 0) THEN
            ! Store the final ensembles after temperature and saturation constraints.
            tref_ens_a(:,:) = tref_ens
            qref_ens_a(:,:) = qref_ens
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
         IF (use_2mt) CALL synop_t%save(.true.)
         IF (use_2mq) THEN
            CALL SYNOP_calc_rh()
            CALL synop_q%save(.true.)
         ENDIF
      ENDIF
      IF (use_2mt) CALL synop_t%clear()
      IF (use_2mq) CALL synop_q%clear()

   END SUBROUTINE run_assim_2MTQ

!-----------------------------------------------------------------------------

   SUBROUTINE end_assim_2MTQ ()

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!-----------------------------------------------------------------------------
      IF (use_2mt) THEN
         CALL synop_t%clear()
         IF (allocated(synop_t%cfgs))    deallocate(synop_t%cfgs)
         IF (allocated(synop_t%sensors)) deallocate(synop_t%sensors)
      ENDIF
      IF (use_2mq) THEN
         CALL synop_q%clear()
         IF (allocated(synop_q%cfgs))    deallocate(synop_q%cfgs)
         IF (allocated(synop_q%sensors)) deallocate(synop_q%sensors)
      ENDIF

      use_2mt = .false.
      use_2mq = .false.

   END SUBROUTINE end_assim_2MTQ

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Assim_2MTQ
#endif
