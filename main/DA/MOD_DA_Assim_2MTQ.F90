#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Assim_2MTQ
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Main 2 m temperature and humidity assimilation driver for SYNOP
!    observations.
!
!    The observation vector combines SYNOP 2 m air temperature [K] and 2 m
!    relative humidity [0-1]. The model humidity state is specific humidity,
!    so H(x) maps qref_ens to relative humidity with tref_ens and surface
!    pressure before the LETKF update (see MOD_DA_Obs_SYNOP). The analysis
!    variables are the model 2 m air temperature and specific humidity,
!    tref_ens and qref_ens.
!
!    A single LETKF transform is computed for each patch and applied directly
!    to the ensemble of tref_ens and qref_ens. The analyzed humidity is
!    constrained to [0, qsat] and the temperature increment is limited to
!    +/-20 K to keep the analysis physical.
!
!    Because tref_ens is itself an analysis variable, the analysis H(x)
!    diagnostic is available immediately after the update.
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

   type(SYNOP_set_type) :: synop

   PUBLIC :: init_assim_2MTQ
   PUBLIC :: run_assim_2MTQ
   PUBLIC :: end_assim_2MTQ

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE init_assim_2MTQ ()

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      integer :: is, i, nsensor
      type(SYNOP_config_type), allocatable :: configs(:)

!-----------------------------------------------------------------------------
      nsensor = 0
      DO i = 1, size(DEF_DA_OBS_TARGET)
         IF (trim(DEF_DA_OBS_TARGET(i)) == '2MTQ') THEN
            nsensor = nsensor + 1
         ENDIF
      ENDDO
      IF (nsensor == 0) RETURN

      allocate(configs(nsensor))
      is = 0
      DO i = 1, size(DEF_DA_OBS_TARGET)
         IF (trim(DEF_DA_OBS_TARGET(i)) == '2MTQ') THEN
            is = is + 1
            configs(is)%source_name = DEF_DA_OBS_SOURCE(i)
            configs(is)%sensor_name = DEF_DA_OBS_SENSOR(i)
            configs(is)%var_name    = DEF_DA_OBS_VAR(i)
            configs(is)%istream     = i
         ENDIF
      ENDDO

      CALL synop%init(landpatch, configs)

      deallocate(configs)

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
      real(r8), allocatable :: tref_ens_f(:,:)
      real(r8), allocatable :: tref_ens_a(:,:)
      real(r8), allocatable :: qref_ens_f(:,:)
      real(r8), allocatable :: qref_ens_a(:,:)

!-----------------------------------------------------------------------------
      ! read SYNOP temperature/humidity and calculate H(x) at station sites
      CALL synop%calcg(idate, deltim)
      CALL synop%save(.false.)

      IF (p_is_worker .and. numpatch > 0) THEN
         tref_f(:) = sum(tref_ens(1:DEF_DA_ENS_NUM,:), dim=1) / DEF_DA_ENS_NUM
         qref_f(:) = sum(qref_ens(1:DEF_DA_ENS_NUM,:), dim=1) / DEF_DA_ENS_NUM
      ENDIF

      ! loop over patches and perform assimilation for each patch
      IF (p_is_worker .and. DEF_DA_ENS_NUM > 1) THEN
         DO np = 1, numpatch
            IF (patchtype(np) >= 3) CYCLE

            ! prepare localized y, R, and H(x) for the target patch
            CALL synop%prepare(np, nobs_p, y_p, hx_p, r_p)
            IF (nobs_p == 0) CYCLE

            ! perform data assimilation
            allocate(trans(DEF_DA_ENS_NUM, DEF_DA_ENS_NUM))
            CALL letkf(DEF_DA_ENS_NUM, nobs_p, hx_p, y_p, r_p, synop%infl, trans)

            ! apply the LETKF transform to the 2 m temperature and specific
            ! humidity ensembles
            allocate(tref_ens_f(1, DEF_DA_ENS_NUM))
            allocate(tref_ens_a(1, DEF_DA_ENS_NUM))
            allocate(qref_ens_f(1, DEF_DA_ENS_NUM))
            allocate(qref_ens_a(1, DEF_DA_ENS_NUM))
            DO iens = 1, DEF_DA_ENS_NUM
               tref_ens_f(1,iens) = tref_ens(iens,np)
               qref_ens_f(1,iens) = qref_ens(iens,np)
            ENDDO

            CALL dgemm('N', 'N', 1, DEF_DA_ENS_NUM, DEF_DA_ENS_NUM, &
                        1.0_r8, tref_ens_f, 1, trans, DEF_DA_ENS_NUM, &
                        0.0_r8, tref_ens_a, 1)
            CALL dgemm('N', 'N', 1, DEF_DA_ENS_NUM, DEF_DA_ENS_NUM, &
                        1.0_r8, qref_ens_f, 1, trans, DEF_DA_ENS_NUM, &
                        0.0_r8, qref_ens_a, 1)

            ! TODO: confirm the post-analysis constraints before finalizing.
            ! The +/-20 K temperature increment is a legacy heuristic from the
            ! soil-temperature DA, not a LETKF or 2 m temperature requirement.
            ! The q in [0, qsat] clipping enforces physical consistency, but it
            ! makes the analysis depart from the direct LETKF solution. Decide
            ! whether to keep them, make them namelist options, or let the
            ! model's own saturation adjustment handle q.
            DO iens = 1, DEF_DA_ENS_NUM
               ! limit the temperature increment to keep the analysis physical
               tref_ens(iens,np) = tref_ens(iens,np) + &
                  max(-20.0_r8, min(20.0_r8, tref_ens_a(1,iens) - tref_ens(iens,np)))

               ! keep specific humidity within [0, qsat]
               qref_ens(iens,np) = max(0.0_r8, qref_ens_a(1,iens))
               IF (tref_ens(iens,np) /= spval .and. &
                   forc_psrf(np)    /= spval .and. &
                   forc_psrf(np)    >  0.0_r8) THEN
                  CALL qsadv(tref_ens(iens,np), forc_psrf(np), es, esdT, qs, qsdT)
                  qref_ens(iens,np) = min(qref_ens(iens,np), qs)
               ENDIF
            ENDDO

            deallocate(y_p, r_p, hx_p, trans)
            deallocate(tref_ens_f, tref_ens_a, qref_ens_f, qref_ens_a)
         ENDDO
      ENDIF

      ! update H(x) for the analyzed ensemble and save the analysis diagnostics
      IF (DEF_DA_ENS_NUM > 1) THEN
         IF (p_is_worker .and. numpatch > 0) THEN
            tref_a(:) = sum(tref_ens(1:DEF_DA_ENS_NUM,:), dim=1) / DEF_DA_ENS_NUM
            qref_a(:) = sum(qref_ens(1:DEF_DA_ENS_NUM,:), dim=1) / DEF_DA_ENS_NUM
         ENDIF
         CALL SYNOP_calc_rh()
         CALL synop%save(.true.)
      ENDIF
      CALL synop%clear()

   END SUBROUTINE run_assim_2MTQ

!-----------------------------------------------------------------------------

   SUBROUTINE end_assim_2MTQ ()

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!-----------------------------------------------------------------------------
      CALL synop%clear()
      IF (allocated(synop%cfgs))    deallocate(synop%cfgs)
      IF (allocated(synop%sensors)) deallocate(synop%sensors)

   END SUBROUTINE end_assim_2MTQ

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Assim_2MTQ
#endif
