#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Assim_ST
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Main soil-temperature assimilation driver for multiple observation
!    streams.
!
!    The forward operator is the same passive-microwave RTM used by soil
!    moisture: the 'ST' streams reuse PM_set_type, and PM_RTM already takes
!    t_soisno_ens as input. Only the state update differs, so run_Assim_ST
!    updates the top two soil layers of t_soisno_ens instead of wliq_soisno_ens.
!
! HISTORY:
!    Lu Li, 05/2026: First implementation
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_DA_Obs_ST
   USE MOD_DA_EnKF
   USE MOD_DA_Vars_TimeVariables
   USE MOD_Vars_TimeInvariants
   USE MOD_Vars_Global, ONLY: nl_soil
   USE MOD_LandPatch
   USE MOD_Namelist
   USE MOD_SPMD_Task
   IMPLICIT NONE

   type(ST_type) :: st

   PUBLIC :: init_assim_ST
   PUBLIC :: run_assim_ST
   PUBLIC :: end_assim_ST

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE init_assim_ST ()

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

   END SUBROUTINE init_assim_ST

!-----------------------------------------------------------------------------

   SUBROUTINE run_assim_ST (idate, deltim)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      integer,  intent(in) :: idate(3)
      real(r8), intent(in) :: deltim

!------------------------ Local Variables ------------------------------------
      integer :: np
      integer :: iens
      integer :: nobs_p
      real(r8), allocatable :: y_p(:)
      real(r8), allocatable :: hx_p(:,:)
      real(r8), allocatable :: r_p(:)
      real(r8), allocatable :: trans(:,:)
      real(r8), allocatable :: tsoi_ens_f(:,:)
      real(r8), allocatable :: tsoi_ens_a(:,:)

!-----------------------------------------------------------------------------
      ! calculate y and H(x) at y locations from different sources
      CALL st%calcg(idate, deltim)
      CALL st%save(.false.)

      IF (p_is_worker .and. numpatch > 0) THEN
         t_soisno_ol(:,:) = t_soisno_ens(:,0,:)
         t_soisno_f (:,:) = sum(t_soisno_ens(:,1:DEF_DA_ENS_NUM,:), dim=2) / DEF_DA_ENS_NUM
      ENDIF

      ! loop over patches and perform assimilation for each patch
      IF (p_is_worker .and. DEF_DA_ENS_NUM > 1) THEN
         DO np = 1, numpatch
            IF (patchtype(np) >= 3) CYCLE

            ! prepare localized y, R, and H(x) for the target patch
            CALL st%prepare(np, nobs_p, y_p, hx_p, r_p)
            IF (nobs_p == 0) CYCLE

            ! perform data assimilation
            allocate(trans(DEF_DA_ENS_NUM, DEF_DA_ENS_NUM))
            CALL letkf(DEF_DA_ENS_NUM, nobs_p, hx_p, y_p, r_p, st%infl, trans)

            allocate(tsoi_ens_f(nl_soil, DEF_DA_ENS_NUM))
            allocate(tsoi_ens_a(nl_soil, DEF_DA_ENS_NUM))
            DO iens = 1, DEF_DA_ENS_NUM
               tsoi_ens_f(:,iens) = t_soisno_ens(1:,iens,np)
            ENDDO
            CALL dgemm('N', 'N', nl_soil, DEF_DA_ENS_NUM, DEF_DA_ENS_NUM, &
                        1.0_r8, tsoi_ens_f, nl_soil, trans, DEF_DA_ENS_NUM, &
                        0.0_r8, tsoi_ens_a, nl_soil)

            ! limit the temperature increment to keep the analysis physical
            DO iens = 1, DEF_DA_ENS_NUM
               t_soisno_ens(1:2,iens,np) = t_soisno_ens(1:2,iens,np) + &
                  max(-20.0_r8, min(20.0_r8, tsoi_ens_a(1:2,iens) - t_soisno_ens(1:2,iens,np)))
            ENDDO

            deallocate(y_p, r_p, hx_p, trans)
            deallocate(tsoi_ens_f, tsoi_ens_a)
         ENDDO
      ENDIF

      ! Recalculate H(x) after LETKF so history can compare prior and analysis
      IF (DEF_DA_ENS_NUM > 1) THEN
         IF (p_is_worker .and. numpatch > 0) THEN
            t_soisno_a(:,:) = sum(t_soisno_ens(:,1:DEF_DA_ENS_NUM,:), dim=2) / DEF_DA_ENS_NUM
         ENDIF
         CALL st%calcp()
         CALL st%save(.true.)
      ENDIF
      CALL st%clear()

   END SUBROUTINE run_assim_ST

!-----------------------------------------------------------------------------

   SUBROUTINE end_assim_ST()

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!-----------------------------------------------------------------------------
      CALL st%clear()
      IF (allocated(st%cfgs)) deallocate(st%cfgs)
      IF (allocated(st%pm%sensors)) deallocate(st%pm%sensors)
      IF (allocated(st%pm%cfgs)) deallocate(st%pm%cfgs)

      st%pm%nsensor = 0
      st%use_pm = .false.

   END SUBROUTINE end_assim_ST

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Assim_ST
#endif
