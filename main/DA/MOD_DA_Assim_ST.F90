#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Assim_ST
!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Soil temperature data assimilation using observations.
!
! WORKFLOW:
!   1. Init: Select soil temperature observation streams from the namelist and
!      initialize Data_ST with their configurations.
!   2. Run: Prepare observations and forecast H(x), then save diagnostics. For
!      assimilation, select nearby observations for each soil patch and apply
!      LETKF to update temperature in the upper two soil layers. Limit
!      temperature increments to +/-20 K. Recompute H(x) and save analysis
!      diagnostics.
!   3. End: Release observation data and stream configurations.
!
! HISTORY:
!   Lu Li, 05/2026: First implementation
!   Lu Li, Zhilong Fan, 11/2026: DA v1.0
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Vars_TimeInvariants
   USE MOD_Vars_Global
   USE MOD_LandPatch
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_DA_Data_ST
   USE MOD_DA_EnKF
   USE MOD_DA_Vars_TimeVars
   IMPLICIT NONE

   type(ST_type)                        :: st

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE init_DA_ST()

!-----------------------------------------------------------------------------
      integer                           :: is, i, nsensor
      type(ST_config_type), allocatable :: configs(:)
!-----------------------------------------------------------------------------
      nsensor = count(DEF_DA_OBS_TARGET == 'ST')
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

   SUBROUTINE run_DA_ST(idate, deltim)

!-----------------------------------------------------------------------------
      integer,  intent(in)              :: idate(3)
      real(r8), intent(in)              :: deltim

      integer                           :: np, iens      
      integer                           :: nobs_p        
      real(r8), allocatable             :: y_p(:)      
      real(r8), allocatable             :: hx_p(:,:)     
      real(r8), allocatable             :: r_p(:)        
      real(r8)                          :: trans(DEF_DA_ENS_NUM, DEF_DA_ENS_NUM)  
      real(r8)                          :: state_f(nl_soil, DEF_DA_ENS_NUM)       
      real(r8)                          :: state_a(nl_soil, DEF_DA_ENS_NUM)   
!-----------------------------------------------------------------------------
      ! Prepare observations and forecast H(x).
      CALL st%forward(idate, deltim)
      CALL st%save(.false.)

      IF (.not. DEF_DA) RETURN

      ! Apply the localized LETKF update to each soil patch.
      IF (p_is_worker) THEN
         DO np = 1, numpatch
            ! Save open-loop and forecast states before analysis.
            t_soisno_ol(:,np) = t_soisno_ens(:,0,np)
            t_soisno_f (:,np) = sum(t_soisno_ens(:,1:DEF_DA_ENS_NUM,np), dim=2) / DEF_DA_ENS_NUM

            IF (patchtype(np) >= 1) CYCLE

            ! Select observations and localize their error variances.
            CALL st%prepare(np, nobs_p, y_p, hx_p, r_p)
            IF (nobs_p == 0) CYCLE

            ! Compute the ensemble transform.
            CALL letkf(DEF_DA_ENS_NUM, nobs_p, hx_p, y_p, r_p, st%infl, trans)

            DO iens = 1, DEF_DA_ENS_NUM
               state_f(:,iens) = t_soisno_ens(1:,iens,np)
            ENDDO
            CALL dgemm('N', 'N', nl_soil, DEF_DA_ENS_NUM, DEF_DA_ENS_NUM, &
               1.0_r8, state_f, nl_soil, trans, DEF_DA_ENS_NUM,           &
               0.0_r8, state_a, nl_soil)

            ! Limit temperature increments in the upper two layers to +/-20 K.
            DO iens = 1, DEF_DA_ENS_NUM
               t_soisno_ens(1:2,iens,np) = t_soisno_ens(1:2,iens,np) + &
                  max(-20.0_r8, min(20.0_r8, state_a(1:2,iens) - t_soisno_ens(1:2,iens,np)))
            ENDDO

            deallocate(y_p, r_p, hx_p)
         ENDDO
      ENDIF

      ! Recompute model- and observation-space H(x), then save analysis diagnostics.
      IF (p_is_worker .and. numpatch > 0) THEN
         t_soisno_a(:,:) = sum(t_soisno_ens(:,1:DEF_DA_ENS_NUM,:), dim=2) / DEF_DA_ENS_NUM
      ENDIF
      CALL st%forward(idate, deltim)
      CALL st%save(.true.)
      CALL st%clear()

   END SUBROUTINE run_DA_ST

!-----------------------------------------------------------------------------

   SUBROUTINE end_DA_ST()

!-----------------------------------------------------------------------------
      CALL st%clear()
      IF (allocated(st%cfgs      )) deallocate(st%cfgs      )
      IF (allocated(st%pm%sensors)) deallocate(st%pm%sensors)
      IF (allocated(st%pm%cfgs   )) deallocate(st%pm%cfgs   )

      st%pm%nsensor = 0
      st%use_pm = .false.

   END SUBROUTINE end_DA_ST

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Assim_ST
#endif
