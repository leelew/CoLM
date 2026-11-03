#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Vars_1DFluxes
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Variables required for land-atmosphere coupling include restart state
!    variables and land-atmosphere fluxes.
!
!    Ensemble state variables written to DA restart files are defined and
!    managed in MOD_DA_Vars_TimeVariables. This module defines and manages
!    ensemble land-atmosphere fluxes and associated diagnostic variables.
!
! HISTORY:
!    Lu Li, 10/2026
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_Vars_Global
   USE MOD_SPMD_Task
   USE MOD_LandPatch
   IMPLICIT NONE
   SAVE

   PUBLIC :: allocate_1D_DAFluxes
   PUBLIC :: deallocate_1D_DAFluxes

!=============================================================================
! ensemble surface flux variables
!=============================================================================
   real(r8), allocatable :: taux_ens            (:,:)   ! east-west wind stress                                    [kg/(m s2)]
   real(r8), allocatable :: tauy_ens            (:,:)   ! north-south wind stress                                  [kg/(m s2)]
   real(r8), allocatable :: fsena_ens           (:,:)   ! sensible heat flux to the atmosphere                     [W/m2]
   real(r8), allocatable :: lfevpa_ens          (:,:)   ! latent heat flux to the atmosphere                       [W/m2]
   real(r8), allocatable :: fevpa_ens           (:,:)   ! evapotranspiration to the atmosphere                     [mm/s]
   real(r8), allocatable :: fsenl_ens           (:,:)   ! sensible heat flux from leaves                           [W/m2]
   real(r8), allocatable :: fevpl_ens           (:,:)   ! evaporation and transpiration from leaves                [mm/s]
   real(r8), allocatable :: etr_ens             (:,:)   ! transpiration rate                                       [mm/s]
   real(r8), allocatable :: fseng_ens           (:,:)   ! sensible heat flux from the ground                       [W/m2]
   real(r8), allocatable :: fevpg_ens           (:,:)   ! evaporation flux from the ground                         [mm/s]
   real(r8), allocatable :: fgrnd_ens           (:,:)   ! ground heat flux                                         [W/m2]
   real(r8), allocatable :: sabvsun_ens         (:,:)   ! shortwave radiation absorbed by sunlit vegetation        [W/m2]
   real(r8), allocatable :: sabvsha_ens         (:,:)   ! shortwave radiation absorbed by shaded vegetation        [W/m2]
   real(r8), allocatable :: sabg_ens            (:,:)   ! shortwave radiation absorbed by the ground               [W/m2]
   real(r8), allocatable :: sr_ens              (:,:)   ! total reflected shortwave radiation                      [W/m2]
   real(r8), allocatable :: solvd_ens           (:,:)   ! incident direct visible shortwave radiation              [W/m2]
   real(r8), allocatable :: solvi_ens           (:,:)   ! incident diffuse visible shortwave radiation             [W/m2]
   real(r8), allocatable :: solnd_ens           (:,:)   ! incident direct near-infrared shortwave radiation        [W/m2]
   real(r8), allocatable :: solni_ens           (:,:)   ! incident diffuse near-infrared shortwave radiation       [W/m2]
   real(r8), allocatable :: srvd_ens            (:,:)   ! reflected direct visible shortwave radiation             [W/m2]
   real(r8), allocatable :: srvi_ens            (:,:)   ! reflected diffuse visible shortwave radiation            [W/m2]
   real(r8), allocatable :: srnd_ens            (:,:)   ! reflected direct near-infrared shortwave radiation       [W/m2]
   real(r8), allocatable :: srni_ens            (:,:)   ! reflected diffuse near-infrared shortwave radiation      [W/m2]
   real(r8), allocatable :: solvdln_ens         (:,:)   ! incident direct visible radiation at local noon          [W/m2]
   real(r8), allocatable :: solviln_ens         (:,:)   ! incident diffuse visible radiation at local noon         [W/m2]
   real(r8), allocatable :: solndln_ens         (:,:)   ! incident direct near-infrared radiation at local noon    [W/m2]
   real(r8), allocatable :: solniln_ens         (:,:)   ! incident diffuse near-infrared radiation at local noon   [W/m2]
   real(r8), allocatable :: srvdln_ens          (:,:)   ! reflected direct visible radiation at local noon         [W/m2]
   real(r8), allocatable :: srviln_ens          (:,:)   ! reflected diffuse visible radiation at local noon        [W/m2]
   real(r8), allocatable :: srndln_ens          (:,:)   ! reflected direct near-infrared radiation at local noon   [W/m2]
   real(r8), allocatable :: srniln_ens          (:,:)   ! reflected diffuse near-infrared radiation at local noon  [W/m2]
   real(r8), allocatable :: olrg_ens            (:,:)   ! outgoing longwave radiation from canopy and ground       [W/m2]
   real(r8), allocatable :: rnet_ens            (:,:)   ! net surface radiation                                    [W/m2]
   real(r8), allocatable :: xerr_ens            (:,:)   ! water balance error                                      [mm/s]
   real(r8), allocatable :: zerr_ens            (:,:)   ! energy balance error                                     [W/m2]
   real(r8), allocatable :: frcsat_ens          (:,:)   ! fraction of saturated area                               [-]
   real(r8), allocatable :: rsur_ens            (:,:)   ! surface runoff                                           [mm/s]
   real(r8), allocatable :: rsur_se_ens         (:,:)   ! saturation-excess surface runoff                         [mm/s]
   real(r8), allocatable :: rsur_ie_ens         (:,:)   ! infiltration-excess surface runoff                       [mm/s]
   real(r8), allocatable :: rsub_ens            (:,:)   ! subsurface runoff                                        [mm/s]
   real(r8), allocatable :: rnof_ens            (:,:)   ! total runoff                                             [mm/s]
   real(r8), allocatable :: qintr_ens           (:,:)   ! interception flux                                        [mm/s]
   real(r8), allocatable :: qinfl_ens           (:,:)   ! infiltration flux                                        [mm/s]
   real(r8), allocatable :: qdrip_ens           (:,:)   ! throughfall flux                                         [mm/s]
   real(r8), allocatable :: assim_ens           (:,:)   ! canopy carbon assimilation rate                          [mol/(m2 s)]
   real(r8), allocatable :: respc_ens           (:,:)   ! respiration rate                                         [mol/(m2 s)]
   real(r8), allocatable :: qcharge_ens         (:,:)   ! groundwater recharge                                     [mm/s]
   real(r8), allocatable :: qlayer_ens          (:,:,:) ! water flux between soil layers                           [mm/s]
   real(r8), allocatable :: lake_deficit_ens    (:,:)   ! lake water deficit due to evaporation                    [mm/s]

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE allocate_1D_DAFluxes()

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            allocate (taux_ens             (0:DEF_DA_ENS_NUM,numpatch)); taux_ens             (:,:) = spval
            allocate (tauy_ens             (0:DEF_DA_ENS_NUM,numpatch)); tauy_ens             (:,:) = spval
            allocate (fsena_ens            (0:DEF_DA_ENS_NUM,numpatch)); fsena_ens            (:,:) = spval
            allocate (lfevpa_ens           (0:DEF_DA_ENS_NUM,numpatch)); lfevpa_ens           (:,:) = spval
            allocate (fevpa_ens            (0:DEF_DA_ENS_NUM,numpatch)); fevpa_ens            (:,:) = spval
            allocate (fsenl_ens            (0:DEF_DA_ENS_NUM,numpatch)); fsenl_ens            (:,:) = spval
            allocate (fevpl_ens            (0:DEF_DA_ENS_NUM,numpatch)); fevpl_ens            (:,:) = spval
            allocate (etr_ens              (0:DEF_DA_ENS_NUM,numpatch)); etr_ens              (:,:) = spval
            allocate (fseng_ens            (0:DEF_DA_ENS_NUM,numpatch)); fseng_ens            (:,:) = spval
            allocate (fevpg_ens            (0:DEF_DA_ENS_NUM,numpatch)); fevpg_ens            (:,:) = spval
            allocate (fgrnd_ens            (0:DEF_DA_ENS_NUM,numpatch)); fgrnd_ens            (:,:) = spval
            allocate (sabvsun_ens          (0:DEF_DA_ENS_NUM,numpatch)); sabvsun_ens          (:,:) = spval
            allocate (sabvsha_ens          (0:DEF_DA_ENS_NUM,numpatch)); sabvsha_ens          (:,:) = spval
            allocate (sabg_ens             (0:DEF_DA_ENS_NUM,numpatch)); sabg_ens             (:,:) = spval
            allocate (sr_ens               (0:DEF_DA_ENS_NUM,numpatch)); sr_ens               (:,:) = spval
            allocate (solvd_ens            (0:DEF_DA_ENS_NUM,numpatch)); solvd_ens            (:,:) = spval
            allocate (solvi_ens            (0:DEF_DA_ENS_NUM,numpatch)); solvi_ens            (:,:) = spval
            allocate (solnd_ens            (0:DEF_DA_ENS_NUM,numpatch)); solnd_ens            (:,:) = spval
            allocate (solni_ens            (0:DEF_DA_ENS_NUM,numpatch)); solni_ens            (:,:) = spval
            allocate (srvd_ens             (0:DEF_DA_ENS_NUM,numpatch)); srvd_ens             (:,:) = spval
            allocate (srvi_ens             (0:DEF_DA_ENS_NUM,numpatch)); srvi_ens             (:,:) = spval
            allocate (srnd_ens             (0:DEF_DA_ENS_NUM,numpatch)); srnd_ens             (:,:) = spval
            allocate (srni_ens             (0:DEF_DA_ENS_NUM,numpatch)); srni_ens             (:,:) = spval
            allocate (solvdln_ens          (0:DEF_DA_ENS_NUM,numpatch)); solvdln_ens          (:,:) = spval
            allocate (solviln_ens          (0:DEF_DA_ENS_NUM,numpatch)); solviln_ens          (:,:) = spval
            allocate (solndln_ens          (0:DEF_DA_ENS_NUM,numpatch)); solndln_ens          (:,:) = spval
            allocate (solniln_ens          (0:DEF_DA_ENS_NUM,numpatch)); solniln_ens          (:,:) = spval
            allocate (srvdln_ens           (0:DEF_DA_ENS_NUM,numpatch)); srvdln_ens           (:,:) = spval
            allocate (srviln_ens           (0:DEF_DA_ENS_NUM,numpatch)); srviln_ens           (:,:) = spval
            allocate (srndln_ens           (0:DEF_DA_ENS_NUM,numpatch)); srndln_ens           (:,:) = spval
            allocate (srniln_ens           (0:DEF_DA_ENS_NUM,numpatch)); srniln_ens           (:,:) = spval
            allocate (olrg_ens             (0:DEF_DA_ENS_NUM,numpatch)); olrg_ens             (:,:) = spval
            allocate (rnet_ens             (0:DEF_DA_ENS_NUM,numpatch)); rnet_ens             (:,:) = spval
            allocate (xerr_ens             (0:DEF_DA_ENS_NUM,numpatch)); xerr_ens             (:,:) = spval
            allocate (zerr_ens             (0:DEF_DA_ENS_NUM,numpatch)); zerr_ens             (:,:) = spval
            allocate (frcsat_ens           (0:DEF_DA_ENS_NUM,numpatch)); frcsat_ens           (:,:) = spval
            allocate (rsur_ens             (0:DEF_DA_ENS_NUM,numpatch)); rsur_ens             (:,:) = spval
            allocate (rsur_se_ens          (0:DEF_DA_ENS_NUM,numpatch)); rsur_se_ens          (:,:) = spval
            allocate (rsur_ie_ens          (0:DEF_DA_ENS_NUM,numpatch)); rsur_ie_ens          (:,:) = spval
            allocate (rsub_ens             (0:DEF_DA_ENS_NUM,numpatch)); rsub_ens             (:,:) = spval
            allocate (rnof_ens             (0:DEF_DA_ENS_NUM,numpatch)); rnof_ens             (:,:) = spval
            allocate (qintr_ens            (0:DEF_DA_ENS_NUM,numpatch)); qintr_ens            (:,:) = spval
            allocate (qinfl_ens            (0:DEF_DA_ENS_NUM,numpatch)); qinfl_ens            (:,:) = spval
            allocate (qdrip_ens            (0:DEF_DA_ENS_NUM,numpatch)); qdrip_ens            (:,:) = spval
            allocate (assim_ens            (0:DEF_DA_ENS_NUM,numpatch)); assim_ens            (:,:) = spval
            allocate (respc_ens            (0:DEF_DA_ENS_NUM,numpatch)); respc_ens            (:,:) = spval
            allocate (qcharge_ens          (0:DEF_DA_ENS_NUM,numpatch)); qcharge_ens          (:,:) = spval
            allocate (qlayer_ens (0:nl_soil,0:DEF_DA_ENS_NUM,numpatch)); qlayer_ens         (:,:,:) = spval
            allocate (lake_deficit_ens     (0:DEF_DA_ENS_NUM,numpatch)); lake_deficit_ens     (:,:) = spval
         ENDIF
      ENDIF

   END SUBROUTINE allocate_1D_DAFluxes

!-----------------------------------------------------------------------------

   SUBROUTINE deallocate_1D_DAFluxes()

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            deallocate (taux_ens                       )
            deallocate (tauy_ens                       )
            deallocate (fsena_ens                      )
            deallocate (lfevpa_ens                     )
            deallocate (fevpa_ens                      )
            deallocate (fsenl_ens                      )
            deallocate (fevpl_ens                      )
            deallocate (etr_ens                        )
            deallocate (fseng_ens                      )
            deallocate (fevpg_ens                      )
            deallocate (fgrnd_ens                      )
            deallocate (sabvsun_ens                    )
            deallocate (sabvsha_ens                    )
            deallocate (sabg_ens                       )
            deallocate (sr_ens                         )
            deallocate (solvd_ens                      )
            deallocate (solvi_ens                      )
            deallocate (solnd_ens                      )
            deallocate (solni_ens                      )
            deallocate (srvd_ens                       )
            deallocate (srvi_ens                       )
            deallocate (srnd_ens                       )
            deallocate (srni_ens                       )
            deallocate (solvdln_ens                    )
            deallocate (solviln_ens                    )
            deallocate (solndln_ens                    )
            deallocate (solniln_ens                    )
            deallocate (srvdln_ens                     )
            deallocate (srviln_ens                     )
            deallocate (srndln_ens                     )
            deallocate (srniln_ens                     )
            deallocate (olrg_ens                       )
            deallocate (rnet_ens                       )
            deallocate (xerr_ens                       )
            deallocate (zerr_ens                       )
            deallocate (frcsat_ens                     )
            deallocate (rsur_ens                       )
            deallocate (rsur_se_ens                    )
            deallocate (rsur_ie_ens                    )
            deallocate (rsub_ens                       )
            deallocate (rnof_ens                       )
            deallocate (qintr_ens                      )
            deallocate (qinfl_ens                      )
            deallocate (qdrip_ens                      )
            deallocate (assim_ens                      )
            deallocate (respc_ens                      )
            deallocate (qcharge_ens                    )
            deallocate (qlayer_ens                     )
            deallocate (lake_deficit_ens               )
         ENDIF
      ENDIF

   END SUBROUTINE deallocate_1D_DAFluxes

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Vars_1DFluxes
#endif