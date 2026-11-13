#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Vars_TimeVars
!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Defines and manages the time-varying model variables required by the data
!   Assimilation workflow.
!
!   This module stores the DA ensemble copies of prognostic state variables,
!   Land-atmosphere coupling variables, history variables, and perturbed
!   Forcing variables
!
!   Ensemble copies include an unperturbed open-loop trajectory and a set
!   Of perturbed members. The open-loop trajectory is retained in memory
!   And restart files for restoration and diagnostics. Only the perturbed
!   Members are used in ensemble assimilation.
!
! HISTORY:
!   Lu Li, 12/2024, 05/2025: Initial DA ensemble storage
!   Lu Li, 05/2026: Store open-loop member as ensemble index 0
!   Lu Li, Zhilong Fan, 11/2026: DA v1.0
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_Vars_Global
   USE MOD_SPMD_Task
   USE MOD_LandPatch
   USE MOD_NetCDFVector
#ifdef GridRiverLakeFlow
   USE MOD_DataType
   USE MOD_NetCDFSerial
   USE MOD_Vector_ReadWrite
   USE MOD_Grid_RiverLakeNetwork
   USE MOD_Grid_Reservoir
#endif
#ifdef RangeCheck
   USE MOD_RangeCheck
#endif
   IMPLICIT NONE
   SAVE

!=============================================================================
! Ensemble restart state variables
!=============================================================================
   real(r8), allocatable :: z_sno_ens        (:,:,:  )  ! Layer depth                                              [m]
   real(r8), allocatable :: dz_sno_ens       (:,:,:  )  ! Snow layer thickness                                     [m]
   real(r8), allocatable :: t_soisno_ens     (:,:,:  )  ! Snow & soil temperature                                  [K]
   real(r8), allocatable :: wliq_soisno_ens  (:,:,:  )  ! Liquid water in layers                                   [kg/m2]
   real(r8), allocatable :: wice_soisno_ens  (:,:,:  )  ! Ice in layers                                            [kg/m2]
   real(r8), allocatable :: smp_ens          (:,:,:  )  ! Soil matrix potential                                    [mm]
   real(r8), allocatable :: hk_ens           (:,:,:  )  ! Hydraulic conductivity                                   [mm h2o/s]
   real(r8), allocatable :: t_grnd_ens       (:,:    )  ! Ground surface temperature                               [K]
   real(r8), allocatable :: tleaf_ens        (:,:    )  ! Leaf temperature                                         [K]
   real(r8), allocatable :: ldew_ens         (:,:    )  ! Depth of water on foliage                                [mm]
   real(r8), allocatable :: ldew_rain_ens    (:,:    )  ! Depth of rain on foliage                                 [mm]
   real(r8), allocatable :: ldew_snow_ens    (:,:    )  ! Depth of snow on foliage                                 [mm]
   real(r8), allocatable :: fwet_snow_ens    (:,:    )  ! Vegetation snow fractional cover                         [-]
   real(r8), allocatable :: sag_ens          (:,:    )  ! Non dimensional snow age                                 [-]
   real(r8), allocatable :: scv_ens          (:,:    )  ! Snow water equivalent                                    [kg/m2]
   real(r8), allocatable :: snowdp_ens       (:,:    )  ! Snow depth                                               [m]
   real(r8), allocatable :: fveg_ens         (:,:    )  ! Fraction vegetation cover                                [-]
   real(r8), allocatable :: fsno_ens         (:,:    )  ! Fraction snow cover                                      [-]
   real(r8), allocatable :: sigf_ens         (:,:    )  ! Fraction of veg cover, excluding snow-covered veg        [-]
   real(r8), allocatable :: green_ens        (:,:    )  ! Leaf greenness                                           [-]
   real(r8), allocatable :: tlai_ens         (:,:    )  ! Leaf area index                                          [m2/m2]
   real(r8), allocatable :: lai_ens          (:,:    )  ! Leaf area index                                          [m2/m2]
   real(r8), allocatable :: tsai_ens         (:,:    )  ! Stem area index                                          [m2/m2]
   real(r8), allocatable :: sai_ens          (:,:    )  ! Stem area index                                          [m2/m2]
   real(r8), allocatable :: alb_ens          (:,:,:,:)  ! Averaged albedo                                          [-]
   real(r8), allocatable :: ssun_ens         (:,:,:,:)  ! Sunlit canopy absorption for solar radiation (0-1)       [-]
   real(r8), allocatable :: ssha_ens         (:,:,:,:)  ! Shaded canopy absorption for solar radiation (0-1)       [-]
   real(r8), allocatable :: ssoi_ens         (:,:,:,:)  ! Soil absorption for solar radiation (0-1)                [-]
   real(r8), allocatable :: ssno_ens         (:,:,:,:)  ! Snow absorption for solar radiation (0-1)                [-]
   real(r8), allocatable :: thermk_ens       (:,:    )  ! Canopy gap fraction for tir radiation                    [-]
   real(r8), allocatable :: extkb_ens        (:,:    )  ! (k, g(mu)/mu) direct solar extinction coefficient        [-]
   real(r8), allocatable :: extkd_ens        (:,:    )  ! Diffuse and scattered diffuse PAR extinction coefficient [-]
   real(r8), allocatable :: zwt_ens          (:,:    )  ! The depth to water table                                 [m]
   real(r8), allocatable :: wdsrf_ens        (:,:    )  ! Depth of surface water                                   [mm]
   real(r8), allocatable :: wa_ens           (:,:    )  ! Water storage in aquifer                                 [mm]
   real(r8), allocatable :: wetwat_ens       (:,:    )  ! Water storage in wetland                                 [mm]
   real(r8), allocatable :: t_lake_ens       (:,:,:  )  ! Lake layer temperature                                   [K]
   real(r8), allocatable :: lake_icefrac_ens (:,:,:  )  ! Lake mass fraction of lake layer that is frozen          [-]
   real(r8), allocatable :: savedtke1_ens    (:,:    )  ! Top level eddy conductivity                              [W/(m K)]

!=============================================================================
! Ensemble land-atmosphere coupling state variables (for online mode)
!=============================================================================
   real(r8), allocatable :: trad_ens         (:,:    )  ! Surface radiative temperature                            [K]
   real(r8), allocatable :: tref_ens         (:,:    )  ! 2 m air temperature                                      [K]
   real(r8), allocatable :: qref_ens         (:,:    )  ! 2 m specific humidity                                    [kg/kg]
   real(r8), allocatable :: rhref_ens        (:,:    )  ! 2 m relative humidity                                    [-]
   real(r8), allocatable :: rst_ens          (:,:    )  ! Canopy stomatal resistance                               [s/m]
   real(r8), allocatable :: emis_ens         (:,:    )  ! Bulk surface emissivity                                  [-]
   real(r8), allocatable :: z0m_ens          (:,:    )  ! Effective momentum roughness length                      [m]
   real(r8), allocatable :: zol_ens          (:,:    )  ! Dimensionless height in Monin-Obukhov theory (z/L)       [-]
   real(r8), allocatable :: rib_ens          (:,:    )  ! Bulk Richardson number in the surface layer              [-]
   real(r8), allocatable :: ustar_ens        (:,:    )  ! Friction velocity in similarity theory                   [m/s]
   real(r8), allocatable :: qstar_ens        (:,:    )  ! Humidity scale in similarity theory                      [kg/kg]
   real(r8), allocatable :: tstar_ens        (:,:    )  ! Temperature scale in similarity theory                   [K]
   real(r8), allocatable :: fm_ens           (:,:    )  ! Integrated momentum profile function                     [-]
   real(r8), allocatable :: fh_ens           (:,:    )  ! Integrated heat profile function                         [-]
   real(r8), allocatable :: fq_ens           (:,:    )  ! Integrated moisture profile function                     [-]

!=============================================================================
! Ensemble forcing variables
!=============================================================================
   real(r8), allocatable :: forc_t_ens       (:,:    )  ! Temperature                                              [K]
   real(r8), allocatable :: forc_frl_ens     (:,:    )  ! Atmospheric infrared (longwave) radiation                [W/m2]
   real(r8), allocatable :: forc_prc_ens     (:,:    )  ! Convective precipitation                                 [mm/s]
   real(r8), allocatable :: forc_prl_ens     (:,:    )  ! Large scale precipitation                                [mm/s]
   real(r8), allocatable :: forc_sols_ens    (:,:    )  ! Atm vis direct beam solar rad onto srf                   [W/m2]
   real(r8), allocatable :: forc_soll_ens    (:,:    )  ! Atm nir direct beam solar rad onto srf                   [W/m2]
   real(r8), allocatable :: forc_solsd_ens   (:,:    )  ! Atm vis diffuse solar rad onto srf                       [W/m2]
   real(r8), allocatable :: forc_solld_ens   (:,:    )  ! Atm nir diffuse solar rad onto srf                       [W/m2]

!=============================================================================
! Diagnostic variables for output history
!=============================================================================
   integer               :: nstream = 0                 ! Number of configured observation streams
   real(r8), allocatable :: hx_f_ens         (:,:,:  )  ! Prior H(x_f)
   real(r8), allocatable :: hx_a_ens         (:,:,:  )  ! Posterior H(x_a)
   real(r8), allocatable :: hx_ol            (:,:    )  ! Open-loop H(x)
   real(r8), allocatable :: hx_f             (:,:    )  ! Ensemble-mean prior H(x_f)
   real(r8), allocatable :: hx_a             (:,:    )  ! Ensemble-mean posterior H(x_a)
   real(r8), allocatable :: wliq_soisno_ol   (:,:    )  ! Open-loop liquid water in layers                         [kg/m2]
   real(r8), allocatable :: wliq_soisno_f    (:,:    )  ! Ensemble-mean prior liquid water in layers               [kg/m2]
   real(r8), allocatable :: wliq_soisno_a    (:,:    )  ! Ensemble-mean posterior liquid water in layers           [kg/m2]
   real(r8), allocatable :: t_soisno_ol      (:,:    )  ! Open-loop temperature in layers                          [K]
   real(r8), allocatable :: t_soisno_f       (:,:    )  ! Ensemble-mean prior temperature in layers                [K]
   real(r8), allocatable :: t_soisno_a       (:,:    )  ! Ensemble-mean posterior temperature in layers            [K]
#ifdef GridRiverLakeFlow
   real(r8), allocatable :: wdsrf_ucat_ol    (:      )  ! Open-loop river water depth                              [m]
   real(r8), allocatable :: wdsrf_ucat_a     (:      )  ! Analysis river water depth                               [m]
   real(r8), allocatable :: veloc_riv_ol     (:      )  ! Open-loop river velocity                                 [m/s]
   real(r8), allocatable :: veloc_riv_a      (:      )  ! Analysis river velocity                                  [m/s]
   real(r8), allocatable :: volresv_ol       (:      )  ! Open-loop reservoir water volume                         [m^3]
   real(r8), allocatable :: volresv_a        (:      )  ! Analysis reservoir water volume                          [m^3]
#endif

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE allocate_DA_TimeVars()

!-----------------------------------------------------------------------------
      nstream = count(len_trim(DEF_DA_OBS_TARGET) > 0)

      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            allocate(z_sno_ens        (  maxsnl+1:0,      0:DEF_DA_ENS_NUM,numpatch)); z_sno_ens        (:,:,:  ) = spval
            allocate(dz_sno_ens       (  maxsnl+1:0,      0:DEF_DA_ENS_NUM,numpatch)); dz_sno_ens       (:,:,:  ) = spval
            allocate(t_soisno_ens     (  maxsnl+1:nl_soil,0:DEF_DA_ENS_NUM,numpatch)); t_soisno_ens     (:,:,:  ) = spval
            allocate(wliq_soisno_ens  (  maxsnl+1:nl_soil,0:DEF_DA_ENS_NUM,numpatch)); wliq_soisno_ens  (:,:,:  ) = spval
            allocate(wice_soisno_ens  (  maxsnl+1:nl_soil,0:DEF_DA_ENS_NUM,numpatch)); wice_soisno_ens  (:,:,:  ) = spval
            allocate(smp_ens          (  1:nl_soil,       0:DEF_DA_ENS_NUM,numpatch)); smp_ens          (:,:,:  ) = spval
            allocate(hk_ens           (  1:nl_soil,       0:DEF_DA_ENS_NUM,numpatch)); hk_ens           (:,:,:  ) = spval
            allocate(t_grnd_ens       (                   0:DEF_DA_ENS_NUM,numpatch)); t_grnd_ens       (:,:    ) = spval
            allocate(tleaf_ens        (                   0:DEF_DA_ENS_NUM,numpatch)); tleaf_ens        (:,:    ) = spval
            allocate(ldew_ens         (                   0:DEF_DA_ENS_NUM,numpatch)); ldew_ens         (:,:    ) = spval
            allocate(ldew_rain_ens    (                   0:DEF_DA_ENS_NUM,numpatch)); ldew_rain_ens    (:,:    ) = spval
            allocate(ldew_snow_ens    (                   0:DEF_DA_ENS_NUM,numpatch)); ldew_snow_ens    (:,:    ) = spval
            allocate(fwet_snow_ens    (                   0:DEF_DA_ENS_NUM,numpatch)); fwet_snow_ens    (:,:    ) = spval
            allocate(sag_ens          (                   0:DEF_DA_ENS_NUM,numpatch)); sag_ens          (:,:    ) = spval
            allocate(scv_ens          (                   0:DEF_DA_ENS_NUM,numpatch)); scv_ens          (:,:    ) = spval
            allocate(snowdp_ens       (                   0:DEF_DA_ENS_NUM,numpatch)); snowdp_ens       (:,:    ) = spval
            allocate(fveg_ens         (                   0:DEF_DA_ENS_NUM,numpatch)); fveg_ens         (:,:    ) = spval
            allocate(fsno_ens         (                   0:DEF_DA_ENS_NUM,numpatch)); fsno_ens         (:,:    ) = spval
            allocate(sigf_ens         (                   0:DEF_DA_ENS_NUM,numpatch)); sigf_ens         (:,:    ) = spval
            allocate(green_ens        (                   0:DEF_DA_ENS_NUM,numpatch)); green_ens        (:,:    ) = spval
            allocate(tlai_ens         (                   0:DEF_DA_ENS_NUM,numpatch)); tlai_ens         (:,:    ) = spval
            allocate(lai_ens          (                   0:DEF_DA_ENS_NUM,numpatch)); lai_ens          (:,:    ) = spval
            allocate(tsai_ens         (                   0:DEF_DA_ENS_NUM,numpatch)); tsai_ens         (:,:    ) = spval
            allocate(sai_ens          (                   0:DEF_DA_ENS_NUM,numpatch)); sai_ens          (:,:    ) = spval
            allocate(alb_ens          (2,2,               0:DEF_DA_ENS_NUM,numpatch)); alb_ens          (:,:,:,:) = spval
            allocate(ssun_ens         (2,2,               0:DEF_DA_ENS_NUM,numpatch)); ssun_ens         (:,:,:,:) = spval
            allocate(ssha_ens         (2,2,               0:DEF_DA_ENS_NUM,numpatch)); ssha_ens         (:,:,:,:) = spval
            allocate(ssoi_ens         (2,2,               0:DEF_DA_ENS_NUM,numpatch)); ssoi_ens         (:,:,:,:) = spval
            allocate(ssno_ens         (2,2,               0:DEF_DA_ENS_NUM,numpatch)); ssno_ens         (:,:,:,:) = spval
            allocate(thermk_ens       (                   0:DEF_DA_ENS_NUM,numpatch)); thermk_ens       (:,:    ) = spval
            allocate(extkb_ens        (                   0:DEF_DA_ENS_NUM,numpatch)); extkb_ens        (:,:    ) = spval
            allocate(extkd_ens        (                   0:DEF_DA_ENS_NUM,numpatch)); extkd_ens        (:,:    ) = spval
            allocate(zwt_ens          (                   0:DEF_DA_ENS_NUM,numpatch)); zwt_ens          (:,:    ) = spval
            allocate(wdsrf_ens        (                   0:DEF_DA_ENS_NUM,numpatch)); wdsrf_ens        (:,:    ) = spval
            allocate(wa_ens           (                   0:DEF_DA_ENS_NUM,numpatch)); wa_ens           (:,:    ) = spval
            allocate(wetwat_ens       (                   0:DEF_DA_ENS_NUM,numpatch)); wetwat_ens       (:,:    ) = spval
            allocate(t_lake_ens       (  nl_lake,         0:DEF_DA_ENS_NUM,numpatch)); t_lake_ens       (:,:,:  ) = spval
            allocate(lake_icefrac_ens (  nl_lake,         0:DEF_DA_ENS_NUM,numpatch)); lake_icefrac_ens (:,:,:  ) = spval
            allocate(savedtke1_ens    (                   0:DEF_DA_ENS_NUM,numpatch)); savedtke1_ens    (:,:    ) = spval

            allocate(trad_ens         (                   0:DEF_DA_ENS_NUM,numpatch)); trad_ens         (:,:    ) = spval
            allocate(tref_ens         (                   0:DEF_DA_ENS_NUM,numpatch)); tref_ens         (:,:    ) = spval
            allocate(qref_ens         (                   0:DEF_DA_ENS_NUM,numpatch)); qref_ens         (:,:    ) = spval
            allocate(rhref_ens        (                   0:DEF_DA_ENS_NUM,numpatch)); rhref_ens        (:,:    ) = spval
            allocate(rst_ens          (                   0:DEF_DA_ENS_NUM,numpatch)); rst_ens          (:,:    ) = spval
            allocate(emis_ens         (                   0:DEF_DA_ENS_NUM,numpatch)); emis_ens         (:,:    ) = spval
            allocate(z0m_ens          (                   0:DEF_DA_ENS_NUM,numpatch)); z0m_ens          (:,:    ) = spval
            allocate(zol_ens          (                   0:DEF_DA_ENS_NUM,numpatch)); zol_ens          (:,:    ) = spval
            allocate(rib_ens          (                   0:DEF_DA_ENS_NUM,numpatch)); rib_ens          (:,:    ) = spval
            allocate(ustar_ens        (                   0:DEF_DA_ENS_NUM,numpatch)); ustar_ens        (:,:    ) = spval
            allocate(qstar_ens        (                   0:DEF_DA_ENS_NUM,numpatch)); qstar_ens        (:,:    ) = spval
            allocate(tstar_ens        (                   0:DEF_DA_ENS_NUM,numpatch)); tstar_ens        (:,:    ) = spval
            allocate(fm_ens           (                   0:DEF_DA_ENS_NUM,numpatch)); fm_ens           (:,:    ) = spval
            allocate(fh_ens           (                   0:DEF_DA_ENS_NUM,numpatch)); fh_ens           (:,:    ) = spval
            allocate(fq_ens           (                   0:DEF_DA_ENS_NUM,numpatch)); fq_ens           (:,:    ) = spval

            allocate(forc_t_ens       (                   0:DEF_DA_ENS_NUM,numpatch)); forc_t_ens       (:,:    ) = spval
            allocate(forc_frl_ens     (                   0:DEF_DA_ENS_NUM,numpatch)); forc_frl_ens     (:,:    ) = spval
            allocate(forc_prc_ens     (                   0:DEF_DA_ENS_NUM,numpatch)); forc_prc_ens     (:,:    ) = spval
            allocate(forc_prl_ens     (                   0:DEF_DA_ENS_NUM,numpatch)); forc_prl_ens     (:,:    ) = spval
            allocate(forc_sols_ens    (                   0:DEF_DA_ENS_NUM,numpatch)); forc_sols_ens    (:,:    ) = spval
            allocate(forc_soll_ens    (                   0:DEF_DA_ENS_NUM,numpatch)); forc_soll_ens    (:,:    ) = spval
            allocate(forc_solsd_ens   (                   0:DEF_DA_ENS_NUM,numpatch)); forc_solsd_ens   (:,:    ) = spval
            allocate(forc_solld_ens   (                   0:DEF_DA_ENS_NUM,numpatch)); forc_solld_ens   (:,:    ) = spval

            allocate(hx_f_ens         (  nstream,         DEF_DA_ENS_NUM,  numpatch)); hx_f_ens         (:,:,:  ) = spval
            allocate(hx_a_ens         (  nstream,         DEF_DA_ENS_NUM,  numpatch)); hx_a_ens         (:,:,:  ) = spval
            allocate(hx_ol            (                   nstream,         numpatch)); hx_ol            (:,:    ) = spval
            allocate(hx_f             (                   nstream,         numpatch)); hx_f             (:,:    ) = spval
            allocate(hx_a             (                   nstream,         numpatch)); hx_a             (:,:    ) = spval
            allocate(wliq_soisno_ol   (                   maxsnl+1:nl_soil,numpatch)); wliq_soisno_ol   (:,:    ) = spval
            allocate(wliq_soisno_f    (                   maxsnl+1:nl_soil,numpatch)); wliq_soisno_f    (:,:    ) = spval
            allocate(wliq_soisno_a    (                   maxsnl+1:nl_soil,numpatch)); wliq_soisno_a    (:,:    ) = spval
            allocate(t_soisno_ol      (                   maxsnl+1:nl_soil,numpatch)); t_soisno_ol      (:,:    ) = spval
            allocate(t_soisno_f       (                   maxsnl+1:nl_soil,numpatch)); t_soisno_f       (:,:    ) = spval
            allocate(t_soisno_a       (                   maxsnl+1:nl_soil,numpatch)); t_soisno_a       (:,:    ) = spval
         ENDIF
#ifdef GridRiverLakeFlow
         IF (numucat > 0) THEN
            allocate(wdsrf_ucat_ol    (                                    numucat )); wdsrf_ucat_ol    (:      ) = spval
            allocate(wdsrf_ucat_a     (                                    numucat )); wdsrf_ucat_a     (:      ) = spval
            allocate(veloc_riv_ol     (                                    numucat )); veloc_riv_ol     (:      ) = spval
            allocate(veloc_riv_a      (                                    numucat )); veloc_riv_a      (:      ) = spval
         ENDIF
         IF (numresv > 0) THEN
            allocate(volresv_ol       (                                    numresv )); volresv_ol       (:      ) = spval
            allocate(volresv_a        (                                    numresv )); volresv_a        (:      ) = spval
         ENDIF
#endif
      ENDIF

   END SUBROUTINE allocate_DA_TimeVars

!-----------------------------------------------------------------------------

   SUBROUTINE deallocate_DA_TimeVars()

!-----------------------------------------------------------------------------
      IF (allocated(z_sno_ens       )) deallocate(z_sno_ens       )
      IF (allocated(dz_sno_ens      )) deallocate(dz_sno_ens      )
      IF (allocated(t_soisno_ens    )) deallocate(t_soisno_ens    )
      IF (allocated(wliq_soisno_ens )) deallocate(wliq_soisno_ens )
      IF (allocated(wice_soisno_ens )) deallocate(wice_soisno_ens )
      IF (allocated(smp_ens         )) deallocate(smp_ens         )
      IF (allocated(hk_ens          )) deallocate(hk_ens          )
      IF (allocated(t_grnd_ens      )) deallocate(t_grnd_ens      )
      IF (allocated(tleaf_ens       )) deallocate(tleaf_ens       )
      IF (allocated(ldew_ens        )) deallocate(ldew_ens        )
      IF (allocated(ldew_rain_ens   )) deallocate(ldew_rain_ens   )
      IF (allocated(ldew_snow_ens   )) deallocate(ldew_snow_ens   )
      IF (allocated(fwet_snow_ens   )) deallocate(fwet_snow_ens   )
      IF (allocated(sag_ens         )) deallocate(sag_ens         )
      IF (allocated(scv_ens         )) deallocate(scv_ens         )
      IF (allocated(snowdp_ens      )) deallocate(snowdp_ens      )
      IF (allocated(fveg_ens        )) deallocate(fveg_ens        )
      IF (allocated(fsno_ens        )) deallocate(fsno_ens        )
      IF (allocated(sigf_ens        )) deallocate(sigf_ens        )
      IF (allocated(green_ens       )) deallocate(green_ens       )
      IF (allocated(tlai_ens        )) deallocate(tlai_ens        )
      IF (allocated(lai_ens         )) deallocate(lai_ens         )
      IF (allocated(tsai_ens        )) deallocate(tsai_ens        )
      IF (allocated(sai_ens         )) deallocate(sai_ens         )
      IF (allocated(alb_ens         )) deallocate(alb_ens         )
      IF (allocated(ssun_ens        )) deallocate(ssun_ens        )
      IF (allocated(ssha_ens        )) deallocate(ssha_ens        )
      IF (allocated(ssoi_ens        )) deallocate(ssoi_ens        )
      IF (allocated(ssno_ens        )) deallocate(ssno_ens        )
      IF (allocated(thermk_ens      )) deallocate(thermk_ens      )
      IF (allocated(extkb_ens       )) deallocate(extkb_ens       )
      IF (allocated(extkd_ens       )) deallocate(extkd_ens       )
      IF (allocated(zwt_ens         )) deallocate(zwt_ens         )
      IF (allocated(wdsrf_ens       )) deallocate(wdsrf_ens       )
      IF (allocated(wa_ens          )) deallocate(wa_ens          )
      IF (allocated(wetwat_ens      )) deallocate(wetwat_ens      )
      IF (allocated(t_lake_ens      )) deallocate(t_lake_ens      )
      IF (allocated(lake_icefrac_ens)) deallocate(lake_icefrac_ens)
      IF (allocated(savedtke1_ens   )) deallocate(savedtke1_ens   )

      IF (allocated(trad_ens        )) deallocate(trad_ens        )
      IF (allocated(tref_ens        )) deallocate(tref_ens        )
      IF (allocated(qref_ens        )) deallocate(qref_ens        )
      IF (allocated(rhref_ens       )) deallocate(rhref_ens       )
      IF (allocated(rst_ens         )) deallocate(rst_ens         )
      IF (allocated(emis_ens        )) deallocate(emis_ens        )
      IF (allocated(z0m_ens         )) deallocate(z0m_ens         )
      IF (allocated(zol_ens         )) deallocate(zol_ens         )
      IF (allocated(rib_ens         )) deallocate(rib_ens         )
      IF (allocated(ustar_ens       )) deallocate(ustar_ens       )
      IF (allocated(qstar_ens       )) deallocate(qstar_ens       )
      IF (allocated(tstar_ens       )) deallocate(tstar_ens       )
      IF (allocated(fm_ens          )) deallocate(fm_ens          )
      IF (allocated(fh_ens          )) deallocate(fh_ens          )
      IF (allocated(fq_ens          )) deallocate(fq_ens          )

      IF (allocated(forc_t_ens      )) deallocate(forc_t_ens      )
      IF (allocated(forc_frl_ens    )) deallocate(forc_frl_ens    )
      IF (allocated(forc_prc_ens    )) deallocate(forc_prc_ens    )
      IF (allocated(forc_prl_ens    )) deallocate(forc_prl_ens    )
      IF (allocated(forc_sols_ens   )) deallocate(forc_sols_ens   )
      IF (allocated(forc_soll_ens   )) deallocate(forc_soll_ens   )
      IF (allocated(forc_solsd_ens  )) deallocate(forc_solsd_ens  )
      IF (allocated(forc_solld_ens  )) deallocate(forc_solld_ens  )

      IF (allocated(hx_f_ens        )) deallocate(hx_f_ens        )
      IF (allocated(hx_a_ens        )) deallocate(hx_a_ens        )
      IF (allocated(hx_ol           )) deallocate(hx_ol           )
      IF (allocated(hx_f            )) deallocate(hx_f            )
      IF (allocated(hx_a            )) deallocate(hx_a            )
      IF (allocated(wliq_soisno_ol  )) deallocate(wliq_soisno_ol  )
      IF (allocated(wliq_soisno_f   )) deallocate(wliq_soisno_f   )
      IF (allocated(wliq_soisno_a   )) deallocate(wliq_soisno_a   )
      IF (allocated(t_soisno_ol     )) deallocate(t_soisno_ol     )
      IF (allocated(t_soisno_f      )) deallocate(t_soisno_f      )
      IF (allocated(t_soisno_a      )) deallocate(t_soisno_a      )
#ifdef GridRiverLakeFlow
      IF (allocated(wdsrf_ucat_ol   )) deallocate(wdsrf_ucat_ol   )
      IF (allocated(wdsrf_ucat_a    )) deallocate(wdsrf_ucat_a    )
      IF (allocated(veloc_riv_ol    )) deallocate(veloc_riv_ol    )
      IF (allocated(veloc_riv_a     )) deallocate(veloc_riv_a     )
      IF (allocated(volresv_ol      )) deallocate(volresv_ol      )
      IF (allocated(volresv_a       )) deallocate(volresv_a       )
#endif

   END SUBROUTINE deallocate_DA_TimeVars

!-----------------------------------------------------------------------------

   SUBROUTINE WRITE_DA_TimeVars(idate, lc_year, site, dir_restart)

!-----------------------------------------------------------------------------
      integer,          intent(in) :: idate (3)       
      integer,          intent(in) :: lc_year        
      character(len=*), intent(in) :: site            
      character(len=*), intent(in) :: dir_restart     

      character(len=256)           :: file_restart    
      character(len=14)            :: cdate           
      character(len=256)           :: cyear           
      character(len=256)           :: file_gridriver 
!-----------------------------------------------------------------------------
      write(cyear,'(i4.4)') lc_year
      write(cdate,'(i4.4,"-",i3.3,"-",i5.5)') idate(1), idate(2), idate(3)

      IF (p_is_master) CALL system('mkdir -p ' // trim(dir_restart)//'/'//trim(cdate))
#ifdef USEMPI
      CALL mpi_barrier(p_comm_glb, p_err)
#endif

      file_restart = trim(dir_restart)// '/'//trim(cdate)//'/' // trim(site) //'_restart_DA_'//trim(cdate)//'_lc'//trim(cyear)//'.nc'
      CALL ncio_create_file_vector(file_restart, landpatch                              )
      CALL ncio_define_dimension_vector(file_restart, landpatch, 'patch'                     )
      CALL ncio_define_dimension_vector(file_restart, landpatch, 'snow',     -maxsnl         )
      CALL ncio_define_dimension_vector(file_restart, landpatch, 'snowp1',   -maxsnl+1       )
      CALL ncio_define_dimension_vector(file_restart, landpatch, 'soilsnow', nl_soil-maxsnl  )
      CALL ncio_define_dimension_vector(file_restart, landpatch, 'soil',     nl_soil         )
      CALL ncio_define_dimension_vector(file_restart, landpatch, 'lake',     nl_lake         )
      CALL ncio_define_dimension_vector(file_restart, landpatch, 'band',     2               )
      CALL ncio_define_dimension_vector(file_restart, landpatch, 'rtyp',     2               )
      CALL ncio_define_dimension_vector(file_restart, landpatch, 'ens',      DEF_DA_ENS_NUM+1)

      CALL ncio_write_vector(file_restart, 'z_sno',       'snow',     -maxsnl,        'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, z_sno_ens,        DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'dz_sno',      'snow',     -maxsnl,        'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, dz_sno_ens,       DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 't_soisno',    'soilsnow', nl_soil-maxsnl, 'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, t_soisno_ens,     DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'wliq_soisno', 'soilsnow', nl_soil-maxsnl, 'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, wliq_soisno_ens,  DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'wice_soisno', 'soilsnow', nl_soil-maxsnl, 'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, wice_soisno_ens,  DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'smp',         'soil',     nl_soil,        'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, smp_ens,          DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'hk',          'soil',     nl_soil,        'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, hk_ens,           DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 't_grnd',                                  'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, t_grnd_ens,       DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'tleaf',                                   'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, tleaf_ens,        DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'ldew',                                    'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, ldew_ens,         DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'ldew_rain',                               'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, ldew_rain_ens,    DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'ldew_snow',                               'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, ldew_snow_ens,    DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'fwet_snow',                               'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, fwet_snow_ens,    DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'sag',                                     'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, sag_ens,          DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'scv',                                     'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, scv_ens,          DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'snowdp',                                  'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, snowdp_ens,       DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'fveg',                                    'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, fveg_ens,         DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'fsno',                                    'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, fsno_ens,         DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'sigf',                                    'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, sigf_ens,         DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'green',                                   'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, green_ens,        DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'tlai',                                    'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, tlai_ens,         DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'lai',                                     'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, lai_ens,          DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'tsai',                                    'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, tsai_ens,         DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'sai',                                     'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, sai_ens,          DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'alb',         'band',     2, 'rtyp', 2,   'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, alb_ens,          DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'ssun',        'band',     2, 'rtyp', 2,   'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, ssun_ens,         DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'ssha',        'band',     2, 'rtyp', 2,   'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, ssha_ens,         DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'ssoi',        'band',     2, 'rtyp', 2,   'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, ssoi_ens,         DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'ssno',        'band',     2, 'rtyp', 2,   'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, ssno_ens,         DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'thermk',                                  'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, thermk_ens,       DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'extkb',                                   'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, extkb_ens,        DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'extkd',                                   'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, extkd_ens,        DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'zwt',                                     'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, zwt_ens,          DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'wdsrf',                                   'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, wdsrf_ens,        DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'wa',                                      'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, wa_ens,           DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'wetwat',                                  'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, wetwat_ens,       DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 't_lake',      'lake',     nl_lake,        'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, t_lake_ens,       DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'lake_icefrc', 'lake',     nl_lake,        'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, lake_icefrac_ens, DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'savedtke1',                               'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, savedtke1_ens,    DEF_REST_CompressLevel)

      CALL ncio_write_vector(file_restart, 'trad',                                    'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, trad_ens,         DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'tref',                                    'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, tref_ens,         DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'qref',                                    'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, qref_ens,         DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'rst',                                     'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, rst_ens,          DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'emis',                                    'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, emis_ens,         DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'z0m',                                     'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, z0m_ens,          DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'zol',                                     'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, zol_ens,          DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'rib',                                     'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, rib_ens,          DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'ustar',                                   'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, ustar_ens,        DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'qstar',                                   'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, qstar_ens,        DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'tstar',                                   'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, tstar_ens,        DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'fm',                                      'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, fm_ens,           DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'fh',                                      'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, fh_ens,           DEF_REST_CompressLevel)
      CALL ncio_write_vector(file_restart, 'fq',                                      'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, fq_ens,           DEF_REST_CompressLevel)

#ifdef GridRiverLakeFlow
      file_gridriver = trim(dir_restart)// '/'//trim(cdate)//'/' // trim(site) //'_restart_gridriver_'//trim(cdate)//'_lc'//trim(cyear)//'.nc'

      CALL vector_gather_and_write(wdsrf_ucat_a, numucat, totalnumucat, ucat_data_address, file_gridriver, 'wdsrf_ucat_a', 'ucatch', longname = 'ensemble mean river or lake water depth', units = 'm'  )
      CALL vector_gather_and_write(veloc_riv_a, numucat, totalnumucat, ucat_data_address, file_gridriver, 'veloc_riv_a', 'ucatch', longname = 'ensemble mean river velocity', units = 'm/s'             )
      
      IF (DEF_Reservoir_Method > 0 .and. totalnumresv > 0) THEN
         CALL vector_gather_and_write(volresv_a, numresv, totalnumresv, resv_data_address, file_gridriver, 'volresv_a', 'reservoir', longname = 'ensemble mean reservoir water volume', units = 'm^3')
      ENDIF
#endif

   END SUBROUTINE WRITE_DA_TimeVars

!-----------------------------------------------------------------------------

   SUBROUTINE READ_DA_TimeVars(idate, lc_year, site, dir_restart)

!-----------------------------------------------------------------------------
      integer,          intent(in) :: idate (3)
      integer,          intent(in) :: lc_year
      character(len=*), intent(in) :: site
      character(len=*), intent(in) :: dir_restart

      character(len=256)           :: file_restart
      character(len=14)            :: cdate
      character(len=14)            :: cyear
      character(len=256)           :: file_gridriver  
      logical                      :: fexists         
!-----------------------------------------------------------------------------
      write (cyear, '(i4.4)') lc_year
      write (cdate, '(i4.4,"-",i3.3,"-",i5.5)') idate(1), idate(2), idate(3)

      file_restart = trim(dir_restart)//'/'//trim(cdate)//'/'//trim(site)//'_restart_DA_'//trim(cdate)//'_lc'//trim(cyear)//'.nc'

      CALL ncio_read_vector(file_restart, 'z_sno',       -maxsnl,          DEF_DA_ENS_NUM+1, landpatch, z_sno_ens       )
      CALL ncio_read_vector(file_restart, 'dz_sno',      -maxsnl,          DEF_DA_ENS_NUM+1, landpatch, dz_sno_ens      )
      CALL ncio_read_vector(file_restart, 't_soisno',    nl_soil - maxsnl, DEF_DA_ENS_NUM+1, landpatch, t_soisno_ens    )
      CALL ncio_read_vector(file_restart, 'wliq_soisno', nl_soil - maxsnl, DEF_DA_ENS_NUM+1, landpatch, wliq_soisno_ens )
      CALL ncio_read_vector(file_restart, 'wice_soisno', nl_soil - maxsnl, DEF_DA_ENS_NUM+1, landpatch, wice_soisno_ens )
      CALL ncio_read_vector(file_restart, 'smp',         nl_soil,          DEF_DA_ENS_NUM+1, landpatch, smp_ens         )
      CALL ncio_read_vector(file_restart, 'hk',          nl_soil,          DEF_DA_ENS_NUM+1, landpatch, hk_ens          )
      CALL ncio_read_vector(file_restart, 't_grnd',                        DEF_DA_ENS_NUM+1, landpatch, t_grnd_ens      )
      CALL ncio_read_vector(file_restart, 'tleaf',                         DEF_DA_ENS_NUM+1, landpatch, tleaf_ens       )
      CALL ncio_read_vector(file_restart, 'ldew',                          DEF_DA_ENS_NUM+1, landpatch, ldew_ens        )
      CALL ncio_read_vector(file_restart, 'ldew_rain',                     DEF_DA_ENS_NUM+1, landpatch, ldew_rain_ens   )
      CALL ncio_read_vector(file_restart, 'ldew_snow',                     DEF_DA_ENS_NUM+1, landpatch, ldew_snow_ens   )
      CALL ncio_read_vector(file_restart, 'fwet_snow',                     DEF_DA_ENS_NUM+1, landpatch, fwet_snow_ens   )
      CALL ncio_read_vector(file_restart, 'sag',                           DEF_DA_ENS_NUM+1, landpatch, sag_ens         )
      CALL ncio_read_vector(file_restart, 'scv',                           DEF_DA_ENS_NUM+1, landpatch, scv_ens         )
      CALL ncio_read_vector(file_restart, 'snowdp',                        DEF_DA_ENS_NUM+1, landpatch, snowdp_ens      )
      CALL ncio_read_vector(file_restart, 'fveg',                          DEF_DA_ENS_NUM+1, landpatch, fveg_ens        )
      CALL ncio_read_vector(file_restart, 'fsno',                          DEF_DA_ENS_NUM+1, landpatch, fsno_ens        )
      CALL ncio_read_vector(file_restart, 'sigf',                          DEF_DA_ENS_NUM+1, landpatch, sigf_ens        )
      CALL ncio_read_vector(file_restart, 'green',                         DEF_DA_ENS_NUM+1, landpatch, green_ens       )
      CALL ncio_read_vector(file_restart, 'tlai',                          DEF_DA_ENS_NUM+1, landpatch, tlai_ens        )
      CALL ncio_read_vector(file_restart, 'lai',                           DEF_DA_ENS_NUM+1, landpatch, lai_ens         )
      CALL ncio_read_vector(file_restart, 'tsai',                          DEF_DA_ENS_NUM+1, landpatch, tsai_ens        )
      CALL ncio_read_vector(file_restart, 'sai',                           DEF_DA_ENS_NUM+1, landpatch, sai_ens         )
      CALL ncio_read_vector(file_restart, 'alb',         2, 2,             DEF_DA_ENS_NUM+1, landpatch, alb_ens         )
      CALL ncio_read_vector(file_restart, 'ssun',        2, 2,             DEF_DA_ENS_NUM+1, landpatch, ssun_ens        )
      CALL ncio_read_vector(file_restart, 'ssha',        2, 2,             DEF_DA_ENS_NUM+1, landpatch, ssha_ens        )
      CALL ncio_read_vector(file_restart, 'ssoi',        2, 2,             DEF_DA_ENS_NUM+1, landpatch, ssoi_ens        )
      CALL ncio_read_vector(file_restart, 'ssno',        2, 2,             DEF_DA_ENS_NUM+1, landpatch, ssno_ens        )
      CALL ncio_read_vector(file_restart, 'thermk',                        DEF_DA_ENS_NUM+1, landpatch, thermk_ens      )
      CALL ncio_read_vector(file_restart, 'extkb',                         DEF_DA_ENS_NUM+1, landpatch, extkb_ens       )
      CALL ncio_read_vector(file_restart, 'extkd',                         DEF_DA_ENS_NUM+1, landpatch, extkd_ens       )
      CALL ncio_read_vector(file_restart, 'zwt',                           DEF_DA_ENS_NUM+1, landpatch, zwt_ens         )
      CALL ncio_read_vector(file_restart, 'wdsrf',                         DEF_DA_ENS_NUM+1, landpatch, wdsrf_ens       )
      CALL ncio_read_vector(file_restart, 'wa',                            DEF_DA_ENS_NUM+1, landpatch, wa_ens          )
      CALL ncio_read_vector(file_restart, 'wetwat',                        DEF_DA_ENS_NUM+1, landpatch, wetwat_ens      )
      CALL ncio_read_vector(file_restart, 't_lake',      nl_lake,          DEF_DA_ENS_NUM+1, landpatch, t_lake_ens      )
      CALL ncio_read_vector(file_restart, 'lake_icefrc', nl_lake,          DEF_DA_ENS_NUM+1, landpatch, lake_icefrac_ens)
      CALL ncio_read_vector(file_restart, 'savedtke1',                     DEF_DA_ENS_NUM+1, landpatch, savedtke1_ens   )

      CALL ncio_read_vector(file_restart, 'trad',                          DEF_DA_ENS_NUM+1, landpatch, trad_ens        )
      CALL ncio_read_vector(file_restart, 'tref',                          DEF_DA_ENS_NUM+1, landpatch, tref_ens        )
      CALL ncio_read_vector(file_restart, 'qref',                          DEF_DA_ENS_NUM+1, landpatch, qref_ens        )
      CALL ncio_read_vector(file_restart, 'rst',                           DEF_DA_ENS_NUM+1, landpatch, rst_ens         )
      CALL ncio_read_vector(file_restart, 'emis',                          DEF_DA_ENS_NUM+1, landpatch, emis_ens        )
      CALL ncio_read_vector(file_restart, 'z0m',                           DEF_DA_ENS_NUM+1, landpatch, z0m_ens         )
      CALL ncio_read_vector(file_restart, 'zol',                           DEF_DA_ENS_NUM+1, landpatch, zol_ens         )
      CALL ncio_read_vector(file_restart, 'rib',                           DEF_DA_ENS_NUM+1, landpatch, rib_ens         )
      CALL ncio_read_vector(file_restart, 'ustar',                         DEF_DA_ENS_NUM+1, landpatch, ustar_ens       )
      CALL ncio_read_vector(file_restart, 'qstar',                         DEF_DA_ENS_NUM+1, landpatch, qstar_ens       )
      CALL ncio_read_vector(file_restart, 'tstar',                         DEF_DA_ENS_NUM+1, landpatch, tstar_ens       )
      CALL ncio_read_vector(file_restart, 'fm',                            DEF_DA_ENS_NUM+1, landpatch, fm_ens          )
      CALL ncio_read_vector(file_restart, 'fh',                            DEF_DA_ENS_NUM+1, landpatch, fh_ens          )
      CALL ncio_read_vector(file_restart, 'fq',                            DEF_DA_ENS_NUM+1, landpatch, fq_ens          )

#ifdef GridRiverLakeFlow
      file_gridriver = trim(dir_restart)//'/'//trim(cdate)//'/'//trim(site)//'_restart_gridriver_'//trim(cdate)//'_lc'//trim(cyear)//'.nc'

      IF (p_is_master) THEN
         fexists = ncio_var_exist(file_gridriver, 'wdsrf_ucat_a', readflag = .false.) .and.  &
                   ncio_var_exist(file_gridriver, 'veloc_riv_a',  readflag = .false.)
      ENDIF
#ifdef USEMPI
      CALL mpi_bcast(fexists, 1, MPI_LOGICAL, p_address_master, p_comm_glb, p_err)
#endif

      IF (fexists) THEN
         CALL vector_read_and_scatter(file_gridriver, wdsrf_ucat_a, numucat, 'wdsrf_ucat_a', ucat_data_address)
         CALL vector_read_and_scatter(file_gridriver, veloc_riv_a,  numucat, 'veloc_riv_a',  ucat_data_address)
      ENDIF

      IF (DEF_Reservoir_Method > 0 .and. totalnumresv > 0) THEN
         IF (p_is_master) fexists = ncio_var_exist(file_gridriver, 'volresv_a', readflag = .false.)
#ifdef USEMPI
         CALL mpi_bcast(fexists, 1, MPI_LOGICAL, p_address_master, p_comm_glb, p_err)
#endif
         IF (fexists) THEN
            CALL vector_read_and_scatter(file_gridriver, volresv_a, numresv, 'volresv_a', resv_data_address   )
         ENDIF
      ENDIF
#endif

#ifdef RangeCheck
      CALL check_DA_TimeVars
#endif

   END SUBROUTINE READ_DA_TimeVars

!-----------------------------------------------------------------------

   SUBROUTINE check_DA_TimeVars()

!-----------------------------------------------------------------------
      CALL check_vector_data('z_sno       [m]    ', z_sno_ens       )
      CALL check_vector_data('dz_sno      [m]    ', dz_sno_ens      )
      CALL check_vector_data('t_soisno    [K]    ', t_soisno_ens    )
      CALL check_vector_data('wliq_soisno [kg/m2]', wliq_soisno_ens )
      CALL check_vector_data('wice_soisno [kg/m2]', wice_soisno_ens )
      CALL check_vector_data('smp         [mm]   ', smp_ens         )
      CALL check_vector_data('hk          [mm/s] ', hk_ens          )
      CALL check_vector_data('t_grnd      [K]    ', t_grnd_ens      )
      CALL check_vector_data('tleaf       [K]    ', tleaf_ens       )
      CALL check_vector_data('ldew        [mm]   ', ldew_ens        )
      CALL check_vector_data('ldew_rain   [mm]   ', ldew_rain_ens   )
      CALL check_vector_data('ldew_snow   [mm]   ', ldew_snow_ens   )
      CALL check_vector_data('fwet_snow   [-]    ', fwet_snow_ens   )
      CALL check_vector_data('sag         [-]    ', sag_ens         )
      CALL check_vector_data('scv         [mm]   ', scv_ens         )
      CALL check_vector_data('snowdp      [m]    ', snowdp_ens      )
      CALL check_vector_data('fveg        [-]    ', fveg_ens        )
      CALL check_vector_data('fsno        [-]    ', fsno_ens        )
      CALL check_vector_data('sigf        [-]    ', sigf_ens        )
      CALL check_vector_data('green       [-]    ', green_ens       )
      CALL check_vector_data('tlai        [-]    ', tlai_ens        )
      CALL check_vector_data('lai         [-]    ', lai_ens         )
      CALL check_vector_data('tsai        [-]    ', tsai_ens        )
      CALL check_vector_data('sai         [-]    ', sai_ens         )
      CALL check_vector_data('alb         [-]    ', alb_ens         )
      CALL check_vector_data('ssun        [-]    ', ssun_ens        )
      CALL check_vector_data('ssha        [-]    ', ssha_ens        )
      CALL check_vector_data('ssoi        [-]    ', ssoi_ens        )
      CALL check_vector_data('ssno        [-]    ', ssno_ens        )
      CALL check_vector_data('thermk      [-]    ', thermk_ens      )
      CALL check_vector_data('extkb       [-]    ', extkb_ens       )
      CALL check_vector_data('extkd       [-]    ', extkd_ens       )
      CALL check_vector_data('zwt         [m]    ', zwt_ens         )
      CALL check_vector_data('wdsrf       [mm]   ', wdsrf_ens       )
      CALL check_vector_data('wa          [mm]   ', wa_ens          )
      CALL check_vector_data('wetwat      [mm]   ', wetwat_ens      )
      CALL check_vector_data('t_lake      [K]    ', t_lake_ens      )
      CALL check_vector_data('lake_icefrc [-]    ', lake_icefrac_ens)
      CALL check_vector_data('savedtke1   [W/m K]', savedtke1_ens   )

   END SUBROUTINE check_DA_TimeVars

!-----------------------------------------------------------------------
END MODULE MOD_DA_Vars_TimeVars
#endif