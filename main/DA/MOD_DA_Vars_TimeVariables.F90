#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Vars_TimeVariables
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Defines and manages the time-varying model variables required by the data
!    assimilation workflow.
!
!    This module stores the DA ensemble copies of prognostic state variables,
!    selected diagnostic variables, and perturbed forcing variables on local
!    worker patches. The ensemble dimension is always allocated as
!    0:DEF_DA_ENS_NUM:
!
!       member 0                 : open-loop
!       members 1:DEF_DA_ENS_NUM : perturbed ensemble members
!
!    Member 0 is kept in memory and in DA restart files so the model can restore
!    the open-loop trajectory after running ensemble forecasts. Assimilation
!    algorithms should use only members 1:DEF_DA_ENS_NUM.
!
! HISTORY:
!    Lu Li, 12/2024, 05/2025: Initial DA ensemble storage
!    Lu Li, 05/2026: Store open-loop member as ensemble index 0
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_TimeManager
   USE MOD_Namelist
   USE MOD_Precision
   USE MOD_Vars_Global
   USE MOD_SPMD_Task
   USE MOD_LandPatch, only: numpatch, landpatch
   USE MOD_NetCDFVector
#ifdef GridRiverLakeFlow
   USE MOD_DataType,              only: pointer_int32_1d
   USE MOD_NetCDFSerial,          only: ncio_define_dimension, ncio_var_exist, ncio_read_serial, ncio_write_serial, ncio_put_attr
   USE MOD_Vector_ReadWrite,      only: vector_gather_to_master, vector_gather_and_write, vector_read_and_scatter
   USE MOD_Grid_RiverLakeNetwork, only: numucat, totalnumucat, ucat_data_address
   USE MOD_Grid_Reservoir,        only: numresv, totalnumresv, resv_data_address
#endif
#ifdef RangeCheck
      USE MOD_RangeCheck
#endif
   IMPLICIT NONE
   SAVE

   PUBLIC :: allocate_DATimeVariables
   PUBLIC :: deallocate_DATimeVariables
   PUBLIC :: READ_DATimeVariables
   PUBLIC :: WRITE_DATimeVariables
#ifdef RangeCheck
   PUBLIC :: check_DATimeVariables
#endif

!=============================================================================
! DA ensemble state variables
!=============================================================================
   real(r8), allocatable :: z_sno_ens       (:,:,:) ! layer depth [m]
   real(r8), allocatable :: dz_sno_ens      (:,:,:) ! snow layer thickness [m]
   real(r8), allocatable :: t_soisno_ens    (:,:,:) ! snow & soil temperature [K]
   real(r8), allocatable :: wliq_soisno_ens (:,:,:) ! liquid water in layers [kg/m2]
   real(r8), allocatable :: wice_soisno_ens (:,:,:) ! ice in layers [kg/m2]
   real(r8), allocatable :: smp_ens         (:,:,:) ! soil matrix potential [mm]
   real(r8), allocatable :: hk_ens          (:,:,:) ! hydraulic conductivity [mm h2o/s]
   real(r8), allocatable :: t_grnd_ens        (:,:) ! ground surface temperature [K]
   real(r8), allocatable :: tleaf_ens         (:,:) ! leaf temperature [K]
   real(r8), allocatable :: ldew_ens          (:,:) ! depth of water on foliage [mm]
   real(r8), allocatable :: ldew_rain_ens     (:,:) ! depth of rain on foliage [mm]
   real(r8), allocatable :: ldew_snow_ens     (:,:) ! depth of snow on foliage [mm]
   real(r8), allocatable :: fwet_snow_ens     (:,:) ! vegetation snow fractional cover [-]
   real(r8), allocatable :: sag_ens           (:,:) ! non dimensional snow age [-]
   real(r8), allocatable :: scv_ens           (:,:) ! snow water equivalent [kg/m2]
   real(r8), allocatable :: snowdp_ens        (:,:) ! snow depth [m]
   real(r8), allocatable :: fveg_ens          (:,:) ! fraction vegetation cover [-]
   real(r8), allocatable :: fsno_ens          (:,:) ! fraction snow cover [-]
   real(r8), allocatable :: sigf_ens          (:,:) ! fraction of veg cover, excluding snow-covered veg [-]
   real(r8), allocatable :: green_ens         (:,:) ! leaf greenness
   real(r8), allocatable :: tlai_ens          (:,:) ! leaf area index
   real(r8), allocatable :: lai_ens           (:,:) ! leaf area index
   real(r8), allocatable :: tsai_ens          (:,:) ! stem area index
   real(r8), allocatable :: sai_ens           (:,:) ! stem area index
   real(r8), allocatable :: alb_ens       (:,:,:,:) ! averaged albedo [-]
   real(r8), allocatable :: ssun_ens      (:,:,:,:) ! sunlit canopy absorption for solar radiation (0-1)
   real(r8), allocatable :: ssha_ens      (:,:,:,:) ! shaded canopy absorption for solar radiation (0-1)
   real(r8), allocatable :: ssoi_ens      (:,:,:,:) ! soil absorption for solar radiation (0-1)
   real(r8), allocatable :: ssno_ens      (:,:,:,:) ! snow absorption for solar radiation (0-1)
   real(r8), allocatable :: thermk_ens        (:,:) ! canopy gap fraction for tir radiation
   real(r8), allocatable :: extkb_ens         (:,:) ! (k, g(mu)/mu) direct solar extinction coefficient
   real(r8), allocatable :: extkd_ens         (:,:) ! diffuse and scattered diffuse PAR extinction coefficient
   real(r8), allocatable :: zwt_ens           (:,:) ! the depth to water table [m]
   real(r8), allocatable :: wdsrf_ens         (:,:) ! depth of surface water [mm]
   real(r8), allocatable :: wa_ens            (:,:) ! water storage in aquifer [mm]
   real(r8), allocatable :: wetwat_ens        (:,:) ! water storage in wetland [mm]
   real(r8), allocatable :: t_lake_ens      (:,:,:) ! lake layer temperature [K]
   real(r8), allocatable :: lake_icefrac_ens(:,:,:) ! lake mass fraction of lake layer that is frozen
   real(r8), allocatable :: savedtke1_ens     (:,:) ! top level eddy conductivity (W/m K)
#ifdef GridRiverLakeFlow
   real(r8), allocatable :: wdsrf_ucat_mean   (:) ! ensemble-mean river or lake water depth [m]
   real(r8), allocatable :: veloc_riv_mean    (:) ! ensemble-mean river velocity [m/s]
   real(r8), allocatable :: volresv_mean      (:) ! ensemble-mean reservoir water volume [m^3]
   real(r8), allocatable :: wdsrf_ucat_ol     (:) ! open-loop river or lake water depth [m]
   real(r8), allocatable :: veloc_riv_ol      (:) ! open-loop river velocity [m/s]
   real(r8), allocatable :: volresv_ol        (:) ! open-loop reservoir water volume [m^3]

   ! Window sums of the routing accumulators shared by the OL and mean calls.
   real(r8), allocatable :: a_discharge_ol       (:) ! open-loop discharge sum           [m^3/s*s]
   real(r8), allocatable :: acctime_ucat_ol      (:) ! open-loop averaging time          [s]
   real(r8), allocatable :: a_wdsrf_ucat_ol      (:) ! open-loop water depth sum         [m*s]
   real(r8), allocatable :: a_veloc_riv_ol       (:) ! open-loop water velocity sum      [m/s*s]
   real(r8), allocatable :: a_floodarea_ol       (:) ! open-loop flooded area sum        [m^2*s]
   real(r8), allocatable :: a_discharge_mean     (:) ! ensemble-mean discharge sum       [m^3/s*s]
   real(r8), allocatable :: acctime_ucat_mean    (:) ! ensemble-mean averaging time      [s]
   real(r8), allocatable :: a_wdsrf_ucat_mean    (:) ! ensemble-mean water depth sum     [m*s]
   real(r8), allocatable :: a_veloc_riv_mean     (:) ! ensemble-mean water velocity sum  [m/s*s]
   real(r8), allocatable :: a_floodarea_mean     (:) ! ensemble-mean flooded area sum    [m^2*s]

   ! Patch-level counterparts of the sums above, remapped at history output.
   real(r8), allocatable :: a_discharge_mean_pch   (:) ! ensemble-mean discharge on patches     [m^3/s]
   real(r8), allocatable :: a_dis_rmth_mean_pch    (:) ! ensemble-mean river mouth discharge on patches [m^3/s]
   real(r8), allocatable :: a_wdsrf_ucat_mean_pch  (:) ! ensemble-mean water depth on patches   [m]
   real(r8), allocatable :: a_veloc_riv_mean_pch   (:) ! ensemble-mean water velocity on patches[m/s]
   real(r8), allocatable :: a_floodfrc_mean_pch    (:) ! ensemble-mean flooded fraction on patches[-]

   ! Reservoir window sums, handled the same; unrelated to the discharge output.
   real(r8), allocatable :: acctime_resv_ol      (:) ! open-loop averaging time           [s]
   real(r8), allocatable :: a_volresv_ol         (:) ! open-loop reservoir volume sum     [m^3*s]
   real(r8), allocatable :: a_qresv_in_ol        (:) ! open-loop reservoir inflow sum     [m^3/s*s]
   real(r8), allocatable :: a_qresv_out_ol       (:) ! open-loop reservoir outflow sum    [m^3/s*s]
   real(r8), allocatable :: acctime_resv_mean    (:) ! ensemble-mean averaging time       [s]
   real(r8), allocatable :: a_volresv_mean       (:) ! ensemble-mean reservoir volume sum [m^3*s]
   real(r8), allocatable :: a_qresv_in_mean      (:) ! ensemble-mean reservoir inflow sum [m^3/s*s]
   real(r8), allocatable :: a_qresv_out_mean     (:) ! ensemble-mean reservoir outflow sum[m^3/s*s]
#endif

!=============================================================================
! diagnostic variables
!=============================================================================
   integer :: nsource = 0                           ! number of observation sources configured in namelist
   real(r8), allocatable :: hx_f_ens        (:,:,:) ! prior H(x_f) [source, ensemble, patch]
   real(r8), allocatable :: hx_a_ens        (:,:,:) ! posterior H(x_a) [source, ensemble, patch]
   real(r8), allocatable :: hx_ol             (:,:) ! open-loop H(x) [source, patch]
   real(r8), allocatable :: hx_f              (:,:) ! ensemble-mean prior H(x_f) [source, patch]
   real(r8), allocatable :: hx_a              (:,:) ! ensemble-mean posterior H(x_a) [source, patch]
   real(r8), allocatable :: wliq_soisno_ol    (:,:) ! open-loop liquid water in layers [kg/m2]
   real(r8), allocatable :: wliq_soisno_f     (:,:) ! ensemble-mean prior liquid water in layers [kg/m2]
   real(r8), allocatable :: wliq_soisno_a     (:,:) ! ensemble-mean posterior liquid water in layers [kg/m2]
   real(r8), allocatable :: t_soisno_ol       (:,:) ! open-loop temperature in layers [K]
   real(r8), allocatable :: t_soisno_f        (:,:) ! ensemble-mean prior temperature in layers [K]
   real(r8), allocatable :: t_soisno_a        (:,:) ! ensemble-mean posterior temperature in layers [K]

!=============================================================================
! land-atmosphere coupling variables
!=============================================================================
   real(r8), allocatable :: trad_ens          (:,:) ! radiative temperature of surface [K]
   real(r8), allocatable :: tref_ens          (:,:) ! 2 m height air temperature [kelvin]
   real(r8), allocatable :: qref_ens          (:,:) ! 2 m height air specific humidity
   real(r8), allocatable :: rhref_ens         (:,:) ! 2 m height air relative humidity
   real(r8), allocatable :: ustar_ens         (:,:) ! u* in similarity theory [m/s]
   real(r8), allocatable :: qstar_ens         (:,:) ! q* in similarity theory [kg/kg]
   real(r8), allocatable :: tstar_ens         (:,:) ! t* in similarity theory [K]
   real(r8), allocatable :: fm_ens            (:,:) ! integral of profile FUNCTION for momentum
   real(r8), allocatable :: fh_ens            (:,:) ! integral of profile FUNCTION for heat
   real(r8), allocatable :: fq_ens            (:,:) ! integral of profile FUNCTION for moisture

!=============================================================================
! DA ensemble forcing variables
!=============================================================================
   real(r8), allocatable :: forc_t_ens        (:,:) ! temperature [K]
   real(r8), allocatable :: forc_frl_ens      (:,:) ! atmospheric infrared (longwave) radiation [W/m2]
   real(r8), allocatable :: forc_prc_ens      (:,:) ! convective precipitation [mm/s]
   real(r8), allocatable :: forc_prl_ens      (:,:) ! large scale precipitation [mm/s]
   real(r8), allocatable :: forc_sols_ens     (:,:) ! atm vis direct beam solar rad onto srf [W/m2]
   real(r8), allocatable :: forc_soll_ens     (:,:) ! atm nir direct beam solar rad onto srf [W/m2]
   real(r8), allocatable :: forc_solsd_ens    (:,:) ! atm vis diffuse solar rad onto srf [W/m2]
   real(r8), allocatable :: forc_solld_ens    (:,:) ! atm nir diffuse solar rad onto srf [W/m2]

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE allocate_DATimeVariables()

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      integer :: i

!-----------------------------------------------------------------------------
      nsource = 0
      DO i = 1, size(DEF_DA_OBS_TARGET)
         IF (trim(DEF_DA_OBS_TARGET(i)) == 'SM' .or. &
             trim(DEF_DA_OBS_TARGET(i)) == 'ST') THEN
            nsource = nsource + 1
         ENDIF
      ENDDO

      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            allocate (z_sno_ens      (maxsnl+1:0,      0:DEF_DA_ENS_NUM,numpatch)); z_sno_ens       (:,:,:) = spval
            allocate (dz_sno_ens     (maxsnl+1:0,      0:DEF_DA_ENS_NUM,numpatch)); dz_sno_ens      (:,:,:) = spval
            allocate (t_soisno_ens   (maxsnl+1:nl_soil,0:DEF_DA_ENS_NUM,numpatch)); t_soisno_ens    (:,:,:) = spval
            allocate (wliq_soisno_ens(maxsnl+1:nl_soil,0:DEF_DA_ENS_NUM,numpatch)); wliq_soisno_ens (:,:,:) = spval
            allocate (wice_soisno_ens(maxsnl+1:nl_soil,0:DEF_DA_ENS_NUM,numpatch)); wice_soisno_ens (:,:,:) = spval
            allocate (smp_ens               (1:nl_soil,0:DEF_DA_ENS_NUM,numpatch)); smp_ens         (:,:,:) = spval
            allocate (hk_ens                (1:nl_soil,0:DEF_DA_ENS_NUM,numpatch)); hk_ens          (:,:,:) = spval
            allocate (t_grnd_ens                      (0:DEF_DA_ENS_NUM,numpatch)); t_grnd_ens        (:,:) = spval
            allocate (tleaf_ens                       (0:DEF_DA_ENS_NUM,numpatch)); tleaf_ens         (:,:) = spval
            allocate (ldew_ens                        (0:DEF_DA_ENS_NUM,numpatch)); ldew_ens          (:,:) = spval
            allocate (ldew_rain_ens                   (0:DEF_DA_ENS_NUM,numpatch)); ldew_rain_ens     (:,:) = spval
            allocate (ldew_snow_ens                   (0:DEF_DA_ENS_NUM,numpatch)); ldew_snow_ens     (:,:) = spval
            allocate (fwet_snow_ens                   (0:DEF_DA_ENS_NUM,numpatch)); fwet_snow_ens     (:,:) = spval
            allocate (sag_ens                         (0:DEF_DA_ENS_NUM,numpatch)); sag_ens           (:,:) = spval
            allocate (scv_ens                         (0:DEF_DA_ENS_NUM,numpatch)); scv_ens           (:,:) = spval
            allocate (snowdp_ens                      (0:DEF_DA_ENS_NUM,numpatch)); snowdp_ens        (:,:) = spval
            allocate (fveg_ens                        (0:DEF_DA_ENS_NUM,numpatch)); fveg_ens          (:,:) = spval
            allocate (fsno_ens                        (0:DEF_DA_ENS_NUM,numpatch)); fsno_ens          (:,:) = spval
            allocate (sigf_ens                        (0:DEF_DA_ENS_NUM,numpatch)); sigf_ens          (:,:) = spval
            allocate (green_ens                       (0:DEF_DA_ENS_NUM,numpatch)); green_ens         (:,:) = spval
            allocate (tlai_ens                        (0:DEF_DA_ENS_NUM,numpatch)); tlai_ens          (:,:) = spval
            allocate (lai_ens                         (0:DEF_DA_ENS_NUM,numpatch)); lai_ens           (:,:) = spval
            allocate (tsai_ens                        (0:DEF_DA_ENS_NUM,numpatch)); tsai_ens          (:,:) = spval
            allocate (sai_ens                         (0:DEF_DA_ENS_NUM,numpatch)); sai_ens           (:,:) = spval
            allocate (alb_ens                     (2,2,0:DEF_DA_ENS_NUM,numpatch)); alb_ens       (:,:,:,:) = spval
            allocate (ssun_ens                    (2,2,0:DEF_DA_ENS_NUM,numpatch)); ssun_ens      (:,:,:,:) = spval
            allocate (ssha_ens                    (2,2,0:DEF_DA_ENS_NUM,numpatch)); ssha_ens      (:,:,:,:) = spval
            allocate (ssoi_ens                    (2,2,0:DEF_DA_ENS_NUM,numpatch)); ssoi_ens      (:,:,:,:) = spval
            allocate (ssno_ens                    (2,2,0:DEF_DA_ENS_NUM,numpatch)); ssno_ens      (:,:,:,:) = spval
            allocate (thermk_ens                      (0:DEF_DA_ENS_NUM,numpatch)); thermk_ens        (:,:) = spval
            allocate (extkb_ens                       (0:DEF_DA_ENS_NUM,numpatch)); extkb_ens         (:,:) = spval
            allocate (extkd_ens                       (0:DEF_DA_ENS_NUM,numpatch)); extkd_ens         (:,:) = spval
            allocate (zwt_ens                         (0:DEF_DA_ENS_NUM,numpatch)); zwt_ens           (:,:) = spval
            allocate (wdsrf_ens                       (0:DEF_DA_ENS_NUM,numpatch)); wdsrf_ens         (:,:) = spval
            allocate (wa_ens                          (0:DEF_DA_ENS_NUM,numpatch)); wa_ens            (:,:) = spval
            allocate (wetwat_ens                      (0:DEF_DA_ENS_NUM,numpatch)); wetwat_ens        (:,:) = spval
            allocate (t_lake_ens              (nl_lake,0:DEF_DA_ENS_NUM,numpatch)); t_lake_ens      (:,:,:) = spval
            allocate (lake_icefrac_ens        (nl_lake,0:DEF_DA_ENS_NUM,numpatch)); lake_icefrac_ens(:,:,:) = spval
            allocate (savedtke1_ens                   (0:DEF_DA_ENS_NUM,numpatch)); savedtke1_ens     (:,:) = spval

            IF (nsource > 0) THEN
               allocate (hx_ol                  (nsource,               numpatch)); hx_ol             (:,:) = spval
               allocate (hx_f_ens               (nsource,DEF_DA_ENS_NUM,numpatch)); hx_f_ens     (:,:,:) = spval
               allocate (hx_a_ens               (nsource,DEF_DA_ENS_NUM,numpatch)); hx_a_ens        (:,:,:) = spval
               allocate (hx_f                   (nsource,               numpatch)); hx_f              (:,:) = spval
               allocate (hx_a                   (nsource,               numpatch)); hx_a              (:,:) = spval
            ENDIF
            allocate (wliq_soisno_ol            (maxsnl+1:nl_soil,      numpatch)); wliq_soisno_ol    (:,:) = spval
            allocate (wliq_soisno_f             (maxsnl+1:nl_soil,      numpatch)); wliq_soisno_f     (:,:) = spval
            allocate (wliq_soisno_a             (maxsnl+1:nl_soil,      numpatch)); wliq_soisno_a     (:,:) = spval
            allocate (t_soisno_ol               (maxsnl+1:nl_soil,      numpatch)); t_soisno_ol       (:,:) = spval
            allocate (t_soisno_f                (maxsnl+1:nl_soil,      numpatch)); t_soisno_f        (:,:) = spval
            allocate (t_soisno_a                (maxsnl+1:nl_soil,      numpatch)); t_soisno_a        (:,:) = spval

            allocate (trad_ens                        (0:DEF_DA_ENS_NUM,numpatch)); trad_ens          (:,:) = spval
            allocate (tref_ens                        (0:DEF_DA_ENS_NUM,numpatch)); tref_ens          (:,:) = spval
            allocate (qref_ens                        (0:DEF_DA_ENS_NUM,numpatch)); qref_ens          (:,:) = spval
            allocate (rhref_ens                       (0:DEF_DA_ENS_NUM,numpatch)); rhref_ens         (:,:) = spval
            allocate (ustar_ens                       (0:DEF_DA_ENS_NUM,numpatch)); ustar_ens         (:,:) = spval
            allocate (qstar_ens                       (0:DEF_DA_ENS_NUM,numpatch)); qstar_ens         (:,:) = spval
            allocate (tstar_ens                       (0:DEF_DA_ENS_NUM,numpatch)); tstar_ens         (:,:) = spval
            allocate (fm_ens                          (0:DEF_DA_ENS_NUM,numpatch)); fm_ens            (:,:) = spval
            allocate (fh_ens                          (0:DEF_DA_ENS_NUM,numpatch)); fh_ens            (:,:) = spval
            allocate (fq_ens                          (0:DEF_DA_ENS_NUM,numpatch)); fq_ens            (:,:) = spval

            allocate (forc_t_ens                      (0:DEF_DA_ENS_NUM,numpatch)); forc_t_ens        (:,:) = spval
            allocate (forc_frl_ens                    (0:DEF_DA_ENS_NUM,numpatch)); forc_frl_ens      (:,:) = spval
            allocate (forc_prc_ens                    (0:DEF_DA_ENS_NUM,numpatch)); forc_prc_ens      (:,:) = spval
            allocate (forc_prl_ens                    (0:DEF_DA_ENS_NUM,numpatch)); forc_prl_ens      (:,:) = spval
            allocate (forc_sols_ens                   (0:DEF_DA_ENS_NUM,numpatch)); forc_sols_ens     (:,:) = spval
            allocate (forc_soll_ens                   (0:DEF_DA_ENS_NUM,numpatch)); forc_soll_ens     (:,:) = spval
            allocate (forc_solsd_ens                  (0:DEF_DA_ENS_NUM,numpatch)); forc_solsd_ens    (:,:) = spval
            allocate (forc_solld_ens                  (0:DEF_DA_ENS_NUM,numpatch)); forc_solld_ens    (:,:) = spval
         ENDIF

#ifdef GridRiverLakeFlow
         IF (numucat > 0) THEN
            allocate (wdsrf_ucat_mean     (numucat)); wdsrf_ucat_mean (:)     = spval
            allocate (veloc_riv_mean      (numucat)); veloc_riv_mean  (:)     = spval
            allocate (wdsrf_ucat_ol       (numucat)); wdsrf_ucat_ol   (:)     = spval
            allocate (veloc_riv_ol        (numucat)); veloc_riv_ol    (:)     = spval

            allocate (a_discharge_ol      (numucat)); a_discharge_ol      (:) = spval
            allocate (acctime_ucat_ol     (numucat)); acctime_ucat_ol     (:) = spval
            allocate (a_wdsrf_ucat_ol     (numucat)); a_wdsrf_ucat_ol     (:) = spval
            allocate (a_veloc_riv_ol      (numucat)); a_veloc_riv_ol      (:) = spval
            allocate (a_floodarea_ol      (numucat)); a_floodarea_ol      (:) = spval
            allocate (a_discharge_mean    (numucat)); a_discharge_mean    (:) = 0.
            allocate (acctime_ucat_mean   (numucat)); acctime_ucat_mean   (:) = 0.
            allocate (a_wdsrf_ucat_mean   (numucat)); a_wdsrf_ucat_mean   (:) = 0.
            allocate (a_veloc_riv_mean    (numucat)); a_veloc_riv_mean    (:) = 0.
            allocate (a_floodarea_mean    (numucat)); a_floodarea_mean    (:) = 0.
         ENDIF
         IF (numpatch > 0) THEN
            allocate (a_discharge_mean_pch  (numpatch)); a_discharge_mean_pch  (:) = spval
            allocate (a_dis_rmth_mean_pch   (numpatch)); a_dis_rmth_mean_pch   (:) = spval
            allocate (a_wdsrf_ucat_mean_pch (numpatch)); a_wdsrf_ucat_mean_pch (:) = spval
            allocate (a_veloc_riv_mean_pch  (numpatch)); a_veloc_riv_mean_pch  (:) = spval
            allocate (a_floodfrc_mean_pch   (numpatch)); a_floodfrc_mean_pch   (:) = spval
         ENDIF
         IF (numresv > 0) THEN
            allocate (volresv_mean        (numresv)); volresv_mean    (:)     = spval
            allocate (volresv_ol          (numresv)); volresv_ol      (:)     = spval

            allocate (acctime_resv_ol     (numresv)); acctime_resv_ol     (:) = spval
            allocate (a_volresv_ol        (numresv)); a_volresv_ol        (:) = spval
            allocate (a_qresv_in_ol       (numresv)); a_qresv_in_ol       (:) = spval
            allocate (a_qresv_out_ol      (numresv)); a_qresv_out_ol      (:) = spval
            allocate (acctime_resv_mean   (numresv)); acctime_resv_mean   (:) = 0.
            allocate (a_volresv_mean      (numresv)); a_volresv_mean      (:) = 0.
            allocate (a_qresv_in_mean     (numresv)); a_qresv_in_mean     (:) = 0.
            allocate (a_qresv_out_mean    (numresv)); a_qresv_out_mean    (:) = 0.
         ENDIF
#endif
      ENDIF

   END SUBROUTINE allocate_DATimeVariables

!-----------------------------------------------------------------------------

   SUBROUTINE deallocate_DATimeVariables()

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!-----------------------------------------------------------------------------
      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            deallocate (z_sno_ens                        )
            deallocate (dz_sno_ens                       )
            deallocate (t_soisno_ens                     )
            deallocate (wliq_soisno_ens                  )
            deallocate (wice_soisno_ens                  )
            deallocate (smp_ens                          )
            deallocate (hk_ens                           )
            deallocate (t_grnd_ens                       )
            deallocate (tleaf_ens                        )
            deallocate (ldew_ens                         )
            deallocate (ldew_rain_ens                    )
            deallocate (ldew_snow_ens                    )
            deallocate (fwet_snow_ens                    )
            deallocate (sag_ens                          )
            deallocate (scv_ens                          )
            deallocate (snowdp_ens                       )
            deallocate (fveg_ens                         )
            deallocate (fsno_ens                         )
            deallocate (sigf_ens                         )
            deallocate (green_ens                        )
            deallocate (tlai_ens                         )
            deallocate (lai_ens                          )
            deallocate (tsai_ens                         )
            deallocate (sai_ens                          )
            deallocate (alb_ens                          )
            deallocate (ssun_ens                         )
            deallocate (ssha_ens                         )
            deallocate (ssoi_ens                         )
            deallocate (ssno_ens                         )
            deallocate (thermk_ens                       )
            deallocate (extkb_ens                        )
            deallocate (extkd_ens                        )
            deallocate (zwt_ens                          )
            deallocate (wdsrf_ens                        )
            deallocate (wa_ens                           )
            deallocate (wetwat_ens                       )
            deallocate (t_lake_ens                       )
            deallocate (lake_icefrac_ens                 )
            deallocate (savedtke1_ens                    )

            IF (allocated(hx_ol))    deallocate (hx_ol   )
            IF (allocated(hx_f_ens)) deallocate (hx_f_ens)
            IF (allocated(hx_a_ens)) deallocate (hx_a_ens)
            IF (allocated(hx_f))     deallocate (hx_f    )
            IF (allocated(hx_a))     deallocate (hx_a    )
            deallocate (wliq_soisno_ol                   )
            deallocate (wliq_soisno_f                    )
            deallocate (wliq_soisno_a                    )
            deallocate (t_soisno_ol                      )
            deallocate (t_soisno_f                       )
            deallocate (t_soisno_a                       )

            deallocate (trad_ens                         )
            deallocate (tref_ens                         )
            deallocate (qref_ens                         )
            deallocate (rhref_ens                        )
            deallocate (ustar_ens                        )
            deallocate (qstar_ens                        )
            deallocate (tstar_ens                        )
            deallocate (fm_ens                           )
            deallocate (fh_ens                           )
            deallocate (fq_ens                           )

            deallocate (forc_t_ens                       )
            deallocate (forc_frl_ens                     )
            deallocate (forc_prc_ens                     )
            deallocate (forc_prl_ens                     )
            deallocate (forc_sols_ens                    )
            deallocate (forc_soll_ens                    )
            deallocate (forc_solsd_ens                   )
            deallocate (forc_solld_ens                   )
         ENDIF

#ifdef GridRiverLakeFlow
         IF (allocated(wdsrf_ucat_mean     )) deallocate (wdsrf_ucat_mean     )
         IF (allocated(veloc_riv_mean      )) deallocate (veloc_riv_mean      )
         IF (allocated(volresv_mean        )) deallocate (volresv_mean        )
         IF (allocated(wdsrf_ucat_ol       )) deallocate (wdsrf_ucat_ol       )
         IF (allocated(veloc_riv_ol        )) deallocate (veloc_riv_ol        )
         IF (allocated(volresv_ol          )) deallocate (volresv_ol          )
         IF (allocated(a_discharge_ol      )) deallocate (a_discharge_ol      )
         IF (allocated(acctime_ucat_ol     )) deallocate (acctime_ucat_ol     )
         IF (allocated(a_wdsrf_ucat_ol     )) deallocate (a_wdsrf_ucat_ol     )
         IF (allocated(a_veloc_riv_ol      )) deallocate (a_veloc_riv_ol      )
         IF (allocated(a_floodarea_ol      )) deallocate (a_floodarea_ol      )
         IF (allocated(a_discharge_mean    )) deallocate (a_discharge_mean    )
         IF (allocated(acctime_ucat_mean   )) deallocate (acctime_ucat_mean   )
         IF (allocated(a_wdsrf_ucat_mean   )) deallocate (a_wdsrf_ucat_mean   )
         IF (allocated(a_veloc_riv_mean    )) deallocate (a_veloc_riv_mean    )
         IF (allocated(a_floodarea_mean    )) deallocate (a_floodarea_mean    )
         IF (allocated(a_discharge_mean_pch  )) deallocate (a_discharge_mean_pch  )
         IF (allocated(a_dis_rmth_mean_pch   )) deallocate (a_dis_rmth_mean_pch   )
         IF (allocated(a_wdsrf_ucat_mean_pch )) deallocate (a_wdsrf_ucat_mean_pch )
         IF (allocated(a_veloc_riv_mean_pch  )) deallocate (a_veloc_riv_mean_pch  )
         IF (allocated(a_floodfrc_mean_pch   )) deallocate (a_floodfrc_mean_pch   )
         IF (allocated(acctime_resv_ol     )) deallocate (acctime_resv_ol     )
         IF (allocated(a_volresv_ol        )) deallocate (a_volresv_ol        )
         IF (allocated(a_qresv_in_ol       )) deallocate (a_qresv_in_ol       )
         IF (allocated(a_qresv_out_ol      )) deallocate (a_qresv_out_ol      )
         IF (allocated(acctime_resv_mean   )) deallocate (acctime_resv_mean   )
         IF (allocated(a_volresv_mean      )) deallocate (a_volresv_mean      )
         IF (allocated(a_qresv_in_mean     )) deallocate (a_qresv_in_mean     )
         IF (allocated(a_qresv_out_mean    )) deallocate (a_qresv_out_mean    )
#endif
      ENDIF

   END SUBROUTINE deallocate_DATimeVariables

!-----------------------------------------------------------------------------

   SUBROUTINE WRITE_DATimeVariables (idate, lc_year, site, dir_restart)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      integer,          intent(in) :: idate(3)    ! current model date [year, day, seconds]
      integer,          intent(in) :: lc_year     ! land-cover year
      character(len=*), intent(in) :: site        ! case/site name used in restart filename
      character(len=*), intent(in) :: dir_restart ! restart output directory

!------------------------ Local Variables ------------------------------------
      character(len=256) :: file_restart          ! full DA restart filename
      character(len=14)  :: cdate                 ! date string used in restart path
      character(len=256) :: cyear                 ! land-cover year string
      integer :: compress                         ! NetCDF compression level
      character(len=256) :: file_gridriver        ! gridriver restart filename

!-----------------------------------------------------------------------------
      compress = DEF_REST_CompressLevel

      write(cyear,'(i4.4)') lc_year
      write(cdate,'(i4.4,"-",i3.3,"-",i5.5)') idate(1), idate(2), idate(3)

      IF (p_is_master) THEN
         CALL system('mkdir -p ' // trim(dir_restart)//'/'//trim(cdate))
      ENDIF
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      file_restart = trim(dir_restart)// '/'//trim(cdate)//'/' // trim(site) //'_restart_DA_'//trim(cdate)//'_lc'//trim(cyear)//'.nc'

      CALL ncio_create_file_vector      (file_restart, landpatch)

      CALL ncio_define_dimension_vector (file_restart, landpatch, 'patch'                   )
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'snow',     -maxsnl       )
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'snowp1',   -maxsnl+1     )
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'soilsnow', nl_soil-maxsnl)
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'soil',     nl_soil       )
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'lake',     nl_lake       )
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'band',     2             )
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'rtyp',     2             )
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'ens',    DEF_DA_ENS_NUM+1)

      CALL ncio_write_vector (file_restart, 'z_sno   '   , 'snow',     -maxsnl,        'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, z_sno_ens,        compress)
      CALL ncio_write_vector (file_restart, 'dz_sno  '   , 'snow',     -maxsnl,        'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, dz_sno_ens,       compress)
      CALL ncio_write_vector (file_restart, 't_soisno'   , 'soilsnow', nl_soil-maxsnl, 'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, t_soisno_ens,     compress)
      CALL ncio_write_vector (file_restart, 'wliq_soisno', 'soilsnow', nl_soil-maxsnl, 'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, wliq_soisno_ens,  compress)
      CALL ncio_write_vector (file_restart, 'wice_soisno', 'soilsnow', nl_soil-maxsnl, 'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, wice_soisno_ens,  compress)
      CALL ncio_write_vector (file_restart, 'smp',         'soil',     nl_soil,        'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, smp_ens,          compress)
      CALL ncio_write_vector (file_restart, 'hk',          'soil',     nl_soil,        'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, hk_ens,           compress)
      CALL ncio_write_vector (file_restart, 't_grnd',                                  'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, t_grnd_ens,       compress)
      CALL ncio_write_vector (file_restart, 'tleaf',                                   'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, tleaf_ens,        compress)
      CALL ncio_write_vector (file_restart, 'ldew',                                    'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, ldew_ens,         compress)
      CALL ncio_write_vector (file_restart, 'ldew_rain',                               'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, ldew_rain_ens,    compress)
      CALL ncio_write_vector (file_restart, 'ldew_snow',                               'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, ldew_snow_ens,    compress)
      CALL ncio_write_vector (file_restart, 'fwet_snow',                               'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, fwet_snow_ens,    compress)
      CALL ncio_write_vector (file_restart, 'sag',                                     'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, sag_ens,          compress)
      CALL ncio_write_vector (file_restart, 'scv',                                     'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, scv_ens,          compress)
      CALL ncio_write_vector (file_restart, 'snowdp',                                  'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, snowdp_ens,       compress)
      CALL ncio_write_vector (file_restart, 'fveg',                                    'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, fveg_ens,         compress)
      CALL ncio_write_vector (file_restart, 'fsno',                                    'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, fsno_ens,         compress)
      CALL ncio_write_vector (file_restart, 'sigf',                                    'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, sigf_ens,         compress)
      CALL ncio_write_vector (file_restart, 'green',                                   'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, green_ens,        compress)
      CALL ncio_write_vector (file_restart, 'tlai',                                    'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, tlai_ens,         compress)
      CALL ncio_write_vector (file_restart, 'lai',                                     'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, lai_ens,          compress)
      CALL ncio_write_vector (file_restart, 'tsai',                                    'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, tsai_ens,         compress)
      CALL ncio_write_vector (file_restart, 'sai',                                     'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, sai_ens,          compress)
      CALL ncio_write_vector (file_restart, 'alb',         'band',     2, 'rtyp', 2,   'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, alb_ens,          compress)
      CALL ncio_write_vector (file_restart, 'ssun',        'band',     2, 'rtyp', 2,   'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, ssun_ens,         compress)
      CALL ncio_write_vector (file_restart, 'ssha',        'band',     2, 'rtyp', 2,   'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, ssha_ens,         compress)
      CALL ncio_write_vector (file_restart, 'ssoi',        'band',     2, 'rtyp', 2,   'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, ssoi_ens,         compress)
      CALL ncio_write_vector (file_restart, 'ssno',        'band',     2, 'rtyp', 2,   'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, ssno_ens,         compress)
      CALL ncio_write_vector (file_restart, 'thermk',                                  'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, thermk_ens,       compress)
      CALL ncio_write_vector (file_restart, 'extkb',                                   'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, extkb_ens,        compress)
      CALL ncio_write_vector (file_restart, 'extkd',                                   'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, extkd_ens,        compress)
      CALL ncio_write_vector (file_restart, 'zwt',                                     'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, zwt_ens,          compress)
      CALL ncio_write_vector (file_restart, 'wdsrf',                                   'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, wdsrf_ens,        compress)
      CALL ncio_write_vector (file_restart, 'wa',                                      'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, wa_ens,           compress)
      CALL ncio_write_vector (file_restart, 'wetwat',                                  'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, wetwat_ens,       compress)
      CALL ncio_write_vector (file_restart, 't_lake',      'lake',     nl_lake,        'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, t_lake_ens,       compress)
      CALL ncio_write_vector (file_restart, 'lake_icefrc', 'lake',     nl_lake,        'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, lake_icefrac_ens, compress)
      CALL ncio_write_vector (file_restart, 'savedtke1',                               'ens',      DEF_DA_ENS_NUM+1, 'patch', landpatch, savedtke1_ens,    compress)

#ifdef GridRiverLakeFlow
      file_gridriver = trim(dir_restart)// '/'//trim(cdate)//'/' // trim(site) //'_restart_gridriver_'//trim(cdate)//'_lc'//trim(cyear)//'.nc'

      CALL vector_gather_and_write (wdsrf_ucat_mean, numucat, totalnumucat, ucat_data_address, &
         file_gridriver, 'wdsrf_ucat_mean', 'ucatch', &
         longname = 'ensemble mean river or lake water depth', units = 'm')

      CALL vector_gather_and_write (veloc_riv_mean, numucat, totalnumucat, ucat_data_address, &
         file_gridriver, 'veloc_riv_mean', 'ucatch', &
         longname = 'ensemble mean river velocity', units = 'm/s')

      IF (DEF_Reservoir_Method > 0) THEN
         IF (totalnumresv > 0) THEN
            CALL vector_gather_and_write (volresv_mean, numresv, totalnumresv, resv_data_address, &
               file_gridriver, 'volresv_mean', 'reservoir', &
               longname = 'ensemble mean reservoir water volume', units = 'm^3')
         ENDIF
      ENDIF
#endif

   END SUBROUTINE WRITE_DATimeVariables

!-----------------------------------------------------------------------------

   SUBROUTINE READ_DATimeVariables (idate, lc_year, site, dir_restart)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      integer,          intent(in) :: idate(3)
      integer,          intent(in) :: lc_year
      character(len=*), intent(in) :: site
      character(len=*), intent(in) :: dir_restart

!------------------------ Local Variables ------------------------------------
      character(len=256) :: file_restart
      character(len=14)  :: cdate
      character(len=14)  :: cyear
      character(len=256) :: file_gridriver        ! gridriver restart filename
      logical            :: fexists               ! ensemble river states present in restart file

!-----------------------------------------------------------------------------
#ifdef USEMPI
      CALL mpi_barrier(p_comm_glb, p_err)
#endif

      IF (p_is_master) THEN
         write (*, *) 'Loading DA Time Variables ...'
      END IF

      write (cyear, '(i4.4)') lc_year
      write (cdate, '(i4.4,"-",i3.3,"-",i5.5)') idate(1), idate(2), idate(3)

      file_restart = trim(dir_restart)//'/'//trim(cdate)//'/'//trim(site)//'_restart_DA_'//trim(cdate)//'_lc'//trim(cyear)//'.nc'

      CALL ncio_read_vector(file_restart, 'z_sno   '   , -maxsnl,          DEF_DA_ENS_NUM+1, landpatch, z_sno_ens       )
      CALL ncio_read_vector(file_restart, 'dz_sno  '   , -maxsnl,          DEF_DA_ENS_NUM+1, landpatch, dz_sno_ens      )
      CALL ncio_read_vector(file_restart, 't_soisno'   , nl_soil - maxsnl, DEF_DA_ENS_NUM+1, landpatch, t_soisno_ens    )
      CALL ncio_read_vector(file_restart, 'wliq_soisno', nl_soil - maxsnl, DEF_DA_ENS_NUM+1, landpatch, wliq_soisno_ens )
      CALL ncio_read_vector(file_restart, 'wice_soisno', nl_soil - maxsnl, DEF_DA_ENS_NUM+1, landpatch, wice_soisno_ens )
      CALL ncio_read_vector(file_restart, 'smp'        , nl_soil,          DEF_DA_ENS_NUM+1, landpatch, smp_ens         )
      CALL ncio_read_vector(file_restart, 'hk'         , nl_soil,          DEF_DA_ENS_NUM+1, landpatch, hk_ens          )
      CALL ncio_read_vector(file_restart, 't_grnd  '   ,                   DEF_DA_ENS_NUM+1, landpatch, t_grnd_ens      )
      CALL ncio_read_vector(file_restart, 'tleaf   '   ,                   DEF_DA_ENS_NUM+1, landpatch, tleaf_ens       )
      CALL ncio_read_vector(file_restart, 'ldew    '   ,                   DEF_DA_ENS_NUM+1, landpatch, ldew_ens        )
      CALL ncio_read_vector(file_restart, 'ldew_rain'  ,                   DEF_DA_ENS_NUM+1, landpatch, ldew_rain_ens   )
      CALL ncio_read_vector(file_restart, 'ldew_snow'  ,                   DEF_DA_ENS_NUM+1, landpatch, ldew_snow_ens   )
      CALL ncio_read_vector(file_restart, 'fwet_snow'  ,                   DEF_DA_ENS_NUM+1, landpatch, fwet_snow_ens   )
      CALL ncio_read_vector(file_restart, 'sag     '   ,                   DEF_DA_ENS_NUM+1, landpatch, sag_ens         )
      CALL ncio_read_vector(file_restart, 'scv     '   ,                   DEF_DA_ENS_NUM+1, landpatch, scv_ens         )
      CALL ncio_read_vector(file_restart, 'snowdp  '   ,                   DEF_DA_ENS_NUM+1, landpatch, snowdp_ens      )
      CALL ncio_read_vector(file_restart, 'fveg    '   ,                   DEF_DA_ENS_NUM+1, landpatch, fveg_ens        )
      CALL ncio_read_vector(file_restart, 'fsno    '   ,                   DEF_DA_ENS_NUM+1, landpatch, fsno_ens        )
      CALL ncio_read_vector(file_restart, 'sigf    '   ,                   DEF_DA_ENS_NUM+1, landpatch, sigf_ens        )
      CALL ncio_read_vector(file_restart, 'green   '   ,                   DEF_DA_ENS_NUM+1, landpatch, green_ens       )
      CALL ncio_read_vector(file_restart, 'lai     '   ,                   DEF_DA_ENS_NUM+1, landpatch, lai_ens         )
      CALL ncio_read_vector(file_restart, 'tlai    '   ,                   DEF_DA_ENS_NUM+1, landpatch, tlai_ens        )
      CALL ncio_read_vector(file_restart, 'sai     '   ,                   DEF_DA_ENS_NUM+1, landpatch, sai_ens         )
      CALL ncio_read_vector(file_restart, 'tsai    '   ,                   DEF_DA_ENS_NUM+1, landpatch, tsai_ens        )
      CALL ncio_read_vector(file_restart, 'alb     '   , 2, 2,             DEF_DA_ENS_NUM+1, landpatch, alb_ens         )
      CALL ncio_read_vector(file_restart, 'ssun    '   , 2, 2,             DEF_DA_ENS_NUM+1, landpatch, ssun_ens        )
      CALL ncio_read_vector(file_restart, 'ssha    '   , 2, 2,             DEF_DA_ENS_NUM+1, landpatch, ssha_ens        )
      CALL ncio_read_vector(file_restart, 'ssoi    '   , 2, 2,             DEF_DA_ENS_NUM+1, landpatch, ssoi_ens        )
      CALL ncio_read_vector(file_restart, 'ssno    '   , 2, 2,             DEF_DA_ENS_NUM+1, landpatch, ssno_ens        )
      CALL ncio_read_vector(file_restart, 'thermk  '   ,                   DEF_DA_ENS_NUM+1, landpatch, thermk_ens      )
      CALL ncio_read_vector(file_restart, 'extkb   '   ,                   DEF_DA_ENS_NUM+1, landpatch, extkb_ens       )
      CALL ncio_read_vector(file_restart, 'extkd   '   ,                   DEF_DA_ENS_NUM+1, landpatch, extkd_ens       )
      CALL ncio_read_vector(file_restart, 'zwt     '   ,                   DEF_DA_ENS_NUM+1, landpatch, zwt_ens         )
      CALL ncio_read_vector(file_restart, 'wdsrf   '   ,                   DEF_DA_ENS_NUM+1, landpatch, wdsrf_ens       )
      CALL ncio_read_vector(file_restart, 'wa      '   ,                   DEF_DA_ENS_NUM+1, landpatch, wa_ens          )
      CALL ncio_read_vector(file_restart, 'wetwat  '   ,                   DEF_DA_ENS_NUM+1, landpatch, wetwat_ens      )
      CALL ncio_read_vector(file_restart, 't_lake  '   , nl_lake,          DEF_DA_ENS_NUM+1, landpatch, t_lake_ens      )
      CALL ncio_read_vector(file_restart, 'lake_icefrc', nl_lake,          DEF_DA_ENS_NUM+1, landpatch, lake_icefrac_ens)
      CALL ncio_read_vector(file_restart, 'savedtke1  ',                   DEF_DA_ENS_NUM+1, landpatch, savedtke1_ens   )

#ifdef GridRiverLakeFlow
      file_gridriver = trim(dir_restart)//'/'//trim(cdate)//'/'//trim(site)//'_restart_gridriver_'//trim(cdate)//'_lc'//trim(cyear)//'.nc'

      fexists = .false.
      IF (p_is_master) fexists = ncio_var_exist (file_gridriver, 'wdsrf_ucat_mean', readflag = .false.)
#ifdef USEMPI
      CALL mpi_bcast (fexists, 1, MPI_LOGICAL, p_address_master, p_comm_glb, p_err)
#endif

      IF (fexists) THEN
         CALL vector_read_and_scatter (file_gridriver, wdsrf_ucat_mean, numucat, 'wdsrf_ucat_mean', ucat_data_address)
         CALL vector_read_and_scatter (file_gridriver, veloc_riv_mean,  numucat, 'veloc_riv_mean',  ucat_data_address)
      ENDIF

      IF (DEF_Reservoir_Method > 0) THEN
         fexists = .false.
         IF (p_is_master) fexists = ncio_var_exist (file_gridriver, 'volresv_mean', readflag = .false.)
#ifdef USEMPI
         CALL mpi_bcast (fexists, 1, MPI_LOGICAL, p_address_master, p_comm_glb, p_err)
#endif
         IF (fexists) THEN
            CALL vector_read_and_scatter (file_gridriver, volresv_mean, numresv, 'volresv_mean', resv_data_address)
         ENDIF
      ENDIF
#endif

#ifdef RangeCheck
      CALL check_DATimeVariables
#endif

      IF (p_is_master) THEN
         write (*, *) 'Loading DA Time Variables done.'
      END IF

   END SUBROUTINE READ_DATimeVariables

#ifdef RangeCheck
!-----------------------------------------------------------------------

   SUBROUTINE check_DATimeVariables ()

!-----------------------------------------------------------------------
      IMPLICIT NONE

#ifdef USEMPI
      CALL mpi_barrier(p_comm_glb, p_err)
#endif
      IF (p_is_master) THEN
         write (*, *) 'Checking DA Time Variables ...'
      END IF

      CALL check_vector_data ('z_sno       [m]    ', z_sno_ens       )
      CALL check_vector_data ('dz_sno      [m]    ', dz_sno_ens      )
      CALL check_vector_data ('t_soisno    [K]    ', t_soisno_ens    )
      CALL check_vector_data ('wliq_soisno [kg/m2]', wliq_soisno_ens )
      CALL check_vector_data ('wice_soisno [kg/m2]', wice_soisno_ens )
      CALL check_vector_data ('smp         [mm]   ', smp_ens         )
      CALL check_vector_data ('hk          [mm/s] ', hk_ens          )
      CALL check_vector_data ('t_grnd      [K]    ', t_grnd_ens      )
      CALL check_vector_data ('tleaf       [K]    ', tleaf_ens       )
      CALL check_vector_data ('ldew        [mm]   ', ldew_ens        )
      CALL check_vector_data ('ldew_rain   [mm]   ', ldew_rain_ens   )
      CALL check_vector_data ('ldew_snow   [mm]   ', ldew_snow_ens   )
      CALL check_vector_data ('fwet_snow   [-]    ', fwet_snow_ens   )
      CALL check_vector_data ('sag         [-]    ', sag_ens         )
      CALL check_vector_data ('scv         [mm]   ', scv_ens         )
      CALL check_vector_data ('snowdp      [m]    ', snowdp_ens      )
      CALL check_vector_data ('fveg        [-]    ', fveg_ens        )
      CALL check_vector_data ('fsno        [-]    ', fsno_ens        )
      CALL check_vector_data ('sigf        [-]    ', sigf_ens        )
      CALL check_vector_data ('green       [-]    ', green_ens       )
      CALL check_vector_data ('lai         [-]    ', lai_ens         )
      CALL check_vector_data ('tlai        [-]    ', tlai_ens        )
      CALL check_vector_data ('sai         [-]    ', sai_ens         )
      CALL check_vector_data ('tsai        [-]    ', tsai_ens        )
      CALL check_vector_data ('alb         [-]    ', alb_ens         )
      CALL check_vector_data ('ssun        [-]    ', ssun_ens        )
      CALL check_vector_data ('ssha        [-]    ', ssha_ens        )
      CALL check_vector_data ('ssoi        [-]    ', ssoi_ens        )
      CALL check_vector_data ('ssno        [-]    ', ssno_ens        )
      CALL check_vector_data ('thermk      [-]    ', thermk_ens      )
      CALL check_vector_data ('extkb       [-]    ', extkb_ens       )
      CALL check_vector_data ('extkd       [-]    ', extkd_ens       )
      CALL check_vector_data ('zwt         [m]    ', zwt_ens         )
      CALL check_vector_data ('wdsrf       [mm]   ', wdsrf_ens       )
      CALL check_vector_data ('wa          [mm]   ', wa_ens          )
      CALL check_vector_data ('wetwat      [mm]   ', wetwat_ens      )
      CALL check_vector_data ('t_lake      [K]    ', t_lake_ens      )
      CALL check_vector_data ('lake_icefrc [-]    ', lake_icefrac_ens)
      CALL check_vector_data ('savedtke1   [W/m K]', savedtke1_ens   )
#ifdef GridRiverLakeFlow
      CALL check_vector_data ('wdsrf_ucat_mean [m]  ', wdsrf_ucat_mean)
      CALL check_vector_data ('veloc_riv_mean  [m/s]', veloc_riv_mean )
      CALL check_vector_data ('volresv_mean    [m^3]', volresv_mean   )
#endif

#ifdef USEMPI
      CALL mpi_barrier(p_comm_glb, p_err)
#endif

   END SUBROUTINE check_DATimeVariables
#endif

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Vars_TimeVariables
#endif
