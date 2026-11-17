#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Operator_EM
!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Microwave emission model adapted from the Community Microwave Emission
!   Modelling Platform (CMEM) v6.1. Simulates horizontally and vertically
!   polarized top-of-atmosphere brightness temperatures from land surface
!   states, microwave frequency, and incidence angle.
!   The surface layers follow soil -> snow -> vegetation, with no burial
!   branch. This configuration targets thin snow and low vegetation.
!
! CONVENTION:
!   Complex permittivities use negative imaginary parts and literature 
!   formulas using positive imaginary parts are complex-conjugated.
!
! REFERENCES:
!   [1] Ulaby, F. T., Moore, R. K., & Fung, A. K. (1986). Microwave remote
!       sensing: Active and passive, Vol. III: From theory to applications.
!       Artech House, Norwood, MA.
!
!   [2] Klein, L. A., & Swift, C. T. (1977). An improved model for the
!       dielectric constant of sea water at microwave frequencies. IEEE
!       Transactions on Antennas and Propagation, 25(1), 104-111.
!       https://doi.org/10.1109/TAP.1977.1141539
!
!   [3] Mätzler, C. (Ed.). (2006). Thermal microwave radiation: Applications
!       for remote sensing. IET Electromagnetic Waves Series 52. Institution
!       of Engineering and Technology, London.
!       https://doi.org/10.1049/PBEW052E
!
!   [4] Wang, J. R., & Schmugge, T. J. (1980). An empirical model for the
!       complex dielectric permittivity of soils as a function of water
!       content. IEEE Transactions on Geoscience and Remote Sensing, GE-18(4),
!       288-295.
!       https://doi.org/10.1109/TGRS.1980.350304
!
!   [5] Dobson, M. C., Ulaby, F. T., Hallikainen, M. T., & El-Rayes, M. A.
!       (1985). Microwave dielectric behavior of wet soil-Part II: Dielectric
!       mixing models. IEEE Transactions on Geoscience and Remote Sensing,
!       GE-23(1), 35-46.
!       https://doi.org/10.1109/TGRS.1985.289498
!
!   [6] Peplinski, N. R., Ulaby, F. T., & Dobson, M. C. (1995a). Dielectric
!       properties of soils in the 0.3-1.3-GHz range. IEEE Transactions on
!       Geoscience and Remote Sensing, 33(3), 803-807.
!       https://doi.org/10.1109/36.387598
!
!   [7] Peplinski, N. R., Ulaby, F. T., & Dobson, M. C. (1995b). Corrections
!       to "Dielectric properties of soils in the 0.3-1.3-GHz range". IEEE
!       Transactions on Geoscience and Remote Sensing, 33(6), 1340.
!       https://doi.org/10.1109/TGRS.1995.477193
!
!   [8] Mironov, V. L., Dobson, M. C., Kaupp, V. H., Komarov, S. A., &
!       Kleshchenko, V. N. (2004). Generalized refractive mixing dielectric
!       model for moist soils. IEEE Transactions on Geoscience and Remote
!       Sensing, 42(4), 773-785.
!       https://doi.org/10.1109/TGRS.2003.823288
!
!   [9] Mironov, V. L., Kosolapova, L. G., & Fomin, S. V. (2009). Physically
!       and mineralogically based spectroscopic dielectric model for moist
!       soils. IEEE Transactions on Geoscience and Remote Sensing, 47(7),
!       2059-2070.
!       https://doi.org/10.1109/TGRS.2008.2011631
!
!   [10] Mironov, V. L., & Fomin, S. V. (2009). Temperature and mineralogy
!       dependable model for microwave dielectric spectra of moist soils.
!       PIERS Online, 5(5), 411-415.
!       https://doi.org/10.2529/PIERS090220054025
!
!   [11] Pellarin, T., Wigneron, J.-P., Calvet, J.-C., Berger, M., Douville,
!       H., Ferrazzoli, P., Kerr, Y. H., Lopez-Baeza, E., Pulliainen, J.,
!       Simmonds, L. P., & Waldteufel, P. (2003). Two-year global simulation
!       of L-band brightness temperatures over land. IEEE Transactions on
!       Geoscience and Remote Sensing, 41(9), 2135-2139.
!       https://doi.org/10.1109/TGRS.2003.815417
!
!   [12] Wigneron, J.-P., Kerr, Y. H., Waldteufel, P., Saleh, K., Escorihuela,
!       M.-J., Richaume, P., Ferrazzoli, P., de Rosnay, P., Gurney, R.,
!       Calvet, J.-C., Grant, J. P., Guglielmetti, M., Hornbuckle, B.,
!       Mätzler, C., Pellarin, T., & Schwank, M. (2007). L-band microwave
!       emission of the biosphere (L-MEB) model: Description and calibration
!       against experimental data sets over crop fields. Remote Sensing of
!       Environment, 107(4), 639-655.
!       https://doi.org/10.1016/j.rse.2006.10.014
!
!   [13] Mätzler, C. (1998). Microwave permittivity of dry sand. IEEE
!       Transactions on Geoscience and Remote Sensing, 36(1), 317-319.
!       https://doi.org/10.1109/36.655342
!
!   [14] Wigneron, J.-P., Jackson, T. J., O'Neill, P., De Lannoy, G., de
!       Rosnay, P., Walker, J. P., Ferrazzoli, P., Mironov, V., Bircher, S.,
!       Grant, J. P., Kurum, M., Schwank, M., Muñoz-Sabater, J., Das, N.,
!       Royer, A., Al-Yaari, A., Al Bitar, A., Fernandez-Moran, R., Lawrence,
!       H., ... Kerr, Y. (2017). Modelling the passive microwave signature
!       from land surfaces: A review of recent results and application to the
!       L-band SMOS & SMAP soil moisture retrieval algorithms. Remote Sensing
!       of Environment, 192, 238-262.
!       https://doi.org/10.1016/j.rse.2017.01.024
!
!   [15] Njoku, E. G., & Kong, J.-A. (1977). Theory for passive microwave
!       remote sensing of near-surface soil moisture. Journal of Geophysical
!       Research, 82(20), 3108-3118.
!       https://doi.org/10.1029/JB082i020p03108
!
!   [16] Kerr, Y. H., & Njoku, E. G. (1990). A semiempirical model for
!       interpreting microwave emission from semiarid land surfaces as seen
!       from space. IEEE Transactions on Geoscience and Remote Sensing, 28(3),
!       384-393.
!       https://doi.org/10.1109/36.54364
!
!   [17] Grody, N. C., & Weng, F. (2008). Microwave emission and scattering
!       from deserts: Theory compared with satellite measurements. IEEE
!       Transactions on Geoscience and Remote Sensing, 46(2), 361-375.
!       https://doi.org/10.1109/TGRS.2007.909920
!
!   [18] Mätzler, C. (1987). Applications of the interaction of microwaves
!       with the natural snow cover. Remote Sensing Reviews, 2(2), 259-387.
!       https://doi.org/10.1080/02757258709532086
!
!   [19] Anderson, E. A. (1976). A point energy and mass balance model of a
!       snow cover. NOAA Technical Report NWS 19. National Weather Service,
!       Silver Spring, MD.
!
!   [20] Ulaby, F. T., Moore, R. K., & Fung, A. K. (1981). Microwave remote
!       sensing: Active and passive, Vol. I: Microwave remote sensing
!       fundamentals and radiometry. Addison-Wesley, Reading, MA.
!
!   [21] Wigneron, J.-P., Laguerre, L., & Kerr, Y. H. (2001). A simple
!       parameterization of the L-band microwave emission from rough
!       agricultural soils. IEEE Transactions on Geoscience and Remote
!       Sensing, 39(8), 1697-1707.
!       https://doi.org/10.1109/36.942548
!
!   [22] Pulliainen, J. T., Grandell, J., & Hallikainen, M. T. (1999). HUT
!       snow emission model and its applicability to snow water equivalent
!       retrieval. IEEE Transactions on Geoscience and Remote Sensing, 37(3),
!       1378-1390.
!       https://doi.org/10.1109/36.763302
!
!   [23] Hallikainen, M. T., Ulaby, F. T., & Van Deventer, T. E. (1987).
!       Extinction behavior of dry snow in the 18- to 90-GHz range. IEEE
!       Transactions on Geoscience and Remote Sensing, GE-25(6), 737-745.
!       https://doi.org/10.1109/TGRS.1987.289743
!
!   [24] Orfanidis, S. J. (2016). Electromagnetic waves and antennas.
!       Online textbook. Rutgers University, New Brunswick, NJ.
!       https://eceweb1.rutgers.edu/~orfanidi/ewa/
!
! HISTORY:
!   Lu Li, 12/2024: Initial version
!   Zhilong Fan, Lu Li, 03/2024: Debug and clean codes
!   Lu Li, 10/2025: Debug and clean codes
!   Lu Li, Zhilong Fan, 11/2026: DA v1.0
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Const_Physical
   USE MOD_Vars_Global, only: pi, nl_soil, N_land_classification, maxsnl, spval, dz_soi, zi_soi
   USE MOD_Namelist
   USE MOD_DA_Operator_Const
   IMPLICIT NONE

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE EM(                               &
      patchtype, patchclass, dz_sno,            &
      forc_topo, htop,                          &
      tref, t_soisno, tleaf,                    &
      wliq_soisno, wice_soisno,                 &
      snowdp, lai, sai,                         &
      wf_clay, wf_sand, wf_silt, BD_all, porsl, &
      sat_theta, sat_fghz,                      &
      tb_toa_h, tb_toa_v)

!-----------------------------------------------------------------------------
      integer, intent(in)      :: patchtype                         ! Land cover type
      integer, intent(in)      :: patchclass                        ! Land cover class
      real(r8), intent(in)     :: dz_sno(maxsnl + 1:0)              ! Layer thickness (m)
      real(r8), intent(in)     :: forc_topo                         ! Topography [m]
      real(r8), intent(in)     :: htop                              ! Retained for caller compatibility [m]
      real(r8), intent(in)     :: tref                              ! 2 m height air temperature [kelvin]
      real(r8), intent(in)     :: tleaf                             ! Leaf temperature [K]
      real(r8), intent(in)     :: t_soisno(maxsnl + 1:nl_soil)      ! Snow and soil layer temperature [K]
      real(r8), intent(in)     :: wliq_soisno(maxsnl + 1:nl_soil)   ! Liquid water mass per unit area in snow and soil layers [kg/m2]
      real(r8), intent(in)     :: wice_soisno(maxsnl + 1:nl_soil)   ! Ice mass per unit area in snow and soil layers [kg/m2]
      real(r8), intent(in)     :: snowdp                            ! Snow depth [meter]
      real(r8), intent(in)     :: lai                               ! Leaf area index
      real(r8), intent(in)     :: sai                               ! Stem area index
      real(r8), intent(in)     :: wf_clay(nl_soil)                  ! Gravimetric fraction of clay
      real(r8), intent(in)     :: wf_sand(nl_soil)                  ! Gravimetric fraction of sand
      real(r8), intent(in)     :: wf_silt(nl_soil)                  ! Gravimetric fraction of silt
      real(r8), intent(in)     :: BD_all(nl_soil)                   ! Bulk density of soil (GRAVELS + ORGANIC MATTER + Mineral Soils,kg/m3)
      real(r8), intent(in)     :: porsl(nl_soil)                    ! Retained for caller compatibility [-]
      real(r8), intent(in)     :: sat_theta                         ! Incidence angle of satellite (rad)
      real(r8), intent(in)     :: sat_fghz                          ! Frequency of satellite (GHz)
      real(r8), intent(out)    :: tb_toa_h                          ! Brightness temperature of top-of-atmosphere for H- polarization
      real(r8), intent(out)    :: tb_toa_v                          ! Brightness temperature of top-of-atmosphere for V- polarization

      real(r8)                 :: dz_soisno(maxsnl+1:nl_soil)       ! Snow and soil layer thickness [m]
      integer                  :: lb                                ! Lower bound of arrays
      real(r8)                 :: tau_atm                           ! Atmospheric optical depth
      real(r8)                 :: r_soil_rough(2)                   ! Rough surface reflectivity for H and V polarizations
      real(r8)                 :: r_snow(2)                         ! Reflectivity of the snow for H and V polarizations
      real(r8)                 :: tb_soil(2)                        ! Brightness temperature of soil for H and V polarizations
      real(r8)                 :: tb_tos(2)                         ! Brightness temperature of snow-covered ground for H and V polarizations
      real(r8)                 :: tb_tov(2)                         ! Brightness temperature of vegetation (consider snow) for H and V polarizations
      real(r8)                 :: tb_atm_up(2)                      ! Upwelling radiation (brightness temperature) of atmosphere
      real(r8)                 :: tb_atm_down(2)                    ! Downwelling radiation (brightness temperature) of atmosphere
      real(r8)                 :: rho_snow                          ! Snow density [kg/m3]
      real(r8)                 :: vliq_snow                         ! Snow liquid water volume fraction [0-1]
      real(r8)                 :: gamma_veg(2)                      ! Canopy transmissivity for H and V polarizations
      real(r8)                 :: tb_veg(2)                         ! Brightness temperature of vegetation for H- and V- polarization
      real(r8)                 :: tb_toa(2)                         ! Brightness temperature of top-of-atmosphere for H- and V- polarization
      real(r8)                 :: wf_total(nl_soil)                 ! Total gravimetric
      real(r8)                 :: BD_all_surf                       ! Surface soil bulk density [kg/m3]
      real(r8)                 :: t_soil_surf                       ! Soil temperature at surface [K]
      real(r8)                 :: t_soil_deep                       ! Soil temperature at deep layer [K]
      real(r8)                 :: vliq_soil                         ! Liquid volumetric water content at surface (m3/m3)
      real(r8)                 :: vice_soil                         ! Ice volumetric water content at surface (m3/m3)
      real(r8)                 :: wf_clay_surf                      ! Clay mass fraction at surface [0-1]
      real(r8)                 :: wf_sand_surf                      ! Sand mass fraction at surface [0-1]
      integer                  :: i
!-----------------------------------------------------------------------------
!#############################################################################
! Prepare parameters & states used in the operator.
!#############################################################################
      ! Get depth of soil and snow layers
      dz_soisno(maxsnl+1:0) = dz_sno(maxsnl+1:0)
      dz_soisno(1:nl_soil) = dz_soi(1:nl_soil)

      ! Calculate thickness-weighted parameters over the first two soil layers.
      wf_total     = wf_clay + wf_sand + wf_silt
      wf_clay_surf = (wf_clay(1)/wf_total(1)*dz_soi(1)+wf_clay(2)/wf_total(2)*dz_soi(2)) / sum(dz_soi(1:2))
      wf_sand_surf = (wf_sand(1)/wf_total(1)*dz_soi(1)+wf_sand(2)/wf_total(2)*dz_soi(2)) / sum(dz_soi(1:2))
      BD_all_surf  = (BD_all(1)*dz_soi(1)+BD_all(2)*dz_soi(2)) / sum(dz_soi(1:2))

      ! Average surface temperature over the first two layers and deep temperature over 0.5-1 m.
      t_soil_surf = (t_soisno(1)*dz_soi(1)+t_soisno(2)*dz_soi(2)) / sum(dz_soi(1:2))
      t_soil_deep = (t_soisno(7)*(zi_soi(7)-0.5)+t_soisno(8)*(1.0-zi_soi(7))) / 0.5

      ! Convert liquid and ice masses to volume fractions over the first two layers.
      vliq_soil_surf = max((wliq_soisno(1)+wliq_soisno(2))/(sum(dz_soi(1:2))*denh2o),0)
      vice_soil_surf = max((wice_soisno(1)+wice_soisno(2))/(sum(dz_soi(1:2))*denice),0)

      ! Calculate lower bound of snow
      lb = 0
      DO i = maxsnl+1, 0
         IF (wliq_soisno(i) + wice_soisno(i) > 0.0) THEN
            lb = i
            EXIT
         ENDIF
      ENDDO

!#############################################################################
! Run the forward operator.
!#############################################################################
      IF (patchtype == 0) THEN

!-----------------------------------------------------------------------------
! Atmosphere module.
!-----------------------------------------------------------------------------
         CALL atm(sat_theta, forc_topo, tref, tau_atm, tb_atm_up, tb_atm_down)
!-----------------------------------------------------------------------------
! Soil module.
!-----------------------------------------------------------------------------
         CALL soil( &
            sat_fghz, sat_theta,                                 &
            patchclass,                                          &
            t_soil_surf, t_soil_deep,                            &
            vliq_soil_surf, vice_soil_surf,                      &
            wf_sand_surf, wf_clay_surf, BD_all_surf,             &
            teff_soil, r_soil_rough, tb_soil)
!-----------------------------------------------------------------------------
! Snow module.
!-----------------------------------------------------------------------------
         IF (snowdp > 0.01) THEN
            rho_snow = (wliq_soisno(lb)+wice_soisno(lb)) / dz_soisno(lb)
            vliq_snow = wliq_soisno(lb) / (dz_soisno(lb)*denh2o)
            CALL snow(sat_fghz, sat_theta, t_soisno(lb), teff_soil, snowdp, rho_snow, vliq_snow, r_soil_rough, r_snow, tb_tos)
         ELSE
            r_snow = r_soil_rough
            tb_tos = tb_soil
         ENDIF
!-----------------------------------------------------------------------------
! Vegetation module.
!-----------------------------------------------------------------------------
         IF (lai + sai > 1e-6) THEN
            CALL veg(sat_theta, patchclass, lai, tleaf, tb_veg, gamma_veg)
            tb_tov = tb_veg + tb_veg * gamma_veg * r_snow + tb_tos * gamma_veg + tb_atm_down * r_snow * gamma_veg**2
         ELSE
            tb_tov = tb_tos + tb_atm_down*r_snow
         ENDIF

         ! Calculate brightness temperature of top-of-atmosphere.
         tb_toa = tb_tov*exp(-tau_atm/cos(sat_theta)) + tb_atm_up
      ELSE
         tb_toa = spval
      ENDIF

      tb_toa_h = tb_toa(1)
      tb_toa_v = tb_toa(2)

   END SUBROUTINE EM

!-----------------------------------------------------------------------------

   SUBROUTINE atm(theta, z, tref, tau_atm, tb_atm_up, tb_atm_down)

!-----------------------------------------------------------------------------
      real(r8), intent(in)     :: theta                             ! Incidence angle [rad]
      real(r8), intent(in)     :: z                                 ! Altitude (m)
      real(r8), intent(in)     :: tref                              ! 2m air temperature (K)
      real(r8), intent(out)    :: tau_atm                           ! Atmospheric optical depth
      real(r8), intent(out)    :: tb_atm_up(2)                      ! Upwelling radiation (brightness temperature) of atmosphere
      real(r8), intent(out)    :: tb_atm_down(2)                    ! Downwelling radiation (brightness temperature) of atmosphere

      real(r8)                 :: tb_cosmic = 2.7                   ! Cosmic background brightness temperature [K]
      real(r8)                 :: t_atm_eq                          ! Equivalent layer temperature [K]
      real(r8)                 :: gamma_atm                         ! Atmospheric transmissity [0-1]
!-----------------------------------------------------------------------------
      ! Atmospheric optical depth ([11], Eq.A1)
      tau_atm = exp(-3.9262-0.2211*z/1000-0.00369*tref)
      gamma_atm = exp(-tau_atm/cos(theta))

      ! Calculate equivalent layer temperature ([11], Eq.A1)
      t_atm_eq = exp(4.9274+0.002195*tref)

      ! Upwelling radiation (brightness temperature) of atmosphere ([11], Eq.A2)
      tb_atm_up(:) = t_atm_eq * (1.-gamma_atm)

      ! Downwelling atmospheric and attenuated cosmic emission ([11], Eq.A2)
      tb_atm_down(:) = t_atm_eq*(1.-gamma_atm) + tb_cosmic*gamma_atm

   END SUBROUTINE atm

!-----------------------------------------------------------------------------

   SUBROUTINE soil( &
      fghz, theta,                                         &
      patchclass,                                          &
      t_soil_surf, t_soil_deep,                            &
      vliq_soil_surf, vice_soil_surf,                      &
      wf_sand_surf, wf_clay_surf, BD_all_surf,             &
      teff_soil, r_soil_rough, tb_soil)

!-----------------------------------------------------------------------------
      real(r8), intent(in)     :: fghz                              ! Microwave frequency [GHz]
      real(r8), intent(in)     :: theta                             ! Incidence angle [rad]
      integer, intent(in)      :: patchclass                        ! Land cover class
      real(r8), intent(in)     :: t_soil_surf                       ! Soil temperature at surface [K]
      real(r8), intent(in)     :: t_soil_deep                       ! Soil temperature at deep layer [K]
      real(r8), intent(in)     :: vliq_soil_surf                    ! Liquid volumetric water content at surface (m3/m3)
      real(r8), intent(in)     :: vice_soil_surf                    ! Ice volumetric water content at surface (m3/m3)
      real(r8), intent(in)     :: wf_sand_surf                      ! Sand mass fraction at surface [0-1]
      real(r8), intent(in)     :: wf_clay_surf                      ! Clay mass fraction at surface [0-1]
      real(r8), intent(in)     :: BD_all_surf                       ! Surface soil bulk density [kg/m3]
      real(r8), intent(out)    :: teff_soil(2)                      ! Effective temperature for H and V polarizations, [K]
      real(r8), intent(out)    :: r_soil_rough(2)                   ! Rough surface reflectivity for H and V polarizations
      real(r8), intent(out)    :: tb_soil(2)                        ! Brightness temperature of soil

      complex(r8)              :: eps_soil                          ! Dielectric constant of soil for H and V polarizations
      real(r8)                 :: r_soil_smooth(2)                  ! Smooth surface reflectivity for H and V polarizations
      logical                  :: is_desert                         ! Flag for desert soil
      real(r8)                 :: sal_soil = 0.0                    ! Soil-water salt mass fraction [g/kg]
      complex(r8)              :: eps_ice = (3.2, -0.1)             ! Dielectric constant of ice
!-----------------------------------------------------------------------------
      ! Whether this patch is desert
      is_desert = .false.
      IF (vliq_soil_surf < 0.02 .and. wf_sand_surf > 0.90) THEN
         is_desert = .true.
      ENDIF

      ! Calculate effective temperature
      CALL cal_teff_soil(vliq_soil_surf, t_soil_surf, t_soil_deep, teff_soil)

      ! Calculate dielectric constant of soil (mixture medium)
      IF (is_desert) THEN
         ! Dry-sand permittivity ([13], Eq.1)
         eps_soil = 2.53 + (2.79-2.53)/(1-jj*(fghz/0.27)) + jj*0.002
         eps_soil = conjg(eps_soil)
      ELSE
         ! Caculate dielectric constant in mixed soil (only liquid water)
         IF (DEF_DA_OPERATOR_EM_DIEL == 1) THEN
            CALL diel_soil_D85(fghz, t_soil_surf, vliq_soil_surf/(1-vice_soil_surf), wf_sand_surf, wf_clay_surf, BD_all_surf, sal_soil, eps_soil)
         ELSEIF (DEF_DA_OPERATOR_EM_DIEL == 2) THEN
            CALL diel_soil_M04(fghz,  vliq_soil_surf/(1-vice_soil_surf), wf_clay_surf, eps_soil)
         ELSEIF (DEF_DA_OPERATOR_EM_DIEL == 3) THEN
            CALL diel_soil_M09(fghz,  vliq_soil_surf/(1-vice_soil_surf), t_soil_surf, wf_clay_surf, eps_soil)
         ENDIF
      ENDIF

      ! Mix dielectric constant of frozen and non-frozen soil
      eps_soil = eps_soil*(1.-vice_soil_surf) + eps_ice*vice_soil_surf

      ! Calculate smooth surface reflectivity
      CALL cal_smooth_refl(theta, eps_soil, r_soil_smooth)

      ! Calculate surface reflectivity and emissivity
      CALL cal_rough_refl(fghz, theta, is_desert, patchclass, eps_soil, r_soil_smooth, r_soil_rough)

      ! Calculate brightness temperature
      tb_soil = teff_soil * (1-r_soil_rough)

   END SUBROUTINE soil

!-----------------------------------------------------------------------------

   SUBROUTINE cal_teff_soil(vliq_soil_surf, t_soil_surf, t_soil_deep, teff_soil)

!-----------------------------------------------------------------------------
      real(r8), intent(in)     :: vliq_soil_surf                    ! Liquid volumetric soil water content [m3/m3]
      real(r8), intent(in)     :: t_soil_surf                       ! Surface soil temperature [K]
      real(r8), intent(in)     :: t_soil_deep                       ! Deep soil temperature [K]
      real(r8), intent(out)    :: teff_soil(2)                      ! Effective temperature for H and V polarizations, [K]

      real(r8)                 :: C                                 ! Moisture-dependent surface-temperature weight
      real(r8)                 :: w0 = 0.30                         ! Reference soil moisture [m3/m3]
      real(r8)                 :: bw = 0.30                         ! Moisture-weighting exponent
!-----------------------------------------------------------------------------
      ! Calculate moisture-dependent weight ([21], Eq.7)
      C = min(1.0, (vliq_soil_surf/w0)**bw)

      ! Calculate effective soil temperature ([21], Eq.8)
      teff_soil(:) = t_soil_deep + (t_soil_surf-t_soil_deep)*C

   END SUBROUTINE cal_teff_soil

!-----------------------------------------------------------------------------

   SUBROUTINE cal_diel_ice(fghz, t_ice, eps_ice)

!-----------------------------------------------------------------------------
      real(r8), intent(in)     :: fghz                              ! Microwave frequency [GHz]
      real(r8), intent(in)     :: t_ice                             ! Ice temperature [K]
      complex(r8), intent(out) :: eps_ice                           ! Complex dielectric permittivity of ice

      real(r8)                 :: t_ice_c                           ! Ice temperature for empirical fits [degC]
      real(r8)                 :: eps_ice_r                         ! Real part of pure ice dielectric constant
      real(r8)                 :: alpha                             ! Alpha parameter
      real(r8)                 :: beta                              ! Beta parameter
      real(r8)                 :: eps_ice_i                         ! Positive loss factor of pure ice
!-----------------------------------------------------------------------------
      t_ice_c = t_ice - tfrz

      ! Calculate real part of pure ice dielectric constant ([3], Eq.5.Q.1)
      eps_ice_r = 3.1884 + 9.1e-4*t_ice_c

      ! Calculate beta parameter ([3], Eq.5.Q.4 - Eq.5.Q.6)
      beta = (0.0207/t_ice)*(exp(335./t_ice)/((exp(335./t_ice)-1.)**2.)) + 1.16e-11*(fghz**2.) + exp(-9.963+0.0372*t_ice_c)

      ! Calculate alpha parameter ([3], Eq.5.Q.3)
      alpha = (0.00504+0.0062*(300./t_ice-1)) * exp(-22.1*(300./t_ice-1))

      ! Imaginary part of pure-ice permittivity ([3], Eq.5.Q.2)
      eps_ice_i = alpha/fghz + beta*fghz

      ! Calculate dielectric constant of pure ice ([3]).
      eps_ice = eps_ice_r - jj*eps_ice_i

   END SUBROUTINE cal_diel_ice

!-----------------------------------------------------------------------------

   SUBROUTINE cal_diel_water(fghz, vliq_soil, t_water, wf_sand, wf_clay, BD_all, sal_water, eps_water)

!-----------------------------------------------------------------------------
      real(r8), intent(in)     :: fghz                              ! Microwave frequency [GHz]
      real(r8), intent(in)     :: vliq_soil                         ! Liquid volumetric soil water content [m3/m3]
      real(r8), intent(in)     :: t_water                           ! Liquid-water temperature [K]
      real(r8), intent(in)     :: wf_sand                           ! Sand mass fraction [0-1]
      real(r8), intent(in)     :: wf_clay                           ! Clay mass fraction [0-1]
      real(r8), intent(in)     :: BD_all                            ! Dry soil bulk density [kg/m3]
      real(r8), intent(in)     :: sal_water                         ! Water salt mass fraction [g/kg]
      complex(r8), intent(out) :: eps_water                         ! Dielectric constant of soil water

      real(r8)                 :: sigma_water                       ! Ionic conductivity (S/m)
      real(r8)                 :: a, b                              ! Parameters
      real(r8)                 :: eps0_water                        ! Static dielectric constant of salinity water
      real(r8)                 :: tau_water                         ! Relaxation time of salinity water
      real(r8)                 :: omega                             ! Angular frequency [rad/s]
      real(r8)                 :: t_water_c                         ! Liquid-water temperature for empirical fits [degC]
      real(r8)                 :: sand_pct                          ! Sand content by mass [%]
      real(r8)                 :: clay_pct                          ! Clay content by mass [%]
      real(r8)                 :: rho_bulk_gcm3                     ! Dry soil bulk density for empirical fits [g/cm3]
!-----------------------------------------------------------------------------
      omega = 2.0 * pi * (fghz*1e9)

      ! Convert interface units for the empirical formulas.
      t_water_c = t_water - tfrz
      sand_pct = wf_sand * 100.0_r8
      clay_pct = wf_clay * 100.0_r8
      rho_bulk_gcm3 = BD_all / 1000.0_r8

      ! Salinity-dependent static permittivity ([2], Eqs.13-15).
      eps0_water = 87.134 - 1.949e-1*t_water_c - 1.276e-2*t_water_c**2 + 2.491e-4*t_water_c**3
      a = 1.000 + 1.613e-5*sal_water*t_water_c - 3.656e-3*sal_water + 3.210e-5*sal_water**2 - 4.232e-7*sal_water**3
      eps0_water = eps0_water * a

      ! Salinity-dependent relaxation time ([2], Eqs.16-18).
      tau_water = 1.768e-11 - 6.086e-13*t_water_c + 1.104e-14*t_water_c**2 - 8.111e-17*t_water_c**3
      b = 1.000 + 2.282e-5*sal_water*t_water_c - 7.638e-4*sal_water - 7.760e-6*sal_water**2 + 1.105e-8*sal_water**3
      tau_water = tau_water * b

      ! Soil-water conductivity ([5] Eq.32; [6] Eq.8).
      sigma_water = max(-1.645+1.939*rho_bulk_gcm3-0.0225622*sand_pct+0.01594*clay_pct,0)

      ! Calculate dielectric constant of soil-water by modified Debye expression ([5] Eq.29)
      eps_water = epsinf_water + (eps0_water-epsinf_water)/(1+jj*omega*tau_water)         &
         - jj*sigma_water/(omega*eps0)*(rho_soil_gcm3-rho_bulk_gcm3)/(rho_soil_gcm3*max(0.001,vliq_soil))

   END SUBROUTINE cal_diel_water

!-----------------------------------------------------------------------------

   SUBROUTINE cal_diel_soil_D85(fghz, t_soil, vliq_soil, wf_sand, wf_clay, BD_all, sal_soil, eps_soil)

!-----------------------------------------------------------------------------
      real(r8), intent(in)     :: fghz                              ! Microwave frequency [GHz]
      real(r8), intent(in)     :: t_soil                            ! Soil temperature [K]
      real(r8), intent(in)     :: vliq_soil                         ! Liquid volumetric soil water content [m3/m3]
      real(r8), intent(in)     :: wf_sand                           ! Sand mass fraction [0-1]
      real(r8), intent(in)     :: wf_clay                           ! Clay mass fraction [0-1]
      real(r8), intent(in)     :: BD_all                            ! Dry soil bulk density [kg/m3]
      real(r8), intent(in)     :: sal_soil                          ! Soil-water salt mass fraction [g/kg]
      complex(r8), intent(out) :: eps_soil                          ! Complex moist-soil relative permittivity [-]

      complex(r8)              :: eps_water                         ! Liquid-water relative permittivity
      real(r8)                 :: eps_bulk                          ! Solid soil-particle relative permittivity [-]
      real(r8), parameter      :: alphas = 0.65_r8                  ! Empirical dielectric mixing exponent [-]
      real(r8)                 :: beta_r                            ! Texture-dependent exponent for the real part [-]
      real(r8)                 :: beta_i                            ! Texture-dependent exponent for dielectric loss [-]
      real(r8)                 :: eps_soil_r                        ! Real part of moist-soil relative permittivity [-]
      real(r8)                 :: eps_soil_i                        ! Nonnegative moist-soil dielectric loss factor [-]
      real(r8)                 :: sand_pct                          ! Sand content by mass [%]
      real(r8)                 :: clay_pct                          ! Clay content by mass [%]
      real(r8)                 :: rho_bulk_gcm3                     ! Dry soil bulk density for empirical mixing [g/cm3]
!-----------------------------------------------------------------------------
      ! Convert interface units for the empirical formulas.
      sand_pct = wf_sand * 100.0_r8
      clay_pct = wf_clay * 100.0_r8
      rho_bulk_gcm3 = BD_all / 1000.0_r8

      ! Calculate liquid-water permittivity using the soil temperature
      CALL cal_diel_water(fghz, vliq_soil, t_soil, wf_sand, wf_clay, BD_all, sal_soil, eps_water)

      ! Calculate solid-particle permittivity ([5], Eq.22) and the real-part moisture exponent and mixing relation ([6], Eqs.4, 2)
      eps_bulk = (1.01_r8+0.44_r8*rho_bulk_gcm3)**2.0_r8 - 0.062_r8
      beta_r = (127.48_r8-0.519_r8*sand_pct-0.152_r8*clay_pct) / 100.0_r8
      eps_soil_r = (1.0_r8+(rho_bulk_gcm3/rho_soil_gcm3)*(eps_bulk**alphas-1.0_r8)         &
         +(max(vliq_soil,0.001)**beta_r)*(real(eps_water)**alphas)                     &
         -max(vliq_soil,0.001))**(1.0_r8/alphas) 

      ! Calculate the moisture exponent and mixing relation for dielectric loss ([6], Eqs.5, 3)
      beta_i = (133.797_r8-0.603_r8*sand_pct-0.166_r8*clay_pct) / 100.0_r8
      eps_soil_i = ((max(vliq_soil,0.001)**beta_i)*((abs(aimag(eps_water)))**alphas))**(1.0_r8/alphas)

      eps_soil = eps_soil_r - jj*eps_soil_i

   END SUBROUTINE cal_diel_soil_D85

!-----------------------------------------------------------------------------

   SUBROUTINE cal_diel_soil_M04(fghz, vliq_soil, wf_clay, eps_soil)

!-----------------------------------------------------------------------------
      real(r8), intent(in)     :: fghz                              ! Microwave frequency [GHz]
      real(r8), intent(in)     :: vliq_soil                         ! Liquid volumetric soil water content [m3/m3]
      real(r8), intent(in)     :: wf_clay                           ! Clay mass fraction [0-1]
      complex(r8), intent(out) :: eps_soil                          ! Complex moist-soil relative permittivity [-]

      real(r8)                 :: n_dsoil                           ! Dry-soil refractive index [-]
      real(r8)                 :: k_dsoil                           ! Dry-soil normalized attenuation coefficient [-]
      real(r8)                 :: vliq_bwater_max                   ! Maximum bound-water volume fraction [m3/m3]
      real(r8)                 :: eps0_bwater                       ! Static bound-water relative permittivity [-]
      real(r8)                 :: tau_bwater                        ! Bound-water relaxation time [s]
      real(r8)                 :: sigma_bwater                      ! Bound-water conductivity [S/m]
      real(r8)                 :: eps0_fwater                       ! Static free-water relative permittivity [-]
      real(r8)                 :: tau_fwater                        ! Free-water relaxation time [s]
      real(r8)                 :: sigma_fwater                      ! Free-water conductivity [S/m]
      real(r8)                 :: eps_bwater_r                      ! Real part of bound-water relative permittivity [-]
      real(r8)                 :: eps_bwater_i                      ! Bound-water dielectric loss factor [-]
      real(r8)                 :: eps_fwater_r                      ! Real part of free-water relative permittivity [-]
      real(r8)                 :: eps_fwater_i                      ! Free-water dielectric loss factor [-]
      real(r8)                 :: n_bwater                          ! Bound-water refractive index [-]
      real(r8)                 :: k_bwater                          ! Bound-water normalized attenuation coefficient [-]
      real(r8)                 :: n_fwater                          ! Free-water refractive index [-]
      real(r8)                 :: k_fwater                          ! Free-water normalized attenuation coefficient [-]
      real(r8)                 :: n_soil                            ! Moist-soil refractive index [-]
      real(r8)                 :: k_soil                            ! Moist-soil normalized attenuation coefficient [-]
      real(r8)                 :: eps_soil_r                        ! Real part of moist-soil relative permittivity [-]
      real(r8)                 :: eps_soil_i                        ! Moist-soil dielectric loss factor [-]
      real(r8)                 :: f_hz                              ! Microwave frequency [Hz]
!-----------------------------------------------------------------------------
      f_hz = fghz * 1e9

      ! Dry-soil refractive parameters ([9], Eq.17-18).
      n_dsoil = 1.634 - 0.539*wf_clay + 0.2748*wf_clay**2
      k_dsoil = 0.03952 - 0.04038*wf_clay

      ! Calculate the maximum bound-water volume fraction ([9], Eq.19).
      vliq_bwater_max = 0.02863 + 0.30673*wf_clay

      ! Calculate bound-water Debye parameters ([9], Eqs.20-22).
      eps0_bwater = 79.8 - 85.4*wf_clay + 32.7*wf_clay*wf_clay
      tau_bwater = 1.062e-11 + 3.450e-12*wf_clay
      sigma_bwater = 0.3112 + 0.467*wf_clay

      ! Set free-water Debye parameters ([9], Eqs.23-25).
      eps0_fwater = 100
      tau_fwater = 8.5e-12
      sigma_fwater = 0.3631 + 1.217*wf_clay

      ! Calculate bound- and free-water permittivities ([9], Eq.16).
      eps_bwater_r = epsinf_water + (eps0_bwater-epsinf_water)/(1.+(2.*pi*f_hz*tau_bwater)**2)
      eps_bwater_i = (eps0_bwater-epsinf_water)/(1.+(2.*pi*f_hz*tau_bwater)**2)*(2.*pi*f_hz*tau_bwater) + &
         sigma_bwater/(2.*pi*eps0*f_hz)
      eps_fwater_r = epsinf_water + (eps0_fwater-epsinf_water)/(1+(2*pi*f_hz*tau_fwater)**2)
      eps_fwater_i = (eps0_fwater-epsinf_water)/(1+(2*pi*f_hz*tau_fwater)**2)*(2.*pi*f_hz*tau_fwater) + &
         sigma_fwater/(2.*pi*eps0*f_hz)

      ! Calculate water refractive indices and normalized attenuation coefficients ([9], Eqs.14-15).
      n_bwater = sqrt(sqrt(eps_bwater_r**2+eps_bwater_i**2)+eps_bwater_r) / sqrt(2.0)
      k_bwater = sqrt(sqrt(eps_bwater_r**2+eps_bwater_i**2)-eps_bwater_r) / sqrt(2.0)
      n_fwater = sqrt(sqrt(eps_fwater_r**2+eps_fwater_i**2)+eps_fwater_r) / sqrt(2.0)
      k_fwater = sqrt(sqrt(eps_fwater_r**2+eps_fwater_i**2)-eps_fwater_r) / sqrt(2.0)

      ! Mix dry soil, bound water, and excess free water ([9], Eqs.12-13).
      IF (vliq_soil <= vliq_bwater_max) THEN
         n_soil = n_dsoil + (n_bwater-1)*vliq_soil
         k_soil = k_dsoil + k_bwater*vliq_soil
      ELSE
         n_soil = n_dsoil + (n_bwater-1)*vliq_bwater_max + (n_fwater-1)*(vliq_soil-vliq_bwater_max)
         k_soil = k_dsoil + k_bwater*vliq_bwater_max + k_fwater*(vliq_soil-vliq_bwater_max)
      ENDIF

      ! Computation of soil dielectric constant ([9], Eq.11).
      eps_soil_r = n_soil**2 - k_soil**2
      eps_soil_i = n_soil * k_soil * 2
      eps_soil = eps_soil_r - jj*eps_soil_i

   END SUBROUTINE cal_diel_soil_M04

!-----------------------------------------------------------------------------

   SUBROUTINE cal_diel_soil_M09(fghz, vliq_soil, t_soil, wf_clay, eps_soil)

!-----------------------------------------------------------------------------
      real(r8), intent(in)     :: fghz                              ! Microwave frequency [GHz]
      real(r8), intent(in)     :: vliq_soil                         ! Liquid volumetric soil water content [m3/m3]
      real(r8), intent(in)     :: t_soil                            ! Soil temperature [K]
      real(r8), intent(in)     :: wf_clay                           ! Clay mass fraction [0-1]
      complex(r8), intent(out) :: eps_soil                          ! Complex moist-soil relative permittivity [-]

      real(r8)                 :: n_dsoil                           ! Dry-soil refractive index [-]
      real(r8)                 :: k_dsoil                           ! Dry-soil normalized attenuation coefficient [-]
      real(r8)                 :: vliq_bwater_max                   ! Maximum bound-water volume fraction [m3/m3]
      real(r8)                 :: ts                                ! Reference temperature for parameter fits [degC]
      real(r8)                 :: eps0_bwater_ts                    ! Static bound-water relative permittivity at ts [-]
      real(r8)                 :: beta_bwater                       ! Bound-water permittivity temperature factor [1/K]
      real(r8)                 :: eps0_bwater                       ! Static bound-water relative permittivity at t_soil [-]
      real(r8)                 :: dHR_bwater                        ! Bound-water activation enthalpy divided by R [K]
      real(r8)                 :: dSR_bwater                        ! Bound-water activation entropy divided by R [-]
      real(r8)                 :: tau_bwater                        ! Bound-water relaxation time at t_soil [s]
      real(r8)                 :: sigma_bwater_ts                   ! Bound-water conductivity at ts [S/m]
      real(r8)                 :: beta_sigma_bwater                 ! Bound-water conductivity temperature slope [S/m/K]
      real(r8)                 :: sigma_bwater                      ! Bound-water conductivity at t_soil [S/m]
      real(r8)                 :: eps0_fwater_ts                    ! Static free-water relative permittivity at ts [-]
      real(r8)                 :: beta_fwater                       ! Free-water permittivity temperature factor [1/K]
      real(r8)                 :: eps0_fwater                       ! Static free-water relative permittivity at t_soil [-]
      real(r8)                 :: dHR_fwater                        ! Free-water activation enthalpy divided by R [K]
      real(r8)                 :: dSR_fwater                        ! Free-water activation entropy divided by R [-]
      real(r8)                 :: tau_fwater                        ! Free-water relaxation time at t_soil [s]
      real(r8)                 :: sigma_fwater_ts                   ! Free-water conductivity at ts [S/m]
      real(r8)                 :: beta_sigma_fwater                 ! Free-water conductivity temperature slope [S/m/K]
      real(r8)                 :: sigma_fwater                      ! Free-water conductivity at t_soil [S/m]
      real(r8)                 :: eps_bwater_r                      ! Real part of bound-water relative permittivity [-]
      real(r8)                 :: eps_bwater_i                      ! Bound-water dielectric loss factor [-]
      real(r8)                 :: eps_fwater_r                      ! Real part of free-water relative permittivity [-]
      real(r8)                 :: eps_fwater_i                      ! Free-water dielectric loss factor [-]
      real(r8)                 :: n_bwater                          ! Bound-water refractive index [-]
      real(r8)                 :: k_bwater                          ! Bound-water normalized attenuation coefficient [-]
      real(r8)                 :: n_fwater                          ! Free-water refractive index [-]
      real(r8)                 :: k_fwater                          ! Free-water normalized attenuation coefficient [-]
      real(r8)                 :: n_soil                            ! Moist-soil refractive index [-]
      real(r8)                 :: k_soil                            ! Moist-soil normalized attenuation coefficient [-]
      real(r8)                 :: eps_soil_r                        ! Real part of moist-soil relative permittivity [-]
      real(r8)                 :: eps_soil_i                        ! Moist-soil dielectric loss factor [-]
      real(r8)                 :: f_hz                              ! Microwave frequency [Hz]
      real(r8)                 :: t_soil_c                          ! Soil temperature for empirical fits [degC]
!-----------------------------------------------------------------------------
      ! Convert interface units for the empirical formulas.
      t_soil_c = t_soil - tfrz
      f_hz = fghz * 1e9

      ! Dry-soil refractive parameters ([10], Eqs.11-12).
      n_dsoil = 1.634 - 0.539*wf_clay + 0.2748*wf_clay**2
      k_dsoil = 0.03952 - 0.04038*wf_clay

      ! Maximum bound water fraction ([10], Eq.13).
      vliq_bwater_max = 0.02863 + 0.30673*wf_clay

      ! Reference temperature for parameter fits ([10], p.413).
      ts = 20.

      ! bound water parameters ([10], Eqs.14-15, 23, 7-8, 9, 18-19, 10, 22).
      eps0_bwater_ts = 79.8 - 85.4*wf_clay + 32.7*wf_clay**2
      beta_bwater = 8.67e-19 - 0.00126*wf_clay + 0.00184*wf_clay**2 - 0.000977*wf_clay**3 - 1.39e-7*wf_clay**4
      eps0_bwater  = (1+2*exp(log((eps0_bwater_ts-1)/(eps0_bwater_ts+2))-beta_bwater*(t_soil_c-ts))) / &
         (1-exp(log((eps0_bwater_ts-1)/(eps0_bwater_ts+2))-beta_bwater*(t_soil_c-ts)))
      dHR_bwater = 1467 + 2697*wf_clay - 980*wf_clay**2 + 0.0001368*wf_clay**3 - 0.0000861*wf_clay**4
      dSR_bwater = 0.888 + 9.7*wf_clay - 4.262*wf_clay**2 + 6.79e-15*wf_clay**3 + 4.263e-14*wf_clay**4
      tau_bwater = 48e-12 * exp(dHR_bwater/t_soil-dSR_bwater) / t_soil
      sigma_bwater_ts = 0.3112 + 0.467*wf_clay
      beta_sigma_bwater = 0.0028 + 0.02094*wf_clay - 0.01229*wf_clay**2 - 5.03e-16*wf_clay**3 + 4.163e-16*wf_clay**4
      sigma_bwater = sigma_bwater_ts + beta_sigma_bwater*(t_soil_c-ts)

      ! Unbound (free) water parameters ([10], Eqs. 7-8, 16-17, 25, 9, 20-21, 10, 24).
      eps0_fwater_ts = 100.
      beta_fwater = 1.11e-4 - 0.00001603*wf_clay + 0.00001239*wf_clay**2 + 8.33e-7*wf_clay**3 - 0.000001007*wf_clay**4
      eps0_fwater  = (1+2*exp(log((eps0_fwater_ts-1)/(eps0_fwater_ts+2))-beta_fwater*(t_soil_c-ts))) / &
         (1-exp(log((eps0_fwater_ts-1)/(eps0_fwater_ts+2))-beta_fwater*(t_soil_c-ts)))
      dHR_fwater = 2231 - 143.1*wf_clay + 223.2*wf_clay**2 - 142.1*wf_clay**3 + 27.14*wf_clay**4
      dSR_fwater = 3.649 - 0.4894*wf_clay + 0.763*wf_clay**2 - 0.4859*wf_clay**3 + 0.0928*wf_clay**4
      tau_fwater = 48e-12 * exp(dHR_fwater/t_soil-dSR_fwater) / t_soil
      sigma_fwater_ts = 0.05 + 1.4*(1.0-(1.0-wf_clay)**4.664)
      beta_sigma_fwater = 0.00108 + 0.1413*wf_clay - 0.2555*wf_clay**2 + 0.2147*wf_clay**3 - 0.0711*wf_clay**4
      sigma_fwater  = sigma_fwater_ts + beta_sigma_fwater*(t_soil_c-ts)

      ! Bound/free-water permittivity ([10], Eq.6; [9], Eq.16).
      eps_bwater_r = epsinf_water + (eps0_bwater-epsinf_water)/(1+(2*pi*f_hz*tau_bwater)**2)
      eps_bwater_i = (eps0_bwater-epsinf_water)/(1+(2*pi*f_hz*tau_bwater)**2)*(2*pi*f_hz*tau_bwater) + &
         sigma_bwater/(2*pi*eps0*f_hz)
      eps_fwater_r = epsinf_water + (eps0_fwater-epsinf_water)/(1+(2*pi*f_hz*tau_fwater)**2)
      eps_fwater_i = (eps0_fwater-epsinf_water)/(1+(2*pi*f_hz*tau_fwater)**2)*(2*pi*f_hz*tau_fwater) + &
         sigma_fwater/(2*pi*eps0*f_hz)

      ! Bound/free-water refractive indices ([10], Eq.5).
      n_bwater = sqrt(sqrt(eps_bwater_r**2+eps_bwater_i**2)+eps_bwater_r) / sqrt(2.0)
      k_bwater = sqrt(sqrt(eps_bwater_r**2+eps_bwater_i**2)-eps_bwater_r) / sqrt(2.0)
      n_fwater = sqrt(sqrt(eps_fwater_r**2+eps_fwater_i**2)+eps_fwater_r) / sqrt(2.0)
      k_fwater = sqrt(sqrt(eps_fwater_r**2+eps_fwater_i**2)-eps_fwater_r) / sqrt(2.0)

      ! Computation of soil refractive index (n_soil & k_soil) ([9], Eqs.12-13).
      IF (vliq_soil <= vliq_bwater_max) THEN
         n_soil = n_dsoil + (n_bwater-1)*vliq_soil
         k_soil = k_dsoil + k_bwater*vliq_soil
      ELSE
         n_soil = n_dsoil + (n_bwater-1)*vliq_bwater_max + (n_fwater-1)*(vliq_soil-vliq_bwater_max)
         k_soil = k_dsoil + k_bwater*vliq_bwater_max + k_fwater*(vliq_soil-vliq_bwater_max)
      ENDIF

      ! Computation of soil dielectric constant ([10], Eq.4).
      eps_soil_r = n_soil**2 - k_soil**2
      eps_soil_i = 2 * n_soil * k_soil
      eps_soil = eps_soil_r - jj*eps_soil_i

   END SUBROUTINE cal_diel_soil_M09

!-----------------------------------------------------------------------------

   SUBROUTINE cal_smooth_refl(theta, eps_soil, r_soil_smooth)

!-----------------------------------------------------------------------------
      real(r8), intent(in)     :: theta                             ! Incidence angle [rad]
      complex(r8), intent(in)  :: eps_soil                          ! Dielectric constant of the soil
      real(r8), intent(out)    :: r_soil_smooth(2)                  ! Reflectivities of flat surfaces for H and V polarizations

      complex(r8)              :: g                                 ! Parameter in Fresnel Law
!-----------------------------------------------------------------------------
      g = sqrt(eps_soil-sin(theta)**2)
      r_soil_smooth(1) = abs((cos(theta)-g)/(cos(theta)+g))**2.
      r_soil_smooth(2) = abs((cos(theta)*eps_soil-g)/(cos(theta)*eps_soil+g))**2.

   END SUBROUTINE cal_smooth_refl

!-----------------------------------------------------------------------------

   SUBROUTINE cal_rough_refl(fghz, theta, is_desert, patchclass, eps_soil, r_soil_smooth, r_soil_rough)

!-----------------------------------------------------------------------------
      real(r8), intent(in)     :: fghz                              ! Microwave frequency [GHz]
      real(r8), intent(in)     :: theta                             ! Incidence angle [rad]
      logical, intent(in)      :: is_desert                         ! Flag for desert soil
      integer, intent(in)      :: patchclass                        ! Patch class
      complex(r8), intent(in)  :: eps_soil                          ! Complex soil permittivity
      real(r8), intent(in)     :: r_soil_smooth(2)                  ! Reflectivities of flat surfaces for H and V polarizations
      real(r8), intent(out)    :: r_soil_rough(2)                   ! Reflectivities of rough surfaces for H and V polarizations

      real(r8)                 :: Q                                 ! Parameter for polarization mixing
      real(r8)                 :: hr(N_land_classification)         ! Roughness parameter
      real(r8)                 :: nrh(N_land_classification)        ! Parameter for H polarization
      real(r8)                 :: nrv(N_land_classification)        ! Parameter for V polarization
      real(r8)                 :: kcm                               ! Wave number [rad/cm]
      real(r8)                 :: f0 = 0.7                          ! Volume fraction of spherical particles
      real(r8)                 :: w                                 ! Single-particle albedo
      real(r8)                 :: g                                 ! Asymmetry parameter
      real(r8)                 :: a                                 ! Similarity parameter
      real(r8)                 :: y_r                               ! Real part of y-parameters
      real(r8)                 :: y_i                               ! Imaginary part of y-parameters
      real(r8)                 :: kr                                ! Size parameter for single-particle scattering albedo
!-----------------------------------------------------------------------------
      IF (is_desert) THEN
         ! Size parameters ([17], p.374)
         kr = (2*pi/(C/(fghz*1e9))) * (0.5*1e-3)

         ! Complex polarizability factors ([17], Eq.A15).
         y_r = (real(eps_soil)-1) / (real(eps_soil)+2)
         y_i = 3 * abs(aimag(eps_soil)) / (real(eps_soil)+2)**2

         ! Single-particle scattering albedo ([17], Eq. A16).
         w = (1-f0)**4 * kr**3 * y_r**2 / &
            ((1-f0)**4*kr**3*y_r**2+1.5*(1+2*f0)**2*y_i)

         ! Asymmetry parameter ([17], p.374).
         g = 0.23 * kr**2

         ! Similarity parameter ([17], Eq.3b).
         a = sqrt((1-w)/(1-w*g))

         ! Desert emissivity ([17], Eq.13).
         r_soil_rough = 1 - ((1-r_soil_smooth)*(2*a/((1+a)-(1-a)*r_soil_smooth)))
      ELSE
         kcm = (2*pi/(C/(fghz*1e9))) / 100.0

         ! Polarization mixing due to surface roughness ([16], Eq.16).
         IF (fghz < 2.) THEN
            Q = 0.
         ELSE
            Q = 0.35 * (1.0-exp(-0.6*rgh_surf**2*fghz))
         ENDIF

         ! Roughness parameters ([14], Table 2).
         IF (DEF_DA_OPERATOR_EM_ROUGH == 0) THEN
            hr(:) = (2.0*kcm*rgh_surf)**2.0
            nrh(:) = 0.0
            nrv(:) = 0.0
         ELSEIF (DEF_DA_OPERATOR_EM_ROUGH == 1) THEN
            hr(:) = hr_SMOS
            nrh(:) = 2.0
            nrv(:) = 0.0
         ELSEIF (DEF_DA_OPERATOR_EM_ROUGH == 2) THEN
            hr(:) = hr_SMAP
            nrh(:) = 2.0
            nrv(:) = 2.0
         ELSEIF (DEF_DA_OPERATOR_EM_ROUGH == 3) THEN
            hr(:) = hr_P16
            nrh(:) = -1.0
            nrv(:) = -1.0
         ENDIF

         ! Rough surface reflectivity for H and V polarizations ([14], Eq.5).
         r_soil_rough(1) = (Q*r_soil_smooth(2)+(1.-Q)*r_soil_smooth(1)) * exp(-hr(patchclass)*cos(theta)**nrh(patchclass))
         r_soil_rough(2) = (Q*r_soil_smooth(1)+(1.-Q)*r_soil_smooth(2)) * exp(-hr(patchclass)*cos(theta)**nrv(patchclass))
      ENDIF

   END SUBROUTINE cal_rough_refl

!-----------------------------------------------------------------------------

   SUBROUTINE veg(theta, patchclass, lai, tleaf, tb_veg, gamma_veg)

!-----------------------------------------------------------------------------
      real(r8), intent(in)     :: theta                             ! Incidence angle [rad]
      integer, intent(in)      :: patchclass                        ! Land cover class
      real(r8), intent(in)     :: lai                               ! Leaf area index
      real(r8), intent(in)     :: tleaf                             ! Leaf temperature (K)
      real(r8), intent(out)    :: tb_veg(2)                         ! Brightness temperature of vegetation for H- and V- polarization
      real(r8), intent(out)    :: gamma_veg(2)                      ! Canopy transmissivity for H and V polarizations

      real(r8)                 :: tau_veg_nadir                     ! Vegetation opacity at nadir
      real(r8)                 :: tau_veg(2)                        ! Vegetation opacity for H- and V- polarization
      integer                  :: ipol
!-----------------------------------------------------------------------------
      ! Calculate vegetation opacity (optical depth) at nadir b*VWC ([12], Eq.22)
      tau_veg_nadir = b1(patchclass)*lai + b2(patchclass)

      ! Calculate vegetation optical depth at H- and V- polarizations ([12], Eqs.23-24).
      tau_veg(1) = tau_veg_nadir * (cos(theta)**2+tth(patchclass)*sin(theta)**2)
      tau_veg(2) = tau_veg_nadir * (cos(theta)**2+ttv(patchclass)*sin(theta)**2)

      ! Calculate brightness temperature of vegetation ([12], Eq. 15).
      DO ipol = 1, 2
         gamma_veg(ipol) = exp(-tau_veg(ipol)/cos(theta))
         tb_veg(ipol) = (1.-w_CMEM(patchclass)) * (1.-gamma_veg(ipol)) * tleaf
      ENDDO

   END SUBROUTINE veg

!-----------------------------------------------------------------------------

   SUBROUTINE snow(fghz, theta, t_snow, t_below, snowdp, rho_snow, vliq_snow, r_below, r_snow, tb_tos)

!-----------------------------------------------------------------------------
      real(r8), intent(in)     :: fghz                              ! Microwave frequency [GHz]
      real(r8), intent(in)     :: theta                             ! Incidence angle in air from the normal [rad]
      real(r8), intent(in)     :: t_snow                            ! Uniform snow-layer temperature [K]
      real(r8), intent(in)     :: t_below                           ! Underlying surface temperature (soil or canopy) [K]
      real(r8), intent(in)     :: snowdp                            ! Geometric snow-layer thickness [m]
      real(r8), intent(in)     :: rho_snow                          ! Bulk snow density including liquid water [kg/m3]
      real(r8), intent(in)     :: vliq_snow                         ! Liquid water volume per bulk snow volume [m3/m3]
      real(r8), intent(in)     :: r_below(2)                        ! Power reflectivity at the lower snow boundary [-]
      real(r8), intent(out)    :: r_snow(2)                         ! Effective reflectivity of snow and underlying surface [-]
      real(r8), intent(out)    :: tb_tos(2)                         ! Emitted brightness temperature above the snow [K]

      real(r8)                 :: sal_snow = 0.0                    ! Snow salinity for the ice-loss correction [ppm]
      real(r8)                 :: rho_snow_gcm3                     ! Bulk snow density for empirical formulas [g/cm3]
      real(r8)                 :: rho_dsnow_gcm3                    ! Density of the dry air/ice matrix [g/cm3]
      real(r8)                 :: grain_snow_mm                     ! Empirical snow grain diameter [mm]

      complex(r8)              :: eps_ice                           ! Complex relative permittivity of pure ice [-]
      real(r8)                 :: eps_ice_r                         ! Real part of ice relative permittivity [-]
      real(r8)                 :: eps_ice_i                         ! Ice loss factor after salinity correction [-]
      real(r8)                 :: eps_ice_i_impure                  ! Reference impure-ice loss factor at -5 degC [-]
      real(r8)                 :: eps_ice_i_pure                    ! Reference pure-ice loss factor at -5 degC [-]
      real(r8)                 :: eps_dsnow_r                       ! Real part of dry-snow relative permittivity [-]
      real(r8)                 :: eps_dsnow_i                       ! Dry-snow dielectric loss factor [-]
      real(r8)                 :: eps_snow_r                        ! Real part of effective snow relative permittivity [-]
      real(r8)                 :: eps_snow_i                        ! Effective snow dielectric loss factor [-]
      complex(r8)              :: eps_wsnow                         ! Complex relative permittivity of wet snow [-]

      real(r8)                 :: eps0_water = 88.                  ! Static liquid-water relative permittivity [-]
      real(r8)                 :: f0_water_a                        ! a-axis water-inclusion relaxation frequency [GHz]
      real(r8)                 :: f0_water_b                        ! b-axis water-inclusion relaxation frequency [GHz]
      real(r8)                 :: f0_water_c                        ! c-axis water-inclusion relaxation frequency [GHz]
      real(r8)                 :: epsinf_water_a                    ! High-frequency a-axis dielectric contribution [-]
      real(r8)                 :: epsinf_water_b                    ! High-frequency b-axis dielectric contribution [-]
      real(r8)                 :: epsinf_water_c                    ! High-frequency c-axis dielectric contribution [-]
      real(r8)                 :: eps0_water_a                      ! Static a-axis dielectric contribution [-]
      real(r8)                 :: eps0_water_b                      ! Static b-axis dielectric contribution [-]
      real(r8)                 :: eps0_water_c                      ! Static c-axis dielectric contribution [-]
      complex(r8)              :: eps_water_a                       ! Complex a-axis dielectric contribution [-]
      complex(r8)              :: eps_water_b                       ! Complex b-axis dielectric contribution [-]
      complex(r8)              :: eps_water_c                       ! Complex c-axis dielectric contribution [-]

      real(r8)                 :: omega                             ! Angular frequency [rad/s]
      real(r8)                 :: k0                                ! Free-space wave number [rad/m]
      real(r8)                 :: alpha_snow                        ! Amplitude attenuation constant [1/m]
      real(r8)                 :: beta_snow                         ! Phase constant [rad/m]
      real(r8)                 :: D_I                               ! Loss term of squared normal wave number [1/m2]
      real(r8)                 :: D_R                               ! Real part of squared normal wave number [1/m2]
      real(r8)                 :: theta_snow                        ! Phase refraction angle from the snow normal [rad]
      complex(r8)              :: z_snow                            ! Intrinsic impedance of snow [ohm]

      real(r8)                 :: b_dsnow                           ! Squared dry-snow dielectric loss tangent [-]
      real(r8)                 :: b_wsnow                           ! Squared effective snow dielectric loss tangent [-]
      real(r8)                 :: ka_dsnow                          ! Dry-snow power absorption coefficient [1/m]
      real(r8)                 :: ka_snow                           ! Effective snow power absorption coefficient [1/m]
      real(r8)                 :: ke_dsnow                          ! Dry-snow power extinction coefficient [1/m]
      real(r8)                 :: ke_snow                           ! Effective snow power extinction coefficient [1/m]
      real(r8)                 :: ks_snow                           ! Snow power scattering coefficient [1/m]
      real(r8)                 :: gamma_snow                        ! One-way effective snow-layer transmissivity [-]

      real(r8)                 :: r_snow_air                        ! Power reflectivity at the upper snow boundary [-]
      real(r8)                 :: tb_snow                           ! Snow-layer emission reaching the upper air side [K]
      real(r8)                 :: tb_below                          ! Underlying emission reaching the upper air side [K]
      integer                  :: ipol                              ! Polarization index (1=H, 2=V)
!-----------------------------------------------------------------------------
      omega = 2.0 * pi * (fghz*1e9)
      k0 = 2 * pi / (C/(fghz*1e9))

      ! Convert interface units for the empirical formulas.
      rho_snow_gcm3 = rho_snow / 1000.0_r8

      ! Pure-ice relative permittivity and positive dielectric loss.
      CALL cal_diel_ice(fghz, t_snow, eps_ice)
      eps_ice_r = real(eps_ice)
      eps_ice_i = abs(aimag(eps_ice))

      ! Salinity correction using reference ice losses at -5 degC.
      eps_ice_i_impure = 0.0026/fghz + 0.00023*(fghz**0.87)         
      eps_ice_i_pure = 6.e-4/fghz + 6.5e-5*(fghz**1.07)                
      eps_ice_i = eps_ice_i + (eps_ice_i_impure-eps_ice_i_pure)*sal_snow/13.0d0    

      ! Dry-snow effective permittivity ([18], Eq.2.19; [23], Eq.21).
      rho_dsnow_gcm3 = denice/1000.0_r8
      IF (vliq_snow < 1.0_r8) THEN
         rho_dsnow_gcm3 = max(0.0_r8, min(denice/1000.0_r8, (rho_snow_gcm3-vliq_snow)/(1.0_r8-vliq_snow)))
      ENDIF
      eps_dsnow_r = 1.0 + 1.58*rho_dsnow_gcm3/(1.0-0.365*rho_dsnow_gcm3)
      eps_dsnow_i = 3.0 * (rho_dsnow_gcm3/(denice/1000.0_r8)) * eps_ice_i * (eps_dsnow_r**2) * &
         (2*eps_dsnow_r+1) / ((eps_ice_r+2*eps_dsnow_r)*(eps_ice_r+2*eps_dsnow_r**2))
 
      ! Wet-snow dielectric mixing ([18], Eqs.2.24-2.26)
      IF (vliq_snow > 0.) THEN                                   
         ! Relaxation frequencies of the three water-inclusion contributions.
         ! Principal-axis depolarization factors: 0.005, 0.4975, 0.4975.
         f0_water_a = f0_water * (1+(0.005*(eps0_water-epsinf_water)/(eps_dsnow_r+(0.005*(epsinf_water-eps_dsnow_r)))))
         f0_water_b = f0_water * (1+(0.4975*(eps0_water-epsinf_water)/(eps_dsnow_r+(0.4975*(epsinf_water-eps_dsnow_r)))))
         f0_water_c = f0_water * (1+(0.4975*(eps0_water-epsinf_water)/(eps_dsnow_r+(0.4975*(epsinf_water-eps_dsnow_r)))))

         ! High-frequency dielectric contributions along the three principal axes.
         epsinf_water_a = (vliq_snow*(epsinf_water-eps_dsnow_r)/3.) / (1.+0.005*((epsinf_water/eps_dsnow_r)-1.))
         epsinf_water_b = (vliq_snow*(epsinf_water-eps_dsnow_r)/3.) / (1.+0.4975*((epsinf_water/eps_dsnow_r)-1.))
         epsinf_water_c = (vliq_snow*(epsinf_water-eps_dsnow_r)/3.) / (1.+0.4975*((epsinf_water/eps_dsnow_r)-1.))

         ! Static dielectric contributions along the three principal axes.
         eps0_water_a = (vliq_snow/3.) * (eps0_water-eps_dsnow_r) / (1.+0.005*((eps0_water/eps_dsnow_r)-1.))
         eps0_water_b = (vliq_snow/3.) * (eps0_water-eps_dsnow_r) / (1.+0.4975*((eps0_water/eps_dsnow_r)-1.))
         eps0_water_c = (vliq_snow/3.) * (eps0_water-eps_dsnow_r) / (1.+0.4975*((eps0_water/eps_dsnow_r)-1.))

         ! Principal-axis Debye contributions with negative imaginary parts.
         eps_water_a = epsinf_water_a + (eps0_water_a-epsinf_water_a)/(1+jj*fghz/f0_water_a)
         eps_water_b = epsinf_water_b + (eps0_water_b-epsinf_water_b)/(1+jj*fghz/f0_water_b)
         eps_water_c = epsinf_water_c + (eps0_water_c-epsinf_water_c)/(1+jj*fghz/f0_water_c)

         ! Add the water-inclusion contributions to the dry air/ice matrix.
         eps_wsnow = eps_water_a + eps_water_b + eps_water_c + (eps_dsnow_r-jj*eps_dsnow_i)
         eps_snow_r = real(eps_wsnow)
         eps_snow_i = abs(aimag(eps_wsnow))
      ELSE
         eps_snow_r = eps_dsnow_r
         eps_snow_i = eps_dsnow_i
      ENDIF

      ! Attenuation and phase constants ([24], Eqs.2.6.12-13).
      alpha_snow = k0 * abs(aimag(sqrt(eps_snow_r-jj*eps_snow_i)))
      beta_snow = k0 * real(sqrt(eps_snow_r-jj*eps_snow_i))

      ! Real and loss terms of the squared normal wave number ([24], Eq.7.9.4).
      D_I = 2. * alpha_snow * beta_snow
      D_R = beta_snow**2 - alpha_snow**2 - (k0*k0)*(sin(theta)**2)

      ! Phase refraction angle ([24], Eqs.7.9.5-6).
      theta_snow = atan(k0*sin(theta)/((1./sqrt(2.)) *sqrt(sqrt(D_I**2.+D_R**2.)+D_R)))

      ! Intrinsic impedance of nonmagnetic snow ([24], Eq.2.6.5).
      z_snow = z0 / sqrt(eps_snow_r-jj*eps_snow_i)

      ! Calculate H and V brightness temperatures above the snow.
      DO ipol = 1, 2
         ! Power reflectivity at the upper snow boundary.
         IF (ipol == 1) THEN
            r_snow_air = abs((z_snow*cos(theta)-z0*cos(theta_snow))/ (z_snow*cos(theta)+z0*cos(theta_snow)))**2
         ELSE
            r_snow_air = abs((z0*cos(theta)-z_snow*cos(theta_snow))/ (z0*cos(theta)+z_snow*cos(theta_snow)))**2
         ENDIF

         ! Snow grain size (mm)([19], Eq.5.1).
         grain_snow_mm = min(1000*(1.6e-4+1.1e-13*(rho_snow**4)),3.0)

         ! Extinction coefficient of dry snow ([23], Eq.14)
         ! NOTE: L-band use is an extrapolation of higher-frequency fit.
         ke_dsnow = 0.0018 * (fghz**2.8) * (grain_snow_mm**2) / 4.3429

         ! Dry-snow power absorption coefficient.
         b_dsnow = (eps_dsnow_i/eps_dsnow_r)**2
         ka_dsnow = 2. * omega * sqrt(mu0*eps0*eps_dsnow_r) * sqrt(b_dsnow/(2.*(sqrt(1.+b_dsnow)+1.)))
         IF (ke_dsnow < ka_dsnow) ke_dsnow = ka_dsnow

         ! Effective snow power absorption coefficient.
         b_wsnow = (eps_snow_i/eps_snow_r)**2
         ka_snow = 2. * omega * sqrt(mu0*eps0*eps_snow_r) * sqrt(b_ws/(2.*(sqrt(1.+b_wsnow)+1)))

         ! Total extinction (assuming scattering is the same for dry and wet snow)
         ke_snow = (ke_dsnow-ka_dsnow) + ka_snow

         ! Scattering coefficient shared by dry and wet snow.
         ks_snow = ke_dsnow - ka_dsnow

         ! Effective snow-layer transmissivity ([22], Eq.4).
         gamma_snow = exp(-(ke_snow-0.96*ks_snow)*snowdp/cos(theta_snow))

         ! Snow-layer emission transmitted through the upper interface ([20], Eq.4.161).
         tb_snow = (1.0_r8+r_below(ipol)*gamma_snow) * (1.0_r8-r_snow_air) * t_snow * &
            (ka_snow/(ke_snow-q*ks_snow)) * (1.0_r8-gamma_snow) / (1.0_r8-r_below(ipol)*r_snow_air*gamma_snow**2)

         ! Underlying-layer emission transmitted through the snow ([20], Eq.4.162).
         tb_below = (1.0_r8-r_below(ipol)) * (1.0_r8-r_snow_air) * t_below * &
            gamma_snow / (1.0_r8-r_below(ipol)*r_snow_air*gamma_snow**2)

         ! Sum snow and underlying-surface emission above the snow.
         tb_tos(ipol) = tb_snow + tb_below

         ! Effective reflectivity including both snow-layer boundaries.
         r_snow(ipol) = 1 - (tb_snow/t_snow+tb_below/t_below)
      ENDDO

   END SUBROUTINE snow
!-----------------------------------------------------------------------------
END MODULE MOD_DA_Operator_EM
#endif
