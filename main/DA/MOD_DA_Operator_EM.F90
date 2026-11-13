#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Operator_EM
!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Microwave emission model adapted from the Community Microwave Emission
!   Modelling Platform (CMEM) v6.1. Simulates horizontally and vertically
!   polarized top-of-atmosphere brightness temperatures from land surface
!   states, microwave frequency, and incidence angle.
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
! HISTORY:
!   Lu Li, 12/2024: Initial version
!   Zhilong Fan, Lu Li, 03/2024: Debug and clean codes
!   Lu Li, 10/2025: Debug and clean codes
!   Lu Li, Zhilong Fan, 11/2026: DA v1.0
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Const_Physical
   USE MOD_Vars_1DForcing
   USE MOD_SPMD_Task
   USE MOD_Vars_Global, only: nl_soil, nl_lake, N_land_classification, maxsnl, spval, dz_soi
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
      real(r8), intent(in)     :: htop                              ! Upper height of vegetation [m]
      real(r8), intent(in)     :: tref                              ! 2 m height air temperature [kelvin]
      real(r8), intent(in)     :: tleaf                             ! Leaf temperature [K]
      real(r8), intent(in)     :: t_soisno(maxsnl + 1:nl_soil)      ! Soil temperature [K]
      real(r8), intent(in)     :: wliq_soisno(maxsnl + 1:nl_soil)   ! Liquid water in layers [kg/m2]
      real(r8), intent(in)     :: wice_soisno(maxsnl + 1:nl_soil)   ! Ice lens in layers [kg/m2]
      real(r8), intent(in)     :: snowdp                            ! Snow depth [meter]
      real(r8), intent(in)     :: lai                               ! Leaf area index
      real(r8), intent(in)     :: sai                               ! Stem area index
      real(r8), intent(in)     :: wf_clay(nl_soil)                  ! Gravimetric fraction of clay
      real(r8), intent(in)     :: wf_sand(nl_soil)                  ! Gravimetric fraction of sand
      real(r8), intent(in)     :: wf_silt(nl_soil)                  ! Gravimetric fraction of silt
      real(r8), intent(in)     :: BD_all(nl_soil)                   ! Bulk density of soil (GRAVELS + ORGANIC MATTER + Mineral Soils,kg/m3)
      real(r8), intent(in)     :: porsl(nl_soil)                    ! Fraction of soil that is voids [-]
      real(r8), intent(in)     :: sat_theta                         ! Incidence angle of satellite (rad)
      real(r8), intent(in)     :: sat_fghz                          ! Frequency of satellite (GHz)
      real(r8), intent(out)    :: tb_toa_h                          ! Brightness temperature of top-of-atmosphere for H- polarization
      real(r8), intent(out)    :: tb_toa_v                          ! Brightness temperature of top-of-atmosphere for V- polarization

      logical                  :: is_low_veg                        ! Flag for low vegetation
      real(r8)                 :: dz_soisno(maxsnl+1:nl_soil)       ! Snow and soil layer thickness [m]
      integer                  :: lb                                ! Lower bound of arrays
      real(r8)                 :: tau_atm                           ! Atmospheric optical depth
      real(r8)                 :: r_r(2)                            ! Rough surface reflectivity for H and V polarizations
      real(r8)                 :: r_sn(2)                           ! Reflectivity between the snow and ground for H and V polarizations
      real(r8)                 :: r_snow(2)                         ! Reflectivity of the snow for H and V polarizations
      real(r8)                 :: tb_soil(2)                        ! Brightness temperature of soil for H and V polarizations
      real(r8)                 :: tb_tos(2)                         ! Brightness temperature of snow-covered ground for H and V polarizations
      real(r8)                 :: tb_tov(2)                         ! Brightness temperature of vegetation (consider snow) for H and V polarizations
      real(r8)                 :: tb_tov_noad(2)                    ! Brightness temperature of vegetation (no downwelling radiation) for H and V polarizations
      real(r8)                 :: tb_au(2)                          ! Upwelling radiation (brightness temperature) of atmosphere
      real(r8)                 :: tb_ad(2)                          ! Downwelling radiation (brightness temperature) of atmosphere
      real(r8)                 :: rho_snow                          ! Snow density (g/cm3)
      real(r8)                 :: liq_snow                          ! Snow liquid water content (cm3/cm3)
      real(r8)                 :: gamma_p(2)                        ! Canopy transmissivity for H and V polarizations
      real(r8)                 :: tb_veg(2)                         ! Brightness temperature of vegetation for H- and V- polarization
      real(r8)                 :: tb_2(2)                           ! The downwelling vegetation emission reflected by the soil and attenuated by the canopy layer
      real(r8)                 :: tb_3(2)                           ! Upwelling soil emission attenuated by the canopy
      real(r8)                 :: tb_4(2)                           ! The downwelling cosmic ray reflected by the soil and attenuated by the canopy layer
      real(r8)                 :: tb_toa(2)                         ! Brightness temperature of top-of-atmosphere for H- and V- polarization
      real(r8)                 :: wf_total(nl_soil)                 ! Total gravimetric
      real(r8)                 :: BD_all_surf                       ! Surface soil bulk density [g/cm3]
      real(r8)                 :: porsl_surf                        ! Soil porosity at surface
      real(r8)                 :: t_surf                            ! Soil temperature at surface (C)
      real(r8)                 :: t_deep                            ! Soil temperature at deep layer (C)
      real(r8)                 :: liq_surf                          ! Liquid volumetric water content at surface (m3/m3)
      real(r8)                 :: ice_surf                          ! Ice volumetric water content at surface (m3/m3)
      real(r8)                 :: wf_clay_surf                      ! Gravimetric clay percent fraction(%) at surface
      real(r8)                 :: wf_sand_surf                      ! Gravimetric sand percent fraction(%) at surface
      integer                  :: i
!-----------------------------------------------------------------------------
!#############################################################################
! Prepare parameters & states used in the operator.
!#############################################################################
      ! Get depth of soil and snow layers
      dz_soisno(maxsnl+1:0) = dz_sno(maxsnl+1:0)
      dz_soisno(1:nl_soil) = dz_soi(1:nl_soil)

      ! Calculate weighted parameters
      wf_total     = wf_clay + wf_sand + wf_silt
      wf_clay_surf = (wf_clay(1)/wf_total(1)*0.0175 + wf_clay(2)/wf_total(2)*0.0276)/0.0451*100
      wf_sand_surf = (wf_sand(1)/wf_total(1)*0.0175 + wf_sand(2)/wf_total(2)*0.0276)/0.0451*100
      BD_all_surf  = (BD_all(1)*0.0175 + BD_all(2)*0.0276)/0.0451/1000
      porsl_surf   = (porsl(1)*0.0175 + porsl(2)*0.0276)/0.0451

      ! Caculate temperature (C) at surface and deep soil layers
      t_surf = ((t_soisno(1)*(0.0175) + t_soisno(2)*(0.0451 - 0.0175))/0.0451) - tfrz
      t_deep = ((t_soisno(7)*(0.8289-0.5) + t_soisno(8)*(1.0 - 0.8289))/0.5) - tfrz

      ! Calculate liquid/ice volumetric water (first two layers)
      liq_surf = (wliq_soisno(1) + wliq_soisno(2))/(0.0451*denh2o)
      ice_surf = (wice_soisno(1) + wice_soisno(2))/(0.0451*denice)

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
         ! Soil module.
!-----------------------------------------------------------------------------
         CALL soil(sat_fghz, sat_theta,                          &
            patchclass,                                          &
            t_surf, t_deep,                                      &
            liq_surf, ice_surf,                                  &
            wf_sand_surf, wf_clay_surf, BD_all_surf, porsl_surf, &
            r_r, tb_soil)
!-----------------------------------------------------------------------------
         ! Vegetation and snow module.
!-----------------------------------------------------------------------------
         ! Roughly judge low or high vegetation (only for IGBP)
         is_low_veg = .true.
         IF (patchclass >= 1 .and. patchclass <= 5) THEN
            is_low_veg = .false.
         ENDIF

         ! Ensure snow density to <= 1 g/cm3
         IF (snowdp > 0.01) THEN
            rho_snow = (wliq_soisno(lb) + wice_soisno(lb))/(dz_soisno(lb)*1e3)
            liq_snow = wliq_soisno(lb)*rho_snow/(dz_soisno(lb)*1e3)

            IF (liq_snow > 1.0 .or. rho_snow > 1.0) THEN
               rho_snow = 1.0
               liq_snow = wliq_soisno(lb)*rho_snow/(dz_soisno(lb)*1e3)
            ENDIF
         ENDIF

         ! 1) no veg and no snow
         IF ((lai + sai < 1e-6) .and. (snowdp < 0.01)) THEN
            tb_tov = tb_soil + tb_ad*r_r
            tb_tov_noad = tb_soil

         ! 2) no veg and has snow
         ELSEIF ((lai + sai < 1e-6) .and. (snowdp > 0.01)) THEN
            CALL snow(sat_fghz, sat_theta, t_soisno(1), t_soisno(1), snowdp, rho_snow, liq_snow, r_r, r_snow, tb_tos)

            tb_tov = tb_tos + tb_ad*r_snow
            tb_tov_noad = tb_tos

         ! 3) has veg and no snow
         ELSEIF ((lai + sai > 1e-6) .and. (snowdp < 0.01)) THEN
            CALL veg(sat_theta, patchclass, lai, htop, 0.0, tleaf, tb_veg, gamma_p)

            DO i = 1, 2
               tb_2(i) = tb_veg(i)*gamma_p(i)*r_r(i)
               tb_3(i) = tb_soil(i)*gamma_p(i)
               tb_4(i) = tb_ad(i)*r_r(i)*(gamma_p(i)**2)
               tb_tov(i) = tb_veg(i) + tb_2(i) + tb_3(i) + tb_4(i)
               tb_tov_noad(i) = tb_veg(i) + tb_2(i) + tb_3(i)
            ENDDO

         ! 4) has veg and has snow
         ELSEIF ((lai + sai > 1e-6) .and. (snowdp > 0.01)) THEN
            IF (htop < snowdp) THEN
               CALL veg(sat_theta, patchclass, lai, htop, snowdp, tleaf, tb_veg, gamma_p)

               DO i = 1, 2
                  tb_2(i) = tb_veg(i)*gamma_p(i)*r_r(i)
                  tb_3(i) = tb_soil(i)*gamma_p(i)
                  tb_4(i) = tb_ad(i)*r_r(i)*(gamma_p(i)**2)
                  tb_tov(i) = tb_veg(i) + tb_2(i) + tb_3(i) + tb_4(i)
                  tb_tov_noad(i) = tb_veg(i) + tb_2(i) + tb_3(i)
               ENDDO

               ! Calculate reflectivity between the snow and low veg (adopted from CMEM)
               r_sn(:) = 1.0 - tb_tov_noad(:)/t_soisno(1)

               ! Calculate brightness temperature of snow-covered ground
               CALL snow(sat_fghz, sat_theta, t_soisno(1), t_soisno(1), snowdp, rho_snow, liq_snow, r_sn, r_snow, tb_tos)

               ! Calculate brightness temperature of top of snow
               tb_tov = tb_tos + tb_ad*r_snow
               tb_tov_noad = tb_tos

            ELSE
               ! Calculate brightness temperature of snow-covered ground
               CALL snow(sat_fghz, sat_theta, t_soisno(1), t_soisno(1), snowdp, rho_snow, liq_snow, r_r, r_snow, tb_tos)

               ! Calculate brightness temperature of high vegetation
               CALL veg(sat_theta, patchclass, lai, htop, snowdp, tleaf, tb_veg, gamma_p)

               ! Calculate brightness temperature of top of vegetation
               DO i = 1, 2
                  tb_2(i) = tb_veg(i)*gamma_p(i)*r_snow(i)
                  tb_3(i) = tb_tos(i)*gamma_p(i)
                  tb_4(i) = tb_ad(i)*r_snow(i)*(gamma_p(i)**2)
                  tb_tov(i) = tb_veg(i) + tb_2(i) + tb_3(i) + tb_4(i)
                  tb_tov_noad(i) = tb_veg(i) + tb_2(i) + tb_3(i)
               ENDDO
            ENDIF
         ENDIF

         ! Atmosphere module.
         CALL atm(sat_theta, forc_topo, tref, tau_atm, tb_au, tb_ad)

         ! Calculate brightness temperature of top-of-atmosphere.
         tb_toa = tb_tov*exp(-tau_atm) + tb_au

      ENDIF

      tb_toa_h = tb_toa(1)
      tb_toa_v = tb_toa(2)

   END SUBROUTINE EM

!-----------------------------------------------------------------------------

   SUBROUTINE atm(theta, z, tref, tau_atm, tb_au, tb_ad)

!-----------------------------------------------------------------------------
      real(r8), intent(in)     :: theta                             ! Incidence angle [rad]
      real(r8), intent(in)     :: z                                 ! Altitude (m)
      real(r8), intent(in)     :: tref                              ! 2m air temperature (K)
      real(r8), intent(out)    :: tau_atm                           ! Atmospheric optical depth
      real(r8), intent(out)    :: tb_au(2)                          ! Upwelling radiation (brightness temperature) of atmosphere
      real(r8), intent(out)    :: tb_ad(2)                          ! Downwelling radiation (brightness temperature) of atmosphere

      real(r8)                 :: t_sky = 2.7                       ! Cosmic background brightness temperature [K]
      real(r8)                 :: t_eq                              ! Equivalent layer temperature
      real(r8)                 :: gossat
!-----------------------------------------------------------------------------
      ! Atmospheric optical depth ([11])
      tau_atm = exp(-3.9262 - 0.2211*z/1000 - 0.00369*tref)/cos(theta)
      gossat = exp(-tau_atm)

      ! Calculate equivalent layer temperature ([11])
      t_eq = exp(4.9274 + 0.002195*tref)

      ! Upwelling radiation (brightness temperature) of atmosphere ([11])
      tb_au(:) = t_eq*(1.-gossat)

      ! Downwelling atmospheric and attenuated cosmic emission ([11])
      tb_ad(:) = t_eq*(1.-gossat) + t_sky*gossat

   END SUBROUTINE atm

!-----------------------------------------------------------------------------

   SUBROUTINE soil(fghz, theta,                            &
      patchclass,                                          &
      t_surf, t_deep,                                      &
      liq_surf, ice_surf,                                  &
      wf_sand_surf, wf_clay_surf, BD_all_surf, porsl_surf, &
      r_r, tb_soil)

!-----------------------------------------------------------------------------
      real(r8), intent(in)     :: fghz                              ! Microwave frequency [GHz]
      real(r8), intent(in)     :: theta                             ! Incidence angle [rad]
      integer, intent(in)      :: patchclass                        ! Land cover class
      real(r8), intent(in)     :: t_surf                            ! Soil temperature at surface (C)
      real(r8), intent(in)     :: t_deep                            ! Soil temperature at deep layer (C)
      real(r8), intent(in)     :: liq_surf                          ! Liquid volumetric water content at surface (m3/m3)
      real(r8), intent(in)     :: ice_surf                          ! Ice volumetric water content at surface (m3/m3)
      real(r8), intent(in)     :: wf_sand_surf                      ! Gravimetric sand percent fraction(%) at surface
      real(r8), intent(in)     :: wf_clay_surf                      ! Gravimetric clay percent fraction(%) at surface
      real(r8), intent(in)     :: BD_all_surf                       ! Surface soil bulk density [g/cm3]
      real(r8), intent(in)     :: porsl_surf                        ! Soil porosity at surface
      real(r8), intent(out)    :: r_r(2)                            ! Rough surface reflectivity for H and V polarizations
      real(r8), intent(out)    :: tb_soil(2)                        ! Brightness temperature of soil

      real(r8)                 :: t_eff(2)                          ! Effective temperature for H and V polarizations, [K]
      complex(r8)              :: eps_soil                          ! Dielectric constant of soil for H and V polarizations
      real(r8)                 :: r_s(2)                            ! Smooth surface reflectivity for H and V polarizations
      complex(r8)              :: ew                                ! Dielectric constant of water
      logical                  :: is_desert                         ! Flag for desert soil
      real(r8)                 :: ffrz                              ! Fraction of frozen soil
      complex(r8)              :: eps_f = (5.0, 0.5)                ! Dielectric constant of frozen soil
      real(r8)                 :: sal_soil = 0.0                    ! Soil salinity (psu)
!-----------------------------------------------------------------------------
      ! Whether this patch is desert
      is_desert = .false.
      IF (liq_surf < 0.02 .and. wf_sand_surf > 90) THEN
         is_desert = .true.
      ENDIF

      ! Calculate ratio of freezed soil
      IF (liq_surf + ice_surf <= 0.0d0) THEN
         ffrz = 0.0d0
      ELSE
         ffrz = ice_surf / (liq_surf + ice_surf)
      ENDIF

      ! Calculate effective temperature
      CALL eff_tsoi(liq_surf, t_surf, t_deep, t_eff)

      ! Calculate dielectric constant of soil (mixture medium) ([13], Eq. 1).
      IF (is_desert) THEN
         ! Dry-sand permittivity.
         eps_soil = 2.53 + (2.79 - 2.53)/(1 - jj*(fghz/0.27)) + jj*0.002
      ELSE
         ! Caculate ice or water dielectric constant
         IF (ffrz > 0.95) THEN
            CALL diel_ice(fghz, t_surf, ew)
         ELSE
            CALL diel_water(fghz, 2, liq_surf, t_surf, wf_sand_surf, wf_clay_surf, BD_all_surf, sal_soil, ew)
         ENDIF

         ! Caculate dielectric constant in mixed soil
         IF (DEF_DA_OPERATOR_EM_DIEL == 0) THEN
            CALL diel_soil_W80(ew, t_surf, liq_surf, wf_sand_surf, wf_clay_surf, porsl_surf, eps_soil)
         ELSEIF (DEF_DA_OPERATOR_EM_DIEL == 1) THEN
            CALL diel_soil_D85(ew, liq_surf, wf_sand_surf, wf_clay_surf, BD_all_surf, eps_soil)
         ELSEIF (DEF_DA_OPERATOR_EM_DIEL == 2) THEN
            CALL diel_soil_M04(fghz, liq_surf, wf_clay_surf, eps_soil)
         ELSEIF (DEF_DA_OPERATOR_EM_DIEL == 3) THEN
            CALL diel_soil_M09(fghz, liq_surf, t_surf, wf_clay_surf, eps_soil)
         ENDIF
      END IF

      ! Mix dielectric constant of frozen and non-frozen soil
      eps_soil = eps_soil*(1.-ffrz) + eps_f*ffrz

      ! Calculate smooth surface reflectivity
      CALL smooth_reflectivity(theta, eps_soil, r_s)

      ! Calculate rough surface reflectivity
      CALL rough_reflectivity(fghz, theta, is_desert, patchclass, r_s, r_r)

      ! Calculate brightness temperature
      IF (is_desert) THEN
         CALL desert(fghz, t_eff, r_r, eps_soil, tb_soil)
      ELSE
         tb_soil = t_eff * (1 - r_r)
      ENDIF

   END SUBROUTINE soil

!-----------------------------------------------------------------------------

   SUBROUTINE eff_tsoi(wc_surf, t_surf, t_deep, t_eff)

!-----------------------------------------------------------------------------
      real(r8), intent(in)     :: wc_surf                           ! Soil moisture at surface (m3/m3)
      real(r8), intent(in)     :: t_surf                            ! Soil temperature (C) at surface
      real(r8), intent(in)     :: t_deep                            ! Soil temperature (C) at deep layer
      real(r8), intent(out)    :: t_eff(2)                          ! Effective temperature for H and V polarizations, [K]

      real(r8)                 :: C                                 ! Moisture-dependent surface-temperature weight
      real(r8)                 :: w0 = 0.30                         ! Reference soil moisture [m3/m3]
      real(r8)                 :: bw = 0.30                         ! Moisture-weighting exponent
!-----------------------------------------------------------------------------
      ! Calculate moisture-dependent weight ([21])
      IF (wc_surf < 0.0) THEN
         C = 0.001
      ELSE
         C = max(0.001, (wc_surf/w0)**bw)
      ENDIF

      ! Calculate effective soil temperature ([14], Eq. 11)
      t_eff(:) = t_deep + (t_surf - t_deep)*C + tfrz

   END SUBROUTINE eff_tsoi

!-----------------------------------------------------------------------------

   SUBROUTINE diel_ice(fghz, t, eps_i)

!-----------------------------------------------------------------------------
      real(r8), intent(in)     :: fghz                              ! Microwave frequency [GHz]
      real(r8), intent(in)     :: t                                 ! Temperature (C)
      complex(r8), intent(out) :: eps_i                             ! Complex dielectric permittivity of ice

      real(r8)                 :: betam                             ! Beta parameter by Mishima et al. (1983)
      real(r8)                 :: dbeta                             ! Corrected delta beta parameter
      real(r8)                 :: beta                              ! Beta parameter
      real(r8)                 :: t_inv                             ! Modified inverse temperature
      real(r8)                 :: tk                                ! Temperature (K)
      real(r8)                 :: alpha                             ! Alpha parameter
      real(r8)                 :: eps_i_r                           ! Real part of pure ice dielectric constant
      real(r8)                 :: eps_i_i                           ! Imaginary part of pure ice dielectric constant
!-----------------------------------------------------------------------------
      ! C to K
      tk = t + tfrz

      ! Calculate beta parameter by Mishima et al. (1983) ([3], Eq. 5.33)
      betam = (0.0207/tk)*(exp(335./tk)/((exp(335./tk) - 1.)**2.)) + 1.16e-11*(fghz**2.)

      ! Calculate delta beta parameter ([3], Eq. 5.35)
      dbeta = exp(-10.02 + 0.0364*t)

      ! Calculate beta parameter ([3], Eq. 5.34)
      beta = betam + dbeta

      ! Calculate alpha parameter ([3], Eq. 5.32)
      t_inv = 300./tk - 1
      alpha = (0.00504 + 0.0062*t_inv)*exp(-22.1*t_inv)

      ! Calculate real part of pure ice dielectric constant ([3], Eq. 5.30)
      eps_i_r = 3.1884 + 9.1e-4*t

      ! Imaginary part of pure-ice permittivity ([3], Eq. 5.31)
      eps_i_i = alpha/fghz + beta*fghz

      ! Calculate dielectric constant of pure ice ([3]).
      eps_i = eps_i_r - jj*eps_i_i

   END SUBROUTINE diel_ice

!-----------------------------------------------------------------------------

   SUBROUTINE diel_water(fghz, type, swc, t, wf_sand, wf_clay, BD_all, sal, eps_w)

!-----------------------------------------------------------------------------
      real(r8), intent(in)     :: fghz                              ! Microwave frequency [GHz]
      integer, intent(in)      :: type                              ! Type of water, 0: pure water, 1: sea water, 2: soil water
      real(r8), intent(in)     :: swc                               ! Soil water content (m3/m3)
      real(r8), intent(in)     :: t                                 ! Soil temperature (C)
      real(r8), intent(in)     :: wf_sand                           ! Gravimetric sand percent fraction(%)
      real(r8), intent(in)     :: wf_clay                           ! Gravimetric clay percent fraction(%)
      real(r8), intent(in)     :: BD_all                            ! Bulk density(g/cm3)
      real(r8), intent(in)     :: sal                               ! Water salinity (psu)
      complex(r8), intent(out) :: eps_w                             ! Dielectric constant of soil water

      real(r8)                 :: sigma                             ! Ionic conductivity (S/m)
      real(r8)                 :: a, b                              ! Parameters
      real(r8)                 :: tau_w                             ! Relaxation time of pure water
      real(r8)                 :: eps_w0                            ! Static dielectric constant of pure water
      real(r8)                 :: wc
      real(r8)                 :: omega                             ! Angular frequency [rad/s]
!-----------------------------------------------------------------------------
      omega = 2.0*pi*(fghz*1e9)

      ! Salinity-dependent relaxation time ([2], Eqs. 16-18).
      tau_w = 1.768e-11 - 6.086e-13*t + 1.104e-14*t**2 - 8.111e-17*t**3
      b = 1.000 + 2.282e-5*sal*t - 7.638e-4*sal - 7.760e-6*sal**2 + 1.105e-8*sal**3
      tau_w = tau_w*b

      ! Salinity-dependent static permittivity ([2], Eqs. 13-15).
      eps_w0 = 87.134 - 1.949e-1*t - 1.276e-2*t**2 + 2.491e-4*t**3
      a = 1.000 + 1.613e-5*sal*t - 3.656e-3*sal + 3.210e-5*sal**2 - 4.232e-7*sal**3
      eps_w0 = eps_w0*a

      IF (type == 0) THEN                                           
         ! Pure-water static permittivity ([2], Eq. 8).
         eps_w0 = 88.045 - 0.4147*t + 6.295e-4*t**2 + 1.075e-5*t**3
         eps_w = eps_w_inf + (eps_w0 - eps_w_inf)/(1 - jj*omega*tau_w)

      ELSEIF (type == 1) THEN                                       
         ! Ionic conductivity ([2], Eqs. 9-12).
         sigma =                                                                         &
            sal*(0.182521 - 1.46192e-3*sal + 2.09324e-5*sal**2 - 1.28205e-7*sal**3)      &
            *exp(-1.*(25 - t)*                                                           &
            (2.033e-2 + 1.266e-4*(25 - t) + 2.464e-6*(25 - t)**2                         &
            - sal*(1.849e-5 - 2.551e-7*(25 - t) + 2.551e-8*(25 - t)**2)))

         ! Saline-water Debye permittivity ([2], Eq. 5).
         eps_w = eps_w_inf + (eps_w0 - eps_w_inf)/(1 - jj*omega*tau_w) + jj*sigma/(omega*eps_0)

      ELSEIF (type == 2) THEN
         ! Soil-water conductivity ([5]).
         sigma = -1.645 + 1.939*BD_all - 0.02256*wf_sand + 0.01594*wf_clay
         
         ! Negative for very sandy soils with low bulk density
         IF (sigma < 0.) THEN
            sigma = 0.                                             
         ENDIF

         ! Calculate dielectric constant of soil-water by modified Debye expression
         wc = max(0.001, swc)
         eps_w = eps_w_inf + (eps_w0 - eps_w_inf)/(1 - jj*omega*tau_w) &
            + jj*sigma/(omega*eps_0)*(rho_soil - BD_all)/(rho_soil*wc)
      ENDIF

   END SUBROUTINE diel_water

!-----------------------------------------------------------------------------

   SUBROUTINE diel_soil_W80(ew, t, wc, wf_sand, wf_clay, porsl, eps)

!-----------------------------------------------------------------------------
      real(r8), intent(in)     :: t                                 ! Soil temperature (C)
      real(r8), intent(in)     :: wc                                ! Volumetric soil moisture (m3/m3)
      real(r8), intent(in)     :: wf_sand                           ! Gravimetric sand percent fraction(%)
      real(r8), intent(in)     :: wf_clay                           ! Gravimetric clay percent fraction(%)
      real(r8), intent(in)     :: porsl                             ! Soil porosity at surface
      complex(r8), intent(in)  :: ew                                ! Dielectric constant of water
      complex(r8), intent(out) :: eps

      real(r8)                 :: wp                                ! Wilting point
      real(r8)                 :: wt                                ! Transition moisture point (cm3/cm3)
      real(r8)                 :: gamma                             ! Fitting parameter
      real(r8)                 :: ecl                               ! Conductivity loss
      real(r8)                 :: alpha                             ! Conductivity loss parameter
      real(r8)                 :: sal_sea = 32.5                    ! Sea water salinity (psu)
      real(r8)                 :: sal_soil = 0.0                    ! Soil salinity (psu)
      complex(r8)              :: eps_x                             ! Dielectric constant of the initially absorbed water
      complex(r8)              :: eps_a = (1.0, 0.0)                ! Dielectric constant of air
      complex(r8)              :: eps_r = (5.5, 0.2)                ! Dielectric constant of rock
      complex(r8)              :: eps_i = (3.2, 0.1)                ! Dielectric constant of ice
      complex(r8)              :: eps_f = (5.0, 0.5)                ! Dielectric constant of frozen soil
!-----------------------------------------------------------------------------
      ! Calculate wilting point at the soil layer ([4], Eq. 1).
      wp = 0.06774 - 0.00064*wf_sand + 0.00478*wf_clay

      ! Calculate fitting parameters ([4], Eq. 8).
      gamma = -0.57*wp + 0.481

      ! Calculate transition moisture point ([4], Eq. 9).
      wt = 0.49*wp + 0.165

      ! Calculate dielectric constant of wet soil ([4], Eqs. 2-5).
      IF (wc <= wt) THEN
         eps_x = eps_i + (ew - eps_i)*(wc/wt)*gamma
         eps = wc*eps_x + (porsl - wc)*eps_a + (1.-porsl)*eps_r
      ELSE
         eps_x = eps_i + (ew - eps_i)*gamma
         eps = wt*eps_x + (wc - wt)*ew + (porsl - wc)*eps_a + (1.-porsl)*eps_r
      END IF

      ! Add conductivity loss for imaginary part ([4], Eq. 6).
      alpha = min(100.*wp, 26.)
      ecl = alpha*wc**2
      eps = eps + jj*ecl

   END SUBROUTINE diel_soil_W80

!-----------------------------------------------------------------------------

   SUBROUTINE diel_soil_M04(fghz, wc, wf_clay, eps)

!-----------------------------------------------------------------------------
      real(r8), intent(in)     :: fghz                              ! Microwave frequency [GHz]
      real(r8), intent(in)     :: wc                                ! Soil moisture (m3/m3)
      real(r8), intent(in)     :: wf_clay                           ! Gravimetric clay percent fraction(%)
      complex(r8), intent(out) :: eps

      real(r8)                 :: znd                               ! Dry-soil refractive index [-]
      real(r8)                 :: zkd                               ! Dry-soil extinction coefficient [-]
      real(r8)                 :: zxmvt                             ! Maximum bound-water volume fraction [m3/m3]
      real(r8)                 :: zep0b                             ! Static bound-water permittivity [-]
      real(r8)                 :: ztaub                             ! Bound-water relaxation time [s]
      real(r8)                 :: zsigmab                           ! Bound-water conductivity [S/m]
      real(r8)                 :: zep0u                             ! Static free-water permittivity [-]
      real(r8)                 :: ztauu                             ! Free-water relaxation time [s]
      real(r8)                 :: zsigmau                           ! Free-water conductivity [S/m]
      real(r8)                 :: zcxb                              ! Bound-water real Debye contribution [-]
      real(r8)                 :: zepwbx                            ! Real bound-water permittivity [-]
      real(r8)                 :: zepwby                            ! Bound-water dielectric loss factor [-]
      real(r8)                 :: zcxu                              ! Free-water real Debye contribution [-]
      real(r8)                 :: zepwux                            ! Real free-water permittivity [-]
      real(r8)                 :: zepwuy                            ! Free-water dielectric loss factor [-]
      real(r8)                 :: znb                               ! Bound-water refractive index [-]
      real(r8)                 :: zkb                               ! Bound-water extinction coefficient [-]
      real(r8)                 :: znu                               ! Free-water refractive index [-]
      real(r8)                 :: zku                               ! Free-water extinction coefficient [-]
      real(r8)                 :: zxmvt2                            ! Actual bound-water volume fraction [m3/m3]
      real(r8)                 :: znm                               ! Moist-soil refractive index [-]
      real(r8)                 :: zkm                               ! Moist-soil extinction coefficient [-]
      real(r8)                 :: zepmx                             ! Real moist-soil permittivity [-]
      real(r8)                 :: zepmy                             ! Moist-soil dielectric loss factor [-]
      integer                  :: zflag                             ! Free-water flag: 1 when wc >= zxmvt
      real(r8)                 :: f                                 ! Frequency [Hz]
!-----------------------------------------------------------------------------
      f = fghz*1e9

      ! Dry-soil refractive parameters ([9], Eq. 18).
      znd = 1.634 - 0.539 * (wf_clay/100) + 0.2748 * (wf_clay/100) ** 2
      zkd = 0.03952 - 0.04038 * (wf_clay / 100)

      ! Maximum bound water fraction ([9], Eq. 19).
      zxmvt = 0.02863 + 0.30673 * wf_clay / 100

      ! Bound water parameters ([9], Eqs. 20-22).
      zep0b   = 79.8 - 85.4  * (wf_clay / 100) + 32.7  * (wf_clay / 100) &
         *(wf_clay / 100)
      ztaub   = 1.062e-11 + 3.450e-12 * (wf_clay / 100)
      zsigmab = 0.3112 + 0.467 * (wf_clay / 100)

      ! Unbound (free) water parameters ([9], Eqs. 24-25).
      zep0u   = 100
      ztauu   = 8.5e-12
      zsigmau = 0.3631 + 1.217 * (wf_clay / 100)

      ! Computation of epsilon water (bound & unbound) ([9], Eq. 16).
      zcxb   = (zep0b - eps_w_inf) / (1. + (2.*pi*f*ztaub)**2)
      zepwbx = eps_w_inf + zcxb
      zepwby = zcxb * (2.*pi*f*ztaub) + zsigmab / (2.*pi*eps_0*f)
      zcxu   = (zep0u - eps_w_inf) / (1 + (2*pi*f*ztauu)**2)
      zepwux = eps_w_inf + zcxu
      zepwuy = zcxu * (2.*pi*f*ztauu) + zsigmau/(2.*pi*eps_0*f)

      ! Bound/free-water refractive indices ([9], Eqs. 14-15).
      znb = sqrt( sqrt( zepwbx**2 + zepwby**2) + zepwbx ) / sqrt(2.0)
      zkb = sqrt( sqrt( zepwbx**2 + zepwby**2) - zepwbx ) / sqrt(2.0)
      znu = sqrt( sqrt( zepwux**2 + zepwuy**2) + zepwux ) / sqrt(2.0)
      zku = sqrt( sqrt( zepwux**2 + zepwuy**2) - zepwux ) / sqrt(2.0)

      ! Computation of soil refractive index (nm & km) ([9], Eqs. 12-13).
      zxmvt2 = min(wc, zxmvt)
      zflag  = 0
      IF ( wc >= zxmvt ) zflag = 1
      znm = znd + (znb - 1) * zxmvt2 + (znu - 1) * (wc-zxmvt) * zflag
      zkm = zkd + zkb * zxmvt2 + zku * (wc-zxmvt) * zflag

      ! Computation of soil dielectric constant ([9], Eq. 11).
      zepmx = znm ** 2 - zkm ** 2
      zepmy = znm * zkm * 2
      eps   = cmplx(zepmx, zepmy, kind=r8)

   END SUBROUTINE diel_soil_M04

!-----------------------------------------------------------------------------

   SUBROUTINE diel_soil_M09(fghz, wc, t, wf_clay, eps)

!-----------------------------------------------------------------------------
      real(r8), intent(in)     :: fghz                              ! Microwave frequency [GHz]
      real(r8), intent(in)     :: t                                 ! Soil temperature
      real(r8), intent(in)     :: wc                                ! Soil moisture (m3/m3)
      real(r8), intent(in)     :: wf_clay                           ! Weighted fraction (%)
      complex(r8), intent(out) :: eps

      real(r8)                 :: nd                                ! Dry-soil refractive index [-]
      real(r8)                 :: kd                                ! Dry-soil extinction coefficient [-]
      real(r8)                 :: mvt                               ! Maximum bound-water volume fraction [m3/m3]
      real(r8)                 :: ts                                ! Reference temperature for parameter fits [degC]
      real(r8)                 :: e0b                               ! Static bound-water permittivity at ts [-]
      real(r8)                 :: Bb                                ! Bound-water permittivity temperature factor [1/K]
      real(r8)                 :: Bsgb                              ! Bound-water conductivity slope [S/m/K]
      real(r8)                 :: Fb                                ! Bound-water permittivity transform at ts [-]
      real(r8)                 :: eb0                               ! Static bound-water permittivity at t [-]
      real(r8)                 :: dSbR                              ! Bound-water activation entropy divided by R [-]
      real(r8)                 :: taub                              ! Bound-water relaxation time at t [s]
      real(r8)                 :: sigmabt                           ! Bound-water conductivity at ts [S/m]
      real(r8)                 :: sigmab                            ! Bound-water conductivity at t [S/m]
      real(r8)                 :: e0u                               ! Static free-water permittivity at ts [-]
      real(r8)                 :: Bu                                ! Free-water permittivity temperature factor [1/K]
      real(r8)                 :: Bsgu                              ! Free-water conductivity slope [S/m/K]
      real(r8)                 :: Fu                                ! Free-water permittivity transform at ts [-]
      real(r8)                 :: eu0                               ! Static free-water permittivity at t [-]
      real(r8)                 :: dHuR                              ! Free-water activation enthalpy divided by R [K]
      real(r8)                 :: dSuR                              ! Free-water activation entropy divided by R [-]
      real(r8)                 :: tauu                              ! Free-water relaxation time at t [s]
      real(r8)                 :: dHbR                              ! Bound-water activation enthalpy divided by R [K]
      real(r8)                 :: sigmau                            ! Free-water conductivity at t [S/m]
      real(r8)                 :: sigmaut                           ! Free-water conductivity at ts [S/m]
      real(r8)                 :: cxb                               ! Bound-water real Debye contribution [-]
      real(r8)                 :: eb_r                              ! Real bound-water permittivity [-]
      real(r8)                 :: eb_i                              ! Bound-water dielectric loss factor [-]
      real(r8)                 :: cxu                               ! Free-water real Debye contribution [-]
      real(r8)                 :: eu_r                              ! Real free-water permittivity [-]
      real(r8)                 :: eu_i                              ! Free-water dielectric loss factor [-]
      real(r8)                 :: nb                                ! Bound-water refractive index [-]
      real(r8)                 :: kb                                ! Bound-water extinction coefficient [-]
      real(r8)                 :: nu                                ! Free-water refractive index [-]
      real(r8)                 :: ku                                ! Free-water extinction coefficient [-]
      real(r8)                 :: nm                                ! Moist-soil refractive index [-]
      real(r8)                 :: km                                ! Moist-soil extinction coefficient [-]
      real(r8)                 :: eps_r                             ! Real moist-soil permittivity [-]
      real(r8)                 :: eps_i                             ! Moist-soil dielectric loss factor [-]
      real(r8)                 :: f                                 ! Frequency [Hz]
!-----------------------------------------------------------------------------
      f = fghz*1e9

      ! Dry-soil refractive parameters ([10], Eqs. 11-12).
      nd = 1.634 - 0.539e-2 * wf_clay + 0.2748e-4 * (wf_clay ** 2)
      kd = 0.03952 - 0.04038e-2 * wf_clay

      ! Maximum bound water fraction ([10], Eq. 13).
      mvt = 0.02863 + 0.30673e-2 * wf_clay

      ! Reference temperature for parameter fits ([10], p. 413).
      ts = 20.

      ! Eb0 computation ([10], Eqs. 7-8, 14-15, 23).
      e0b  = 79.8 - 85.4e-2 * wf_clay + 32.7e-4 * (wf_clay **2)
      Bb   = 8.67e-19 - 0.00126e-2 * wf_clay + 0.00184e-4 * (wf_clay ** 2)  - 9.77e-10*(wf_clay &
         **3) - 1.39e-15 *(wf_clay**4)
      Bsgb = 0.0028  + 0.02094e-2*wf_clay - 0.01229e-4*(wf_clay**2) - 5.03e-22*(wf_clay**3) + &
         4.163e-24*(wf_clay**4)
      Fb   = log((e0b - 1)/(e0b + 2))
      eb0  = (1 + 2*exp(Fb-Bb*(t-ts))) / (1 - exp(Fb-Bb*(t-ts)))

      ! Taub computation ([10], Eqs. 9, 18-19).
      dHbR = 1467 + 2697e-2*wf_clay - 980e-4 *(wf_clay**2) + 1.368e-10*(wf_clay**3) - 8.61e-13 &
         *(wf_clay**4)
      dSbR = 0.888 + 9.7e-2 *wf_clay - 4.262e-4*(wf_clay**2) + 6.79e-21 *(wf_clay**3) + &
         4.263e-22*(wf_clay**4)
      taub = 48e-12 * exp(dHbR/(t+tfrz)-dSbR)/(t+tfrz)

      ! Sigmab computation ([10], Eqs. 10, 22).
      sigmabt = 0.3112 + 0.467e-2*wf_clay
      sigmab  = sigmabt + Bsgb*(t-ts)

      ! Unbound (free) water parameters ([10], Eqs. 7-8, 16-17, 25).
      ! Eu0 computation.
      e0u  = 100.
      Bu   = 1.11e-4 - 1.603e-7 *wf_clay + 1.239e-9 *(wf_clay**2) + 8.33e-13 *(wf_clay**3) - &
         1.007e-14*(wf_clay**4)
      Bsgu = 0.00108 + 0.1413e-2*wf_clay - 0.2555e-4*(wf_clay**2) + 0.2147e-6*(wf_clay**3) - &
         0.0711e-8*(wf_clay**4)
      Fu   = log((e0u - 1)/(e0u + 2))
      eu0  = (1 + 2*exp(Fu-Bu*(t-ts))) / (1-exp(Fu-Bu*(t-ts)))

      ! Tauu computation ([10], Eqs. 9, 20-21).
      dHuR = 2231 - 143.1e-2 *wf_clay + 223.2e-4*(wf_clay**2) - 142.1e-6*(wf_clay**3) + &
         27.14e-8 *(wf_clay**4)
      dSuR = 3.649 - 0.4894e-2*wf_clay + 0.763e-4*(wf_clay**2) - 0.4859e-6*(wf_clay**3) + &
         0.0928e-8*(wf_clay**4)
      tauu = 48e-12 * exp(dHuR/(t+tfrz)-dSuR)/(t+tfrz)

      ! Sigmau computation ([10], Eqs. 10, 24).
      sigmaut = 0.05 + 1.4*(1.0 - (1.0 - wf_clay*1.e-2)**4.664)
      sigmau  = sigmaut + Bsgu*(t-ts)

      ! Bound/free-water permittivity ([10], Eq. 6), ([9], Eq. 16).
      cxb  = (eb0-eps_w_inf) / (1+(2*pi*f*taub)**2)
      eb_r = eps_w_inf + cxb
      eb_i = cxb*(2*pi*f*taub) + sigmab/(2*pi*eps_0*f)
      cxu  = (eu0-eps_w_inf) / (1+(2*pi*f*tauu)**2)
      eu_r = eps_w_inf + cxu
      eu_i = cxu*(2*pi*f*tauu) + sigmau/(2*pi*eps_0*f)

      ! Bound/free-water refractive indices ([10], Eq. 5).
      nb = sqrt(sqrt(eb_r**2+eb_i**2)+eb_r) / sqrt(2.0)
      kb = sqrt(sqrt(eb_r**2+eb_i**2)-eb_r) / sqrt(2.0)
      nu = sqrt(sqrt(eu_r**2+eu_i**2)+eu_r) / sqrt(2.0)
      ku = sqrt(sqrt(eu_r**2+eu_i**2)-eu_r) / sqrt(2.0)

      ! Computation of soil refractive index (nm & km) ([9], Eqs. 12-13).
      IF (wc <= mvt) THEN
         nm = nd + (nb-1)*wc
         km = kd + kb*wc
      ELSE
         nm = nd + (nb-1)*mvt + (nu-1)*(wc-mvt)
         km = kd + kb*mvt + ku*(wc-mvt)
      ENDIF

      ! Computation of soil dielectric constant ([10], Eq. 4).
      eps_r = nm**2 - km**2
      eps_i = 2* nm * km
      eps   = cmplx(eps_r, eps_i, kind=r8)

   END SUBROUTINE diel_soil_M09

!-----------------------------------------------------------------------------

   SUBROUTINE diel_soil_D85(ew, swc, wf_sand, wf_clay, BD_all, eps)

!-----------------------------------------------------------------------------
      complex(r8), intent(in)  :: ew                                ! Complex water permittivity [-]
      real(r8), intent(in)     :: swc                               ! Volumetric soil moisture [m3/m3]
      real(r8), intent(in)     :: wf_sand                           ! Sand mass fraction [%]
      real(r8), intent(in)     :: wf_clay                           ! Clay mass fraction [%]
      real(r8), intent(in)     :: BD_all                            ! Soil bulk density [g/cm3]
      complex(r8), intent(out) :: eps                               ! Complex moist-soil permittivity [-]

      real(r8)                 :: alphas = 0.65_r8                  ! Empirical dielectric mixing exponent [-]
      real(r8)                 :: beta                              ! Texture-dependent moisture exponent [-]
      real(r8)                 :: eaa                               ! Intermediate power-law mixture [-]
      real(r8)                 :: eps_s                             ! Solid soil-particle permittivity [-]
      real(r8)                 :: epsi                              ! Moist-soil dielectric loss factor [-]
      real(r8)                 :: epsr                              ! Real moist-soil permittivity [-]
      real(r8)                 :: wc                                ! Soil moisture used in mixing [m3/m3]
!-----------------------------------------------------------------------------
      ! Apply the existing lower bound to volumetric soil moisture.
      wc    = max(swc, 0.001_r8)

      ! Solid soil-particle permittivity from particle density ([5], Eq. 22).
      eps_s = (1.01_r8 + 0.44_r8 * rho_soil)**2.0_r8 - 0.062_r8

      ! Real-part dielectric mixing ([5], Eqs. 28, 30), ([6], Eq. 2), ([7]).
      ! Derive the moisture exponent from sand and clay mass percentages.
      beta  = (127.48_r8 - 0.519_r8 * wf_sand - 0.152_r8 * wf_clay) / 100.0_r8

      ! Combine air, soil particles and liquid water in power-law space.
      eaa   = 1.0_r8 + (BD_all / rho_soil) * (eps_s ** alphas - 1.0_r8) &
         &      + (wc ** beta) * (real(ew) ** alphas) - wc

      ! Invert the power-law transform to obtain real soil permittivity.
      epsr  = eaa ** (1.0_r8/alphas)

      ! Dielectric-loss mixing ([5], Eq. 31), ([6], Eqs. 3, 5).
      ! Derive the loss-part moisture exponent using percentages.
      beta  = (133.797_r8 - 0.603_r8 * wf_sand - 0.166_r8 * wf_clay) / &
         100.0_r8

      ! Weight the water dielectric loss by volumetric soil moisture.
      eaa   = (wc ** beta) * (abs(aimag(ew)) ** alphas)

      ! Invert the power-law transform to obtain the soil loss factor.
      epsi  = eaa ** (1.0_r8/alphas)

      ! Assemble complex soil permittivity with a positive imaginary part.
      eps   = cmplx(epsr, epsi, kind=r8)

   END SUBROUTINE diel_soil_D85

!-----------------------------------------------------------------------------

   SUBROUTINE smooth_reflectivity(theta, eps, r_s)

!-----------------------------------------------------------------------------
      real(r8), intent(in)     :: theta                             ! Incidence angle [rad]
      complex(r8), intent(in)  :: eps                               ! Dielectric constant of the surface
      real(r8), intent(out)    :: r_s(2)                            ! Reflectivities of flat surfaces for H and V polarizations

      complex(r8)              :: g                                 ! Parameter in Fresnel Law
!-----------------------------------------------------------------------------
      g = sqrt(eps - sin(theta)**2)
      r_s(1) = abs((cos(theta) - g)/(cos(theta) + g))**2.
      r_s(2) = abs((cos(theta)*eps - g)/(cos(theta)*eps + g))**2.

   END SUBROUTINE smooth_reflectivity

!-----------------------------------------------------------------------------

   SUBROUTINE rough_reflectivity(fghz, theta, is_desert, patchclass, r_s, r_r)

!-----------------------------------------------------------------------------
      real(r8), intent(in)     :: fghz                              ! Microwave frequency [GHz]
      real(r8), intent(in)     :: theta                             ! Incidence angle [rad]
      logical, intent(in)      :: is_desert                         ! Flag for desert soil
      integer, intent(in)      :: patchclass                        ! Patch class
      real(r8), intent(in)     :: r_s(2)                            ! Reflectivities of flat surfaces for H and V polarizations
      real(r8), intent(out)    :: r_r(2)                            ! Reflectivities of rough surfaces for H and V polarizations

      real(r8)                 :: Q                                 ! Parameter for polarization mixing
      real(r8)                 :: hr(N_land_classification)         ! Roughness parameter
      real(r8)                 :: nrh(N_land_classification)        ! Parameter for H polarization
      real(r8)                 :: nrv(N_land_classification)        ! Parameter for V polarization
      real(r8)                 :: kcm                               ! Wave number [rad/cm]
!-----------------------------------------------------------------------------
      kcm = (2*pi/(C/(fghz*1e9)))/100.0

      ! Calculate parameter for polarization mixing due to surface roughness ([16], Eq. 16).
      IF (is_desert) THEN
         r_r = r_s
      ELSE
         IF (fghz < 2.) THEN
            Q = 0.                                                
         ELSE
            Q = 0.35*(1.0 - exp(-0.6*rgh_surf**2*fghz))
         ENDIF

         ! Roughness parameters summarized in ([14], Table 2).
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

         ! Rough surface reflectivity for H and V polarizations ([14], Eq. 5).
         r_r(1) = (Q*r_s(2) + (1.-Q)*r_s(1))*exp(-hr(patchclass)*cos(theta) &
            **nrh(patchclass))
         r_r(2) = (Q*r_s(1) + (1.-Q)*r_s(2))*exp(-hr(patchclass)*cos(theta) &
            **nrv(patchclass))
      ENDIF

   END SUBROUTINE rough_reflectivity

!-----------------------------------------------------------------------------

   SUBROUTINE desert(fghz, t_soil, r_r, eps, tb_desert)

!-----------------------------------------------------------------------------
      real(r8), intent(in)     :: fghz                              ! Microwave frequency [GHz]
      real(r8), intent(in)     :: t_soil(2)                         ! Desert surface temperature (K) (h-pol and v-pol)
      real(r8), intent(in)     :: r_r(2)                            ! Reflectivity of rough surface  (h-pol and v-pol)
      complex(r8), intent(in)  :: eps                               ! Diel. const of desert
      real(r8), intent(out)    :: tb_desert(2)                      ! Brightness temperature of soil for H- and V- polarization

      real(r8)                 :: f0 = 0.7                          ! The fractional volume of spherical particles (f = (4/3)*pi*r^3*n0),
      real(r8)                 :: w                                 ! Single-particle albedo
      real(r8)                 :: g                                 ! Asymmetry parameter
      real(r8)                 :: a                                 ! Similarity parameter
      real(r8)                 :: em(2)                             ! Desert soil emissivity
      real(r8)                 :: y_r                               ! Real part of y-parameters
      real(r8)                 :: y_i                               ! Imaginary part of y-parameters
      real(r8)                 :: kr                                ! Size parameter for single-particle scattering albedo
!-----------------------------------------------------------------------------
      kr = (2*pi/(C/(fghz*1e9)))*(0.5*1e-3)

      ! Complex polarizability factors ([17], Eq. A15).
      y_r = (real(eps) - 1)/(real(eps) + 2)
      y_i = 3*aimag(eps)/(real(eps) + 2)**2

      ! Single-particle scattering albedo ([17], Eq. A16).
      w = (1 - f0)**4*kr**3*y_r**2/ &
         ((1 - f0)**4*kr**3*y_r**2 + 1.5*(1 + 2*f0)**2*y_i)

      ! Asymmetry parameter ([17], p. 374).
      g = 0.23*kr**2

      ! Similarity parameter ([17], Eq. 3b).
      a = sqrt((1 - w)/(1 - w*g))

      ! Desert emissivity ([17], Eq. 13).
      em = (1 - r_r)*(2*a/((1 + a) - (1 - a)*r_r))

      ! Calculate brightness temperature of desert
      tb_desert = t_soil*em

   END SUBROUTINE desert

!-----------------------------------------------------------------------------

   SUBROUTINE veg(theta, patchclass, lai, htop, snowdp, tleaf, tb_veg, gamma_p)

!-----------------------------------------------------------------------------
      real(r8), intent(in)     :: theta                             ! Incidence angle [rad]
      integer, intent(in)      :: patchclass                        ! Land cover class
      real(r8), intent(in)     :: lai                               ! Leaf area index
      real(r8), intent(in)     :: tleaf                             ! Leaf temperature (K)
      real(r8), intent(out)    :: tb_veg(2)                         ! Brightness temperature of vegetation for H- and V- polarization
      real(r8), intent(out)    :: gamma_p(2)                        ! Vegetation opacity for H- and V- polarization
      real(r8), intent(in)     :: htop, snowdp

      real(r8)                 :: tau_nadir                         ! Vegetation opacity at nadir
      real(r8)                 :: tau_veg(2)                        ! Vegetation opacity for H- and V- polarization
      integer                  :: i
!-----------------------------------------------------------------------------
      ! Calculate vegetation opacity (optical depth) at nadir b*VWC ([12], Eq. 22)
      IF (htop < snowdp) THEN
         tau_nadir = b1(patchclass)*lai + b2(patchclass)
      ELSE
         tau_nadir = b3(patchclass)                               
      ENDIF

      ! Calculate vegetation optical depth at H- and V- polarizations ([12], Eqs. 23-24).
      tau_veg(1) = tau_nadir*(cos(theta)**2 + tth(patchclass)*sin(theta)**2)
      tau_veg(2) = tau_nadir*(cos(theta)**2 + ttv(patchclass)*sin(theta)**2)

      ! Calculate brightness temperature of vegetation ([12], Eq. 15).
      DO i = 1, 2
         gamma_p(i) = exp(-tau_veg(i)/cos(theta))
         tb_veg(i) = (1.-w_CMEM(patchclass))*(1.-gamma_p(i))*tleaf
      ENDDO

   END SUBROUTINE veg

!-----------------------------------------------------------------------------

   SUBROUTINE snow(fghz, theta, t_snow, t, snowdp, rho_snow, liq_snow, r_sn, r_snow, tb_tos)

!-----------------------------------------------------------------------------
      real(r8), intent(in)     :: fghz                              ! Microwave frequency [GHz]
      real(r8), intent(in)     :: theta                             ! Incidence angle [rad]
      real(r8), intent(in)     :: t_snow                            ! Average snow temperature (K)
      real(r8), intent(in)     :: t                                 ! Temperature at bottom of snow (K), i.e., soil or leaf
      real(r8), intent(in)     :: snowdp                            ! Snow depth (m)
      real(r8), intent(in)     :: rho_snow                          ! Snow density (g/cm3)
      real(r8), intent(in)     :: liq_snow                          ! Snow liquid water content (cm3/cm3)
      real(r8), intent(in)     :: r_sn(2)                           ! Reflectivity between the snow and ground at (1, H-POL. 2, V.)
      real(r8), intent(out)    :: r_snow(2)                         ! Reflectivity between the snow and air for H- and V- polarization
      real(r8), intent(out)    :: tb_tos(2)                         ! Brightness temperature of snow-cover ground for H- and V- polarization

      real(r8)                 :: sal_snow = 0.0                    ! Snow salinity (pmm)
      real(r8)                 :: eps_i_r                           ! Real part of dielectric constant of ice
      real(r8)                 :: eps_i_i                           ! Imaginary part of dielectric constant of ice
      real(r8)                 :: eps_i_is                          ! Imaginary part of dielectric constant of impure ice -5(C)
      real(r8)                 :: eps_i_ip                          ! Imaginary part of dielectric constant of pure ice -5(C)
      real(r8)                 :: eps_ds_r                          ! Real part of dielectric constant of dry snow
      real(r8)                 :: eps_ds_i                          ! Imaginary part of dielectric constant of dry snow
      real(r8)                 :: eps_ws_i                          ! Imaginary part of dielectric constant of wet snow
      real(r8)                 :: eps_ws_r                          ! Real part of dielectric constant of wet snow
      real(r8)                 :: eps_w_s = 88.                     ! Dielectric constant of static water
      real(r8)                 :: eps_a_inf, eps_b_inf, eps_c_inf   ! Infinite frequency dielectric constant of three parts
      real(r8)                 :: eps_a_s, eps_b_s, eps_c_s         ! Static dielectric constant of three parts
      complex(r8)              :: eps_i                             ! Dielectric constant of ice
      complex(r8)              :: eps_a
      complex(r8)              :: eps_b
      complex(r8)              :: eps_c                             ! Dielectric constant of three parts
      complex(r8)              :: eps                               ! Dielectric constant of wet snow
      real(r8)                 :: rho_ds                            ! Density of dry snow  (g/cm3)
      real(r8)                 :: rho_i = 0.916                     ! Density of ice (g/cm3)
      real(r8)                 :: aa = 0.005
      real(r8)                 :: bb = 0.4975
      real(r8)                 :: cc = 0.4975                       ! Fitting parameters
      real(r8)                 :: fa, fb, fc                        ! Relaxation frequency of wet snow
      real(r8)                 :: d                                 ! Snow grain size (mm)
      real(r8)                 :: alpha, beta, pp, qq               ! Parameter used to calculate propagation angle in snow
      real(r8)                 :: theta_s                           ! Propagation angle in snow
      complex(r8)              :: z_s                               ! Wave impedance in snow
      real(r8)                 :: r_sa                              ! Reflectivity between the snow and air for H- and V- polarization
      real(r8)                 :: tb_2                              ! The net apparent temperature contributions due to emission by layers 2 (snow)
      real(r8)                 :: tb_3                              ! The net apparent temperature contributions due to emission by layers 3 (soil)
      real(r8)                 :: l2_apu                            ! Extinction coefficient of snow (Beer's Law)
      real(r8)                 :: q = 0.96                          ! Forward-scattering fraction
      real(r8)                 :: ka_ws                             ! Absorption coefficient of wet snow
      real(r8)                 :: ka_ds                             ! Absorption coefficient of dry snow
      real(r8)                 :: ke                                ! Extinction coefficient of wet snow
      real(r8)                 :: ke_ds                             ! Extinction coefficient of dry snow
      real(r8)                 :: ks                                ! Scattering coefficient of snow
      real(r8)                 :: b_ds, b_ws
      real(r8)                 :: wk_h                              ! [m], equal to cmem cmem_snow_set_var.F90:476
      integer                  :: i
      real(r8)                 :: omega                             ! Angular frequency [rad/s]
      real(r8)                 :: k                                 ! Wave number [rad/m]
!-----------------------------------------------------------------------------
      omega = 2.0*pi*(fghz*1e9)
      k = 2*pi/(C/(fghz*1e9))

      IF (snowdp > 0.01) THEN                                    
         ! Calculate dielectric constant of ice
         CALL diel_ice(fghz, t_snow - tfrz, eps_i)
         eps_i_r = real(eps_i)
         eps_i_i = aimag(eps_i)

         ! Consider the effect of salinity on the dielectric constant of ice
         eps_i_is = 0.0026/fghz + 0.00023*(fghz**0.87)              ! Impure ice -5(C)
         eps_i_ip = 6.e-4/fghz + 6.5e-5*(fghz**1.07)                ! Pure ice -5(C)
         eps_i_i = eps_i_i + (eps_i_is - eps_i_ip)*sal_snow/13.0d0  ! Corrected imaginary part of diel cons of ice

         ! Dry-snow air/ice dielectric mixing ([18]).
         rho_ds = (rho_snow - liq_snow)/(1.0 - liq_snow)            ! Calculate density of dry snow,
         eps_ds_r = 1.0 + 1.58*rho_ds/(1.0 - 0.365*rho_ds)
         eps_ds_i = 3.0*(rho_ds/rho_i)*eps_i_i*(eps_ds_r**2)*(2*eps_ds_r + 1)/ &
            ((eps_i_r + 2*eps_ds_r)*(eps_i_r + 2*eps_ds_r**2))      ! Negative imaginary part of diel cons of dry snow

         ! Wet-snow dielectric mixing ([18]).
         IF (liq_snow > 0.) THEN                                    ! Wet snow
            ! Three-component wet-snow relaxation frequencies.
            fa = f0w*(1 + (aa*(eps_w_s - eps_w_inf)/(eps_ds_r + (aa*(eps_w_inf - eps_ds_r)))))
            fb = f0w*(1 + (bb*(eps_w_s - eps_w_inf)/(eps_ds_r + (bb*(eps_w_inf - eps_ds_r)))))
            fc = f0w*(1 + (cc*(eps_w_s - eps_w_inf)/(eps_ds_r + (cc*(eps_w_inf - eps_ds_r)))))

            ! Calculate infinite frequency dielectric constant of three parts
            eps_a_inf = (liq_snow*(eps_w_inf - eps_ds_r)/3.)/ &
               (1.+aa*((eps_w_inf/eps_ds_r) - 1.))
            eps_b_inf = (liq_snow*(eps_w_inf - eps_ds_r)/3.)/ &
               (1.+bb*((eps_w_inf/eps_ds_r) - 1.))
            eps_c_inf = (liq_snow*(eps_w_inf - eps_ds_r)/3.)/ &
               (1.+cc*((eps_w_inf/eps_ds_r) - 1.))

            ! Calculate static dielectric constant of three parts
            eps_a_s = (liq_snow/3.)*(eps_w_s - eps_ds_r)/ &
               (1.+aa*((eps_w_s/eps_ds_r) - 1.))
            eps_b_s = (liq_snow/3.)*(eps_w_s - eps_ds_r)/ &
               (1.+bb*((eps_w_s/eps_ds_r) - 1.))
            eps_c_s = (liq_snow/3.)*(eps_w_s - eps_ds_r)/ &
               (1.+cc*((eps_w_s/eps_ds_r) - 1.))

            ! Three-component Debye relaxation ([18]).
            eps_a = eps_a_inf + (eps_a_s - eps_a_inf)/(1 + jj*fghz/fa)
            eps_b = eps_b_inf + (eps_b_s - eps_b_inf)/(1 + jj*fghz/fb)
            eps_c = eps_c_inf + (eps_c_s - eps_c_inf)/(1 + jj*fghz/fc)

            ! Calculate dielectric constant of wet snow
            eps = eps_a + eps_b + eps_c + (eps_ds_r - jj*eps_ds_i)
            eps_ws_r = real(eps)
            eps_ws_i = -1.*aimag(eps)
         ELSE
            eps_ws_r = eps_ds_r
            eps_ws_i = eps_ds_i
         ENDIF

         ! Calculate propagation angle in snow (change medium from air to snow)
         alpha = k*abs(aimag(sqrt(eps_ws_r - jj*eps_ws_i)))
         beta = k*real(sqrt(eps_ws_r - jj*eps_ws_i))
         pp = 2.*alpha*beta
         qq = beta**2 - alpha**2 - (k*k)*(sin(theta)**2)
         theta_s = atan(k*sin(theta)/((1./sqrt(2.)) &
            *sqrt(sqrt(pp**2.+qq**2.) + qq)))

         ! Calculate wave impedance in snow
         z_s = z0/sqrt(eps_ws_r - jj*eps_ws_i)

         ! Calculate brightness temperature above snow for H- V- polarization
         DO i = 1, 2
            ! Fresnel reflection coefficient between snow and air
            IF (i == 1) THEN
               r_sa = abs((z_s*cos(theta) - z0*cos(theta_s))/ &
                  (z_s*cos(theta) + z0*cos(theta_s)))**2
            ELSE
               r_sa = abs((z0*cos(theta) - z_s*cos(theta_s))/ &
                  (z0*cos(theta) + z_s*cos(theta_s)))**2
            ENDIF

            ! Density-based snow grain diameter; capped at 3 mm here ([19]).
            d = min(1000*(1.6e-4 + 1.1e-13*((rho_snow*1000.0)**4)), 3.0)

            ! Extinction coefficient of dry snow ([23], Eq. 14).
            ! Higher-frequency extinction fit ; L-band use is an extrapolation.
            ke_ds = 0.0018*(fghz**2.8)*(d**2)/4.3429

            ! Absorption coefficient of dry snow
            b_ds = (eps_ds_i/eps_ds_r)**2
            ka_ds = 2.*omega*sqrt(mu0*eps_0*eps_ds_r)* &
               sqrt(b_ds/(2.*(sqrt(1.+b_ds) + 1.)))

            IF (ke_ds < ka_ds) ke_ds = ka_ds

            ! Absorption coefficient of wet snow
            b_ws = (eps_ws_i/eps_ws_r)**2
            ka_ws = 2.*omega*sqrt(mu0*eps_0*eps_ws_r)* &
               sqrt(b_ws/(2.*(sqrt(1.+b_ws) + 1)))

            ! Total extinction (assuming scattering is the same for dry and wet snow)
            ke = (ke_ds - ka_ds) + ka_ws

            ! Scattering coefficient of dry and wet snow
            ks = ke_ds - ka_ds

            wk_h = snowdp/rho_snow
            IF (((ke - q*ks)*(1./cos(theta_s))*wk_h) > (log(HUGE(1.)) - 1.)) THEN
               l2_apu = sqrt(HUGE(1.))
            ELSE
               l2_apu = min(sqrt(HUGE(1.)), exp((ke - q*ks)*(1./cos(theta_s))*wk_h))
            ENDIF

            ! Snow-layer emission transmitted through the upper interface ([20]).
            tb_2 = (1.+r_sn(i)/l2_apu)*(1.-r_sa)*t_snow* &
               (ka_ws/(ke - q*ks))*(1.-1./l2_apu)/       &
               (1.-r_sn(i)*r_sa/l2_apu**2.)

            ! Underlying-layer emission transmitted through the snow ([20]).
            tb_3 = ((1.-r_sn(i))*(1.-r_sa)*t)/ &
               (l2_apu*(1.-r_sn(i)*r_sa/l2_apu**2.))

            ! Brightness temperature of snow-cover ground
            tb_tos(i) = tb_2 + tb_3

            ! Emission by layers 2 (snow) and 3 (soil)
            ! Reflectivity by layers 2 (snow) and 3 (soil)
            r_snow(i) = 1 - (tb_2/t_snow + tb_3/t)
         ENDDO
      ENDIF

   END SUBROUTINE snow
!-----------------------------------------------------------------------------
END MODULE MOD_DA_Operator_EM
#endif
