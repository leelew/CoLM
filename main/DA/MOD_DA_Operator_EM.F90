#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Operator_EM
!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Microwave emission model for top-of-atmosphere brightness temperatures.
!
!   EM prepares the land-surface inputs and combines the atmospheric, soil,
!   vegetation, and snow contributions of the emission model.
!
! WORKFLOW:
!   1. Prepare soil properties, temperatures, water contents, and snow states.
!   2. Evaluate the component processes for the supplied frequency and angle.
!   3. Combine emission, attenuation, and reflection to obtain H- and
!      V-polarized top-of-atmosphere brightness temperatures.
!
! NOTES:
!   Only soil patches (patchtype = 0) are evaluated; all other patch types
!   return a missing value.
!
!   Frequency, incidence angle, and the configured dielectric and roughness
!   schemes are passed explicitly to the component processes.
!   Effective soil temperature, surface reflection, vegetation, snow, and
!   atmosphere parameterizations are kept here with the emission workflow.
!   Water, ice, and soil dielectric calculations use PhysicsDiel, which can
!   also be called by other microwave observation operators.
!
! REFERENCES:
!   【1】Wigneron, J.-P., Kerr, Y., Waldteufel, P., Saleh, K., Escorihuela, M.-J.,
!       et al. (2007). L-band Microwave Emission of the Biosphere (L-MEB) Model:
!       Description and calibration against experimental data sets over crop
!       fields. Remote Sensing of Environment, 107(4), 639-655.
!       DOI: 10.1016/j.rse.2006.10.014
!   【2】Wigneron, J.-P., Jackson, T. J., O'Neill, P., De Lannoy, G., de Rosnay, P.,
!       et al. (2017). Modelling the passive microwave signature from land
!       surfaces: A review of recent results and application to the L-band SMOS &
!       SMAP soil moisture retrieval algorithms. Remote Sensing of Environment,
!       192, 238-262. DOI: 10.1016/j.rse.2017.01.024
!   【3】Njoku, E. G., & Kong, J.-A. (1977). Theory for passive microwave remote
!       sensing of near-surface soil moisture. Journal of Geophysical Research,
!       82(20), 3108-3118. DOI: 10.1029/JB082i020p03108
!   【4】Lv, S., Wen, J., Zeng, Y., Tian, H., & Su, Z. (2014). An improved
!       two-layer algorithm for estimating effective soil temperature in
!       microwave radiometry using in situ temperature and soil moisture
!       measurements. Remote Sensing of Environment, 152, 356-363.
!       DOI: 10.1016/j.rse.2014.07.007
!   【5】Kerr, Y. H., & Njoku, E. G. (1990). A semiempirical model for
!       interpreting microwave emission from semiarid land surfaces as seen from
!       space. IEEE Transactions on Geoscience and Remote Sensing, 28(3),
!       384-393. DOI: 10.1109/36.54364
!   【6】Grody, N. C., & Weng, F. (2008). Microwave emission and scattering from
!       deserts: Theory compared with satellite measurements. IEEE Transactions
!       on Geoscience and Remote Sensing, 46(2), 361-375.
!       DOI: 10.1109/TGRS.2007.909920
!   【7】Mätzler, C. (1998). Microwave permittivity of dry sand. IEEE
!       Transactions on Geoscience and Remote Sensing, 36(1), 317-319.
!       DOI: 10.1109/36.655342
!   【8】Mätzler, C. (1987). Applications of the interaction of microwaves with
!       the natural snow cover. Remote Sensing Reviews, 2(2), 259-387.
!       DOI: 10.1080/02757258709532086
!   【9】Anderson, E. A. (1976). A point energy and mass balance model of a snow
!       cover. NOAA Technical Report NWS 19, U.S. Department of Commerce,
!       Washington, D.C.
!   【10】Hallikainen, M. T., Ulaby, F. T., & Van Deventer, T. (1987). Extinction
!       behavior of dry snow in the 18- to 90-GHz range. IEEE Transactions on
!       Geoscience and Remote Sensing, GE-25(6), 737-745.
!       DOI: 10.1109/TGRS.1987.289743
!   【11】Ulaby, F. T., Moore, R. K., & Fung, A. K. (1986). Microwave Remote
!       Sensing: Active and Passive, Vol. III: From Theory to Applications.
!       Artech House, Norwood, MA.
!   【12】Pellarin, T., Wigneron, J.-P., Calvet, J.-C., Berger, M., Douville, H.,
!       et al. (2003). Two-year global simulation of L-band brightness
!       temperatures over land. IEEE Transactions on Geoscience and Remote
!       Sensing, 41(9), 2135-2139. DOI: 10.1109/TGRS.2003.815417
!
! HISTORY:
!   Lu Li, 12/2024: Initial version
!   Zhilong Fan, Lu Li, 03/2024: Debug and clean codes
!   Lu Li, 10/2025: Debug and clean codes
!   Lu Li, Zhilong Fan, 11/2026: DA v1.0
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Const_Physical
   USE MOD_Vars_Global
   USE MOD_Namelist
   USE MOD_DA_Operator_Const
   USE MOD_DA_Operator_PhysicsDiel
   IMPLICIT NONE
   SAVE

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE EM(patchtype, patchclass, dz_sno, forc_topo, htop, tref, t_soisno, tleaf, wliq_soisno, &
      wice_soisno, snowdp, lai, sai, wf_clay, wf_sand, wf_silt, BD_all, porsl, sat_theta, sat_fghz,  &
      tb_toa_h, tb_toa_v)

!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Compute top-of-atmosphere brightness temperatures for one soil patch.
!-----------------------------------------------------------------------------
!------------------------------ Dummy Arguments ------------------------------
      integer, intent(in)   :: patchtype                        ! Patch type (0 = soil)
      integer, intent(in)   :: patchclass                       ! Land cover class
      real(r8), intent(in)  :: dz_sno(maxsnl + 1:0)             ! Layer thickness (m)
      real(r8), intent(in)  :: forc_topo                        ! Topography [m]
      real(r8), intent(in)  :: htop                             ! Upper height of vegetation [m]
      real(r8), intent(in)  :: tref                             ! 2 m height air temperature [kelvin]
      real(r8), intent(in)  :: tleaf                            ! Leaf temperature [K]
      real(r8), intent(in)  :: t_soisno(maxsnl + 1:nl_soil)     ! Soil temperature [K]
      real(r8), intent(in)  :: wliq_soisno(maxsnl + 1:nl_soil)  ! Soil and snow liquid water content [kg/m2]
      real(r8), intent(in)  :: wice_soisno(maxsnl + 1:nl_soil)  ! Ice lens in layers [kg/m2]
      real(r8), intent(in)  :: snowdp                           ! Snow depth [meter]
      real(r8), intent(in)  :: lai                              ! Leaf area index
      real(r8), intent(in)  :: sai                              ! Stem area index
      real(r8), intent(in)  :: wf_clay(nl_soil)                 ! Gravimetric fraction of clay
      real(r8), intent(in)  :: wf_sand(nl_soil)                 ! Gravimetric fraction of sand
      real(r8), intent(in)  :: wf_silt(nl_soil)                 ! Gravimetric fraction of silt
      real(r8), intent(in)  :: BD_all(nl_soil)                  ! Bulk density of soil [kg/m3]
      real(r8), intent(in)  :: porsl(nl_soil)                   ! Fraction of soil that is voids [-]
      real(r8), intent(in)  :: sat_theta                        ! Incidence angle of satellite (rad)
      real(r8), intent(in)  :: sat_fghz                         ! Frequency of satellite (GHz)
      real(r8), intent(out) :: tb_toa_h                         ! Brightness temperature of top-of-atmosphere for H- polarization
      real(r8), intent(out) :: tb_toa_v                         ! Brightness temperature of top-of-atmosphere for V- polarization

!------------------------------ Local Variables ------------------------------
      logical               :: is_low_veg                       ! Flag for low vegetation
      real(r8)              :: dz_soisno(maxsnl+1:nl_soil)      ! Soil and snow layer thickness [m]
      integer               :: lb                               ! Lower bound of arrays
      real(r8)              :: tau_atm                          ! Atmospheric optical depth
      real(r8)              :: r_r(2)                           ! Rough surface reflectivity for H and V polarizations
      real(r8)              :: r_sn(2)                          ! Reflectivity between the snow and the ground
      real(r8)              :: r_snow(2)                        ! Reflectivity between the snow and the air
      real(r8)              :: tb_soil(2)                       ! Brightness temperature of soil for H and V polarizations
      real(r8)              :: tb_tos(2)                        ! Brightness temperature of snow-covered ground
      real(r8)              :: tb_tov(2)                        ! Brightness temperature of the vegetation
      real(r8)              :: tb_tov_noad(2)                   ! Vegetation brightness temperature without downwelling
      real(r8)              :: tb_au(2)                         ! Upwelling brightness temperature of the atmosphere 【12】Eq. (A2)
      real(r8)              :: tb_ad(2)                         ! Downwelling radiation (brightness temperature) of atmosphere
      real(r8)              :: rho_snow                         ! Snow density (g/cm3)
      real(r8)              :: liq_snow                         ! Snow liquid water content (cm3/cm3)
      real(r8)              :: gamma_p(2)                       ! Vegetation opacity for H- and V- polarization
      real(r8)              :: tb_veg(2)                        ! Brightness temperature of vegetation for H- and V- polarization
      real(r8)              :: tb_2(2)                          ! Downwelling vegetation emission reflected by the soil
      real(r8)              :: tb_3(2)                          ! Upwelling soil emission attenuated by the canopy
      real(r8)              :: tb_4(2)                          ! Downwelling sky radiation reflected by the soil
      real(r8)              :: tb_toa(2)                        ! Top-of-atmosphere brightness temperature
      real(r8)              :: wf_total(nl_soil)                ! Total gravimetric fraction
      real(r8)              :: BD_all_surf                      ! Bulk density of soil (g/cm3) at surface
      real(r8)              :: porsl_surf                       ! Soil porosity at surface
      real(r8)              :: t_surf                           ! Soil temperature at surface (C)
      real(r8)              :: t_deep                           ! Soil temperature at deep layer (C)
      real(r8)              :: liq_surf                         ! Liquid volumetric water content at surface (m3/m3)
      real(r8)              :: ice_surf                         ! Ice volumetric water content at surface (m3/m3)
      real(r8)              :: wf_clay_surf                     ! Gravimetric clay percent fraction(%) at surface
      real(r8)              :: wf_sand_surf                     ! Gravimetric sand percent fraction(%) at surface
      integer               :: i                                ! Loop index
!-----------------------------------------------------------------------------
! Prepare the parameters and states used by the operator
      ! Get the depths of the soil and snow layers
      dz_soisno(maxsnl+1:0) = dz_sno(maxsnl+1:0)
      dz_soisno(1:nl_soil) = dz_soi(1:nl_soil)

      ! Calculate the weighted surface parameters
      wf_total     = wf_clay + wf_sand + wf_silt
      wf_clay_surf = (wf_clay(1)/wf_total(1)*0.0175 + wf_clay(2)/wf_total(2)*0.0276)/0.0451*100
      wf_sand_surf = (wf_sand(1)/wf_total(1)*0.0175 + wf_sand(2)/wf_total(2)*0.0276)/0.0451*100
      BD_all_surf  = (BD_all(1)*0.0175 + BD_all(2)*0.0276)/0.0451/1000
      porsl_surf   = (porsl(1)*0.0175 + porsl(2)*0.0276)/0.0451

      ! Calculate the surface and deep soil temperatures [C]
      t_surf = ((t_soisno(1)*(0.0175) + t_soisno(2)*(0.0451 - 0.0175))/0.0451) - tfrz
      t_deep = ((t_soisno(7)*(0.8289-0.5) + t_soisno(8)*(1.0 - 0.8289))/0.5) - tfrz

      ! Calculate the liquid and ice volumetric water contents of the top two layers
      liq_surf = (wliq_soisno(1) + wliq_soisno(2))/(0.0451*denh2o)
      ice_surf = (wice_soisno(1) + wice_soisno(2))/(0.0451*denice)

      ! Calculate the lower bound of the snow layers
      lb = 0
      DO i = maxsnl+1, 0
         IF (wliq_soisno(i) + wice_soisno(i) > 0.0) THEN
            lb = i
            EXIT
         ENDIF
      ENDDO

!-----------------------------------------------------------------------------
! Run the forward operator
!-----------------------------------------------------------------------------
      ! Only soil patches are evaluated; other patch types are not computed.
      IF (patchtype /= 0) THEN
         tb_toa = spval
      ELSE

!-----------------------------------------------------------------------------
! Atmosphere module
!-----------------------------------------------------------------------------
         CALL atm(sat_theta, forc_topo, tref, tau_atm, tb_au, tb_ad)

!-----------------------------------------------------------------------------
! Soil module
!-----------------------------------------------------------------------------
         CALL EM_soil(sat_fghz, sat_theta, DEF_DA_OPERATOR_RTM_DIEL, DEF_DA_OPERATOR_RTM_ROUGH, patchclass, &
            t_surf, t_deep, liq_surf, ice_surf, wf_sand_surf, wf_clay_surf, BD_all_surf, porsl_surf, r_r,   &
            tb_soil)

!-----------------------------------------------------------------------------
! Vegetation and snow module
!   Four cases are distinguished for the calculations:
!   1) no vegetation and no snow
!   2) no vegetation with snow
!   3) vegetation without snow
!   4) vegetation with snow
!-----------------------------------------------------------------------------
         ! Roughly distinguish low and high vegetation (IGBP only)
         is_low_veg = .true.
         IF (patchclass >= 1 .and. patchclass <= 5) THEN
            is_low_veg = .false.
         END IF

         ! Limit the snow density to 1 g/cm3
         IF (snowdp > 0.01) THEN
            rho_snow = (wliq_soisno(lb) + wice_soisno(lb))/(dz_soisno(lb)*1e3)
            liq_snow = wliq_soisno(lb)*rho_snow/(dz_soisno(lb)*1e3)

            IF (liq_snow > 1.0 .or. rho_snow > 1.0) then
               rho_snow = 1.0
               liq_snow = wliq_soisno(lb)*rho_snow/(dz_soisno(lb)*1e3)
            END IF
         END IF

         ! Main procedure
         !-----------------------------------------------------------------------------
         ! 1) no veg and no snow
         !   two components:
         !      (1) brightness temperature of soil
         !      (2) the downwelling radiation reflected by the soil
         !-----------------------------------------------------------------------------
         IF ((lai + sai < 1e-6) .and. (snowdp < 0.01)) THEN
            tb_tov = tb_soil + tb_ad*r_r
            tb_tov_noad = tb_soil

         !-----------------------------------------------------------------------------
         ! 2) no veg and has snow
         !   two components:
         !      (1) brightness temperature of snow
         !      (2) the downwelling radiation reflected by the snow
         !-----------------------------------------------------------------------------
         ELSE IF ((lai + sai < 1e-6) .and. (snowdp > 0.01)) THEN
            ! Calculate the brightness temperature of snow-covered ground
            CALL snow(sat_fghz, sat_theta, t_soisno(1), t_soisno(1), snowdp, rho_snow, liq_snow, r_r, &
               r_snow, tb_tos)

            tb_tov = tb_tos + tb_ad*r_snow
            tb_tov_noad = tb_tos

         !-----------------------------------------------------------------------------
         ! 3) has veg and no snow
         !   four components:
         !      (1) the direct upwelling vegetation emission,
         !      (2) the downwelling vegetation emission reflected by the soil and attenuated by the canopy layer
         !      (3) upwelling soil emission attenuated by the canopy
         !      (4) the downwelling reflected by the soil and attenuated by the canopy layer
         !-----------------------------------------------------------------------------
         ELSE IF ((lai + sai > 1e-6) .and. (snowdp < 0.01)) THEN
            ! Calculate the brightness temperature of vegetation
            CALL veg(sat_theta, patchclass, lai, htop, 0.0, tleaf, tb_veg, gamma_p)

            DO i = 1, 2
               tb_2(i) = tb_veg(i)*gamma_p(i)*r_r(i)
               tb_3(i) = tb_soil(i)*gamma_p(i)
               tb_4(i) = tb_ad(i)*r_r(i)*(gamma_p(i)**2)
               tb_tov(i) = tb_veg(i) + tb_2(i) + tb_3(i) + tb_4(i)
               tb_tov_noad(i) = tb_veg(i) + tb_2(i) + tb_3(i)
            ENDDO

         !-----------------------------------------------------------------------------
         ! 4) has veg and has snow
         !   Determine the positional relationship between the vegetation and the snow.
         !
         !   If vegetation is higher than snow,
         !   first calculate the brightness temperature of the snow (soil boundary), then
         !   combine four components to derive the brightness temperature at the top of the
         !   vegetation:
         !      (1) the direct upwelling vegetation emission,
         !      (2) the downwelling vegetation emission reflected by the snow and attenuated by the canopy layer
         !      (3) upwelling snow emission attenuated by the canopy
         !      (4) the downwelling reflected by the snow and attenuated by the canopy layer
         !
         !   If vegetation is lower than snow
         !   first calculate the brightness temperature at the top of the vegetation (soil
         !   boundary), then combine four components to derive the brightness temperature at
         !   the top of the snow:
         !      (1) the direct upwelling vegetation emission,
         !      (2) the downwelling vegetation emission reflected by the soil and attenuated by the canopy layer
         !      (3) upwelling soil emission attenuated by the canopy
         !      (4) the downwelling reflected by the soil and attenuated by the canopy layer
         !-----------------------------------------------------------------------------
         ELSE IF ((lai + sai > 1e-6) .and. (snowdp > 0.01)) THEN
            IF (htop < snowdp) THEN
               ! Calculate the brightness temperature of low vegetation
               CALL veg(sat_theta, patchclass, lai, htop, snowdp, tleaf, tb_veg, gamma_p)

               ! Calculate the brightness temperature at the top of the vegetation
               DO i = 1, 2
                  tb_2(i) = tb_veg(i)*gamma_p(i)*r_r(i)
                  tb_3(i) = tb_soil(i)*gamma_p(i)
                  tb_4(i) = tb_ad(i)*r_r(i)*(gamma_p(i)**2)
                  tb_tov(i) = tb_veg(i) + tb_2(i) + tb_3(i) + tb_4(i)
                  tb_tov_noad(i) = tb_veg(i) + tb_2(i) + tb_3(i)
               ENDDO

               ! Calculate the reflectivity between the snow and low vegetation (adopted from CMEM)
               r_sn(:) = 1.0 - tb_tov_noad(:)/t_soisno(1)

               ! Calculate the brightness temperature of snow-covered ground
               CALL snow(sat_fghz, sat_theta, t_soisno(1), t_soisno(1), snowdp, rho_snow, liq_snow, r_sn, &
                  r_snow, tb_tos)

               ! Calculate the brightness temperature at the top of the snow
               tb_tov = tb_tos + tb_ad*r_snow
               tb_tov_noad = tb_tos

            ELSE
               ! Calculate the brightness temperature of snow-covered ground
               CALL snow(sat_fghz, sat_theta, t_soisno(1), t_soisno(1), snowdp, rho_snow, liq_snow, r_r, &
                  r_snow, tb_tos)

               ! Calculate the brightness temperature of high vegetation
               CALL veg(sat_theta, patchclass, lai, htop, snowdp, tleaf, tb_veg, gamma_p)

               ! Calculate the brightness temperature at the top of the vegetation
               DO i = 1, 2
                  tb_2(i) = tb_veg(i)*gamma_p(i)*r_snow(i)
                  tb_3(i) = tb_tos(i)*gamma_p(i)
                  tb_4(i) = tb_ad(i)*r_snow(i)*(gamma_p(i)**2)
                  tb_tov(i) = tb_veg(i) + tb_2(i) + tb_3(i) + tb_4(i)
                  tb_tov_noad(i) = tb_veg(i) + tb_2(i) + tb_3(i)
               ENDDO
            END IF
         END IF

!-----------------------------------------------------------------------------
! Calculate the top-of-atmosphere brightness temperature
!-----------------------------------------------------------------------------
         tb_toa = tb_tov*exp(-tau_atm) + tb_au

      END IF

      tb_toa_h = tb_toa(1)
      tb_toa_v = tb_toa(2)

   END SUBROUTINE EM

!-----------------------------------------------------------------------------

   SUBROUTINE EM_soil(fghz, theta, diel_scheme, rough_scheme, patchclass, t_surf, t_deep, liq_surf, ice_surf, &
      wf_sand_surf, wf_clay_surf, BD_all_surf, porsl_surf, r_r, tb_soil)

!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the brightness temperature of the soil surface.
!-----------------------------------------------------------------------------
!------------------------------ Dummy Arguments ------------------------------
      real(r8), intent(in)  :: fghz                           ! Microwave frequency [GHz]
      real(r8), intent(in)  :: theta                          ! Incidence angle [rad]
      integer, intent(in)   :: diel_scheme                    ! Soil dielectric mixing scheme
      integer, intent(in)   :: rough_scheme                   ! Surface roughness scheme
      integer, intent(in)   :: patchclass                     ! Land cover class
      real(r8), intent(in)  :: t_surf                         ! Soil temperature at surface (C)
      real(r8), intent(in)  :: t_deep                         ! Soil temperature at deep layer (C)
      real(r8), intent(in)  :: liq_surf                       ! Liquid volumetric water content at surface (m3/m3)
      real(r8), intent(in)  :: ice_surf                       ! Ice volumetric water content at surface (m3/m3)
      real(r8), intent(in)  :: wf_sand_surf                   ! Gravimetric sand percent fraction(%) at surface
      real(r8), intent(in)  :: wf_clay_surf                   ! Gravimetric clay percent fraction(%) at surface
      real(r8), intent(in)  :: BD_all_surf                    ! Bulk density of soil (g/cm3) at surface
      real(r8), intent(in)  :: porsl_surf                     ! Soil porosity at surface
      real(r8), intent(out) :: r_r(2)                         ! Rough surface reflectivity for H and V polarizations
      real(r8), intent(out) :: tb_soil(2)                     ! Brightness temperature of soil

!------------------------------ Local Variables ------------------------------
      real(r8)              :: t_eff(2)                       ! Effective temperature for H and V polarizations, [K]
      complex(r8)           :: eps_soil                       ! Dielectric constant of soil for H and V polarizations
      real(r8)              :: r_s(2)                         ! Smooth surface reflectivity for H and V polarizations
      complex(r8)           :: ew                             ! Dielectric constant of water
      logical               :: is_desert                      ! Flag for desert soil
      real(r8)              :: ffrz                           ! Fraction of frozen soil
      complex(r8)           :: eps_f = (5.0, 0.5)             ! Dielectric constant of frozen soil
      real(r8)              :: sal_soil = 0.0                 ! Soil salinity (psu)
!-----------------------------------------------------------------------------
      ! Check whether this patch is desert
      is_desert = .false.
      IF (liq_surf < 0.02 .and. wf_sand_surf > 90) THEN
         is_desert = .true.
      END IF

      ! Calculate the fraction of frozen soil
      IF (liq_surf + ice_surf <= 0.0d0) THEN
         ffrz = 0.0d0
      ELSE
         ffrz = ice_surf / (liq_surf + ice_surf)
      ENDIF

      ! Calculate the effective temperature 【4】
      CALL eff_soil_temp(liq_surf, t_surf, t_deep, t_eff)

      ! Calculate the dielectric constant of the soil mixture
      IF (is_desert) THEN
         ! Microwave permittivity of dry sand at 1-10 GHz 【7】Eq. (1)
         eps_soil = 2.53 + (2.79 - 2.53)/(1 - jj*(fghz/0.27)) + jj*0.002
      ELSE
         ! Define bulk density and porosity (CMEM)
         ! BD_all_surf = (wf_sand_surf*1.60d0 + wf_clay_surf*1.10d0 + (100.0d0 - wf_sand_surf - wf_clay_surf)*1.20d0)/100.0d0
         ! porsl_surf = 1.0d0 - BD_all_surf/2.660d0

         ! Calculate the dielectric constant of ice or water
         IF (ffrz > 0.95) THEN
            CALL diel_ice(fghz, t_surf, ew)
         ELSE
            CALL diel_water(fghz, 2, liq_surf, t_surf, wf_sand_surf, wf_clay_surf, BD_all_surf, sal_soil, &
               ew)
         END IF

         ! Calculate the dielectric constant of the mixed soil
         IF (diel_scheme == 0) THEN
            CALL diel_soil_W80(ew, t_surf, liq_surf, wf_sand_surf, wf_clay_surf, porsl_surf, eps_soil)
         ELSE IF (diel_scheme == 1) THEN
            CALL diel_soil_D85(ew, liq_surf, wf_sand_surf, wf_clay_surf, BD_all_surf, eps_soil)
         ELSE IF (diel_scheme == 2) THEN
            CALL diel_soil_M04(fghz, liq_surf, wf_clay_surf, eps_soil)
         ELSE IF (diel_scheme == 3) THEN
            CALL diel_soil_M09(fghz, liq_surf, t_surf, wf_clay_surf, eps_soil)
         ENDIF
      END IF

      ! Mix the dielectric constants of frozen and unfrozen soil
      eps_soil = eps_soil*(1.-ffrz) + eps_f*ffrz

      ! Calculate the smooth surface reflectivity
      CALL smooth_reflectivity(theta, eps_soil, r_s)

      ! Calculate the rough surface reflectivity
      CALL rough_reflectivity(fghz, theta, rough_scheme, is_desert, patchclass, r_s, r_r)

      ! Calculate the brightness temperature
      IF (is_desert) THEN
         CALL desert(fghz, t_eff, r_r, eps_soil, tb_soil)
      ELSE
         tb_soil = t_eff * (1 - r_r) ! 【1】
      END IF

   END SUBROUTINE EM_soil

!-----------------------------------------------------------------------------

   SUBROUTINE eff_soil_temp(wc_surf, t_surf, t_deep, t_eff)

!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the effective soil temperature with the two-layer algorithm of
!   Lv et al. (2014).
!-----------------------------------------------------------------------------
!------------------------------ Dummy Arguments ------------------------------
      real(r8), intent(in)  :: wc_surf                        ! Soil moisture at surface (m3/m3)
      real(r8), intent(in)  :: t_surf                         ! Soil temperature (C) at surface
      real(r8), intent(in)  :: t_deep                         ! Soil temperature (C) at deep layer
      real(r8), intent(out) :: t_eff(2)                       ! Effective temperature for H and V polarizations, [K]

!------------------------------ Local Variables ------------------------------
      ! Coefficient depending mainly on frequency and soil moisture 【4】; a larger
      ! C brings the effective temperature closer to the surface temperature
      real(r8)              :: C
      real(r8)              :: w0 = 0.30                      ! Reference soil moisture
      real(r8)              :: bw = 0.30                      ! Soil moisture exponent
!-----------------------------------------------------------------------------
      IF (wc_surf < 0.0) THEN
         C = 0.001
      ELSE
         C = max(0.001, (wc_surf/w0)**bw)
      ENDIF
      t_eff(:) = t_deep + (t_surf - t_deep)*C + tfrz

   END SUBROUTINE eff_soil_temp

!-----------------------------------------------------------------------------

   SUBROUTINE smooth_reflectivity(theta, eps, r_s)

!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the smooth surface reflectivity with the Fresnel equations
!   (Njoku and Kong, 1977).
!-----------------------------------------------------------------------------
!------------------------------ Dummy Arguments ------------------------------
      real(r8), intent(in)    :: theta                        ! Incidence angle [rad]
      complex(r8), intent(in) :: eps                          ! Dielectric constant of the surface
      real(r8), intent(out)   :: r_s(2)                       ! Reflectivities of flat surfaces for H and V polarizations

!------------------------------ Local Variables ------------------------------
      complex(r8)             :: g                            ! Intermediate Fresnel parameter
!-----------------------------------------------------------------------------
      ! Fresnel reflection coefficients for H and V polarizations 【3】
      g = sqrt(eps - sin(theta)**2)
      r_s(1) = abs((cos(theta) - g)/(cos(theta) + g))**2.
      r_s(2) = abs((cos(theta)*eps - g)/(cos(theta)*eps + g))**2.

   END SUBROUTINE smooth_reflectivity

!-----------------------------------------------------------------------------

   SUBROUTINE rough_reflectivity(fghz, theta, rough_scheme, is_desert, patchclass, r_s, r_r)

!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the rough surface reflectivity from the Q/h model of Kerr and
!   Njoku (1990) with the land-cover parameters compiled by Wigneron et al.
!   (2017).
!-----------------------------------------------------------------------------
!------------------------------ Dummy Arguments ------------------------------
      real(r8), intent(in)  :: fghz                           ! Microwave frequency [GHz]
      real(r8), intent(in)  :: theta                          ! Incidence angle [rad]
      integer, intent(in)   :: rough_scheme                   ! Surface roughness scheme
      logical, intent(in)   :: is_desert                      ! Flag for desert soil
      integer, intent(in)   :: patchclass                     ! Land cover class
      real(r8), intent(in)  :: r_s(2)                         ! Reflectivities of flat surfaces for H and V polarizations
      real(r8), intent(out) :: r_r(2)                         ! Reflectivities of rough surfaces for H and V polarizations

!------------------------------ Local Variables ------------------------------
      real(r8)              :: f                              ! Microwave frequency [Hz]
      real(r8)              :: lam                            ! Wavelength [m]
      real(r8)              :: k                              ! Wave number [rad/m]
      real(r8)              :: kcm                            ! Wave number [rad/cm]
      real(r8)              :: Q                              ! Parameter for polarization mixing
      real(r8)              :: hr(N_land_classification)      ! Roughness parameter
      real(r8)              :: nrh(N_land_classification)     ! Parameter for H polarization
      real(r8)              :: nrv(N_land_classification)     ! Parameter for V polarization
!-----------------------------------------------------------------------------
      f = fghz*1e9
      lam = C/f
      k = 2*pi/lam
      kcm = k/100.0

      IF (is_desert) THEN
         r_r = r_s
      ELSE
         ! Calculate the polarization mixing parameter
         IF (fghz < 2.) THEN
            Q = 0. ! Q is assumed zero at low frequency
         ELSE
            Q = 0.35*(1.0 - exp(-0.6*rgh_surf**2*fghz)) ! 【5】Eq. (16)
         END IF

         ! Calculate the rough surface reflectivity (default settings of 【2】)
         IF (rough_scheme == 0) THEN
            hr(:) = (2.0*kcm*rgh_surf)**2.0
            nrh(:) = 0.0
            nrv(:) = 0.0
         ELSE IF (rough_scheme == 1) THEN
            hr(:) = hr_SMOS
            nrh(:) = 2.0
            nrv(:) = 0.0
         ELSE IF (rough_scheme == 2) THEN
            hr(:) = hr_SMAP
            nrh(:) = 2.0
            nrv(:) = 2.0
         ELSE IF (rough_scheme == 3) THEN
            hr(:) = hr_P16
            nrh(:) = -1.0
            nrv(:) = -1.0
         END IF

         ! Calculate the rough surface reflectivity for H and V polarizations 【5】Eqs. (16)-(17)
         r_r(1) = (Q*r_s(2) + (1.-Q)*r_s(1))*exp(-hr(patchclass)*cos(theta)**nrh(patchclass))
         r_r(2) = (Q*r_s(1) + (1.-Q)*r_s(2))*exp(-hr(patchclass)*cos(theta)**nrv(patchclass))
      END IF

   END SUBROUTINE rough_reflectivity

!-----------------------------------------------------------------------------

   SUBROUTINE desert(fghz, t_soil, r_r, eps, tb_desert)

!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Calculate desert emissivity with the Grody and Weng (2008) model.
!-----------------------------------------------------------------------------
!------------------------------ Dummy Arguments ------------------------------
      real(r8), intent(in)    :: fghz                         ! Microwave frequency [GHz]
      real(r8), intent(in)    :: t_soil(2)                    ! Desert surface temperature (K) (h-pol and v-pol)
      real(r8), intent(in)    :: r_r(2)                       ! Reflectivity of rough surface (h-pol and v-pol)
      complex(r8), intent(in) :: eps                          ! Dielectric constant of desert
      real(r8), intent(out)   :: tb_desert(2)                 ! Brightness temperature of soil for H- and V- polarization

!------------------------------ Local Variables ------------------------------
      real(r8)                :: f                            ! Microwave frequency [Hz]
      real(r8)                :: lam                          ! Wavelength [m]
      real(r8)                :: k                            ! Wave number [rad/m]
      real(r8)                :: kr                           ! Particle size parameter
      ! Fractional volume of spherical particles: f = (4/3)*pi*r^3*n0, where
      ! r = 0.5 mm is the particle radius and n0 is the number density
      real(r8)                :: f0 = 0.7                     ! Fractional volume of particles
      real(r8)                :: w                            ! Single-particle albedo
      real(r8)                :: g                            ! Asymmetry parameter
      real(r8)                :: a                            ! Similarity parameter
      real(r8)                :: em(2)                        ! Desert soil emissivity
      real(r8)                :: y_r                          ! Real part of the y parameter
      real(r8)                :: y_i                          ! Imaginary part of the y parameter
!-----------------------------------------------------------------------------
      ! Calculate the y parameters
      f = fghz*1e9
      lam = C/f
      k = 2*pi/lam
      kr = k*(0.5*1e-3)

      y_r = (real(eps) - 1)/(real(eps) + 2)  ! 【6】Eq. (A15)
      y_i = 3*aimag(eps)/(real(eps) + 2)**2  ! 【6】Eq. (A15)

      ! Calculate the single-particle albedo
      w = (1 - f0)**4*kr**3*y_r**2/ &
         ((1 - f0)**4*kr**3*y_r**2 + 1.5*(1 + 2*f0)**2*y_i) ! 【6】Eq. (A16)

      ! Calculate the asymmetry parameter
      g = 0.23*kr**2 ! 【6】p. 374

      ! Calculate the similarity parameter
      a = sqrt((1 - w)/(1 - w*g)) ! 【6】Eq. (3b)

      ! Calculate the desert soil emissivity
      em = (1 - r_r)*(2*a/((1 + a) - (1 - a)*r_r)) ! 【6】Eq. (A13)

      ! Calculate the brightness temperature of desert
      tb_desert = t_soil*em

   END SUBROUTINE desert

!-----------------------------------------------------------------------------

   SUBROUTINE veg(theta, patchclass, lai, htop, snowdp, tleaf, tb_veg, gamma_p)

!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the vegetation opacity and the brightness temperature at the
!   top of the vegetation canopy (L-MEB).
!-----------------------------------------------------------------------------
!------------------------------ Dummy Arguments ------------------------------
      real(r8), intent(in)  :: theta                          ! Incidence angle [rad]
      integer, intent(in)   :: patchclass                     ! Land cover class
      real(r8), intent(in)  :: lai                            ! Leaf area index
      real(r8), intent(in)  :: tleaf                          ! Leaf temperature (K)
      real(r8), intent(out) :: tb_veg(2)                      ! Brightness temperature of vegetation for H- and V- polarization
      real(r8), intent(in)  :: htop                           ! Canopy top height [m]
      real(r8), intent(in)  :: snowdp                         ! Snow depth [m]
      real(r8), intent(out) :: gamma_p(2)                     ! Vegetation opacity for H- and V- polarization

!------------------------------ Local Variables ------------------------------
      real(r8)              :: tau_nadir                      ! Vegetation opacity at nadir
      real(r8)              :: tau_veg(2)                     ! Vegetation opacity for H- and V- polarization
      integer               :: i                              ! Loop index
!-----------------------------------------------------------------------------
      ! Calculate the vegetation opacity (optical depth) at nadir
      IF (htop < snowdp) THEN
         tau_nadir = b1(patchclass)*lai + b2(patchclass) ! Low vegetation ! 【1】Eq. (22)
      ELSE
         tau_nadir = b3(patchclass) ! High vegetation
      END IF

      ! Calculate the vegetation optical depth at H and V polarizations
      tau_veg(1) = tau_nadir*(cos(theta)**2 + tth(patchclass)*sin(theta)**2)  ! 【1】Eq. (23)
      tau_veg(2) = tau_nadir*(cos(theta)**2 + ttv(patchclass)*sin(theta)**2)  ! 【1】Eq. (24)

      ! Calculate the brightness temperature of vegetation
      DO i = 1, 2
         gamma_p(i) = exp(-tau_veg(i)/cos(theta)) ! 【1】Eq. (15)
         tb_veg(i) = (1.-w_CMEM(patchclass))*(1.-gamma_p(i))*tleaf
      ENDDO

   END SUBROUTINE veg

!-----------------------------------------------------------------------------

   SUBROUTINE snow(fghz, theta, t_snow, t, snowdp, rho_snow, liq_snow, r_sn, r_snow, tb_tos)

!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the brightness temperature of snow-covered ground.
!-----------------------------------------------------------------------------
!------------------------------ Dummy Arguments ------------------------------
      real(r8), intent(in)  :: fghz                             ! Microwave frequency [GHz]
      real(r8), intent(in)  :: theta                            ! Incidence angle [rad]
      real(r8), intent(in)  :: t_snow                           ! Average snow temperature (K)
      real(r8), intent(in)  :: t                                ! Temperature at bottom of snow (K), i.e., soil or leaf
      real(r8), intent(in)  :: snowdp                           ! Snow depth (m)
      real(r8), intent(in)  :: rho_snow                         ! Snow density (g/cm3)
      real(r8), intent(in)  :: liq_snow                         ! Snow liquid water content (cm3/cm3)
      real(r8), intent(in)  :: r_sn(2)                          ! Reflectivity between the snow and the ground
      real(r8), intent(out) :: r_snow(2)                        ! Reflectivity between the snow and air for H- and V- polarization
      real(r8), intent(out) :: tb_tos(2)                        ! Brightness temperature of snow-covered ground

!------------------------------ Local Variables ------------------------------
      real(r8)              :: f                                ! Microwave frequency [Hz]
      real(r8)              :: omega                            ! Angular frequency [rad/s]
      real(r8)              :: lam                              ! Wavelength [m]
      real(r8)              :: k                                ! Wave number [rad/m]
      real(r8)              :: sal_snow = 0.0                   ! Snow salinity [ppm]
      real(r8)              :: eps_i_r                          ! Real part of dielectric constant of ice
      real(r8)              :: eps_i_i                          ! Imaginary part of dielectric constant of ice
      real(r8)              :: eps_i_is                         ! Imaginary part of dielectric constant of impure ice -5(C)
      real(r8)              :: eps_i_ip                         ! Imaginary part of dielectric constant of pure ice -5(C)
      real(r8)              :: eps_ds_r                         ! Real part of dielectric constant of dry snow
      real(r8)              :: eps_ds_i                         ! Imaginary part of dielectric constant of dry snow
      real(r8)              :: eps_ws_i                         ! Imaginary part of dielectric constant of wet snow
      real(r8)              :: eps_ws_r                         ! Real part of dielectric constant of wet snow
      real(r8)              :: eps_w_s = 88.                    ! Static dielectric constant of water
      real(r8)              :: eps_a_inf, eps_b_inf, eps_c_inf  ! Infinite frequency dielectric constant of three parts
      real(r8)              :: eps_a_s, eps_b_s, eps_c_s        ! Static dielectric constant of three parts
      complex(r8)           :: eps_i                            ! Dielectric constant of ice
      complex(r8)           :: eps_a                            ! Dielectric constant of part a
      complex(r8)           :: eps_b                            ! Dielectric constant of part b
      complex(r8)           :: eps_c                            ! Dielectric constant of part c
      complex(r8)           :: eps                              ! Dielectric constant of wet snow
      real(r8)              :: rho_ds                           ! Density of dry snow (g/cm3)
      real(r8)              :: rho_i = 0.916                    ! Density of ice (g/cm3)
      real(r8)              :: aa = 0.005                       ! Fitting parameter a
      real(r8)              :: bb = 0.4975                      ! Fitting parameter b
      real(r8)              :: cc = 0.4975                      ! Fitting parameter c
      real(r8)              :: fa, fb, fc                       ! Relaxation frequencies of wet-snow parts
      real(r8)              :: d                                ! Snow grain size [mm]
      real(r8)              :: alpha, beta, pp, qq              ! Parameters for the propagation angle in snow
      real(r8)              :: theta_s                          ! Propagation angle in snow
      complex(r8)           :: z_s                              ! Wave impedance in snow
      real(r8)              :: r_sa                             ! Reflectivity between the snow and air for H- and V- polarization
      real(r8)              :: tb_2                             ! Net emission of the snow layer
      real(r8)              :: tb_3                             ! Net emission of the soil layer
      real(r8)              :: l2_apu                           ! Two-way attenuation factor of the snow layer
      real(r8)              :: q = 0.96                         ! Scattering efficiency parameter
      real(r8)              :: ka_ws                            ! Absorption coefficient of wet snow
      real(r8)              :: ka_ds                            ! Absorption coefficient of dry snow
      real(r8)              :: ke                               ! Extinction coefficient of wet snow
      real(r8)              :: ke_ds                            ! Extinction coefficient of dry snow
      real(r8)              :: ks                               ! Scattering coefficient of snow
      real(r8)              :: b_ds, b_ws                       ! Loss-tangent terms of dry and wet snow
      real(r8)              :: wk_h                             ! Snow water path [m] (cf. CMEM)
      integer               :: i                                ! Loop index
!-----------------------------------------------------------------------------
      f = fghz*1e9
      omega = 2.0*pi*f
      lam = C/f
      k = 2*pi/lam

      IF (snowdp > 0.01) THEN ! Snow depth larger than 1 cm
         ! Calculate the dielectric constant of ice
         CALL diel_ice(fghz, t_snow - tfrz, eps_i)
         eps_i_r = real(eps_i)
         eps_i_i = aimag(eps_i)

         ! Consider the effect of salinity on the dielectric constant of ice
         eps_i_is = 0.0026/fghz + 0.00023*(fghz**0.87)              ! Impure ice at -5 C
         eps_i_ip = 6.e-4/fghz + 6.5e-5*(fghz**1.07)                ! Pure ice at -5 C
         eps_i_i = eps_i_i + (eps_i_is - eps_i_ip)*sal_snow/13.0d0  ! Corrected imaginary part of the ice dielectric constant

         ! Calculate the dielectric constant of dry snow with the Polder-van Santen
         ! mixing model 【8】
         rho_ds = (rho_snow - liq_snow)/(1.0 - liq_snow) ! Calculate the density of dry snow
         eps_ds_r = 1.0 + 1.58*rho_ds/(1.0 - 0.365*rho_ds)
         eps_ds_i = 3.0*(rho_ds/rho_i)*eps_i_i*(eps_ds_r**2)*(2*eps_ds_r + 1)/ &
            ((eps_i_r + 2*eps_ds_r)*(eps_i_r + 2*eps_ds_r**2)) ! Negative imaginary part of the dry-snow dielectric constant

         ! Calculate the dielectric constant of wet snow (Mätzler 1987) 【8】
         IF (liq_snow > 0.) THEN ! Wet snow
            ! Calculate the relaxation frequency of the three parts 【8】Eq. (2.26)
            fa = f0w*(1 + (aa*(eps_w_s - eps_w_inf)/(eps_ds_r + (aa*(eps_w_inf - eps_ds_r)))))
            fb = f0w*(1 + (bb*(eps_w_s - eps_w_inf)/(eps_ds_r + (bb*(eps_w_inf - eps_ds_r)))))
            fc = f0w*(1 + (cc*(eps_w_s - eps_w_inf)/(eps_ds_r + (cc*(eps_w_inf - eps_ds_r)))))

            ! Calculate the infinite-frequency dielectric constant of the three parts
            eps_a_inf = (liq_snow*(eps_w_inf - eps_ds_r)/3.)/ &
               (1.+aa*((eps_w_inf/eps_ds_r) - 1.))
            eps_b_inf = (liq_snow*(eps_w_inf - eps_ds_r)/3.)/ &
               (1.+bb*((eps_w_inf/eps_ds_r) - 1.))
            eps_c_inf = (liq_snow*(eps_w_inf - eps_ds_r)/3.)/ &
               (1.+cc*((eps_w_inf/eps_ds_r) - 1.))

            ! Calculate the static dielectric constant of the three parts
            eps_a_s = (liq_snow/3.)*(eps_w_s - eps_ds_r)/ &
               (1.+aa*((eps_w_s/eps_ds_r) - 1.))
            eps_b_s = (liq_snow/3.)*(eps_w_s - eps_ds_r)/ &
               (1.+bb*((eps_w_s/eps_ds_r) - 1.))
            eps_c_s = (liq_snow/3.)*(eps_w_s - eps_ds_r)/ &
               (1.+cc*((eps_w_s/eps_ds_r) - 1.))

            ! Debye equations
            eps_a = eps_a_inf + (eps_a_s - eps_a_inf)/(1 + jj*fghz/fa)
            eps_b = eps_b_inf + (eps_b_s - eps_b_inf)/(1 + jj*fghz/fb)
            eps_c = eps_c_inf + (eps_c_s - eps_c_inf)/(1 + jj*fghz/fc)

            ! Calculate the dielectric constant of wet snow
            eps = eps_a + eps_b + eps_c + (eps_ds_r - jj*eps_ds_i)
            eps_ws_r = real(eps)
            eps_ws_i = -1.*aimag(eps)
         ELSE
            eps_ws_r = eps_ds_r
            eps_ws_i = eps_ds_i
         END IF

         ! Calculate the propagation angle in snow
         alpha = k*abs(aimag(sqrt(eps_ws_r - jj*eps_ws_i)))
         beta = k*real(sqrt(eps_ws_r - jj*eps_ws_i))
         pp = 2.*alpha*beta
         qq = beta**2 - alpha**2 - (k*k)*(sin(theta)**2)
         theta_s = atan(k*sin(theta)/((1./sqrt(2.)) &
            *sqrt(sqrt(pp**2.+qq**2.) + qq)))

         ! Calculate the wave impedance in snow
         z_s = z0/sqrt(eps_ws_r - jj*eps_ws_i)

         ! Calculate the brightness temperature above snow for H and V polarizations
         DO i = 1, 2
            ! Fresnel reflection coefficient between snow and air 【8】
            IF (i == 1) THEN
               r_sa = abs((z_s*cos(theta) - z0*cos(theta_s))/ &
                  (z_s*cos(theta) + z0*cos(theta_s)))**2
            ELSE
               r_sa = abs((z0*cos(theta) - z_s*cos(theta_s))/ &
                  (z0*cos(theta) + z_s*cos(theta_s)))**2
            END IF

            ! Calculate the snow grain size 【9】Eq. (5.1)
            d = min(1000*(1.6e-4 + 1.1e-13*((rho_snow*1000.0)**4)), 3.0)

            ! Extinction coefficient of dry snow
            ! TODO: The formula was not derived for L-band and may not be suitable.
            ke_ds = 0.0018*(fghz**2.8)*(d**2)/4.3429 ! 【10】Eq. (14)

            ! Absorption coefficient of dry snow
            b_ds = (eps_ds_i/eps_ds_r)**2
            ka_ds = 2.*omega*sqrt(mu0*eps_0*eps_ds_r)* &
               sqrt(b_ds/(2.*(sqrt(1.+b_ds) + 1.)))

            IF (ke_ds < ka_ds) ke_ds = ka_ds

            ! Absorption coefficient of wet snow
            b_ws = (eps_ws_i/eps_ws_r)**2
            ka_ws = 2.*omega*sqrt(mu0*eps_0*eps_ws_r)* &
               sqrt(b_ws/(2.*(sqrt(1.+b_ws) + 1)))

            ! Total extinction (assuming the same scattering for dry and wet snow)
            ke = (ke_ds - ka_ds) + ka_ws

            ! Scattering coefficient of dry and wet snow
            ks = ke_ds - ka_ds

            wk_h = snowdp/rho_snow
            IF (((ke - q*ks)*(1./cos(theta_s))*wk_h) > (log(HUGE(1.)) - 1.)) THEN
               l2_apu = sqrt(HUGE(1.))
            ELSE
               l2_apu = min(sqrt(HUGE(1.)), exp((ke - q*ks)*(1./cos(theta_s))*wk_h))
            END IF

            ! Brightness temperature through the snow-air layer 【11】Eq. (4.161)
            tb_2 = (1.+r_sn(i)/l2_apu)*(1.-r_sa)*t_snow* &
               (ka_ws/(ke - q*ks))*(1.-1./l2_apu)/       &
               (1.-r_sn(i)*r_sa/l2_apu**2.)

            ! Brightness temperature through the soil-snow layer 【11】Eq. (4.162)
            tb_3 = ((1.-r_sn(i))*(1.-r_sa)*t)/ &
               (l2_apu*(1.-r_sn(i)*r_sa/l2_apu**2.))

            ! Brightness temperature of snow-covered ground
            tb_tos(i) = tb_2 + tb_3

            ! Emission by layers 2 (snow) and 3 (soil)
            ! Reflectivity of layers 2 (snow) and 3 (soil)
            r_snow(i) = 1 - (tb_2/t_snow + tb_3/t)
         ENDDO
      END IF

   END SUBROUTINE snow

!-----------------------------------------------------------------------------

   SUBROUTINE atm(theta, z, tref, tau_atm, tb_au, tb_ad)

!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the atmospheric optical depth and the upwelling and downwelling
!   brightness temperatures of the atmosphere (Pellarin et al., 2003).
!-----------------------------------------------------------------------------
!------------------------------ Dummy Arguments ------------------------------
      real(r8), intent(in)  :: theta                          ! Incidence angle [rad]
      real(r8), intent(in)  :: z                              ! Altitude [m]
      real(r8), intent(in)  :: tref                           ! 2 m air temperature [K]
      real(r8), intent(out) :: tau_atm                        ! Atmospheric optical depth
      real(r8), intent(out) :: tb_au(2)                       ! Upwelling brightness temperature of the atmosphere 【12】Eq. (A2)
      real(r8), intent(out) :: tb_ad(2)                       ! Downwelling radiation (brightness temperature) of atmosphere

!------------------------------ Local Variables ------------------------------
      real(r8)              :: t_sky = 2.7                    ! Cosmic background temperature [K]
      real(r8)              :: t_eq                           ! Equivalent layer temperature
      real(r8)              :: gossat                         ! Atmospheric transmissivity
!-----------------------------------------------------------------------------
      ! Calculate the optical depth of the atmosphere 【12】Eq. (A1)
      tau_atm = exp(-3.9262 - 0.2211*z/1000 - 0.00369*tref)/cos(theta)
      gossat = exp(-tau_atm)

      ! Calculate the equivalent layer temperature 【12】
      t_eq = exp(4.9274 + 0.002195*tref)

      ! Upwelling brightness temperature of the atmosphere 【12】Eq. (A2)
      tb_au(:) = t_eq*(1.-gossat)

      ! Downwelling radiation (brightness temperature) of the atmosphere 【12】Eq. (A2)
      tb_ad(:) = t_eq*(1.-gossat) + t_sky*gossat

   END SUBROUTINE atm

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Operator_EM
#endif

