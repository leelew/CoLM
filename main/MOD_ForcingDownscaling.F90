#include <define.h>

MODULE MOD_ForcingDownscaling

!-----------------------------------------------------------------------------
! !DESCRIPTION:
!  Downscaling meteorological forcings
!
! !INITIAL:
!  The Community Land Model version 5.0 (CLM5.0)
!
! !REVISIONS:
!  Updated by Yongjiu Dai, January 16, 2023
!  Revised by Lu Li, January 24, 2024
!  Revised by Sisi Chen, Lu Li, June, 2024
!-----------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Qsadv
   USE MOD_Namelist
   USE MOD_Const_Physical
   USE MOD_Vars_Global
   USE MOD_UserDefFun
   USE MOD_SpatialMapping
   USE MOD_Vars_TimeInvariants
   USE MOD_Vars_TimeVariables
   USE MOD_Vars_1DForcing
   USE MOD_OrbCoszen
   USE MOD_OrbCosazi
   USE MOD_DataType
   USE MOD_SpatialMapping
   USE MOD_SPMD_Task

   IMPLICIT NONE

   SAVE

   PUBLIC :: downscale_forcings 

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE downscale_forcings (&
         mapping,                           &
         julian_day,                        &
         forc_elv_part,  forc_maxelv_part,  &
         forc_t_part,    forc_q_part,       &
         forc_pbot_part, forc_rhoair_part,  &
         forc_prc_part,  forc_prl_part,     &
         forc_frl_part,  forc_swrad_part,   &
         forc_us_part,   forc_vs_part
         )

!-----------------------------------------------------------------------------
! DESCRIPTION:
!  Downscale atmospheric forcings using elevation differences 
!  and subgrid topographic characteristics and deep generative models.
!
! REFERENCES:
!  Chen,S., Li,L., & Dai, Y (2024). 
!  Exploring topography downscaling methods for hyper‐resolution land surface modeling.
!  Journal of Geophysical Research: Atmospheres, 129, e2024JD041338.
!
!  Chen,S., Li,L., & Dai, Y (2026).
!  Deep generative downscaling of precipitation for hyper-Resolution land surface modeling.
!  Journal of Advances in Modeling Earth Systems.
!-----------------------------------------------------------------------------

   IMPLICIT NONE

!------------------------ Dummy Argument -------------------------------------
   type(spatial_mapping_type), intent(in):: mapping                ! mapping of forcing grids to patch (also part)
   real(r8), intent(in)                  :: julian_day             ! day of year

   type(pointer_real8_1d), intent(in)    :: forc_elv_part(:)       ! map surface height of coarse forcing grid at part[m]
   type(pointer_real8_1d), intent(in)    :: forc_maxelv_part(:)    ! map max surface height of coarse forcing grid at part [m]

   type(pointer_real8_1d), intent(inout) :: forc_t_part(:)         ! atmospheric temperature [Kelvin]
   type(pointer_real8_1d), intent(inout) :: forc_q_part(:)         ! atmospheric specific humidity [kg/kg]
   type(pointer_real8_1d), intent(inout) :: forc_pbot_part(:)      ! atmospheric pressure [Pa]
   type(pointer_real8_1d), intent(inout) :: forc_rhoair_part(:)    ! atmospheric density [kg/m**3]
   type(pointer_real8_1d), intent(inout) :: forc_prc_part(:)       ! column convective precipitation [mm/s]
   type(pointer_real8_1d), intent(inout) :: forc_prl_part(:)       ! column large-scale precipitation [mm/s]
   type(pointer_real8_1d), intent(inout) :: forc_frl_part(:)       ! column downward longwave [W/m**2]
   type(pointer_real8_1d), intent(inout) :: forc_swrad_part(:)     ! column downward shortwave [W/m**2]
   type(pointer_real8_1d), intent(inout) :: forc_us_part(:)        ! column eastward wind [m/s]
   type(pointer_real8_1d), intent(inout) :: forc_vs_part(:)        ! column northward wind [m/s]

!------------------------ Local Variables ------------------------------------
   ! forcing fields before downscaling
   real(r8) :: forc_t_g      
   real(r8) :: forc_q_g      
   real(r8) :: forc_pbot_g   
   real(r8) :: forc_rhoair_g 
   real(r8) :: forc_prc_g   
   real(r8) :: forc_prl_g    
   real(r8) :: forc_frl_g    
   real(r8) :: forc_swrad_g  
   real(r8) :: forc_hgt_g    
   real(r8) :: forc_us_g    
   real(r8) :: forc_vs_g     

   ! forcing fields after downscaling
   real(r8) :: forc_t_c     
   real(r8) :: forc_q_c     
   real(r8) :: forc_pbot_c  
   real(r8) :: forc_rhoair_c
   real(r8) :: forc_prc_c   
   real(r8) :: forc_prl_c   
   real(r8) :: forc_frl_c   
   real(r8) :: forc_swrad_c 
   real(r8) :: forc_us_c    
   real(r8) :: forc_vs_c   

   ! terrain parameters
   real(r8) :: elv_g
   real(r8) :: elv_c
!-----------------------------------------------------------------------------

   IF (p_is_worker) THEN
      DO np = 1, mapping%npset
         DO ipart = 1, mapping%npart(np)

            ! Prepare terrain factors
#ifdef SinglePoint
            elv_g = elvf
#else
            elv_g = forc_elv_part(np)%val(ipart)
#endif
            elv_c = forc_topo(np)

            slp_c = slp(:,np)
            area_c = area(:,np)
            cur_c = cur(np)
            IF (DEF_DS_MODE == 'C') THEN
               asp_c = asp(:,np)
               svf_c = svf(np)
               sf_c = sf(:,:,np)
            ELSE
               asp_c(1) = 0.0_r8*PI/180     ! north
               asp_c(2) = 45.0_r8*PI/180    ! northeast
               asp_c(3) = 90.0_r8*PI/180    ! east
               asp_c(4) = 135.0_r8*PI/180   ! southeast
               asp_c(5) = 180.0_r8*PI/180   ! south
               asp_c(6) = 225.0_r8*PI/180   ! southwest
               asp_c(7) = 270.0_r8*PI/180   ! west
               asp_c(8) = 315.0_r8*PI/180   ! northwest
               asp_c(9) = -9999.0_r8        ! flat
               svf_c = -9999.
               sf_c = -9999.
            ENDIF

            balb = alb
            IF (forc_sols(np)+forc_solsd(np)+forc_soll(np)+forc_solld(np) == 0.) THEN
               balb = 0
            ELSE
               balb = ( alb(1,1,np)*forc_sols (np) + alb(1,2,np)*forc_solsd(np)   &
                     + alb(2,1,np)*forc_soll (np) + alb(2,2,np)*forc_solld(np) ) &
                  / (forc_sols(np)+forc_solsd(np)+forc_soll(np)+forc_solld(np))
            ENDIF

            ! Prepare 
            coszen = orb_coszen(julian_day, patchlonr(np), patchlatr(np))
            cosazi = orb_cosazi(julian_day, patchlonr(np), patchlatr(np), coszen)

            ! Prepare 
            forc_t_g = forc_t_part(np)%val(ipart)
            forc_pbot_g = forc_pbot_part(np)%val(ipart)
            forc_q_g = forc_q_part(np)%val(ipart)
            forc_frl_g = forc_frl_part(np)%val(ipart)

            ! Downscaling
            CALL downscale_t (&
                     elv_g, forc_t_g, &
                     elv_c, forc_t_c)

            CALL downscale_pbot (&
                     elv_g, forc_t_g, forc_pbot_g, &
                     elv_c, forc_t_c, forc_pbot_c)

            CALL downscale_q (&
                     forc_t_g, forc_pbot_g, forc_q_g, &
                     forc_t_c, forc_pbot_c, forc_q_c)

            CALL downscale_rhoair (&
                     forc_t_c, forc_pbot_c, forc_q_c, forc_rhoair_c)

            CALL downscale_frl (&
                     elv_g, forc_t_g, forc_q_g, forc_pbot_g, forc_frl_g, &
                     elv_c, forc_t_c, forc_q_c, forc_pbot_c, forc_frl_c)

            CALL downscale_swrad (&
                     julian_day, coszen, cosazi,              &
                     slp_c, asp_c, area_c, svf_c, sf_c, balb, &
                     forc_pbot_g, forc_swrad_g,               &
                     forc_pbot_c, forc_swrad_c)

            CALL downscale_ws (&
                     )

            
         ENDDO
      ENDDO
   ENDIF




   END SUBROUTINE downscale_forcings

!-----------------------------------------------------------------------------

   SUBROUTINE downscale_t (elv_g, forc_t_g, elv_c, forc_t_c)

!-----------------------------------------------------------------------------
! DESCRIPTION:
!  Adjust air temperature for the elevation difference between the
!  forcing grid and the target patch using a fixed surface lapse rate:
!
! NOTES:
!  - The lapse rate is set to 6 K km-1.
!  - Surface lapse rates vary with location, season, and time of day.
!    The fixed value is a simplification.
!  - Minder et al. (2010) report annual mean lapse rates of
!    3.9-5.2 K km-1 on the windward side of the Cascade Mountains,
!    below the commonly assumed value of 6.5 K km-1. Rates range
!    from 2.5-3.5 K km-1 for late-summer minimum temperatures to
!    6.5-7.5 K km-1 for spring maximum temperatures. Windward and
!    leeward slopes also show substantial differences.
!  - Kunkel (1989) reports monthly lapse rates (January-December)
!    of 4.4, 5.9, 7.1, 7.8, 8.1, 8.2, 8.1, 8.1, 7.7, 6.8, 5.5,
!    and 4.7 K km-1, illustrating seasonal variability.
!-----------------------------------------------------------------------------

   IMPLICIT NONE

   real(r8), intent(in)  :: elv_g                    ! forcing-grid elevation [m]
   real(r8), intent(in)  :: elv_c                    ! target column elevation [m]
   real(r8), intent(in)  :: forc_t_g                 ! forcing-grid air temperature [K]
   real(r8), intent(out) :: forc_t_c                 ! column temperature [K]

   real(r8), parameter   :: lapse_rate = 0.006_r8    ! temperature lapse rate [K/m]
!-----------------------------------------------------------------------------

      forc_t_c = forc_t_g - lapse_rate * (elv_c - elv_g)

   END SUBROUTINE downscale_t

!-----------------------------------------------------------------------------

   SUBROUTINE downscale_pbot (elv_g, forc_t_g, forc_pbot_g, elv_c, forc_t_c, forc_pbot_c)

!-----------------------------------------------------------------------------
! DESCRIPTION:
!  Adjust air pressure from the forcing-grid elevation to the target
!  column elevation using an exponential pressure-height relation.
!-----------------------------------------------------------------------------

   IMPLICIT NONE

   real(r8), intent(in)  :: elv_g        ! forcing-grid elevation [m]
   real(r8), intent(in)  :: elv_c        ! target column elevation [m]
   real(r8), intent(in)  :: forc_t_g     ! forcing-grid air temperature [K]
   real(r8), intent(in)  :: forc_t_c     ! column air temperature [K]
   real(r8), intent(in)  :: forc_pbot_g  ! forcing-grid air pressure [Pa]
   real(r8), intent(out) :: forc_pbot_c  ! column air pressure [Pa]

   real(r8) :: shbot                     ! pressure scale height [m]
!-----------------------------------------------------------------------------

      ! Scale height based on the mean air temperature [m].
      shbot = rgas * 0.5_r8 * (forc_t_g + forc_t_c) / grav

      ! Pressure at column elevation from the hydrostatic relation [Pa].
      forc_pbot_c = forc_pbot_g * exp(-(elv_c-elv_g)/shbot)

   END SUBROUTINE downscale_pbot

!-----------------------------------------------------------------------------

   SUBROUTINE downscale_q (forc_t_g, forc_pbot_g, forc_q_g, forc_t_c, forc_pbot_c, forc_q_c)

!-----------------------------------------------------------------------------
! DESCRIPTION:
!  Adjust specific humidity using the ratio of saturation specific
!  humidity at the column and forcing-grid temperature and pressure.
!-----------------------------------------------------------------------------

   IMPLICIT NONE

   real(r8), intent(in)  :: forc_t_g     ! forcing-grid air temperature [K]
   real(r8), intent(in)  :: forc_pbot_g  ! forcing-grid air pressure [Pa]
   real(r8), intent(in)  :: forc_q_g     ! forcing-grid specific humidity [kg/kg]
   real(r8), intent(in)  :: forc_t_c     ! column air temperature [K]
   real(r8), intent(in)  :: forc_pbot_c  ! column air pressure [Pa]
   real(r8), intent(out) :: forc_q_c     ! column specific humidity [kg/kg]

   real(r8) :: es_g, es_c                ! saturation vapor pressure [Pa]
   real(r8) :: qs_g, qs_c                ! saturation specific humidity [kg/kg]
   real(r8) :: dum1, dum2
!-----------------------------------------------------------------------------

      ! Saturation specific humidity at forcing-grid conditions.
      CALL Qsadv(forc_t_g, forc_pbot_g, es_g, dum1, qs_g, dum2)

      ! Saturation specific humidity at column conditions.
      CALL Qsadv(forc_t_c, forc_pbot_c, es_c, dum1, qs_c, dum2)

      ! Preserve the grid relative humidity while adjusting specific humidity.
      forc_q_c = forc_q_g * (qs_c / qs_g)

   END SUBROUTINE downscale_q

!-----------------------------------------------------------------------------

   SUBROUTINE downscale_rhoair (forc_t_c, forc_q_c, forc_pbot_c, forc_rhoair_c)

!-----------------------------------------------------------------------------
! DESCRIPTION:
!  Calculate column air density from the adjusted specific humidity,
!  air pressure, and air temperature.
!-----------------------------------------------------------------------------

   IMPLICIT NONE

   real(r8), intent(in)  :: forc_q_c      ! column specific humidity [kg/kg]
   real(r8), intent(in)  :: forc_pbot_c   ! column air pressure [Pa]
   real(r8), intent(in)  :: forc_t_c      ! column air temperature [K]
   real(r8), intent(out) :: forc_rhoair_c ! column moist air density [kg/m3]

   real(r8) :: e_c                        ! water-vapor partial pressure [Pa]
   real(r8) :: pd_c                       ! dry-air partial pressure [Pa]
   real(r8) :: rho_dryair                 ! dry-air density [kg/m3]
   real(r8) :: rho_vapor                  ! water-vapor density [kg/m3]
!-----------------------------------------------------------------------------

      ! Convert specific humidity to water-vapor partial pressure.
      e_c = forc_q_c * forc_pbot_c / (0.622_r8 + 0.378_r8 * forc_q_c)

      ! The total pressure is the sum of dry-air and water-vapor pressures.
      pd_c = forc_pbot_c - e_c

      ! Calculate the density of each gas using its own gas constant.
      rho_dryair = pd_c / (rgas * forc_t_c)
      rho_vapor = e_c / (rwat * forc_t_c)

      ! Moist-air density is the sum of dry-air and water-vapor densities.
      forc_rhoair_c = rho_dryair + rho_vapor

   END SUBROUTINE downscale_rhoair

!-----------------------------------------------------------------------------

   SUBROUTINE downscale_frl (&
      elv_g, forc_t_g, forc_q_g, forc_pbot_g, forc_frl_g, &
      elv_c, forc_t_c, forc_q_c, forc_pbot_c, forc_frl_c)

!-----------------------------------------------------------------------------
! DESCRIPTION:
!  Downscale longwave radiation

! REFERENCE:
!  Fiddes and Gruber, 2014
!  TopoSCALE v.1.0: downscaling gridded climate data in complex terrain. 
!  Geosci. Model Dev., 7, 387-405.  
!-----------------------------------------------------------------------------

   IMPLICIT NONE

   real(r8), intent(in) :: elv_g           ! atmospheric surface height (m)
   real(r8), intent(in) :: forc_t_g        ! forcing-grid atmospheric temperature [Kelvin]
   real(r8), intent(in) :: forc_q_g        ! forcing-grid atmospheric specific humidity [kg/kg]
   real(r8), intent(in) :: forc_pbot_g     ! forcing-grid atmospheric pressure [Pa]
   real(r8), intent(in) :: forc_frl_g      ! forcing-grid downward longwave (W/m**2)

   real(r8), intent(in) :: elv_c           ! column surface height (m)
   real(r8), intent(in) :: forc_t_c        ! column atmospheric temperature [Kelvin]
   real(r8), intent(in) :: forc_q_c        ! column atmospheric specific humidity [kg/kg]
   real(r8), intent(in) :: forc_pbot_c     ! column atmospheric pressure [Pa]
   real(r8), intent(out):: forc_frl_c      ! column downward longwave (W/m**2)

   real(r8) :: e_g                         ! forcing-grid water vapor pressure [Pa]
   real(r8) :: e_c                         ! column water vapor pressure [Pa]
   real(r8) :: emis_clearsky_g             ! forcing-grid clear-sky emissivity [-]
   real(r8) :: emis_clearsky_c             ! column clear-sky emissivity [-]
   real(r8) :: emis_allsky_g               ! forcing-grid all-sky emissivity [-]

   real(r8), parameter :: frl_bd = 0.5_r8  ! boundary factor of adjustment
!--------------------------------------------------------------------------

      ! Cal vapor pressure [Pa] at forcing-grid & column
      e_g = forc_q_g * forc_pbot_g / (0.622_r8 + 0.378_r8 * forc_q_g)
      e_c = forc_q_c * forc_pbot_c / (0.622_r8 + 0.378_r8 * forc_q_c)

      ! Cal clear-sky emissivity from vapor pressure [hPa].
      emis_clearsky_g = 0.23_r8 + 0.43_r8*((e_g/100)*forc_t_g)**(1._r8/5.7_r8)
      emis_clearsky_c = 0.23_r8 + 0.43_r8*((e_c/100)/forc_t_c)**(1._r8/5.7_r8)

      ! Infer grid-cell all-sky emissivity from the incoming longwave forcing.
      emis_allsky_g = forc_frl_g / (5.67e-8_r8 * forc_t_g**4)

      ! Apply the column clear-sky emissivity while retaining the grid cloudy-sky emissivity.
      forc_frl_c = (emis_clearsky_c + (emis_allsky_g - emis_clearsky_g)) * 5.67e-8_r8*forc_t_c**4
   
      ! But ensure that we don't depart too far from the atmospheric forcing value:
      ! negative values of lwrad are certainly bad, but small positive values might
      ! also be bad. We can especially run into trouble due to the normalization: a
      ! small lwrad value in one column can lead to a big normalization factor,
      ! leading to huge lwrad values in other columns.
      forc_frl_c = min(forc_frl_c, forc_frl_g * (1._r8 + frl_bd))
      forc_frl_c = max(forc_frl_c, forc_frl_g * (1._r8 - frl_bd))

   END SUBROUTINE downscale_frl

!-----------------------------------------------------------------------------

   SUBROUTINE downscale_swrad (&
      julian_day, coszen, cosazi,   &
      slp, asp, area, svf, sf, balb, &
      forc_pbot_g, forc_swrad_g,    &
      forc_pbot_c, forc_swrad_c)

!-----------------------------------------------------------------------------
! DESCRIPTION:
!  Downscale shortwave radiation

! REFERENCE:
!  Rouf, T., Mei, Y., Maggioni, V., Houser, P., & Noonan, M. (2020). 
!  A Physically Based Atmospheric Variables Downscaling Technique. 
!  Journal of Hydrometeorology, 21(1), 93-108.
!
!  Chen,S., Li,L., & (2024).
!  Exploring topography downscaling methods for hyper-resolution land surface modeling.
!  Journal of Geophysical Research: Atmospheres, 129, e2024JD041338.
!-----------------------------------------------------------------------------

   IMPLICIT NONE

   integer,  intent(in) :: julian_day
   real(r8), intent(in) :: coszen
   real(r8), intent(in) :: cosazi

   real(r8), intent(in) :: slp(:)
   real(r8), intent(in) :: asp(:)
   real(r8), intent(in) :: area(:)
   real(r8), intent(in), optional :: svf
   real(r8), intent(in), optional :: sf(:,:)
   real(r8), intent(in), optional :: balb

   real(r8), intent(in) :: forc_pbot_g       ! forcing-grid atmospheric pressure [Pa]
   real(r8), intent(in) :: forc_swrad_g      ! forcing-grid downward shortwave (W/m**2)
   real(r8), intent(in) :: forc_pbot_c       ! column atmospheric pressure [Pa]
   real(r8), intent(out):: forc_swrad_c      ! column downward shortwave (W/m**2)

   integer  :: i                             ! Terrain-type index
   integer  :: idx_azi                       ! Solar-azimuth sector index

   real(r8) :: cloud                         ! Cloudiness index
   real(r8) :: difrat                        ! Diffuse fraction of incoming shortwave radiation
   real(r8) :: beam_swrad_g                  ! Direct shortwave radiation at the forcing grid
   real(r8) :: dif_swrad_g                   ! Diffuse shortwave radiation at the forcing grid
   real(r8) :: beam_swrad_c                  ! Area-weighted direct shortwave radiation for the column
   real(r8) :: dif_swrad_c                   ! Diffuse shortwave radiation for the column
   real(r8) :: refl_swrad_c                  ! Terrain-reflected shortwave radiation for the column
   real(r8) :: toa_swrad                     ! Top-of-atmosphere incident shortwave radiation
   real(r8) :: clr_idx                       ! Atmospheric clearness index
   real(r8) :: k_c                           ! Broadband atmospheric attenuation coefficient
   real(r8) :: opt_factor                    ! Optical-path correction factor
   real(r8) :: rt_R                          ! Earth-Sun distance correction factor
   real(r8) :: zen_rad                       ! Solar zenith angle (rad)
   real(r8) :: zenith_segment                ! Zenith-angle threshold for the shadow function (rad)
   real(r8) :: a1                            ! First fitted parameter of the shadow function
   real(r8) :: a2                            ! Second fitted parameter of the shadow function
   real(r8) :: swrad_bd                      ! Lower bound for the shortwave scaling factor
   real(r8) :: cosill(size(slp))             ! Illumination correction for each terrain type
   real(r8) :: beam_swrad_c_sub(size(slp))   ! Direct shortwave contribution from each terrain type
   real(r8) :: tcf(size(slp))                ! Terrain configuration factor for each terrain type
   real(r8) :: refl_swrad_c_sub(size(slp))   ! Reflected shortwave contribution from each terrain type
!--------------------------------------------------------------------------

      ! Participate diffuse and beam radiation
      IF (coszen .eq. 0)THEN
         cloud = 0.
      ELSE
         cloud = (1160.*coszen-forc_swrad_g)/(963.*coszen)
      ENDIF
      cloud = max(cloud,0.)
      cloud = min(cloud,1.)
      cloud = max(0.58,cloud)

      difrat = 0.0604/(coszen-0.0223)+0.0683
      IF(difrat.lt.0.) difrat = 0.
      IF(difrat.gt.1.) difrat = 1.

      difrat = difrat+(1.0-difrat)*cloud

      diff_swrad_g = forc_swrad_g * difrat
      beam_swrad_g = forc_swrad_g * (1-difrat)

      ! Cal top-of-atmosphere shortwave radiation, accounting for
      ! the seasonal variation in Earth-Sun distance and solar zenith angle.
      rt_R = 1 - 0.01672 * cos(0.9856*(julian_day-4)*PI/180._r8)
      toa_swrad = S*(rt_R**2)*coszen

      ! Derive the broadband attenuation coefficient k_c [Pa-1] from the
      ! ratio of surface to top-of-atmosphere shortwave radiation. Following
      ! Rouf et al. (2020), pressure represents the atmospheric optical path.
      clr_idx = max(min(forc_swrad_g/toa_swrad, 1), 0)

      IF (clr_idx.le.0) THEN
         k_c = 0
      ELSE
         k_c = log(clr_idx)/forc_pbot_c
      ENDIF
      
      ! Account for the optical-path difference caused by the pressure
      ! difference between the forcing grid and target column.
      opt_factor = exp(k_c * (forc_pbot_g - forc_pbot_c))

      ! Adjust beam radiation according to illustration angle and optical-path difference
      zen_rad = min(acos(coszen), 85*PI/180.)
      azi_rad = acos(cosazi)

      DO i = 1, size(slp)
         IF (slp(i) == spval) THEN
            beam_swrad_c_sub(i) = spval
         ELSE
            IF (asp(i) >= 0._r8) THEN 
               cosill(i) = cos(atan(slp(i)))+tan(zen_rad)*sin(atan(slp(i)))*cos(azi_rad-asp(i))
               cosill(i) = max(cosill(i), 0)
            ELSE 
               cosill(i) = 1.
            ENDIF

            beam_swrad_c_sub(i) = beam_swrad_g * cosill(i) * opt_factor * area(i)
         ENDIF
      ENDDO

      beam_swrad_c = sum(beam_swrad_c_sub, mask=beam_swrad_c_sub/=spval)
      diff_swrad_c = diff_swrad_g
      refl_swrad_c = 0

      ! Adjust shortwave radiation using other factors in complex mode
      IF (DEF_DS_MODE == 'C') THEN

         ! Cal shadow factor using shadow factor parameters
         ! sf = exp(-1 * exp(a1*zenith + a2))
         idx_azi = int((azi_rad*180/PI)/(360/16))+1
         zen_segment = sf(idx_azi, 1)      
         a1 = sf(idx_azi, 2)                   
         a2 = sf(idx_azi, 3)                

         IF (zen_rad <= zen_segment) THEN
            sf_c = 1.
         ELSEIF (a1<=1e-10) THEN
            sf_c = 1.
         ELSE
            sf_c = exp(-1*exp(min(a1*zen_rad+a2,3.5)))
         ENDIF

         ! Adjust beam radiation using shadow factor
         beam_swrad_c = beam_swrad_c * sf_c

         ! Adjust diffuse radiation using sky view factor
         diff_swrad_c = diff_swrad_g * svf

         ! Cal reflect radiation by surronding terrains by configuration factor
         DO i = 1, size(slp)
            tcf(i) = max((1+cos(slp(i)))/2-svf, 0)

            IF (slp(i) == spval) THEN
               refl_swrad_c_sub(i) = spval
            ELSE
               refl_swrad_c_sub(i) = balb * tcf(i) * (beam_swrad_c*coszen+(1-svf)*diff_swrad_c)
            ENDIF
         ENDDO
         refl_swrad_c = sum(refl_swrad_c_sub, mask=refl_swrad_c_sub/=spval)
      ENDIF

      ! bound the downscaled shortwave radiation to ensure that the denominator 
      ! is not 0 during shortwave normalization
      forc_swrad_c = beam_swrad_c + diff_swrad_c + refl_swrad_c
      forc_swrad_c = min(forc_swrad_c, forc_swrad_g * (1._r8 + swrad_bd)) 
      forc_swrad_c = max(forc_swrad_c, forc_swrad_g * (1._r8 - swrad_bd))
      IF (forc_swrad_c < 1.e-4) forc_swrad_c = 0.0001

   END SUBROUTINE downscale_swrad

!-----------------------------------------------------------------------------

   SUBROUTINE downscale_wind ( &
         slp_c, asp_c, cur_c, area_c, &
         forc_us_g, forc_vs_g,        &
         forc_us_c, forc_vs_c)

!-----------------------------------------------------------------------------
! DESCRIPTION:
!  Adjust wind speed using slope along the wind direction and curvature.
!  Average terrain classes by area and retain the original vector direction.
!
! REFERENCE:
!  Liston, G. E. and Elder, K. (2006)
!  A meteorological distribution system for high-resolution terrestrial modeling (MicroMet)
!  Journal of Hydrometeorology, 7, 217-234.
!-----------------------------------------------------------------------------

   IMPLICIT NONE

   real(r8), intent(in)  :: forc_us_g, forc_vs_g                 ! grid eastward/northward wind [m/s]
   real(r8), intent(in)  :: slp_c(:)                             ! slope angle [rad]; -1.e36: missing
   real(r8), intent(in)  :: asp_c(:)                             ! clockwise from north [rad]; -9999: flat
   real(r8), intent(in)  :: cur_c                                ! normalized terrain curvature [-]
   real(r8), intent(in)  :: area_c(:)                            ! terrain-class area fraction [-]
   real(r8), intent(out) :: forc_us_c, forc_vs_c                 ! column eastward/northward wind [m/s]

   integer               :: i
   real(r8)              :: ws_dir, slp_ws_dir, scale_factor, weighted_scale
!-----------------------------------------------------------------------------

      ! Cal wind speed and direction and convert angles to match the aspect-data convention.
      ! Before: east = 0, [-PI, PI].
      !
      !                         N (PI/2) 
      !                            |
      !             W (+/-PI) -----+----- E (0)
      !                            |
      !                         S (-PI/2)
      !
      ! After: north = 0, [0, 2*PI).
      !
      !                         N (0)
      !                            |
      !             W (3*PI/2) ----+----- E (PI/2)
      !                            |
      !                         S (PI)
      !
      ws_dir = modulo(PI/2.0_r8 - atan2(forc_vs_g, forc_us_g), 2.0_r8*PI)
      ws_g = sqrt(forc_us_g**2 + forc_vs_g**2)

      !
      DO i = 1, size(slp_c)
         IF (slp_c(i) == spval) THEN  ! no this terrain type
            ws_c_sub(i) = spval
         ELSE
            ws_dir_slp = slp_c(i) * cos(ws_dir - asp_c(i))
            IF (asp_c(i) == -9999) THEN  ! flatten area
               scale_factor = 1.0
            ELSE
               scale_factor = min(1.5_r8, abs(1.0_r8 + 0.58_r8*slp_ws_dir + 0.42_r8*cur_c))  ! We let scale factor larger than 0
            ENDIF
            ws_c_sub(i) = ws_g * scale_factor * area_c(i)
         ENDIF
      ENDDO
      ws_c = sum(ws_c_sub(:), mask=ws_c_sub(:)/=spval)

      ! NOTE: we use the wind direction after coverting to the aspect-data convention to calculate the eastward/northward wind components.
      forc_us_c = ws_c * sin(ws_dir)
      forc_vs_c = ws_c * cos(ws_dir)

   END SUBROUTINE downscale_wind

!-----------------------------------------------------------------------------
END MODULE MOD_ForcingDownscaling