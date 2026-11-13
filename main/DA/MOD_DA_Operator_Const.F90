#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Operator_Const
!-----------------------------------------------------------------------
! DESCRIPTION:
!   1. Define constants (do not rely on satellite parameters) used in EM.
!   2. Define parameters (IGBP) used in EM.
!
! REFERENCES:
!   [1] Wigneron, J.-P., Kerr, Y. H., Waldteufel, P., Saleh, K.,
!       Escorihuela, M.-J., Richaume, P., Ferrazzoli, P., de Rosnay, P.,
!       Gurney, R., Calvet, J.-C., Grant, J. P., Guglielmetti, M.,
!       Hornbuckle, B., Mätzler, C., Pellarin, T., & Schwank, M. (2007).
!       L-band microwave emission of the biosphere (L-MEB) model:
!       Description and calibration against experimental data sets over
!       crop fields. Remote Sensing of Environment, 107(4), 639-655.
!       https://doi.org/10.1016/j.rse.2006.10.014
!
!   [2] Wigneron, J.-P., Jackson, T. J., O'Neill, P., De Lannoy, G., de
!       Rosnay, P., Walker, J. P., Ferrazzoli, P., Mironov, V., Bircher,
!       S., Grant, J. P., Kurum, M., Schwank, M., Muñoz-Sabater, J.,
!       Das, N., Royer, A., Al-Yaari, A., Al Bitar, A., Fernandez-Moran,
!       R., Lawrence, H., ... Kerr, Y. (2017). Modelling the passive
!       microwave signature from land surfaces: A review of recent
!       results and application to the L-band SMOS & SMAP soil moisture
!       retrieval algorithms. Remote Sensing of Environment, 192,
!       238-262.
!       https://doi.org/10.1016/j.rse.2017.01.024
!
! HISTORY:
!   Lu Li, 12/2024: Initial version
!   Lu Li, Zhilong Fan, 11/2026: DA v1.0
!-----------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Vars_Global, only: pi, N_land_classification
   IMPLICIT NONE
   SAVE

   ! Constant variables
   real(r8), parameter    :: C = 2.998e8           ! Speed of light (m/s)
   real(r8), parameter    :: mu0 = 4.*pi*1e-7      ! Vacuum permeability (H/m)
   real(r8), parameter    :: eps_0 = 8.854e-12     ! Dielectric constant of free space (Klein and Swift 1977) [Farads/meter]
   real(r8), parameter    :: z0 = sqrt(mu0/eps_0)  ! Impedance of free space (Ohm)
   real(r8), parameter    :: eps_w_inf = 4.9       ! Dielectric constant at infinite frequency (Stogryn 1971),
   real(r8), parameter    :: rho_soil = 2.66       ! Soil specific density (g/cm3)
   real(r8), parameter    :: f0w = 9.              ! Relaxation frequency of liquid water (GHz)
   real(r8), parameter    :: rgh_surf = 2.2        ! Soil surface roughness (cm)
   complex(r8), parameter :: jj = (0., 1.)         ! Imaginary unit for complex number

#ifdef LULC_IGBP
   ! 0  Ocean
   ! 1  Evergreen Needleleaf Forests
   ! 2  Evergreen Broadleaf Forests
   ! 3  Deciduous Needleleaf Forests
   ! 4  Deciduous Broadleaf Forests
   ! 5  Mixed Forests
   ! 6  Closed Shrublands
   ! 7  Open Shrublands
   ! 8  Woody Savannas
   ! 9  Savannas
   !10  Grasslands
   !11  Permanent Wetlands
   !12  Croplands
   !13  Urban and Built-up Lands
   !14  Cropland/Natural Vegetation Mosaics
   !15  Permanent Snow and Ice
   !16  Barren
   !17  Water Bodies

   ! Empirical parameters to account for the dependence of optical depth on incidence angle [1]
   real(r8), parameter, dimension(N_land_classification) :: tth &
      = (/0.80, 1.00, 0.80, 0.49, 0.49, &
      1.00, 1.00, 1.00, 1.00, 1.00,  &
      1.00, 1.00, 1.00, 1.00, 1.00,  &
      1.00, 2.00/)
   real(r8), parameter, dimension(N_land_classification) :: ttv &
      = (/0.80, 1.00, 0.80, 0.46, 0.46, &
      1.00, 1.00, 1.00, 1.00, 1.00,  &
      1.00, 2.00, 1.00, 2.00, 1.00,  &
      1.00, 1.00/)

   ! Empirical roughness parameters (Table 2 in [2])
   real(r8), parameter, dimension(N_land_classification) :: hr_SMAP &
      = (/0.160, 0.160, 0.160, 0.160, 0.160, &
      0.110, 0.110, 0.125, 0.156, 0.156, &
      0.100, 0.108, 0.000, 0.130, 0.000, &
      0.150, 0.000/)

   real(r8), parameter, dimension(N_land_classification) :: hr_SMOS &
      = (/0.300, 0.300, 0.300, 0.300, 0.300, &
      0.100, 0.100, 0.100, 0.100, 0.100, &
      0.100, 0.100, 0.100, 0.100, 0.000, &
      0.100, 0.000/)

   real(r8), parameter, dimension(N_land_classification) :: hr_P16  &
      = (/0.350, 0.460, 0.430, 0.450, 0.410, &
      0.260, 0.170, 0.350, 0.230, 0.130, &
      0.020, 0.170, 0.190, 0.220, 0.000, &
      0.020, 0.000/)

   ! b parameters for Wigneron vegetation model (ref?)
   real(r8), parameter, dimension(N_land_classification) :: b1 &
      = (/0.2600, 0.2260, 0.2600, 0.2260, 0.2260, &
      0.0375, 0.0375, 0.0375, 0.0375, 0.0375, &
      0.0000, 0.0500, 0.0000, 0.0500, 0.0000, &
      0.0000, 0.0500/)

   real(r8), parameter, dimension(N_land_classification) :: b2 &
      = (/0.0060, 0.0010, 0.0060, 0.0010, 0.0010, &
      0.0500, 0.0500, 0.0500, 0.0500, 0.0500, &
      0.0000, 0.0000, 0.0000, 0.0000, 0.0000, &
      0.0000, 0.0000/)

   real(r8), parameter, dimension(N_land_classification) :: b3 &
      = (/0.6900, 0.7000, 0.6900, 0.7000, 0.7000, &
      0.0000, 0.0000, 0.0000, 0.0000, 0.0000, &
      0.0000, 0.0000, 0.0000, 0.0000, 0.0000, &
      0.0000, 0.0000/)

   ! Effective diffusion albedo (Table 3 in [2])
   real(r8), parameter, dimension(N_land_classification) :: w_SMAPL2 &
      = (/0.050, 0.050, 0.050, 0.050, 0.050, &
      0.050, 0.050, 0.050, 0.080, 0.050, &
      0.050, 0.000, 0.065, 0.000, 0.000, &
      0.000, 0.000/)

   real(r8), parameter, dimension(N_land_classification) :: w_CMEM &
      = (/0.080, 0.095, 0.080, 0.070, 0.070, &
      0.050, 0.050, 0.050, 0.050, 0.050, &
      0.000, 0.000, 0.000, 0.000, 0.000, &
      0.000, 0.000/)

   real(r8), parameter, dimension(N_land_classification) :: w_K16 &
      = (/0.050, 0.050, 0.060, 0.030, 0.050, &
      0.030, 0.050, 0.040, 0.020, 0.030, &
      0.000, 0.040, 0.000, 0.020, 0.000, &
      0.000, 0.000/)

   real(r8), parameter, dimension(N_land_classification) :: w_SMAPL4 &
      = (/0.120, 0.080, 0.120, 0.100, 0.120, &
      0.140, 0.110, 0.130, 0.120, 0.070, &
      0.000, 0.120, 0.000, 0.150, 0.000, &
      0.000, 0.000/)
#endif

END MODULE MOD_DA_Operator_Const

!-----------------------------------------------------------------------------
#endif
