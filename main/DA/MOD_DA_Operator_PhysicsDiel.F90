#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Operator_PhysicsDiel
!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Dielectric properties used by the microwave observation operators.
!
!   Provides liquid-water and ice permittivity and the soil dielectric mixing
!   schemes W80, D85, M04, and M09. These functions are shared by the soil,
!   snow, and other electromagnetic calculations.
!
! REFERENCES:
!   [1] Ulaby, F. T., Moore, R. K., & Fung, A. K. (1986). Microwave Remote
!       Sensing: Active and Passive, Vol. III: From Theory to Applications.
!       Artech House, Norwood, MA.
!
!   [2] Klein, L. A., & Swift, C. T. (1977). An improved model for the
!       dielectric constant of sea water at microwave frequencies. IEEE
!       Transactions on Antennas and Propagation, 25(1), 104-111.
!       DOI: 10.1109/TAP.1977.1141539
!
!   【3】Mätzler, C. (2006). Thermal Microwave Radiation: Applications for
!       Remote Sensing. IET Electromagnetic Waves Series 52, pp. 456-461.
!       DOI: 10.1049/PBEW052E
!   【4】Wang, J. R., & Schmugge, T. J. (1980). An empirical model for the
!       complex dielectric permittivity of soils as a function of water
!       content. IEEE Transactions on Geoscience and Remote Sensing,
!       GE-18(4), 288-295. DOI: 10.1109/TGRS.1980.350304
!   【5】Dobson, M. C., Ulaby, F. T., Hallikainen, M. T., & El-Rayes, M. A.
!       (1985). Microwave dielectric behavior of wet soil - Part II:
!       Dielectric mixing models. IEEE Transactions on Geoscience and Remote
!       Sensing, GE-23(1), 35-46. DOI: 10.1109/TGRS.1985.289498
!   【6】Peplinski, N. R., Ulaby, F. T., & Dobson, M. C. (1995). Dielectric
!       properties of soils in the 0.3-1.3-GHz range. IEEE Transactions on
!       Geoscience and Remote Sensing, 33(3), 803-807.
!       DOI: 10.1109/36.387598
!   【7】Peplinski, N. R., Ulaby, F. T., & Dobson, M. C. (1995). Corrections
!       to "Dielectric properties of soils in the 0.3-1.3-GHz range". IEEE
!       Transactions on Geoscience and Remote Sensing, 33(6), 1340.
!       DOI: 10.1109/TGRS.1995.477193
!   【8】Mironov, V. L., Dobson, M. C., Kaupp, V. H., Komarov, S. A., &
!       Kleshchenko, V. N. (2004). Generalized refractive mixing dielectric
!       model for moist soils. IEEE Transactions on Geoscience and Remote
!       Sensing, 42(4), 773-785. DOI: 10.1109/TGRS.2003.823288
!   【9】Mironov, V. L., Kosolapova, L. G., & Fomin, S. V. (2009). Physically
!       and mineralogically based spectroscopic dielectric model for moist
!       soils. IEEE Transactions on Geoscience and Remote Sensing, 47(7),
!       2059-2070. DOI: 10.1109/TGRS.2008.2011631
!   【10】Mironov, V. L., & Fomin, S. V. (2009). Temperature and mineralogy
!       dependable model for microwave dielectric spectra of moist soils.
!       PIERS Online, 5(5), 411-415. DOI: 10.2529/PIERS090220054025
!
! HISTORY:
!   Lu Li, 12/2024: Initial version
!   Zhilong Fan, Lu Li, 03/2024: Debug and clean codes
!   Lu Li, 10/2025: Debug and clean codes
!   Lu Li, Zhilong Fan, 11/2026: DA v1.0
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Const_Physical
   USE MOD_DA_Operator_Const
   IMPLICIT NONE
   SAVE

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE diel_ice(fghz, t, eps_i)

!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the dielectric constant of pure ice.
!-----------------------------------------------------------------------------
!------------------------------ Dummy Arguments ------------------------------
      real(r8), intent(in)     :: fghz                        ! Microwave frequency [GHz]
      real(r8), intent(in)     :: t                           ! Temperature (C)
      complex(r8), intent(out) :: eps_i                       ! Dielectric constant of ice water

!------------------------------ Local Variables ------------------------------
      real(r8)                 :: betam                       ! Beta parameter by Mishima et al. (1983)
      real(r8)                 :: dbeta                       ! Corrected delta beta parameter
      real(r8)                 :: beta                        ! Beta parameter
      real(r8)                 :: t_inv                       ! Modified inverse temperature
      real(r8)                 :: tk                          ! Temperature (K)
      real(r8)                 :: alpha                       ! Alpha parameter
      real(r8)                 :: eps_i_r                     ! Real part of pure ice dielectric constant
      real(r8)                 :: eps_i_i                     ! Imaginary part of pure ice dielectric constant
!-----------------------------------------------------------------------------
      ! Convert the temperature from Celsius to Kelvin
      tk = t + tfrz

      ! Calculate the beta parameter by Mishima et al. (1983)
      betam = (0.0207/tk)*(exp(335./tk)/((exp(335./tk) - 1.)**2.)) + 1.16e-11*(fghz**2.) ! 【3】Eq. (5.33)

      ! Calculate the corrected delta beta parameter
      dbeta = exp(-10.02 + 0.0364*t) ! 【3】Eq. (5.35)

      ! Calculate the beta parameter
      beta = betam + dbeta ! 【3】Eq. (5.34)

      ! Calculate the alpha parameter
      t_inv = 300./tk - 1                                ! 【3】p. 457
      alpha = (0.00504 + 0.0062*t_inv)*exp(-22.1*t_inv)  ! 【3】Eq. (5.32)

      ! Calculate the real part of the pure ice dielectric constant
      eps_i_r = 3.1884 + 9.1e-4*t ! 【3】Eq. (5.30)

      ! Calculate the imaginary part of the pure ice dielectric constant
      eps_i_i = alpha/fghz + beta*fghz ! 【3】Eq. (5.31)

      ! Calculate the dielectric constant of pure ice
      eps_i = eps_i_r - jj*eps_i_i ! 【3】Eq. (5.31)

   END SUBROUTINE diel_ice

!-----------------------------------------------------------------------------

   SUBROUTINE diel_water(fghz, type, swc, t, wf_sand, wf_clay, BD_all, sal, eps_w)

!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the dielectric constant of pure water, sea water, or soil water.
!-----------------------------------------------------------------------------
!------------------------------ Dummy Arguments ------------------------------
      real(r8), intent(in)     :: fghz                        ! Microwave frequency [GHz]
      integer, intent(in)      :: type                        ! Type of water, 0: pure water, 1: sea water, 2: soil water
      real(r8), intent(in)     :: swc                         ! Soil water content (m3/m3)
      real(r8), intent(in)     :: t                           ! Soil temperature (C)
      real(r8), intent(in)     :: wf_sand                     ! Gravimetric sand percent fraction(%)
      real(r8), intent(in)     :: wf_clay                     ! Gravimetric clay percent fraction(%)
      real(r8), intent(in)     :: BD_all                      ! Bulk density(g/cm3)
      real(r8), intent(in)     :: sal                         ! Water salinity (psu)
      complex(r8), intent(out) :: eps_w                       ! Dielectric constant of soil water

!------------------------------ Local Variables ------------------------------
      real(r8)                 :: f                           ! Microwave frequency [Hz]
      real(r8)                 :: omega                       ! Angular frequency [rad/s]
      real(r8)                 :: sigma                       ! Ionic conductivity (S/m)
      real(r8)                 :: a, b                        ! Parameters
      real(r8)                 :: tau_w                       ! Relaxation time of pure water
      real(r8)                 :: eps_w0                      ! Static dielectric constant of pure water 【2】Eq. (14)
      real(r8)                 :: wc
!-----------------------------------------------------------------------------
      ! 【2】Eq. (16): tau(T, sal) = tau_w(T)*b(sal, T)
      ! Calculate the relaxation time of pure water
      f = fghz*1e9
      omega = 2.0*pi*f

      tau_w = 1.768e-11 - 6.068e-13*t + 1.104e-14*t**2 - 8.111e-17*t**3              ! 【2】Eq. (17)
      b = 1.000 + 2.282e-5*sal*t - 7.638e-4*sal - 7.760e-6*sal**2 + 1.105e-8*sal**3  ! 【2】Eq. (18)
      tau_w = tau_w*b                                                                ! 【2】Eq. (16)

      ! 【2】Eq. (13): eps_w0(sal, T) = eps_w0(T)*a(sal, T)
      ! Static dielectric constant of pure water 【2】Eq. (14)
      eps_w0 = 87.134 - 1.949e-1*t - 1.276e-2*t**2 + 2.491e-4*t**3                   ! 【2】Eq. (14)
      a = 1.000 + 1.613e-5*sal*t - 3.656e-3*sal + 3.210e-5*sal**2 - 4.232e-7*sal**3  ! 【2】Eq. (15)
      eps_w0 = eps_w0*a                                                              ! 【2】Eq. (13)

      IF (type == 0) THEN ! Pure water
         ! 【1】Eq. (19)
         eps_w0 = 88.045 - 0.4147*t + 6.295e-4*t**2 + 1.075e-5*t**3
         eps_w = eps_w_inf + (eps_w0 - eps_w_inf)/(1 - jj*omega*tau_w)

      ELSEIF (type == 1) THEN ! Sea water
         ! Calculate the ionic conductivity 【1】Eqs. (27)-(28)
         sigma = sal*(0.182521 - 1.46192e-3*sal + 2.09324e-5*sal**2 - 1.28205e-7*sal**3) &
            *exp(-1.*(25 - t)*                                                           &
            (2.033e-2 + 1.266e-4*(25 - t) + 2.464e-6*(25 - t)**2                         &
            - sal*(1.849e-5 - 2.551e-7*(25 - t) + 2.551e-8*(25 - t)**2)))

         ! Dielectric constant of sea water 【1】Eq. (21)
         eps_w = eps_w_inf + (eps_w0 - eps_w_inf)/(1 - jj*omega*tau_w) + jj*sigma/(omega*eps_0)
      ELSEIF (type == 2) THEN
         ! Calculate the soil conductivity 【1】
         sigma = -1.645 + 1.939*BD_all - 0.02256*wf_sand + 0.01594*wf_clay
         IF (sigma < 0.) THEN
            sigma = 0. ! Negative for very sandy soils with a low bulk density
         END IF

         ! Calculate the dielectric constant of soil water with a modified
         ! Debye expression 【2】Eq. (16)
         wc = max(0.001, swc)
         eps_w = eps_w_inf + (eps_w0 - eps_w_inf)/(1 - jj*omega*tau_w) &
            + jj*sigma/(omega*eps_0)*(rho_soil - BD_all)/(rho_soil*wc)
      END IF

   END SUBROUTINE diel_water

!-----------------------------------------------------------------------------

   SUBROUTINE diel_soil_W80(ew, t, wc, wf_sand, wf_clay, porsl, eps)

!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the dielectric constant of wet soil with the Wang and Schmugge
!   (1980) empirical mixing model.
!-----------------------------------------------------------------------------
!------------------------------ Dummy Arguments ------------------------------
      real(r8), intent(in)     :: t                           ! Soil temperature (C)
      real(r8), intent(in)     :: wc                          ! Volumetric soil moisture (m3/m3)
      real(r8), intent(in)     :: wf_sand                     ! Gravimetric sand percent fraction(%)
      real(r8), intent(in)     :: wf_clay                     ! Gravimetric clay percent fraction(%)
      real(r8), intent(in)     :: porsl                       ! Soil porosity at surface
      complex(r8), intent(in)  :: ew                          ! Dielectric constant of water
      complex(r8), intent(out) :: eps                         ! Dielectric constant

!------------------------------ Local Variables ------------------------------
      real(r8)                 :: wp                          ! Wilting point
      real(r8)                 :: wt                          ! Transition moisture point (cm3/cm3)
      real(r8)                 :: gamma                       ! Fitting parameter
      real(r8)                 :: ecl                         ! Conductivity loss
      real(r8)                 :: alpha                       ! Conductivity loss parameter
      real(r8)                 :: sal_sea = 32.5              ! Sea water salinity (psu)
      real(r8)                 :: sal_soil = 0.0              ! Soil salinity (psu)
      complex(r8)              :: eps_x                       ! Dielectric constant of the initially absorbed water
      complex(r8)              :: eps_a = (1.0, 0.0)          ! Dielectric constant of air 【4】Sec. IV
      complex(r8)              :: eps_r = (5.5, 0.2)          ! Dielectric constant of rock 【4】Sec. IV
      complex(r8)              :: eps_i = (3.2, 0.1)          ! Dielectric constant of ice 【4】Sec. IV
      complex(r8)              :: eps_f = (5.0, 0.5)          ! Dielectric constant of frozen soil
!-----------------------------------------------------------------------------
      ! Calculate the wilting point
      wp = 0.06774 - 0.00064*wf_sand + 0.00478*wf_clay ! 【4】Eq. (1)

      ! Calculate the fitting parameter
      gamma = -0.57*wp + 0.481 ! 【4】Eq. (8)

      ! Calculate the transition moisture
      wt = 0.49*wp + 0.165 ! 【4】Eq. (9)

      ! Calculate the dielectric constant of wet soil; when all soil is
      ! frozen, eps_x = eps_i
      IF (wc <= wt) THEN
         eps_x = eps_i + (ew - eps_i)*(wc/wt)*gamma              ! 【4】Eq. (3)
         eps = wc*eps_x + (porsl - wc)*eps_a + (1.-porsl)*eps_r  ! 【4】Eq. (2)
      ELSE
         eps_x = eps_i + (ew - eps_i)*gamma                                     ! 【4】Eq. (5)
         eps = wt*eps_x + (wc - wt)*ew + (porsl - wc)*eps_a + (1.-porsl)*eps_r  ! 【4】Eq. (4)
      END IF

      ! Add the conductivity loss to the imaginary part
      alpha = min(100.*wp, 26.)
      ecl = alpha*wc**2                 ! 【4】Eq. (6)
      eps = eps + jj*ecl                ! 【4】Eq. (6)

   END SUBROUTINE diel_soil_W80

!-----------------------------------------------------------------------------

   SUBROUTINE diel_soil_D85(ew, swc, wf_sand, wf_clay, BD_all, eps)

!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the dielectric constant of wet soil with the Dobson et al.
!   (1985) mixing model.
!-----------------------------------------------------------------------------
!------------------------------ Dummy Arguments ------------------------------
      complex(r8), intent(in)  :: ew                          ! Dielectric constant of water
      real(r8), intent(in)     :: swc                         ! Soil moisture
      real(r8), intent(in)     :: wf_sand                     ! Gravimetric sand percent fraction(%)
      real(r8), intent(in)     :: wf_clay                     ! Gravimetric clay percent fraction(%)
      real(r8), intent(in)     :: BD_all                      ! Soil bulk density (g/cm3)
      complex(r8), intent(out) :: eps                         ! Dielectric constant

!------------------------------ Local Variables ------------------------------
      real(r8)                 :: alphas = 0.65_r8
      real(r8)                 :: beta, eaa, eps_s, epsi, epsr, wc
!-----------------------------------------------------------------------------
      wc    = max(swc, 0.001_r8)
      eps_s = (1.01_r8 + 0.44_r8 * rho_soil)**2.0_r8 - 0.062_r8                 ! 【5】Eq. (22)
      beta  = (127.48_r8 - 0.519_r8 * wf_sand - 0.152_r8 * wf_clay) / 100.0_r8  ! 【5】Eq. (30)
      eaa   = 1.0_r8 + (BD_all / rho_soil) * (eps_s ** alphas - 1.0_r8) &
         + (wc ** beta) * (real(ew) ** alphas) - wc                              ! 【5】Eq. (28)
      epsr  = eaa ** (1.0_r8/alphas)                                             ! 【6】Eq. (2), 【7】
      beta  = (133.797_r8 - 0.603_r8 * wf_sand - 0.166_r8 * wf_clay) / 100.0_r8  ! 【5】Eq. (31) 1.33797 -> 133.797, 【6】Eq. (5)
      eaa   = (wc ** beta) * (abs(aimag(ew)) ** alphas)                          ! 【6】Eq. (3)
      epsi  = eaa ** (1.0_r8/alphas)                                             ! 【6】Eq. (3)
      eps   = cmplx(epsr, epsi, kind=r8)

   END SUBROUTINE diel_soil_D85

!-----------------------------------------------------------------------------

   SUBROUTINE diel_soil_M04(fghz, wc, wf_clay, eps)

!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the dielectric constant of wet soil with the generalized
!   refractive mixing dielectric model (GRMDM) of Mironov et al. (2004) 【8】.
!-----------------------------------------------------------------------------
!------------------------------ Dummy Arguments ------------------------------
      real(r8), intent(in)     :: fghz                        ! Microwave frequency [GHz]
      real(r8), intent(in)     :: wc                          ! Soil moisture (m3/m3)
      real(r8), intent(in)     :: wf_clay                     ! Gravimetric clay percent fraction(%)
      complex(r8), intent(out) :: eps                         ! Dielectric constant

!------------------------------ Local Variables ------------------------------
      real(r8)                 :: f                           ! Microwave frequency [Hz]
      real(r8)                 :: znd, zkd, zxmvt, zep0b, ztaub, zsigmab, zep0u, ztauu
      real(r8)                 :: zsigmau, zcxb, zepwbx, zepwby, zcxu, zepwux, zepwuy
      real(r8)                 :: znb, zkb, znu, zku, zxmvt2, znm, zkm, zepmx, zepmy
      integer                  :: zflag
!-----------------------------------------------------------------------------
! Initialize the GRMDM spectroscopic parameters from the clay fraction
      ! Refractive index and normalized attenuation coefficient of dry soil
      f = fghz*1e9

      znd = 1.634 - 0.539 * (wf_clay/100) + 0.2748 * (wf_clay/100) ** 2
      zkd = 0.03952 - 0.04038 * (wf_clay / 100) ! 【9】Eq. (18)

      ! Maximum bound water fraction
      zxmvt = 0.02863 + 0.30673 * wf_clay / 100 ! 【9】Eq. (19)

      ! Bound water parameters
      zep0b   = 79.8 - 85.4  * (wf_clay / 100) + 32.7  * (wf_clay / 100)*(wf_clay / 100)  ! 【9】Eq. (20)
      ztaub   = 1.062e-11 + 3.450e-12 * (wf_clay / 100)                                   ! 【9】Eq. (21)
      zsigmab = 0.3112 + 0.467 * (wf_clay / 100)                                          ! 【9】Eq. (22)

      ! Unbound (free) water parameters
      zep0u   = 100                     ! 【9】Eq. (24)
      ztauu   = 8.5e-12                 ! 【9】Eq. (25)
      zsigmau = 0.3631 + 1.217 * (wf_clay / 100)

      ! Compute the dielectric constant of bound and unbound water
      zcxb   = (zep0b - eps_w_inf) / (1. + (2.*pi*f*ztaub)**2)     ! 【9】Eq. (16)
      zepwbx = eps_w_inf + zcxb                                    ! 【9】Eq. (16)
      zepwby = zcxb * (2.*pi*f*ztaub) + zsigmab / (2.*pi*eps_0*f)  ! 【9】Eq. (16)
      zcxu   = (zep0u - eps_w_inf) / (1 + (2*pi*f*ztauu)**2)       ! 【9】Eq. (16)
      zepwux = eps_w_inf + zcxu                                    ! 【9】Eq. (16)
      zepwuy = zcxu * (2.*pi*f*ztauu) + zsigmau/(2.*pi*eps_0*f)

      ! Compute the refractive index of bound and unbound water
      znb = sqrt( sqrt( zepwbx**2 + zepwby**2) + zepwbx ) / sqrt(2.0)  ! 【9】Eq. (14)
      zkb = sqrt( sqrt( zepwbx**2 + zepwby**2) - zepwbx ) / sqrt(2.0)  ! 【9】Eq. (15)
      znu = sqrt( sqrt( zepwux**2 + zepwuy**2) + zepwux ) / sqrt(2.0)  ! 【9】Eq. (14)
      zku = sqrt( sqrt( zepwux**2 + zepwuy**2) - zepwux ) / sqrt(2.0)  ! 【9】Eq. (15)

      ! Compute the soil refractive index
      zxmvt2 = min(wc, zxmvt)
      zflag  = 0
      IF ( wc >= zxmvt ) zflag = 1
      znm = znd + (znb - 1) * zxmvt2 + (znu - 1) * (wc-zxmvt) * zflag  ! 【9】Eq. (12)
      zkm = zkd + zkb * zxmvt2 + zku * (wc-zxmvt) * zflag              ! 【9】Eq. (13)

      ! Compute the soil dielectric constant
      zepmx = znm ** 2 - zkm ** 2       ! 【9】Eq. (11)
      zepmy = znm * zkm * 2             ! 【9】Eq. (11)
      eps   = cmplx(zepmx, zepmy, kind=r8)

   END SUBROUTINE diel_soil_M04

!-----------------------------------------------------------------------------

   SUBROUTINE diel_soil_M09(fghz, wc, t, wf_clay, eps)

!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the dielectric constant of wet soil with the temperature-
!   dependent model of Mironov and Fomin (2009).
!-----------------------------------------------------------------------------
!------------------------------ Dummy Arguments ------------------------------
      real(r8), intent(in)     :: fghz                        ! Microwave frequency [GHz]
      real(r8), intent(in)     :: t                           ! Soil temperature
      real(r8), intent(in)     :: wc                          ! Soil moisture (m3/m3)
      real(r8), intent(in)     :: wf_clay                     ! Weighted fraction (%)
      complex(r8), intent(out) :: eps                         ! Dielectric constant

!------------------------------ Local Variables ------------------------------
      real(r8)                 :: f                           ! Microwave frequency [Hz]
      real(r8)                 :: nd, kd, mvt, ts, e0b, Bb, Bsgb, Fb, eb0
      real(r8)                 :: dSbR, taub, sigmabt, sigmab, e0u, Bu, Bsgu, Fu, eu0, dHuR, dSuR, tauu, dHbR
      real(r8)                 :: sigmau, sigmaut, cxb, eb_r, eb_i, cxu, eu_r, eu_i, nb, kb, nu, ku, nm, km, eps_r, eps_i
!-----------------------------------------------------------------------------
! Initialize the GRMDM spectroscopic parameters from the clay fraction
      ! Refractive index and normalized attenuation coefficient of dry soil
      f = fghz*1e9

      nd = 1.634 - 0.539e-2 * wf_clay + 0.2748e-4 * (wf_clay ** 2)  ! 【10】Eq. (11)
      kd = 0.03952 - 0.04038e-2 * wf_clay                           ! 【10】Eq. (12)

      ! Maximum bound water fraction
      mvt = 0.02863 + 0.30673e-2 * wf_clay ! 【10】Eq. (13)

      ! Starting temperature of the parameter fit (【10】p. 413)
      ts = 20.

      ! Compute the static dielectric constant of bound water
      ! 【10】Eq. (14)
      e0b  = 79.8 - 85.4e-2 * wf_clay + 32.7e-4 * (wf_clay **2)
      ! 【10】Eq. (15)
      Bb   = 8.67e-19 - 0.00126e-2 * wf_clay + 0.00184e-4 * (wf_clay ** 2)  - 9.77e-10*(wf_clay**3) - 1.39e-15 *(wf_clay**4)
      ! 【10】Eq. (23)
      Bsgb = 0.0028  + 0.02094e-2*wf_clay - 0.01229e-4*(wf_clay**2) - 5.03e-22*(wf_clay**3) + 4.163e-24*(wf_clay**4)
      ! 【10】Eq. (8), ep0 -> e0p
      Fb   = log((e0b - 1)/(e0b + 2))
      ! 【10】Eq. (7), e0p -> ep0
      eb0  = (1 + 2*exp(Fb-Bb*(t-ts))) / (1 - exp(Fb-Bb*(t-ts)))

      ! Compute the relaxation time of bound water
      dHbR = 1467 + 2697e-2*wf_clay - 980e-4 *(wf_clay**2) + 1.368e-10*(wf_clay**3) - 8.61e-13 *(wf_clay**4)    ! 【10】Eq. (18)
      dSbR = 0.888 + 9.7e-2 *wf_clay - 4.262e-4*(wf_clay**2) + 6.79e-21 *(wf_clay**3) + 4.263e-22*(wf_clay**4)  ! 【10】Eq. (19)
      taub = 48e-12 * exp(dHbR/(t+tfrz)-dSbR)/(t+tfrz)                                                          ! 【10】Eq. (9)

      ! Compute the conductivity of bound water
      sigmabt = 0.3112 + 0.467e-2*wf_clay  ! 【10】Eq. (22)
      sigmab  = sigmabt + Bsgb*(t-ts)      ! 【10】Eq. (10)

      ! Unbound (free) water parameters
!-----------------------------------------------------------------------------
! Compute the static dielectric constant of unbound water
!-----------------------------------------------------------------------------
      e0u  = 100.                                                                                                    ! 【10】Eq. (16)
      Bu   = 1.11e-4 - 1.603e-7 *wf_clay + 1.239e-9 *(wf_clay**2) + 8.33e-13 *(wf_clay**3) - 1.007e-14*(wf_clay**4)  ! 【10】Eq. (17)
      Bsgu = 0.00108 + 0.1413e-2*wf_clay - 0.2555e-4*(wf_clay**2) + 0.2147e-6*(wf_clay**3) - 0.0711e-8*(wf_clay**4)  ! 【10】Eq. (25)
      ! 【10】Eq. (8), ep0 -> e0p
      Fu   = log((e0u - 1)/(e0u + 2))
      ! 【10】Eq. (7), e0p -> ep0
      eu0  = (1 + 2*exp(Fu-Bu*(t-ts))) / (1-exp(Fu-Bu*(t-ts)))

!-----------------------------------------------------------------------------
! Compute the relaxation time of unbound water
!-----------------------------------------------------------------------------
      dHuR = 2231 - 143.1e-2 *wf_clay + 223.2e-4*(wf_clay**2) - 142.1e-6*(wf_clay**3) + 27.14e-8 *(wf_clay**4)    ! 【10】Eq. (20)
      dSuR = 3.649 - 0.4894e-2*wf_clay + 0.763e-4*(wf_clay**2) - 0.4859e-6*(wf_clay**3) + 0.0928e-8*(wf_clay**4)  ! 【10】Eq. (21)
      tauu = 48e-12 * exp(dHuR/(t+tfrz)-dSuR)/(t+tfrz)                                                            ! 【10】Eq. (9)

!-----------------------------------------------------------------------------
! Compute the conductivity of unbound water
!-----------------------------------------------------------------------------
      sigmaut = 0.05 + 1.4*(1.0 - (1.0 - wf_clay*1.e-2)**4.664)  ! 【10】Eq. (24)
      sigmau  = sigmaut + Bsgu*(t-ts)                            ! 【10】Eq. (10)

!-----------------------------------------------------------------------------
! Compute the dielectric constant of bound and unbound water
!-----------------------------------------------------------------------------
      cxb  = (eb0-eps_w_inf) / (1+(2*pi*f*taub)**2)     ! 【10】Eq. (6), 【9】Eq. (16)
      eb_r = eps_w_inf + cxb                            ! 【10】Eq. (6), 【9】Eq. (16)
      eb_i = cxb*(2*pi*f*taub) + sigmab/(2*pi*eps_0*f)  ! 【10】Eq. (6), 【9】Eq. (16)
      cxu  = (eu0-eps_w_inf) / (1+(2*pi*f*tauu)**2)     ! 【10】Eq. (6), 【9】Eq. (16)
      eu_r = eps_w_inf + cxu                            ! 【10】Eq. (6), 【9】Eq. (16)
      eu_i = cxu*(2*pi*f*tauu) + sigmau/(2*pi*eps_0*f)  ! 【10】Eq. (6), 【9】Eq. (16)

!-----------------------------------------------------------------------------
! Compute the refractive index of bound and unbound water
!-----------------------------------------------------------------------------
      nb = sqrt(sqrt(eb_r**2+eb_i**2)+eb_r) / sqrt(2.0)  ! 【10】Eq. (5)
      kb = sqrt(sqrt(eb_r**2+eb_i**2)-eb_r) / sqrt(2.0)  ! 【10】Eq. (5)
      nu = sqrt(sqrt(eu_r**2+eu_i**2)+eu_r) / sqrt(2.0)  ! 【10】Eq. (5)
      ku = sqrt(sqrt(eu_r**2+eu_i**2)-eu_r) / sqrt(2.0)  ! 【10】Eq. (5)

!-----------------------------------------------------------------------------
! Compute the soil refractive index
!-----------------------------------------------------------------------------
      IF (wc <= mvt) THEN
         nm = nd + (nb-1)*wc            ! 【9】Eq. (12)
         km = kd + kb*wc                ! 【9】Eq. (13)
      ELSE
         nm = nd + (nb-1)*mvt + (nu-1)*(wc-mvt)  ! 【9】Eq. (12)
         km = kd + kb*mvt + ku*(wc-mvt)          ! 【9】Eq. (13)
      ENDIF

!-----------------------------------------------------------------------------
! Compute the soil dielectric constant
!-----------------------------------------------------------------------------
      eps_r = nm**2 - km**2             ! 【10】Eq. (4)
      eps_i = 2* nm * km                ! 【10】Eq. (4)
      eps   = cmplx(eps_r, eps_i, kind=r8)

   END SUBROUTINE diel_soil_M09

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Operator_PhysicsDiel
#endif

