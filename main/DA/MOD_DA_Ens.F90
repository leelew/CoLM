#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Ens
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Provide procedures for ensemble generation and management in data
!   Assimilation, supporting perturbations to atmospheric forcing,
!   Model parameters, and model states.
!
! REFERENCES:
!   [1] Reichle, R. H., Koster, R. D., De Lannoy, G. J. M., Crow, W. T.,
!       And Kimball, J. S. (2014). Algorithm Theoretical Basis Document:
!       Level 4 Surface and Root Zone Soil Moisture (L4_SM) Data Product.
!       Revision A, JPL D-66483.
!       https://nsidc.org/sites/default/files/272_l4_sm_reva_web_0.pdf
!
!   [2] Seo, E., Lee, M.-I., and Reichle, R. H. (2021). Assimilation of
!       SMAP and ASCAT soil moisture retrievals into the JULES land surface
!       Model using the Local Ensemble Transform Kalman Filter.
!       Remote Sensing of Environment, 253, 112222.
!       https://doi.org/10.1016/j.rse.2020.112222
!
!   [3] Ewerdwalbesloh, Y., Montzka, C., Kusche, J., and Springer, A. (2026).
!       GRACE TWSA data assimilation in land surface models-The role of
!       State vector and model uncertainty representation.
!       Water Resources Research, 62, e2025WR042324.
!       https://doi.org/10.1029/2025WR042324
!
! AUTHOR:
!   Lu Li, 12/2024: Initial version
!   Lu Li, 10/2025: Consider correlation & AR(1) process
!   Lu Li, Zhilong Fan, 11/2026: DA v1.0
!-----------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_LandPatch
   USE MOD_SPMD_Task
   USE MOD_Vars_Global
   USE MOD_Vars_1DForcing
   USE MOD_Vars_1DFluxes
   USE MOD_Vars_TimeVariables
   USE MOD_DA_Vars_TimeVars
   IMPLICIT NONE
   SAVE

   ! Parameters
   real(r8)                                  :: dt                                    ! Time step in hours

   ! Forcing perturbations (Table 4 in [1])
   integer,  parameter                       :: nvar_forc = 4                         ! Number of pertutated forcing variables (prcp, sw, lw, t)
   real(r8), parameter                       :: tau_forc = 24.0                       ! Correlation time scale in hours
   real(r8)                                  :: phi_forc                              ! AR(1) autocorrelation coefficient consider time scale
   real(r8)                                  :: sigma_eps_forc                        ! Standard deviation of noise in AR(1) process
   real(r8), dimension(nvar_forc)            :: sigma_forc = (/0.5, 0.3, 20.0, 1.0/)  ! Standard deviation of perturbed forcing variables
   real(r8), dimension(nvar_forc, nvar_forc) :: C = reshape([ &
      1.0, -0.8,  0.5, 0.0, &
      -0.8,  1.0, -0.5, 0.4, &
      0.5, -0.5,  1.0, 0.4, &
      0.0,  0.4,  0.4, 1.0], shape=[nvar_forc, nvar_forc])                            ! Cross-correlation matrix between perturbed forcing variables
   real(r8), allocatable                     :: r_prev_forc (:,:,:)                   ! Previous step forcing perturbations
   real(r8), allocatable                     :: r_curr_forc (:,:,:)                   ! Current step forcing perturbations
   logical                                   :: initialized = .false.                 ! Flag to indicate if is initialized

   ! Soil moisture perturbations
   integer,  parameter                       :: nvar_sm = 2                           ! Number of pertutated soil moisture layers
   real(r8), parameter                       :: tau_sm = 3.0                          ! Correlation time scale in hours
   real(r8)                                  :: phi_sm                                ! AR(1) autocorrelation coefficient consider time scale
   real(r8)                                  :: sigma_eps_sm                          ! Standard deviation of noise in AR(1) process
   real(r8), dimension(nvar_sm)              :: sigma_sm = (/0.035, 0.0552/)          ! Standard deviation of perturbed soil moisture (equal to 0.002 m3/m3)
   real(r8), allocatable                     :: r_prev_sm   (:,:,:)                   ! Previous step soil moisture perturbation
   real(r8), allocatable                     :: r_curr_sm   (:,:,:)                   ! Current step soil moisture perturbation

!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE ensemble(deltim)

!-----------------------------------------------------------------------
! DESCRIPTION:
!   Generate ensemble perturbations to atmospheric forcing and soil
!   Moisture using temporally correlated AR(1) processes. Forcing
!   Perturbations also account for cross-variable error correlations.
!-----------------------------------------------------------------------
      IMPLICIT NONE

      real(r8), intent(in) :: deltim

      integer              :: np, i, j, info
      real(r8)             :: u1(DEF_DA_ENS_NUM/2), u2(DEF_DA_ENS_NUM/2)  ! Uniform random variables

      real(r8)             :: cov_matrix(nvar_forc,nvar_forc)             ! Covariance matrix between perturbed forcing variables
      real(r8)             :: L(nvar_forc,nvar_forc)                      ! Cholesky factor of the forcing covariance matrix
      real(r8)             :: z_forc(nvar_forc,DEF_DA_ENS_NUM)            ! Standard normal random variables for forcing
      real(r8)             :: mean_z_forc(nvar_forc)                      ! Ensemble means of raw forcing normal samples
      real(r8)             :: std_z_forc(nvar_forc)                       ! Sample standard deviations of raw forcing normal samples
      real(r8)             :: mean_r_forc(nvar_forc)                      ! Ensemble means of forcing perturbations after AR(1)
      real(r8)             :: zxL(nvar_forc,DEF_DA_ENS_NUM)               ! Correlated forcing innovations

      real(r8)             :: z_sm(DEF_DA_ENS_NUM)                        ! Standard normal random variables for soil moisture
      real(r8)             :: mean_z_sm                                   ! Ensemble means of soil moisture
      real(r8)             :: std_z_sm                                    ! Sample standard deviations of soil moisture
      real(r8)             :: mean_r_sm(nvar_sm)                          ! Ensemble means of soil moisture perturbations after AR(1)

!-----------------------------------------------------------------------
      ! Initialize persistent variables
      IF (.not. initialized) THEN
         CALL random_seed()
         allocate (r_prev_forc (numpatch,nvar_forc,DEF_DA_ENS_NUM))
         allocate (r_curr_forc (numpatch,nvar_forc,DEF_DA_ENS_NUM))
         allocate (r_prev_sm   (numpatch,nvar_sm,  DEF_DA_ENS_NUM))
         allocate (r_curr_sm   (numpatch,nvar_sm,  DEF_DA_ENS_NUM))
         r_prev_forc = 0.0_r8
         r_curr_forc = 0.0_r8
         r_prev_sm   = 0.0_r8
         r_curr_sm   = 0.0_r8
         initialized = .true.
      ENDIF

      ! Calculate AR(1) parameters
      dt = deltim/3600
      phi_forc = exp(-dt/tau_forc)
      sigma_eps_forc = sqrt(1.0 - phi_forc**2)
      phi_sm = exp(-dt/tau_sm)
      sigma_eps_sm = sqrt(1.0 - phi_sm**2)

      ! Calculate covariance matrix by cross-correlation matrix and standard deviation
      cov_matrix = 0.0_r8
      DO i = 1, nvar_forc
         DO j = 1, nvar_forc
            cov_matrix(i,j) = C(i,j) * sigma_forc(i) * sigma_forc(j)
         ENDDO
      ENDDO

      ! Perform Cholesky decomposition of covariance matrix
      L = cov_matrix
      CALL dpotrf('L', nvar_forc, L, nvar_forc, info)
      DO i = 1, nvar_forc
         DO j = i+1, nvar_forc
            L(i,j) = 0.0_r8
         ENDDO
      ENDDO

      ! Generate ensemble samples for forcing variables
      DO np = 1, numpatch
         ! Generate disturbance ensemble samples ~ N(0, I)
         CALL random_number(u1)
         CALL random_number(u2)
         DO i = 1, DEF_DA_ENS_NUM/2
            u1(i) = max(u1(i), 1e-10)  ! Ensure u1 is not zero
            z_forc(1,i*2-1) = sqrt(-2.0*log(u1(i))) * cos(2.0*pi*u2(i))
            z_forc(1,i*2)   = sqrt(-2.0*log(u1(i))) * sin(2.0*pi*u2(i))
         ENDDO
         CALL random_number(u1)
         CALL random_number(u2)
         DO i = 1, DEF_DA_ENS_NUM/2
            u1(i) = max(u1(i), 1e-10)
            z_forc(2,i*2-1) = sqrt(-2.0*log(u1(i))) * cos(2.0*pi*u2(i))
            z_forc(2,i*2)   = sqrt(-2.0*log(u1(i))) * sin(2.0*pi*u2(i))
         ENDDO
         CALL random_number(u1)
         CALL random_number(u2)
         DO i = 1, DEF_DA_ENS_NUM/2
            u1(i) = max(u1(i), 1e-10)
            z_forc(3,i*2-1) = sqrt(-2.0*log(u1(i))) * cos(2.0*pi*u2(i))
            z_forc(3,i*2)   = sqrt(-2.0*log(u1(i))) * sin(2.0*pi*u2(i))
         ENDDO
         CALL random_number(u1)
         CALL random_number(u2)
         DO i = 1, DEF_DA_ENS_NUM/2
            u1(i) = max(u1(i), 1e-10)
            z_forc(4,i*2-1) = sqrt(-2.0*log(u1(i))) * cos(2.0*pi*u2(i))
            z_forc(4,i*2)   = sqrt(-2.0*log(u1(i))) * sin(2.0*pi*u2(i))
         ENDDO

         ! Normalize z to mean 0 and std 1
         DO i = 1, nvar_forc
            mean_z_forc(i) = sum(z_forc(i, :))/DEF_DA_ENS_NUM
            std_z_forc(i) = sqrt(sum((z_forc(i, :) - mean_z_forc(i))**2)/(DEF_DA_ENS_NUM - 1))
         ENDDO
         DO i = 1, nvar_forc
            z_forc(i,:) = (z_forc(i,:)-mean_z_forc(i))/std_z_forc(i)
         ENDDO

         ! Apply the Cholesky factor to introduce forcing cross-correlations.
         CALL dgemm('N', 'N', nvar_forc, DEF_DA_ENS_NUM, nvar_forc, &
            1.0_r8, L, nvar_forc, z_forc, nvar_forc, &
            0.0_r8, zxL, nvar_forc)

         ! Introduce correlation using AR(1) process
         DO i = 1, nvar_forc
            DO j = 1, DEF_DA_ENS_NUM
               r_curr_forc(np,i,j) = phi_forc * r_prev_forc(np,i,j) + sigma_eps_forc * zxL(i,j)
            ENDDO
         ENDDO
         ! No AR(1) process, directly use correlated random variables
         ! R_curr_forc(np,:,:) = zxL

         ! Normalize the disturbance ensemble samples to mean 0
         mean_r_forc = sum(r_curr_forc(np,:,:), dim=2)/DEF_DA_ENS_NUM
         DO i = 1, nvar_forc
            DO j = 1, DEF_DA_ENS_NUM
               r_curr_forc(np,i,j) = r_curr_forc(np,i,j) - mean_r_forc(i)
            ENDDO
         ENDDO

         ! Save current perturbation as previous perturbation for next time step
         r_prev_forc = r_curr_forc

         ! Generate ensemble samples according different types
         DO j = 1, DEF_DA_ENS_NUM
            forc_prc_ens  (j,np) = forc_prc  (np) * exp(r_curr_forc(np,1,j) - 0.5 * sigma_forc(1)**2)
            forc_prl_ens  (j,np) = forc_prl  (np) * exp(r_curr_forc(np,1,j) - 0.5 * sigma_forc(1)**2)
            forc_sols_ens (j,np) = forc_sols (np) * exp(r_curr_forc(np,2,j) - 0.5 * sigma_forc(2)**2)
            forc_soll_ens (j,np) = forc_soll (np) * exp(r_curr_forc(np,2,j) - 0.5 * sigma_forc(2)**2)
            forc_solsd_ens(j,np) = forc_solsd(np) * exp(r_curr_forc(np,2,j) - 0.5 * sigma_forc(2)**2)
            forc_solld_ens(j,np) = forc_solld(np) * exp(r_curr_forc(np,2,j) - 0.5 * sigma_forc(2)**2)
            forc_frl_ens  (j,np) = forc_frl  (np) + r_curr_forc(np,3,j)
            forc_t_ens    (j,np) = forc_t    (np) + r_curr_forc(np,4,j)
         ENDDO

         IF (DEF_DA_ENS_SM) THEN
            ! Generate ensemble samples (0, I) for soil moisture
            CALL random_number(u1)
            CALL random_number(u2)
            DO i = 1, DEF_DA_ENS_NUM/2
               u1(i) = max(u1(i), 1e-10)
               z_sm(i*2-1) = sqrt(-2.0*log(u1(i))) * cos(2.0*pi*u2(i))
               z_sm(i*2)   = sqrt(-2.0*log(u1(i))) * sin(2.0*pi*u2(i))
            ENDDO
            mean_z_sm = sum(z_sm)/DEF_DA_ENS_NUM
            std_z_sm = sqrt(sum((z_sm - mean_z_sm)**2)/(DEF_DA_ENS_NUM - 1))
            z_sm = (z_sm - mean_z_sm)/std_z_sm

            ! Introduce correlation using AR(1) process
            DO i = 1, nvar_sm
               DO j = 1, DEF_DA_ENS_NUM
                  r_curr_sm(np,i,j) = phi_sm * r_prev_sm(np,i,j) + sigma_eps_sm * sigma_sm(i) * z_sm(j)
               ENDDO
            ENDDO

            ! Normalize the disturbance ensemble samples to mean 0
            mean_r_sm = sum(r_curr_sm(np,:,:), dim=2)/DEF_DA_ENS_NUM
            DO i = 1, nvar_sm
               DO j = 1, DEF_DA_ENS_NUM
                  r_curr_sm(np,i,j) = r_curr_sm(np,i,j) - mean_r_sm(i)
               ENDDO
            ENDDO

            ! Save current perturbation as previous perturbation for next time step
            r_prev_sm   = r_curr_sm

            ! Apply layer perturbations and the existing groundwater adjustment.
            wa_ens(1:DEF_DA_ENS_NUM,np) = wa_ens(1:DEF_DA_ENS_NUM,np) + &
               sum(max(1e-10, wliq_soisno_ens(1:nvar_sm,1:DEF_DA_ENS_NUM,np) + r_curr_sm(np,:,:)) - &
               wliq_soisno_ens(1:nvar_sm,1:DEF_DA_ENS_NUM,np), dim=1)
            wliq_soisno_ens(1:nvar_sm,1:DEF_DA_ENS_NUM,np) = &
               max(1e-10, wliq_soisno_ens(1:nvar_sm,1:DEF_DA_ENS_NUM,np) + r_curr_sm(np,:,:))
         ENDIF
      ENDDO

   END SUBROUTINE ensemble

!-----------------------------------------------------------------------

   SUBROUTINE save_member(iens)

!-----------------------------------------------------------------------
! DESCRIPTION:
!   Save the current model states, forcing, and diagnostic and flux
!   Outputs into ensemble member.
!-----------------------------------------------------------------------
      IMPLICIT NONE
      integer, intent(in) :: iens

      ! Restart states
      z_sno_ens       (  :,iens,:) = z_sno
      dz_sno_ens      (  :,iens,:) = dz_sno
      t_soisno_ens    (  :,iens,:) = t_soisno
      wliq_soisno_ens (  :,iens,:) = wliq_soisno
      wice_soisno_ens (  :,iens,:) = wice_soisno
      smp_ens         (  :,iens,:) = smp
      hk_ens          (  :,iens,:) = hk
      t_grnd_ens      (    iens,:) = t_grnd
      tleaf_ens       (    iens,:) = tleaf
      ldew_ens        (    iens,:) = ldew
      ldew_rain_ens   (    iens,:) = ldew_rain
      ldew_snow_ens   (    iens,:) = ldew_snow
      fwet_snow_ens   (    iens,:) = fwet_snow
      sag_ens         (    iens,:) = sag
      scv_ens         (    iens,:) = scv
      snowdp_ens      (    iens,:) = snowdp
      fveg_ens        (    iens,:) = fveg
      fsno_ens        (    iens,:) = fsno
      sigf_ens        (    iens,:) = sigf
      green_ens       (    iens,:) = green
      tlai_ens        (    iens,:) = tlai
      lai_ens         (    iens,:) = lai
      tsai_ens        (    iens,:) = tsai
      sai_ens         (    iens,:) = sai
      alb_ens         (:,:,iens,:) = alb
      ssun_ens        (:,:,iens,:) = ssun
      ssha_ens        (:,:,iens,:) = ssha
      ssoi_ens        (:,:,iens,:) = ssoi
      ssno_ens        (:,:,iens,:) = ssno
      thermk_ens      (    iens,:) = thermk
      extkb_ens       (    iens,:) = extkb
      extkd_ens       (    iens,:) = extkd
      zwt_ens         (    iens,:) = zwt
      wdsrf_ens       (    iens,:) = wdsrf
      wa_ens          (    iens,:) = wa
      wetwat_ens      (    iens,:) = wetwat
      t_lake_ens      (  :,iens,:) = t_lake
      lake_icefrac_ens(  :,iens,:) = lake_icefrac
      savedtke1_ens   (    iens,:) = savedtke1

      ! Land-atmosphere coupling states
      trad_ens        (    iens,:) = trad
      tref_ens        (    iens,:) = tref
      qref_ens        (    iens,:) = qref
      rst_ens         (    iens,:) = rst
      emis_ens        (    iens,:) = emis
      z0m_ens         (    iens,:) = z0m
      zol_ens         (    iens,:) = zol
      rib_ens         (    iens,:) = rib
      ustar_ens       (    iens,:) = ustar
      qstar_ens       (    iens,:) = qstar
      tstar_ens       (    iens,:) = tstar
      fm_ens          (    iens,:) = fm
      fh_ens          (    iens,:) = fh
      fq_ens          (    iens,:) = fq

      ! Forcings
      forc_t_ens      (    iens,:) = forc_t
      forc_frl_ens    (    iens,:) = forc_frl
      forc_prc_ens    (    iens,:) = forc_prc
      forc_prl_ens    (    iens,:) = forc_prl
      forc_sols_ens   (    iens,:) = forc_sols
      forc_soll_ens   (    iens,:) = forc_soll
      forc_solsd_ens  (    iens,:) = forc_solsd
      forc_solld_ens  (    iens,:) = forc_solld

      ! Surface fluxes
      taux_ens        (    iens,:) = taux
      tauy_ens        (    iens,:) = tauy
      fsena_ens       (    iens,:) = fsena
      lfevpa_ens      (    iens,:) = lfevpa
      fevpa_ens       (    iens,:) = fevpa
      fsenl_ens       (    iens,:) = fsenl
      fevpl_ens       (    iens,:) = fevpl
      etr_ens         (    iens,:) = etr
      fseng_ens       (    iens,:) = fseng
      fevpg_ens       (    iens,:) = fevpg
      fgrnd_ens       (    iens,:) = fgrnd
      sabvsun_ens     (    iens,:) = sabvsun
      sabvsha_ens     (    iens,:) = sabvsha
      sabg_ens        (    iens,:) = sabg
      sr_ens          (    iens,:) = sr
      solvd_ens       (    iens,:) = solvd
      solvi_ens       (    iens,:) = solvi
      solnd_ens       (    iens,:) = solnd
      solni_ens       (    iens,:) = solni
      srvd_ens        (    iens,:) = srvd
      srvi_ens        (    iens,:) = srvi
      srnd_ens        (    iens,:) = srnd
      srni_ens        (    iens,:) = srni
      solvdln_ens     (    iens,:) = solvdln
      solviln_ens     (    iens,:) = solviln
      solndln_ens     (    iens,:) = solndln
      solniln_ens     (    iens,:) = solniln
      srvdln_ens      (    iens,:) = srvdln
      srviln_ens      (    iens,:) = srviln
      srndln_ens      (    iens,:) = srndln
      srniln_ens      (    iens,:) = srniln
      olrg_ens        (    iens,:) = olrg
      rnet_ens        (    iens,:) = rnet
      xerr_ens        (    iens,:) = xerr
      zerr_ens        (    iens,:) = zerr
      frcsat_ens      (    iens,:) = frcsat
      rsur_ens        (    iens,:) = rsur
      rsur_se_ens     (    iens,:) = rsur_se
      rsur_ie_ens     (    iens,:) = rsur_ie
      rsub_ens        (    iens,:) = rsub
      rnof_ens        (    iens,:) = rnof
      qintr_ens       (    iens,:) = qintr
      qinfl_ens       (    iens,:) = qinfl
      qdrip_ens       (    iens,:) = qdrip
      assim_ens       (    iens,:) = assim
      respc_ens       (    iens,:) = respc
      qcharge_ens     (    iens,:) = qcharge
      qlayer_ens      (  :,iens,:) = qlayer
      lake_deficit_ens(    iens,:) = lake_deficit

   END SUBROUTINE save_member

!-----------------------------------------------------------------------------

   SUBROUTINE load_member(iens)

!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Load the states, forcing, and flux outputs of an ensemble member
!   Into the model variables.
!-----------------------------------------------------------------------------
      IMPLICIT NONE
      integer, intent(in) :: iens

      ! Restart states
      z_sno        = z_sno_ens       (  :,iens,:)
      dz_sno       = dz_sno_ens      (  :,iens,:)
      t_soisno     = t_soisno_ens    (  :,iens,:)
      wliq_soisno  = wliq_soisno_ens (  :,iens,:)
      wice_soisno  = wice_soisno_ens (  :,iens,:)
      smp          = smp_ens         (  :,iens,:)
      hk           = hk_ens          (  :,iens,:)
      t_grnd       = t_grnd_ens      (    iens,:)
      tleaf        = tleaf_ens       (    iens,:)
      ldew         = ldew_ens        (    iens,:)
      ldew_rain    = ldew_rain_ens   (    iens,:)
      ldew_snow    = ldew_snow_ens   (    iens,:)
      fwet_snow    = fwet_snow_ens   (    iens,:)
      sag          = sag_ens         (    iens,:)
      scv          = scv_ens         (    iens,:)
      snowdp       = snowdp_ens      (    iens,:)
      fveg         = fveg_ens        (    iens,:)
      fsno         = fsno_ens        (    iens,:)
      sigf         = sigf_ens        (    iens,:)
      green        = green_ens       (    iens,:)
      tlai         = tlai_ens        (    iens,:)
      lai          = lai_ens         (    iens,:)
      tsai         = tsai_ens        (    iens,:)
      sai          = sai_ens         (    iens,:)
      alb          = alb_ens         (:,:,iens,:)
      ssun         = ssun_ens        (:,:,iens,:)
      ssha         = ssha_ens        (:,:,iens,:)
      ssoi         = ssoi_ens        (:,:,iens,:)
      ssno         = ssno_ens        (:,:,iens,:)
      thermk       = thermk_ens      (    iens,:)
      extkb        = extkb_ens       (    iens,:)
      extkd        = extkd_ens       (    iens,:)
      zwt          = zwt_ens         (    iens,:)
      wdsrf        = wdsrf_ens       (    iens,:)
      wa           = wa_ens          (    iens,:)
      wetwat       = wetwat_ens      (    iens,:)
      t_lake       = t_lake_ens      (  :,iens,:)
      lake_icefrac = lake_icefrac_ens(  :,iens,:)
      savedtke1    = savedtke1_ens   (    iens,:)

      ! Land-atmosphere coupling states
      trad         = trad_ens        (    iens,:)
      tref         = tref_ens        (    iens,:)
      qref         = qref_ens        (    iens,:)
      rst          = rst_ens         (    iens,:)
      emis         = emis_ens        (    iens,:)
      z0m          = z0m_ens         (    iens,:)
      zol          = zol_ens         (    iens,:)
      rib          = rib_ens         (    iens,:)
      ustar        = ustar_ens       (    iens,:)
      qstar        = qstar_ens       (    iens,:)
      tstar        = tstar_ens       (    iens,:)
      fm           = fm_ens          (    iens,:)
      fh           = fh_ens          (    iens,:)
      fq           = fq_ens          (    iens,:)

      ! Forcings
      forc_t       = forc_t_ens      (    iens,:)
      forc_frl     = forc_frl_ens    (    iens,:)
      forc_prc     = forc_prc_ens    (    iens,:)
      forc_prl     = forc_prl_ens    (    iens,:)
      forc_sols    = forc_sols_ens   (    iens,:)
      forc_soll    = forc_soll_ens   (    iens,:)
      forc_solsd   = forc_solsd_ens  (    iens,:)
      forc_solld   = forc_solld_ens  (    iens,:)

      ! Surface fluxes
      taux         = taux_ens        (    iens,:)
      tauy         = tauy_ens        (    iens,:)
      fsena        = fsena_ens       (    iens,:)
      lfevpa       = lfevpa_ens      (    iens,:)
      fevpa        = fevpa_ens       (    iens,:)
      fsenl        = fsenl_ens       (    iens,:)
      fevpl        = fevpl_ens       (    iens,:)
      etr          = etr_ens         (    iens,:)
      fseng        = fseng_ens       (    iens,:)
      fevpg        = fevpg_ens       (    iens,:)
      fgrnd        = fgrnd_ens       (    iens,:)
      sabvsun      = sabvsun_ens     (    iens,:)
      sabvsha      = sabvsha_ens     (    iens,:)
      sabg         = sabg_ens        (    iens,:)
      sr           = sr_ens          (    iens,:)
      solvd        = solvd_ens       (    iens,:)
      solvi        = solvi_ens       (    iens,:)
      solnd        = solnd_ens       (    iens,:)
      solni        = solni_ens       (    iens,:)
      srvd         = srvd_ens        (    iens,:)
      srvi         = srvi_ens        (    iens,:)
      srnd         = srnd_ens        (    iens,:)
      srni         = srni_ens        (    iens,:)
      solvdln      = solvdln_ens     (    iens,:)
      solviln      = solviln_ens     (    iens,:)
      solndln      = solndln_ens     (    iens,:)
      solniln      = solniln_ens     (    iens,:)
      srvdln       = srvdln_ens      (    iens,:)
      srviln       = srviln_ens      (    iens,:)
      srndln       = srndln_ens      (    iens,:)
      srniln       = srniln_ens      (    iens,:)
      olrg         = olrg_ens        (    iens,:)
      rnet         = rnet_ens        (    iens,:)
      xerr         = xerr_ens        (    iens,:)
      zerr         = zerr_ens        (    iens,:)
      frcsat       = frcsat_ens      (    iens,:)
      rsur         = rsur_ens        (    iens,:)
      rsur_se      = rsur_se_ens     (    iens,:)
      rsur_ie      = rsur_ie_ens     (    iens,:)
      rsub         = rsub_ens        (    iens,:)
      rnof         = rnof_ens        (    iens,:)
      qintr        = qintr_ens       (    iens,:)
      qinfl        = qinfl_ens       (    iens,:)
      qdrip        = qdrip_ens       (    iens,:)
      assim        = assim_ens       (    iens,:)
      respc        = respc_ens       (    iens,:)
      qcharge      = qcharge_ens     (    iens,:)
      qlayer       = qlayer_ens      (  :,iens,:)
      lake_deficit = lake_deficit_ens(    iens,:)

   END SUBROUTINE load_member

!-----------------------------------------------------------------------------

   SUBROUTINE postprocess()

!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Set the final model states and outputs after data assimilation.
!
!   When land-atmosphere feedback is enabled, apply the analysis
!   Ensemble means for restart states and fluxes.
!
!   Otherwise, restore the open-loop trajectory for diagnostic.
!-----------------------------------------------------------------------------
      IMPLICIT NONE

      CALL load_member(0)

      IF (DEF_DA_ONLINE) THEN
         ! Restart states
         z_sno        = sum(z_sno_ens(:,1:DEF_DA_ENS_NUM,:), dim=2)        / DEF_DA_ENS_NUM
         dz_sno       = sum(dz_sno_ens(:,1:DEF_DA_ENS_NUM,:), dim=2)       / DEF_DA_ENS_NUM
         t_soisno     = sum(t_soisno_ens(:,1:DEF_DA_ENS_NUM,:), dim=2)     / DEF_DA_ENS_NUM
         wliq_soisno  = sum(wliq_soisno_ens(:,1:DEF_DA_ENS_NUM,:), dim=2)  / DEF_DA_ENS_NUM
         wice_soisno  = sum(wice_soisno_ens(:,1:DEF_DA_ENS_NUM,:), dim=2)  / DEF_DA_ENS_NUM
         smp          = sum(smp_ens(:,1:DEF_DA_ENS_NUM,:), dim=2)          / DEF_DA_ENS_NUM
         hk           = sum(hk_ens(:,1:DEF_DA_ENS_NUM,:), dim=2)           / DEF_DA_ENS_NUM
         t_grnd       = sum(t_grnd_ens(1:DEF_DA_ENS_NUM,:), dim=1)         / DEF_DA_ENS_NUM
         tleaf        = sum(tleaf_ens(1:DEF_DA_ENS_NUM,:), dim=1)          / DEF_DA_ENS_NUM
         ldew         = sum(ldew_ens(1:DEF_DA_ENS_NUM,:), dim=1)           / DEF_DA_ENS_NUM
         ldew_rain    = sum(ldew_rain_ens(1:DEF_DA_ENS_NUM,:), dim=1)      / DEF_DA_ENS_NUM
         ldew_snow    = sum(ldew_snow_ens(1:DEF_DA_ENS_NUM,:), dim=1)      / DEF_DA_ENS_NUM
         fwet_snow    = sum(fwet_snow_ens(1:DEF_DA_ENS_NUM,:), dim=1)      / DEF_DA_ENS_NUM
         sag          = sum(sag_ens(1:DEF_DA_ENS_NUM,:), dim=1)            / DEF_DA_ENS_NUM
         scv          = sum(scv_ens(1:DEF_DA_ENS_NUM,:), dim=1)            / DEF_DA_ENS_NUM
         snowdp       = sum(snowdp_ens(1:DEF_DA_ENS_NUM,:), dim=1)         / DEF_DA_ENS_NUM
         fveg         = sum(fveg_ens(1:DEF_DA_ENS_NUM,:), dim=1)           / DEF_DA_ENS_NUM
         fsno         = sum(fsno_ens(1:DEF_DA_ENS_NUM,:), dim=1)           / DEF_DA_ENS_NUM
         sigf         = sum(sigf_ens(1:DEF_DA_ENS_NUM,:), dim=1)           / DEF_DA_ENS_NUM
         green        = sum(green_ens(1:DEF_DA_ENS_NUM,:), dim=1)          / DEF_DA_ENS_NUM
         tlai         = sum(tlai_ens(1:DEF_DA_ENS_NUM,:), dim=1)           / DEF_DA_ENS_NUM
         lai          = sum(lai_ens(1:DEF_DA_ENS_NUM,:), dim=1)            / DEF_DA_ENS_NUM
         tsai         = sum(tsai_ens(1:DEF_DA_ENS_NUM,:), dim=1)           / DEF_DA_ENS_NUM
         sai          = sum(sai_ens(1:DEF_DA_ENS_NUM,:), dim=1)            / DEF_DA_ENS_NUM
         alb          = sum(alb_ens(:,:,1:DEF_DA_ENS_NUM,:), dim=3)        / DEF_DA_ENS_NUM
         ssun         = sum(ssun_ens(:,:,1:DEF_DA_ENS_NUM,:), dim=3)       / DEF_DA_ENS_NUM
         ssha         = sum(ssha_ens(:,:,1:DEF_DA_ENS_NUM,:), dim=3)       / DEF_DA_ENS_NUM
         ssoi         = sum(ssoi_ens(:,:,1:DEF_DA_ENS_NUM,:), dim=3)       / DEF_DA_ENS_NUM
         ssno         = sum(ssno_ens(:,:,1:DEF_DA_ENS_NUM,:), dim=3)       / DEF_DA_ENS_NUM
         thermk       = sum(thermk_ens(1:DEF_DA_ENS_NUM,:), dim=1)         / DEF_DA_ENS_NUM
         extkb        = sum(extkb_ens(1:DEF_DA_ENS_NUM,:), dim=1)          / DEF_DA_ENS_NUM
         extkd        = sum(extkd_ens(1:DEF_DA_ENS_NUM,:), dim=1)          / DEF_DA_ENS_NUM
         zwt          = sum(zwt_ens(1:DEF_DA_ENS_NUM,:), dim=1)            / DEF_DA_ENS_NUM
         wdsrf        = sum(wdsrf_ens(1:DEF_DA_ENS_NUM,:), dim=1)          / DEF_DA_ENS_NUM
         wa           = sum(wa_ens(1:DEF_DA_ENS_NUM,:), dim=1)             / DEF_DA_ENS_NUM
         wetwat       = sum(wetwat_ens(1:DEF_DA_ENS_NUM,:), dim=1)         / DEF_DA_ENS_NUM
         t_lake       = sum(t_lake_ens(:,1:DEF_DA_ENS_NUM,:), dim=2)       / DEF_DA_ENS_NUM
         lake_icefrac = sum(lake_icefrac_ens(:,1:DEF_DA_ENS_NUM,:), dim=2) / DEF_DA_ENS_NUM
         savedtke1    = sum(savedtke1_ens(1:DEF_DA_ENS_NUM,:), dim=1)      / DEF_DA_ENS_NUM

         ! Land-atmosphere coupling states
         trad         = sum(trad_ens(1:DEF_DA_ENS_NUM,:), dim=1)           / DEF_DA_ENS_NUM
         tref         = sum(tref_ens(1:DEF_DA_ENS_NUM,:), dim=1)           / DEF_DA_ENS_NUM
         qref         = sum(qref_ens(1:DEF_DA_ENS_NUM,:), dim=1)           / DEF_DA_ENS_NUM
         rst          = sum(rst_ens(1:DEF_DA_ENS_NUM,:), dim=1)            / DEF_DA_ENS_NUM
         emis         = sum(emis_ens(1:DEF_DA_ENS_NUM,:), dim=1)           / DEF_DA_ENS_NUM
         z0m          = sum(z0m_ens(1:DEF_DA_ENS_NUM,:), dim=1)            / DEF_DA_ENS_NUM
         zol          = sum(zol_ens(1:DEF_DA_ENS_NUM,:), dim=1)            / DEF_DA_ENS_NUM
         rib          = sum(rib_ens(1:DEF_DA_ENS_NUM,:), dim=1)            / DEF_DA_ENS_NUM
         ustar        = sum(ustar_ens(1:DEF_DA_ENS_NUM,:), dim=1)          / DEF_DA_ENS_NUM
         qstar        = sum(qstar_ens(1:DEF_DA_ENS_NUM,:), dim=1)          / DEF_DA_ENS_NUM
         tstar        = sum(tstar_ens(1:DEF_DA_ENS_NUM,:), dim=1)          / DEF_DA_ENS_NUM
         fm           = sum(fm_ens(1:DEF_DA_ENS_NUM,:), dim=1)             / DEF_DA_ENS_NUM
         fh           = sum(fh_ens(1:DEF_DA_ENS_NUM,:), dim=1)             / DEF_DA_ENS_NUM
         fq           = sum(fq_ens(1:DEF_DA_ENS_NUM,:), dim=1)             / DEF_DA_ENS_NUM

         ! Surface fluxes
         taux         = sum(taux_ens(1:DEF_DA_ENS_NUM,:), dim=1)           / DEF_DA_ENS_NUM
         tauy         = sum(tauy_ens(1:DEF_DA_ENS_NUM,:), dim=1)           / DEF_DA_ENS_NUM
         fsena        = sum(fsena_ens(1:DEF_DA_ENS_NUM,:), dim=1)          / DEF_DA_ENS_NUM
         lfevpa       = sum(lfevpa_ens(1:DEF_DA_ENS_NUM,:), dim=1)         / DEF_DA_ENS_NUM
         fevpa        = sum(fevpa_ens(1:DEF_DA_ENS_NUM,:), dim=1)          / DEF_DA_ENS_NUM
         fsenl        = sum(fsenl_ens(1:DEF_DA_ENS_NUM,:), dim=1)          / DEF_DA_ENS_NUM
         fevpl        = sum(fevpl_ens(1:DEF_DA_ENS_NUM,:), dim=1)          / DEF_DA_ENS_NUM
         etr          = sum(etr_ens(1:DEF_DA_ENS_NUM,:), dim=1)            / DEF_DA_ENS_NUM
         fseng        = sum(fseng_ens(1:DEF_DA_ENS_NUM,:), dim=1)          / DEF_DA_ENS_NUM
         fevpg        = sum(fevpg_ens(1:DEF_DA_ENS_NUM,:), dim=1)          / DEF_DA_ENS_NUM
         fgrnd        = sum(fgrnd_ens(1:DEF_DA_ENS_NUM,:), dim=1)          / DEF_DA_ENS_NUM
         sabvsun      = sum(sabvsun_ens(1:DEF_DA_ENS_NUM,:), dim=1)        / DEF_DA_ENS_NUM
         sabvsha      = sum(sabvsha_ens(1:DEF_DA_ENS_NUM,:), dim=1)        / DEF_DA_ENS_NUM
         sabg         = sum(sabg_ens(1:DEF_DA_ENS_NUM,:), dim=1)           / DEF_DA_ENS_NUM
         sr           = sum(sr_ens(1:DEF_DA_ENS_NUM,:), dim=1)             / DEF_DA_ENS_NUM
         solvd        = sum(solvd_ens(1:DEF_DA_ENS_NUM,:), dim=1)          / DEF_DA_ENS_NUM
         solvi        = sum(solvi_ens(1:DEF_DA_ENS_NUM,:), dim=1)          / DEF_DA_ENS_NUM
         solnd        = sum(solnd_ens(1:DEF_DA_ENS_NUM,:), dim=1)          / DEF_DA_ENS_NUM
         solni        = sum(solni_ens(1:DEF_DA_ENS_NUM,:), dim=1)          / DEF_DA_ENS_NUM
         srvd         = sum(srvd_ens(1:DEF_DA_ENS_NUM,:), dim=1)           / DEF_DA_ENS_NUM
         srvi         = sum(srvi_ens(1:DEF_DA_ENS_NUM,:), dim=1)           / DEF_DA_ENS_NUM
         srnd         = sum(srnd_ens(1:DEF_DA_ENS_NUM,:), dim=1)           / DEF_DA_ENS_NUM
         srni         = sum(srni_ens(1:DEF_DA_ENS_NUM,:), dim=1)           / DEF_DA_ENS_NUM
         solvdln      = sum(solvdln_ens(1:DEF_DA_ENS_NUM,:), dim=1)        / DEF_DA_ENS_NUM
         solviln      = sum(solviln_ens(1:DEF_DA_ENS_NUM,:), dim=1)        / DEF_DA_ENS_NUM
         solndln      = sum(solndln_ens(1:DEF_DA_ENS_NUM,:), dim=1)        / DEF_DA_ENS_NUM
         solniln      = sum(solniln_ens(1:DEF_DA_ENS_NUM,:), dim=1)        / DEF_DA_ENS_NUM
         srvdln       = sum(srvdln_ens(1:DEF_DA_ENS_NUM,:), dim=1)         / DEF_DA_ENS_NUM
         srviln       = sum(srviln_ens(1:DEF_DA_ENS_NUM,:), dim=1)         / DEF_DA_ENS_NUM
         srndln       = sum(srndln_ens(1:DEF_DA_ENS_NUM,:), dim=1)         / DEF_DA_ENS_NUM
         srniln       = sum(srniln_ens(1:DEF_DA_ENS_NUM,:), dim=1)         / DEF_DA_ENS_NUM
         olrg         = sum(olrg_ens(1:DEF_DA_ENS_NUM,:), dim=1)           / DEF_DA_ENS_NUM
         rnet         = sum(rnet_ens(1:DEF_DA_ENS_NUM,:), dim=1)           / DEF_DA_ENS_NUM
         xerr         = sum(xerr_ens(1:DEF_DA_ENS_NUM,:), dim=1)           / DEF_DA_ENS_NUM
         zerr         = sum(zerr_ens(1:DEF_DA_ENS_NUM,:), dim=1)           / DEF_DA_ENS_NUM
         frcsat       = sum(frcsat_ens(1:DEF_DA_ENS_NUM,:), dim=1)         / DEF_DA_ENS_NUM
         rsur         = sum(rsur_ens(1:DEF_DA_ENS_NUM,:), dim=1)           / DEF_DA_ENS_NUM
         rsur_se      = sum(rsur_se_ens(1:DEF_DA_ENS_NUM,:), dim=1)        / DEF_DA_ENS_NUM
         rsur_ie      = sum(rsur_ie_ens(1:DEF_DA_ENS_NUM,:), dim=1)        / DEF_DA_ENS_NUM
         rsub         = sum(rsub_ens(1:DEF_DA_ENS_NUM,:), dim=1)           / DEF_DA_ENS_NUM
         rnof         = sum(rnof_ens(1:DEF_DA_ENS_NUM,:), dim=1)           / DEF_DA_ENS_NUM
         qintr        = sum(qintr_ens(1:DEF_DA_ENS_NUM,:), dim=1)          / DEF_DA_ENS_NUM
         qinfl        = sum(qinfl_ens(1:DEF_DA_ENS_NUM,:), dim=1)          / DEF_DA_ENS_NUM
         qdrip        = sum(qdrip_ens(1:DEF_DA_ENS_NUM,:), dim=1)          / DEF_DA_ENS_NUM
         assim        = sum(assim_ens(1:DEF_DA_ENS_NUM,:), dim=1)          / DEF_DA_ENS_NUM
         respc        = sum(respc_ens(1:DEF_DA_ENS_NUM,:), dim=1)          / DEF_DA_ENS_NUM
         qcharge      = sum(qcharge_ens(1:DEF_DA_ENS_NUM,:), dim=1)        / DEF_DA_ENS_NUM
         qlayer       = sum(qlayer_ens(:,1:DEF_DA_ENS_NUM,:), dim=2)       / DEF_DA_ENS_NUM
         lake_deficit = sum(lake_deficit_ens(1:DEF_DA_ENS_NUM,:), dim=1)   / DEF_DA_ENS_NUM
      ENDIF

   END SUBROUTINE postprocess

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Ens
#endif