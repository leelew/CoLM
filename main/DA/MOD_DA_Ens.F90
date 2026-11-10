#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Ens
!-----------------------------------------------------------------------
! DESCRIPTION:
!    Provide procedures for ensemble generation and management in data
!    assimilation, including perturbation generation, member state exchange,
!    and analysis feedback to the model.
!
! REFERENCES:
!    [1] Algorithm Theoretical Basis Document Level 4 Surface and Root
!        Zone Soil Moisture (L4_SM) Data Product
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
   USE MOD_Vars_TimeInvariants, only: porsl, hksati, &
      vf_sand, vf_clay, patchtype, z_soil_ptb
   USE MOD_Vars_1DForcing
   USE MOD_Vars_1DFluxes
   USE MOD_Vars_TimeVariables
   USE MOD_DA_Vars_Fluxes
   USE MOD_DA_Vars_TimeVars
   IMPLICIT NONE
   SAVE

   ! parameters
   real(r8) :: dt                                                          ! time step in hours

   ! forcing [parameters used here is consistent with SMAP L4 (Table 4 in [1])]
   integer,  parameter :: nvar_forc = 4                                    ! number of pertutated forcing variables
   real(r8), parameter :: tau_forc = 24.0                                  ! correlation time scale in hours

   real(r8) :: phi_forc                                                    ! AR(1) autocorrelation coefficient consider time scale
   real(r8) :: sigma_eps_forc                                              ! standard deviation of noise in AR(1) process
   real(r8), dimension(nvar_forc) :: sigma_forc = (/0.5, 0.3, 20.0, 1.0/)  ! standard deviation of perturbed forcing variables (prcp, sw, lw, t)
   real(r8), dimension(nvar_forc, nvar_forc) :: C = reshape([ &
         1.0, -0.8,  0.5, 0.0, &
        -0.8,  1.0, -0.5, 0.4, &
         0.5, -0.5,  1.0, 0.4, &
         0.0,  0.4,  0.4, 1.0], shape=[nvar_forc, nvar_forc])              ! cross-correlation matrix between perturbed forcing variables
   real(r8), allocatable :: r_prev_forc(:,:,:)                             ! previous perturbation (numpatch, nvar_forc, DEF_DA_ENS_NUM)
   real(r8), allocatable :: r_curr_forc(:,:,:)                             ! current perturbation (numpatch, nvar_forc, DEF_DA_ENS_NUM)
   logical :: initialized = .false.                                        ! flag to indicate if is initialized

   ! soil moisture [10-layer scheme from the legacy MOD_DA_Ensemble]
   integer,  parameter :: nvar_sm = nl_soil                                ! number of perturbed soil moisture layers
   real(r8), parameter :: tau_sm = 3.0                                     ! correlation time scale in hours
   real(r8) :: phi_sm                                                      ! AR(1) autocorrelation coefficient consider time scale
   real(r8) :: sigma_eps_sm                                                ! standard deviation of noise in AR(1) process
   real(r8), dimension(nvar_sm) :: sigma_sm = (/ &
         0.105_r8, 0.166_r8, 0.102_r8, 6.61e-3_r8, 3.69e-3_r8, &
         1.22e-3_r8, 1.83e-3_r8, 2.61e-3_r8, 3.65e-3_r8, 4.55e-3_r8 /) ! kg/m2
   real(r8), dimension(nvar_sm,nvar_sm) :: C_sm = reshape([ &
         1.0_r8, 0.9_r8, 0.8_r8, 0.7_r8, 0.6_r8, 0.5_r8, 0.4_r8, 0.3_r8, 0.2_r8, 0.1_r8, &
         0.9_r8, 1.0_r8, 0.9_r8, 0.8_r8, 0.7_r8, 0.6_r8, 0.5_r8, 0.4_r8, 0.3_r8, 0.2_r8, &
         0.8_r8, 0.9_r8, 1.0_r8, 0.9_r8, 0.8_r8, 0.7_r8, 0.6_r8, 0.5_r8, 0.4_r8, 0.3_r8, &
         0.7_r8, 0.8_r8, 0.9_r8, 1.0_r8, 0.9_r8, 0.8_r8, 0.7_r8, 0.6_r8, 0.5_r8, 0.4_r8, &
         0.6_r8, 0.7_r8, 0.8_r8, 0.9_r8, 1.0_r8, 0.9_r8, 0.8_r8, 0.7_r8, 0.6_r8, 0.5_r8, &
         0.5_r8, 0.6_r8, 0.7_r8, 0.8_r8, 0.9_r8, 1.0_r8, 0.9_r8, 0.8_r8, 0.7_r8, 0.6_r8, &
         0.4_r8, 0.5_r8, 0.6_r8, 0.7_r8, 0.8_r8, 0.9_r8, 1.0_r8, 0.9_r8, 0.8_r8, 0.7_r8, &
         0.3_r8, 0.4_r8, 0.5_r8, 0.6_r8, 0.7_r8, 0.8_r8, 0.9_r8, 1.0_r8, 0.9_r8, 0.8_r8, &
         0.2_r8, 0.3_r8, 0.4_r8, 0.5_r8, 0.6_r8, 0.7_r8, 0.8_r8, 0.9_r8, 1.0_r8, 0.9_r8, &
         0.1_r8, 0.2_r8, 0.3_r8, 0.4_r8, 0.5_r8, 0.6_r8, 0.7_r8, 0.8_r8, 0.9_r8, 1.0_r8], &
         shape=[nvar_sm, nvar_sm])                                         ! cross-layer correlation
   real(r8), allocatable :: r_prev_sm(:,:,:)                               ! previous perturbation (numpatch, nvar_sm, DEF_DA_ENS_NUM)
   real(r8), allocatable :: r_curr_sm(:,:,:)                               ! current perturbation (numpatch, nvar_sm, DEF_DA_ENS_NUM)

   ! perturbed soil hydraulic parameters
   real(r8), allocatable :: porsl_ens (:,:,:)                              ! fraction of soil that is voids (nl_soil, 0:N, numpatch) [-]
   real(r8), allocatable :: hksati_ens(:,:,:)                              ! saturated hydraulic conductivity (nl_soil, 0:N, numpatch) [mm/s]

!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE ensemble(deltim)

!-----------------------------------------------------------------------
      IMPLICIT NONE

      real(r8), intent(in) :: deltim

!------------------------ Local Variables ------------------------------
      integer  ::  np, i, j

      real(r8) :: cov_matrix(nvar_forc, nvar_forc)                         ! covariance matrix between perturbed forcing variables
      real(r8) :: L(nvar_forc, nvar_forc)                                  ! Cholesky decomposition of correlation matrix
      integer  :: info                                                     ! info flag for Cholesky decomposition
      real(r8) :: u1(DEF_DA_ENS_NUM/2), u2(DEF_DA_ENS_NUM/2)               ! uniform random variables
      real(r8) :: z_forc(nvar_forc, DEF_DA_ENS_NUM)                        ! standard normal random variables
      real(r8) :: mean_z_forc(nvar_forc)                                   ! mean of perturbation (nvar_forc)
      real(r8) :: std_z_forc(nvar_forc)                                    ! std of perturbation (nvar_forc)
      real(r8) :: mean_r_forc(nvar_forc)                                  ! Ensemble means of forcing perturbations after AR(1)
      real(r8) :: zxL(nvar_forc, DEF_DA_ENS_NUM)                           ! correlated random variables (nvar, DEF_DA_ENS_NUM)

      real(r8) :: cov_matrix_sm(nvar_sm, nvar_sm)                          ! covariance matrix between soil moisture layers
      real(r8) :: L_sm(nvar_sm, nvar_sm)                                   ! Cholesky decomposition of soil moisture covariance
      real(r8) :: z_sm(nvar_sm, DEF_DA_ENS_NUM)                            ! standard normal random variables (layer, member)
      real(r8) :: mean_z_sm(nvar_sm)                                       ! mean of perturbation (nvar_sm)
      real(r8) :: std_z_sm(nvar_sm)                                        ! std of perturbation (nvar_sm)
      real(r8) :: zxL_sm(nvar_sm, DEF_DA_ENS_NUM)                          ! correlated random variables (nvar_sm, DEF_DA_ENS_NUM)
      real(r8) :: mean_r_sm(nvar_sm)                                       ! mean of perturbation for soil moisture (nvar_sm)
      real(r8) :: a_sm(nvar_sm, DEF_DA_ENS_NUM)                            ! disturbed soil moisture (layer, member)

!-----------------------------------------------------------------------
      IF (.not. p_is_worker) RETURN

      ! initialize persistent variables
      IF (.not. initialized) THEN
         CALL random_seed()
         allocate(r_prev_forc(numpatch, nvar_forc, DEF_DA_ENS_NUM))
         allocate(r_curr_forc(numpatch, nvar_forc, DEF_DA_ENS_NUM))
         allocate(r_prev_sm  (numpatch, nvar_sm,   DEF_DA_ENS_NUM))
         allocate(r_curr_sm  (numpatch, nvar_sm,   DEF_DA_ENS_NUM))
         r_prev_forc = 0.0_r8
         r_curr_forc = 0.0_r8
         r_prev_sm   = 0.0_r8
         r_curr_sm   = 0.0_r8
         initialized = .true.
      ENDIF

      ! calculate AR(1) parameters
      dt = deltim/3600
      phi_forc = exp(-dt/tau_forc)
      sigma_eps_forc = sqrt(1.0 - phi_forc**2)
      phi_sm = exp(-dt/tau_sm)
      sigma_eps_sm = sqrt(1.0 - phi_sm**2)

      ! calculate covariance matrix by cross-correlation matrix and standard deviation
      cov_matrix = 0.0_r8
      DO i = 1, nvar_forc
         DO j = 1, nvar_forc
            cov_matrix(i,j) = C(i,j) * sigma_forc(i) * sigma_forc(j)
         ENDDO
      ENDDO

      ! perform Cholesky decomposition of covariance matrix
      L = cov_matrix
      CALL dpotrf('L', nvar_forc, L, nvar_forc, info)
      DO i = 1, nvar_forc
         DO j = i+1, nvar_forc
            L(i,j) = 0.0_r8
         ENDDO
      ENDDO

      ! calculate covariance matrix and its Cholesky decomposition for soil moisture
      cov_matrix_sm = 0.0_r8
      DO i = 1, nvar_sm
         DO j = 1, nvar_sm
            cov_matrix_sm(i,j) = C_sm(i,j) * sigma_sm(i) * sigma_sm(j)
         ENDDO
      ENDDO
      L_sm = cov_matrix_sm
      CALL dpotrf('L', nvar_sm, L_sm, nvar_sm, info)
      DO i = 1, nvar_sm
         DO j = i+1, nvar_sm
            L_sm(i,j) = 0.0_r8
         ENDDO
      ENDDO

      ! Generate ensemble samples for forcing variables
      DO np = 1, numpatch
         ! generate disturbance ensemble samples ~ N(0, I)
         CALL random_number(u1)
         CALL random_number(u2)
         DO i = 1, DEF_DA_ENS_NUM/2
            u1(i) = max(u1(i), 1e-10)  ! ensure u1 is not zero
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

         ! normalize z to mean 0 and std 1
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

         ! introduce correlation using AR(1) process
         DO i = 1, nvar_forc
            DO j = 1, DEF_DA_ENS_NUM
               r_curr_forc(np,i,j) = phi_forc * r_prev_forc(np,i,j) + sigma_eps_forc * zxL(i,j)
            ENDDO
         ENDDO
         ! no AR(1) process, directly use correlated random variables
         ! r_curr_forc(np,:,:) = zxL

         ! normalize the disturbance ensemble samples to mean 0
         mean_r_forc = sum(r_curr_forc(np,:,:), dim=2)/DEF_DA_ENS_NUM
         DO i = 1, nvar_forc
            DO j = 1, DEF_DA_ENS_NUM
               r_curr_forc(np,i,j) = r_curr_forc(np,i,j) - mean_r_forc(i)
            ENDDO
         ENDDO

         ! save current perturbation as previous perturbation for next time step
         r_prev_forc(np,:,:) = r_curr_forc(np,:,:)

         ! generate ensemble samples according different types
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
            ! generate disturbance ensemble samples ~ N(0, I) for each soil layer
            DO i = 1, nvar_sm
               CALL random_number(u1)
               CALL random_number(u2)
               DO j = 1, DEF_DA_ENS_NUM/2
                  u1(j) = max(u1(j), 1e-10)  ! ensure u1 is not zero
                  z_sm(i,j*2-1) = sqrt(-2.0*log(u1(j))) * cos(2.0*pi*u2(j))
                  z_sm(i,j*2)   = sqrt(-2.0*log(u1(j))) * sin(2.0*pi*u2(j))
               ENDDO
            ENDDO

            ! normalize z to mean 0 and std 1
            DO i = 1, nvar_sm
               mean_z_sm(i) = sum(z_sm(i,:))/DEF_DA_ENS_NUM
               std_z_sm(i) = sqrt(sum((z_sm(i,:) - mean_z_sm(i))**2)/(DEF_DA_ENS_NUM - 1))
               z_sm(i,:) = (z_sm(i,:) - mean_z_sm(i))/std_z_sm(i)
            ENDDO

            ! multiply by Cholesky factor to introduce cross-layer correlation (L*z)
            CALL dgemm('N', 'N', nvar_sm, DEF_DA_ENS_NUM, nvar_sm, &
               1.0_r8, L_sm, nvar_sm, z_sm, nvar_sm, 0.0_r8, zxL_sm, nvar_sm)

            ! introduce correlation using AR(1) process
            DO i = 1, nvar_sm
               DO j = 1, DEF_DA_ENS_NUM
                  r_curr_sm(np,i,j) = phi_sm * r_prev_sm(np,i,j) + sigma_eps_sm * zxL_sm(i,j)
               ENDDO
            ENDDO

            ! normalize the disturbance ensemble samples to mean 0
            mean_r_sm = sum(r_curr_sm(np,:,:), dim=2)/DEF_DA_ENS_NUM
            DO i = 1, nvar_sm
               DO j = 1, DEF_DA_ENS_NUM
                  r_curr_sm(np,i,j) = r_curr_sm(np,i,j) - mean_r_sm(i)
               ENDDO
            ENDDO

            ! save current perturbation as previous perturbation for next time step
            r_prev_sm(np,:,:) = r_curr_sm(np,:,:)

            ! generate perturbed soil moisture samples
            DO i = 1, nvar_sm
               DO j = 1, DEF_DA_ENS_NUM
                  a_sm(i,j) = max(1.0e-10_r8, wliq_soisno_ens(i,j,np) + r_curr_sm(np,i,j))
               ENDDO
            ENDDO

            ! move residual water to water table
            DO j = 1, DEF_DA_ENS_NUM
               wa_ens(j, np) = wa_ens(j, np) - sum(wliq_soisno_ens(1:nvar_sm, j, np) - a_sm(:,j))
               wliq_soisno_ens(1:nvar_sm, j, np) = a_sm(:,j)
            ENDDO
         ENDIF
      ENDDO

   END SUBROUTINE ensemble

!-----------------------------------------------------------------------

   SUBROUTINE save_member(iens)

!-----------------------------------------------------------------------
! DESCRIPTION:
!    Save the current model states, forcing, and diagnostic and flux
!    outputs into ensemble member.
!-----------------------------------------------------------------------
      IMPLICIT NONE
      integer, intent(in) :: iens

      IF (p_is_worker) THEN
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
      ENDIF

   END SUBROUTINE save_member

!-----------------------------------------------------------------------------

   SUBROUTINE load_member(iens)

!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Load the states, forcing, and flux outputs of an ensemble member
!    into the model variables.
!-----------------------------------------------------------------------------
      IMPLICIT NONE
      integer, intent(in) :: iens

      IF (p_is_worker) THEN
         ! Soil hydraulic parameters
         porsl        = porsl_ens       (  :,iens,:)
         hksati       = hksati_ens      (  :,iens,:)

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
      ENDIF

   END SUBROUTINE load_member

!-----------------------------------------------------------------------------

   SUBROUTINE postprocess()

!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Set the final model states and outputs after data assimilation.
!
!    When land-atmosphere feedback is enabled, apply the analysis
!    ensemble means for restart states and fluxes.
!
!    Otherwise, restore the open-loop trajectory for diagnostic.
!-----------------------------------------------------------------------------
      IMPLICIT NONE

      IF (p_is_worker) THEN
         CALL load_member(0)

         IF (DEF_DA_FEEDBACK) THEN
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
      ENDIF

   END SUBROUTINE postprocess

!-----------------------------------------------------------------------------

   SUBROUTINE init_soil_ptb()

!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Initialize the soil hydraulic parameter perturbation.
!
!    The standard normal draws are generated by mkinidata and stored in the
!    const file. Here they are only read, normalized across the members in
!    use, and expanded into one parameter set per member: no random number
!    is generated at run time.
!
!    Each member applies one global adjustment, shared by all patches and
!    layers, which represents the overall mean bias of the soil parameter
!    dataset.
!
!    Only porsl and hksati are perturbed. The standard deviations follow
!    Ewerdwalbesloh et al. (2026), Eqs. (6) and (8), which are the Cosby et al.
!    (1984) pedotransfer uncertainty regressions: sigma_hksati is the standard
!    deviation of log10(hksati), sigma_porsl is the standard deviation of
!    porsl itself. Both are damped by DEF_DA_ENS_SOIL_PTB_SCALE.
!-----------------------------------------------------------------------------
      IMPLICIT NONE

      integer  :: i, iens, ipatch, il
      real(r8) :: mean_z, std_z
      real(r8) :: z_porsl, z_hksati
      real(r8) :: sand_pct, clay_pct, silt_pct
      real(r8) :: sigma_porsl, sigma_hksati

      IF (.not. DEF_DA .or. DEF_DA_ENS_NUM < 1 .or. &
          .not. p_is_worker .or. numpatch <= 0) RETURN

      ! porsl_ens is used by the observation operator and the SM analysis even
      ! when the parameter perturbation is switched off, so allocate it here
      ! for every DA run; member 0 always keeps the original parameter.
      IF (.not. allocated(porsl_ens)) THEN
         allocate (porsl_ens  (nl_soil,0:DEF_DA_ENS_NUM,numpatch))
         allocate (hksati_ens (nl_soil,0:DEF_DA_ENS_NUM,numpatch))
         porsl_ens  = spread(porsl,  2, DEF_DA_ENS_NUM + 1)
         hksati_ens = spread(hksati, 2, DEF_DA_ENS_NUM + 1)
      ENDIF

      IF (.not. DEF_DA_ENS_SOIL_PTB) RETURN

      ! normalize each parameter across the members in use to mean 0 and std 1
      ! z_soil_ptb row 1: porsl, row 2: hksati
      DO i = 1, 2
         mean_z = sum(z_soil_ptb(i,1:DEF_DA_ENS_NUM))/DEF_DA_ENS_NUM
         std_z  = sqrt(sum((z_soil_ptb(i,1:DEF_DA_ENS_NUM) - mean_z)**2)/(DEF_DA_ENS_NUM - 1))
         std_z  = max(std_z, 1.0e-10_r8)
         z_soil_ptb(i,1:DEF_DA_ENS_NUM) = (z_soil_ptb(i,1:DEF_DA_ENS_NUM) - mean_z)/std_z
      ENDDO

      ! expand each member's global adjustment over all patches and layers
      DO iens = 1, DEF_DA_ENS_NUM
         z_porsl  = z_soil_ptb(1, iens)
         z_hksati = z_soil_ptb(2, iens)

         DO ipatch = 1, numpatch
            IF (patchtype(ipatch) > 2) CYCLE
            DO il = 1, nl_soil
               sand_pct = vf_sand(il,ipatch) * 100.0_r8
               clay_pct = vf_clay(il,ipatch) * 100.0_r8
               silt_pct = 100.0_r8 - sand_pct - clay_pct

               sigma_porsl  = DEF_DA_ENS_SOIL_PTB_SCALE * (7.73_r8 - 0.073_r8*clay_pct)/100.0_r8
               sigma_hksati = DEF_DA_ENS_SOIL_PTB_SCALE * (0.459_r8 + 0.00321_r8*silt_pct)

               ! porosity: additive in raw space
               porsl_ens(il,iens,ipatch) = max(0.01_r8, min(0.99_r8, porsl_ens(il,0,ipatch) + sigma_porsl * z_porsl))

               ! saturated hydraulic conductivity: perturbation of log10
               hksati_ens(il,iens,ipatch) = 10.0_r8**(log10(hksati_ens(il,0,ipatch)) + sigma_hksati * z_hksati)
               hksati_ens(il,iens,ipatch) = max(hksati_ens(il,iens,ipatch), 1.0e-10_r8)
            ENDDO
         ENDDO
      ENDDO

   END SUBROUTINE init_soil_ptb

!-----------------------------------------------------------------------------

   SUBROUTINE end_DA_Ens()

!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Release the ensemble perturbation buffers.
!-----------------------------------------------------------------------------
      IMPLICIT NONE

      IF (allocated(r_prev_forc))   deallocate(r_prev_forc)
      IF (allocated(r_curr_forc))   deallocate(r_curr_forc)
      IF (allocated(r_prev_sm  ))   deallocate(r_prev_sm  )
      IF (allocated(r_curr_sm  ))   deallocate(r_curr_sm  )
      IF (allocated(porsl_ens  ))   deallocate(porsl_ens  )
      IF (allocated(hksati_ens ))   deallocate(hksati_ens )
      initialized = .false.

   END SUBROUTINE end_DA_Ens

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Ens
#endif
