#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Vars_1DFluxes
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Stores ensemble surface fluxes and land-atmosphere coupling variables
!    specifically for land-atmosphere coupling data assimilation.
!
! HISTORY:
!    Lu Li, 10/2026
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_Vars_Global
   USE MOD_SPMD_Task
   USE MOD_LandPatch
   IMPLICIT NONE
   SAVE

   PUBLIC :: allocate_1D_DAFluxes
   PUBLIC :: deallocate_1D_DAFluxes

!-----------------------------------------------------------------------------
! Surface flux variables
!-----------------------------------------------------------------------------
   real(r8), allocatable :: fsena_ens        (:,:)      ! sensible heat flux to the atmosphere                     [W/m2]
   real(r8), allocatable :: lfevpa_ens       (:,:)      ! latent heat flux to the atmosphere                       [W/m2]
   real(r8), allocatable :: fevpa_ens        (:,:)      ! evapotranspiration to the atmosphere                     [mm/s]
   real(r8), allocatable :: rsur_ens         (:,:)      ! surface runoff                                           [mm/s]
   real(r8), allocatable :: rnof_ens         (:,:)      ! total runoff                                             [mm/s]

!-----------------------------------------------------------------------------
! Land-atmosphere coupling variables
!-----------------------------------------------------------------------------
   real(r8), allocatable :: trad_ens         (:,:)      ! surface radiative temperature                            [K]
   real(r8), allocatable :: tref_ens         (:,:)      ! 2 m air temperature                                      [K]
   real(r8), allocatable :: tref_ens_f       (:,:)      ! prior 2 m air temperature                                [K]
   real(r8), allocatable :: tref_ens_a       (:,:)      ! posterior 2 m air temperature                            [K]
   real(r8), allocatable :: tref_f           (:)        ! ensemble-mean prior 2 m air temperature                  [K]
   real(r8), allocatable :: tref_a           (:)        ! ensemble-mean posterior 2 m air temperature              [K]
   real(r8), allocatable :: qref_ens         (:,:)      ! 2 m specific humidity                                    [kg/kg]
   real(r8), allocatable :: qref_ens_f       (:,:)      ! prior 2 m specific humidity                              [kg/kg]
   real(r8), allocatable :: qref_ens_a       (:,:)      ! posterior 2 m specific humidity                          [kg/kg]
   real(r8), allocatable :: qref_f           (:)        ! ensemble-mean prior 2 m specific humidity                [kg/kg]
   real(r8), allocatable :: qref_a           (:)        ! ensemble-mean posterior 2 m specific humidity            [kg/kg]
   real(r8), allocatable :: rhref_ens        (:,:)      ! 2 m relative humidity                                    [-]
   real(r8), allocatable :: rhref_ens_f      (:,:)      ! prior 2 m relative humidity                              [-]
   real(r8), allocatable :: rhref_ens_a      (:,:)      ! posterior 2 m relative humidity                          [-]
   real(r8), allocatable :: ustar_ens        (:,:)      ! friction velocity in similarity theory                   [m/s]
   real(r8), allocatable :: qstar_ens        (:,:)      ! humidity scale in similarity theory                      [kg/kg]
   real(r8), allocatable :: tstar_ens        (:,:)      ! temperature scale in similarity theory                   [K]
   real(r8), allocatable :: fm_ens           (:,:)      ! integrated momentum profile function                     [-]
   real(r8), allocatable :: fh_ens           (:,:)      ! integrated heat profile function                         [-]
   real(r8), allocatable :: fq_ens           (:,:)      ! integrated moisture profile function                     [-]

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE allocate_1D_DAFluxes()

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            allocate (fsena_ens        (0:DEF_DA_ENS_NUM,numpatch)); fsena_ens        (:,:) = spval
            allocate (lfevpa_ens       (0:DEF_DA_ENS_NUM,numpatch)); lfevpa_ens       (:,:) = spval
            allocate (fevpa_ens        (0:DEF_DA_ENS_NUM,numpatch)); fevpa_ens        (:,:) = spval
            allocate (rsur_ens         (0:DEF_DA_ENS_NUM,numpatch)); rsur_ens         (:,:) = spval
            allocate (rnof_ens         (0:DEF_DA_ENS_NUM,numpatch)); rnof_ens         (:,:) = spval
            allocate (trad_ens         (0:DEF_DA_ENS_NUM,numpatch)); trad_ens         (:,:) = spval
            allocate (tref_ens         (0:DEF_DA_ENS_NUM,numpatch)); tref_ens         (:,:) = spval
            allocate (tref_ens_f       (0:DEF_DA_ENS_NUM,numpatch)); tref_ens_f       (:,:) = spval
            allocate (tref_ens_a       (0:DEF_DA_ENS_NUM,numpatch)); tref_ens_a       (:,:) = spval
            allocate (tref_f           (numpatch)); tref_f           (:) = spval
            allocate (tref_a           (numpatch)); tref_a           (:) = spval
            allocate (qref_ens         (0:DEF_DA_ENS_NUM,numpatch)); qref_ens         (:,:) = spval
            allocate (qref_ens_f       (0:DEF_DA_ENS_NUM,numpatch)); qref_ens_f       (:,:) = spval
            allocate (qref_ens_a       (0:DEF_DA_ENS_NUM,numpatch)); qref_ens_a       (:,:) = spval
            allocate (qref_f           (numpatch)); qref_f           (:) = spval
            allocate (qref_a           (numpatch)); qref_a           (:) = spval
            allocate (rhref_ens        (0:DEF_DA_ENS_NUM,numpatch)); rhref_ens        (:,:) = spval
            allocate (rhref_ens_f      (0:DEF_DA_ENS_NUM,numpatch)); rhref_ens_f      (:,:) = spval
            allocate (rhref_ens_a      (0:DEF_DA_ENS_NUM,numpatch)); rhref_ens_a      (:,:) = spval
            allocate (ustar_ens        (0:DEF_DA_ENS_NUM,numpatch)); ustar_ens        (:,:) = spval
            allocate (qstar_ens        (0:DEF_DA_ENS_NUM,numpatch)); qstar_ens        (:,:) = spval
            allocate (tstar_ens        (0:DEF_DA_ENS_NUM,numpatch)); tstar_ens        (:,:) = spval
            allocate (fm_ens           (0:DEF_DA_ENS_NUM,numpatch)); fm_ens           (:,:) = spval
            allocate (fh_ens           (0:DEF_DA_ENS_NUM,numpatch)); fh_ens           (:,:) = spval
            allocate (fq_ens           (0:DEF_DA_ENS_NUM,numpatch)); fq_ens           (:,:) = spval
         ENDIF
      ENDIF

   END SUBROUTINE allocate_1D_DAFluxes

!-----------------------------------------------------------------------------

   SUBROUTINE deallocate_1D_DAFluxes()

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            deallocate (fsena_ens       )
            deallocate (lfevpa_ens      )
            deallocate (fevpa_ens       )
            deallocate (rsur_ens        )
            deallocate (rnof_ens        )
            deallocate (trad_ens        )
            deallocate (tref_ens        )
            deallocate (tref_ens_f      )
            deallocate (tref_ens_a      )
            deallocate (tref_f          )
            deallocate (tref_a          )
            deallocate (qref_ens        )
            deallocate (qref_ens_f      )
            deallocate (qref_ens_a      )
            deallocate (qref_f          )
            deallocate (qref_a          )
            deallocate (rhref_ens       )
            deallocate (rhref_ens_f     )
            deallocate (rhref_ens_a     )
            deallocate (ustar_ens       )
            deallocate (qstar_ens       )
            deallocate (tstar_ens       )
            deallocate (fm_ens          )
            deallocate (fh_ens          )
            deallocate (fq_ens          )
         ENDIF
      ENDIF

   END SUBROUTINE deallocate_1D_DAFluxes

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Vars_1DFluxes
#endif
