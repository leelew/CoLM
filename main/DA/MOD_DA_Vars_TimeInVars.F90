#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Vars_TimeInVars
!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Store the time-invariant soil hydraulic parameters for ensemble members.
!   Member 0 retains the original parameters; perturbed members remain fixed
!   throughout the integration.
!
!   Generate separate standard normal samples for porsl and hksati when writing
!   the const file. Read and broadcast these samples during initialization;
!   MOD_DA_Ens applies one adjustment per parameter and member to all soil
!   patches and layers on the first ensemble forecast.
!
! HISTORY:
!   Zhilong Fan, 11/2026: Initial version
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_Vars_Global
   USE MOD_SPMD_Task
   USE MOD_LandPatch, only: numpatch
   USE MOD_NetCDFSerial, only: ncio_define_dimension, &
      ncio_write_serial_real8_1d, ncio_read_bcast_serial_real8_1d
   IMPLICIT NONE
   SAVE

   ! Ensemble soil hydraulic parameters
   real(r8), allocatable            :: porsl_ens (:,:,:)    ! Soil porosity [-]
   real(r8), allocatable            :: hksati_ens(:,:,:)    ! Saturated hydraulic conductivity [mm/s]

   ! Ensemble soil hydraulic parameters
   real(r8), allocatable            :: z_porsl   (:)        ! Porosity adjustment per member [-]
   real(r8), allocatable            :: z_hksati  (:)        ! Log10 conductivity adjustment per member [-]

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE allocate_DA_TimeInVars()

!-----------------------------------------------------------------------------
      IF (p_is_worker .and. numpatch > 0) THEN
         allocate(porsl_ens (nl_soil,0:DEF_DA_ENS_NUM,numpatch))
         allocate(hksati_ens(nl_soil,0:DEF_DA_ENS_NUM,numpatch))
      ENDIF

   END SUBROUTINE allocate_DA_TimeInVars

!-----------------------------------------------------------------------------

   SUBROUTINE deallocate_DA_TimeInVars()

!-----------------------------------------------------------------------------
      IF (allocated(porsl_ens )) deallocate(porsl_ens )
      IF (allocated(hksati_ens)) deallocate(hksati_ens)
      IF (allocated(z_porsl   )) deallocate(z_porsl   )
      IF (allocated(z_hksati  )) deallocate(z_hksati  )

   END SUBROUTINE deallocate_DA_TimeInVars

!-----------------------------------------------------------------------------

   SUBROUTINE WRITE_DA_TimeInVars(file_restart)

!-----------------------------------------------------------------------------
      character(len=*), intent(in)  :: file_restart

      integer                       :: i, k
      real(r8)                      :: u1(DEF_DA_ENS_NUM/2), u2(DEF_DA_ENS_NUM/2)
      real(r8)                      :: z(DEF_DA_ENS_NUM,2)
      real(r8)                      :: mean_z, std_z
!-----------------------------------------------------------------------------
      IF (.not. DEF_DA_ENS_SOILPARAM) RETURN

      CALL random_seed()
      DO k = 1, 2
         CALL random_number(u1)
         CALL random_number(u2)
            DO i = 1, DEF_DA_ENS_NUM/2
               u1(i) = max(u1(i), 1.0e-10_r8)
            z(2*i-1,k) = sqrt(-2.0_r8*log(u1(i))) * cos(2.0_r8*pi*u2(i))
            z(2*i,  k) = sqrt(-2.0_r8*log(u1(i))) * sin(2.0_r8*pi*u2(i))
         ENDDO
         mean_z = sum(z(:,k))/DEF_DA_ENS_NUM
         std_z  = sqrt(sum((z(:,k)-mean_z)**2)/(DEF_DA_ENS_NUM-1))
         z(:,k) = (z(:,k)-mean_z) / max(std_z, 1.0e-10_r8)
      ENDDO

      CALL ncio_define_dimension(file_restart, 'ens', DEF_DA_ENS_NUM)
      CALL ncio_write_serial_real8_1d(file_restart, 'z_porsl',  z(:,1), 'ens')
      CALL ncio_write_serial_real8_1d(file_restart, 'z_hksati', z(:,2), 'ens')

   END SUBROUTINE WRITE_DA_TimeInVars

!-----------------------------------------------------------------------------

   SUBROUTINE READ_DA_TimeInVars(file_restart)

!-----------------------------------------------------------------------------
      character(len=*), intent(in)  :: file_restart
!-----------------------------------------------------------------------------
      IF (.not. DEF_DA_ENS_SOILPARAM) RETURN

      CALL ncio_read_bcast_serial_real8_1d(file_restart, 'z_porsl',  z_porsl)
      CALL ncio_read_bcast_serial_real8_1d(file_restart, 'z_hksati', z_hksati)

   END SUBROUTINE READ_DA_TimeInVars

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Vars_TimeInVars
#endif
