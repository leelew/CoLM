#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Data_ANATQ
!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Analysis-derived 2m temperature and relative humidity pseudo-observation
!   handling for land data assimilation.
!
! WORKFLOW:
!   1. Forward: Select patches with changed forecast and analysis ensemble
!      means independently for each variable. Use analysis ensemble means as
!      pseudo-observations and forecast members as matching H(x). Gather
!      observations through PseudoType, then gather and broadcast H(x) in the
!      same record order. Combine temperature records first, followed by
!      relative humidity records.
!   2. Clear: Release current pseudo-observations and H(x).
!
! SUPPORTED:
!   2MT:
!      Pseudo-observations from 2 m temperature analyses [K].
!      The observation error standard deviation is 1 K.
!   2MRH:
!      Pseudo-observations from 2 m relative humidity analyses [0-1].
!      The observation error standard deviation is 0.04.
!
! HISTORY:
!   Lu Li, 11/2026: DA v1.0
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_LandPatch
   USE MOD_SPMD_Task
   USE MOD_Vars_Global
   USE MOD_Vars_TimeInvariants
   USE MOD_DA_Data_PseudoType
   USE MOD_DA_Assim_TQ
   IMPLICIT NONE

!-----------------------------------------------------------------------------
! Analysis-derived pseudo-observation derived type
!-----------------------------------------------------------------------------
   type :: ANATQ_type

      type(pseudo_obs_type)                :: y                             ! Pseudo-observations and patch locations
      real(r8), allocatable                :: hx_obs(:,:)                   ! Observation-space forecast ensemble H(x)

   CONTAINS

      procedure                            :: forward => ANATQ_forward      ! Construct y and matching observation-space H(x)
      procedure                            :: clear   => ANATQ_clear        ! Clear current pseudo-observations and H(x)

   END type ANATQ_type

!-----------------------------------------------------------------------------
! Collection of temperature and relative humidity pseudo-observations
!-----------------------------------------------------------------------------
   type :: ANATQ_set_type

      type(ANATQ_type)                     :: tref                          ! 2 m temperature pseudo-observations [K]
      type(ANATQ_type)                     :: rhref                         ! 2 m relative humidity pseudo-observations [0-1]

      integer                              :: nobs = 0                      ! Total number of temperature and relative humidity records
      real(r8), allocatable                :: lat(:)                        ! Observation latitudes [deg]
      real(r8), allocatable                :: lon(:)                        ! Observation longitudes [deg]
      real(r8), allocatable                :: y(:)                          ! Analysis ensemble means of the selected variable
      real(r8), allocatable                :: sigma(:)                      ! Observation error standard deviations
      real(r8), allocatable                :: hx_obs(:,:)                   ! Observation-space forecast ensemble H(x)

   CONTAINS

      procedure                            :: forward => ANATQ_set_forward  ! Prepare both variables and combine their records
      procedure                            :: concat  => ANATQ_set_concat   ! Combine temperature records followed by humidity
      procedure                            :: clear   => ANATQ_set_clear    ! Clear both variables and their combined arrays

   END type ANATQ_set_type

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE ANATQ_forward(this, var_name)

!-----------------------------------------------------------------------------
      class(ANATQ_type), intent(inout)     :: this
      character(len=*), intent(in)         :: var_name

      real(r8)                             :: sigma                         ! Observation error standard deviation
      integer                              :: iens
      integer                              :: iworker
      integer                              :: nlocal                        ! Number of selected patches on one worker
      integer                              :: nobs
      integer                              :: offset                        ! Number of records already received on master
      logical,  allocatable                :: mask(:)                       ! Local patches with valid, changed analysis means
      real(r8), allocatable                :: y(:)                          ! Analysis means on local patches
      real(r8), allocatable                :: hx_obs(:,:)                   ! Forecast H(x) at selected patches on each worker
      real(r8), allocatable                :: hx_recv(:,:)                  ! Forecast H(x) gathered on master
!-----------------------------------------------------------------------------
      CALL this%clear()

      ! Select valid patches with changed analysis ensemble means.
      IF (p_is_worker) THEN
         allocate(mask(numpatch), y(numpatch))
         mask = .false.

         SELECT CASE (trim(var_name))
         CASE ('2MT')
            sigma = 1.0_r8
            y = sum(tref_ens_a(1:DEF_DA_ENS_NUM,:), dim=1) / DEF_DA_ENS_NUM
            mask = y /= sum(tref_ens_f(1:DEF_DA_ENS_NUM,:), dim=1) / DEF_DA_ENS_NUM
         CASE ('2MRH')
            sigma = 0.04_r8
            y = sum(rhref_ens_a(1:DEF_DA_ENS_NUM,:), dim=1) / DEF_DA_ENS_NUM
            mask = y /= sum(rhref_ens_f(1:DEF_DA_ENS_NUM,:), dim=1) / DEF_DA_ENS_NUM
         END SELECT
      ENDIF

      ! Gather observations through PseudoType; all processes participate.
      CALL this%y%read(mask, y, sigma)
      nobs = this%y%nobs
      IF (nobs == 0) RETURN

      ! Workers send H(x) in the same local patch order used for observations.
      IF (p_is_worker) THEN
         nlocal = count(mask)
         CALL mpi_send(nlocal, 1, MPI_INTEGER, p_address_master, mpi_tag_size, p_comm_glb, p_err)
         IF (nlocal > 0) THEN
            allocate(hx_obs(DEF_DA_ENS_NUM,nlocal))
            SELECT CASE (trim(var_name))
            CASE ('2MT')
               DO iens = 1, DEF_DA_ENS_NUM
                  hx_obs(iens,:) = pack(tref_ens_f(iens,:), mask)
               ENDDO
            CASE ('2MRH')
               DO iens = 1, DEF_DA_ENS_NUM
                  hx_obs(iens,:) = pack(rhref_ens_f(iens,:), mask)
               ENDDO
            END SELECT
            CALL mpi_send(hx_obs, DEF_DA_ENS_NUM*nlocal, MPI_REAL8, &
               p_address_master, mpi_tag_data, p_comm_glb, p_err)
         ENDIF
      ENDIF

      ! Master receives H(x) in the same worker order used by PseudoType.
      IF (p_is_master) THEN
         allocate(hx_recv(DEF_DA_ENS_NUM,nobs))
         offset = 0
         DO iworker = 0, p_np_worker-1
            CALL mpi_recv(nlocal, 1, MPI_INTEGER, p_address_worker(iworker), &
               mpi_tag_size, p_comm_glb, p_stat, p_err)
            IF (nlocal > 0) THEN
               CALL mpi_recv(hx_recv(:,offset+1:offset+nlocal), DEF_DA_ENS_NUM*nlocal, MPI_REAL8, &
                  p_address_worker(iworker), mpi_tag_data, p_comm_glb, p_stat, p_err)
               offset = offset + nlocal
            ENDIF
         ENDDO
      ENDIF

      ! Broadcast complete observation-space H(x) to all processes.
      allocate(this%hx_obs(nobs,DEF_DA_ENS_NUM))
      IF (p_is_master) this%hx_obs = transpose(hx_recv)
      CALL mpi_bcast(this%hx_obs, nobs*DEF_DA_ENS_NUM, MPI_REAL8, &
         p_address_master, p_comm_glb, p_err)

   END SUBROUTINE ANATQ_forward

!-----------------------------------------------------------------------------

   SUBROUTINE ANATQ_clear(this)

!-----------------------------------------------------------------------------
      class(ANATQ_type), intent(inout)     :: this
!-----------------------------------------------------------------------------
      CALL this%y%clear()
      IF (allocated(this%hx_obs)) deallocate(this%hx_obs)

   END SUBROUTINE ANATQ_clear

!-----------------------------------------------------------------------------

   SUBROUTINE ANATQ_set_forward(this)

!-----------------------------------------------------------------------------
      class(ANATQ_set_type), intent(inout) :: this
!-----------------------------------------------------------------------------
      CALL this%clear()
      CALL this%tref%forward('2MT')
      CALL this%rhref%forward('2MRH')
      CALL this%concat()

   END SUBROUTINE ANATQ_set_forward

!-----------------------------------------------------------------------------

   SUBROUTINE ANATQ_set_concat(this)

!-----------------------------------------------------------------------------
      class(ANATQ_set_type), intent(inout) :: this
      integer                              :: nt, nq
!-----------------------------------------------------------------------------
      nt = this%tref%y%nobs
      nq = this%rhref%y%nobs
      this%nobs = nt + nq
      IF (this%nobs == 0) RETURN

      allocate(this%lat(this%nobs), this%lon(this%nobs), this%y(this%nobs), this%sigma(this%nobs))
      allocate(this%hx_obs(this%nobs,DEF_DA_ENS_NUM))
      IF (nt > 0) THEN
         this%lat(1:nt) = this%tref%y%lat
         this%lon(1:nt) = this%tref%y%lon
         this%y(1:nt) = this%tref%y%y
         this%sigma(1:nt) = this%tref%y%sigma
         this%hx_obs(1:nt,:) = this%tref%hx_obs
      ENDIF
      IF (nq > 0) THEN
         this%lat(nt+1:nt+nq) = this%rhref%y%lat
         this%lon(nt+1:nt+nq) = this%rhref%y%lon
         this%y(nt+1:nt+nq) = this%rhref%y%y
         this%sigma(nt+1:nt+nq) = this%rhref%y%sigma
         this%hx_obs(nt+1:nt+nq,:) = this%rhref%hx_obs
      ENDIF

   END SUBROUTINE ANATQ_set_concat

!-----------------------------------------------------------------------------

   SUBROUTINE ANATQ_set_clear(this)

!-----------------------------------------------------------------------------
      class(ANATQ_set_type), intent(inout) :: this
!-----------------------------------------------------------------------------
      CALL this%tref%clear()
      CALL this%rhref%clear()
      this%nobs = 0
      IF (allocated(this%lat   )) deallocate(this%lat   )
      IF (allocated(this%lon   )) deallocate(this%lon   )
      IF (allocated(this%y     )) deallocate(this%y     )
      IF (allocated(this%sigma )) deallocate(this%sigma )
      IF (allocated(this%hx_obs)) deallocate(this%hx_obs)

   END SUBROUTINE ANATQ_set_clear

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Data_ANATQ
#endif
