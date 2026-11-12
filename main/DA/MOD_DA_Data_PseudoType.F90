#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Data_PseudoType
!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Handle one patch-based pseudo-observation variable derived from model
!   analyses.
!
! NOTES:
!   In operational forecasting systems, screen-level temperature and humidity
!   analyses can provide observational constraints on land-surface variables,
!   supporting land-atmosphere coupled data assimilation.
!
!   Following SYNOP 2 m temperature and humidity assimilation, patch-level
!   analysis values with nonzero analysis increments are treated as
!   pseudo-observations for subsequent land-surface data assimilation.
!
! WORKFLOW:
!   1. Read: Pack selected records on each worker, gather them on master,
!      and broadcast the combined observation vectors to all processes.
!   2. Clear: Release the observation vectors and reset the record count.
!
! HISTORY:
!   Lu Li, 11/2026: DA v1.0
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_LandPatch
   USE MOD_SPMD_Task
   USE MOD_Vars_Global
   USE MOD_Vars_TimeInvariants
   IMPLICIT NONE

!-----------------------------------------------------------------------------
! Pseudo-observation derived type
!-----------------------------------------------------------------------------
   type :: pseudo_obs_type

      integer :: nobs = 0                                       ! Number of selected patches
      real(r8), allocatable :: lat(:)                           ! Patch latitude
      real(r8), allocatable :: lon(:)                           ! Patch longitude
      real(r8), allocatable :: y(:)                             ! Analysis ensemble means of the pseudo-observation variable
      real(r8), allocatable :: r(:)                             ! Observation error standard deviations

   CONTAINS

      procedure, PUBLIC :: read  => pseudo_obs_read             ! Gather and broadcast selected patch observations
      procedure, PUBLIC :: clear => pseudo_obs_clear            ! Deallocate observation arrays and reset the record count

   END type pseudo_obs_type

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE pseudo_obs_read(this, mask, y, r)

!-----------------------------------------------------------------------------
      class(pseudo_obs_type), intent(inout) :: this     ! Pseudo-observation collection
      logical, allocatable, intent(in)      :: mask(:)  ! Selected local patches; allocated only on workers
      real(r8), allocatable, intent(in)     :: y(:)     ! Analysis means on local patches; allocated only on workers
      real(r8), intent(in)                  :: r        ! Observation error standard deviation

      integer               :: np               ! Local patch index
      integer               :: nlocal           ! Number of selected patches on one worker
      integer               :: iobs             ! Record index in the worker send buffer
      integer               :: iworker          ! Worker index in the global worker list
      integer               :: offset           ! Number of records already received on master
      real(r8), allocatable :: obs(:,:)         ! Selected latitude, longitude, y and r on each worker
      real(r8), allocatable :: obs_global(:,:)  ! Observation records gathered on master
!-----------------------------------------------------------------------------
      CALL this%clear()

      ! Count selected patches on workers and obtain the global record count.
      nlocal = 0
      IF (p_is_worker) nlocal = count(mask)
      CALL mpi_allreduce(nlocal, this%nobs, 1, MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
      IF (this%nobs == 0) RETURN

      ! Workers send selected records in local patch order.
      IF (p_is_worker) THEN
         CALL mpi_send(nlocal, 1, MPI_INTEGER, p_address_master, mpi_tag_size, p_comm_glb, p_err)
         IF (nlocal > 0) THEN
            allocate(obs(4,nlocal))
            iobs = 0
            DO np = 1, numpatch
               IF (.not. mask(np)) CYCLE
               iobs = iobs + 1
               obs(1,iobs) = patchlatr(np)*180.0_r8/pi
               obs(2,iobs) = (modulo(patchlonr(np)+pi,2.0_r8*pi)-pi)*180.0_r8/pi
               obs(3,iobs) = y(np)
               obs(4,iobs) = r
            ENDDO
            CALL mpi_send(obs, 4*nlocal, MPI_REAL8, p_address_master, mpi_tag_data, p_comm_glb, p_err)
         ENDIF
      ENDIF

      ! Master concatenates records in worker order, matching ANATQ H(x).
      IF (p_is_master) THEN
         allocate(obs_global(4,this%nobs))
         offset = 0
         DO iworker = 0, p_np_worker-1
            CALL mpi_recv(nlocal, 1, MPI_INTEGER, p_address_worker(iworker), &
               mpi_tag_size, p_comm_glb, p_stat, p_err)
            IF (nlocal > 0) THEN
               CALL mpi_recv(obs_global(:,offset+1:offset+nlocal), 4*nlocal, MPI_REAL8, &
                  p_address_worker(iworker), mpi_tag_data, p_comm_glb, p_stat, p_err)
               offset = offset + nlocal
            ENDIF
         ENDDO
      ENDIF

      ! Broadcast complete observation vectors to all processes.
      allocate(this%lat(this%nobs), this%lon(this%nobs), this%y(this%nobs), this%r(this%nobs))
      IF (p_is_master) THEN
         this%lat = obs_global(1,:)
         this%lon = obs_global(2,:)
         this%y   = obs_global(3,:)
         this%r   = obs_global(4,:)
      ENDIF
      CALL mpi_bcast(this%lat, this%nobs, MPI_REAL8, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast(this%lon, this%nobs, MPI_REAL8, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast(this%y, this%nobs, MPI_REAL8, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast(this%r, this%nobs, MPI_REAL8, p_address_master, p_comm_glb, p_err)

   END SUBROUTINE pseudo_obs_read

!-----------------------------------------------------------------------------

   SUBROUTINE pseudo_obs_clear (this)

!-----------------------------------------------------------------------------
      class(pseudo_obs_type), intent(inout) :: this
!-----------------------------------------------------------------------------
      this%nobs = 0
      IF (allocated(this%lat)) deallocate(this%lat)
      IF (allocated(this%lon)) deallocate(this%lon)
      IF (allocated(this%y  )) deallocate(this%y  )
      IF (allocated(this%r  )) deallocate(this%r  )

   END SUBROUTINE pseudo_obs_clear

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Data_PseudoType
#endif

