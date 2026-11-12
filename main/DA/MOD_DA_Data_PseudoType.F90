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

   SUBROUTINE pseudo_obs_read (this, mask, y, r)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------------ Dummy Arguments ------------------------------
      class(pseudo_obs_type), intent(inout) :: this             ! Pseudo-observations collections
      logical, intent(in)                   :: mask(:)          ! Selected local patches 
      real(r8), intent(in)                  :: y(:)             ! Analysis means at patch observation locations
      real(r8), intent(in)                  :: r                ! Observation error standard deviation

!------------------------------ Local Variables ------------------------------
      integer                               :: np               ! Local patch index
      integer                               :: nlocal           ! Number of selected patches on this rank
      integer                               :: iobs             ! Selected-patch index in send buffers
      integer                               :: iproc            ! Rank index in p_comm_glb
      real(r8), allocatable                 :: obs(:,:)         ! Send buffer: lat, lon, y, r 
      real(r8), allocatable                 :: obs_global(:,:)  ! Gathered lat, lon, y, r on master 
      integer, allocatable                  :: counts(:)        ! Selected-patch count per rank 
      integer, allocatable                  :: displs(:)        ! Receive offset per rank, in patch units 

!-----------------------------------------------------------------------------
      CALL this%clear()

      ! Pack selected patches in local patch order.
      IF (p_is_worker) THEN 
         nlocal = count(mask)
      ELSE
         nlocal = 0
      allocate(obs(4,nlocal))

      IF (p_is_worker) THEN
         iobs = 0
         DO np = 1, numpatch
            IF (.not. mask(np)) CYCLE
            iobs = iobs + 1
            obs(1,iobs) = patchlatr(np)*180.0_r8/pi
            obs(2,iobs) = (modulo(patchlonr(np)+pi,2.0_r8*pi)-pi)*180.0_r8/pi
            obs(3,iobs) = y(np)
            obs(4,iobs) = r
         ENDDO
      ENDIF

#ifdef USEMPI
      ! Gather selected-patch counts on master, then broadcast the total count.
      allocate(counts(0:p_np_glb-1), displs(0:p_np_glb-1))
      counts = 0
      displs = 0
      CALL mpi_gather(nlocal, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
      IF (p_is_master) this%nobs = sum(counts)
      CALL mpi_bcast(this%nobs, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
      IF (this%nobs == 0) RETURN

      ! Compute receive offsets in patch units and allocate receive buffers.
      IF (p_is_master) THEN
         DO iproc = 1, p_np_glb-1
            displs(iproc) = displs(iproc-1) + counts(iproc-1)
         ENDDO
         allocate(obs_global(4,this%nobs))
      ELSE
         allocate(obs_global(4,1))
      ENDIF
      CALL mpi_gatherv(obs, 4*nlocal, MPI_REAL8, obs_global, 4*counts, 4*displs, MPI_REAL8, p_address_master, p_comm_glb, p_err)

      ! Store the gathered data in object members on master and broadcast to workers
      allocate(this%lat(this%nobs), this%lon(this%nobs), this%y(this%nobs), this%r(this%nobs))
      IF (p_is_master) THEN
         this%lat = obs_global(1,:)
         this%lon = obs_global(2,:)
         this%y   = obs_global(3,:)
         this%r   = obs_global(4,:)
      ENDIF
      CALL mpi_bcast(this%lat, this%nobs, MPI_REAL8, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast(this%lon, this%nobs, MPI_REAL8, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast(this%y,   this%nobs, MPI_REAL8, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast(this%r,   this%nobs, MPI_REAL8, p_address_master, p_comm_glb, p_err)
#endif

   END SUBROUTINE pseudo_obs_read

!-----------------------------------------------------------------------------

   SUBROUTINE pseudo_obs_clear (this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(pseudo_obs_type), intent(inout) :: this

      this%nobs = 0
      IF (allocated(this%lat)) deallocate(this%lat)
      IF (allocated(this%lon)) deallocate(this%lon)
      IF (allocated(this%y  )) deallocate(this%y  )
      IF (allocated(this%r  )) deallocate(this%r  )

   END SUBROUTINE pseudo_obs_clear

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Data_PseudoType
#endif