#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Obs_PseudoType
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Handle one patch-based pseudo-observation variable derived from model
!    analyses.
!
!    In operational forecasting systems, screen-level temperature and humidity
!    analyses can provide observational constraints on land-surface variables,
!    supporting land-atmosphere coupled data assimilation.
!
!    Following SYNOP 2 m temperature and humidity assimilation, patch-level
!    analysis values with nonzero analysis increments are treated as
!    pseudo-observations for subsequent land-surface data assimilation.
!
! HISTORY:
!    Lu Li, 10/2026: First operational version
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_LandPatch
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_Vars_Global
   USE MOD_Vars_TimeInvariants
   IMPLICIT NONE

!-----------------------------------------------------------------------------
! Pseudo-observation derived type
!
! Both y and hx are derived from model ensembles and processed together:
! y contains analysis ensemble means, while hx contains forecast member values.
! Unlike GridType, which reads external observations and leaves hx calculation
! to the source-specific module, this type constructs and stores both y and hx.
!-----------------------------------------------------------------------------
   type :: pseudo_obs_type

      integer :: nobs = 0                            ! Number of changed patches across all workers
      real(r8), allocatable :: lat(:)                ! Patch latitude [nobs], degrees
      real(r8), allocatable :: lon(:)                ! Patch longitude [nobs], degrees
      real(r8), allocatable :: y(:)                  ! Analysis ensemble means [nobs]
      real(r8), allocatable :: r(:)                  ! Observation error variances [nobs]
      real(r8), allocatable :: hx(:,:)               ! Forecast ensemble values [nobs, ens]

   CONTAINS

      procedure, PUBLIC :: init  => pseudo_obs_init
      procedure, PUBLIC :: read  => pseudo_obs_read
      procedure, PUBLIC :: clear => pseudo_obs_clear

   END type pseudo_obs_type

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE pseudo_obs_init (this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(pseudo_obs_type), intent(inout) :: this

!-----------------------------------------------------------------------------
      CALL this%clear()

   END SUBROUTINE pseudo_obs_init

!-----------------------------------------------------------------------------

   SUBROUTINE pseudo_obs_read (this, hx_ens_f, hx_ens_a, r)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      class(pseudo_obs_type), intent(inout) :: this
      real(r8), allocatable, intent(in) :: hx_ens_f(:,:) ! forecast [0:ens, local patch]
      real(r8), allocatable, intent(in) :: hx_ens_a(:,:) ! analysis [0:ens, local patch]
      real(r8), intent(in) :: r                         ! pseudo-observation error variance

!------------------------ Local Variables ------------------------------------
      integer :: np, nlocal, iobs
      logical, allocatable :: valid(:)
      real(r8), allocatable :: y_f(:), y_a(:)
      real(r8), allocatable :: obs(:,:), obs_global(:,:)
      real(r8), allocatable :: hx(:,:), hx_global(:,:)
#ifdef USEMPI
      integer :: iproc
      integer, allocatable :: counts(:), displs(:)
#endif

!-----------------------------------------------------------------------------
      CALL this%clear()

      ! Select local patches with nonzero increments in the ensemble mean.
      nlocal = 0
      IF (p_is_worker) THEN
         allocate(valid(numpatch), y_f(numpatch), y_a(numpatch))
         valid = .false.

         DO np = 1, numpatch
            IF (patchtype(np) >= 3) CYCLE

            y_f(np) = sum(hx_ens_f(1:DEF_DA_ENS_NUM,np)) / DEF_DA_ENS_NUM
            y_a(np) = sum(hx_ens_a(1:DEF_DA_ENS_NUM,np)) / DEF_DA_ENS_NUM
            valid(np) = y_a(np) /= y_f(np)
         ENDDO

         nlocal = count(valid)
      ENDIF

      ! Zero-length send buffers are valid when this rank has no selected patches.
      allocate(obs(4,nlocal), hx(DEF_DA_ENS_NUM,nlocal))

      ! Pack obs (lat, lon, y, r) and H(x) in the same selected-patch order.
      IF (nlocal > 0) THEN
         iobs = 0
         DO np = 1, numpatch
            IF (.not. valid(np)) CYCLE
            iobs = iobs + 1
            obs(1,iobs) = patchlatr(np)*180.0_r8/pi
            obs(2,iobs) = (modulo(patchlonr(np)+pi,2.0_r8*pi)-pi)*180.0_r8/pi
            obs(3,iobs) = y_a(np)
            obs(4,iobs) = r
            hx(:,iobs)  = hx_ens_f(1:DEF_DA_ENS_NUM,np)
         ENDDO
      ENDIF

#ifdef USEMPI
      ! All ranks participate; non-worker ranks contribute zero records.
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
         allocate(obs_global(4,this%nobs), hx_global(DEF_DA_ENS_NUM,this%nobs))
      ELSE
         ! MPI_Gatherv ignores receive buffers on non-root ranks.
         allocate(obs_global(4,1), hx_global(DEF_DA_ENS_NUM,1))
      ENDIF
      ! Gather obs and H(x) separately, preserving the same patch ordering.
      ! Scale patch counts and offsets to numbers of real elements.
      CALL mpi_gatherv(obs, 4*nlocal, MPI_REAL8, obs_global, 4*counts, 4*displs, &
         MPI_REAL8, p_address_master, p_comm_glb, p_err)
      CALL mpi_gatherv(hx, DEF_DA_ENS_NUM*nlocal, MPI_REAL8, hx_global, &
         DEF_DA_ENS_NUM*counts, DEF_DA_ENS_NUM*displs, &
         MPI_REAL8, p_address_master, p_comm_glb, p_err)

      ! Store the gathered data in object members on master.
      allocate(this%lat(this%nobs), this%lon(this%nobs))
      allocate(this%y(this%nobs), this%r(this%nobs), this%hx(this%nobs,DEF_DA_ENS_NUM))
      IF (p_is_master) THEN
         this%lat = obs_global(1,:)
         this%lon = obs_global(2,:)
         this%y   = obs_global(3,:)
         this%r   = obs_global(4,:)
         this%hx  = transpose(hx_global) ! [member, observation] -> [observation, member]
      ENDIF
      ! Broadcast the complete observation vectors and H(x) to all ranks.
      CALL mpi_bcast(this%lat, this%nobs, MPI_REAL8, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast(this%lon, this%nobs, MPI_REAL8, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast(this%y,   this%nobs, MPI_REAL8, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast(this%r,   this%nobs, MPI_REAL8, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast(this%hx,  this%nobs*DEF_DA_ENS_NUM, MPI_REAL8, p_address_master, p_comm_glb, p_err)
#endif

   END SUBROUTINE pseudo_obs_read

!-----------------------------------------------------------------------------

   SUBROUTINE pseudo_obs_clear (this)

!-----------------------------------------------------------------------------
      IMPLICIT NONE

      class(pseudo_obs_type), intent(inout) :: this

!-----------------------------------------------------------------------------
      this%nobs = 0
      IF (allocated(this%lat)) deallocate(this%lat)
      IF (allocated(this%lon)) deallocate(this%lon)
      IF (allocated(this%y  )) deallocate(this%y  )
      IF (allocated(this%r  )) deallocate(this%r  )
      IF (allocated(this%hx )) deallocate(this%hx )

   END SUBROUTINE pseudo_obs_clear

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Obs_PseudoType
#endif
