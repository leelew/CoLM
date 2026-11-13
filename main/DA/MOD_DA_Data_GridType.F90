#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Data_GridType
!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Grid observation derived type for data assimilation.
!
!   A grid_obs_type represents one gridded observation stream: a satellite or
!   gridded product that comes with its own spatial grid. It stores the
!   observation grid, the grid-to-patch spatial mapping, and
!   timestep-dependent observation vectors.
!
!   Source-specific modules for passive microwave (PM), GNSS reflectometry
!   (GNSS-R), and moderate-resolution optical/thermal remote sensing (MR)
!   build on this type to share storage and search logic. Observation
!   variables include brightness temperature (TB), retrieved soil moisture
!   (SM), delay-Doppler maps (DDM), land surface temperature (LST), and
!   fractional snow-covered area (fSCA).
!
! WORKFLOW:
!   1. Init: Read the global reference grid and build the area-weighted
!      mapping between observation cells and model patches once.
!   2. Read: Read and broadcast the orbit/swath observation vectors, select
!      records within the current time step and model domain, and retain their
!      values, metadata, and reference-grid indices on each process.
!   3. Clear: Release the current observation vectors while retaining the
!      reference grid and spatial mapping for subsequent time steps.
!
! FILE CONVENTIONS:
!   Files are located under <DEF_DA_OBS_DIR>/<sensor>/.
!
!   grid.nc defines a complete global reference grid using cell-center
!   coordinates 'latitude' and 'longitude'. This static grid is used to build
!   the spatial mapping between model patches and observation cells. Its cell
!   size should approximate the satellite observation footprint.
!
!   <var_name>_YYYY_MM_DD_HH.nc stores observations for each time period as
!   vectors following the satellite orbit or swath, rather than as a full
!   global field. Each record contains 'latitude', 'longitude', 'time',
!   <var_name>, 'r', 'theta', 'ii', and 'jj'. The 1-based latitude index ii
!   and longitude index jj locate the record in grid.nc, linking each
!   observation to the corresponding cell and its mapped model patches.
!
! SUPPORTED:
!   +--------+----------+----------------------+----------------------------+
!   | Source | Variable | Assimilated state    | Example datasets           |
!   +--------+----------+----------------------+----------------------------+
!   | PM     | TB_H     | Soil moisture, soil  | SMAP L1C, FY3D MWRI,       |
!   |        | TB_V     | temperature, snow    | AMSR-E.                    |
!   |        |          | variables            |                            |
!   +--------+----------+----------------------+----------------------------+
!   | PM     | SM       | Soil moisture        | SMAP L2 SM, SMOS L2 SM,    |
!   |        |          |                      | AMSR2 SM                   |
!   +--------+----------+----------------------+----------------------------+
!   | GNSS   | DDM      | Soil moisture        | FY3E, CYGNSS               |
!   +--------+----------+----------------------+----------------------------+
!   | MR     | LST      | Land surface         | MODIS LST, VIIRS LST, FY3F |
!   |        |          | temperature, soil    |                            |
!   |        |          | temperature          |                            |
!   +--------+----------+----------------------+----------------------------+
!   | MR     | fSCA     | Snow variables       | MODIS, VIIRS, FY3F         |
!   +--------+----------+----------------------+----------------------------+
!
! HISTORY:
!   Lu Li, 05/2026: Initial version
!   Lu Li, Zhilong Fan, 11/2026: DA v1.0
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Grid
   USE MOD_Namelist
   USE MOD_NetCDFSerial
   USE MOD_Pixelset
   USE MOD_SpatialMapping
   USE MOD_TimeManager
   IMPLICIT NONE

!-----------------------------------------------------------------------------
! Grid observation derived type
!-----------------------------------------------------------------------------
   type :: grid_obs_type

      character(len=16)                     :: source_name = ''         ! E.g., "PM"
      character(len=32)                     :: sensor_name = ''         ! E.g., "SMAP_L1C_D"
      character(len=16)                     :: var_name    = ''         ! E.g., "TB_H"

      type(grid_type)                       :: grid                     ! Grid info for obs
      type(spatial_mapping_type)            :: mg2p                     ! Mapping from obs grid to patch

      integer                               :: nobs                     ! Total avaliable obs
      real(r8), allocatable                 :: lat(:)                   ! Latitude of each obs point
      real(r8), allocatable                 :: lon(:)                   ! Longitude of each obs point
      real(r8), allocatable                 :: time(:)                  ! UTC time of each obs point (seconds since begin of current day)
      real(r8), allocatable                 :: y(:)                     ! Observed value
      real(r8), allocatable                 :: r(:)                     ! Observation error standard deviation
      real(r8), allocatable                 :: theta(:)                 ! Incidence angle for each obs point [rad]
      integer,  allocatable                 :: ii(:)                    ! Latitude index of each observation in the global reference grid
      integer,  allocatable                 :: jj(:)                    ! Longitude index of each observation in the global reference grid

   CONTAINS

      procedure, PUBLIC                     :: init  => grid_obs_init   ! Initialize the obs type (grid info, mapping, etc.)
      procedure, PUBLIC                     :: read  => grid_obs_read   ! Read obs data for this time step (from file)
      procedure, PUBLIC                     :: clear => grid_obs_clear  ! Clear obs data for this time step (keep init info)

   END type grid_obs_type

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE grid_obs_init(this, pixelset, source_name, sensor_name, var_name)

!-----------------------------------------------------------------------------
      class(grid_obs_type),  intent(inout)  :: this
      type(pixelset_type),   intent(in)     :: pixelset
      character(len=*),      intent(in)     :: source_name
      character(len=*),      intent(in)     :: sensor_name
      character(len=*),      intent(in)     :: var_name

      character(len=256)                    :: grid_file_name
!-----------------------------------------------------------------------------
      this%source_name = source_name
      this%sensor_name = sensor_name
      this%var_name = var_name
      this%nobs = 0

      grid_file_name = trim(DEF_DA_OBS_DIR)//'/'//trim(this%sensor_name)//'/grid.nc'
      CALL this%grid%define_from_file(grid_file_name, 'latitude', 'longitude')
      CALL this%mg2p%build_arealweighted(this%grid, pixelset)

   END SUBROUTINE grid_obs_init

!-----------------------------------------------------------------------------

   SUBROUTINE grid_obs_clear(this)

!-----------------------------------------------------------------------------
      class(grid_obs_type),   intent(inout) :: this
!-----------------------------------------------------------------------------
      IF (allocated(this%lat  ))  deallocate(this%lat  )
      IF (allocated(this%lon  ))  deallocate(this%lon  )
      IF (allocated(this%time ))  deallocate(this%time )
      IF (allocated(this%y    ))  deallocate(this%y    )
      IF (allocated(this%r    ))  deallocate(this%r    )
      IF (allocated(this%theta))  deallocate(this%theta)
      IF (allocated(this%ii   ))  deallocate(this%ii   )
      IF (allocated(this%jj   ))  deallocate(this%jj   )
      this%nobs = 0

   END SUBROUTINE grid_obs_clear

!-----------------------------------------------------------------------------

   SUBROUTINE grid_obs_read(this, idate, deltim)

!-----------------------------------------------------------------------------
      class(grid_obs_type), intent(inout)   :: this
      integer,              intent(in)      :: idate(3)
      real(r8),             intent(in)      :: deltim

      integer                               :: sdate(3)                 ! Timestep start time adjusted to the beginning of step
      integer                               :: month, mday, hour        ! Calendar month, day of month, and UTC hour
      character(len=4)                      :: yearstr                  ! YYYY string for file name
      character(len=2)                      :: monthstr                 ! MM string for file name
      character(len=2)                      :: daystr                   ! DD string for file name
      character(len=2)                      :: hourstr                  ! HH string for file name
      character(len=256)                    :: obs_file_name            ! Full path of current timestep obs file
      logical                               :: exists_file              ! Whether the expected obs file exists
      real(r8), allocatable                 :: lat(:)                   ! Raw obs latitudes
      real(r8), allocatable                 :: lon(:)                   ! Raw obs longitudes 
      real(r8), allocatable                 :: time(:)                  ! Raw obs times [s since current UTC day]
      real(r8), allocatable                 :: y(:)                     ! Raw obs values 
      real(r8), allocatable                 :: r(:)                     ! Raw obs error standard deviation 
      real(r8), allocatable                 :: theta(:)                 ! Raw incidence angles [rad]
      integer,  allocatable                 :: ii(:)                    ! Raw obs ii index 
      integer,  allocatable                 :: jj(:)                    ! Raw obs jj index 
      logical,  allocatable                 :: valid(:)                 ! Full-length mask for obs inside time window and domain
!-----------------------------------------------------------------------------
      CALL this%clear()

      ! UTC begin time of the current DA window
      sdate = idate
      CALL backtime(deltim, sdate)
      CALL adj2begin(sdate)
      CALL julian2monthday(sdate(1), sdate(2), month, mday)
      hour = int(sdate(3)/3600)
      write (yearstr,  '(I4.4)') sdate(1)
      write (monthstr, '(I2.2)') month
      write (daystr,   '(I2.2)') mday
      write (hourstr,  '(I2.2)') hour

      ! Whether exists obs for this timestep at domain
      obs_file_name = trim(DEF_DA_OBS_DIR)//'/'//trim(this%sensor_name)//'/'//trim(this%var_name)//'_'// &
         trim(yearstr)//'_'//trim(monthstr)//'_'//trim(daystr)//'_'//trim(hourstr)//'.nc'
      inquire (file=trim(obs_file_name), exist=exists_file)

      IF (exists_file) THEN
         CALL ncio_read_bcast_serial(obs_file_name, 'latitude',          lat  )
         CALL ncio_read_bcast_serial(obs_file_name, 'longitude',         lon  )
         CALL ncio_read_bcast_serial(obs_file_name, 'time',              time )
         CALL ncio_read_bcast_serial(obs_file_name, trim(this%var_name), y    )
         CALL ncio_read_bcast_serial(obs_file_name, 'r',                 r    )
         CALL ncio_read_bcast_serial(obs_file_name, 'theta',             theta)
         CALL ncio_read_bcast_serial(obs_file_name, 'ii',                ii   )
         CALL ncio_read_bcast_serial(obs_file_name, 'jj',                jj   )

         allocate(valid(size(time)))
         valid = time >= sdate(3)          .and.  &
            time <= sdate(3) + deltim      .and.  &
            lat  >= DEF_domain%edges       .and.  &
            lat  <= DEF_domain%edgen       .and.  &
            lon  >= DEF_domain%edgew       .and.  &
            lon  <= DEF_domain%edgee
         this%nobs = count(valid)

         IF (this%nobs > 0) THEN
            allocate(this%lat   (this%nobs))
            allocate(this%lon   (this%nobs))
            allocate(this%time  (this%nobs))
            allocate(this%y     (this%nobs))
            allocate(this%r     (this%nobs))
            allocate(this%theta (this%nobs))
            allocate(this%ii    (this%nobs))
            allocate(this%jj    (this%nobs))

            this%lat    = pack(lat,   valid)
            this%lon    = pack(lon,   valid)
            this%time   = pack(time,  valid)
            this%y      = pack(y,     valid)
            this%r      = pack(r,     valid)
            this%theta  = pack(theta, valid)
            this%ii     = pack(ii,    valid)
            this%jj     = pack(jj,    valid)
         ENDIF

         IF (allocated(valid)) deallocate(valid)
         IF (allocated(lat  )) deallocate(lat  )
         IF (allocated(lon  )) deallocate(lon  )
         IF (allocated(time )) deallocate(time )
         IF (allocated(y    )) deallocate(y    )
         IF (allocated(r    )) deallocate(r    )
         IF (allocated(theta)) deallocate(theta)
         IF (allocated(ii   )) deallocate(ii   )
         IF (allocated(jj   )) deallocate(jj   )
      ENDIF

      IF (p_is_master) THEN
         IF (exists_file) THEN
            WRITE (*,'(A,I0,2A)') '[DA] Selected ', this%nobs, 'observations from ', trim(obs_file_name)
         ELSE
            WRITE (*,'(2A)') '[DA] Observation file not found: ', trim(obs_file_name)
         ENDIF
      ENDIF

   END SUBROUTINE grid_obs_read

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Data_GridType
#endif
