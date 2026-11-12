#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Data_SiteType
!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Site observation derived type for data assimilation.
!
!   A site_obs_type represents one station-network observation stream.
!   It stores the station table, the site-to-patch spatial mapping, and
!   Per-timestep observation vectors.
!
!   Source-specific modules for SYNOP, ISMN, SNOTEL, etc build on
!   this type to share storage and search logic. Observation variables
!   include screen-level temperature and humidity, soil moisture, soil
!   temperature, and snow depth.
!
! WORKFLOW:
!   1. Init: Read and broadcast the station table, locate each site in the
!      model pixel grid, and retain the worker and first patch found in
!      the containing pixel. Build this site-to-patch mapping once.
!   2. Read: Read and broadcast observation vectors, select records within
!      the current time step and model domain, and discard records with
!      invalid site indices or unmapped sites.
!   3. Clear: Release the current observation vectors while retaining the
!      station table and spatial mapping for subsequent time steps.
!
! FILE CONVENTIONS:
!   Files are located under <DEF_DA_OBS_DIR>/<sensor>/.
!
!   sites.nc defines the station table using 'latitude' and 'longitude'.
!   Its row order provides a fixed index for each site. This table is used
!   to build the mapping from station locations to model patches.
!
!   <var_name>_YYYY_MM_DD_HH.nc stores observations for each time period as
!   vectors of station records. Each record contains 'latitude', 'longitude',
!   'time', 'id', <var_name>, and 'r'. The 1-based row index id locates the
!   station in sites.nc, linking each observation to its mapped model patch.
!
! SUPPORTED:
!   +--------+----------+----------------------+--------------------------+
!   | Source | Variable | Assimilated state    | Example datasets         |
!   +--------+----------+----------------------+--------------------------+
!   | site   | tref     | 2 m air temperature  | SYNOP                    |
!   +--------+----------+----------------------+--------------------------+
!   | site   | qref     | 2 m specific humidity| SYNOP                    |
!   +--------+----------+----------------------+--------------------------+
!   | site   | SM       | Soil moisture        | ISMN                     |
!   +--------+----------+----------------------+--------------------------+
!   | site   | SD       | snow depth           | ISMN                     |
!   +--------+----------+----------------------+--------------------------+
!
! HISTORY:
!   Lu Li, 05/2026: Initial version
!   Lu Li, Zhilong Fan, 11/2026: DA v1.0
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_NetCDFSerial
   USE MOD_Pixelset
   USE MOD_Pixel
   USE MOD_Mesh
   USE MOD_Vars_Global
   USE MOD_TimeManager
   IMPLICIT NONE

!-----------------------------------------------------------------------------
! Site observation derived type
!-----------------------------------------------------------------------------
   type :: site_obs_type

      character(len=16)     :: source_name = ''                  ! "site"
      character(len=32)     :: sensor_name = ''                  ! E.g., "SYNOP", "SNOTEL"
      character(len=16)     :: var_name    = ''                  ! NetCDF variable name, e.g., "tref", "qref"

      integer               :: nsite                             ! Number of sites in sites.nc
      real(r8), allocatable :: site_lat(:)                       ! Site latitude 
      real(r8), allocatable :: site_lon(:)                       ! Site longitude 
      integer,  allocatable :: site_iworker(:)                   ! Global rank of worker owning the site
      integer,  allocatable :: site_iworker_ipatch(:)            ! Local patch index on the owning worker

      integer               :: nobs                              ! Number of selected observations
      integer,  allocatable :: isite(:)                          ! Row index in sites.nc for each observation
      real(r8), allocatable :: lat(:)                            ! Latitude of each observation 
      real(r8), allocatable :: lon(:)                            ! Longitude of each observation 
      real(r8), allocatable :: time(:)                           ! UTC time of each observation [s since begin of current day]
      real(r8), allocatable :: y(:)                              ! Observed value
      real(r8), allocatable :: r(:)                              ! Observation error standard deviation

   CONTAINS

      procedure, PUBLIC     :: init   => site_obs_init           ! Read sites and build site-to-patch map
      procedure, PUBLIC     :: read   => site_obs_read           ! Read observations for this time step
      procedure, PUBLIC     :: clear  => site_obs_clear          ! Clear per-timestep data (keep site table)

   END type site_obs_type

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE site_obs_init (this, pixelset, source_name, sensor_name, var_name)

!-----------------------------------------------------------------------------
      class(site_obs_type), intent(inout) :: this                ! Site observation collection
      type(pixelset_type),  intent(in)    :: pixelset            ! Local pixel sets used for site-to-patch mapping
      character(len=*),     intent(in)    :: source_name         ! Observation source name
      character(len=*),     intent(in)    :: sensor_name         ! Sensor or station network name
      character(len=*),     intent(in)    :: var_name            ! Observation variable name

      character(len=256)                  :: site_file_name      ! Path to sites.nc
      integer                             :: i                   ! Site index
      integer                             :: np                  ! Local patch index
      integer                             :: ie                  ! Local mesh element index
      integer                             :: ipxl                ! Pixel index within the mesh element
      integer                             :: ipxstt              ! First pixel index of the patch in the mesh element
      integer                             :: ipxend              ! Last pixel index of the patch in the mesh element
      integer                             :: ilat                ! Latitude index in the shared pixel grid
      integer                             :: ilon                ! Longitude index in the shared pixel grid
      integer, allocatable                :: site_pixel_ilat(:)  ! Site pixel latitude indices; -1 if not found
      integer, allocatable                :: site_pixel_ilon(:)  ! Site pixel longitude indices; -1 if not found
!-----------------------------------------------------------------------------
      ! Initialize the observation stream.
      this%source_name = source_name
      this%sensor_name = sensor_name
      this%var_name = var_name
      this%nobs = 0

      ! Read and broadcast site coordinates.
      site_file_name = trim(DEF_DA_OBS_DIR)//'/'//trim(this%sensor_name)//'/sites.nc'
      CALL ncio_read_bcast_serial(site_file_name, 'latitude',  this%site_lat)
      CALL ncio_read_bcast_serial(site_file_name, 'longitude', this%site_lon)
      this%nsite = size(this%site_lat)
      IF (this%nsite == 0) RETURN

      ! Locate each site in the shared pixel grid.
      allocate (site_pixel_ilat(this%nsite))
      allocate (site_pixel_ilon(this%nsite))
      site_pixel_ilat(:) = -1
      site_pixel_ilon(:) = -1

      IF (p_is_master) THEN
         DO i = 1, this%nsite
            DO ilat = 1, pixel%nlat
               IF (pixel%lat_s(ilat) <= this%site_lat(i) .and. &
                  pixel%lat_n(ilat) >  this%site_lat(i)) THEN
                  site_pixel_ilat(i) = ilat
               ENDIF
            ENDDO
            DO ilon = 1, pixel%nlon
               IF (pixel%lon_w(ilon) <= this%site_lon(i) .and. &
                  pixel%lon_e(ilon) >  this%site_lon(i)) THEN
                  site_pixel_ilon(i) = ilon
               ENDIF
            ENDDO
         ENDDO
      ENDIF
#ifdef USEMPI
      CALL mpi_bcast (site_pixel_ilat, this%nsite, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (site_pixel_ilon, this%nsite, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif

      ! Find the patch containing each site's pixel.
      allocate (this%site_iworker        (this%nsite))
      allocate (this%site_iworker_ipatch (this%nsite))
      this%site_iworker       (:) = -1
      this%site_iworker_ipatch(:) = -1

      IF (p_is_worker) THEN
         DO i = 1, this%nsite
            IF (site_pixel_ilat(i) == -1 .or. site_pixel_ilon(i) == -1) CYCLE

            patch_loop: DO np = 1, pixelset%nset
               ie = pixelset%ielm(np)
               ipxstt = pixelset%ipxstt(np)
               ipxend = pixelset%ipxend(np)

               DO ipxl = ipxstt, ipxend
                  IF (mesh(ie)%ilat(ipxl) /= site_pixel_ilat(i)) CYCLE
                  IF (mesh(ie)%ilon(ipxl) /= site_pixel_ilon(i)) CYCLE

                  this%site_iworker(i)        = p_iam_glb
                  this%site_iworker_ipatch(i) = np
                  EXIT patch_loop
               ENDDO
            ENDDO patch_loop
         ENDDO
      ENDIF

      ! Share the site-to-patch mapping across all ranks.
#ifdef USEMPI
      CALL mpi_allreduce (MPI_IN_PLACE, this%site_iworker,        this%nsite, MPI_INTEGER, MPI_MAX, p_comm_glb, p_err)
      CALL mpi_allreduce (MPI_IN_PLACE, this%site_iworker_ipatch, this%nsite, MPI_INTEGER, MPI_MAX, p_comm_glb, p_err)
#endif
      deallocate (site_pixel_ilat, site_pixel_ilon)

   END SUBROUTINE site_obs_init

!-----------------------------------------------------------------------------

   SUBROUTINE site_obs_clear (this)

!-----------------------------------------------------------------------------
      class(site_obs_type), intent(inout) :: this
!-----------------------------------------------------------------------------
      IF (allocated(this%isite)) deallocate (this%isite)
      IF (allocated(this%lat  )) deallocate (this%lat  )
      IF (allocated(this%lon  )) deallocate (this%lon  )
      IF (allocated(this%time )) deallocate (this%time )
      IF (allocated(this%y    )) deallocate (this%y    )
      IF (allocated(this%r    )) deallocate (this%r    )
      this%nobs = 0

   END SUBROUTINE site_obs_clear

!-----------------------------------------------------------------------------

   SUBROUTINE site_obs_read (this, idate, deltim)

!-----------------------------------------------------------------------------
      class(site_obs_type), intent(inout) :: this
      integer,              intent(in)    :: idate(3)
      real(r8),             intent(in)    :: deltim

      integer                             :: i                 
      integer                             :: sdate(3)            ! Timestep start time adjusted to the beginning of step
      integer                             :: month, mday, hour   ! Calendar month, day of month, and UTC hour
      character(len=4)                    :: yearstr             ! YYYY string for file name
      character(len=2)                    :: monthstr            ! MM string for file name
      character(len=2)                    :: daystr              ! DD string for file name
      character(len=2)                    :: hourstr             ! HH string for file name
      character(len=256)                  :: obs_file_name       ! Full path of current timestep obs file
      logical                             :: exists_file         ! Whether the expected obs file exists
      real(r8), allocatable               :: lat(:)              ! Raw obs latitudes
      real(r8), allocatable               :: lon(:)              ! Raw obs longitudes
      real(r8), allocatable               :: time(:)             ! Raw obs times [s since current UTC day]
      real(r8), allocatable               :: y(:)                ! Raw obs values
      real(r8), allocatable               :: r(:)                ! Raw obs error standard deviation
      integer,  allocatable               :: isite(:)            ! Raw obs site index in sites.nc
      logical,  allocatable               :: valid(:)            ! Mask for obs inside time window and domain with mapped sites
!-----------------------------------------------------------------------------
      CALL this%clear()

      sdate = idate
      CALL backtime(deltim, sdate)
      CALL adj2begin(sdate)
      CALL julian2monthday(sdate(1), sdate(2), month, mday)
      hour = int(sdate(3)/3600)
      write (yearstr,  '(I4.4)') sdate(1)
      write (monthstr, '(I2.2)') month
      write (daystr,   '(I2.2)') mday
      write (hourstr,  '(I2.2)') hour

      obs_file_name = trim(DEF_DA_OBS_DIR)//'/'//trim(this%sensor_name)//'/'//trim(this%var_name)//'_'// &
         trim(yearstr)//'_'//trim(monthstr)//'_'//trim(daystr)//'_'//trim(hourstr)//'.nc'
      inquire (file=trim(obs_file_name), exist=exists_file)

      IF (exists_file) THEN
         CALL ncio_read_bcast_serial (obs_file_name, 'latitude',          lat  )
         CALL ncio_read_bcast_serial (obs_file_name, 'longitude',         lon  )
         CALL ncio_read_bcast_serial (obs_file_name, 'time',              time )
         CALL ncio_read_bcast_serial (obs_file_name, 'isite',             isite)
         CALL ncio_read_bcast_serial (obs_file_name, trim(this%var_name), y    )
         CALL ncio_read_bcast_serial (obs_file_name, 'r',                 r    )

         allocate (valid(size(time)))
         valid = time >= sdate(3)          .and.  &
            time <= sdate(3) + deltim      .and.  &
            lat  >= DEF_domain%edges       .and.  &
            lat  <= DEF_domain%edgen       .and.  &
            lon  >= DEF_domain%edgew       .and.  &
            lon  <= DEF_domain%edgee       
         DO i = 1, size(valid)
            IF (valid(i)) valid(i) = this%site_iworker(id(i)) /= -1
         ENDDO
         this%nobs = count(valid)

         IF (this%nobs > 0) THEN
            allocate (this%isite (this%nobs))
            allocate (this%lat   (this%nobs))
            allocate (this%lon   (this%nobs))
            allocate (this%time  (this%nobs))
            allocate (this%y     (this%nobs))
            allocate (this%r     (this%nobs))

            this%isite = pack(isite, valid)
            this%lat   = pack(lat,   valid)
            this%lon   = pack(lon,   valid)
            this%time  = pack(time,  valid)
            this%y     = pack(y,     valid)
            this%r     = pack(r,     valid)
         ENDIF

         IF (allocated(valid)) deallocate (valid)
         IF (allocated(lat  )) deallocate (lat  )
         IF (allocated(lon  )) deallocate (lon  )
         IF (allocated(time )) deallocate (time )
         IF (allocated(y    )) deallocate (y    )
         IF (allocated(r    )) deallocate (r    )
         IF (allocated(isite)) deallocate (isite)
      ENDIF

      IF (p_is_master) THEN
         IF (exists_file) THEN
            WRITE (*,'(A,I0,2A)') '[DA] Selected', this%nobs, 'observations from ', trim(obs_file_name)
         ELSE
            WRITE (*,'(2A)') '[DA] Observation file not found: ', trim(obs_file_name)
         ENDIF
      ENDIF

   END SUBROUTINE site_obs_read

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Data_SiteType
#endif