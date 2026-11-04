#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Obs_SiteType
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Site observation derived type for data assimilation.
!
!    A site_obs_type represents one in-situ observation stream, such as
!    SYNOP, ISMN, or other station networks. It owns:
!      1. the site table read from <DEF_DA_OBS_DIR>/<sensor>/sites.nc;
!      2. the site-to-patch mapping, built once from pixel containment and
!         the first patch found in the containing pixel;
!      3. the per-timestep observation vectors read from
!         <DEF_DA_OBS_DIR>/<sensor>/<var_name>_YYYY_MM_DD_HH.nc.
!
!    Source-specific modules compute H(x), extract values at the mapped
!    patches, and share the results across MPI ranks.
!
!    File conventions:
!      sites.nc                        station table with 'latitude' and
!                                       'longitude'
!      <var_name>_YYYY_MM_DD_HH.nc     per-timestep observations with
!                                       'latitude', 'longitude', 'time',
!                                       'id', <var_name>, and 'r'
!
!    The 'id' variable is the 1-based row index in sites.nc. Observations
!    whose id is out of range or whose site cannot be mapped to a local
!    patch are dropped.
!
! SUPPORTED:
! +--------+----------+----------------------+-------------------------------+
! | Source | Variable | Assimilated state    | Example datasets              |
! +--------+----------+----------------------+-------------------------------+
! | site   | tref     | 2 m air temperature  | SYNOP 2 m air temperature     |
! +--------+----------+----------------------+-------------------------------+
! | site   | qref     | 2 m specific humidity| SYNOP 2 m relative humidity   |
! +--------+----------+----------------------+-------------------------------+
! | site   | SM       | Soil moisture        | ISMN soil moisture            |
! +--------+----------+----------------------+-------------------------------+
! | site   | SD       | snow depth           | ISMN soil temperature         |
! +--------+----------+----------------------+-------------------------------+
!
! HISTORY:
!    Lu Li, 05/2026: Initial version
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

        ! initialize once
        character(len=16) :: source_name = ''          ! "site"
        character(len=32) :: sensor_name = ''          ! e.g., "SYNOP", "SNOTEL"
        character(len=16) :: var_name    = ''          ! NetCDF variable name, e.g., "tref", "qref"

        ! site table (initialize once)
        integer :: nsite = 0                           ! number of sites in sites.nc
        real(r8), allocatable :: site_lat(:)            ! site latitude [deg]
        real(r8), allocatable :: site_lon(:)            ! site longitude [deg]
        integer,  allocatable :: site_iworker(:)        ! global rank of worker owning the site, -1 if unmapped
        integer,  allocatable :: site_iworker_ipatch(:) ! local patch index on the owning worker, -1 if unmapped

        ! observations at this timestep
        integer :: nobs = 0                            ! number of selected observations
        integer,  allocatable :: isite(:)              ! row index in sites.nc for each observation
        real(r8), allocatable :: lat(:)                ! latitude of each observation [deg]
        real(r8), allocatable :: lon(:)                ! longitude of each observation [deg]
        real(r8), allocatable :: time(:)               ! UTC time of each observation [s since begin of current day]
        real(r8), allocatable :: y(:)                  ! observed value
        real(r8), allocatable :: r(:)                  ! observation error standard deviation

    CONTAINS

        procedure, PUBLIC :: init   => site_obs_init   ! read sites and build site-to-patch map
        procedure, PUBLIC :: read   => site_obs_read   ! read observations for this time step
        procedure, PUBLIC :: clear  => site_obs_clear  ! clear per-timestep data (keep site table)

    END type site_obs_type

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

    SUBROUTINE site_obs_init (this, pixelset, source_name, sensor_name, var_name)

!-----------------------------------------------------------------------------
        IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
        class(site_obs_type), intent(inout) :: this        ! Site observation collection
        type(pixelset_type),  intent(in)    :: pixelset    ! Local pixel sets used for site-to-patch mapping
        character(len=*),     intent(in)    :: source_name ! Observation source name
        character(len=*),     intent(in)    :: sensor_name ! Sensor or station network name
        character(len=*),     intent(in)    :: var_name    ! Observation variable name

!------------------------ Local Variables ------------------------------------
        character(len=256) :: site_file_name               ! Path to sites.nc
        integer :: i                                       ! Site index
        integer :: np                                      ! Local patch index
        integer :: ie                                      ! Local mesh element index
        integer :: ipxl                                    ! Pixel index within the mesh element
        integer :: ipxstt                                  ! First pixel index of the patch in the mesh element
        integer :: ipxend                                  ! Last pixel index of the patch in the mesh element
        integer :: ilat                                    ! Latitude index in the shared pixel grid
        integer :: ilon                                    ! Longitude index in the shared pixel grid
        integer, allocatable :: site_pixel_ilat(:)         ! Site pixel latitude indices; -1 if not found
        integer, allocatable :: site_pixel_ilon(:)         ! Site pixel longitude indices; -1 if not found

!-----------------------------------------------------------------------------
        ! Reset the site table and patch mapping.
        this%source_name = source_name
        this%sensor_name = sensor_name
        this%var_name = var_name
        this%nsite = 0
        this%nobs = 0
        IF (allocated(this%site_lat           )) deallocate(this%site_lat           )
        IF (allocated(this%site_lon           )) deallocate(this%site_lon           )
        IF (allocated(this%site_iworker       )) deallocate(this%site_iworker       )
        IF (allocated(this%site_iworker_ipatch)) deallocate(this%site_iworker_ipatch)

        ! Read and broadcast site coordinates.
        site_file_name = trim(DEF_DA_OBS_DIR)//'/'//trim(this%sensor_name)//'/sites.nc'
        CALL ncio_read_bcast_serial(site_file_name, 'latitude',  this%site_lat)
        CALL ncio_read_bcast_serial(site_file_name, 'longitude', this%site_lon)
        this%nsite = size(this%site_lat)
        IF (this%nsite == 0) RETURN

        allocate (this%site_iworker        (this%nsite))
        allocate (this%site_iworker_ipatch (this%nsite))
        this%site_iworker        (:) = -1
        this%site_iworker_ipatch (:) = -1

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

        ! Find the first local patch containing each site's pixel.
        IF (p_is_worker .and. pixelset%nset > 0) THEN
            DO i = 1, this%nsite
                IF (site_pixel_ilat(i) < 1 .or. site_pixel_ilon(i) < 1) CYCLE

                patch_loop: DO np = 1, pixelset%nset
                    ie = pixelset%ielm(np)
                    ipxstt = pixelset%ipxstt(np)
                    ipxend = pixelset%ipxend(np)

                    IF (ipxstt == -1) THEN
                        ipxstt = 1
                        ipxend = mesh(ie)%npxl
                    ENDIF

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
        CALL mpi_allreduce (MPI_IN_PLACE, this%site_iworker, this%nsite, MPI_INTEGER, &
            MPI_MAX, p_comm_glb, p_err)
        CALL mpi_allreduce (MPI_IN_PLACE, this%site_iworker_ipatch, this%nsite, MPI_INTEGER, &
            MPI_MAX, p_comm_glb, p_err)
#endif
        deallocate (site_pixel_ilat, site_pixel_ilon)


    END SUBROUTINE site_obs_init

!-----------------------------------------------------------------------------

    SUBROUTINE site_obs_clear (this)

!-----------------------------------------------------------------------------
        IMPLICIT NONE

        class(site_obs_type), intent(inout) :: this

        IF (allocated(this%isite)) deallocate(this%isite)
        IF (allocated(this%lat))   deallocate(this%lat)
        IF (allocated(this%lon))   deallocate(this%lon)
        IF (allocated(this%time))  deallocate(this%time)
        IF (allocated(this%y))     deallocate(this%y)
        IF (allocated(this%r))     deallocate(this%r)
        this%nobs = 0

    END SUBROUTINE site_obs_clear

!-----------------------------------------------------------------------------

    SUBROUTINE site_obs_read (this, idate, deltim)

!-----------------------------------------------------------------------------
        IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
        class(site_obs_type), intent(inout) :: this
        integer,              intent(in)    :: idate(3)
        real(r8),             intent(in)    :: deltim

!------------------------ Local Variables ------------------------------------
        logical :: exists_file                 ! whether the expected obs file exists
        integer :: i                           ! loop index
        integer :: sdate(3)                    ! timestep start time adjusted to the beginning of step
        integer :: month, mday, hour           ! calendar month, day of month, and UTC hour
        character(len=4)   :: yearstr          ! YYYY string for file name
        character(len=2)   :: monthstr         ! MM string for file name
        character(len=2)   :: daystr           ! DD string for file name
        character(len=2)   :: hourstr          ! HH string for file name
        character(len=256) :: obs_file_name    ! full path of current timestep obs file
        real(r8), allocatable :: lat(:)        ! raw obs latitudes read from file [deg]
        real(r8), allocatable :: lon(:)        ! raw obs longitudes read from file [deg]
        real(r8), allocatable :: time(:)       ! raw obs times read from file [s since current UTC day]
        real(r8), allocatable :: y(:)          ! raw obs values for var
        real(r8), allocatable :: r(:)          ! raw obs error standard deviation for var
        integer,  allocatable :: id(:)         ! raw site row indices in sites.nc
        logical,  allocatable :: valid(:)      ! full-length mask for selected observations

!-----------------------------------------------------------------------------

        CALL this%clear()

        IF (this%nsite == 0) RETURN

        sdate = idate
        CALL backtime(deltim, sdate)
        CALL adj2begin(sdate)
        CALL julian2monthday(sdate(1), sdate(2), month, mday)
        hour = int(sdate(3)/3600)
        write (yearstr, '(I4.4)') sdate(1)
        write (monthstr, '(I2.2)') month
        write (daystr, '(I2.2)') mday
        write (hourstr, '(I2.2)') hour

        obs_file_name = trim(DEF_DA_OBS_DIR)//'/'//trim(this%sensor_name)//'/'//trim(this%var_name)//'_'// &
            trim(yearstr)//'_'//trim(monthstr)//'_'//trim(daystr)//'_'//trim(hourstr)//'.nc'
        inquire (file=trim(obs_file_name), exist=exists_file)

        IF (exists_file) THEN
            CALL ncio_read_bcast_serial(obs_file_name, 'latitude', lat)
            CALL ncio_read_bcast_serial(obs_file_name, 'longitude', lon)
            CALL ncio_read_bcast_serial(obs_file_name, 'time', time)
            CALL ncio_read_bcast_serial(obs_file_name, 'id', id)
            CALL ncio_read_bcast_serial(obs_file_name, trim(this%var_name), y)
            CALL ncio_read_bcast_serial(obs_file_name, 'r', r)

            allocate (valid(size(time)))
            valid = time >= sdate(3)          .and. &
                    time <= sdate(3) + deltim .and. &
                    lat  >= DEF_domain%edges  .and. &
                    lat  <= DEF_domain%edgen  .and. &
                    lon  >= DEF_domain%edgew  .and. &
                    lon  <= DEF_domain%edgee  .and. &
                    id   >= 1                 .and. &
                    id   <= this%nsite

            DO i = 1, size(valid)
                IF (valid(i)) THEN
                    valid(i) = this%site_iworker(id(i)) >= 0
                ENDIF
            ENDDO

            this%nobs = count(valid)
            IF (this%nobs > 0) THEN
                allocate (this%isite(this%nobs))
                allocate (this%lat(this%nobs))
                allocate (this%lon(this%nobs))
                allocate (this%time(this%nobs))
                allocate (this%y(this%nobs))
                allocate (this%r(this%nobs))

                this%isite = pack(id  , valid)
                this%lat   = pack(lat , valid)
                this%lon   = pack(lon , valid)
                this%time  = pack(time, valid)
                this%y     = pack(y   , valid)
                this%r     = pack(r   , valid)
            ENDIF

            IF (allocated(valid)) deallocate(valid)
            IF (allocated(lat  )) deallocate(lat  )
            IF (allocated(lon  )) deallocate(lon  )
            IF (allocated(time )) deallocate(time )
            IF (allocated(y    )) deallocate(y    )
            IF (allocated(r    )) deallocate(r    )
            IF (allocated(id   )) deallocate(id   )
        ENDIF

        IF (p_is_master) THEN
            IF (exists_file) THEN
                WRITE (*,'(A,I0,2A)') '[DA] Selected ', this%nobs, &
                    ' observations from ', trim(obs_file_name)
            ELSE
                WRITE (*,'(2A)') '[DA] Observation file not found: ', trim(obs_file_name)
            ENDIF
        ENDIF

    END SUBROUTINE site_obs_read

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Obs_SiteType
#endif