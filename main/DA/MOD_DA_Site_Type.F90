#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Site_Type
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
!         <DEF_DA_OBS_DIR>/<sensor>/<var_name>_YYYY_MM_DD_HH.nc;
!      4. a generic gather that maps a patch-space field, for example
!         ensemble 2 m temperature, to the observation locations.
!
!    The source-specific observation operator, for example the SYNOP
!    observation operator, stays in its own module. It supplies the
!    patch-space H(x) field and calls gather to obtain H(x) at sites.
!
!    File conventions:
!      sites.nc                         station table with 'latitude' and
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
! | situ   | tref     | 2 m air temperature  | SYNOP 2 m air temperature     |
! +--------+----------+----------------------+-------------------------------+
! | situ   | qref     | 2 m specific humidity| SYNOP 2 m relative humidity   |
! +--------+----------+----------------------+-------------------------------+
! | situ   | SM       | Soil moisture        | ISMN soil moisture            |
! +--------+----------+----------------------+-------------------------------+
! | situ   | ST       | Soil temperature     | ISMN soil temperature         |
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
    USE MOD_LandPatch, ONLY: numpatch
    USE MOD_Vars_Global, ONLY: spval
    USE MOD_TimeManager
    IMPLICIT NONE

    type :: site_obs_type

        ! initialize once
        character(len=16) :: source_name = ''        ! "situ"
        character(len=32) :: sensor_name = ''        ! sensor/subdirectory name, e.g., "SYNOP"
        character(len=16) :: var_name    = ''        ! NetCDF variable name, e.g., "tref", "qref", "SM", "ST"

        ! site table (initialize once)
        integer :: nsite = 0                         ! number of sites in sites.nc
        real(r8), allocatable :: site_lat(:)         ! site latitude [deg]
        real(r8), allocatable :: site_lon(:)         ! site longitude [deg]
        integer,  allocatable :: site_worker(:)      ! global rank owning the site, -1 if unmapped
        integer,  allocatable :: site_patch(:)       ! local patch index owning the site, -1 if unmapped

        ! observations at this timestep
        integer :: nobs = 0                          ! number of selected observations
        integer,  allocatable :: isite(:)            ! row index in sites.nc for each observation
        real(r8), allocatable :: lat(:)              ! latitude of each observation [deg]
        real(r8), allocatable :: lon(:)              ! longitude of each observation [deg]
        real(r8), allocatable :: time(:)             ! UTC time of each observation [s since begin of current day]
        real(r8), allocatable :: y(:)                ! observed value
        real(r8), allocatable :: r(:)                ! observation error variance

    CONTAINS

        procedure, PUBLIC :: init   => site_obs_init   ! read sites and build site-to-patch map
        procedure, PUBLIC :: read   => site_obs_read   ! read observations for this time step
        procedure, PUBLIC :: gather => site_obs_gather ! map a patch-space field to site locations
        procedure, PUBLIC :: clear  => site_obs_clear  ! clear per-timestep data (keep site table)

    END type site_obs_type

    ! Per-pixel-latitude lists of pixel longitude indices and the patch that
    ! contains each pixel. Used only while building the site-to-patch map.
    type :: site_pixel_index_type
        integer, allocatable :: ilon(:)
        integer, allocatable :: ipatch(:)
    END type site_pixel_index_type


!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

    SUBROUTINE site_obs_init (this, pixelset, source_name, sensor_name, var_name)

!-----------------------------------------------------------------------------
        IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
        class(site_obs_type), intent(inout) :: this
        type(pixelset_type),  intent(in)    :: pixelset
        character(len=*),     intent(in)    :: source_name
        character(len=*),     intent(in)    :: sensor_name
        character(len=*),     intent(in)    :: var_name

!------------------------ Local Variables ------------------------------------
        character(len=256) :: site_file_name   ! full path of the site table
        logical :: exists_site                 ! whether the site table exists
        integer :: i, j, np, ie, ipxl          ! loop indices
        integer :: ipxstt, ipxend              ! local pixel range of a patch
        integer :: ilat, ilon                  ! pixel indices of a site or patch
        integer :: numpxl_lat, numpxl_lon      ! number of containing pixel rows and columns
        integer, allocatable :: ilat_site(:)   ! pixel latitude index of each site, -1 if invalid
        integer, allocatable :: ilon_site(:)   ! pixel longitude index of each site, -1 if invalid
        integer, allocatable :: counter(:)     ! number of patches attached to each pixel row
        integer, allocatable :: lut_worker(:)  ! global rank owning each site
        integer, allocatable :: lut_patch(:)   ! local patch owning each site
        type(site_pixel_index_type), allocatable :: idx(:)

!-----------------------------------------------------------------------------

        this%source_name = source_name
        this%sensor_name = sensor_name
        this%var_name = var_name

        IF (allocated(this%site_lat))    deallocate(this%site_lat)
        IF (allocated(this%site_lon))    deallocate(this%site_lon)
        IF (allocated(this%site_worker)) deallocate(this%site_worker)
        IF (allocated(this%site_patch))  deallocate(this%site_patch)
        this%nsite = 0

        site_file_name = trim(DEF_DA_OBS_DIR)//'/'//trim(this%sensor_name)//'/sites.nc'
        exists_site = .false.
        IF (p_is_master) THEN
            exists_site = ncio_var_exist(site_file_name, 'latitude', .false.)
        ENDIF
#ifdef USEMPI
        CALL mpi_bcast (exists_site, 1, MPI_LOGICAL, p_address_master, p_comm_glb, p_err)
#endif
        IF (.not. exists_site) THEN
            IF (p_is_master) WRITE (*,'(A)') '[CoLM-DA] site table not found: '//trim(site_file_name)
            RETURN
        ENDIF

        CALL ncio_read_bcast_serial(site_file_name, 'latitude', this%site_lat)
        CALL ncio_read_bcast_serial(site_file_name, 'longitude', this%site_lon)
        this%nsite = size(this%site_lat)
        IF (this%nsite == 0) RETURN

        allocate (this%site_worker(this%nsite))
        allocate (this%site_patch(this%nsite))
        this%site_worker(:) = -1
        this%site_patch(:) = -1

        ! Find the pixel row and column that contain each site. A site is
        ! usable only when it falls inside exactly one pixel row and one
        ! pixel column.
        allocate (ilat_site(this%nsite))
        allocate (ilon_site(this%nsite))
        ilat_site(:) = -1
        ilon_site(:) = -1

        IF (p_is_master) THEN
            DO i = 1, this%nsite
                numpxl_lat = 0
                DO ilat = 1, pixel%nlat
                    IF (pixel%lat_s(ilat) <= this%site_lat(i) .and. &
                        pixel%lat_n(ilat) >  this%site_lat(i)) THEN
                        numpxl_lat = numpxl_lat + 1
                        ilat_site(i) = ilat
                    ENDIF
                ENDDO
                IF (numpxl_lat /= 1) ilat_site(i) = -1

                numpxl_lon = 0
                DO ilon = 1, pixel%nlon
                    IF (pixel%lon_w(ilon) <= this%site_lon(i) .and. &
                        pixel%lon_e(ilon) >  this%site_lon(i)) THEN
                        numpxl_lon = numpxl_lon + 1
                        ilon_site(i) = ilon
                    ENDIF
                ENDDO
                IF (numpxl_lon /= 1) ilon_site(i) = -1
            ENDDO
        ENDIF
#ifdef USEMPI
        CALL mpi_bcast (ilat_site, this%nsite, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
        CALL mpi_bcast (ilon_site, this%nsite, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif

        ! Build, on each worker, a list of (pixel longitude, patch) pairs for
        ! every pixel row. This is the same lookup as the legacy SYNOP
        ! implementation and gives the first patch in the containing pixel.
        IF (p_is_worker .and. pixelset%nset > 0) THEN

            allocate (counter(pixel%nlat))
            counter(:) = 0

            DO np = 1, pixelset%nset
                ie = pixelset%ielm(np)
                ipxstt = pixelset%ipxstt(np)
                ipxend = pixelset%ipxend(np)

                ! -1 means the set element covers the whole mesh element
                IF (ipxstt == -1) THEN
                    ipxstt = 1
                    ipxend = mesh(ie)%npxl
                ENDIF

                DO ipxl = ipxstt, ipxend
                    ilat = mesh(ie)%ilat(ipxl)
                    counter(ilat) = counter(ilat) + 1
                ENDDO
            ENDDO

            allocate (idx(pixel%nlat))
            DO ilat = 1, pixel%nlat
                IF (counter(ilat) > 0) THEN
                    allocate (idx(ilat)%ilon(counter(ilat)))
                    allocate (idx(ilat)%ipatch(counter(ilat)))
                ENDIF
            ENDDO

            counter(:) = 0
            DO np = 1, pixelset%nset
                ie = pixelset%ielm(np)
                ipxstt = pixelset%ipxstt(np)
                ipxend = pixelset%ipxend(np)

                IF (ipxstt == -1) THEN
                    ipxstt = 1
                    ipxend = mesh(ie)%npxl
                ENDIF

                DO ipxl = ipxstt, ipxend
                    ilat = mesh(ie)%ilat(ipxl)
                    ilon = mesh(ie)%ilon(ipxl)
                    counter(ilat) = counter(ilat) + 1
                    idx(ilat)%ilon(counter(ilat)) = ilon
                    idx(ilat)%ipatch(counter(ilat)) = np
                ENDDO
            ENDDO

            ! First patch whose pixel longitude matches the site pixel.
            DO i = 1, this%nsite
                ilat = ilat_site(i)
                ilon = ilon_site(i)
                IF (ilat < 1 .or. ilon < 1) CYCLE
                IF (.not. allocated(idx(ilat)%ilon)) CYCLE

                DO j = 1, size(idx(ilat)%ilon)
                    IF (idx(ilat)%ilon(j) == ilon) THEN
                        this%site_worker(i) = p_iam_glb
                        this%site_patch(i)  = idx(ilat)%ipatch(j)
                        EXIT
                    ENDIF
                ENDDO
            ENDDO

            deallocate (counter)
            DO ilat = 1, pixel%nlat
                IF (allocated(idx(ilat)%ilon))   deallocate(idx(ilat)%ilon)
                IF (allocated(idx(ilat)%ipatch)) deallocate(idx(ilat)%ipatch)
            ENDDO
            deallocate (idx)

        ENDIF

        ! Combine worker-local maps into one global site-to-patch table.
        ! A site belongs to at most one worker, so MPI_MAX selects the owner.
        allocate (lut_worker(this%nsite))
        allocate (lut_patch(this%nsite))
        lut_worker(:) = this%site_worker(:)
        lut_patch(:)  = this%site_patch(:)
#ifdef USEMPI
        CALL mpi_allreduce (MPI_IN_PLACE, lut_worker, this%nsite, MPI_INTEGER, &
            MPI_MAX, p_comm_glb, p_err)
        CALL mpi_allreduce (MPI_IN_PLACE, lut_patch, this%nsite, MPI_INTEGER, &
            MPI_MAX, p_comm_glb, p_err)
#endif
        this%site_worker(:) = lut_worker(:)
        this%site_patch(:)  = lut_patch(:)

        deallocate (lut_worker, lut_patch, ilat_site, ilon_site)

        IF (p_is_master) THEN
            WRITE (*,'(A,I0,A,I0,A)') '[CoLM-DA] sites: ', this%nsite, &
                ', mapped: ', count(this%site_worker >= 0), &
                ' ('//trim(site_file_name)//')'
        ENDIF

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
        real(r8), allocatable :: r(:)          ! raw obs error variance for var
        integer,  allocatable :: id(:)         ! raw site row indices in sites.nc
        logical,  allocatable :: valid(:)      ! full-length mask for selected observations
        character(len=512) :: logmsg           ! log message buffer

!-----------------------------------------------------------------------------

        CALL this%clear()

        IF (this%nsite == 0) RETURN

        ! UTC begin time of the current DA window
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

            ! drop sites that are not mapped to any local patch
            DO i = 1, size(valid)
                IF (valid(i)) THEN
                    valid(i) = this%site_worker(id(i)) >= 0
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

        write (logmsg,'(A,A,A,L1,A,I0)') 'site obs file=', trim(obs_file_name), &
            ' exists=', exists_file, ' selected_nobs=', this%nobs
        IF (p_is_master) WRITE (*,'(A)') '[CoLM-DA] '//trim(logmsg)

    END SUBROUTINE site_obs_read

!-----------------------------------------------------------------------------

    SUBROUTINE site_obs_gather (this, hx_pset, hx)

!-----------------------------------------------------------------------------
! Map a patch-space field to the observation locations.
!
! This is a collective operation: every MPI rank must call it. It returns the
! same global (nobs,0:ens) array on all ranks. Observations whose site is not
! owned by any worker, or whose patch falls outside the local patch range,
! are set to spval.
!-----------------------------------------------------------------------------
        IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
        class(site_obs_type),  intent(inout) :: this
        real(r8), allocatable, intent(in)    :: hx_pset(:,:)  ! H(x) on patches [0:ens, patch]
        real(r8), allocatable, intent(out)   :: hx(:,:)       ! H(x) at obs [nobs, 0:ens]

!------------------------ Local Variables ------------------------------------
        integer :: iobs, isite, np
        integer :: nens
        integer, allocatable :: owner(:)       ! number of ranks that filled each obs

!-----------------------------------------------------------------------------

        nens = DEF_DA_ENS_NUM
        allocate (hx(this%nobs, 0:nens))
        hx(:,:) = spval
        IF (this%nobs == 0) RETURN

        allocate (owner(this%nobs))
        owner(:) = 0
        hx(:,:) = 0.0_r8

        IF (p_is_worker .and. numpatch > 0 .and. allocated(hx_pset)) THEN
            IF (size(hx_pset, 2) >= numpatch) THEN
                DO iobs = 1, this%nobs
                    isite = this%isite(iobs)
                    IF (this%site_worker(isite) == p_iam_glb) THEN
                        np = this%site_patch(isite)
                        IF (np >= 1 .and. np <= numpatch) THEN
                            hx(iobs, 0:nens) = hx_pset(0:nens, np)
                            owner(iobs) = 1
                        ENDIF
                    ENDIF
                ENDDO
            ENDIF
        ENDIF

#ifdef USEMPI
        CALL mpi_allreduce (MPI_IN_PLACE, hx, this%nobs*(nens+1), MPI_REAL8, &
            MPI_SUM, p_comm_glb, p_err)
        CALL mpi_allreduce (MPI_IN_PLACE, owner, this%nobs, MPI_INTEGER, &
            MPI_SUM, p_comm_glb, p_err)
#endif

        DO iobs = 1, this%nobs
            IF (owner(iobs) /= 1) hx(iobs, 0:nens) = spval
        ENDDO

        deallocate (owner)

    END SUBROUTINE site_obs_gather

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Site_Type
#endif
