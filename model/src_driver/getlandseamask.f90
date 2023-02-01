MODULE getlandseamask
  !
  !================================================================================================================================
  !! MODULE   : getlandseamask
  !!
  !> BRIEF      This module contains a few methods in order to build a land surface mask which can serve
  !!            to test the routing and its ability to build a river network.
  !!
  !! DESCRIPTION  : None
  !!
  !! RECENT CHANGE(S): None
  !!
  !! MAIN OUTPUT VARIABLE(S):
  !!
  !! REFERENCES   : None
  !!
  !! FLOWCHART    : This module has 2 possible modes for building the land/sea mask :
  !!                TOPOFILE : This configuration option contains a file name of the orography file
  !!                           which will be used in order to build the land/sea mask.
  !!                IIM_LON, JJM_lon : the number of points the grid should have in longitude and
  !!                                   latitude. This is only used when we build the mask based on orography.
  !!                                   Should these parameters not be provided, then the gird of the orography
  !!                                   file is used.
  !!                CONTFRACFILE :: File name in which we will find variable "Contfrac". This land/sea mask will be 
  !!                                used with the associated lat/lon grid.
  !!                WEST_EAST,  SOUTH_NORTH : The zoom to be performed on the grid so as to retain only this window
  !!                                          of the land-sea mask.
  !! \n
  !_================================================================================================================================

  USE ioipsl_para
  USE mod_orchidee_para
  USE control
  USE constantes_soil_var
  USE constantes_var
  USE constantes

  USE netcdf

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: getlandseamask_init, getlandseamask_read

  INTEGER(i_std), SAVE                           :: iim, jjm, nbland
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:) :: lon, lat
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:) :: contfrac, orog

CONTAINS
  !================================================================================================================================
  !! SUBROUTINE   : getlandseamask_init
  !!
  !> BRIEF         Builds the land/sea mask and returns its size to the calling program.
  !!
  !! DESCRIPTION  : None
  !!
  !! RECENT CHANGE(S): None
  !!
  !! MAIN OUTPUT VARIABLE(S):
  !!
  !! REFERENCES   : None
  !!
  !! FLOWCHART    : None
  !! \n
  !_================================================================================================================================
  SUBROUTINE getlandseamask_init(iim_out, jjm_out, nbland_out)
    !
    ! ARGUMENTS
    !
    INTEGER(i_std), INTENT(out) :: iim_out, jjm_out, nbland_out
    !
    ! LOCAL
    !
    CHARACTER(LEN=120)                                  :: filename                !! File name for topography (var: ROSE in "etopo20.nc")
    CHARACTER(LEN=120)                                  :: contfracfile            !! filename for contfrac
    INTEGER(i_std)                                      :: iim_full, jjm_full      !! Lon/Lat dimensions for ROSE in etopo20.nc 
    REAL(r_std), ALLOCATABLE, DIMENSION (:,:)           :: lon_full, lat_full      !! Lon/Lat array for ROSE
    REAL(r_std), ALLOCATABLE, DIMENSION (:,:)           :: etopo, contfrac_full    !! Array for reading ROSE in etopo20.nc
    INTEGER(i_std), DIMENSION(1)                        :: ibeg, iend
    INTEGER(i_std), DIMENSION(1)                        :: jbeg, jend
    INTEGER(i_std)                                      :: fid, i, j
    REAL(r_std)                                         :: dx, dy
    !
    REAL(r_std), DIMENSION(2)                           :: zoom_lon, zoom_lat               !! For crop region when reading
    !
    LOGICAL                                             :: buildcontfrac=.FALSE.
    !
    !  Read the full topography which we will use to get the land-see mask for all resolutions.
    !
    !Config Key   = TOPOFILE
    !Config Desc  = File containing high resolution topographic informations
    !Config If    = 
    !Config Def   = etopo20.nc
    !Config Help  = This needs to be a netCDF file.
    !Config Units = [-] 
    !- 
    filename = 'NONE'
    CALL getin('TOPOFILE', filename)
    WRITE(*,*) "TOPOFILE = ", filename
    !
    !Config Key   = CONTFRACFILE
    !Config Desc  = File containing a contfrac variable as defined in ORCHIDEE.
    !Config If    = 
    !Config Def   = NONE
    !Config Help  = This needs to be a netCDF file and can be any history file ORCHIDEE
    !               has produced ina previous simulation.
    !Config Units = [-] 
    !- 
    contfracfile = 'NONE'
    CALL getin('CONTFRACFILE', contfracfile)
    WRITE(*,*) "CONTFRACFILE = ", contfracfile
    !
    ! Initialize Orchidee parameter and calendar
    !
    IF ( INDEX(filename,"NONE") <= 0 ) THEN
       IF ( is_root_prc) THEN
          CALL topo_getsize(filename, iim_full, jjm_full)
          !
          ALLOCATE(lon_full(iim_full,jjm_full),lat_full(iim_full,jjm_full))
          ALLOCATE(etopo(iim_full,jjm_full))
          !
          CALL topo_getvar(filename, iim_full, jjm_full, lon_full, lat_full, etopo)
          !
          WRITE(*,*) MINVAL(lon_full), ' < LON_FULL < ', MAXVAL(lon_full)
          WRITE(*,*) MINVAL(lat_full), ' < LAT_FULL < ', MAXVAL(lat_full)
          WRITE(*,*) 'ETOPO :', MINVAL(etopo), MAXVAL(etopo)
          !
          buildcontfrac=.TRUE.
          !
       ENDIF
    ELSE IF ( INDEX( contfracfile,"NONE") <= 0 ) THEN
       IF ( is_root_prc) THEN
          !
          ! Read lon, lat and contfrac from a netCDF file generated by ORCHIDEE
          !
          CALL opencontfrac(contfracfile, iim_full, jjm_full, fid)
          !
          ALLOCATE(lon_full(iim_full,jjm_full), lat_full(iim_full,jjm_full), contfrac_full(iim_full,jjm_full))
          CALL readcontfrac(fid, contfracfile, iim_full, jjm_full, lon_full, lat_full, contfrac_full)
          !
          buildcontfrac=.FALSE.
          !
       ENDIF
    ELSE
       WRITE(*,*) "Neither a orography nore a contfrac file was provided"
       WRITE(*,*) "CONTFRACFILE = ", contfracfile
       write(*,*) "TOPOFILE = ", filename
    ENDIF
    !
    !-
    !- Define the zoom
    !-
    zoom_lon=(/-180,180/)
    zoom_lat=(/-90,90/)
    !
    !Config Key   = LIMIT_WEST
    !Config Desc  = Western limit of region
    !Config If    = [-]
    !Config Def   = -180.
    !Config Help  = Western limit of the region we are 
    !Config         interested in. Between -180 and +180 degrees
    !Config         The model will use the smalest regions from
    !Config         region specified here and the one of the forcing file.
    !Config Units = [Degrees] 
    !- 
    CALL getin_p('LIMIT_WEST',zoom_lon(1))
    !-
    !Config Key   = LIMIT_EAST
    !Config Desc  = Eastern limit of region
    !Config If    = [-]
    !Config Def   = 180.
    !Config Help  = Eastern limit of the region we are
    !Config         interested in. Between -180 and +180 degrees
    !Config         The model will use the smalest regions from
    !Config         region specified here and the one of the forcing file.
    !Config Units = [Degrees] 
    !-
    CALL getin_p('LIMIT_EAST',zoom_lon(2))
    !-
    !Config Key   = LIMIT_NORTH
    !Config Desc  = Northern limit of region
    !Config If    = [-]
    !Config Def   = 90.
    !Config Help  = Northern limit of the region we are
    !Config         interested in. Between +90 and -90 degrees
    !Config         The model will use the smalest regions from
    !Config         region specified here and the one of the forcing file.
    !Config Units = [Degrees]
    !-
    CALL getin_p('LIMIT_NORTH',zoom_lat(2))
    !-
    !Config Key   = LIMIT_SOUTH
    !Config Desc  = Southern limit of region
    !Config If    = [-]
    !Config Def   = -90.
    !Config Help  = Southern limit of the region we are
    !Config         interested in. Between 90 and -90 degrees
    !Config         The model will use the smalest regions from
    !Config         region specified here and the one of the forcing file.
    !Config Units = [Degrees]
    !-
    CALL getin_p('LIMIT_SOUTH',zoom_lat(1))
    IF ( (zoom_lon(1)+180 < EPSILON(zoom_lon(1))) .AND. (zoom_lon(2)-180 < EPSILON(zoom_lon(2))) .AND.&
        &(zoom_lat(1)+90 < EPSILON(zoom_lat(1))) .AND. (zoom_lat(2)-90 < EPSILON(zoom_lat(2))) ) THEN
       !
       !Config Key   = WEST_EAST
       !Config Desc  = Longitude interval to use from the forcing data
       !Config If    = [-]
       !Config Def   = -180, 180
       !Config Help  = This function allows to zoom into the forcing data
       !Config Units = [degrees east] 
       !- 
       CALL getin('WEST_EAST', zoom_lon)
       !
       !Config Key   = SOUTH_NORTH
       !Config Desc  = Latitude interval to use from the forcing data
       !Config If    = [-]
       !Config Def   = -90, 90
       !Config Help  = This function allows to zoom into the forcing data
       !Config Units = [degrees north] 
       !- 
       CALL getin('SOUTH_NORTH', zoom_lat)
    ENDIF
    !
    ! Get the finer coarser grid
    !
    !
    ! We see if the user has specified the resolution, in number of
    ! points he wishes to use.
    !
    IF ( buildcontfrac ) THEN
       IF ( is_root_prc) THEN
          iim=-1
          jjm=-1
          !Config Key   = IIM_LON
          !Config Desc  = Number of points in longitude
          !Config If    = [-]
          !Config Def   = -1
          !Config Help  = This allows to select the resolution at which the land/sea mask should be built.
          !Config Units = - 
          CALL getin('IIM_LON', iim)
          !Config Key   = IIM_LAT
          !Config Desc  = Number of points in latitude
          !Config If    = [-]
          !Config Def   = -1
          !Config Help  = This allows to select the resolution at which the land/sea mask should be built.
          !Config Units = - 
          CALL getin('JJM_LAT', jjm)
          !
          IF (iim > 0 .AND. jjm > 0 ) THEN
             !
             ALLOCATE(lon(iim,jjm), lat(iim,jjm), orog(iim,jjm))
             ! Trung: Calculate resolution
             dx = (MAXVAL(zoom_lon)-MINVAL(zoom_lon))/REAL(iim, r_std)
             dy = (MAXVAL(zoom_lat)-MINVAL(zoom_lat))/REAL(jjm, r_std)
             WRITE (*,*) "Resolution in longitude [degrees] dx = ",dx
             WRITE (*,*) "Resolution in laltiude [degrees] dy = ",dy
             !
             ! Generate lon and lat
             !
             DO i=1,iim
                lon(i,:) = MINVAL(zoom_lon) + (i-1)*dx + dx/2.0 
             ENDDO
             DO j=1,jjm
                lat(:,j) = MINVAL(zoom_lat) + (j-1)*dy + dy/2.0
             ENDDO
             !
             ! Interpolate to new orography
             !
             CALL interpol(iim_full, jjm_full, lon_full, lat_full, etopo, iim, jjm, lon, lat, orog)
             !
             ! Write orog to netcdf for checking (deleted)
             !
          ELSE
             ! Just transfer the data over the zoomed area
             ibeg=MINLOC(ABS(lon_full(:,1)-MINVAL(zoom_lon)))
             iend=MINLOC(ABS(lon_full(:,1)-MAXVAL(zoom_lon)))
             jbeg=MINLOC(ABS(lat_full(1,:)-MINVAL(zoom_lat)))
             jend=MINLOC(ABS(lat_full(1,:)-MAXVAL(zoom_lat)))
             iim = (iend(1)-ibeg(1))+1
             jjm = (jend(1)-jbeg(1))+1
             !
             !
             ALLOCATE(lon(iim,jjm), lat(iim,jjm), orog(iim,jjm), contfrac(iim,jjm))
             !
             DO j=1,jjm
                lon(:,j) = lon_full(ibeg(1):iend(1),j)
             ENDDO
             DO i=1,iim
                lat(i,:) = lat_full(i,jbeg(1):jend(1))
             ENDDO
             ! Some treatment for the orogography
             orog(:,:) = 0.0
             DO i=ibeg(1),iend(1)
                DO j=jbeg(1),jend(1)
                   orog(i-ibeg(1)+1,j-jbeg(1)+1) = MAX(etopo(i,j), 0.0)
                ENDDO
             ENDDO
             !
             ! Generate the contfrac field
             !
             contfrac(:,:)= 0.0
             DO i=1,iim
                DO j=1,jjm
                   IF ( orog(i,j) > 0 ) THEN
                      contfrac(i,j) = 1.0
                   ENDIF
                ENDDO
             ENDDO
          ENDIF
          !
       ENDIF
       !
    ELSE
       !
       ! Just zoom into the contfrac if needed and build a dummy orog
       !
       IF ( is_root_prc) THEN
          IF ( MINVAL(zoom_lon) > -180 .OR. MAXVAL(zoom_lon) < 180 .OR.&
               & MINVAL(zoom_lat) > -90 .OR. MAXVAL(zoom_lat) < 90 ) THEN
             ! Just transfer the data over the zoomed area
             ibeg=MINLOC(ABS(lon_full(:,1)-MINVAL(zoom_lon)))
             iend=MINLOC(ABS(lon_full(:,1)-MAXVAL(zoom_lon)))
             jbeg=MINLOC(ABS(lat_full(1,:)-MINVAL(zoom_lat)))
             jend=MINLOC(ABS(lat_full(1,:)-MAXVAL(zoom_lat)))
             iim = (iend(1)-ibeg(1))+1
             jjm = (jend(1)-jbeg(1))+1
             !
             !
             ALLOCATE(lon(iim,jjm), lat(iim,jjm), orog(iim,jjm), contfrac(iim,jjm))
             !
             DO j=1,jjm
                lon(:,j) = lon_full(ibeg(1):iend(1),j)
             ENDDO
             DO i=1,iim
                lat(i,:) = lat_full(i,jbeg(1):jend(1))
             ENDDO
             ! Some treatment for the orogography
             orog(:,:) = 0.0
             contfrac(:,:) = 0.0
             DO i=ibeg(1),iend(1)
                DO j=jbeg(1),jend(1)
                   orog(i-ibeg(1)+1,j-jbeg(1)+1) = 1.0
                   contfrac(i-ibeg(1)+1,j-jbeg(1)+1) = contfrac_full(i,j)
                ENDDO
             ENDDO
             !
          ELSE
             iim = iim_full
             jjm = jjm_full
             ALLOCATE(lon(iim,jjm), lat(iim,jjm), orog(iim,jjm), contfrac(iim,jjm))
             lon(:,:) = lon_full(:,:)
             lat(:,:) = lat_full(:,:)
             contfrac(:,:) = contfrac_full(:,:)
             orog(:,:) = contfrac_full(:,:)+100
          ENDIF
       ENDIF
       !
    ENDIF
    !
    IF ( is_root_prc) THEN
       WRITE(*,*) "Dimensions : ", iim, jjm
       WRITE(*,*) MINVAL(lon), ' < LON < ', MAXVAL(lon)
       WRITE(*,*) MINVAL(lat), ' < LAT < ', MAXVAL(lat)  
       WRITE(*,*) MINVAL(orog), " < OROG  < ", MAXVAL(orog)
       WRITE(*,*) MINVAL(contfrac), " < contfrac  < ", MAXVAL(contfrac)
    ENDIF
    !
    ! Free some memory
    !
    IF (ALLOCATED(lon_full)) DEALLOCATE(lon_full)
    IF (ALLOCATED(lat_full)) DEALLOCATE(lat_full)
    IF (ALLOCATED(etopo)) DEALLOCATE(etopo)
    IF (ALLOCATED(contfrac_full)) DEALLOCATE(contfrac_full)
    !
    ! Distribute the information to all processors
    !
    CALL bcast(iim)
    CALL bcast(jjm)
    IF ( .NOT. ALLOCATED(lon)) ALLOCATE(lon(iim,jjm))
    IF ( .NOT. ALLOCATED(lat)) ALLOCATE(lat(iim,jjm))
    IF ( .NOT. ALLOCATED(orog)) ALLOCATE(orog(iim,jjm))
    IF ( .NOT. ALLOCATED(contfrac)) ALLOCATE(contfrac(iim,jjm))
    CALL bcast(lon)
    CALL bcast(lat)
    CALL bcast(orog)
    CALL bcast(contfrac)
    !
    ! Count the number of land points
    !
    nbland = COUNT(contfrac > 0.0)
    !
    iim_out = iim
    jjm_out = jjm
    nbland_out = nbland
    !
  END SUBROUTINE getlandseamask_init
  !!
  !================================================================================================================================
  !! SUBROUTINE   : getlandseamask_read
  !!
  !> BRIEF         Returns the land/sea mask
  !!
  !! DESCRIPTION  : None
  !!
  !! RECENT CHANGE(S): None
  !!
  !! MAIN OUTPUT VARIABLE(S):
  !!
  !! REFERENCES   : None
  !!
  !! FLOWCHART    : None
  !! \n
  !_================================================================================================================================
  !!
  SUBROUTINE getlandseamask_read(lon_out, lat_out, mask_out, orog_out)
    !
    ! ARGUMENTS
    !
    REAL(r_std), DIMENSION(:,:), INTENT(out) :: lon_out, lat_out
    REAL(r_std), DIMENSION(:,:), INTENT(out) :: mask_out, orog_out
    !
    ! LOCAL
    !
    lon_out(:,:) = lon(:,:)
    lat_out(:,:) = lat(:,:)
    mask_out(:,:) = contfrac(:,:)
    orog_out(:,:) = orog(:,:)
    !
  END SUBROUTINE getlandseamask_read
  !!
  !!
  !
  !================================================================================================================================
  !! SUBROUTINE   : opencontfrac
  !!
  !> BRIEF         This subroutine gets the dimensions of the fields in the filename. Lat and lon in particular.
  !!
  !! DESCRIPTION  : This routines caters for ORCHIDEE as well as the WRF geo files.
  !!
  !! RECENT CHANGE(S): None
  !!
  !! MAIN OUTPUT VARIABLE(S):
  !!
  !! REFERENCES   : None
  !!
  !! FLOWCHART    : None
  !! \n
  !_================================================================================================================================
  SUBROUTINE opencontfrac(filename, iim_full, jjm_full, fid)
    !
    ! ARGUMENTS
    !
    CHARACTER(LEN=*), INTENT(in) :: filename
    INTEGER(i_std), INTENT(out) :: iim_full, jjm_full, fid
    !
    INTEGER(i_std) :: iret, i, lll, ndims, nvars
    CHARACTER(LEN=20) :: dimname
    !
    iret = NF90_OPEN(filename, NF90_NOWRITE, fid)
    IF (iret /= NF90_NOERR) THEN
       WRITE(*,*) "==>",trim(nf90_strerror(iret))
       WRITE(*,*) "Error opening ", filename
       STOP
    ENDIF
    iret = NF90_INQUIRE (fid, nDimensions=ndims, nVariables=nvars)
    DO i=1,ndims
       iret = NF90_INQUIRE_DIMENSION (fid, i, name=dimname, len=lll)
       CALL change_to_lowercase(dimname)
       IF ( INDEX(dimname,"lon") > 0 ) THEN
          iim_full = lll
       ELSE IF (INDEX(dimname,"west_east") > 0 .AND. LEN_TRIM(dimname) == LEN_TRIM("west_east")) THEN
          iim_full = lll
       ELSE IF ( INDEX(dimname,"lat") > 0 ) THEN
          jjm_full = lll
       ELSE IF (INDEX(dimname,"south_north") > 0 .AND. LEN_TRIM(dimname) == LEN_TRIM("south_north")) THEN
          jjm_full = lll
       ENDIF
    ENDDO
    !
    WRITE(*,*) "XXX iimf, jjmf = ", iim_full, jjm_full
    !
  END SUBROUTINE opencontfrac

  !================================================================================================================================
  !! SUBROUTINE   : readcontfrac
  !!
  !> BRIEF         This subroutine gets the topography out of the file
  !!
  !! DESCRIPTION  : This routines caters for ORCHIDEE as well as the WRF geo files.
  !!
  !! RECENT CHANGE(S): None
  !!
  !! MAIN OUTPUT VARIABLE(S):
  !!
  !! REFERENCES   : None
  !!
  !! FLOWCHART    : None
  !! \n
  !_================================================================================================================================
  SUBROUTINE readcontfrac(fid, filename, iimf, jjmf, lon, lat, contf)
    !
    ! ARGUMENTS
    !
    INTEGER(i_std), INTENT(in)               :: fid
    CHARACTER(LEN=*), INTENT(in)             :: filename
    INTEGER(i_std), INTENT(in)               :: iimf, jjmf
    REAL(r_std), DIMENSION(:,:), INTENT(out) :: lon, lat
    REAL(r_std), DIMENSION(:,:), INTENT(out) :: contf
    !
    ! LOCAL
    !
    INTEGER(i_std) :: iret, ndims, nvars, xid, yid, vid, wid
    INTEGER(i_std) :: xndims, yndims, vndims, wndims
    INTEGER(i_std) :: iv, ndimsvar, i, j
    INTEGER(i_std), DIMENSION(4) :: dimids
    CHARACTER(LEN=60) :: varname, units
    CHARACTER(LEN=3)  :: model='ORC'
    INTEGER(i_std), DIMENSION(3) :: start, count
    REAL(r_std), ALLOCATABLE, DIMENSION(:) :: vtmp
    !
    iret = NF90_INQUIRE (fid, nDimensions=ndims, nVariables=nvars)
    !
    xid = -1
    yid = -1
    vid = -1
    wid = -1
    !
    DO iv=1,nvars
       !
       iret = NF90_INQUIRE_VARIABLE(fid, iv, name=varname, ndims=ndimsvar, dimids=dimids)
       iret = NF90_GET_ATT(fid, iv, 'units', units)
       CALL change_to_lowercase(units)
       CALL change_to_lowercase(varname)
       !
       IF ( INDEX(units,"east") > 0 .OR. INDEX(varname, "xlong_m") > 0 ) THEN
          xid = iv
          xndims = ndimsvar
          IF ( INDEX(varname, "xlong_m") > 0 ) THEN
             model='WRF'
          ENDIF
       ELSE IF (INDEX(units,"north") > 0 .OR. INDEX(varname, "xlat_m") > 0 ) THEN
          yid = iv
          yndims = ndimsvar
          IF ( INDEX(varname, "xlat_m") > 0 ) THEN
             model='WRF'
          ENDIF
       ENDIF
       !
       IF ( INDEX(varname,"contfrac") > 0 ) THEN
          vid = iv
          vndims = ndimsvar
       ENDIF
       IF ( INDEX(varname,"landmask") > 0 ) THEN
          wid = iv
          wndims = ndimsvar
       ENDIF
       !
    ENDDO
    !
    ! Get longitude
    !
    IF ( xid > 0 ) THEN
       IF ( xndims == 1 ) THEN
          ALLOCATE(vtmp(iimf))
          iret = NF90_GET_VAR(fid, xid, vtmp)
          DO j=1,jjmf
             lon(:,j) = vtmp(:)
          ENDDO
          DEALLOCATE(vtmp)
       ELSE IF ( xndims == 2 ) THEN
          iret = NF90_GET_VAR(fid, xid, lon)
       ELSE IF ( xndims == 3 ) THEN
          start = (/1,1,1/)
          count = (/iimf,jjmf,1/)
          iret = NF90_GET_VAR(fid, xid, lon, start, count)
       ELSE
          WRITE(*,*) "Unforeseen number of dimensions for longitude ", xndims
          STOP
       ENDIF
    ELSE
       WRITE(*,*) "could not find the longitude in file ", filename
       STOP
    ENDIF
    !
    ! Get latitude
    !
    IF ( yid > 0 ) THEN
       IF ( yndims == 1 ) THEN
          ALLOCATE(vtmp(jjmf))
          iret = NF90_GET_VAR(fid, xid, vtmp)
          DO i=1,iimf
             lat(i,:) = vtmp(:)
          ENDDO
          DEALLOCATE(vtmp)
       ELSE IF ( yndims == 2 ) THEN
          iret = NF90_GET_VAR(fid, yid, lat)
       ELSE IF ( yndims == 3 ) THEN
          start = (/1,1,1/)
          count = (/iimf,jjmf,1/)
          iret = NF90_GET_VAR(fid, yid, lat, start, count)
       ELSE
          WRITE(*,*) "Unforeseen number of dimensions for latitude ", yndims
          STOP
       ENDIF
    ELSE
       WRITE(*,*) "could not find the latitude in file ", filename
       STOP
    ENDIF
    !
    ! Get the variable
    !
    contf(:,:) = 0.0
    IF ( model == "ORC" ) THEN
       IF ( vid > 0 ) THEN
          iret = NF90_GET_VAR(fid, vid, contf)
          IF (iret /= NF90_NOERR) THEN
             WRITE(*,*) "==>",trim(nf90_strerror(iret))
             WRITE(*,*) "Cannot read variable contfrac"
             STOP
          ENDIF
       ELSE
          WRITE(*,*) "could not find the variable contfrac in file ", filename
          STOP
       ENDIF
    ELSE IF ( model == "WRF" ) THEN
       !
       ! This is a time varying variable so start and count needs to be specified
       !
       start = (/1,1,1/)
       count = (/iimf,jjmf,1/)
       IF ( wid > 0 ) THEN
          iret = NF90_GET_VAR(fid, wid, contf, start, count)
          IF (iret /= NF90_NOERR) THEN
             WRITE(*,*) "==>",trim(nf90_strerror(iret))
             WRITE(*,*) "Cannot read variable landmask"
             STOP
          ENDIF
       ELSE
          WRITE(*,*) "could not find the variable contfrac in file ", filename
          STOP
       ENDIF
       WRITE(*,*) "XXX", MINVAL(contf), " < contf < ", MAXVAL(contf)
    ELSE
       WRITE(*,*) "Unknown model = ", model
       STOP
    ENDIF
    !
    ! Delete undef values
    !
    DO i=1,iimf
       DO j=1,jjmf
          IF ( contf(i,j) > 1.0 .OR. contf(i,j) < 0.0 ) THEN
             contf(i,j) = 0.0
          ENDIF
       ENDDO
    ENDDO
    !
    iret = NF90_CLOSE(fid)
    !
  END SUBROUTINE readcontfrac
  !!
  !
  !================================================================================================================================
  !! SUBROUTINE   : topo_getsize
  !!
  !> BRIEF         This subroutine gets the size of the topography field
  !!
  !! DESCRIPTION  : None
  !!
  !! RECENT CHANGE(S): None
  !!
  !! MAIN OUTPUT VARIABLE(S):
  !!
  !! REFERENCES   : None
  !!
  !! FLOWCHART    : None
  !! \n
  !_================================================================================================================================
  SUBROUTINE topo_getsize(filename, iim_full, jjm_full)
    !
    USE defprec
    USE netcdf
    !
    ! ARGUMENTS
    !
    CHARACTER(LEN=*), INTENT(in) :: filename
    INTEGER(i_std), INTENT(out)  :: iim_full, jjm_full
    !
    ! LOCAL
    !
    INTEGER(i_std) :: fid, iret, xid, yid, ndims, nvars, i, lll
    CHARACTER(LEN=20) :: dimname
    !
    iret = NF90_OPEN(filename, NF90_NOWRITE, fid)
    IF (iret /= NF90_NOERR) THEN
       WRITE(*,*) "Error opening ", filename
       STOP
    ENDIF
    iret = NF90_INQUIRE (fid, nDimensions=ndims, nVariables=nvars)
    DO i=1,ndims
       iret = NF90_INQUIRE_DIMENSION (fid, i, name=dimname, len=lll)
       CALL change_to_lowercase(dimname)
       IF ( INDEX(dimname,"x") > 0 ) THEN
          iim_full = lll
       ELSE IF ( INDEX(dimname,"y") > 0 ) THEN
          jjm_full = lll
       ENDIF
    ENDDO
    !
    iret = NF90_CLOSE(fid)
    !
  END SUBROUTINE topo_getsize
  !================================================================================================================================
  !! SUBROUTINE   : topo_getvar
  !!
  !> BRIEF         This subroutine gets the topography out of the file
  !!
  !! DESCRIPTION  : None
  !!
  !! RECENT CHANGE(S): None
  !!
  !! MAIN OUTPUT VARIABLE(S):
  !!
  !! REFERENCES   : None
  !!
  !! FLOWCHART    : None
  !! \n
  !_================================================================================================================================
  SUBROUTINE topo_getvar(filename, iim, jjm, lon, lat, etopo)
    !
    USE defprec
    USE netcdf
    !
    ! ARGUMENTS
    !
    CHARACTER(LEN=*), INTENT(in)                 :: filename
    INTEGER(i_std), INTENT(in)                   :: iim, jjm
    REAL(r_std), DIMENSION(iim,jjm), INTENT(out) :: lon 
    REAL(r_std), DIMENSION(iim,jjm), INTENT(out) :: lat
    REAL(r_std), DIMENSION(iim,jjm), INTENT(out) :: etopo
    !
    ! LOCAL
    !
    INTEGER(i_std) :: fid, iret, xid, yid, vid
    INTEGER(i_std) :: i, j, iv, ndimsvar, ndims, nvars
    CHARACTER(LEN=20) :: varname, units
    INTEGER(i_std), DIMENSION(1) :: il
    REAL(r_std), ALLOCATABLE, DIMENSION(:)   :: lon_loc, lon_read, lat_read
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: etopo_loc
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:) :: dimids, londim_ids, latdim_ids
    !
    ALLOCATE(lon_loc(iim))
    ALLOCATE(lon_read(iim))
    ALLOCATE(lat_read(jjm))
    ALLOCATE(etopo_loc(iim,jjm))
    ALLOCATE(dimids(NF90_MAX_VAR_DIMS), londim_ids(NF90_MAX_VAR_DIMS), latdim_ids(NF90_MAX_VAR_DIMS))
    !
    iret = NF90_OPEN(filename, NF90_NOWRITE, fid)
    IF (iret /= NF90_NOERR) THEN
       WRITE(*,*) "Error opening ", filename
       STOP
    ENDIF
    !
    iret = NF90_INQUIRE (fid, nDimensions=ndims, nVariables=nvars)
    !
    ! Identify the variables
    !
    xid = -1
    yid = -1
    vid = -1
    !
    DO iv=1,nvars
       iret = NF90_INQUIRE_VARIABLE(fid, iv, name=varname, ndims=ndimsvar, dimids=dimids)
       iret = NF90_GET_ATT(fid, iv, 'units', units)
       CALL change_to_lowercase(units)
       !
       IF ( INDEX(units,"east") > 0 ) THEN
          xid = iv
       ELSE IF (INDEX(units,"north") > 0 ) THEN
          yid = iv
       ELSE IF ( INDEX(units,"meters") > 0 ) THEN
          vid = iv
       ENDIF
    ENDDO
    !
    ! Get variables now that they have been identified
    !
    IF ( xid > 0 ) THEN
       iret = NF90_GET_VAR(fid, xid, lon_read)
    ELSE
       WRITE(*,*) "cound not find the longitude in file ", filename
       STOP
    ENDIF
    IF ( yid > 0 ) THEN
       iret = NF90_GET_VAR(fid, yid, lat_read)
    ELSE
       WRITE(*,*) "cound not find the latitude in file ", filename
       STOP
    ENDIF
    IF ( vid > 0 ) THEN
       iret = NF90_GET_VAR(fid, vid, etopo_loc)
    ELSE
       WRITE(*,*) "cound not find the topography in file ", filename
       STOP
    ENDIF
    iret = NF90_CLOSE(fid)
    !
    ! Shift the topography around the longitude to get a -180:180 range 
    ! Only if needed !
    !
    IF ( MAXVAL(lon_read) > 180. ) THEN
       ! 
       lon_loc(:) = lon_read(:)
       !
       DO i=1,iim
          IF (lon_loc(i) > 180. ) THEN
             lon_loc(i) = lon_loc(i)-360.
          ENDIF
       ENDDO
       !
       DO i=1,iim
          il = MINLOC(lon_loc)
          lon_read(i) = lon_loc(il(1))
          etopo(i,:) = etopo_loc(il(1),:)
          lon_loc(il(1)) = 9999999.999999
       ENDDO
       !
       DO j=1,jjm
          lon(:,j) = lon_read(:)
       ENDDO
       DO i=1,iim
          lat(i,:) = lat_read(:)
       ENDDO
    ELSE
       DO j=1,jjm
          lon(:,j) = lon_read(:)
       ENDDO
       DO i=1,iim
          lat(i,:) = lat_read(:)
       ENDDO
       etopo(:,:) = etopo_loc(:,:)
    ENDIF
    !
    !
    DEALLOCATE(lon_read)
    DEALLOCATE(lon_loc)
    DEALLOCATE(lat_read)
    DEALLOCATE(etopo_loc)
    !
  END SUBROUTINE topo_getvar
  !
  !================================================================================================================================
  !! SUBROUTINE   : interpol
  !!
  !> BRIEF         This subroutine interpolates topography to new resolution.
  !!
  !! DESCRIPTION  : None
  !!
  !! RECENT CHANGE(S): None
  !!
  !! MAIN OUTPUT VARIABLE(S):
  !!
  !! REFERENCES   : None
  !!
  !! FLOWCHART    : None
  !! \n
  !_================================================================================================================================

  SUBROUTINE interpol(iim_fine, jjm_fine, lon_fine, lat_fine, orog, iim, jjm, lon, lat, neworog)
    !
    USE constantes_var
    USE ioipsl_para
    !
    IMPLICIT NONE
    ! Parameter
    !INTEGER(i_std), PARAMETER                                  :: nbvmax = 800
    INTEGER(i_std), PARAMETER                                  :: nbvmax = 30000
    !
    INTEGER(i_std)                                             :: iim_fine, jjm_fine, iim, jjm
    REAL(r_std), DIMENSION (iim_fine,jjm_fine)                 :: lon_fine
    REAL(r_std), DIMENSION (iim_fine,jjm_fine)                 :: lat_fine
    REAL(r_std), DIMENSION (iim_fine,jjm_fine)                 :: orog
    REAL(r_std), DIMENSION (iim,jjm)                           :: lon
    REAL(r_std), DIMENSION (jjm,iim)                           :: lat
    REAL(r_std), DIMENSION (iim,jjm)                           :: neworog
    !
    REAL(r_std), DIMENSION (nbvmax)                            :: area, topo
    !
    INTEGER(i_std)                                             :: ip, jp, ig, jg, i, fopt, lastjp
    REAL(r_std)                                                :: dx, dy, ax, ay, ox, lon_up, lon_low, lat_up, lat_low, sgn
    REAL(r_std)                                                :: totarea, landarea, height
    REAL(r_std), ALLOCATABLE, DIMENSION (:,:)                  :: lon_ful, lat_ful, laup_rel, loup_rel, lalow_rel, lolow_rel
    REAL(r_std)                                                :: lonrel, lolowrel, louprel, coslat
    !
    ! Allocate arrays
    !
    ALLOCATE (laup_rel(iim_fine,jjm_fine))
    ALLOCATE (loup_rel(iim_fine,jjm_fine))
    ALLOCATE (lalow_rel(iim_fine,jjm_fine))
    ALLOCATE (lolow_rel(iim_fine,jjm_fine))
    ALLOCATE (lat_ful(iim_fine+2,jjm_fine+2))
    ALLOCATE (lon_ful(iim_fine+2,jjm_fine+2))
    !
    ! Duplicate the border assuming we have a global grid going from west to east
    !
    DO jp=2,jjm_fine+1
       lon_ful(2:iim_fine+1,jp) = lon_fine(1:iim_fine,jp)
    ENDDO
    DO ip=2,iim_fine+1
       lat_ful(ip,2:jjm_fine+1) = lat_fine(ip,1:jjm_fine)
    ENDDO
    !
    IF ( lon_fine(iim_fine,1) .LT. lon_ful(2,2)) THEN
       lon_ful(1,2:jjm_fine+1) = lon_fine(iim_fine,:)
       lat_ful(1,2:jjm_fine+1) = lat_fine(1,1:jjm_fine)
    ELSE
       lon_ful(1,2:jjm_fine+1) = lon_fine(iim_fine,:)-360.
       lat_ful(1,2:jjm_fine+1) = lat_fine(1,1:jjm_fine)
    ENDIF

    IF ( lon_fine(1,1) .GT. lon_ful(iim_fine+1,2)) THEN
       lon_ful(iim_fine+2,2:jjm_fine+1) = lon_fine(1,1:jjm_fine)
       lat_ful(iim_fine+2,2:jjm_fine+1) = lat_fine(iim_fine,1:jjm_fine)
    ELSE
       lon_ful(iim_fine+2,2:jjm_fine+1) = lon_fine(1,1:jjm_fine)+360
       lat_ful(iim_fine+2,2:jjm_fine+1) = lat_fine(iim_fine,1:jjm_fine)
    ENDIF
    !
    sgn = lat_fine(1,1)/ABS(lat_fine(1,1))
    lat_ful(2:iim_fine+1,1) = sgn*180 - lat_fine(1,1)
    sgn = lat_fine(iim_fine,jjm_fine)/ABS(lat_fine(iim_fine,jjm_fine))
    lat_ful(2:iim_fine+1,jjm_fine+2) = sgn*180 - lat_fine(iim_fine,jjm_fine)
    lat_ful(1,1) = lat_ful(iim_fine+1,1)
    lat_ful(iim_fine+2,1) = lat_ful(2,1)
    lat_ful(1,jjm_fine+2) = lat_ful(iim_fine+1,jjm_fine+2)
    lat_ful(iim_fine+2,jjm_fine+2) = lat_ful(2,jjm_fine+2)
    !
    ! Add the longitude lines to the top and bottom
    !
    lon_ful(:,1) = lon_ful(:,2) 
    lon_ful(:,jjm_fine+2) = lon_ful(:,jjm_fine+1) 
    !
    ! Get the upper and lower limits of each grid box
    !
    DO ip=1,iim_fine
       DO jp=1,jjm_fine
          loup_rel(ip,jp) =MAX(0.5*(lon_ful(ip,jp+1)+lon_ful(ip+1,jp+1)), 0.5*(lon_ful(ip+1,jp+1)+lon_ful(ip+2,jp+1)))
          lolow_rel(ip,jp) =MIN(0.5*(lon_ful(ip,jp+1)+lon_ful(ip+1,jp+1)), 0.5*(lon_ful(ip+1,jp+1)+lon_ful(ip+2,jp+1)))
          laup_rel(ip,jp) =MAX(0.5*(lat_ful(ip+1,jp)+lat_ful(ip+1,jp+1)), 0.5*(lat_ful(ip+1,jp+1)+lat_ful(ip+1,jp+2)))
          lalow_rel(ip,jp) =MIN(0.5*(lat_ful(ip+1,jp)+lat_ful(ip+1,jp+1)), 0.5*(lat_ful(ip+1,jp+1)+lat_ful(ip+1,jp+2)))
       ENDDO
    ENDDO
    !
    ! Now we take each grid point and find out which values from the forcing we need to average
    !
    dx = ABS(lon(2,1) - lon(1,1))
    dy = ABS(lat(1,2) - lat(1,1))
    !
    DO ig = 1,iim
       DO jg = 1,jjm
          !
          !
          lon_up = lon(ig,jg) + dx/2.0
          lon_low = lon(ig,jg) - dx/2.0
          !
          !
          lat_up = lat(ig,jg) + dy/2.0
          lat_low = lat(ig,jg) - dy/2.0
          !
          !  Find the grid boxes from the data that go into the model's boxes
          !  We still work as if we had a regular grid ! Well it needs to be localy regular so
          !  so that the longitude at the latitude of the last found point is close to the one of the next point.
          !
          totarea = 0.0
          landarea = 0.0
          height = 0.0
          !
          fopt = 0
          lastjp = 1
          DO ip=1,iim_fine
             !
             DO jp = 1, jjm_fine
                !  Either the center of the data grid point is in the interval of the model grid or
                !  the East and West limits of the data grid point are on either sides of the border of
                !  the data grid.
                !
                !  To do that correctly we have to check if the grid box sits on the date-line.
                !
                IF ( lon_low < -180.0 ) THEN
                   lonrel = MOD( lon_fine(ip,jp) - 360.0, 360.0)
                   lolowrel = MOD( lolow_rel(ip,jp) - 360.0, 360.0)
                   louprel = MOD( loup_rel(ip,jp) - 360.0, 360.0)
                   !
                ELSE IF ( lon_up > 180.0 ) THEN
                   lonrel = MOD( 360. - lon_fine(ip,jp), 360.0)
                   lolowrel = MOD( 360. - lolow_rel(ip,jp), 360.0)
                   louprel = MOD( 360. - loup_rel(ip,jp), 360.0)
                ELSE
                   lonrel = lon_fine(ip,jp)
                   lolowrel = lolow_rel(ip,jp)
                   louprel = loup_rel(ip,jp)
                ENDIF
                !
                !
                IF ( lonrel > lon_low .AND. lonrel < lon_up .OR. &
                     & lolowrel < lon_low .AND.  louprel > lon_low .OR. &
                     & lolowrel < lon_up  .AND.  louprel > lon_up ) THEN
                   !  
                   !
                   ! Now that we have the longitude let us find the latitude
                   !
                   IF ( lat_fine(ip,jp) > lat_low .AND. lat_fine(ip,jp) < lat_up .OR. &
                        & lalow_rel(ip,jp) < lat_low .AND. laup_rel(ip,jp) > lat_low .OR.&
                        & lalow_rel(ip,jp) < lat_up .AND. laup_rel(ip,jp) > lat_up) THEN
                      !
                      lastjp = jp
                      !
                      fopt = fopt + 1
                      IF ( fopt .GT. nbvmax) THEN
                         WRITE(*,*) 'Please increase nbvmax in subroutine interpol', ig, jg
                         WRITE(*,*) "nbvmax = ", nbvmax
                         STOP
                      ELSE
                         !
                         ! If we sit on the date line we need to do the same transformations as above.
                         ! 
                         IF ( lon_low < -180.0 ) THEN
                            lolowrel = MOD( lolow_rel(ip,jp) - 360.0, 360.0)
                            louprel = MOD( loup_rel(ip,jp) - 360.0, 360.0)
                            !
                         ELSE IF ( lon_up > 180.0 ) THEN
                            lolowrel = MOD( 360. - lolow_rel(ip,jp), 360.0)
                            louprel = MOD( 360. - loup_rel(ip,jp), 360.0)
                         ELSE
                            lolowrel = lolow_rel(ip,jp)
                            louprel = loup_rel(ip,jp)
                         ENDIF
                         !
                         ! Get the area of the fine grid in the model grid
                         !
                         coslat = MAX( COS( lat_fine(ip,jp) * pi/180. ), 0.001 )
                         ax = (MIN(lon_up,louprel)-MAX(lon_low, lolowrel))*pi/180. * R_Earth * coslat
                         ay = (MIN(lat_up, laup_rel(ip,jp))-MAX(lat_low,lalow_rel(ip,jp)))*pi/180. * R_Earth
                         !
                         totarea = totarea + ax*ay
                         IF ( orog(ip, jp) > 0.0 ) THEN
                            landarea = landarea + ax*ay
                            height = height + ax*ay*orog(ip,jp)
                         ENDIF
                         !
                      ENDIF
                   ENDIF
                   !
                ENDIF
                !
             ENDDO
             !
          ENDDO
          !
          IF ( landarea > 0 .AND. landarea > 0.1*totarea ) THEN
             neworog(ig,jg) = height/landarea
          ELSE
             neworog(ig,jg) = 0.0
          ENDIF
          !
       ENDDO
       !
    ENDDO
  END SUBROUTINE interpol

  SUBROUTINE change_to_lowercase (str)
    !---------------------------------------------------------------------
    !- Converts a string into lower case letters
    !---------------------------------------------------------------------
    IMPLICIT NONE
    !-
    CHARACTER(LEN=*) :: str
    !-
    INTEGER :: i,ic
    !---------------------------------------------------------------------
    DO i=1,LEN_TRIM(str)
       ic = IACHAR(str(i:i))
       IF ( (ic >= 65).AND.(ic <= 90) )  str(i:i) = ACHAR(ic+32)
    ENDDO
    !--------------------------
  END SUBROUTINE change_to_lowercase

END MODULE getlandseamask
