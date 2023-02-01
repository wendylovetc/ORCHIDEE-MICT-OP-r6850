!  ==============================================================================================================================\n
!  MODULE globgrd : This module is dedicated to managing the spatial grid of the forcing. It can either read a file
!                 containing the grid information, as is the case for WRF forcing, or obtain the grid from the forcing files.
!                 The module has also the possibility to create a grid description files for certain applications like
!                 for instance in a coupling of ORCHIDEE through OASIS. 
!                 For this purpose the module provides 4 subroutines :
!                 globgrd_getdomsz : This routine allows to get the domain size of the forcing based on a file it will explore.
!                 globgrd_getgrid : This routine extracts the coordinates and land/sea mask from the domain files.
!                 globgrd_writevar : Writes a variables into a netCDF file which can then be analysed to verify that the
!                                    forcing grid was well read and interpreted by this module.
!                 globgrd_writegrid : Write a grid description file which has the WRF flavor. It allows to exchange grid information
!                                     between an atmospheric model (mostly a driver !) and ORCHIDEE which are coupled through OASIS.
!
!  CONTACT      : jan.polcher@lmd.jussieu.fr
!
!  LICENCE      : IPSL (2016)
!  This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        
!! 
!! RECENT CHANGE(S): None 
!! 
!! REFERENCE(S) : None
!! 
!_ ================================================================================================================================
MODULE globgrd
  !
  !
  USE defprec
  USE netcdf
  !
  USE ioipsl
  !
  USE grid
  USE forcing_tools
  !
  IMPLICIT NONE
  !
  PRIVATE
  PUBLIC :: globgrd_getdomsz, globgrd_getgrid, globgrd_writevar, globgrd_writegrid
  !
  !
  LOGICAL, SAVE  :: is_forcing_file=.FALSE.
  !
CONTAINS
!!
!!  =============================================================================================================================
!! SUBROUTINE: globgrd_getdomsz
!!
!>\BRIEF  This routine allows to get the domain size of the forcing based on a file it will explore.
!!
!! DESCRIPTION: The routine opens the file and explores it. It can either be a forcing file or a grid description
!!              file from WRF. Progressively this should be opened to other ways of describing the grid over which
!!              the forcing is provided.
!!              The routing will return the sizes in I and J and the number of land points.
!!              The zooming interval is also provided so that only the dimensions over the domain used can be computed.
!!
!! \n
!_ ==============================================================================================================================
!!
  !---------------------------------------------------------------------
  !-
  !- 
  !---------------------------------------------------------------------
  SUBROUTINE globgrd_getdomsz(filename, iim, jjm, nbland, model_guess, fid, forcingfile, zoom_lon, zoom_lat)
    !
    ! INPUT
    !
    CHARACTER(LEN=*), INTENT(in)  :: filename
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: forcingfile(:)
    REAL(r_std), DIMENSION(2), INTENT(in), OPTIONAL :: zoom_lon, zoom_lat
    !
    ! OUTPUT
    !
    INTEGER(i_std), INTENT(out)    :: fid
    INTEGER(i_std), INTENT(out)    :: iim, jjm, nbland
    CHARACTER(LEN=*), INTENT(out)  :: model_guess
    !
    ! LOCAL
    !
    INTEGER(i_std) :: iret, ndims, nvars, nb_atts, id_unlim, iv, lll
    INTEGER(i_std) :: iim_full, jjm_full, nbland_full
    CHARACTER(LEN=20) :: axname, varname
    CHARACTER(LEN=120) :: tmpfile
    REAL(r_std), DIMENSION(2) :: loczoom_lon, loczoom_lat
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: mask
    !
    ! Set default values against which we can test
    !
    iim = -1 
    jjm = -1
    !
    ! Verify the grid file name
    ! 
    IF ( INDEX(filename,"NONE") >= 1 ) THEN
       is_forcing_file=.TRUE.
       IF ( PRESENT(forcingfile) ) THEN
          tmpfile=forcingfile(1)
       ELSE
          CALL ipslerr (3,'globgrd_getdomsz',"Error No grid file provided :",tmpfile, " ")
       ENDIF
    ELSE
       is_forcing_file=.FALSE.
       tmpfile=filename
    ENDIF
    !
    ! Verify that the zomm is provided. Else choose the entire globe
    !
    IF ( PRESENT(zoom_lon) .AND. PRESENT(zoom_lat) ) THEN
       loczoom_lon = zoom_lon
       loczoom_lat = zoom_lat
    ELSE
       loczoom_lon(1) = -180.0
       loczoom_lon(2) = 180.0
       loczoom_lat(1) = -90.0
       loczoom_lat(2) = 90.0
    ENDIF
    !
    ! Open the correct file
    !
    iret = NF90_OPEN (tmpfile, NF90_NOWRITE, fid)
    IF (iret /= NF90_NOERR) THEN
       CALL ipslerr (3,'globgrd_getdomsz',"Error opening the grid file :",tmpfile, " ")
    ENDIF
    !
    !
    iret = NF90_INQUIRE (fid, nDimensions=ndims, nVariables=nvars, &
         nAttributes=nb_atts, unlimitedDimId=id_unlim)
    IF (iret /= NF90_NOERR) THEN
       CALL ipslerr (3,'globgrd_getdomsz',"Error in NF90_INQUIRE :",tmpfile, " ")
    ENDIF
    !
    DO iv=1,ndims
       !
       iret = NF90_INQUIRE_DIMENSION (fid, iv, name=axname, len=lll)
       IF (iret /= NF90_NOERR) THEN
          CALL ipslerr (3,'globgrd_getdomsz',"Could not get size of dimension :"," "," ")
       ENDIF
       !
       ! This can be refined by testing the actual grid found in the file.
       !
       SELECT CASE(axname)
          !
          !! Coordinate variables used by WRF.
       CASE("west_east")
          iim = lll
          model_guess = "WRF"
       CASE("south_north")
          jjm = lll
          model_guess = "WRF"
          !
          !! Variables used in WFDEI
       CASE("lon")
          iim = lll
          model_guess = "regular"
       CASE("lat")
          jjm = lll
          model_guess = "regular"
       CASE("nbland")
          nbland = lll
          !
          !! Variables used by CRU-NCEP
       CASE("nav_lon")
          iim = lll
          model_guess = "regular"
       CASE("nav_lat")
          jjm = lll
          model_guess = "regular"
       CASE("land")
          nbland = lll
       END SELECT
    ENDDO
    !
    ! If we have a WRF file we need to count the number of land points.
    !
    IF (  model_guess == "WRF" ) THEN

       ALLOCATE(mask(iim,jjm))

       varname = "LANDMASK"
       iret = NF90_INQ_VARID (fid, varname, iv)
       IF (iret /= NF90_NOERR) THEN
          CALL ipslerr (3,'globgrd_getdomsz',"Could not find variable ", varname," ")
       ELSE
          iret = NF90_GET_VAR (fid,iv,mask)
       ENDIF

       nbland = COUNT(mask > 0.5)

    ENDIF
    !
    !
    ! If we are in the case of a forcing file a few functions from forcing_tools need to be called
    ! so that the file is analysed with the tools of the forcing module.
    !
    IF ( is_forcing_file ) THEN
       !
       ! Because we are re-using routines from the forcing module, we have to
       ! close the file. It will be opened again by the forcing module.
       !
       iret = NF90_CLOSE(fid)
       IF (iret /= NF90_NOERR) THEN
          CALL ipslerr (3,'globgrd_getdomzz',"Error closing the output file :",filename, " ")
       ENDIF
       !
       ! This has to be a regular grid. A more clever clasification of files will be needed.
       model_guess = "regular"
       !
       CALL forcing_getglogrid (1, forcingfile, iim_full, jjm_full, nbland_full, .TRUE.)
       WRITE(*,*) forcingfile, "Forcing file with dimensions : ", iim_full, jjm_full, nbland_full
       !
       CALL forcing_zoomgrid (loczoom_lon, loczoom_lat, forcingfile(1), .TRUE.)
       !
       CALL forcing_givegridsize (iim, jjm, nbland)
       !
    ENDIF
    !
    ! Do a final test to see if we got the information needed.
    !
    IF ( iim < 0 .OR. jjm < 0 ) THEN
       CALL ipslerr (3,'globgrd_getdomsz',"Could not get the horizontal size of the domaine out of the file",&
            & filename,"Are you sure that the case for this type of file is foreseen ? ")
    ENDIF
    !
    !
  END SUBROUTINE globgrd_getdomsz
!!
!!  =============================================================================================================================
!! SUBROUTINE: globgrd_getgrid
!!
!>\BRIEF        This routine extracts the coordinates and land/sea mask from the domain files.     
!!
!! DESCRIPTION: The domain size is provided together with the netCDF file ID so that the main information can
!!              be extracted from the file. We will read the longitude, latitude, land/sea mask and calendar.
!!              This allows to set-up ORCHIDEE. We also provide the corners of the grid-boxes as this is needed
!!              for setting-up OASIS but is computed more correctly in grid.f90 for ORCHIDEE. 
!!              This routine is only an interface to globgrd_getwrf, globgrd_getregular and forcing_givegrid.
!!              forcing_givegrid is an interface to the forcing_tools.f90 module so that we are certain to have
!!              the same grid information between both modules.
!!
!! \n
!_ ==============================================================================================================================
!!
  !---------------------------------------------------------------------
  !-
  !- 
  !---------------------------------------------------------------------
  SUBROUTINE globgrd_getgrid(fid, iim, jjm, nbland, model_guess, lon, lat, mask, area, corners, &
       &                     lindex, contfrac, calendar)
    !
    !
    ! This subroutine only switched between routines to extract and compte the grid data needed for 
    ! ORCHIDEE.
    !
    !
    ! INPUT
    !
    INTEGER(i_std), INTENT(in)   :: fid
    INTEGER(i_std), INTENT(in)   :: iim, jjm, nbland
    CHARACTER(LEN=*), INTENT(in) :: model_guess
    !
    ! OUTPUT
    !
    REAL(r_std),DIMENSION(iim,jjm), INTENT(out)     :: lon, lat, mask, area
    REAL(r_std),DIMENSION(iim,jjm,4,2), INTENT(out) :: corners
    INTEGER(i_std), DIMENSION(nbland), INTENT(out)  :: lindex
    REAL(r_std),DIMENSION(nbland), INTENT(out)      :: contfrac
    CHARACTER(LEN=20), INTENT(out)                  :: calendar
    !
    SELECT CASE(model_guess)

    CASE("WRF")
       CALL globgrd_getwrf(fid, iim, jjm, nbland, lon, lat, mask, area, corners, &
            &               lindex, contfrac, calendar)
    CASE("regular")
       IF ( .NOT. is_forcing_file ) THEN
          CALL globgrd_getregular(fid, iim, jjm, nbland, lon, lat, mask, area, corners, &
               &                   lindex, contfrac, calendar)
       ELSE
          CALL forcing_givegrid(lon, lat, mask, area, corners, lindex, contfrac, calendar)
          CALL forcing_close()
       ENDIF
    CASE DEFAULT
       CALL ipslerr (3,'globgrd_getgrid',"The model/grid type we guessed is not recognized here. model_guess =",&
            & model_guess,"Have you used the right file and are you sure that this case is foreseen ? ")
    END SELECT
    !
    !
  END SUBROUTINE globgrd_getgrid
!!
!!  =============================================================================================================================
!! SUBROUTINE: globgrd_regular
!!
!>\BRIEF       The routine to obtain regular grids from the file.     
!!
!! DESCRIPTION:	  Read the regular grid and its information from the opened file (fid).
!!
!! \n
!_ ==============================================================================================================================
!!
  SUBROUTINE globgrd_getregular(fid, iim, jjm, nbland, lon, lat, mask, area, corners, &
       &                     lindex, contfrac, calendar)
    !
    USE defprec
    USE netcdf
    !
    ! INPUT
    !
    INTEGER(i_std), INTENT(in)   :: fid
    INTEGER(i_std), INTENT(in)   :: iim, jjm, nbland
    !
    ! OUTPUT
    !
    REAL(r_std),DIMENSION(iim,jjm), INTENT(out)     :: lon, lat, mask, area
    REAL(r_std),DIMENSION(iim,jjm,4,2), INTENT(out) :: corners
    INTEGER(i_std), DIMENSION(nbland), INTENT(out)  :: lindex
    REAL(r_std),DIMENSION(nbland), INTENT(out)      :: contfrac
    CHARACTER(LEN=20), INTENT(out)                  :: calendar
    !
    ! LOCAL
    !
    INTEGER(i_std) :: iret, iv, nvars, varndim
    INTEGER(i_std) :: i, j
    CHARACTER(LEN=20) :: varname
    INTEGER(i_std), DIMENSION(4) :: vardims
    REAL(r_std) :: dx
    !
    ! Set some default values agains which we can check 
    !
    lon(:,:) = val_exp
    lat(:,:) = val_exp
    mask(:,:) = val_exp
    area(:,:) = val_exp
    corners(:,:,:,:) = val_exp
    !
    lindex(:) = INT(val_exp)
    contfrac(:) = val_exp
    !
    ! Get the global attributes from grid file
    !
    iret = NF90_GET_ATT(fid, NF90_GLOBAL, 'calendar', calendar)
    IF (iret /= NF90_NOERR) THEN
       CALL ipslerr (3,'globgrd_getregular',"Could not read the calendar in grid file.", " ", " ")
    ENDIF
    !
    iret = NF90_INQUIRE (fid, nVariables=nvars)
    !
    DO iv = 1,nvars
       !
       iret = NF90_INQUIRE_VARIABLE(fid, iv, name=varname, ndims=varndim, dimids=vardims)
       !
       !
       SELECT CASE(varname)
          !
       CASE("longitude")
          IF (varndim == 1 ) THEN
             DO j=1,jjm
                iret = NF90_GET_VAR(fid, iv, lon(:,j))
             ENDDO
          ELSE IF (varndim == 2 ) THEN
             iret = NF90_GET_VAR(fid, iv, lon)
          ELSE
             CALL ipslerr (3,'globgrd_getregular',"Longitude cannot have more than 2 dimensions","","")
          ENDIF

       CASE ("latitude")
          IF (varndim == 1 ) THEN
             DO i=1,iim
                iret = NF90_GET_VAR(fid, iv, lat(i,:))
             ENDDO
          ELSE IF (varndim == 2 ) THEN
             iret = NF90_GET_VAR(fid, iv, lon)
          ELSE
             CALL ipslerr (3,'globgrd_getregular',"Latitude cannot have more than 2 dimensions","","")
          ENDIF

       CASE ("mask")
          IF (varndim /= 2 ) THEN
             CALL ipslerr (3,'globgrd_getregular',"mask needs to have 2 dimensions","","")
          ELSE
             iret = NF90_GET_VAR (fid, iv, mask)
          ENDIF

       CASE ("areas")
          IF (varndim /= 2 ) THEN
             CALL ipslerr (3,'globgrd_getregular',"Areas needs to have 2 dimensions","","")
          ELSE
             iret = NF90_GET_VAR (fid, iv, area)
          ENDIF

       CASE ("corners")
          IF (varndim /= 4 ) THEN
             CALL ipslerr (3,'globgrd_getregular',"corners needs to have 4 dimensions","","")
          ELSE
             iret = NF90_GET_VAR (fid, iv, corners)
          ENDIF

       CASE ("landindex")
          IF (varndim /= 1 ) THEN
             CALL ipslerr (3,'globgrd_getregular',"landindex is the list of continental points to be gathered", &
                  &          "Thus it can only have 1 dimensions","")
          ELSE
             iret = NF90_GET_VAR (fid, iv, lindex)
          ENDIF

       CASE ("contfrac")
          IF (varndim /= 1 ) THEN
             CALL ipslerr (3,'globgrd_getregular',"Contfrac needs to be a gathered variable", &
                  &          "thus it needs only 1 dimensions","")
          ELSE
             iret = NF90_GET_VAR (fid, iv, contfrac)
          ENDIF

       END SELECT
       !
    ENDDO
    !
    !
    iret = NF90_CLOSE(fid)
    !
    ! Verify that we have al the variables needed to describe the ORCHIDEE grid
    !
    IF ( ANY( lon(:,:) == val_exp ) ) THEN
       CALL ipslerr (3,'globgrd_getregular',"The longitude of the ORCHIDEE grid could not be extracted from the",&
            & "grid definition file","")
    ENDIF
    !
    IF ( ANY( lat(:,:) == val_exp ) ) THEN
       CALL ipslerr (3,'globgrd_getregular',"The latitude of the ORCHIDEE grid could not be extracted from the",&
            & "grid definition file","")
    ENDIF
    !
    IF ( ANY( lindex(:) == INT(val_exp) ) ) THEN
       CALL ipslerr (3,'globgrd_getregular',"The lindex of the ORCHIDEE grid could not be extracted from the",&
            & "grid definition file","")
    ENDIF
    !
    IF ( ALL( mask(:,:) == val_exp ) ) THEN
       CALL ipslerr (3,'globgrd_getregular',"The land mask of the ORCHIDEE grid could not be extracted from the",&
            & "grid definition file","")
    ELSE IF (MAXVAL(mask) > 1 ) THEN
       CALL ipslerr (2,'globgrd_getregular',"We have a special case for the mask which needs to be treated.",&
            & "The field contains the indices of the land points on a compressed grid.","So we replace them with 1 or 0.")
       DO i=1,iim
          DO j=1,jjm
             IF ( mask(i,j) > iim*jjm ) THEN
                mask(i,j) = 0
             ELSE
                mask(i,j) = 1
             ENDIF
          ENDDO
       ENDDO
    ENDIF
    !
    IF ( ANY( contfrac(:) == val_exp ) ) THEN
       CALL ipslerr (2,'globgrd_getregular',"The continental fraction of the ORCHIDEE grid could not be extracted from the",&
            & "grid definition file","Thus on all land points it is set to 1.")
       contfrac(:) = 1.
    ENDIF
    !
    IF ( ANY( corners(:,:,:,:) == val_exp ) ) THEN
       CALL ipslerr (3,'globgrd_getregular',"The corners for the ORCHIDEE grid could not be extracted from the",&
            & "grid definition file","As we have to assume a very general grid we cannot do anything !")
    ENDIF
    !
    !
  END SUBROUTINE globgrd_getregular
!!
!!  =============================================================================================================================
!! SUBROUTINE: globgrd_getwrf
!!
!>\BRIEF       Routine to read the WRF grid description file.
!!
!! DESCRIPTION:	Read the WRF grid and its information from the opened file (fid) and convert
!!              it to the variables needed by ORCHIDEE.
!!
!! \n
!_ ==============================================================================================================================
!!
  SUBROUTINE globgrd_getwrf(fid, iim, jjm, nbland, lon, lat, mask, area, corners, &
       &                     lindex, contfrac, calendar)
    !
    USE defprec
    USE netcdf
    !
    ! INPUT
    !
    INTEGER(i_std), INTENT(in)   :: fid
    INTEGER(i_std), INTENT(in)   :: iim, jjm, nbland
    !
    ! OUTPUT
    !
    REAL(r_std),DIMENSION(iim,jjm), INTENT(out)     :: lon, lat, mask, area
    REAL(r_std),DIMENSION(iim,jjm,4,2), INTENT(out) :: corners
    INTEGER(i_std), DIMENSION(nbland), INTENT(out)  :: lindex
    REAL(r_std),DIMENSION(nbland), INTENT(out)      :: contfrac
    CHARACTER(LEN=20), INTENT(out)                  :: calendar
    !
    ! LOCAL
    !
    INTEGER(i_std) :: i, ip, jp, k, iret, iv, nvars, varndim
    INTEGER(i_std) :: imm1, imp1
    CHARACTER(LEN=20) :: varname
    REAL(r_std),DIMENSION(iim,jjm) :: mapfac_x, mapfac_y
    INTEGER(i_std), DIMENSION(4) :: vardims
    INTEGER(i_std), DIMENSION(8) :: rose
    !
    !
    ! Set some default values agains which we can check 
    !
    lon(:,:) = val_exp
    lat(:,:) = val_exp
    mask(:,:) = val_exp
    area(:,:) = val_exp
    corners(:,:,:,:) = val_exp
    !
    lindex(:) = INT(val_exp)
    contfrac(:) = val_exp
    !
    calendar = 'gregorian'
    !
    !
    !  Init projection in grid.f90 so that it can be used later for projections.
    !
    CALL grid_initproj(fid, iim, jjm)
    !
    iret = NF90_INQUIRE (fid, nVariables=nvars)
    !
    IF (iret /= NF90_NOERR) THEN
       CALL ipslerr (3,'globgrd_getwrf',"Error inquiering variables from WRF grid file."," ", " ")
    ENDIF
    !
    DO iv = 1,nvars
       !
       iret = NF90_INQUIRE_VARIABLE(fid, iv, name=varname, ndims=varndim, dimids=vardims)
       !
       SELECT CASE(varname)
       !
       CASE("XLONG_M")
          iret = NF90_GET_VAR(fid, iv, lon)
          IF (iret /= NF90_NOERR) THEN
             CALL ipslerr (3,'globgrd_getwrf',"Could not read the longitude from the WRF grid file.", " ", " ")
          ENDIF
       CASE("XLAT_M")
          iret = NF90_GET_VAR(fid, iv, lat)
          IF (iret /= NF90_NOERR) THEN
             CALL ipslerr (3,'globgrd_getwrf',"Could not read the latitude from the WRF grid file.", " ", " ")
          ENDIF
       CASE("LANDMASK")
          iret = NF90_GET_VAR(fid, iv, mask)
          IF (iret /= NF90_NOERR) THEN
             CALL ipslerr (3,'globgrd_getwrf',"Could not read the land mask from the WRF grid file.", " ", " ")
          ENDIF
       CASE("MAPFAC_MX")
          iret = NF90_GET_VAR(fid, iv, mapfac_x)
          IF (iret /= NF90_NOERR) THEN
             CALL ipslerr (3,'globgrd_getwrf',"Could not read the land mask from the WRF grid file.", " ", " ")
          ENDIF
       CASE("MAPFAC_MY")
          iret = NF90_GET_VAR(fid, iv, mapfac_y)
          IF (iret /= NF90_NOERR) THEN
             CALL ipslerr (3,'globgrd_getwrf',"Could not read the land mask from the WRF grid file.", " ", " ")
          ENDIF
          !
       END SELECT
    ENDDO
    !
    ! Compute corners and area on the iimxjjm grid
    !
    DO ip=1,iim
       DO jp=1,jjm
          ! Corners
          CALL grid_tolola(ip+0.5, jp+0.5, corners(ip,jp,1,1), corners(ip,jp,1,2))
          CALL grid_tolola(ip+0.5, jp-0.5, corners(ip,jp,2,1), corners(ip,jp,2,2))
          CALL grid_tolola(ip-0.5, jp-0.5, corners(ip,jp,3,1), corners(ip,jp,3,2))
          CALL grid_tolola(ip-0.5, jp-0.5, corners(ip,jp,4,1), corners(ip,jp,4,2))
          !
       ENDDO
    ENDDO
    !
    ! Compute resolution and area on the gathered grid
    !
    k=0
    !
    DO jp=1,jjm
       DO ip=1,iim
          IF ( mask(ip,jp) > 0.5 ) THEN
             !
             k=k+1
             lindex(k) = (jp-1)*iim+ip 
             contfrac(k) = 1.0
             !
          ENDIF
       ENDDO
    ENDDO
    !
    iret = NF90_CLOSE (fid)
    IF (iret /= NF90_NOERR) THEN
       CALL ipslerr (3,'globgrd_getwrf',"Error closing the WRF grid file :", " ", " ")
    ENDIF
    !
  END SUBROUTINE globgrd_getwrf
!!
!!  =============================================================================================================================
!! SUBROUTINE: globgrd_writegrid
!!
!>\BRIEF      Allows to write the grid to a netDF file for later usage by the glogrid module.
!!
!! DESCRIPTION: This routine will write a grid description to a netCDF file. mask is on the iimxjjm grid while other
!! variables are on the gathered grid. 
!!
!! \n
!_ ==============================================================================================================================
!!
!
!
  SUBROUTINE globgrd_writegrid (gridfilename)
    !
    ! This routine will write a grid description to a netCDF file. mask is on the iimxjjm grid while other
    ! variables are on the gathered grid.
    !
    ! ARGUMENTS
    !
    CHARACTER(LEN=*), INTENT(in) :: gridfilename
    !
    ! LOCAL Grid description
    !
    INTEGER(i_std) :: iim, jjm, nblindex
    !
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)    :: lon, lat
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)    :: mask
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)    :: area
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:,:):: corners
    REAL(r_std), ALLOCATABLE, DIMENSION(:)      :: contfrac
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:)   :: lindex
    CHARACTER(LEN=20) :: calendar
    !
    ! LOCAL netCDF and helping variables
    !
    INTEGER(i_std) :: iret, fid, i
    INTEGER(i_std) :: lonid, latid, landdimid, resid, neighid, maskid, nbcornersid
    INTEGER(i_std) :: londimid, latdimid, contfracid, resolutionid, neighbourid
    INTEGER(i_std) :: landindexid, areaid, cornerid
    !
    ! Get the grid size from the forcing module
    !
    CALL forcing_givegridsize (iim, jjm, nblindex)
    WRITE(*,*) "Dimension of grid for forcing (iim,jjm,nblindex):", iim,jjm,nblindex
    !
    ! Allocate fields
    !
    ALLOCATE(lon(iim,jjm), lat(iim,jjm))
    ALLOCATE(mask(iim,jjm))
    ALLOCATE(area(iim,jjm))
    ALLOCATE(corners(iim,jjm,4,2))
    ALLOCATE(lindex(nblindex))
    ALLOCATE(contfrac(nblindex))
    !
    ! Get the actual grid from the forcing module
    !
    CALL forcing_givegrid(lon, lat, mask, area, corners, lindex, contfrac, calendar)
    ! 
    !
    iret = NF90_CREATE(gridfilename, NF90_WRITE, fid)
    IF (iret /= NF90_NOERR) THEN
       CALL ipslerr (3,'globgrd_writegrid',"Error opening the output file :",gridfilename, " ")
    ENDIF
    !
    ! Define dimensions
    !
    iret = NF90_DEF_DIM(fid,'lon',iim,londimid)
    iret = NF90_DEF_DIM(fid,'lat',jjm,latdimid)
    iret = NF90_DEF_DIM(fid,'nbland',nblindex,landdimid)
    iret = NF90_DEF_DIM(fid,'nbres',2,resid)
    iret = NF90_DEF_DIM(fid,'nbcorners',4,nbcornersid)
    !
    !
    ! We need to verify here that we have a regulat grid befor deciding if we write lon and lat in 1D or 2D !
    !
    !
    iret = NF90_DEF_VAR(fid,"longitude",NF90_REAL4,londimid,lonid)
    iret = NF90_PUT_ATT(fid,lonid,'standard_name',"longitude")
    iret = NF90_PUT_ATT(fid,lonid,'units',"degrees_east")
    iret = NF90_PUT_ATT(fid,lonid,'valid_min',MINVAL(lon))
    iret = NF90_PUT_ATT(fid,lonid,'valid_max',MAXVAL(lon))
    iret = NF90_PUT_ATT(fid,lonid,'long_name',"Longitude")
    !
    iret = NF90_DEF_VAR(fid,"latitude",NF90_REAL4,latdimid,latid)
    iret = NF90_PUT_ATT(fid,latid,'standard_name',"latitude")
    iret = NF90_PUT_ATT(fid,latid,'units',"degrees_north")
    iret = NF90_PUT_ATT(fid,latid,'valid_min',MINVAL(lat))
    iret = NF90_PUT_ATT(fid,latid,'valid_max',MAXVAL(lat))
    iret = NF90_PUT_ATT(fid,latid,'long_name',"Latitude")
    !
    iret = NF90_DEF_VAR(fid,"mask",NF90_REAL4,(/lonid,latid/),maskid)
    iret = NF90_PUT_ATT(fid,maskid,'standard_name',"mask")
    iret = NF90_PUT_ATT(fid,maskid,'units',"-")
    iret = NF90_PUT_ATT(fid,maskid,'valid_min',MINVAL(mask))
    iret = NF90_PUT_ATT(fid,maskid,'valid_max',MAXVAL(mask))
    iret = NF90_PUT_ATT(fid,maskid,'long_name',"Land surface mask")
    !
    iret = NF90_DEF_VAR(fid,"area",NF90_REAL4,(/lonid,latid/), areaid)
    iret = NF90_PUT_ATT(fid,areaid,'standard_name',"area")
    iret = NF90_PUT_ATT(fid,areaid,'units',"m*m")
    iret = NF90_PUT_ATT(fid,areaid,'valid_min',MINVAL(area))
    iret = NF90_PUT_ATT(fid,areaid,'valid_max',MAXVAL(area))
    iret = NF90_PUT_ATT(fid,areaid,'long_name',"Area of grid box")
    !
    iret = NF90_DEF_VAR(fid,"corners",NF90_REAL4,(/lonid,latid,nbcornersid,resid/), cornerid)
    iret = NF90_PUT_ATT(fid,cornerid,'standard_name',"gridcorners")
    iret = NF90_PUT_ATT(fid,cornerid,'units',"m*m")
    iret = NF90_PUT_ATT(fid,cornerid,'valid_min',MINVAL(corners))
    iret = NF90_PUT_ATT(fid,cornerid,'valid_max',MAXVAL(corners))
    iret = NF90_PUT_ATT(fid,cornerid,'long_name',"corners of grid boxes")
    !
    iret = NF90_DEF_VAR(fid,"landindex",NF90_INT, landdimid, landindexid)
    iret = NF90_PUT_ATT(fid,landindexid,'standard_name',"landindex")
    iret = NF90_PUT_ATT(fid,landindexid,'units',"-")
    iret = NF90_PUT_ATT(fid,landindexid,'valid_min',MINVAL(lindex))
    iret = NF90_PUT_ATT(fid,landindexid,'valid_max',MAXVAL(lindex))
    iret = NF90_PUT_ATT(fid,landindexid,'long_name',"Land index on global grid (FORTRAN convention)")
    !
    iret = NF90_DEF_VAR(fid,"contfrac",NF90_INT,(/landdimid/), contfracid)
    iret = NF90_PUT_ATT(fid,contfracid,'standard_name',"contfrac")
    iret = NF90_PUT_ATT(fid,contfracid,'units',"-")
    iret = NF90_PUT_ATT(fid,contfracid,'valid_min',MINVAL(contfrac))
    iret = NF90_PUT_ATT(fid,contfracid,'valid_max',MAXVAL(contfrac))
    iret = NF90_PUT_ATT(fid,contfracid,'long_name',"Fraction of continent in grid box")
    !
    ! Global attributes
    ! 
    iret = NF90_PUT_ATT(fid, NF90_GLOBAL,'calendar', calendar)
    !
    iret = NF90_ENDDEF (fid)
    IF (iret /= NF90_NOERR) THEN
       CALL ipslerr (3,'globgrd_writegrid',"Error ending definitions in file :",gridfilename, " ")
    ENDIF
    !
    ! Write variables
    !
    iret = NF90_PUT_VAR(fid, lonid, lon(:,1))
    iret = NF90_PUT_VAR(fid, latid, lat(1,:))
    iret = NF90_PUT_VAR(fid, maskid, mask)
    iret = NF90_PUT_VAR(fid, areaid, area)
    iret = NF90_PUT_VAR(fid, cornerid, corners)
    !
    iret = NF90_PUT_VAR(fid, landindexid,lindex)
    iret = NF90_PUT_VAR(fid, contfracid, contfrac)
    !
    ! Close file
    !
    iret = NF90_CLOSE(fid)
    IF (iret /= NF90_NOERR) THEN
       CALL ipslerr (3,'globgrd_writegrid',"Error closing the output file :",gridfilename, " ")
    ENDIF
    !
  END SUBROUTINE globgrd_writegrid
!!
!!  =============================================================================================================================
!! SUBROUTINE: globgrd_writevar
!!
!>\BRIEF      Writes the grid and a variable to a netCDF file to check if the grid was correctly interpreted by the module.
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================
!!
!
  SUBROUTINE globgrd_writevar(iim, jjm, lon, lat, nbpt, lalo, var, varname, filename)
    !
    ! Subroutine used to dump a compressed variable into a full lat/lon grid of a netCDF file
    !
    USE netcdf
    !
    ! ARGUMENTS
    !
    INTEGER(i_std), INTENT(in) :: iim, jjm, nbpt
    REAL(r_std), INTENT(in)    :: lon(iim,jjm), lat(iim,jjm)
    REAL(r_std), INTENT(in)    :: lalo(nbpt,2)
    REAL(r_std), INTENT(in)    :: var(nbpt)
    CHARACTER(LEN=*), INTENT(in) :: varname
    CHARACTER(LEN=*), INTENT(in) :: filename
    !
    ! LOCAL
    !
    INTEGER(i_std) :: iret, fid, i, ii, jj, nlonid, nlatid, varid
    INTEGER(i_std) :: ip1, im1, jp1, jm1, di, dj
    REAL(r_std) :: limlon, limlat
    INTEGER(i_std), DIMENSION(2) :: lolaid 
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: varfull, dist
    INTEGER(i_std), DIMENSION(2)             :: closest
    REAL(r_std), PARAMETER :: epsilon=0.001
    !
    !
    WRITE(*,*) "globgrd_writevar WRITE ", TRIM(varname), " into file ", TRIM(filename)
    !
    ALLOCATE(varfull(iim,jjm), dist(iim,jjm))
    varfull(:,:) = nf90_fill_real
    !
    ! Locate each point on the global grid
    !
    DO i=1,nbpt
       closest(1) = 99999999
       closest(2) = 99999999
       DO ii=1,iim
          DO jj=1,jjm
             ! Get neighbours
             ip1=ii+1
             im1=ii-1
             jp1=jj+1
             jm1=jj-1
             di=2
             dj=2
             ! Treat exceptions
             IF (ip1 > iim) THEN
                ip1=iim
                di=1
             ENDIF
             IF (im1 < 1) THEN
                im1=1
                di=1
             ENDIF
             IF ( jp1 > jjm) THEN
                jp1=jjm
                dj=1
             ENDIF
             IF ( jm1 < 1) THEN
                jm1=1
                dj=1
             ENDIF
             ! Calculate limits
             limlon=ABS(lon(ip1,jj)-lon(im1,jj))/di-epsilon
             limlat=ABS(lat(ii,jp1)-lat(ii,jm1))/dj-epsilon
             !
             IF ( ABS(lalo(i,1)-lat(ii,jj)) < limlat .AND. ABS(lalo(i,2)-lon(ii,jj)) < limlon ) THEN
                closest(1) = ii
                closest(2) = jj
             ENDIF
          ENDDO
       ENDDO
       IF ( closest(1) >  99999998 .OR. closest(2) >  99999998 ) THEN
          WRITE(*,*) "LALO closest : ", closest
          STOP "ERROR in globgrd_writevar"
       ELSE
          varfull(closest(1),closest(2)) = var(i)
       ENDIF
    ENDDO
    !
    ! Write the full variable into a NETCDF file
    !
    iret = NF90_CREATE(filename, NF90_WRITE, fid)
    IF (iret /= NF90_NOERR) THEN
       CALL ipslerr (3,'globgrd_writevar',"Error opening the output file :",filename, " ")
    ENDIF
    !
    ! Define dimensions
    !
    iret = NF90_DEF_DIM(fid,'Longitude',iim,lolaid(1))
    iret = NF90_DEF_DIM(fid,'Latitude',jjm,lolaid(2))
    !
    iret = NF90_DEF_VAR(fid,"Longitude",NF90_REAL4,lolaid,nlonid)
    iret = NF90_PUT_ATT(fid,nlonid,'standard_name',"longitude")
    iret = NF90_PUT_ATT(fid,nlonid,'units',"degrees_east")
    iret = NF90_PUT_ATT(fid,nlonid,'valid_min',MINVAL(lon))
    iret = NF90_PUT_ATT(fid,nlonid,'valid_max',MAXVAL(lon))
    iret = NF90_PUT_ATT(fid,nlonid,'long_name',"Longitude")
    !
    iret = NF90_DEF_VAR(fid,"Latitude",NF90_REAL4,lolaid,nlatid)
    iret = NF90_PUT_ATT(fid,nlatid,'standard_name',"latitude")
    iret = NF90_PUT_ATT(fid,nlatid,'units',"degrees_north")
    iret = NF90_PUT_ATT(fid,nlatid,'valid_min',MINVAL(lat))
    iret = NF90_PUT_ATT(fid,nlatid,'valid_max',MAXVAL(lat))
    iret = NF90_PUT_ATT(fid,nlatid,'long_name',"Latitude")
    !
    iret = NF90_DEF_VAR(fid,varname,NF90_REAL4,lolaid,varid)
    iret = NF90_PUT_ATT(fid,varid,'_FillValue',NF90_FILL_REAL)
    !
    iret = NF90_ENDDEF (fid)
    IF (iret /= NF90_NOERR) THEN
       CALL ipslerr (3,'globgrd_writevar',"Error ending definitions in file :",filename, " ")
    ENDIF
    !
    ! Write variables
    !
    iret = NF90_PUT_VAR(fid,nlonid,lon)
    iret = NF90_PUT_VAR(fid,nlatid,lat)
    iret = NF90_PUT_VAR(fid,varid,varfull)
    !
    ! Close file
    !
    iret = NF90_CLOSE(fid)
    IF (iret /= NF90_NOERR) THEN
       CALL ipslerr (3,'globgrd_writevar',"Error closing the output file :",filename, " ")
    ENDIF
    !
    DEALLOCATE(varfull,dist)
    !
    WRITE(*,*) "globgrd_writevar CLOSE file ", TRIM(filename)
    !
  END SUBROUTINE globgrd_writevar
!
END MODULE globgrd
