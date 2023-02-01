MODULE interpweight
!=================================================================================================================================
! MODULE       : interpweight
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       ! Specific module for the interpolation in ORCHIDEE
!!
!!\n DESCRIPTION: L. Fita, October 2016
!!   This module is used in the re-organization done in the section of the code
!! responsible of the interpolation of the morphology datat (PFT, albedo, LAI,
!! ...) from different files done by L. Fita in the mark of the HEAT project.
!! Previous version of the interpolation of these files, proceed by three 
!! basic steps common for all the files:
!!   1.: Modification of values from file
!!   2.: Creation of land/sea mask using values from file
!!   3.: Interpolate values
!!   This module servs to generalize all these steps and construct a common 
!! interface independently of the file. Each step is generalize and the 
!! different processes almost specific for each file in the previous form
!! are introduced as different methodologies. Thus now it provides:
!!   1. Modification of initial values: Removing that wrong values from the file
!|     `interpweight_modifying_input[1/2/3/4]D': subroutines to modify 1D, 2D, 3D
!!     or 4D input data from file, using different methods
!!   2. Masking input: Compute ’on fly’ the values of the land/sea mask taking
!!    the values in the file
!!     `interpweight_masking_input[1/2/3/4]D': subroutines to compute the mask
!!     using data from input file using different methods: nomask, mbelow, 
!!     mabove, msumrange, var
!!   3. Area weights: Get coincident areas from input file to the target 
!!    projection using aggregate_p
!!      No changes on it
!!   4. Interpolate: Use obtained areas and interpolate values at each grid
!!    point at the same format as it will be provided by XIOS
!!      `interpweight_provide_fractions[1/2/3/4]D': Perform interpolation for
!!      fraction/category data
!!      `interpweight_provide_interpolation[2/4]D': Perform interpolation for
!!      continuous data
!!      `variableusetypes': Variable to provide the values along the additional
!!      dimension to perform the interpolation. (e.g.: in case the 13 pfts
!!      are: 1,3,6,18,23,34,35,39,48,...)
!!   5. A new variable is added which explains the ‘availability’ of data to 
!!    perform the interpolation (0, 1). When it is negative it means that there
!!    was no data for that grid point
!!            availability = SUM(area_source)/Area_target_grid_point
!!   This module will disappear once interpolation will be done throuhguot XIOS
!! 
!! More details at:
!!   https://forge.ipsl.jussieu.fr/orchidee/wiki/DevelopmentActivities/ORCHIDEE-DYNAMICO
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S): None
!!
!! SVN         :
!! $HeadURL: $
!! $Date: $
!! $Revision: $
!! \n
!_================================================================================================================================

  USE ioipsl_para
  USE interpol_help
  USE grid

  IMPLICIT NONE

  PRIVATE
  PUBLIC interpweight_1D, interpweight_2D, interpweight_3D, interpweight_4D, &
    interpweight_2Dcont, interpweight_3Dcont, interpweight_4Dcont, &
    interpweight_2Dcont_tstep, interpweight_3Dcont_tstep, & ! With timsteps
    interpweight_get_varNdims_file, interpweight_RangeR,                         & 
    interpweight_get_var2dims_file, interpweight_get_var3dims_file, interpweight_get_var4dims_file,   &
    interpweight_Index2DLonLat, interpweight_ValVecR

  CONTAINS

!! ================================================================================================================================
!! SUBROUTINE   : interpweight_1D
!!
!>\BRIEF        ! Interpolate any input file to the grid of the model for a 1D-variable  
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): outvar1D, aoutvar
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE interpweight_1D(nbpt, Nvariabletypes, variabletypes, lalo, resolution, neighbours,       &
    contfrac, filename, varname, inlonname, inlatname, varmin, varmax, noneg, masktype,               &
    maskvalues, maskvarname, dim1, dim2, initime, typefrac,                                           &
    maxresollon, maxresollat, outvar1D, aoutvar)

    USE ioipsl
    USE constantes_var
    USE mod_orchidee_para

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !
    !! 0.1 Input variables
    !
    INTEGER(i_std), INTENT(in)                 :: nbpt                    !! Number of points for which the data 
                                                                          !!   needs to be interpolated
    INTEGER(i_std), INTENT(in)                 :: Nvariabletypes          !! Number of types of the variable
    REAL(r_std), DIMENSION(Nvariabletypes),    &                          !! Vector of values of the types
      INTENT(in)    :: variabletypes                                      !!   (-1, not used)
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in) :: lalo                    !! Vector of latitude and longitudes 
                                                                          !! (beware of the order = 1 : latitude, 
                                                                          !!   2 : longitude)
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in) :: resolution              !! The size in km of each grid-box in X and Y
    !
    INTEGER(i_std), DIMENSION(nbpt,NbNeighb), INTENT(in) :: neighbours    !! Vector of neighbours for each grid point 
                                                                          !!   1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW) 
                            
    REAL(r_std), DIMENSION(nbpt),  INTENT(in)  :: contfrac                !! Fraction of land in each grid box.
    CHARACTER(LEN=80), INTENT(in)              :: filename                !! name of the input map to read

    CHARACTER(LEN=80), INTENT(in)              :: varname                 !! name of the variable to interpolate
    CHARACTER(LEN=80), INTENT(in)              :: inlonname, inlatname    !! names of the longitude and latitude
    REAL(r_std), INTENT(inout)                 :: varmin, varmax          !! min/max values to use for the
                                                                          !!   renormalization
! Transformations to apply to the original read values from the file
    LOGICAL, INTENT(in)                        :: noneg                   !! no negative values
    CHARACTER(LEN=50), INTENT(in)              :: masktype                !! Type of masking
                                                                          !!   'nomask': no-mask is applied
                                                                          !!   'mbelow': take values below maskvalues(1)
                                                                          !!   'mabove': take values above maskvalues(1)
                                                                          !!   'msumrange': take values within 2 ranges;
                                                                          !!      maskvalues(2) <= SUM(vals(k)) <= maskvalues(1)
                                                                          !!      maskvalues(1) < SUM(vals(k)) <= maskvalues(3)
                                                                          !!        (normalized by maskvalues(3))
                                                                          !!   'var': mask values are taken from a 
                                                                          !!     variable (>0)
    REAL(r_std), DIMENSION(3), INTENT(in)      :: maskvalues              !! values to use to mask (according to 
                                                                          !!   `masktype') 
    CHARACTER(LEN=250), INTENT(in)             :: maskvarname             !! name of the variable to use to mask 
    INTEGER(i_std), INTENT(in)                 :: dim1, dim2              !! 4D complementary dimension sizes
    CHARACTER(LEN=50), INTENT(in)              :: typefrac                !! Type of fraction retrieval:
                                                                          !!   'default': standard

    INTEGER(i_std), INTENT(in)                 :: initime                 !! Initial time to take (for input data 
                                                                          !!   with time values, -1 not used)
    REAL(r_std), INTENT(in)                    :: maxresollon,maxresollat !! lon,lat maximum resolutions (in m)
                                                                          !!   (-un, not used)
    !! 0.2 Modified variables
    !
    !! 0.3 output variables
    !
    REAL(r_std), DIMENSION(nbpt), INTENT(out) :: outvar1D                 !! 1D output variable and re-dimensioned
    REAL(r_std), DIMENSION(nbpt), INTENT(out)  :: aoutvar                 !! availability of input data to
                                                                          !!   interpolate output variable 
                                                                          !!   (on the nbpt space)
    !
    !! 0.4 Local variables
    !
    INTEGER(i_std)                             :: iml, jml, lml, tml, fid, ib, ip, jp, it, jj, jv
    INTEGER                                    :: inNdims                 !! Number of dimensions of the input data
    INTEGER                                    :: inNLldims               !! Number of dimensions of the lon/lat input data
    REAL(r_std), ALLOCATABLE, DIMENSION(:)     :: lat_in, lon_in          !! latitude and
                                                                          !! longitude, extract from input map
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)   :: sub_area                !! the area of the fine grid in the 
                                                                          !!   model grid ???
                                                                          !! cf src_global/interpol_help.f90, 
                                                                          !!   line 377, called "areaoverlap"
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:)  :: sub_index             !! the indexes from the grid boxes from the 
                                                                          !!   data that go into the 
                                                                          !!   model's boxes  
                                                                          !! cf src_global/interpol_help.f90,
                                                                          !!   line 300, called "ip"

    REAL(r_std), ALLOCATABLE, DIMENSION(:)       :: invar1D               !! 1D input values
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)     :: invar2D               !! 2D input values
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)   :: invar3D               !! 3D input values
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:,:) :: invar4D               !! 4D input values

    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)   :: resol_in
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:)    :: mask
    REAL(r_std), ALLOCATABLE, DIMENSION(:)       :: maskvar
    INTEGER(i_std)                             :: nix, njx
    INTEGER(i_std)                             :: idi, nbvmax             !! nbvmax : number of maximum vegetation map
                                                                          !! points in the GCM grid ; idi : its counter
    CHARACTER(LEN=30)                          :: callsign                !! Allows to specify which variable is beeing
                                                                          !! treated
    LOGICAL                                    :: ok_interpol             !! optionnal return of aggregate_2d
    INTEGER                                    :: ALLOC_ERR 
    CHARACTER(LEN=50)                          :: S1
    REAL(r_std), ALLOCATABLE, DIMENSION(:)     :: variableusetypes
    CHARACTER(LEN=256)                         :: msg
    INTEGER, DIMENSION(:), ALLOCATABLE         :: indims

! Debug vars
    INTEGER                                    :: nb_coord, nb_var, nb_gat


!_ ================================================================================================================================
! Allocating input data

    IF (printlev >= 3) WRITE(numout,*)"In interpweight_1D filename: ",TRIM(filename),          &
      " varname:'" // TRIM(varname) // "'"
    
    IF (is_root_prc) THEN
       IF (printlev >=5) THEN
          WRITE(numout,*) "Entering 'interpweight_1D'. Debug mode."
          WRITE(numout,*)'(/," --> fliodmpf")'
          CALL fliodmpf (TRIM(filename))
          WRITE(numout,*)'(/," --> flioopfd")'
       ENDIF
       CALL flioopfd (TRIM(filename),fid,nb_dim=nb_coord,nb_var=nb_var,nb_gat=nb_gat)
       IF (printlev >=5) THEN
          WRITE (numout,'(" Number of coordinate        in the file : ",I2)') nb_coord
          WRITE (numout,'(" Number of variables         in the file : ",I2)') nb_var
          WRITE (numout,'(" Number of global attributes in the file : ",I2)') nb_gat
       ENDIF
    ENDIF

! Getting shape of input variable from input file
    inNdims = interpweight_get_varNdims_file(filename, varname)
    IF (is_root_prc) CALL flininfo(filename, iml, jml, lml, tml, fid)

    InputDims: SELECT CASE (inNdims)
      CASE (1)
        CALL bcast(iml)

        ALLOCATE(invar1D(iml), STAT=ALLOC_ERR)
        IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_1D',"Problem in allocation of variable '" //     &
          TRIM(varname) ,'','')
        IF (is_root_prc) CALL flinget(fid, TRIM(varname), iml, 0, 0, 0, 1, 1, invar1D)
        IF (printlev >= 5) THEN
          WRITE(numout,*)'  interpweight_1D input data dimensions _______'
          WRITE(numout,*)'    iml:',iml
        END IF

      CASE DEFAULT
        WRITE (S1,'(I3)') inNdims
        CALL ipslerr_p(3,'interpweight_1D',"Input number of dimensions: '" // TRIM(S1) // "' not ready !!", &
          '','')

    END SELECT InputDims
    CALL bcast(invar1D)

    IF (printlev >= 5) WRITE(numout,*)'  interpweight_1D getting variable ended!'

! Getting longitudes/latitudes from input file
    inNLldims = interpweight_get_varNdims_file(filename, inlonname)
    
    IF (inNLldims == 1) THEN
      ! Allocate memory for latitudes
      ALLOCATE(lon_in(iml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_1D','Problem in allocation of variable lon_in','','')

      ALLOCATE(lat_in(iml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_1D','Problem in allocation of variable lat_in','','')

      IF (is_root_prc) THEN
        CALL flinget(fid, TRIM(inlonname), iml, 0, 0, 0, 1, 1, lon_in)
        CALL flinget(fid, TRIM(inlatname), iml, 0, 0, 0, 1, 1, lat_in)
      ENDIF

    END IF
    CALL bcast(lon_in)
    CALL bcast(lat_in)

    IF (noneg) THEN
      IF (is_root_prc) THEN
        SELECT CASE (inNdims)
          CASE(1)
            CALL interpweight_modifying_input1D(iml, 0, 0, 0, 1, 1, zero, invar1D)
        END SELECT       
      END IF
    END IF
    IF (printlev >= 5) WRITE(numout,*)'  interpweight_1D modifying input ended!'

    IF (is_root_prc) CALL flinclo(fid)

    !
    ! Consider all points a priori
    !
    ALLOCATE(mask(iml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_1D','Problem in allocation of variable mask','','')
    mask = zero

    IF (TRIM(masktype) == 'var') THEN
      ALLOCATE(maskvar(iml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_1D','Problem in allocation of variable maskvar','','')
      maskvar = zero
! Reads the mask from the file
      IF (is_root_prc) THEN
! ATTTENTION!! One can not use flinget in this way. It always spect to have a REAL variable!
        CALL flinget(fid, TRIM(maskvarname), iml, 0, 0, 0, 1, 1, maskvar)
        WHERE (maskvar > zero)
          mask = un 
        END WHERE
      END IF
      CALL bcast(mask)
    ELSE
      SELECT CASE (inNdims)
        CASE(1)
          CALL interpweight_masking_input1D(iml, 0, 0, 0, 0, 0, masktype, maskvalues, invar1D,        &
            mask)
      END SELECT
    END IF
    IF (printlev >= 5) WRITE(numout,*)'  interpweight_1D masking input ended!'

    ! Computation of nbvmax
    ! In 'condveg_soilalb, slowproc_update' nbvmax=200
    !
    ! The number of maximum vegetation map points in the GCM grid is estimated.
    ! Some lmargin is taken.
    !
    IF (is_root_prc) THEN
      IF ( (maxresollon /= -un) .AND. (maxresollat /= -un) ) THEN
!       Assuming regular lon/lat projection for the input data but it is not caluclated assuming a maximum  &
!         lon,lat resolution given as input paramater
         nix=INT(MAXVAL(resolution(:,1))/maxresollon)+2
         njx=INT(MAXVAL(resolution(:,2))/maxresollat)+2
      ELSE
!       Assuming regular lon/lat projection for the input data
        CALL interpweight_calc_resolution_in(lon_in, lat_in, iml, jml, mincos, R_Earth, pi, resol_in)
        nix=INT(MAXVAL(resolution(:,1))/MAXVAL(resol_in(:,:,1)))+2
        njx=INT(MAXVAL(resolution(:,2))/MAXVAL(resol_in(:,:,2)))+2
      END IF
      nbvmax = MAX(nix*njx, 200)
    END IF  
    CALL bcast(nbvmax)
    CALL bcast(nix)
    CALL bcast(njx)

    !
    callsign = TRIM(varname) // ' map'
    !
    ok_interpol = .FALSE.
    DO WHILE ( .NOT. ok_interpol )
       IF (printlev>=3) WRITE(numout,*) "Projection arrays for ",callsign," : "
       IF (printlev>=3) WRITE(numout,*) "nbvmax = ",nbvmax, nix, njx


       ALLOCATE(sub_index(nbpt, nbvmax), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_1D','Problem in allocation of variable sub_index','','')

       sub_index(:,:)=zero

       ALLOCATE(sub_area(nbpt, nbvmax), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_1D','Problem in allocation of variable sub_area','','')

       sub_area(:,:)=zero
       IF (printlev>=3) WRITE(numout,*) 'Input data range LON:', MINVAL(lon_in), MAXVAL(lon_in)
       IF (printlev>=3) WRITE(numout,*) 'Input data range LAT:', MINVAL(lat_in), MAXVAL(lat_in)

       !
       CALL aggregate_p(nbpt, lalo, neighbours, resolution, contfrac,  &
            &                iml, lon_in, lat_in, maxresollon, maxresollat, callsign, &
            &                nbvmax, sub_index, sub_area, ok_interpol)
       
       !
       IF ( .NOT. ok_interpol ) THEN
          DEALLOCATE(sub_area)
          DEALLOCATE(sub_index)
          nbvmax = nbvmax * 2
       ENDIF
    ENDDO
    IF (printlev>=4) WRITE (numout,*) '  interpweight_1D: read/allocate OK'

! Getting variables thresholds
    ALLOCATE(variableusetypes(Nvariabletypes), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_1D','Problem in allocation of variableusetypes','','')

    IF (variabletypes(1) == -un) THEN
      variableusetypes = (/ (ip*1., ip=1,Nvariabletypes) /)
      IF (printlev>=4) WRITE(numout,*) 'interpweight_1D: Note: Using default ',&
           'equivalence between dimension in the file and dimension in the variable, for file ', filename
      varmin = 1._r_std
      varmax = Nvariabletypes*1._r_std
    ELSE
      variableusetypes = variabletypes    
    END IF

    outvar1D = zero

! Providing the fractions of each type for each interpolated grid point
    CALL interpweight_provide_fractions1D(typefrac, nbpt, Nvariabletypes, variableusetypes,           &
      iml, jml, lml, tml, nbvmax, zero, invar1D, sub_area, sub_index, varmax, varmin,                 &
      deux, lalo, outvar1D, aoutvar)
    IF (printlev>=5) WRITE (numout,*)'  interpweight_1D end of interp_provide_fractions1D'

    IF (printlev>=3) THEN
       IF (COUNT(aoutvar == -1) >= AINT(nbpt*0.4)) THEN
          WRITE(numout,*) '  interpweight_1D total number of points: ', nbpt
          WRITE(numout,*) '  interpweight_1D interpolated: ', COUNT(aoutvar /= -1),                 &
               ' non-interpolated: ', COUNT(aoutvar == -1)
          msg = "Less than 60% of points '" // TRIM(filename) // "' variable '" // TRIM(varname) //       & 
               "' have been succesfully interpolated !!"
          WRITE(numout,*) '  interpweight_1D: ' // TRIM(msg)
       ELSE
          msg = TRIM(filename) // "' variable: '" // TRIM(varname) // "' interpolation Done _______"
          WRITE(numout,*) "  interpweight_1D: '" // TRIM(msg)
          WRITE(numout,*) '  interpweight_1D total number of points: ', nbpt
          WRITE(numout,*) '  interpweight_1D interpolated: ', COUNT(aoutvar /= -1),               &
               ' non-interpolated: ', COUNT(aoutvar == -1)
       END IF
    END IF

! Scaling to 1 the quality variable
    DO ip=1, nbpt
      IF (aoutvar(ip) /= -1) aoutvar(ip) = aoutvar(ip)/(resolution(ip,1)*resolution(ip,2))/contfrac(ip)
    END DO

    IF (ALLOCATED(lon_in)) DEALLOCATE(lon_in)
    IF (ALLOCATED(lat_in)) DEALLOCATE(lat_in)
    IF (ALLOCATED(invar1D)) DEALLOCATE(invar1D)
    DEALLOCATE(mask)
    IF (ALLOCATED(maskvar)) DEALLOCATE(maskvar)
    DEALLOCATE(sub_area)
    DEALLOCATE(sub_index)

    RETURN

  END SUBROUTINE interpweight_1D

!! ================================================================================================================================
!! SUBROUTINE   : interpweight_2D
!!
!>\BRIEF        ! Interpolate any input file to the grid of the model for a 2D-variable  
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): outvar2D, aoutvar
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE interpweight_2D(nbpt, Nvariabletypes, variabletypes, lalo, resolution, neighbours,       &
    contfrac, filename, varname, inlonname, inlatname, varmin, varmax, noneg, masktype,               &
    maskvalues, maskvarname, dim1, dim2, initime, typefrac,                                           & 
    maxresollon, maxresollat, outvar2D, aoutvar)

    USE ioipsl
    USE constantes_var
    USE mod_orchidee_para

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !
    !! 0.1 Input variables
    !
    INTEGER(i_std), INTENT(in)                 :: nbpt                    !! Number of points for which the data 
                                                                          !!   needs to be interpolated
    INTEGER(i_std), INTENT(in)                 :: Nvariabletypes          !! Number of types of the variable
    REAL(r_std), DIMENSION(Nvariabletypes),    &                          !! Vector of values of the types
      INTENT(in)    :: variabletypes                                      !!   (-1, not used)
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in) :: lalo                    !! Vector of latitude and longitudes 
                                                                          !! (beware of the order = 1 : latitude, 
                                                                          !!   2 : longitude)
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in) :: resolution              !! The size in km of each grid-box in X and Y
    !
    INTEGER(i_std), DIMENSION(nbpt,NbNeighb), INTENT(in) :: neighbours    !! Vector of neighbours for each grid point 
                                                                          !!   1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW) 
                            
    REAL(r_std), DIMENSION(nbpt),  INTENT(in)  :: contfrac                !! Fraction of land in each grid box.
    CHARACTER(LEN=80), INTENT(in)              :: filename                !! name of the input map to read

    CHARACTER(LEN=80), INTENT(in)              :: varname                 !! name of the variable to interpolate
    CHARACTER(LEN=80), INTENT(in)              :: inlonname, inlatname    !! names of the longitude and latitude
    REAL(r_std), INTENT(inout)                 :: varmin, varmax          !! min/max values to use for the
                                                                          !!   renormalization
! Transformations to apply to the original read values from the file
    LOGICAL, INTENT(in)                        :: noneg                   !! no negative values
    CHARACTER(LEN=50), INTENT(in)              :: masktype                !! Type of masking
                                                                          !!   'nomask': no-mask is applied
                                                                          !!   'mbelow': take values below maskvalues(1)
                                                                          !!   'mabove': take values above maskvalues(1)
                                                                          !!   'msumrange': take values within 2 ranges;
                                                                          !!      maskvalues(2) <= SUM(vals(k)) <= maskvalues(1)
                                                                          !!      maskvalues(1) < SUM(vals(k)) <= maskvalues(3)
                                                                          !!        (normalized by maskvalues(3))
                                                                          !!   'var': mask values are taken from a 
                                                                          !!     variable (>0)
    REAL(r_std), DIMENSION(3), INTENT(in)      :: maskvalues              !! values to use to mask (according to 
                                                                          !!   `masktype') 
    CHARACTER(LEN=250), INTENT(in)             :: maskvarname             !! name of the variable to use to mask 
    CHARACTER(LEN=50), INTENT(in)              :: typefrac                !! Type of fraction retrieval:
                                                                          !!   'default': standard
    INTEGER(i_std), INTENT(in)                 :: dim1, dim2              !! 4D complementary dimension sizes
    INTEGER(i_std), INTENT(in)                 :: initime                 !! Initial time to take (for input data 
                                                                          !!   with time values, -1 not used)
    REAL(r_std), INTENT(in)                    :: maxresollon,maxresollat !! lon,lat maximum resolutions (in m)
    !! 0.2 Modified variables
    !
    !  0.3 OUTPUT
    !
!    REAL(r_std), INTENT(out)                   ::  laimap(nbpt,nvm,12)    !! lai read variable and re-dimensioned
    REAL(r_std), DIMENSION(nbpt,Nvariabletypes), INTENT(out) :: outvar2D  !! 2D output variable and re-dimensioned
    REAL(r_std), DIMENSION(nbpt), INTENT(out)  :: aoutvar                 !! availability of input data to
                                                                          !!   interpolate output variable 
                                                                          !!   (on the nbpt space)
    !
    !  0.4 LOCAL
    !
    INTEGER(i_std)                             :: iml, jml, lml, tml, fid, ib, ip, jp, it, jj, jv
    INTEGER                                    :: inNdims                 !! Number of dimensions of the input data
    INTEGER                                    :: inNLldims               !! Number of dimensions of the lon/lat input data
    REAL(r_std), ALLOCATABLE, DIMENSION(:)     :: lat_in, lon_in          !! latitude and
                                                                          !! longitude, extract from input map
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)   :: lat, lon                !! en 2D ???
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)   :: sub_area                !! the area of the fine grid in the 
                                                                          !!   model grid ???
                                                                          !! cf src_global/interpol_help.f90, 
                                                                          !!   line 377, called "areaoverlap"
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:,:)  :: sub_index           !! the indexes from the grid boxes from the 
                                                                          !!   data that go into the 
                                                                          !!   model's boxes  
                                                                          !! cf src_global/interpol_help.f90,
                                                                          !!   line 300, called "ip"

!    REAL(r_std), ALLOCATABLE, DIMENSION(:)       :: invar1D               !! 1D input values
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)     :: invar2D               !! 2D input values
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)   :: invar3D               !! 3D input values
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:,:) :: invar4D               !! 4D input values

    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)   :: resol_in
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:)  :: mask
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)     :: maskvar
    INTEGER(i_std)                             :: nix, njx
    INTEGER(i_std)                             :: idi, nbvmax             !! nbvmax : number of maximum vegetation map
                                                                          !! points in the GCM grid ; idi : its counter
    CHARACTER(LEN=30)                          :: callsign                !! Allows to specify which variable is beeing
                                                                          !! treated
    LOGICAL                                    :: ok_interpol             !! optionnal return of aggregate_2d
    !
    INTEGER                                    :: ALLOC_ERR 
    CHARACTER(LEN=50)                          :: S1
    REAL(r_std), ALLOCATABLE, DIMENSION(:)     :: variableusetypes
    CHARACTER(LEN=256)                         :: msg
    INTEGER, DIMENSION(:), ALLOCATABLE         :: indims

! Debug vars
    INTEGER                                    :: nb_coord, nb_var, nb_gat


!_ ================================================================================================================================
! Allocating input data

    IF (printlev >= 3) WRITE(numout,*)"In interpweight_2D filename: ",TRIM(filename),          &
      " varname:'" // TRIM(varname) // "'"
    !
    IF (is_root_prc) THEN
       IF (printlev >=5) THEN
          WRITE(numout,*) "Entering 'interpweight_2D'. Debug mode."
          WRITE(numout,*)'(/," --> fliodmpf")'
          CALL fliodmpf (TRIM(filename))
          WRITE(numout,*)'(/," --> flioopfd")'
       ENDIF
       CALL flioopfd (TRIM(filename),fid,nb_dim=nb_coord,nb_var=nb_var,nb_gat=nb_gat)
       IF (printlev >=5) THEN
          WRITE (numout,'(" Number of coordinate        in the file : ",I2)') nb_coord
          WRITE (numout,'(" Number of variables         in the file : ",I2)') nb_var
          WRITE (numout,'(" Number of global attributes in the file : ",I2)') nb_gat
       ENDIF
    ENDIF

! Getting shape of input variable from input file
    inNdims = interpweight_get_varNdims_file(filename, varname)
    IF (is_root_prc) CALL flininfo(filename, iml, jml, lml, tml, fid)
    
    ! Read variable (varname) from input file
    InputDims: SELECT CASE (inNdims)
      CASE (2)
        CALL bcast(iml)
        CALL bcast(jml)

        ALLOCATE(invar2D(iml,jml), STAT=ALLOC_ERR)
        IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2D',"Problem in allocation of variable '" //     &
          TRIM(varname) ,'','')
        IF (is_root_prc) CALL flinget(fid, TRIM(varname), iml, jml, 0, 0, 1, 1, invar2D)
        IF (printlev >= 5) THEN
          WRITE(numout,*)'  interpweight_2D input data dimensions: '
          WRITE(numout,*)'    iml, jml:',iml, jml
        END IF

      CASE (3)
        ALLOCATE(indims(3), STAT=ALLOC_ERR)
        IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2D','Problem in allocation of indims','','')
        indims = interpweight_get_var3dims_file(filename, varname)

        lml = indims(3)
        CALL bcast(iml)
        CALL bcast(jml)
        CALL bcast(lml)

        ALLOCATE(invar3D(iml,jml, lml), STAT=ALLOC_ERR)
        IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2D',"Problem in allocation of variable '" //     &
          TRIM(varname) ,'','')
        IF (is_root_prc) CALL flinget(fid, TRIM(varname), iml, jml, lml, 0, 1, 1, invar3D)
        IF (printlev >= 5) THEN
          WRITE(numout,*)'  interpweight_2D input data dimensions:'
          WRITE(numout,*)'    iml, jml, lml:',iml, jml, lml
        END IF

      CASE (4)
        ALLOCATE(indims(4), STAT=ALLOC_ERR)
        IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2D','Problem in allocation of indims','','')
        indims = interpweight_get_var4dims_file(filename, varname)
        lml = indims(3)
        CALL bcast(iml)
        CALL bcast(jml)
        CALL bcast(lml)
!       Taken from `slowproc_update'
        IF (initime /= -1) THEN
          inNdims = 3
          IF (printlev >= 3)  WRITE(numout,*) &
               'interpweight_2D: taking into account the initial reading time:',initime
          IF (initime <= 0) tml = 0

          CALL bcast(tml)
          ALLOCATE(invar3D(iml,jml,lml), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2D',"Problem in allocation of variable '" //   &
            TRIM(varname) ,'','')
          IF ( initime > 0 ) THEN
            IF (is_root_prc) THEN
              IF (initime <= tml) THEN
                 CALL flinget(fid, TRIM(varname), iml, jml, lml, 1, initime, 1, invar3D)
              ELSE
                 CALL flinget(fid, TRIM(varname), iml, jml, lml, 1, tml, 1, invar3D)
              ENDIF
            ENDIF
          ELSE
            IF (is_root_prc) THEN
              CALL flinget(fid, TRIM(varname), iml, jml, lml, 1, 1, 1, invar3D)
            ENDIF
          ENDIF
        ELSE
          CALL bcast(tml)
          ALLOCATE(invar4D(iml,jml,lml,tml), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2D',"Problem in allocation of variable '" //   &
            TRIM(varname) ,'','')
          IF (is_root_prc) CALL flinget(fid, TRIM(varname), iml, jml, lml, tml, 1, tml, invar4D)
        END IF

        IF (printlev >= 5) THEN
          WRITE(numout,*)'  interpweight_2D input data dimensions'
          WRITE(numout,*)'    iml, jml, lml, tml :',iml, jml, lml, tml
        END IF

      CASE DEFAULT
        WRITE (S1,'(I3)') inNdims
        CALL ipslerr_p(3,'interpweight_2D',"Input number of dimensions: '" // TRIM(S1) // "' not ready !!", &
          '','')

    END SELECT InputDims
    CALL bcast(invar2D)

    IF (printlev >= 5) WRITE(numout,*)'  interpweight_2D getting variable ended!'

! Getting longitudes/latitudes from input file
    inNLldims = interpweight_get_varNdims_file(filename, inlonname)
    
    IF (inNLldims == 1) THEN
      ! Allocate memory for latitudes
      ALLOCATE(lon_in(iml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2D','Problem in allocation of variable lon_in','','')

      ALLOCATE(lat_in(jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2D','Problem in allocation of variable lat_in','','')

      IF (is_root_prc) THEN
        CALL flinget(fid, TRIM(inlonname), iml, 0, 0, 0, 1, 1, lon_in)
        CALL flinget(fid, TRIM(inlatname), jml, 0, 0, 0, 1, 1, lat_in)
      ENDIF
      CALL bcast(lon_in)
      CALL bcast(lat_in)

      ALLOCATE(lon(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2D','Problem in allocation of variable lon','','')
      ALLOCATE(lat(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2D','Problem in allocation of variable lat','','')

      DO ip=1,iml
         lat(ip,:) = lat_in(:)
      ENDDO
      DO jp=1,jml
         lon(:,jp) = lon_in(:)
      ENDDO
    ELSE
      ! Allocate memory for longitude
      ALLOCATE(lon(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2D','Pb in allocation for lon','','')
      ! Allocate memory for latitudes
      ALLOCATE(lat(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2D','Pb in allocation for lat','','')

      IF (is_root_prc) THEN
        CALL flinget(fid, TRIM(inlonname), iml, jml, 0, 0, 1, 1, lon)
        CALL flinget(fid, TRIM(inlatname), iml, jml, 0, 0, 1, 1, lat)
      ENDIF
    END IF
    CALL bcast(lon)
    CALL bcast(lat)

    IF (is_root_prc) THEN
      ALLOCATE(resol_in(iml,jml,2), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2D','Problem in allocation of variable resol_in','','')
    END IF

    IF (noneg) THEN
      IF (is_root_prc) THEN
        SELECT CASE (inNdims)
          CASE(2)
            CALL interpweight_modifying_input2D(iml, jml, 0, 0, 1, 1, zero, invar2D)
          CASE(3)
            CALL interpweight_modifying_input3D(iml, jml, lml, 0, 1, 1, zero, invar3D)
          CASE(4)
            CALL interpweight_modifying_input4D(iml, jml, lml, tml, 1, tml, zero, invar4D)
        END SELECT       
      END IF
    END IF
    IF (printlev >= 5) WRITE(numout,*)'  interpweight_2D modifying input ended!'

    !
    ! Consider all points a priori
    !
    ALLOCATE(mask(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2D','Problem in allocation of variable mask','','')
    mask = zero

    IF (TRIM(masktype) == 'var') THEN
      ALLOCATE(maskvar(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2D','Problem in allocation of variable maskvar','','')
      maskvar(:,:) = 0
! Reads the mask from the file
      IF (is_root_prc) THEN
! ATTTENTION!! One can not use flinget in this way. It always spect to have a REAL variable!
        CALL flinget(fid, TRIM(maskvarname), iml, jml, 0, 0, 1, 1, maskvar)
        WHERE (maskvar > zero)
          mask = un 
        END WHERE
      END IF
      CALL bcast(mask)
    ELSE
      ! Create land-sea mask based on values from the variable according to maskvalues and masktype
      SELECT CASE (inNdims)
        CASE(2)
          CALL interpweight_masking_input2D(iml, jml, 0, 0, 0, 0, masktype, maskvalues, invar2D,      &
            mask)
        CASE(3)
          CALL interpweight_masking_input3D(iml, jml, 0, 0, lml, 0, masktype, maskvalues, invar3D,    &
            mask)
        CASE(4)
          CALL interpweight_masking_input4D(iml, jml, 0, 0, lml, tml, masktype, maskvalues, invar4D,  &
            mask)
      END SELECT
    END IF
    IF (printlev >= 5) WRITE(numout,*)'  interpweight_2D masking input ended!'

    IF (is_root_prc) CALL flinclo(fid)

! Assuming regular lon/lat projection for the input data
    IF (is_root_prc) THEN
      CALL interpweight_calc_resolution_in(lon, lat, iml, jml, mincos, R_Earth, pi, resol_in)
    END IF

    ! Computation of nbvmax
    ! In 'condveg_soilalb, slowproc_update' nbvmax=200
    !
    ! The number of maximum vegetation map points in the GCM grid is estimated.
    ! Some lmargin is taken.
    !
    IF (is_root_prc) THEN
      IF ( (maxresollon /= -un) .AND. (maxresollat /= -un) ) THEN
!       Assuming regular lon/lat projection for the input data but it is not caluclated assuming a maximum  &
!         lon,lat resolution given as input paramater
         nix=INT(MAXVAL(resolution(:,1))/maxresollon)+2
         njx=INT(MAXVAL(resolution(:,2))/maxresollat)+2
      ELSE
!       Assuming regular lon/lat projection for the input data
        CALL interpweight_calc_resolution_in(lon, lat, iml, jml, mincos, R_Earth, pi, resol_in)
        nix=INT(MAXVAL(resolution(:,1))/MAXVAL(resol_in(:,:,1)))+2
        njx=INT(MAXVAL(resolution(:,2))/MAXVAL(resol_in(:,:,2)))+2
      END IF
      nbvmax = MAX(nix*njx, 200)
    END IF  
    CALL bcast(nbvmax)
    CALL bcast(nix)
    CALL bcast(njx)

    !
    callsign = TRIM(varname) // ' map'
    !
    ok_interpol = .FALSE.
    DO WHILE ( .NOT. ok_interpol )
       IF (printlev>=3) WRITE(numout,*) "Projection arrays for ",callsign," : "
       IF (printlev>=3) WRITE(numout,*) "nbvmax = ",nbvmax, nix, njx


       ALLOCATE(sub_index(nbpt, nbvmax, 2), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2D','Problem in allocation of variable sub_index','','')

       sub_index(:,:,:)=zero

       ALLOCATE(sub_area(nbpt, nbvmax), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2D','Problem in allocation of variable sub_area','','')

       sub_area(:,:)=zero
       IF (printlev>=3) WRITE(numout,*) 'Input data range LON:', MINVAL(lon), MAXVAL(lon)
       IF (printlev>=3) WRITE(numout,*) 'Input data range LAT:', MINVAL(lat), MAXVAL(lat)

       !
       CALL aggregate_p(nbpt, lalo, neighbours, resolution, contfrac, &
            &                iml, jml, lon, lat, mask, callsign, &
            &                nbvmax, sub_index, sub_area, ok_interpol)
       
       !
       IF ( .NOT. ok_interpol ) THEN
          DEALLOCATE(sub_area)
          DEALLOCATE(sub_index)
          nbvmax = nbvmax * 2
       ENDIF
    ENDDO
    IF (printlev>=4) WRITE (numout,*) '  interpweight_2D: read/allocate OK'

! Getting variables thresholds
    ALLOCATE(variableusetypes(Nvariabletypes), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2D','Problem in allocation of variableusetypes','','')
    

    IF (variabletypes(1) == -un) THEN
      variableusetypes = (/ (ip*1., ip=1,Nvariabletypes) /)
      IF (printlev>=4) THEN
         WRITE(numout,*)'  interpweight_2D: Note: Using default equivalence between dimension ',&
              'in the file and dimension in the variable, for file ', filename
      END IF
      varmin = 1._r_std
      varmax = Nvariabletypes*1._r_std
    ELSE
      variableusetypes = variabletypes    
    END IF

    outvar2D = zero

! Providing the fractions of each type for each interpolated grid point
    CALL interpweight_provide_fractions2D(typefrac, nbpt, Nvariabletypes, variableusetypes,           &
      iml, jml, lml, tml, nbvmax, zero, invar2D, sub_area, sub_index, varmax, varmin,                 &
      deux, lalo, outvar2D, aoutvar)
    IF (printlev>=5) WRITE (numout,*)'  interpweight_2D end of provide_fractions2D'

    IF (printlev>=3) THEN
       IF (COUNT(aoutvar == -1) >= AINT(nbpt*0.4)) THEN
          WRITE(numout,*) '  interpweight_2D total number of points: ', nbpt
          WRITE(numout,*) '  interpweight_2D interpolated: ', COUNT(aoutvar /= -1),                   &
               ' non-interpolated: ', COUNT(aoutvar == -1)
          msg = "Less than 60% of points '" // TRIM(filename) // "' variable '" // TRIM(varname) //       & 
               "' have been succesfully interpolated !!"
          WRITE(numout,*) '  interpweight_2D: ' // TRIM(msg)
       ELSE
          msg = TRIM(filename) // "' variable: '" // TRIM(varname) // "' interpolation Done _______"
          WRITE(numout,*) "interpweight_2D: '" // TRIM(msg)
          WRITE(numout,*) '  Total number of points: ', nbpt
          WRITE(numout,*) '  Interpolated: ', COUNT(aoutvar /= -1),                 &
               ' Non-interpolated: ', COUNT(aoutvar == -1)
       END IF
    END IF

! Scaling to 1 the quality variable
    DO ip=1, nbpt
      IF (aoutvar(ip) /= -1) aoutvar(ip) = aoutvar(ip)/(resolution(ip,1)*resolution(ip,2))/contfrac(ip)
    END DO

    IF (ALLOCATED(lon_in)) DEALLOCATE(lon_in)
    IF (ALLOCATED(lat_in)) DEALLOCATE(lat_in)
    DEALLOCATE(lon)
    DEALLOCATE(lat)
    IF (ALLOCATED(invar2D)) DEALLOCATE(invar2D)
    IF (ALLOCATED(invar3D)) DEALLOCATE(invar3D)
    IF (ALLOCATED(invar4D)) DEALLOCATE(invar4D)
    DEALLOCATE(mask)
    IF (ALLOCATED(maskvar)) DEALLOCATE(maskvar)
    DEALLOCATE(sub_area)
    DEALLOCATE(sub_index)
    IF (ALLOCATED(resol_in)) DEALLOCATE(resol_in)

    RETURN

  END SUBROUTINE interpweight_2D

!! ================================================================================================================================
!! SUBROUTINE   : interpweight_3D
!!
!>\BRIEF        ! Interpolate any input file to the grid of the model for a 3D-variable  
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): outvar3D, aoutvar
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE interpweight_3D(nbpt, Nvariabletypes, variabletypes, lalo, resolution, neighbours,       &
    contfrac, filename, varname, inlonname, inlatname, varmin, varmax, noneg, masktype,               &
    maskvalues, maskvarname, dim1, dim2, initime, typefrac,                                           &
    maxresollon, maxresollat, outvar3D, aoutvar)

    USE ioipsl
    USE constantes_var
    USE mod_orchidee_para

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !
    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)                 :: nbpt                    !! Number of points for which the data 
                                                                          !!   needs to be interpolated
    INTEGER(i_std), INTENT(in)                 :: Nvariabletypes          !! Number of types of the variable
    REAL(r_std), DIMENSION(Nvariabletypes),    &                          !! Vector of values of the types
      INTENT(in)    :: variabletypes                                      !!   (-1, not used)
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in) :: lalo                    !! Vector of latitude and longitudes 
                                                                          !! (beware of the order = 1 : latitude, 
                                                                          !!   2 : longitude)
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in) :: resolution              !! The size in km of each grid-box in X and Y
    INTEGER(i_std), DIMENSION(nbpt,NbNeighb), INTENT(in) :: neighbours    !! Vector of neighbours for each grid point 
                                                                          !!   1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW) 
                            
    REAL(r_std), DIMENSION(nbpt), INTENT(in)   :: contfrac                !! Fraction of land in each grid box.
    CHARACTER(LEN=80), INTENT(in)              :: filename                !! name of the input map to read

    CHARACTER(LEN=80), INTENT(in)              :: varname                 !! name of the variable to interpolate
    CHARACTER(LEN=80), INTENT(in)              :: inlonname, inlatname    !! names of the longitude and latitude
                                                                          !!   in the input file
! Transformations to apply to the original read values from the file
    LOGICAL, INTENT(in)                        :: noneg                   !! no negative values
    CHARACTER(LEN=50), INTENT(in)              :: masktype                !! Type of masking
                                                                          !!   'nomask': no-mask is applied
                                                                          !!   'mbelow': take values below maskvalues(1)
                                                                          !!   'mabove': take values above maskvalues(1)
                                                                          !!   'msumrange': take values within 2 ranges;
                                                                          !!      maskvalues(2) <= SUM(vals(k)) <= maskvalues(1)
                                                                          !!      maskvalues(1) < SUM(vals(k)) <= maskvalues(3)
                                                                          !!        (normalized by maskvalues(3))
                                                                          !!   'var': mask values are taken from a 
                                                                          !!     variable (>0)
    REAL(r_std), DIMENSION(3), INTENT(in)      :: maskvalues              !! values to use to mask (according to 
                                                                          !!   `masktype') 
    CHARACTER(LEN=250), INTENT(in)             :: maskvarname             !! name of the variable to use to mask 
    INTEGER, INTENT(in)                        :: dim1, dim2              !! 3/4D dimensions of the output variable
    REAL(r_std), DIMENSION(dim1), INTENT(inout) :: varmin, varmax          !! min/max values to use for the renormalization
    CHARACTER(LEN=50), INTENT(in)              :: typefrac                !! Type of fraction retrieval:
                                                                          !!   'default': standard
    INTEGER(i_std), INTENT(in)                 :: initime                 !! Initial time to take (for input data 
                                                                          !!   with time values, -1 not used)
    REAL(r_std), INTENT(in)                    :: maxresollon,maxresollat !! lon,lat maximum resolutions (in m)
                                                                          !!   (-un, not used)
    !! 0.2 Modified variables
    !
    !  0.3 OUTPUT
    !
    REAL(r_std), DIMENSION(nbpt,dim1), INTENT(out) :: outvar3D            !! 3D output variable and re-dimensioned
! aovar: availability of input data to interpolate output variable (on the nbpt space)
    REAL(r_std), DIMENSION(nbpt), INTENT(out)  :: aoutvar                 !! availability of input data to
                                                                          !!   interpolate output variable 
                                                                          !!   (on the nbpt space)
    !
    !  0.4 Local variables
    !
    INTEGER(i_std)                             :: iml, jml, lml, tml, fid, ib, ip, jp, zp, it, jj, jv
    INTEGER                                    :: inNdims                 !! Number of dimensions of the input data
    INTEGER                                    :: inNLldims               !! Number of dimensions of the lon/lat input data
    REAL(r_std), ALLOCATABLE, DIMENSION(:)     :: lat_in, lon_in          !! latitude and
                                                                          !! longitude, extract from input map
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)   :: lat, lon                !! en 2D ???
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)   :: sub_area                !! the area of the fine grid in the 
                                                                          !!   model grid ???
                                                                          !! cf src_global/interpol_help.f90, 
                                                                          !!   line 377, called "areaoverlap"
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:,:)  :: sub_index           !! the indexes from the grid boxes from the 
                                                                          !!   data that go into the 
                                                                          !!   model's boxes  
                                                                          !! cf src_global/interpol_help.f90,
                                                                          !!   line 300, called "ip"

!    REAL(r_std), ALLOCATABLE, DIMENSION(:)       :: invar1D               !! 1D input values
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)     :: invar2D               !! 2D input values
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)   :: invar3D               !! 3D input values
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:,:) :: invar4D               !! 4D input values

    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)   :: resol_in
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:)  :: mask
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)     :: maskvar
    INTEGER(i_std)                             :: nix, njx
    INTEGER(i_std)                             :: idi, nbvmax             !! nbvmax : number of maximum vegetation map
                                                                          !! points in the GCM grid ; idi : its counter
    CHARACTER(LEN=30)                          :: callsign                !! Allows to specify which variable is beeing
                                                                          !! treated
    LOGICAL                                    :: ok_interpol             !! optionnal return of aggregate_2d
    INTEGER                                    :: ALLOC_ERR 

    CHARACTER(LEN=50)                          :: S1
    REAL(r_std), ALLOCATABLE, DIMENSION(:)     :: variableusetypes
    CHARACTER(LEN=256)                         :: msg

! Debug vars
    INTEGER                                    :: nb_coord, nb_var, nb_gat
    CHARACTER(LEN=5)                           :: iS
    CHARACTER(LEN=50)                          :: suminS

    INTEGER, ALLOCATABLE, DIMENSION(:)         :: indims


!_ ================================================================================================================================
    IF (printlev >= 3) WRITE(numout, *)"In interpweight_3D filename: ",TRIM(filename),         &
      " varname:'" // TRIM(varname) // "'"
! Allocating input data
    IF (is_root_prc) THEN
       IF (printlev >=5) THEN
          WRITE(numout,*) "Entering 'interpweight_3D'. Debug mode."
          WRITE(numout,*)'(/," --> fliodmpf")'
          CALL fliodmpf (TRIM(filename))
          WRITE(numout,*)'(/," --> flioopfd")'
       ENDIF
       CALL flioopfd (TRIM(filename),fid,nb_dim=nb_coord,nb_var=nb_var,nb_gat=nb_gat)
       IF (printlev >=5) THEN
          WRITE (numout,'(" Number of coordinate        in the file : ",I2)') nb_coord
          WRITE (numout,'(" Number of variables         in the file : ",I2)') nb_var
          WRITE (numout,'(" Number of global attributes in the file : ",I2)') nb_gat
       ENDIF
    ENDIF

! Getting shape of input variable from input file
    inNdims = interpweight_get_varNdims_file(filename, varname)
    IF (is_root_prc) CALL flininfo(filename, iml, jml, lml, tml, fid)

    InputDims: SELECT CASE (inNdims)
      CASE (2)
        CALL bcast(iml)
        CALL bcast(jml)

        ALLOCATE(invar2D(iml,jml), STAT=ALLOC_ERR)
        IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_3D',"Problem in allocation of variable '" //     &
          TRIM(varname) ,'','')
        IF (is_root_prc) CALL flinget(fid, TRIM(varname), iml, jml, 0, 0, 1, 1, invar2D)
        IF (printlev >= 5) THEN
          WRITE(numout,*)'  interpweight_3D input data dimensions'
          WRITE(numout,*)'    iml, jml:',iml, jml
        END IF

      CASE (3)
        ALLOCATE(indims(3), STAT=ALLOC_ERR)
        IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_3D','Problem in allocation of indims','','')
        indims = interpweight_get_var3dims_file(filename, varname)

        lml = indims(3)
        CALL bcast(iml)
        CALL bcast(jml)
        CALL bcast(lml)

        ALLOCATE(invar3D(iml,jml, lml), STAT=ALLOC_ERR)
        IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_3D',"Problem in allocation of variable '" //     &
          TRIM(varname) ,'','') 
        IF (is_root_prc) CALL flinget(fid, TRIM(varname), iml, jml, lml, 0, 2, 1, invar3D)
        IF (printlev >= 5) THEN
          WRITE(numout,*)'  interpweight_3D input data dimensions'
          WRITE(numout,*)'    iml, jml, lml :',iml, jml, lml
        END IF

      CASE (4)
        ALLOCATE(indims(4), STAT=ALLOC_ERR)
        IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_3D','Problem in allocation of indims','','')
        
        indims = interpweight_get_var4dims_file(filename, varname)
        lml = indims(3)
        CALL bcast(iml)
        CALL bcast(jml)
        CALL bcast(lml)
!       Taken from `slowproc_update'
        IF (initime /= -1) THEN
          inNdims = 3
          IF (printlev >= 3) WRITE (numout,*) &
               'interpweight_3D: taking into account the initial reading time:',initime
          IF (initime <= 0) tml = 0

          CALL bcast(tml)
          ALLOCATE(invar3D(iml,jml,lml), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_3D','Problem in allocation of variable invar3D','','')
          IF ( initime > 0 ) THEN
            IF (is_root_prc) THEN
              IF (initime <= tml) THEN
                CALL flinget(fid, TRIM(varname), iml, jml, lml, 1, initime, 1, invar3D)
              ELSE
                CALL flinget(fid, TRIM(varname), iml, jml, lml, 1, tml, 1, invar3D)
              ENDIF
            ENDIF
          ELSE
            IF (is_root_prc) THEN
              CALL flinget(fid, TRIM(varname), iml, jml, lml, 1, 1, 1, invar3D)
            ENDIF
          ENDIF
        ELSE
          CALL bcast(tml)
          ALLOCATE(invar4D(iml,jml,lml,tml), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_3D','Problem in allocation of variable invar4D','','')
          IF (is_root_prc) CALL flinget(fid, TRIM(varname), iml, jml, lml, tml, 1, tml, invar4D)
        END IF

        IF (printlev >= 5) THEN
          WRITE(numout,*)'  interpweight_3D input data dimensions'
          WRITE(numout,*)'    iml, jml, lml, tml:',iml, jml, lml, tml
        END IF

      CASE DEFAULT
        WRITE (S1,'(I3)') inNdims
        CALL ipslerr_p(3,'interpweight_3D',"Input number of dimensions: '" // TRIM(S1) // "' not ready !!", &
          '','')

    END SELECT InputDims
    CALL bcast(invar3D)
  
    IF (printlev >= 5) WRITE(numout,*)'  interpweight_3D getting variable ended!'
! Nan values are not allowed
    DO ip=1,SIZE(invar3d, DIM=1)
      DO jp=1,SIZE(invar3d, DIM=2)
        DO zp=1,SIZE(invar3d, DIM=3)
          IF ( invar3D(ip, jp, zp) /= invar3d(ip, jp, zp) ) THEN
            CALL ipslerr_p(3, 'interpweight_3D', 'Nan values are not allowed in', 'netcdf orchidee files', 'missing_values have to be defined as 1.e+20f')
          ENDIF
        ENDDO
      ENDDO
    ENDDO

! Test lecture
    DO ip=1,lml
      IF ( SUM(invar3D(:,:,ip)) < 0.00001 ) THEN
        WRITE(iS, '(i5)')ip 
        WRITE(suminS, '(f20.10)')SUM(invar3d(:,:,ip)) 
        msg = 'Wrong lecture of data for type ' // TRIM(iS) // ' all values are zero !!'
        CALL ipslerr_p(3,'interpweight_3D',TRIM(msg),'','')
      END IF
    END DO

! Getting longitudes/latitudes from input file
    inNLldims = interpweight_get_varNdims_file(filename, inlonname)
    IF (inNLldims == 1) THEN
      ! Allocate memory for latitudes
      ALLOCATE(lon_in(iml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_3D','Problem in allocation of variable lon_in','','')

      ALLOCATE(lat_in(jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_3D','Problem in allocation of variable lat_in','','')

      IF (is_root_prc) THEN
        CALL flininspect(fid)
        CALL flinget(fid, TRIM(inlonname), iml, 0, 0, 0, 1, 1, lon_in)
        CALL flinget(fid, TRIM(inlatname), jml, 0, 0, 0, 1, 1, lat_in)
      ENDIF
      CALL bcast(lon_in)
      CALL bcast(lat_in)

      ALLOCATE(lon(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_3D','Problem in allocation of variable lon','','')
      ALLOCATE(lat(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_3D','Problem in allocation of variable lat','','')

      DO ip=1,iml
         lat(ip,:) = lat_in(:)
      ENDDO
      DO jp=1,jml
         lon(:,jp) = lon_in(:)
      ENDDO
    ELSE
      ! Allocate memory for longitude
      ALLOCATE(lon(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_3D','Pb in allocation for lon','','')
      ! Allocate memory for latitudes
      ALLOCATE(lat(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_3D','Pb in allocation for lat','','')

      IF (is_root_prc) THEN
        CALL flinget(fid, TRIM(inlonname), iml, jml, 0, 0, 1, 1, lon)
        CALL flinget(fid, TRIM(inlatname), iml, jml, 0, 0, 1, 1, lat)
      ENDIF
    END IF
    CALL bcast(lon)
    CALL bcast(lat)

    ALLOCATE(resol_in(iml,jml,2), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_3D','Problem in allocation of variable resol_in','','')

    IF (noneg) THEN
      IF (is_root_prc) THEN
        SELECT CASE (inNdims)
          CASE(2)
            CALL interpweight_modifying_input2D(iml, jml, 0, 0, 1, 1, zero, invar2D)
          CASE(3)
            CALL interpweight_modifying_input3D(iml, jml, lml, 0, 1, 1, zero, invar3D)
          CASE(4)
            CALL interpweight_modifying_input4D(iml, jml, lml, tml, 1, tml, zero, invar4D)
        END SELECT       
      END IF
    END IF
    IF (printlev >= 5) WRITE(numout,*)'  interpweight_3D modifying input ended!'

    ALLOCATE(mask(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_3D','Problem in allocation of variable mask','','')
    mask(:,:) = 0

    IF (TRIM(masktype) == 'var') THEN
      ALLOCATE(maskvar(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_3D','Problem in allocation of variable maskvar','','')
      maskvar(:,:) = 0.
! Reads the mask from the file
      IF (is_root_prc) THEN
! ATTTENTION!! One can not use flinget in this way. It always spect to have a REAL variable!
        CALL flinget(fid, TRIM(maskvarname), iml, jml, 0, 0, 1, 1, maskvar)
        WHERE (maskvar > zero)
          mask = un 
        END WHERE
      END IF
      CALL bcast(mask)
    ELSE
      CALL interpweight_masking_input3D(iml, jml, 0, 0, lml, 0, masktype, maskvalues, invar3D,    &
            mask)
    END IF
    IF (printlev >= 5) WRITE(numout,*)'  interpweight_3D masking input ended!'

    IF (is_root_prc) CALL flinclo(fid)

    ! Computation of nbvmax
    ! In 'condveg_soilalb, slowproc_update' nbvmax=200
    !
    ! The number of maximum vegetation map points in the GCM grid is estimated.
    ! Some lmargin is taken.
    !
    IF (is_root_prc) THEN
      IF ( (maxresollon /= -un) .AND. (maxresollat /= -un) ) THEN
!       Assuming regular lon/lat projection for the input data but it is not caluclated assuming a maximum  &
!         lon,lat resolution given as input paramater
         nix=INT(MAXVAL(resolution(:,1))/maxresollon)+2
         njx=INT(MAXVAL(resolution(:,2))/maxresollat)+2
      ELSE
!       Assuming regular lon/lat projection for the input data
        CALL interpweight_calc_resolution_in(lon, lat, iml, jml, mincos, R_Earth, pi, resol_in)
        nix=INT(MAXVAL(resolution(:,1))/MAXVAL(resol_in(:,:,1)))+2
        njx=INT(MAXVAL(resolution(:,2))/MAXVAL(resol_in(:,:,2)))+2
      END IF
      nbvmax = MAX(nix*njx, 200)
    END IF  
    CALL bcast(nbvmax)
    CALL bcast(nix)
    CALL bcast(njx)

    !
    callsign = TRIM(varname) // ' map'
    !
    ok_interpol = .FALSE.
    DO WHILE ( .NOT. ok_interpol )
       IF (printlev>=3) WRITE(numout,*) "Projection arrays for ",callsign," : "
       IF (printlev>=3) WRITE(numout,*) "nbvmax = ",nbvmax, nix, njx

       ALLOCATE(sub_index(nbpt, nbvmax, 2), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_3D','Problem in allocation of variable sub_index','','')

       sub_index(:,:,:)=zero

       ALLOCATE(sub_area(nbpt, nbvmax), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_3D','Problem in allocation of variable sub_area','','')

       sub_area(:,:)=zero

       IF (printlev >= 3) THEN
         WRITE(numout,*) 'interpweight_3D:'
         WRITE(numout,*) '  Carteveg range LON:', MINVAL(lon), MAXVAL(lon)
         WRITE(numout,*) '  Carteveg range LAT:', MINVAL(lat), MAXVAL(lat)
       END IF
       
       CALL aggregate_p(nbpt, lalo, neighbours, resolution, contfrac, &
            &                iml, jml, lon, lat, mask, callsign, &
            &                nbvmax, sub_index, sub_area, ok_interpol)
       
       IF ( .NOT. ok_interpol ) THEN
          DEALLOCATE(sub_area)
          DEALLOCATE(sub_index)
          nbvmax = nbvmax * 2
       ENDIF
    ENDDO

! Getting variables thresholds
    ALLOCATE(variableusetypes(Nvariabletypes),STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_3D','Problem in allocation of variable variableusetypes','','')

    IF (variabletypes(1) == -un) THEN
      variableusetypes = (/ (ip*1., ip=1,Nvariabletypes) /)
      IF (printlev>=4) WRITE(numout,*)'interpweight_3D: Note: Using default equivalence between ',&
           'dimension in the file and dimension in the variable, for file ', filename
      varmin = 1._r_std
      varmax = Nvariabletypes*1._r_std
    ELSE
      variableusetypes = variabletypes
    END IF
    
    outvar3D = zero

! Providing the fractions of each type for each interpolated grid point
    CALL interpweight_provide_fractions3D(typefrac, nbpt, Nvariabletypes, variableusetypes,           &
      iml, jml, lml, tml, nbvmax, zero, invar3D, sub_area, sub_index, varmax, varmin,                 &
      deux, lalo, outvar3D, aoutvar)
    IF (printlev >= 5) WRITE(numout,*)'interpweight_3D end of interpweight_provide_fractions3D'

    IF (printlev>=2) THEN
       IF (COUNT(aoutvar == -1) >= AINT(nbpt*0.4)) THEN
          WRITE(numout,*) 'interpweight_3D total number of points: ', nbpt
          WRITE(numout,*) '  Interpolated: ', COUNT(aoutvar /= -1),                 &
               ' Non-interpolated: ', COUNT(aoutvar == -1)
          msg = "Less than 60% of points '" // TRIM(filename) // "' variable '" // TRIM(varname) //       & 
               "' have been succesfully interpolated on the current process !!"
          WRITE(numout,*)'   ' // TRIM(msg)
       ELSE
          msg = TRIM(filename) // "' variable: '" // TRIM(varname) // "' interpolation Done"
          WRITE(numout,*) "interpweight_3D: '" // TRIM(msg)
          WRITE(numout,*) '  Total number of points: ', nbpt
          WRITE(numout,*) '  Interpolated: ', COUNT(aoutvar /= -1),                 &
               ' Non-interpolated: ', COUNT(aoutvar == -1)
       END IF
    END IF

! Scaling to 1 the quality variable
    DO ip=1, nbpt
      IF (aoutvar(ip) /= -1) aoutvar(ip) = aoutvar(ip)/(resolution(ip,1)*resolution(ip,2))/contfrac(ip)
    END DO

    IF (ALLOCATED(lat_in)) DEALLOCATE(lat_in)
    IF (ALLOCATED(lon_in)) DEALLOCATE(lon_in)
    DEALLOCATE(lon)
    DEALLOCATE(lat)
    IF (ALLOCATED(invar2D)) DEALLOCATE(invar2D)
    IF (ALLOCATED(invar3D)) DEALLOCATE(invar3D)
    IF (ALLOCATED(invar4D)) DEALLOCATE(invar4D)
    DEALLOCATE(mask)
    IF (ALLOCATED(maskvar)) DEALLOCATE(maskvar)
    DEALLOCATE(sub_area)
    DEALLOCATE(sub_index)
    IF (ALLOCATED(resol_in)) DEALLOCATE(resol_in)

    RETURN

  END SUBROUTINE interpweight_3D

!! ================================================================================================================================
!! SUBROUTINE   : interpweight_4D
!!
!>\BRIEF        ! Interpolate any input file to the grid of the model for a 4D-variable  
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): outvar4D, aoutvar
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE interpweight_4D(nbpt, Nvariabletypes, variabletypes, lalo, resolution, neighbours,       &
    contfrac, filename, varname, inlonname, inlatname, varmin, varmax, noneg, masktype,               &
    maskvalues, maskvarname, dim1, dim2, initime, typefrac,                                           &
    maxresollon, maxresollat, outvar4D, aoutvar)

    USE ioipsl
    USE constantes_var
    USE mod_orchidee_para

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !
    !! 0.1 Input variables
    !
    INTEGER(i_std), INTENT(in)                 :: nbpt                    !! Number of points for which the data 
                                                                          !!   needs to be interpolated
    INTEGER(i_std), INTENT(in)                 :: Nvariabletypes          !! Number of types of the variable
    REAL(r_std), DIMENSION(Nvariabletypes),    &                          !! Vector of values of the types
      INTENT(in)    :: variabletypes                                      !!   (-1, not used)
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in) :: lalo                    !! Vector of latitude and longitudes 
                                                                          !! (beware of the order = 1 : latitude, 
                                                                          !!   2 : longitude)
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in) :: resolution              !! The size in km of each grid-box in X and Y
    INTEGER(i_std), DIMENSION(nbpt,NbNeighb), INTENT(in) :: neighbours    !! Vector of neighbours for each grid point 
                                                                          !!   1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW) 
                            
    REAL(r_std), DIMENSION(nbpt), INTENT(in)   :: contfrac                !! Fraction of land in each grid box.
    CHARACTER(LEN=80), INTENT(in)              :: filename                !! name of the input map to read

    CHARACTER(LEN=80), INTENT(in)              :: varname                 !! name of the variable to interpolate
    CHARACTER(LEN=80), INTENT(in)              :: inlonname, inlatname    !! names of the longitude and latitude
! Transformations to apply to the original read values from the file
    LOGICAL, INTENT(in)                        :: noneg                   !! no negative values
    CHARACTER(LEN=50), INTENT(in)              :: masktype                !! Type of masking
                                                                          !!   'nomask': no-mask is applied
                                                                          !!   'mbelow': take values below maskvalues(1)
                                                                          !!   'mabove': take values above maskvalues(1)
                                                                          !!   'msumrange': take values within 2 ranges;
                                                                          !!      maskvalues(2) <= SUM(vals(k)) <= maskvalues(1)
                                                                          !!      maskvalues(1) < SUM(vals(k)) <= maskvalues(3)
                                                                          !!        (normalized by maskedvalues(3))
                                                                          !!   'var': mask values are taken from a 
                                                                          !!     variable (>0)
    REAL(r_std), DIMENSION(3), INTENT(in)      :: maskvalues              !! values to use to mask (according to 
                                                                          !!   `masktype') 
    CHARACTER(LEN=250), INTENT(in)             :: maskvarname             !! name of the variable to use to mask 
    INTEGER, INTENT(in)                        :: dim1, dim2              !! 3/4D dimensions of the output variable
                                                                          !! dim1 = 1 fpr outvar2D
                                                                          !!   in the input file
    REAL(r_std), DIMENSION(dim1), INTENT(inout) :: varmin, varmax         !! min/max values to use for the renormalization
    CHARACTER(LEN=50), INTENT(in)              :: typefrac                !! Type of fractions retrieval:
                                                                          !!   'default': standard

    INTEGER(i_std), INTENT(in)                 :: initime                 !! Initial time to take (for input data 
                                                                          !!   with time values, -1 not used)
    REAL(r_std), INTENT(in)                    :: maxresollon,maxresollat !! lon,lat maximum resolutions (in m)
                                                                          !!   (-un, not used)
    !! 0.2 Modified variables
    !
    !! 0.3 Output variables
    !
!    REAL(r_std), INTENT(out)                   ::  laimap(nbpt,nvm,12)    !! lai read variable and re-dimensioned
    REAL(r_std), DIMENSION(nbpt,dim1,dim2), INTENT(out) :: outvar4D       !! 4D output variable and re-dimensioned
    REAL(r_std), DIMENSION(nbpt), INTENT(out)  :: aoutvar                 !! availability of input data to
                                                                          !!   interpolate output variable 
                                                                          !!   (on the nbpt space)
    !
    !! 0.4 Local variables
    !
    INTEGER(i_std)                             :: iml, jml, lml, tml, fid, ib, ip, jp, it, jj, jv
    INTEGER                                    :: inNdims                 !! Number of dimensions of the input data
    INTEGER                                    :: inNLldims               !! Number of dimensions of the lon/lat input data
    REAL(r_std), ALLOCATABLE, DIMENSION(:)     :: lat_in, lon_in          !! latitude and
                                                                          !! longitude, extract from input map
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)   :: lat, lon                !! en 2D ???
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)   :: sub_area                !! the area of the fine grid in the 
                                                                          !!   model grid ???
                                                                          !! cf src_global/interpol_help.f90, 
                                                                          !!   line 377, called "areaoverlap"
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:,:)  :: sub_index           !! the indexes from the grid boxes from the 
                                                                          !!   data that go into the 
                                                                          !!   model's boxes  
                                                                          !! cf src_global/interpol_help.f90,
                                                                          !!   line 300, called "ip"

    INTEGER                                      :: ALLOC_ERR
!    REAL(r_std), ALLOCATABLE, DIMENSION(:)       :: invar1D               !! 1D input values
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)     :: invar2D               !! 2D input values
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)   :: invar3D               !! 3D input values
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:,:) :: invar4D               !! 4D input values

    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)   :: resol_in
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:)  :: mask
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)     :: maskvar
    INTEGER(i_std)                             :: nix, njx
    INTEGER(i_std)                             :: idi, nbvmax             !! nbvmax : number of maximum vegetation map
                                                                          !! points in the GCM grid ; idi : its counter
    CHARACTER(LEN=30)                          :: callsign                !! Allows to specify which variable is beeing
                                                                          !! treated
    LOGICAL                                    :: ok_interpol             !! optionnal return of aggregate_2d
    CHARACTER(LEN=50)                          :: S1
    REAL(r_std), ALLOCATABLE, DIMENSION(:)     :: variableusetypes
    CHARACTER(LEN=256)                         :: msg

! Debug vars
    INTEGER                                    :: nb_coord, nb_var, nb_gat

    INTEGER, ALLOCATABLE, DIMENSION(:)         :: indims


!_ ================================================================================================================================
    IF (printlev >= 3) WRITE(numout, *)"In interpweight_4D filename: ",TRIM(filename),         &
      " varname:'" // TRIM(varname) // "'"

! Allocating input data
    IF (is_root_prc) THEN
       IF (printlev >=5) THEN
          WRITE(numout,*) "Entering 'interpweight_4D'. Debug mode."
          WRITE(numout,*)'(/," --> fliodmpf")'
          CALL fliodmpf (TRIM(filename))
          WRITE(numout,*)'(/," --> flioopfd")'
       ENDIF
       CALL flioopfd (TRIM(filename),fid,nb_dim=nb_coord,nb_var=nb_var,nb_gat=nb_gat)
       IF (printlev >=5) THEN
          WRITE (numout,'(" Number of coordinate        in the file : ",I2)') nb_coord
          WRITE (numout,'(" Number of variables         in the file : ",I2)') nb_var
          WRITE (numout,'(" Number of global attributes in the file : ",I2)') nb_gat
       ENDIF
    ENDIF

! Getting shape of input variable from input file
    inNdims = interpweight_get_varNdims_file(filename, varname)
    IF (is_root_prc) CALL flininfo(filename, iml, jml, lml, tml, fid)

    InputDims: SELECT CASE (inNdims)
      CASE (2)
        CALL bcast(iml)
        CALL bcast(jml)

        ALLOCATE(invar2D(iml,jml), STAT=ALLOC_ERR)
        IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4D',"Problem in allocation of variable '" //     &
          TRIM(varname) ,'','')
        IF (is_root_prc) CALL flinget(fid, TRIM(varname), iml, jml, 0, 0, 1, 1, invar2D)
        IF (printlev >= 5) THEN
          WRITE(numout,*)'  interpweight_4D input data dimensions _______'
          WRITE(numout,*)'    iml, jml:',iml, jml
        END IF

      CASE (3)
        ALLOCATE(indims(3), STAT=ALLOC_ERR)
        IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4D','Problem in allocation of indims','','')
        indims = interpweight_get_var3dims_file(filename, varname)
        lml = indims(3)
        CALL bcast(iml)
        CALL bcast(jml)
        CALL bcast(lml)

        ALLOCATE(invar3D(iml,jml, lml), STAT=ALLOC_ERR)
        IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4D',"Problem in allocation of variable '" //     &
          TRIM(varname) ,'','')
        IF (is_root_prc) CALL flinget(fid, TRIM(varname), iml, jml, lml, 0, 1, 1, invar3D)
        IF (printlev >= 5) THEN
          WRITE(numout,*)'  interpweight_4D input data dimensions:'
          WRITE(numout,*)'    iml, jml, lml:',iml, jml, lml
        END IF

      CASE (4)
        ALLOCATE(indims(4), STAT=ALLOC_ERR)
        IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4D','Problem in allocation of indims','','')
        indims = interpweight_get_var4dims_file(filename, varname)
        lml = indims(3)
        CALL bcast(iml)
        CALL bcast(jml)
        CALL bcast(lml)
!       Taken from `slowproc_update'
        IF (initime /= -1) THEN
          inNdims = 3
          IF (printlev >= 3) WRITE(numout,*) &
               'interpweight_4D: taking into account the initial reading time:',initime
          IF (initime <= 0) tml = 0

          CALL bcast(tml)
          ALLOCATE(invar3D(iml,jml,lml), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4D',"Problem in allocation of variable '" //   &
            TRIM(varname) ,'','')
          IF ( initime > 0 ) THEN
            IF (is_root_prc) THEN
              IF (initime <= tml) THEN
                CALL flinget(fid, TRIM(varname), iml, jml, lml, 1, initime, 1, invar3D)
              ELSE
                CALL flinget(fid, TRIM(varname), iml, jml, lml, 1, tml, 1, invar3D)
              ENDIF
            ENDIF
          ELSE
            IF (is_root_prc) THEN
              CALL flinget(fid, TRIM(varname), iml, jml, lml, 1, 1, 1, invar3D)
            ENDIF
          ENDIF
        ELSE
          CALL bcast(tml)
          ALLOCATE(invar4D(iml,jml,lml,tml), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4D',"Problem in allocation of variable '" //   &
            TRIM(varname) ,'','')
          IF (is_root_prc) CALL flinget(fid, TRIM(varname), iml, jml, lml, tml, 1, tml, invar4D)
        END IF

        IF (printlev >= 5) THEN
          WRITE(numout,*)'  interpweight_4D input data dimensions'
          WRITE(numout,*)'    iml, jml, lml, tml:',iml, jml, lml, tml
        END IF

      CASE DEFAULT
        WRITE (S1,'(I3)') inNdims
        CALL ipslerr_p(3,'interpweight_4D',"Input number of dimensions: '" // TRIM(S1) // "' not ready !!", &
          '','')

    END SELECT InputDims
    CALL bcast(invar4D)
  
! Getting longitudes/latitudes from input file
    inNLldims = interpweight_get_varNdims_file(filename, inlonname)
    IF (inNLldims == 1) THEN
      ! Allocate memory for latitudes
      ALLOCATE(lon_in(iml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4D','Problem in allocation of variable lon_in','','')

      ALLOCATE(lat_in(jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4D','Problem in allocation of variable lat_in','','')

      IF (is_root_prc) THEN
        CALL flinget(fid, TRIM(inlonname), iml, 0, 0, 0, 1, 1, lon_in)
        CALL flinget(fid, TRIM(inlatname), jml, 0, 0, 0, 1, 1, lat_in)
      ENDIF
      CALL bcast(lon_in)
      CALL bcast(lat_in)

      ALLOCATE(lon(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4D','Problem in allocation of variable lon','','')
      ALLOCATE(lat(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4D','Problem in allocation of variable lat','','')

      DO ip=1,iml
         lat(ip,:) = lat_in(:)
      ENDDO
      DO jp=1,jml
         lon(:,jp) = lon_in(:)
      ENDDO
    ELSE
      ! Allocate memory for longitude
      ALLOCATE(lon(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4D','Pb in allocation for lon','','')
      ! Allocate memory for latitudes
      ALLOCATE(lat(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4D','Pb in allocation for lat','','')

      IF (is_root_prc) THEN
        CALL flinget(fid, TRIM(inlonname), iml, jml, 0, 0, 1, 1, lon)
        CALL flinget(fid, TRIM(inlatname), iml, jml, 0, 0, 1, 1, lat)
      ENDIF
    END IF
    CALL bcast(lon)
    CALL bcast(lat)

    IF (is_root_prc) THEN
      ALLOCATE(resol_in(iml,jml,2), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4D','Problem in allocation of variable resol_in','','')
    END IF

    IF (noneg) THEN
      IF (is_root_prc) THEN
        SELECT CASE (inNdims)
          CASE(2)
            CALL interpweight_modifying_input2D(iml, jml, 0, 0, 1, 1, zero, invar2D)
          CASE(3)
            CALL interpweight_modifying_input3D(iml, jml, lml, 0, 1, 1, zero, invar3D)
          CASE(4)
            CALL interpweight_modifying_input4D(iml, jml, lml, tml, 1, tml, zero, invar4D)
        END SELECT       
      END IF
    END IF

    IF (printlev >= 5) WRITE(numout,*)'  interpweight_4D modifying input ended!'

    ALLOCATE(mask(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4D','Problem in allocation of variable mask','','')

    IF (TRIM(masktype) == 'var') THEN
      ALLOCATE(maskvar(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4D','Problem in allocation of variable maskvar','','')
      maskvar(:,:) = 0
! Reads the mask from the file
      IF (is_root_prc) THEN
! ATTTENTION!! One can not use flinget in this way. It always spect to have a REAL variable!
        CALL flinget(fid, TRIM(maskvarname), iml, jml, 0, 0, 1, 1, maskvar)
        WHERE (maskvar > zero)
          mask = un 
        END WHERE
      END IF
      CALL bcast(mask)
    ELSE
      SELECT CASE (inNdims)
        CASE(2)
          CALL interpweight_masking_input2D(iml, jml, 0, 0, 0, 0, masktype, maskvalues, invar2D,      &
            mask)
        CASE(3)
          CALL interpweight_masking_input3D(iml, jml, 0, 0, lml, 0, masktype, maskvalues, invar3D,    &
            mask)
        CASE(4)
          CALL interpweight_masking_input4D(iml, jml, 0, 0, lml, tml, masktype, maskvalues, invar4D,  &
            mask)
      END SELECT
    END IF
    IF (printlev >= 5) WRITE(numout,*)'  interpweight_4D masking input ended!'

    IF (is_root_prc) CALL flinclo(fid)

    ! Computation of nbvmax
    ! In 'condveg_soilalb, slowproc_update' nbvmax=200
    !
    ! The number of maximum vegetation map points in the GCM grid is estimated.
    ! Some lmargin is taken.
    !
    IF (is_root_prc) THEN
      IF ( (maxresollon /= -un) .AND. (maxresollat /= -un) ) THEN
!       Assuming regular lon/lat projection for the input data but it is not caluclated assuming a maximum  &
!         lon,lat resolution given as input paramater
         nix=INT(MAXVAL(resolution(:,1))/maxresollon)+2
         njx=INT(MAXVAL(resolution(:,2))/maxresollat)+2
      ELSE
!       Assuming regular lon/lat projection for the input data
        CALL interpweight_calc_resolution_in(lon, lat, iml, jml, mincos, R_Earth, pi, resol_in)
        nix=INT(MAXVAL(resolution(:,1))/MAXVAL(resol_in(:,:,1)))+2
        njx=INT(MAXVAL(resolution(:,2))/MAXVAL(resol_in(:,:,2)))+2
      END IF
      nbvmax = MAX(nix*njx, 200)
    END IF  
    CALL bcast(nbvmax)
    CALL bcast(nix)
    CALL bcast(njx)

    callsign = TRIM(varname) // ' map'

    ok_interpol = .FALSE.
    DO WHILE ( .NOT. ok_interpol )
       IF (printlev>=3) WRITE(numout,*) "Projection arrays for ",callsign," : "
       IF (printlev>=3) WRITE(numout,*) "nbvmax = ",nbvmax, nix, njx

       ALLOCATE(sub_index(nbpt, nbvmax, 2), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4D','Problem in allocation of variable sub_index','','')

       sub_index(:,:,:)=zero

       ALLOCATE(sub_area(nbpt, nbvmax), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4D','Problem in allocation of variable sub_area','','')

       sub_area(:,:)=zero
       CALL aggregate_p(nbpt, lalo, neighbours, resolution, contfrac, &
            &                iml, jml, lon, lat, mask, callsign, &
            &                nbvmax, sub_index, sub_area, ok_interpol)
       
       IF ( .NOT. ok_interpol ) THEN
          DEALLOCATE(sub_area)
          DEALLOCATE(sub_index)
          nbvmax = nbvmax * 2
       ENDIF
    ENDDO

! Getting variables thresholds
    ALLOCATE(variableusetypes(Nvariabletypes), STAT=ALLOC_ERR)
            IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4D','Problem in allocation of variable variableusetypesr','','')

    IF (variabletypes(1) == -un) THEN
      variableusetypes = (/ (ip*1., ip=1,Nvariabletypes) /)
      IF (printlev>=4) WRITE(numout,*)'  interpweight_4D: Note: Using default ',&
           'equivalence between dimension in the file and dimension in the variable, for file ', filename
      varmin = 1._r_std
      varmax = Nvariabletypes*1._r_std
    ELSE
      variableusetypes = variabletypes
    END IF
    
    outvar4D = zero

! Doing the interpolation as function of the type fo the values
    CALL interpweight_provide_fractions4D(typefrac, nbpt, Nvariabletypes, variableusetypes,           &
      iml, jml, lml, tml, nbvmax, zero, invar4D, sub_area, sub_index, varmax, varmin,                 &
      deux, lalo, outvar4D, aoutvar)
    IF (printlev >= 5) WRITE(numout,*)'  interpweight_4D end of interpweight_provide_fractions4D'

    IF (printlev>=2) THEN
       IF (COUNT(aoutvar == -1) >= AINT(nbpt*0.4)) THEN
          WRITE(numout,*) 'interpweight_4D total number of points: ', nbpt
          WRITE(numout,*) '  Interpolated: ', COUNT(aoutvar /= -1),                 &
               ' Non-interpolated: ', COUNT(aoutvar == -1)
          msg = "Less than 60% of points '" // TRIM(filename) // "' variable '" // TRIM(varname) //       & 
               "' have been succesfully interpolated !!"
          WRITE(numout,*)'  ' // TRIM(msg)
       ELSE
          msg = TRIM(filename) // "' variable: '" // TRIM(varname) // "' interpolation Done _______"
          WRITE(numout,*) "interpweight_4D: '" // TRIM(msg)
          WRITE(numout,*) '  Total number of points: ', nbpt
          WRITE(numout,*) '  Interpweight_4D interpolated: ', COUNT(aoutvar /= -1),                 &
               ' Non-interpolated: ', COUNT(aoutvar == -1)
       END IF
    END IF

! Scaling to 1 the quality variable
    DO ip=1, nbpt
      IF (aoutvar(ip) /= -1) aoutvar(ip) = aoutvar(ip)/(resolution(ip,1)*resolution(ip,2))/contfrac(ip)
    END DO

    IF (ALLOCATED(lat_in)) DEALLOCATE(lat_in)
    IF (ALLOCATED(lon_in)) DEALLOCATE(lon_in)
    DEALLOCATE(lon)
    DEALLOCATE(lat)
    IF (ALLOCATED(invar2D)) DEALLOCATE(invar2D)
    IF (ALLOCATED(invar3D)) DEALLOCATE(invar3D)
    IF (ALLOCATED(invar4D)) DEALLOCATE(invar4D)
    DEALLOCATE(mask)
    IF (ALLOCATED(maskvar)) DEALLOCATE(maskvar)
    DEALLOCATE(sub_area)
    DEALLOCATE(sub_index)
    IF (ALLOCATED(resol_in)) DEALLOCATE(resol_in)

    RETURN

  END SUBROUTINE interpweight_4D

!! ================================================================================================================================
!! SUBROUTINE   : interpweight_2Dcont_tstep
!!
!>\BRIEF        ! Interpolate any input file to the grid of the model for a continuos 2D-variable  
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): outvar2D, aoutvar
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE interpweight_2Dcont_tstep(nbpt, lalo, resolution, neighbours,                      &
    contfrac, filename, varname, inlonname, inlatname, varmin, varmax, noneg, masktype,               &
    maskvalues, maskvarname, initime, typefrac, defaultvalue, defaultNOvalue,                         &
    outvar2D, aoutvar)

    USE ioipsl
    USE constantes_var
    USE mod_orchidee_para

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !
    !! 0.1 Input variables
    !
    INTEGER(i_std), INTENT(in)                 :: nbpt                    !! Number of points for which the data 
                                                                          !!   needs to be interpolated
!    INTEGER(i_std), INTENT(in)                 :: dim1, dim2              !! 4D complementary dimension sizes
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in) :: lalo                    !! Vector of latitude and longitudes 
                                                                          !! (beware of the order = 1 : latitude, 
                                                                          !!   2 : longitude)
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in) :: resolution              !! The size in km of each grid-box in X and Y
    INTEGER(i_std), DIMENSION(nbpt,NbNeighb), INTENT(in) :: neighbours    !! Vector of neighbours for each grid point 
                                                                          !!   1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW) 
                            
    REAL(r_std), DIMENSION(nbpt), INTENT(in)   :: contfrac                !! Fraction of land in each grid box.
    CHARACTER(LEN=*), INTENT(in)              :: filename                !! name of the input map to read

    CHARACTER(LEN=*), INTENT(in)              :: varname                 !! name of the variable to interpolate
    CHARACTER(LEN=*), INTENT(in)              :: inlonname, inlatname    !! names of the longitude and latitude
    REAL(r_std), INTENT(in)                    :: varmin, varmax          !! min/max values to use for the
                                                                          !!   renormalization
    REAL(r_std), INTENT(in)                    :: defaultvalue            !! default interpolated value 
    REAL(r_std), INTENT(in)                    :: defaultNOvalue          !! default interpolated NO-value 

! Transformations to apply to the original read values from the file
    LOGICAL, INTENT(in)                        :: noneg                   !! no negative values
    CHARACTER(LEN=50), INTENT(in)              :: masktype                !! Type of masking
                                                                          !!   'nomask': no-mask is applied
                                                                          !!   'mbelow': take values below maskvalues(1)
                                                                          !!   'mabove': take values above maskvalues(1)
                                                                          !!   'msumrange': take values within 2 ranges;
                                                                          !!      maskvalues(2) <= SUM(vals(k)) <= maskvalues(1)
                                                                          !!      maskvalues(1) < SUM(vals(k)) <= maskvalues(3)
                                                                          !!        (normalized by maskedvalues(3))
                                                                          !!   'var': mask values are taken from a 
                                                                          !!     variable (1/0)
    REAL(r_std), DIMENSION(3), INTENT(in)      :: maskvalues              !! values to use to mask (according to 
                                                                          !!   `masktype') 
    CHARACTER(LEN=250), INTENT(in)             :: maskvarname             !! name of the variable to use to mask 
                                                                          !!   (according to `masktype') 
    CHARACTER(LEN=50), INTENT(in)              :: typefrac                !! Type of fraction retrieval:
                                                                          !!   'default': standard

    INTEGER(i_std), INTENT(in)                 :: initime                 !! Initial time to take (for input data 
                                                                          !!   with time values, -1: not used)
    !! 0.2 Modified variables
    !
    !! 0.3 Output variables
    !
    REAL(r_std), DIMENSION(:,:), INTENT(out)  :: outvar2D                 !! 2D output variable and re-dimensioned kjpindex, tsteps
    REAL(r_std), DIMENSION(:), INTENT(out)  :: aoutvar                    !! availability of input data to
                                                                          !!   interpolate output variable 
                                                                          !!   (on the nbpt space)
    !
    !! 0.4 Local variables
    !
    INTEGER(i_std)                             :: iml, jml, lml, tml, fid, ib, ip, jp, it, jj, jv
    INTEGER                                    :: inNdims                 !! Number of dimensions of the input data
    INTEGER                                    :: inNLldims               !! Number of dimensions of the lon/lat input data
    REAL(r_std), ALLOCATABLE, DIMENSION(:)     :: lat_in, lon_in          !! latitude and
                                                                          !! longitude, extract from input map
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)   :: lat, lon                !! en 2D ???
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)   :: sub_area                !! the area of the fine grid in the 
                                                                          !!   model grid ???
                                                                          !! cf src_global/interpol_help.f90, 
                                                                          !!   line 377, called "areaoverlap"
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:,:)  :: sub_index           !! the indexes from the grid boxes from the 
                                                                          !!   data that go into the 
                                                                          !!   model's boxes  
                                                                          !! cf src_global/interpol_help.f90,
                                                                          !!   line 300, called "ip"

!    REAL(r_std), ALLOCATABLE, DIMENSION(:)       :: invar1D               !! 1D input values
!    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)     :: invar2D               !! 2D input values
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)   :: invar3D               !! 3D input values
!    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:,:) :: invar4D               !! 4D input values

    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)   :: resol_in
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:)  :: mask
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)     :: maskvar
    INTEGER(i_std)                             :: nix, njx
    INTEGER(i_std)                             :: idi, nbvmax             !! nbvmax : number of maximum vegetation map
                                                                          !! points in the GCM grid ; idi : its counter
    CHARACTER(LEN=30)                          :: callsign                !! Allows to specify which variable is beeing
                                                                          !! treated
    LOGICAL                                    :: ok_interpol             !! optionnal return of aggregate_2d
    INTEGER                                    :: ALLOC_ERR 

    CHARACTER(LEN=50)                          :: S1
    CHARACTER(LEN=250)                         :: msg


! Debug vars
    INTEGER                                    :: nb_coord, nb_var, nb_gat

!_ ================================================================================================================================
    IF (printlev >= 3) WRITE(numout,*)"In interpweight_2Dcont filename: ",TRIM(filename),          &
      " varname:'" // TRIM(varname) // "'"
! Allocating input data
    IF (is_root_prc) THEN
       IF (printlev >=5) THEN
          WRITE(numout,*) "Entering 'interpweight_2Dcont'. Debug mode."
          WRITE(numout,*)'(/," --> fliodmpf")'
          CALL fliodmpf (TRIM(filename))
          WRITE(numout,*)'(/," --> flioopfd")'
       ENDIF
       CALL flioopfd (TRIM(filename),fid,nb_dim=nb_coord,nb_var=nb_var,nb_gat=nb_gat)
       IF (printlev >=5) THEN
          WRITE (numout,'(" Number of coordinate        in the file : ",I2)') nb_coord
          WRITE (numout,'(" Number of variables         in the file : ",I2)') nb_var
          WRITE (numout,'(" Number of global attributes in the file : ",I2)') nb_gat
       ENDIF
    ENDIF

! Getting shape of input variable from input file
    inNdims = interpweight_get_varNdims_file(filename, varname)
    IF (is_root_prc) CALL flininfo(filename, iml, jml, lml, tml, fid)

    CALL bcast(iml)
    CALL bcast(jml)
    CALL bcast(lml)
    CALL bcast(tml)
    CALL bcast(fid)

    IF (lml .GT. 0) CALL ipslerr_p(3, 'interpweight_2Dcont_tstep','No expected dimension for lml ','But Found:',lml)
    IF (tml .LE. 0) CALL ipslerr_p(3, 'interpweight_2Dcont_tstep','Expected dimension for tml','But Found:',tml)

    ALLOCATE(invar3D(iml,jml, tml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2Dcont',"Problem in allocation of variable '" //     &
          TRIM(varname) ,'','')
    IF (is_root_prc) CALL flinget(fid, TRIM(varname), iml, jml, lml, tml, 1, tml, invar3D)
    IF (printlev >= 5) THEN
          WRITE(numout,*)'  interpweight_2Dcont input data dimensions:'
          WRITE(numout,*)'    iml, jml, lml:',iml, jml, lml
    END IF
    CALL bcast(invar3D)

    IF (printlev >= 5) WRITE(numout,*)'  interpweight_2Dcont getting variable ended!'

! Getting longitudes/latitudes from input file
    inNLldims = interpweight_get_varNdims_file(filename, inlonname)

    IF (inNLldims == 1) THEN
      ! Allocate memory for latitudes
      ALLOCATE(lon_in(iml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2Dcont','Problem in allocation of variable lon_in','','')

      ALLOCATE(lat_in(jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2Dcont','Problem in allocation of variable lat_in','','')

      IF (is_root_prc) THEN
        CALL flinget(fid, TRIM(inlonname), iml, 0, 0, 0, 1, 1, lon_in)
        CALL flinget(fid, TRIM(inlatname), jml, 0, 0, 0, 1, 1, lat_in)
      ENDIF
      CALL bcast(lon_in)
      CALL bcast(lat_in)

      ALLOCATE(lon(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2Dcont','Problem in allocation of variable lon','','')
      ALLOCATE(lat(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2Dcont','Problem in allocation of variable lat','','')

      DO ip=1,iml
         lat(ip,:) = lat_in(:)
      ENDDO
      DO jp=1,jml
         lon(:,jp) = lon_in(:)
      ENDDO
    ELSE
      ! Allocate memory for longitude
      ALLOCATE(lon(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2Dcont','Pb in allocation for lon','','')
      ! Allocate memory for latitudes
      ALLOCATE(lat(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2Dcont','Pb in allocation for lat','','')

      IF (is_root_prc) THEN
        CALL flinget(fid, TRIM(inlonname), iml, jml, 0, 0, 1, 1, lon)
        CALL flinget(fid, TRIM(inlatname), iml, jml, 0, 0, 1, 1, lat)
      ENDIF
    END IF
    CALL bcast(lon)
    CALL bcast(lat)

    IF (is_root_prc) THEN
      ALLOCATE(resol_in(iml,jml,2), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2Dcont','Problem in allocation of variable resol_in','','')
    END IF

    IF (noneg) THEN
      IF (is_root_prc) THEN
         CALL interpweight_modifying_input3D(iml, jml, lml, tml, 1, 1, zero, invar3D)
      END IF
    END IF
    IF (printlev >= 5) WRITE(numout,*)'  interpweight_2Dcont modifying input ended!'

    !
    ! Consider all points a priori
    !
    ALLOCATE(mask(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2Dcont','Problem in allocation of variable mask','','')
    mask(:,:) = 0

    IF (TRIM(masktype) == 'var') THEN
      ALLOCATE(maskvar(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2Dcont','Problem in allocation of variable maskvar','','')
      maskvar(:,:) = 0
! Reads the mask from the file
      IF (is_root_prc) THEN
! ATTTENTION!! One can not use flinget in this way. It always spect to have a REAL variable!
        CALL flinget(fid, TRIM(maskvarname), iml, jml, 0, 0, 1, 1, maskvar)
        WHERE (maskvar > zero)
          mask = un 
        END WHERE
      END IF
      CALL bcast(mask)
    ELSE
      CALL interpweight_masking_input2D(iml, jml, 0, 0, 0, 0, masktype, maskvalues, invar3D,    &
          mask)
    END IF
    IF (printlev >= 5) WRITE(numout,*)'  interpweight_2Dcont masking input ended!'
    IF (is_root_prc) CALL flinclo(fid)

    IF (printlev>=5) WRITE(numout,*)"  interpweight_2Dcont '" // TRIM(typefrac) // "' interpolation"

! Assuming regular lon/lat projection for the input data
    IF (is_root_prc) THEN
      CALL interpweight_calc_resolution_in(lon, lat, iml, jml, mincos, R_Earth, pi, resol_in)
    END IF

!   Computation of nbvmax
!   In 'condveg_soilalb, slowproc_update' nbvmax=200

!
!   The number of maximum vegetation map points in the GCM grid is estimated.
!   Some lmargin is taken.
!
    IF (is_root_prc) THEN
       nix=INT(MAXVAL(resolution(:,1))/MAXVAL(resol_in(:,:,1)))+2
       njx=INT(MAXVAL(resolution(:,2))/MAXVAL(resol_in(:,:,2)))+2
       nbvmax = MAX(nix*njx, 200)
    ENDIF
    CALL bcast(nbvmax)
    CALL bcast(nix)
    CALL bcast(njx)

    callsign = TRIM(varname) // ' map'

    ok_interpol = .FALSE.
    DO WHILE ( .NOT. ok_interpol )
       IF (printlev>=3) WRITE(numout,*) "Projection arrays for ",callsign," : "
       IF (printlev>=3) WRITE(numout,*) "nbvmax = ",nbvmax, nix, njx

       ALLOCATE(sub_index(nbpt, nbvmax, 2), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2Dcont','Problem in allocation of variable sub_index','','')

       sub_index(:,:,:)=0

       ALLOCATE(sub_area(nbpt, nbvmax), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2Dcont','Problem in allocation of variable sub_area','','')

       sub_area(:,:)=zero
       IF (printlev>=3) WRITE(numout,*) 'Input data range LON:', MINVAL(lon), MAXVAL(lon)
       IF (printlev>=3) WRITE(numout,*) 'Input data range LAT:', MINVAL(lat), MAXVAL(lat)

       CALL aggregate_p(nbpt, lalo, neighbours, resolution, contfrac, &
            &                iml, jml, lon, lat, mask, callsign, &
            &                nbvmax, sub_index, sub_area, ok_interpol)
       
       IF ( .NOT. ok_interpol ) THEN
          DEALLOCATE(sub_area)
          DEALLOCATE(sub_index)
          nbvmax = nbvmax * 2
       ENDIF
    ENDDO
    IF (printlev>=4) WRITE (numout,*) '  interpweight_2Dcont: read/allocate OK'

    outvar2D = zero
! Providing the interpolated grid point

    CALL interpweight_provide_interpolation2D_tsteps(typefrac, nbpt, 0, 0, iml, jml, lml, tml,               &
      nbvmax, zero, invar3D, sub_area, sub_index, varmax, varmin, un, deux, lalo,                     &
      defaultvalue, defaultNOvalue, outvar2D, aoutvar)
    IF (printlev >= 5) WRITE(numout,*)'  interpweight_2Dcont end of interpweight_provide_interpolation2D_tsteps'

    IF (printlev>=2) THEN
       IF (COUNT(aoutvar == -1) >= AINT(nbpt*0.4)) THEN
          WRITE(numout,*) 'interpweight_2Dcont total number of points: ', nbpt
          WRITE(numout,*) '  Interpolated: ', COUNT(aoutvar /= -1),                 &
               ' non-interpolated: ', COUNT(aoutvar == -1)
          msg = "Less than 60% of points '" // TRIM(filename) // "' variable '" // TRIM(varname) //       &
               "' have been succesfully interpolated !!"
          WRITE(numout,*) '  ' // TRIM(msg)
       ELSE
          msg = TRIM(filename) // "' variable: '" // TRIM(varname) // "' interpolation Done _______"
          WRITE(numout,*) "interpweight_2Dcont: '" // TRIM(msg)
          WRITE(numout,*) '  Total number of points: ', nbpt
          WRITE(numout,*) '  Interpolated: ', COUNT(aoutvar /= -1),                 &
               ' Non-interpolated: ', COUNT(aoutvar == -1)
       END IF
    END IF

! Scaling to 1 the quality variable
    DO ip=1, nbpt
      IF (aoutvar(ip) /= -1) aoutvar(ip) = aoutvar(ip)/(resolution(ip,1)*resolution(ip,2))/contfrac(ip)
    END DO

    IF (printlev>=5) WRITE(numout,*) '  interpweight_2Dcont: Interpolation Done'

    IF (ALLOCATED(lat_in)) DEALLOCATE(lat_in)
    IF (ALLOCATED(lon_in)) DEALLOCATE(lon_in)
    DEALLOCATE(lon)
    DEALLOCATE(lat)
!    IF (ALLOCATED(invar2D)) DEALLOCATE(invar2D)
    IF (ALLOCATED(invar3D)) DEALLOCATE(invar3D)
!    IF (ALLOCATED(invar4D)) DEALLOCATE(invar4D)
    DEALLOCATE(mask)
    IF (ALLOCATED(maskvar)) DEALLOCATE(maskvar)
    IF (ALLOCATED(sub_area)) DEALLOCATE(sub_area)
    IF (ALLOCATED(sub_index)) DEALLOCATE(sub_index)
    IF (ALLOCATED(resol_in)) DEALLOCATE(resol_in)

    RETURN

  END SUBROUTINE interpweight_2Dcont_tstep

!! ================================================================================================================================
!! SUBROUTINE   : interpweight_3Dcont_tstep
!!
!>\BRIEF        ! Interpolate any input file to the grid of the model for a continuos 2D-variable  
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): outvar3D, aoutvar
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE interpweight_3Dcont_tstep(nbpt, lalo, resolution, neighbours,                      &
    contfrac, filename, varname, inlonname, inlatname, varmin, varmax, noneg, masktype,               &
    maskvalues, maskvarname, initime, typefrac, defaultvalue, defaultNOvalue,                         &
    outvar3D, aoutvar)

    USE ioipsl
    USE constantes_var
    USE mod_orchidee_para

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !
    !! 0.1 Input variables
    !
    INTEGER(i_std), INTENT(in)                 :: nbpt                    !! Number of points for which the data 
                                                                          !!   needs to be interpolated
!    INTEGER(i_std), INTENT(in)                 :: dim1, dim2              !! 4D complementary dimension sizes
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in) :: lalo                    !! Vector of latitude and longitudes 
                                                                          !! (beware of the order = 1 : latitude, 
                                                                          !!   2 : longitude)
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in) :: resolution              !! The size in km of each grid-box in X and Y
    INTEGER(i_std), DIMENSION(nbpt,NbNeighb), INTENT(in) :: neighbours    !! Vector of neighbours for each grid point 
                                                                          !!   1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW) 
                            
    REAL(r_std), DIMENSION(nbpt), INTENT(in)   :: contfrac                !! Fraction of land in each grid box.
    CHARACTER(LEN=*), INTENT(in)              :: filename                !! name of the input map to read

    CHARACTER(LEN=*), INTENT(in)              :: varname                 !! name of the variable to interpolate
    CHARACTER(LEN=*), INTENT(in)              :: inlonname, inlatname    !! names of the longitude and latitude
    REAL(r_std), INTENT(in)                    :: varmin, varmax          !! min/max values to use for the
                                                                          !!   renormalization
    REAL(r_std), INTENT(in)                    :: defaultvalue            !! default interpolated value 
    REAL(r_std), INTENT(in)                    :: defaultNOvalue          !! default interpolated NO-value 

! Transformations to apply to the original read values from the file
    LOGICAL, INTENT(in)                        :: noneg                   !! no negative values
    CHARACTER(LEN=50), INTENT(in)              :: masktype                !! Type of masking
                                                                          !!   'nomask': no-mask is applied
                                                                          !!   'mbelow': take values below maskvalues(1)
                                                                          !!   'mabove': take values above maskvalues(1)
                                                                          !!   'msumrange': take values within 2 ranges;
                                                                          !!      maskvalues(2) <= SUM(vals(k)) <= maskvalues(1)
                                                                          !!      maskvalues(1) < SUM(vals(k)) <= maskvalues(3)
                                                                          !!        (normalized by maskedvalues(3))
                                                                          !!   'var': mask values are taken from a 
                                                                          !!     variable (1/0)
    REAL(r_std), DIMENSION(3), INTENT(in)      :: maskvalues              !! values to use to mask (according to 
                                                                          !!   `masktype') 
    CHARACTER(LEN=250), INTENT(in)             :: maskvarname             !! name of the variable to use to mask 
                                                                          !!   (according to `masktype') 
    CHARACTER(LEN=50), INTENT(in)              :: typefrac                !! Type of fraction retrieval:
                                                                          !!   'default': standard

    INTEGER(i_std), INTENT(in)                 :: initime                 !! Initial time to take (for input data 
                                                                          !!   with time values, -1: not used)
    !! 0.2 Modified variables
    !
    !! 0.3 Output variables
    !                   nbpt, dim1, tsteps
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)  :: outvar3D                !! 2D output variable and re-dimensioned
    REAL(r_std), DIMENSION(nbpt), INTENT(out)  :: aoutvar                 !! availability of input data to
                                                                          !!   interpolate output variable 
                                                                          !!   (on the nbpt space)
    !
    !! 0.4 Local variables
    !
    INTEGER(i_std)                             :: iml, jml, lml, tml, fid, ib, ip, jp, it, jj, jv
    INTEGER                                    :: inNdims                 !! Number of dimensions of the input data
    INTEGER                                    :: inNLldims               !! Number of dimensions of the lon/lat input data
    REAL(r_std), ALLOCATABLE, DIMENSION(:)     :: lat_in, lon_in          !! latitude and
                                                                          !! longitude, extract from input map
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)   :: lat, lon                !! en 2D ???
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)   :: sub_area                !! the area of the fine grid in the 
                                                                          !!   model grid ???
                                                                          !! cf src_global/interpol_help.f90, 
                                                                          !!   line 377, called "areaoverlap"
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:,:)  :: sub_index           !! the indexes from the grid boxes from the 
                                                                          !!   data that go into the 
                                                                          !!   model's boxes  
                                                                          !! cf src_global/interpol_help.f90,
                                                                          !!   line 300, called "ip"

!    REAL(r_std), ALLOCATABLE, DIMENSION(:)       :: invar1D               !! 1D input values
!    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)     :: invar2D               !! 2D input values
!    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)   :: invar3D               !! 3D input values
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:,:) :: invar4D               !! 4D input values (lat, lon, dim1, tsteps)

    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)   :: resol_in
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:)  :: mask
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)     :: maskvar
    INTEGER(i_std)                             :: nix, njx
    INTEGER(i_std)                             :: idi, nbvmax             !! nbvmax : number of maximum vegetation map
                                                                          !! points in the GCM grid ; idi : its counter
    CHARACTER(LEN=30)                          :: callsign                !! Allows to specify which variable is beeing
                                                                          !! treated
    LOGICAL                                    :: ok_interpol             !! optionnal return of aggregate_2d
    INTEGER                                    :: ALLOC_ERR 

    CHARACTER(LEN=50)                          :: S1
    CHARACTER(LEN=250)                         :: msg


! Debug vars
    INTEGER                                    :: nb_coord, nb_var, nb_gat

!_ ================================================================================================================================
    IF (printlev >= 3) WRITE(numout,*)"In interpweight_3Dcont filename: ",TRIM(filename),          &
      " varname:'" // TRIM(varname) // "'"
! Allocating input data
    IF (is_root_prc) THEN
       IF (printlev >=5) THEN
          WRITE(numout,*) "Entering 'interpweight_3Dcont'. Debug mode."
          WRITE(numout,*)'(/," --> fliodmpf")'
          CALL fliodmpf (TRIM(filename))
          WRITE(numout,*)'(/," --> flioopfd")'
       ENDIF
       CALL flioopfd (TRIM(filename),fid,nb_dim=nb_coord,nb_var=nb_var,nb_gat=nb_gat)
       IF (printlev >=5) THEN
          WRITE (numout,'(" Number of coordinate        in the file : ",I2)') nb_coord
          WRITE (numout,'(" Number of variables         in the file : ",I2)') nb_var
          WRITE (numout,'(" Number of global attributes in the file : ",I2)') nb_gat
       ENDIF
    ENDIF

! Getting shape of input variable from input file
    inNdims = interpweight_get_varNdims_file(filename, varname)
    IF (is_root_prc) CALL flininfo(filename, iml, jml, lml, tml, fid)

    CALL bcast(iml)
    CALL bcast(jml)
    CALL bcast(lml)
    CALL bcast(tml)
    CALL bcast(fid)

    IF (lml .LE. 0) CALL ipslerr_p(3, 'interpweight_3Dcont_tstep','Expected dimension for lml ','But Found:',lml)
    IF (tml .LE. 0) CALL ipslerr_p(3, 'interpweight_3Dcont_tstep','Expected dimension for tml','But Found:',tml)

    ALLOCATE(invar4D(iml, jml, lml, tml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_3Dcont_tsteps',"Problem in allocation of variable '" //     &
          TRIM(varname) ,'','')
    IF (is_root_prc) CALL flinget(fid, TRIM(varname), iml, jml, lml, tml, 1, tml, invar4D)
    IF (printlev >= 5) THEN
          WRITE(numout,*)'  interpweight_3Dcont_tsteps input data dimensions:'
          WRITE(numout,*)'    iml, jml, lml:',iml, jml, lml
    END IF
    CALL bcast(invar4D)

    IF (printlev >= 5) WRITE(numout,*)'  interpweight_3Dcont_tsteps getting variable ended!'

! Getting longitudes/latitudes from input file
    inNLldims = interpweight_get_varNdims_file(filename, inlonname)

    IF (inNLldims == 1) THEN
      ! Allocate memory for latitudes
      ALLOCATE(lon_in(iml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_3Dcont_tsteps','Problem in allocation of variable lon_in','','')

      ALLOCATE(lat_in(jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_3Dcont_tsteps','Problem in allocation of variable lat_in','','')

      IF (is_root_prc) THEN
        CALL flinget(fid, TRIM(inlonname), iml, 0, 0, 0, 1, 1, lon_in)
        CALL flinget(fid, TRIM(inlatname), jml, 0, 0, 0, 1, 1, lat_in)
      ENDIF
      CALL bcast(lon_in)
      CALL bcast(lat_in)

      ALLOCATE(lon(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_3Dcont_tsteps','Problem in allocation of variable lon','','')
      ALLOCATE(lat(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_3Dcont_tsteps','Problem in allocation of variable lat','','')

      DO ip=1,iml
         lat(ip,:) = lat_in(:)
      ENDDO
      DO jp=1,jml
         lon(:,jp) = lon_in(:)
      ENDDO
    ELSE
      ! Allocate memory for longitude
      ALLOCATE(lon(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_3Dcont_tsteps','Pb in allocation for lon','','')
      ! Allocate memory for latitudes
      ALLOCATE(lat(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_3Dcont_tsteps','Pb in allocation for lat','','')

      IF (is_root_prc) THEN
        CALL flinget(fid, TRIM(inlonname), iml, jml, 0, 0, 1, 1, lon)
        CALL flinget(fid, TRIM(inlatname), iml, jml, 0, 0, 1, 1, lat)
      ENDIF
    END IF
    CALL bcast(lon)
    CALL bcast(lat)

    IF (is_root_prc) THEN
      ALLOCATE(resol_in(iml,jml,2), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_3Dcont_tsteps','Problem in allocation of variable resol_in','','')
    END IF

    IF (noneg) THEN
      IF (is_root_prc) THEN
         CALL interpweight_modifying_input3D(iml, jml, lml, tml, 1, 1, zero, invar4D)
      END IF
    END IF
    IF (printlev >= 5) WRITE(numout,*)'  interpweight_3Dcont_tsteps modifying input ended!'

    !
    ! Consider all points a priori
    !
    ALLOCATE(mask(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_3Dcont_tsteps','Problem in allocation of variable mask','','')
    mask(:,:) = 0

    IF (TRIM(masktype) == 'var') THEN
      ALLOCATE(maskvar(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_3Dcont_tsteps','Problem in allocation of variable maskvar','','')
      maskvar(:,:) = 0
! Reads the mask from the file
      IF (is_root_prc) THEN
! ATTTENTION!! One can not use flinget in this way. It always spect to have a REAL variable!
        CALL flinget(fid, TRIM(maskvarname), iml, jml, 0, 0, 1, 1, maskvar)
        WHERE (maskvar > zero)
          mask = un 
        END WHERE
      END IF
      CALL bcast(mask)
    ELSE
      CALL interpweight_masking_input2D(iml, jml, 0, 0, 0, 0, masktype, maskvalues, invar4D,    &
          mask)
    END IF
    IF (printlev >= 5) WRITE(numout,*)'  interpweight_3Dcont_tsteps masking input ended!'
    IF (is_root_prc) CALL flinclo(fid)

    IF (printlev>=5) WRITE(numout,*)"  interpweight_3Dcont_tsteps '" // TRIM(typefrac) // "' interpolation"

! Assuming regular lon/lat projection for the input data
    IF (is_root_prc) THEN
      CALL interpweight_calc_resolution_in(lon, lat, iml, jml, mincos, R_Earth, pi, resol_in)
    END IF

!   Computation of nbvmax
!   In 'condveg_soilalb, slowproc_update' nbvmax=200

!
!   The number of maximum vegetation map points in the GCM grid is estimated.
!   Some lmargin is taken.
!
    IF (is_root_prc) THEN
       nix=INT(MAXVAL(resolution(:,1))/MAXVAL(resol_in(:,:,1)))+2
       njx=INT(MAXVAL(resolution(:,2))/MAXVAL(resol_in(:,:,2)))+2
       nbvmax = MAX(nix*njx, 200)
    ENDIF
    CALL bcast(nbvmax)
    CALL bcast(nix)
    CALL bcast(njx)

    callsign = TRIM(varname) // ' map'

    ok_interpol = .FALSE.
    DO WHILE ( .NOT. ok_interpol )
       IF (printlev>=3) WRITE(numout,*) "Projection arrays for ",callsign," : "
       IF (printlev>=3) WRITE(numout,*) "nbvmax = ",nbvmax, nix, njx

       ALLOCATE(sub_index(nbpt, nbvmax, 2), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_3Dcont_tsteps','Problem in allocation of variable sub_index','','')

       sub_index(:,:,:)=0

       ALLOCATE(sub_area(nbpt, nbvmax), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_3Dcont_tsteps','Problem in allocation of variable sub_area','','')

       sub_area(:,:)=zero
       IF (printlev>=3) WRITE(numout,*) 'Input data range LON:', MINVAL(lon), MAXVAL(lon)
       IF (printlev>=3) WRITE(numout,*) 'Input data range LAT:', MINVAL(lat), MAXVAL(lat)

       CALL aggregate_p(nbpt, lalo, neighbours, resolution, contfrac, &
            &                iml, jml, lon, lat, mask, callsign, &
            &                nbvmax, sub_index, sub_area, ok_interpol)
       
       IF ( .NOT. ok_interpol ) THEN
          DEALLOCATE(sub_area)
          DEALLOCATE(sub_index)
          nbvmax = nbvmax * 2
       ENDIF
    ENDDO
    IF (printlev>=4) WRITE (numout,*) '  interpweight_3Dcont_tsteps: read/allocate OK'

    outvar3D = zero
! Providing the interpolated grid point

    CALL interpweight_provide_interpolation3D_tsteps(typefrac, nbpt, 0, 0, iml, jml, lml, tml,               &
      nbvmax, zero, invar4D, sub_area, sub_index, varmax, varmin, un, deux, lalo,                     &
      defaultvalue, defaultNOvalue, outvar3D, aoutvar)
    IF (printlev >= 5) WRITE(numout,*)'  interpweight_3Dcont_tsteps end of interpweight_provide_interpolation3D_tsteps'

    IF (printlev>=2) THEN
       IF (COUNT(aoutvar == -1) >= AINT(nbpt*0.4)) THEN
          WRITE(numout,*) 'interpweight_3Dcont_tsteps total number of points: ', nbpt
          WRITE(numout,*) '  Interpolated: ', COUNT(aoutvar /= -1),                 &
               ' non-interpolated: ', COUNT(aoutvar == -1)
          msg = "Less than 60% of points '" // TRIM(filename) // "' variable '" // TRIM(varname) //       &
               "' have been succesfully interpolated !!"
          WRITE(numout,*) '  ' // TRIM(msg)
       ELSE
          msg = TRIM(filename) // "' variable: '" // TRIM(varname) // "' interpolation Done _______"
          WRITE(numout,*) "interpweight_3Dcont_tsteps: '" // TRIM(msg)
          WRITE(numout,*) '  Total number of points: ', nbpt
          WRITE(numout,*) '  Interpolated: ', COUNT(aoutvar /= -1),                 &
               ' Non-interpolated: ', COUNT(aoutvar == -1)
       END IF
    END IF

! Scaling to 1 the quality variable
    DO ip=1, nbpt
      IF (aoutvar(ip) /= -1) aoutvar(ip) = aoutvar(ip)/(resolution(ip,1)*resolution(ip,2))/contfrac(ip)
    END DO

    IF (printlev>=5) WRITE(numout,*) '  interpweight_3Dcont_tsteps: Interpolation Done'

    IF (ALLOCATED(lat_in)) DEALLOCATE(lat_in)
    IF (ALLOCATED(lon_in)) DEALLOCATE(lon_in)
    DEALLOCATE(lon)
    DEALLOCATE(lat)
!    IF (ALLOCATED(invar2D)) DEALLOCATE(invar2D)
!    IF (ALLOCATED(invar3D)) DEALLOCATE(invar3D)
    IF (ALLOCATED(invar4D)) DEALLOCATE(invar4D)
    DEALLOCATE(mask)
    IF (ALLOCATED(maskvar)) DEALLOCATE(maskvar)
    IF (ALLOCATED(sub_area)) DEALLOCATE(sub_area)
    IF (ALLOCATED(sub_index)) DEALLOCATE(sub_index)
    IF (ALLOCATED(resol_in)) DEALLOCATE(resol_in)

    RETURN

  END SUBROUTINE interpweight_3Dcont_tstep

!! ================================================================================================================================
!! SUBROUTINE   : interpweight_2Dcont
!!
!>\BRIEF        ! Interpolate any input file to the grid of the model for a continuos 2D-variable  
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): outvar2D, aoutvar
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE interpweight_2Dcont(nbpt, dim1, dim2, lalo, resolution, neighbours,                      &
    contfrac, filename, varname, inlonname, inlatname, varmin, varmax, noneg, masktype,               &
    maskvalues, maskvarname, initime, typefrac, defaultvalue, defaultNOvalue,                         &
    outvar2D, aoutvar)

    USE ioipsl
    USE constantes_var
    USE mod_orchidee_para

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !
    !! 0.1 Input variables
    !
    INTEGER(i_std), INTENT(in)                 :: nbpt                    !! Number of points for which the data 
                                                                          !!   needs to be interpolated
    INTEGER(i_std), INTENT(in)                 :: dim1, dim2              !! 4D complementary dimension sizes
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in) :: lalo                    !! Vector of latitude and longitudes 
                                                                          !! (beware of the order = 1 : latitude, 
                                                                          !!   2 : longitude)
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in) :: resolution              !! The size in km of each grid-box in X and Y
    INTEGER(i_std), DIMENSION(nbpt,NbNeighb), INTENT(in) :: neighbours    !! Vector of neighbours for each grid point 
                                                                          !!   1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW) 
                            
    REAL(r_std), DIMENSION(nbpt), INTENT(in)   :: contfrac                !! Fraction of land in each grid box.
    CHARACTER(LEN=80), INTENT(in)              :: filename                !! name of the input map to read

    CHARACTER(LEN=80), INTENT(in)              :: varname                 !! name of the variable to interpolate
    CHARACTER(LEN=80), INTENT(in)              :: inlonname, inlatname    !! names of the longitude and latitude
    REAL(r_std), INTENT(in)                    :: varmin, varmax          !! min/max values to use for the
                                                                          !!   renormalization
    REAL(r_std), INTENT(in)                    :: defaultvalue            !! default interpolated value 
    REAL(r_std), INTENT(in)                    :: defaultNOvalue          !! default interpolated NO-value 

! Transformations to apply to the original read values from the file
    LOGICAL, INTENT(in)                        :: noneg                   !! no negative values
    CHARACTER(LEN=50), INTENT(in)              :: masktype                !! Type of masking
                                                                          !!   'nomask': no-mask is applied
                                                                          !!   'mbelow': take values below maskvalues(1)
                                                                          !!   'mabove': take values above maskvalues(1)
                                                                          !!   'msumrange': take values within 2 ranges;
                                                                          !!      maskvalues(2) <= SUM(vals(k)) <= maskvalues(1)
                                                                          !!      maskvalues(1) < SUM(vals(k)) <= maskvalues(3)
                                                                          !!        (normalized by maskedvalues(3))
                                                                          !!   'var': mask values are taken from a 
                                                                          !!     variable (1/0)
    REAL(r_std), DIMENSION(3), INTENT(in)      :: maskvalues              !! values to use to mask (according to 
                                                                          !!   `masktype') 
    CHARACTER(LEN=250), INTENT(in)             :: maskvarname             !! name of the variable to use to mask 
                                                                          !!   (according to `masktype') 
    CHARACTER(LEN=50), INTENT(in)              :: typefrac                !! Type of fraction retrieval:
                                                                          !!   'default': standard

    INTEGER(i_std), INTENT(in)                 :: initime                 !! Initial time to take (for input data 
                                                                          !!   with time values, -1: not used)
    !! 0.2 Modified variables
    !
    !! 0.3 Output variables
    !
    REAL(r_std), DIMENSION(nbpt), INTENT(out)  :: outvar2D                !! 2D output variable and re-dimensioned
    REAL(r_std), DIMENSION(nbpt), INTENT(out)  :: aoutvar                 !! availability of input data to
                                                                          !!   interpolate output variable 
                                                                          !!   (on the nbpt space)
    !
    !! 0.4 Local variables
    !
    INTEGER(i_std)                             :: iml, jml, lml, tml, fid, ib, ip, jp, it, jj, jv
    INTEGER                                    :: inNdims                 !! Number of dimensions of the input data
    INTEGER                                    :: inNLldims               !! Number of dimensions of the lon/lat input data
    REAL(r_std), ALLOCATABLE, DIMENSION(:)     :: lat_in, lon_in          !! latitude and
                                                                          !! longitude, extract from input map
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)   :: lat, lon                !! en 2D ???
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)   :: sub_area                !! the area of the fine grid in the 
                                                                          !!   model grid ???
                                                                          !! cf src_global/interpol_help.f90, 
                                                                          !!   line 377, called "areaoverlap"
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:,:)  :: sub_index           !! the indexes from the grid boxes from the 
                                                                          !!   data that go into the 
                                                                          !!   model's boxes  
                                                                          !! cf src_global/interpol_help.f90,
                                                                          !!   line 300, called "ip"

!    REAL(r_std), ALLOCATABLE, DIMENSION(:)       :: invar1D               !! 1D input values
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)     :: invar2D               !! 2D input values
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)   :: invar3D               !! 3D input values
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:,:) :: invar4D               !! 4D input values

    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)   :: resol_in
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:)  :: mask
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)     :: maskvar
    INTEGER(i_std)                             :: nix, njx
    INTEGER(i_std)                             :: idi, nbvmax             !! nbvmax : number of maximum vegetation map
                                                                          !! points in the GCM grid ; idi : its counter
    CHARACTER(LEN=30)                          :: callsign                !! Allows to specify which variable is beeing
                                                                          !! treated
    LOGICAL                                    :: ok_interpol             !! optionnal return of aggregate_2d
    INTEGER                                    :: ALLOC_ERR 

    CHARACTER(LEN=50)                          :: S1
    CHARACTER(LEN=250)                         :: msg


! Debug vars
    INTEGER                                    :: nb_coord, nb_var, nb_gat

!_ ================================================================================================================================
    IF (printlev >= 3) WRITE(numout,*)"In interpweight_2Dcont filename: ",TRIM(filename),          &
      " varname:'" // TRIM(varname) // "'"
! Allocating input data
    IF (is_root_prc) THEN
       IF (printlev >=5) THEN
          WRITE(numout,*) "Entering 'interpweight_2Dcont'. Debug mode."
          WRITE(numout,*)'(/," --> fliodmpf")'
          CALL fliodmpf (TRIM(filename))
          WRITE(numout,*)'(/," --> flioopfd")'
       ENDIF
       CALL flioopfd (TRIM(filename),fid,nb_dim=nb_coord,nb_var=nb_var,nb_gat=nb_gat)
       IF (printlev >=5) THEN
          WRITE (numout,'(" Number of coordinate        in the file : ",I2)') nb_coord
          WRITE (numout,'(" Number of variables         in the file : ",I2)') nb_var
          WRITE (numout,'(" Number of global attributes in the file : ",I2)') nb_gat
       ENDIF
    ENDIF

! Getting shape of input variable from input file
    inNdims = interpweight_get_varNdims_file(filename, varname)
    IF (is_root_prc) CALL flininfo(filename, iml, jml, lml, tml, fid)

    InputDims: SELECT CASE (inNdims)
      CASE (2)
        CALL bcast(iml)
        CALL bcast(jml)

        ALLOCATE(invar2D(iml,jml), STAT=ALLOC_ERR)
        IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2Dcont',"Problem in allocation of variable '" //     &
          TRIM(varname) ,'','')
        IF (is_root_prc) CALL flinget(fid, TRIM(varname), iml, jml, 0, 0, 1, 1, invar2D)
        IF (printlev >= 5) THEN
          WRITE(numout,*)'  interpweight_2Dcont input data dimensions _______'
          WRITE(numout,*)'    iml, jml:',iml, jml
        END IF

      CASE (3)
        CALL bcast(iml)
        CALL bcast(jml)
        CALL bcast(lml)

        ALLOCATE(invar3D(iml,jml, lml), STAT=ALLOC_ERR)
        IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2Dcont',"Problem in allocation of variable '" //     &
          TRIM(varname) ,'','')
        IF (is_root_prc) CALL flinget(fid, TRIM(varname), iml, jml, lml, 0, 1, 1, invar3D)
        IF (printlev >= 5) THEN
          WRITE(numout,*)'  interpweight_2Dcont input data dimensions:'
          WRITE(numout,*)'    iml, jml, lml:',iml, jml, lml
        END IF

      CASE (4)
        CALL bcast(iml)
        CALL bcast(jml)
        CALL bcast(lml)
!       Taken from `slowproc_update'
        IF (initime /= -1) THEN
          inNdims = 3
          IF (printlev >= 3) WRITE(numout,*) &
               'interpweight_2Dcont: taking into account the initial reading time:',initime
          IF (initime <= 0) tml = 0
          
          CALL bcast(tml)
          ALLOCATE(invar3D(iml,jml,lml), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2Dcont',"Problem in allocation of variable '" //   &
            TRIM(varname) ,'','')
          IF ( initime > 0 ) THEN
            IF (is_root_prc) THEN
              IF (initime <= tml) THEN
                CALL flinget(fid, TRIM(varname), iml, jml, lml, 1, initime, 1, invar3D)
              ELSE
                CALL flinget(fid, TRIM(varname), iml, jml, lml, 1, tml, 1, invar3D)
              ENDIF
            ENDIF
          ELSE
            IF (is_root_prc) THEN
              CALL flinget(fid, TRIM(varname), iml, jml, lml, 1, 1, 1, invar3D)
            ENDIF
          ENDIF
        ELSE
          CALL bcast(tml)
          ALLOCATE(invar4D(iml,jml,lml,tml), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2Dcont',"Problem in allocation of variable '" //   &
            TRIM(varname) ,'','')
          IF (is_root_prc) CALL flinget(fid, TRIM(varname), iml, jml, lml, tml, 1, tml, invar4D)
        END IF

        IF (printlev >= 5) THEN
          WRITE(numout,*)'  interpweight_2Dcont input data dimensions:'
          WRITE(numout,*)'    iml, jml, lml, tml :',iml, jml, lml, tml
        END IF

      CASE DEFAULT
        WRITE (S1,'(I3)') inNdims
        CALL ipslerr_p(3,'interpweight_2Dcont',"Input number of dimensions: '" // TRIM(S1) // "' not ready !!", &
          '','')

    END SELECT InputDims
    CALL bcast(invar2D)

    IF (printlev >= 5) WRITE(numout,*)'  interpweight_2Dcont getting variable ended!'

! Getting longitudes/latitudes from input file
    inNLldims = interpweight_get_varNdims_file(filename, inlonname)

    IF (inNLldims == 1) THEN
      ! Allocate memory for latitudes
      ALLOCATE(lon_in(iml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2Dcont','Problem in allocation of variable lon_in','','')

      ALLOCATE(lat_in(jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2Dcont','Problem in allocation of variable lat_in','','')

      IF (is_root_prc) THEN
        CALL flinget(fid, TRIM(inlonname), iml, 0, 0, 0, 1, 1, lon_in)
        CALL flinget(fid, TRIM(inlatname), jml, 0, 0, 0, 1, 1, lat_in)
      ENDIF
      CALL bcast(lon_in)
      CALL bcast(lat_in)

      ALLOCATE(lon(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2Dcont','Problem in allocation of variable lon','','')
      ALLOCATE(lat(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2Dcont','Problem in allocation of variable lat','','')

      DO ip=1,iml
         lat(ip,:) = lat_in(:)
      ENDDO
      DO jp=1,jml
         lon(:,jp) = lon_in(:)
      ENDDO
    ELSE
      ! Allocate memory for longitude
      ALLOCATE(lon(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2Dcont','Pb in allocation for lon','','')
      ! Allocate memory for latitudes
      ALLOCATE(lat(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2Dcont','Pb in allocation for lat','','')

      IF (is_root_prc) THEN
        CALL flinget(fid, TRIM(inlonname), iml, jml, 0, 0, 1, 1, lon)
        CALL flinget(fid, TRIM(inlatname), iml, jml, 0, 0, 1, 1, lat)
      ENDIF
    END IF
    CALL bcast(lon)
    CALL bcast(lat)

    IF (is_root_prc) THEN
      ALLOCATE(resol_in(iml,jml,2), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2Dcont','Problem in allocation of variable resol_in','','')
    END IF

    IF (noneg) THEN
      IF (is_root_prc) THEN
        SELECT CASE (inNdims)
          CASE(2)
            CALL interpweight_modifying_input2D(iml, jml, 0, 0, 1, 1, zero, invar2D)
          CASE(3)
            CALL interpweight_modifying_input3D(iml, jml, lml, 0, 1, 1, zero, invar3D)
          CASE(4)
            CALL interpweight_modifying_input4D(iml, jml, lml, tml, 1, tml, zero, invar4D)
        END SELECT       
      END IF
    END IF
    IF (printlev >= 5) WRITE(numout,*)'  interpweight_2Dcont modifying input ended!'

    !
    ! Consider all points a priori
    !
    ALLOCATE(mask(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2Dcont','Problem in allocation of variable mask','','')
    mask(:,:) = 0

    IF (TRIM(masktype) == 'var') THEN
      ALLOCATE(maskvar(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2Dcont','Problem in allocation of variable maskvar','','')
      maskvar(:,:) = 0
! Reads the mask from the file
      IF (is_root_prc) THEN
! ATTTENTION!! One can not use flinget in this way. It always spect to have a REAL variable!
        CALL flinget(fid, TRIM(maskvarname), iml, jml, 0, 0, 1, 1, maskvar)
        WHERE (maskvar > zero)
          mask = un 
        END WHERE
      END IF
      CALL bcast(mask)
    ELSE
      SELECT CASE (inNdims)
        CASE(2)
          CALL interpweight_masking_input2D(iml, jml, 0, 0, 0, 0, masktype, maskvalues, invar2D,      &
            mask)
        CASE(3)
          CALL interpweight_masking_input3D(iml, jml, 0, 0, lml, 0, masktype, maskvalues, invar3D,    &
            mask)
        CASE(4)
          CALL interpweight_masking_input4D(iml, jml, 0, 0, lml, tml, masktype, maskvalues, invar4D,  &
            mask)
      END SELECT
    END IF
    IF (printlev >= 5) WRITE(numout,*)'  interpweight_2Dcont masking input ended!'
    IF (is_root_prc) CALL flinclo(fid)

    IF (printlev>=5) WRITE(numout,*)"  interpweight_2Dcont '" // TRIM(typefrac) // "' interpolation"

! Assuming regular lon/lat projection for the input data
    IF (is_root_prc) THEN
      CALL interpweight_calc_resolution_in(lon, lat, iml, jml, mincos, R_Earth, pi, resol_in)
    END IF

!   Computation of nbvmax
!   In 'condveg_soilalb, slowproc_update' nbvmax=200

!
!   The number of maximum vegetation map points in the GCM grid is estimated.
!   Some lmargin is taken.
!
    IF (is_root_prc) THEN
       nix=INT(MAXVAL(resolution(:,1))/MAXVAL(resol_in(:,:,1)))+2
       njx=INT(MAXVAL(resolution(:,2))/MAXVAL(resol_in(:,:,2)))+2
       nbvmax = MAX(nix*njx, 200)
    ENDIF
    CALL bcast(nbvmax)
    CALL bcast(nix)
    CALL bcast(njx)

    callsign = TRIM(varname) // ' map'

    ok_interpol = .FALSE.
    DO WHILE ( .NOT. ok_interpol )
       IF (printlev>=3) WRITE(numout,*) "Projection arrays for ",callsign," : "
       IF (printlev>=3) WRITE(numout,*) "nbvmax = ",nbvmax, nix, njx

       ALLOCATE(sub_index(nbpt, nbvmax, 2), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2Dcont','Problem in allocation of variable sub_index','','')

       sub_index(:,:,:)=0

       ALLOCATE(sub_area(nbpt, nbvmax), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_2Dcont','Problem in allocation of variable sub_area','','')

       sub_area(:,:)=zero
       IF (printlev>=3) WRITE(numout,*) 'Input data range LON:', MINVAL(lon), MAXVAL(lon)
       IF (printlev>=3) WRITE(numout,*) 'Input data range LAT:', MINVAL(lat), MAXVAL(lat)

       CALL aggregate_p(nbpt, lalo, neighbours, resolution, contfrac, &
            &                iml, jml, lon, lat, mask, callsign, &
            &                nbvmax, sub_index, sub_area, ok_interpol)
       
       IF ( .NOT. ok_interpol ) THEN
          DEALLOCATE(sub_area)
          DEALLOCATE(sub_index)
          nbvmax = nbvmax * 2
       ENDIF
    ENDDO
    IF (printlev>=4) WRITE (numout,*) '  interpweight_2Dcont: read/allocate OK'

    outvar2D = zero
! Providing the interpolated grid point

    CALL interpweight_provide_interpolation2D(typefrac, nbpt, 0, 0, iml, jml, lml, tml,               &
      nbvmax, zero, invar2D, sub_area, sub_index, varmax, varmin, un, deux, lalo,                     &
      defaultvalue, defaultNOvalue, outvar2D, aoutvar)
    IF (printlev >= 5) WRITE(numout,*)'  interpweight_2Dcont end of interpweight_provide_interpolation2D'

    IF (printlev>=2) THEN
       IF (COUNT(aoutvar == -1) >= AINT(nbpt*0.4)) THEN
          WRITE(numout,*) 'interpweight_2Dcont total number of points: ', nbpt
          WRITE(numout,*) '  Interpolated: ', COUNT(aoutvar /= -1),                 &
               ' non-interpolated: ', COUNT(aoutvar == -1)
          msg = "Less than 60% of points '" // TRIM(filename) // "' variable '" // TRIM(varname) //       &
               "' have been succesfully interpolated !!"
          WRITE(numout,*) '  ' // TRIM(msg)
       ELSE
          msg = TRIM(filename) // "' variable: '" // TRIM(varname) // "' interpolation Done _______"
          WRITE(numout,*) "interpweight_2Dcont: '" // TRIM(msg)
          WRITE(numout,*) '  Total number of points: ', nbpt
          WRITE(numout,*) '  Interpolated: ', COUNT(aoutvar /= -1),                 &
               ' Non-interpolated: ', COUNT(aoutvar == -1)
       END IF
    END IF

! Scaling to 1 the quality variable
    DO ip=1, nbpt
      IF (aoutvar(ip) /= -1) aoutvar(ip) = aoutvar(ip)/(resolution(ip,1)*resolution(ip,2))/contfrac(ip)
    END DO

    IF (printlev>=5) WRITE(numout,*) '  interpweight_2Dcont: Interpolation Done'

    IF (ALLOCATED(lat_in)) DEALLOCATE(lat_in)
    IF (ALLOCATED(lon_in)) DEALLOCATE(lon_in)
    DEALLOCATE(lon)
    DEALLOCATE(lat)
    IF (ALLOCATED(invar2D)) DEALLOCATE(invar2D)
    IF (ALLOCATED(invar3D)) DEALLOCATE(invar3D)
    IF (ALLOCATED(invar4D)) DEALLOCATE(invar4D)
    DEALLOCATE(mask)
    IF (ALLOCATED(maskvar)) DEALLOCATE(maskvar)
    IF (ALLOCATED(invar4D)) DEALLOCATE(sub_area)
    IF (ALLOCATED(invar4D)) DEALLOCATE(sub_index)
    IF (ALLOCATED(resol_in)) DEALLOCATE(resol_in)

    RETURN

  END SUBROUTINE interpweight_2Dcont

!! ================================================================================================================================
!! SUBROUTINE   : interpweight_3Dcont
!!
!>\BRIEF        ! Interpolate any input file to the grid of the model for a continuos 4D-variable  
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): outvar4D, aoutvar
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE interpweight_3Dcont(nbpt, dim1, lalo, resolution, neighbours,                      &
    contfrac, filename, varname, inlonname, inlatname, varmin, varmax, noneg, masktype,               &
    maskvalues, maskvarname, initime, typefrac, defaultvalue, defaultNOvalue,                         &
    outvar3D, aoutvar)

    USE ioipsl
    USE constantes_var
    USE mod_orchidee_para

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !
    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)                 :: nbpt                    !! Number of points for which the data 
                                                                          !!   needs to be interpolated
    INTEGER(i_std), INTENT(in)                 :: dim1                    !! 4D complementary dimension sizes
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in) :: lalo                    !! Vector of latitude and longitudes 
                                                                          !! (beware of the order = 1 : latitude, 
                                                                          !!   2 : longitude)
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in) :: resolution              !! The size in km of each grid-box in X and Y
    INTEGER(i_std), DIMENSION(nbpt,NbNeighb), INTENT(in) :: neighbours    !! Vector of neighbours for each grid point 
                                                                          !!   1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW) 
                            
    REAL(r_std), DIMENSION(nbpt), INTENT(in)   :: contfrac                !! Fraction of land in each grid box.
    CHARACTER(LEN=*), INTENT(in)              :: filename                !! name of the input map to read

    CHARACTER(LEN=*), INTENT(in)              :: varname                 !! name of the variable to interpolate
    CHARACTER(LEN=*), INTENT(in)              :: inlonname, inlatname    !! names of the longitude and latitude
    REAL(r_std), INTENT(in)                    :: varmin, varmax          !! min/max values to use for the
                                                                          !!   renormalization
    REAL(r_std), DIMENSION(dim1), INTENT(in) :: defaultvalue         !! default interpolated value 
    REAL(r_std), INTENT(in)                    :: defaultNOvalue          !! default interpolated NO-value 

! Transformations to apply to the original read values from the file
    LOGICAL, INTENT(in)                        :: noneg                   !! no negative values
    CHARACTER(LEN=50), INTENT(in)              :: masktype                !! Type of masking
                                                                          !!   'nomask': no-mask is applied
                                                                          !!   'mbelow': take values below maskvalues(1)
                                                                          !!   'mabove': take values above maskvalues(1)
                                                                          !!   'msumrange': take values within 2 ranges;
                                                                          !!      maskvalues(2) <= SUM(vals(k)) <= maskvalues(1)
                                                                          !!      maskvalues(1) < SUM(vals(k)) <= maskvalues(3)
                                                                          !!        (normalized by maskedvalues(3))
                                                                          !!   'var': mask values are taken from a 
                                                                          !!     variable (1/0)
    REAL(r_std), DIMENSION(3), INTENT(in)      :: maskvalues              !! values to use to mask (according to 
                                                                          !!   `masktype') 
    CHARACTER(LEN=250), INTENT(in)             :: maskvarname             !! name of the variable to use to mask 
                                                                          !!   (according to `masktype') 
    CHARACTER(LEN=50), INTENT(in)              :: typefrac                !! Type of fraction retrieval:
                                                                          !!   'default': standard

    INTEGER(i_std), INTENT(in)                 :: initime                 !! Initial time to take (for input data 
                                                                          !!   with time values, -1: not used)
    !! 0.2 Modified variables
    !
    !! 0.3 Output variables
    !
    REAL(r_std), DIMENSION(nbpt,dim1), INTENT(out) :: outvar3D       !! 2D output variable and re-dimensioned
    REAL(r_std), DIMENSION(nbpt), INTENT(out)  :: aoutvar                 !! availability of input data to
                                                                          !!   interpolate output variable 
                                                                          !!   (on the nbpt space)
    !
    !! 0.4 Local variables
    !
    INTEGER(i_std)                             :: iml, jml, lml, tml, fid, ib, ip, jp, it, jj, jv
    INTEGER                                    :: inNdims                 !! Number of dimensions of the input data
    INTEGER                                    :: inNLldims               !! Number of dimensions of the lon/lat input data
    REAL(r_std), ALLOCATABLE, DIMENSION(:)     :: lat_in, lon_in          !! latitude and
                                                                          !! longitude, extract from input map
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)   :: lat, lon                !! en 2D ???
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)   :: sub_area                !! the area of the fine grid in the 
                                                                          !!   model grid ???
                                                                          !! cf src_global/interpol_help.f90, 
                                                                          !!   line 377, called "areaoverlap"
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:,:)  :: sub_index           !! the indexes from the grid boxes from the 
                                                                          !!   data that go into the 
                                                                          !!   model's boxes  
                                                                          !! cf src_global/interpol_help.f90,
                                                                          !!   line 300, called "ip"

!    REAL(r_std), ALLOCATABLE, DIMENSION(:)       :: invar1D               !! 1D input values
!    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)     :: invar2D               !! 2D input values
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)   :: invar3D               !! 3D input values
!    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:,:) :: invar4D               !! 4D input values

    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)   :: resol_in
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)     :: maskvar
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:)  :: mask
    INTEGER(i_std)                             :: nix, njx
    INTEGER(i_std)                             :: idi, nbvmax             !! nbvmax : number of maximum vegetation map
                                                                          !! points in the GCM grid ; idi : its counter
    CHARACTER(LEN=30)                          :: callsign                !! Allows to specify which variable is beeing
                                                                          !! treated
    LOGICAL                                    :: ok_interpol             !! optionnal return of aggregate_2d
    INTEGER                                    :: ALLOC_ERR 

    CHARACTER(LEN=50)                          :: S1
    INTEGER, ALLOCATABLE, DIMENSION(:)         :: indims
    CHARACTER(LEN=250)                         :: msg

! Debug vars
    INTEGER                                    :: nb_coord, nb_var, nb_gat
    INTEGER                                    :: nblL
    INTEGER, DIMENSION(2)                      :: lLpt

!_ ================================================================================================================================
    IF (printlev >= 3) WRITE(numout, *)"In interpweight_3Dcont filename: '",TRIM(filename),         &
      " varname:'" // TRIM(varname) // "'"
! Allocating input data
    IF (is_root_prc) THEN
       IF (printlev >=5) THEN
          WRITE(numout,*) "Entering 'interpweight_3Dcont'. Debug mode."
          WRITE(numout,*)'(/," --> fliodmpf")'
          CALL fliodmpf (TRIM(filename))
          WRITE(numout,*)'(/," --> flioopfd")'
       ENDIF
       CALL flioopfd (TRIM(filename),fid,nb_dim=nb_coord,nb_var=nb_var,nb_gat=nb_gat)
       IF (printlev >=5) THEN
          WRITE (numout,'(" Number of coordinate        in the file : ",I2)') nb_coord
          WRITE (numout,'(" Number of variables         in the file : ",I2)') nb_var
          WRITE (numout,'(" Number of global attributes in the file : ",I2)') nb_gat
       ENDIF
    ENDIF

! Getting shape of input variable from input file
    inNdims = interpweight_get_varNdims_file(filename, varname)
    IF (is_root_prc) CALL flininfo(filename, iml, jml, lml, tml, fid)

    IF (inNdims .NE. 3) CALL ipslerr_p(3,'interpweight_3Dcont','Expected 3 dimenions in the file '//filename, &
                                         'but found:',inNdims)

    ALLOCATE(indims(3),STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_3Dcont','Problem in allocation of variable indims','','')
    indims = interpweight_get_var3dims_file(filename, varname)
    lml = indims(3)
    CALL bcast(iml)
    CALL bcast(jml)
    CALL bcast(lml)

    ALLOCATE(invar3D(iml,jml, lml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_3Dcont',"Problem in allocation of variable '" //     &
      TRIM(varname) ,'','')
    IF (is_root_prc) CALL flinget(fid, TRIM(varname), iml, jml, lml, 0, 1, 1, invar3D)
    IF (printlev >= 5) THEN
      WRITE(numout,*)'  interpweight_3Dcont input data dimensions:'
      WRITE(numout,*)'    iml, jml, lml:',iml, jml, lml
    END IF

    CALL bcast(invar3D)

    IF (printlev >= 5) WRITE(numout,*)'  interpweight_3Dcont getting variable ended!'

! Getting longitudes/latitudes from input file
    inNLldims = interpweight_get_varNdims_file(filename, inlonname)

    IF (inNLldims == 1) THEN
      ! Allocate memory for latitudes
      ALLOCATE(lon_in(iml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_3Dcont','Problem in allocation of variable lon_in','','')

      ALLOCATE(lat_in(jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_3Dcont','Problem in allocation of variable lat_in','','')

      IF (is_root_prc) THEN
        CALL flinget(fid, TRIM(inlonname), iml, 0, 0, 0, 1, 1, lon_in)
        CALL flinget(fid, TRIM(inlatname), jml, 0, 0, 0, 1, 1, lat_in)
      ENDIF
      CALL bcast(lon_in)
      CALL bcast(lat_in)

      ALLOCATE(lon(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4Dcont','Problem in allocation of variable lon','','')
      ALLOCATE(lat(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4Dcont','Problem in allocation of variable lat','','')

      DO ip=1,iml
         lat(ip,:) = lat_in(:)
      ENDDO
      DO jp=1,jml
         lon(:,jp) = lon_in(:)
      ENDDO
    ELSE
      ! Allocate memory for longitude
      ALLOCATE(lon(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4Dcont','Pb in allocation for lon','','')
      ! Allocate memory for latitudes
      ALLOCATE(lat(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4Dcont','Pb in allocation for lat','','')

      IF (is_root_prc) THEN
        CALL flinget(fid, TRIM(inlonname), iml, jml, 0, 0, 1, 1, lon)
        CALL flinget(fid, TRIM(inlatname), iml, jml, 0, 0, 1, 1, lat)
      ENDIF
    END IF
    CALL bcast(lon)
    CALL bcast(lat)

    IF (is_root_prc) THEN
      ALLOCATE(resol_in(iml,jml,2), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4Dcont','Problem in allocation of variable resol_in','','')
    END IF

    IF (noneg) THEN
      IF (is_root_prc) THEN
          CALL interpweight_modifying_input3D(iml, jml, lml, 0, 0, 0, zero, invar3D)
      END IF
    END IF
    CALL bcast(invar3D)
    IF (printlev >= 5) WRITE(numout,*)'  interpweight_4Dcont modifying input ended!'

    !
    ! Consider all points a priori
    !
    ALLOCATE(mask(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4Dcont','Problem in allocation of variable mask','','')
    mask(:,:) = zero

    IF (TRIM(masktype) == 'var') THEN
      ALLOCATE(maskvar(iml,jml), STAT=ALLOC_ERR)
      maskvar(:,:) = zero
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4Dcont','Problem in allocation of variable maskvar','','')
      ! Reads the mask from the file
      IF (is_root_prc) THEN
         ! ATTTENTION!! One can not use flinget in this way. It always spect to have a REAL variable!
        CALL flinget(fid, TRIM(maskvarname), iml, jml, 0, 0, 1, 1, maskvar)
        WHERE (maskvar > zero)
          mask = un 
        ENDWHERE
      END IF
      CALL bcast(mask)
    ELSE
      CALL interpweight_masking_input3D(iml, jml, 0, 0, lml, 0, masktype, maskvalues, invar3D,    &
          mask)
    END IF

    IF (is_root_prc) CALL flinclo(fid)

! Assuming regular lon/lat projection for the input data
    IF (is_root_prc) THEN
      CALL interpweight_calc_resolution_in(lon, lat, iml, jml, mincos, R_Earth, pi, resol_in)
    END IF

    ! Computation of nbvmax
    ! In 'condveg_soilalb, slowproc_update' nbvmax=200

    !
    ! The number of maximum vegetation map points in the GCM grid is estimated.
    ! Some lmargin is taken.
    !
    IF (is_root_prc) THEN
       nix=INT(MAXVAL(resolution(:,1))/MAXVAL(resol_in(:,:,1)))+2
       njx=INT(MAXVAL(resolution(:,2))/MAXVAL(resol_in(:,:,2)))+2
       nbvmax = MAX(nix*njx, 200)
    ENDIF
    CALL bcast(nbvmax)
    CALL bcast(nix)
    CALL bcast(njx)

    !
    callsign = TRIM(varname) // ' map'
    !
    ok_interpol = .FALSE.
    DO WHILE ( .NOT. ok_interpol )
       IF (printlev>=3) WRITE(numout,*) "Projection arrays for ",callsign," : "
       IF (printlev>=3) WRITE(numout,*) "nbvmax = ",nbvmax, nix, njx


       ALLOCATE(sub_index(nbpt, nbvmax, 2), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_3Dcont','Problem in allocation of variable sub_index','','')

       sub_index(:,:,:)=0

       ALLOCATE(sub_area(nbpt, nbvmax), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_3Dcont','Problem in allocation of variable sub_area','','')

       sub_area(:,:)=zero
       IF (printlev>=3) WRITE(numout,*) 'Input data range LON:', MINVAL(lon), MAXVAL(lon)
       IF (printlev>=3) WRITE(numout,*) 'Input data range LAT:', MINVAL(lat), MAXVAL(lat)

       CALL aggregate_p(nbpt, lalo, neighbours, resolution, contfrac, &
            &                iml, jml, lon, lat, mask, callsign, &
            &                nbvmax, sub_index, sub_area, ok_interpol)
       
       IF ( .NOT. ok_interpol ) THEN
          DEALLOCATE(sub_area)
          DEALLOCATE(sub_index)
          nbvmax = nbvmax * 2
       ENDIF
    ENDDO
    IF (printlev>=4) WRITE (numout,*) '  interpweight_3Dcont: read/allocate OK'

    outvar3D = zero

! Providing the interpolated grid point
    CALL interpweight_provide_interpolation3D(typefrac, nbpt, dim1, 0, iml, jml, lml, tml,         &
      nbvmax, zero, invar3D, sub_area, sub_index, varmax, varmin, un, deux, lalo,                     &
      defaultvalue, defaultNOvalue, outvar3D, aoutvar)
    IF (printlev >= 5) WRITE(numout,*)'  interpweight_3Dcont end of interpweight_provide_interpolation3D'
    
    IF (printlev>=2) THEN
       IF (COUNT(aoutvar == -1) >= AINT(nbpt*0.4)) THEN
          WRITE(numout,*) 'interpweight_3Dcont total number of points: ', nbpt
          WRITE(numout,*) '  Interpolated: ', COUNT(aoutvar /= -1),                 &
               ' non-interpolated: ', COUNT(aoutvar == -1)
          msg = "Less than 60% of points '" // TRIM(filename) // "' variable '" // TRIM(varname) //       &
               "' have been succesfully interpolated !!"
          WRITE(numout,*) '  ' // TRIM(msg)
       ELSE
          msg = TRIM(filename) // "' variable: '" // TRIM(varname) // "' interpolation Done _______"
          WRITE(numout,*) "interpweight_3Dcont: '" // TRIM(msg)
          WRITE(numout,*) '  Total number of points: ', nbpt
          WRITE(numout,*) '  Interpolated: ', COUNT(aoutvar /= -1),                 &
               ' Non-interpolated: ', COUNT(aoutvar == -1)
       END IF
    END IF

! Scaling to 1 the quality variable
    DO ip=1, nbpt
      IF (aoutvar(ip) /= -1) aoutvar(ip) = aoutvar(ip)/(resolution(ip,1)*resolution(ip,2))/contfrac(ip)
    END DO

    IF (printlev>=5) WRITE(numout,*) '  interpweight_3Dcont: Interpolation Done'

    IF (ALLOCATED(lat_in)) DEALLOCATE(lat_in)
    IF (ALLOCATED(lon_in)) DEALLOCATE(lon_in)
    DEALLOCATE(lon)
    DEALLOCATE(lat)
    IF (ALLOCATED(invar3D)) DEALLOCATE(invar3D)
    DEALLOCATE(mask)
    IF (ALLOCATED(maskvar)) DEALLOCATE(maskvar)
    DEALLOCATE(sub_area)
    DEALLOCATE(sub_index)
    IF (ALLOCATED(resol_in)) DEALLOCATE(resol_in)

    RETURN

  END SUBROUTINE interpweight_3Dcont

!! ================================================================================================================================
!! SUBROUTINE   : interpweight_4Dcont
!!
!>\BRIEF        ! Interpolate any input file to the grid of the model for a continuos 4D-variable  
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): outvar4D, aoutvar
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE interpweight_4Dcont(nbpt, dim1, dim2, lalo, resolution, neighbours,                      &
    contfrac, filename, varname, inlonname, inlatname, varmin, varmax, noneg, masktype,               &
    maskvalues, maskvarname, initime, typefrac, defaultvalue, defaultNOvalue,                         &
    outvar4D, aoutvar)

    USE ioipsl
    USE constantes_var
    USE mod_orchidee_para

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !
    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)                 :: nbpt                    !! Number of points for which the data 
                                                                          !!   needs to be interpolated
    INTEGER(i_std), INTENT(in)                 :: dim1, dim2              !! 4D complementary dimension sizes
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in) :: lalo                    !! Vector of latitude and longitudes 
                                                                          !! (beware of the order = 1 : latitude, 
                                                                          !!   2 : longitude)
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in) :: resolution              !! The size in km of each grid-box in X and Y
    INTEGER(i_std), DIMENSION(nbpt,NbNeighb), INTENT(in) :: neighbours    !! Vector of neighbours for each grid point 
                                                                          !!   1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW) 
                            
    REAL(r_std), DIMENSION(nbpt), INTENT(in)   :: contfrac                !! Fraction of land in each grid box.
    CHARACTER(LEN=80), INTENT(in)              :: filename                !! name of the input map to read

    CHARACTER(LEN=80), INTENT(in)              :: varname                 !! name of the variable to interpolate
    CHARACTER(LEN=80), INTENT(in)              :: inlonname, inlatname    !! names of the longitude and latitude
    REAL(r_std), INTENT(in)                    :: varmin, varmax          !! min/max values to use for the
                                                                          !!   renormalization
    REAL(r_std), DIMENSION(dim1,dim2), INTENT(in) :: defaultvalue         !! default interpolated value 
    REAL(r_std), INTENT(in)                    :: defaultNOvalue          !! default interpolated NO-value 

! Transformations to apply to the original read values from the file
    LOGICAL, INTENT(in)                        :: noneg                   !! no negative values
    CHARACTER(LEN=50), INTENT(in)              :: masktype                !! Type of masking
                                                                          !!   'nomask': no-mask is applied
                                                                          !!   'mbelow': take values below maskvalues(1)
                                                                          !!   'mabove': take values above maskvalues(1)
                                                                          !!   'msumrange': take values within 2 ranges;
                                                                          !!      maskvalues(2) <= SUM(vals(k)) <= maskvalues(1)
                                                                          !!      maskvalues(1) < SUM(vals(k)) <= maskvalues(3)
                                                                          !!        (normalized by maskedvalues(3))
                                                                          !!   'var': mask values are taken from a 
                                                                          !!     variable (1/0)
    REAL(r_std), DIMENSION(3), INTENT(in)      :: maskvalues              !! values to use to mask (according to 
                                                                          !!   `masktype') 
    CHARACTER(LEN=250), INTENT(in)             :: maskvarname             !! name of the variable to use to mask 
                                                                          !!   (according to `masktype') 
    CHARACTER(LEN=50), INTENT(in)              :: typefrac                !! Type of fraction retrieval:
                                                                          !!   'default': standard

    INTEGER(i_std), INTENT(in)                 :: initime                 !! Initial time to take (for input data 
                                                                          !!   with time values, -1: not used)
    !! 0.2 Modified variables
    !
    !! 0.3 Output variables
    !
    REAL(r_std), DIMENSION(nbpt,dim1,dim2), INTENT(out) :: outvar4D       !! 2D output variable and re-dimensioned
    REAL(r_std), DIMENSION(nbpt), INTENT(out)  :: aoutvar                 !! availability of input data to
                                                                          !!   interpolate output variable 
                                                                          !!   (on the nbpt space)
    !
    !! 0.4 Local variables
    !
    INTEGER(i_std)                             :: iml, jml, lml, tml, fid, ib, ip, jp, it, jj, jv
    INTEGER                                    :: inNdims                 !! Number of dimensions of the input data
    INTEGER                                    :: inNLldims               !! Number of dimensions of the lon/lat input data
    REAL(r_std), ALLOCATABLE, DIMENSION(:)     :: lat_in, lon_in          !! latitude and
                                                                          !! longitude, extract from input map
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)   :: lat, lon                !! en 2D ???
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)   :: sub_area                !! the area of the fine grid in the 
                                                                          !!   model grid ???
                                                                          !! cf src_global/interpol_help.f90, 
                                                                          !!   line 377, called "areaoverlap"
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:,:)  :: sub_index           !! the indexes from the grid boxes from the 
                                                                          !!   data that go into the 
                                                                          !!   model's boxes  
                                                                          !! cf src_global/interpol_help.f90,
                                                                          !!   line 300, called "ip"

!    REAL(r_std), ALLOCATABLE, DIMENSION(:)       :: invar1D               !! 1D input values
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)     :: invar2D               !! 2D input values
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)   :: invar3D               !! 3D input values
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:,:) :: invar4D               !! 4D input values

    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)   :: resol_in
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)     :: maskvar
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:)  :: mask
    INTEGER(i_std)                             :: nix, njx
    INTEGER(i_std)                             :: idi, nbvmax             !! nbvmax : number of maximum vegetation map
                                                                          !! points in the GCM grid ; idi : its counter
    CHARACTER(LEN=30)                          :: callsign                !! Allows to specify which variable is beeing
                                                                          !! treated
    LOGICAL                                    :: ok_interpol             !! optionnal return of aggregate_2d
    INTEGER                                    :: ALLOC_ERR 

    CHARACTER(LEN=50)                          :: S1
    INTEGER, ALLOCATABLE, DIMENSION(:)         :: indims
    CHARACTER(LEN=250)                         :: msg

! Debug vars
    INTEGER                                    :: nb_coord, nb_var, nb_gat
    INTEGER                                    :: nblL
    INTEGER, DIMENSION(2)                      :: lLpt

!_ ================================================================================================================================
    IF (printlev >= 3) WRITE(numout, *)"In interpweight_4Dcont filename: '",TRIM(filename),         &
      " varname:'" // TRIM(varname) // "'"
! Allocating input data
    IF (is_root_prc) THEN
       IF (printlev >=5) THEN
          WRITE(numout,*) "Entering 'interpweight_4Dcont'. Debug mode."
          WRITE(numout,*)'(/," --> fliodmpf")'
          CALL fliodmpf (TRIM(filename))
          WRITE(numout,*)'(/," --> flioopfd")'
       ENDIF
       CALL flioopfd (TRIM(filename),fid,nb_dim=nb_coord,nb_var=nb_var,nb_gat=nb_gat)
       IF (printlev >=5) THEN
          WRITE (numout,'(" Number of coordinate        in the file : ",I2)') nb_coord
          WRITE (numout,'(" Number of variables         in the file : ",I2)') nb_var
          WRITE (numout,'(" Number of global attributes in the file : ",I2)') nb_gat
       ENDIF
    ENDIF

! Getting shape of input variable from input file
    inNdims = interpweight_get_varNdims_file(filename, varname)
    IF (is_root_prc) CALL flininfo(filename, iml, jml, lml, tml, fid)

    InputDims: SELECT CASE (inNdims)
      CASE (2)
        CALL bcast(iml)
        CALL bcast(jml)

        ALLOCATE(invar2D(iml,jml), STAT=ALLOC_ERR)
        IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4Dcont',"Problem in allocation of variable '" //     &
          TRIM(varname) ,'','')
        IF (is_root_prc) CALL flinget(fid, TRIM(varname), iml, jml, 0, 0, 1, 1, invar2D)
        IF (printlev >= 5) THEN
          WRITE(numout,*)'  interpweight_4Dcont input data dimensions:'
          WRITE(numout,*)'    iml, jml:',iml, jml
        END IF

      CASE (3)
        ALLOCATE(indims(3),STAT=ALLOC_ERR)
        IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4Dcont','Problem in allocation of variable indims','','')
        indims = interpweight_get_var3dims_file(filename, varname)
        lml = indims(3)
        CALL bcast(iml)
        CALL bcast(jml)
        CALL bcast(lml)

        ALLOCATE(invar3D(iml,jml, lml), STAT=ALLOC_ERR)
        IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4Dcont',"Problem in allocation of variable '" //     &
          TRIM(varname) ,'','')
        IF (is_root_prc) CALL flinget(fid, TRIM(varname), iml, jml, lml, 0, 1, 1, invar3D)
        IF (printlev >= 5) THEN
          WRITE(numout,*)'  interpweight_4Dcont input data dimensions:'
          WRITE(numout,*)'    iml, jml, lml:',iml, jml, lml
        END IF

      CASE (4)
        ALLOCATE(indims(4), STAT=ALLOC_ERR)
        IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4Dcont','Problem in allocation of indims','','')
        indims = interpweight_get_var4dims_file(filename, varname)
        lml = indims(3)
        CALL bcast(iml)
        CALL bcast(jml)
        CALL bcast(lml)
!       Taken from `slowproc_update'
        IF (initime /= -1) THEN
          inNdims = 3
          IF (printlev >= 3) WRITE(numout,*) &
               'interpweight_4Dcont: taking into account the initial reading time:',initime
          IF (initime <= 0) tml = 0

          CALL bcast(tml)
          ALLOCATE(invar3D(iml,jml,lml), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4Dcont',"Problem in allocation of variable '" //   &
            TRIM(varname) ,'','')
          IF ( initime > 0 ) THEN
            IF (is_root_prc) THEN
              IF (initime <= tml) THEN
                CALL flinget(fid, TRIM(varname), iml, jml, lml, 1, initime, 1, invar3D)
              ELSE
                CALL flinget(fid, TRIM(varname), iml, jml, lml, 1, tml, 1, invar3D)
              ENDIF
            ENDIF
          ELSE
            IF (is_root_prc) THEN
              CALL flinget(fid, TRIM(varname), iml, jml, lml, 1, 1, 1, invar3D)
            ENDIF
          ENDIF
        ELSE
          CALL bcast(tml)
          ALLOCATE(invar4D(iml,jml,lml,tml), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4Dcont',"Problem in allocation of variable '" //   &
            TRIM(varname) ,'','')
          IF (is_root_prc) CALL flinget(fid, TRIM(varname), iml, jml, lml, tml, 1, tml, invar4D)
        END IF

        IF (printlev >= 5) THEN
          WRITE(numout,*)'interpweight_4Dcont input data dimensions:'
          WRITE(numout,*)'    iml, jml, lml, tml:',iml, jml, lml, tml
        END IF

      CASE DEFAULT
        WRITE (S1,'(I3)') inNdims
        CALL ipslerr_p(3,'interpweight_4Dcont',"Input number of dimensions: '" // S1 // "' not ready !!",'','')

    END SELECT InputDims
    CALL bcast(invar4D)

    IF (printlev >= 5) WRITE(numout,*)'  interpweight_4Dcont getting variable ended!'

! Getting longitudes/latitudes from input file
    inNLldims = interpweight_get_varNdims_file(filename, inlonname)

    IF (inNLldims == 1) THEN
      ! Allocate memory for latitudes
      ALLOCATE(lon_in(iml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4Dcont','Problem in allocation of variable lon_in','','')

      ALLOCATE(lat_in(jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4Dcont','Problem in allocation of variable lat_in','','')

      IF (is_root_prc) THEN
        CALL flinget(fid, TRIM(inlonname), iml, 0, 0, 0, 1, 1, lon_in)
        CALL flinget(fid, TRIM(inlatname), jml, 0, 0, 0, 1, 1, lat_in)
      ENDIF
      CALL bcast(lon_in)
      CALL bcast(lat_in)

      ALLOCATE(lon(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4Dcont','Problem in allocation of variable lon','','')
      ALLOCATE(lat(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4Dcont','Problem in allocation of variable lat','','')

      DO ip=1,iml
         lat(ip,:) = lat_in(:)
      ENDDO
      DO jp=1,jml
         lon(:,jp) = lon_in(:)
      ENDDO
    ELSE
      ! Allocate memory for longitude
      ALLOCATE(lon(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4Dcont','Pb in allocation for lon','','')
      ! Allocate memory for latitudes
      ALLOCATE(lat(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4Dcont','Pb in allocation for lat','','')

      IF (is_root_prc) THEN
        CALL flinget(fid, TRIM(inlonname), iml, jml, 0, 0, 1, 1, lon)
        CALL flinget(fid, TRIM(inlatname), iml, jml, 0, 0, 1, 1, lat)
      ENDIF
    END IF
    CALL bcast(lon)
    CALL bcast(lat)

    IF (is_root_prc) THEN
      ALLOCATE(resol_in(iml,jml,2), STAT=ALLOC_ERR)
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4Dcont','Problem in allocation of variable resol_in','','')
    END IF

    IF (noneg) THEN
      IF (is_root_prc) THEN
        SELECT CASE (inNdims)
          CASE(2)
            CALL interpweight_modifying_input2D(iml, jml, 0, 0, 0, 0, zero, invar2D)
          CASE(3)
            CALL interpweight_modifying_input3D(iml, jml, lml, 0, 0, 0, zero, invar3D)
          CASE(4)
            CALL interpweight_modifying_input4D(iml, jml, lml, tml, 1, tml, zero, invar4D)
        END SELECT       
      END IF
    END IF
    CALL bcast(invar4D)
    IF (printlev >= 5) WRITE(numout,*)'  interpweight_4Dcont modifying input ended!'

    !
    ! Consider all points a priori
    !
    ALLOCATE(mask(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4Dcont','Problem in allocation of variable mask','','')
    mask(:,:) = zero

    IF (TRIM(masktype) == 'var') THEN
      ALLOCATE(maskvar(iml,jml), STAT=ALLOC_ERR)
      maskvar(:,:) = zero
      IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4Dcont','Problem in allocation of variable maskvar','','')
      ! Reads the mask from the file
      IF (is_root_prc) THEN
         ! ATTTENTION!! One can not use flinget in this way. It always spect to have a REAL variable!
        CALL flinget(fid, TRIM(maskvarname), iml, jml, 0, 0, 1, 1, maskvar)
        WHERE (maskvar > zero)
          mask = un 
        ENDWHERE
      END IF
      CALL bcast(mask)
    ELSE
      SELECT CASE (inNdims)
        CASE(2)
          CALL interpweight_masking_input2D(iml, jml, 0, 0, 0, 0, masktype, maskvalues, invar2D,      &
            mask)
        CASE(3)
          CALL interpweight_masking_input3D(iml, jml, 0, 0, lml, 0, masktype, maskvalues, invar3D,    &
            mask)
        CASE(4)
          CALL interpweight_masking_input4D(iml, jml, 0, 0, lml, tml, masktype, maskvalues, invar4D,  &
            mask)
      END SELECT
    END IF

    IF (is_root_prc) CALL flinclo(fid)

! Assuming regular lon/lat projection for the input data
    IF (is_root_prc) THEN
      CALL interpweight_calc_resolution_in(lon, lat, iml, jml, mincos, R_Earth, pi, resol_in)
    END IF

    ! Computation of nbvmax
    ! In 'condveg_soilalb, slowproc_update' nbvmax=200

    !
    ! The number of maximum vegetation map points in the GCM grid is estimated.
    ! Some lmargin is taken.
    !
    IF (is_root_prc) THEN
       nix=INT(MAXVAL(resolution(:,1))/MAXVAL(resol_in(:,:,1)))+2
       njx=INT(MAXVAL(resolution(:,2))/MAXVAL(resol_in(:,:,2)))+2
       nbvmax = MAX(nix*njx, 200)
    ENDIF
    CALL bcast(nbvmax)
    CALL bcast(nix)
    CALL bcast(njx)

    !
    callsign = TRIM(varname) // ' map'
    !
    ok_interpol = .FALSE.
    DO WHILE ( .NOT. ok_interpol )
       IF (printlev>=3) WRITE(numout,*) "Projection arrays for ",callsign," : "
       IF (printlev>=3) WRITE(numout,*) "nbvmax = ",nbvmax, nix, njx


       ALLOCATE(sub_index(nbpt, nbvmax, 2), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4Dcont','Problem in allocation of variable sub_index','','')

       sub_index(:,:,:)=0

       ALLOCATE(sub_area(nbpt, nbvmax), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_4Dcont','Problem in allocation of variable sub_area','','')

       sub_area(:,:)=zero
       IF (printlev>=3) WRITE(numout,*) 'Input data range LON:', MINVAL(lon), MAXVAL(lon)
       IF (printlev>=3) WRITE(numout,*) 'Input data range LAT:', MINVAL(lat), MAXVAL(lat)

       CALL aggregate_p(nbpt, lalo, neighbours, resolution, contfrac, &
            &                iml, jml, lon, lat, mask, callsign, &
            &                nbvmax, sub_index, sub_area, ok_interpol)
       
       IF ( .NOT. ok_interpol ) THEN
          DEALLOCATE(sub_area)
          DEALLOCATE(sub_index)
          nbvmax = nbvmax * 2
       ENDIF
    ENDDO
    IF (printlev>=4) WRITE (numout,*) '  interpweight_4Dcont: read/allocate OK'

    outvar4D = zero

! Providing the interpolated grid point
    CALL interpweight_provide_interpolation4D(typefrac, nbpt, dim1, dim2, iml, jml, lml, tml,         &
      nbvmax, zero, invar4D, sub_area, sub_index, varmax, varmin, un, deux, lalo,                     &
      defaultvalue, defaultNOvalue, outvar4D, aoutvar)
    IF (printlev >= 5) WRITE(numout,*)'  interpweight_4Dcont end of interpweight_provide_interpolation4D'
    
    IF (printlev>=2) THEN
       IF (COUNT(aoutvar == -1) >= AINT(nbpt*0.4)) THEN
          WRITE(numout,*) 'interpweight_4Dcont total number of points: ', nbpt
          WRITE(numout,*) '  Interpolated: ', COUNT(aoutvar /= -1),                 &
               ' non-interpolated: ', COUNT(aoutvar == -1)
          msg = "Less than 60% of points '" // TRIM(filename) // "' variable '" // TRIM(varname) //       &
               "' have been succesfully interpolated !!"
          WRITE(numout,*) '  ' // TRIM(msg)
       ELSE
          msg = TRIM(filename) // "' variable: '" // TRIM(varname) // "' interpolation Done _______"
          WRITE(numout,*) "interpweight_4Dcont: '" // TRIM(msg)
          WRITE(numout,*) '  Total number of points: ', nbpt
          WRITE(numout,*) '  Interpolated: ', COUNT(aoutvar /= -1),                 &
               ' Non-interpolated: ', COUNT(aoutvar == -1)
       END IF
    END IF

! Scaling to 1 the quality variable
    DO ip=1, nbpt
      IF (aoutvar(ip) /= -1) aoutvar(ip) = aoutvar(ip)/(resolution(ip,1)*resolution(ip,2))/contfrac(ip)
    END DO

    IF (printlev>=5) WRITE(numout,*) '  interpweight_4Dcont: Interpolation Done'

    IF (ALLOCATED(lat_in)) DEALLOCATE(lat_in)
    IF (ALLOCATED(lon_in)) DEALLOCATE(lon_in)
    DEALLOCATE(lon)
    DEALLOCATE(lat)
    IF (ALLOCATED(invar2D)) DEALLOCATE(invar2D)
    IF (ALLOCATED(invar3D)) DEALLOCATE(invar3D)
    IF (ALLOCATED(invar4D)) DEALLOCATE(invar4D)
    DEALLOCATE(mask)
    IF (ALLOCATED(maskvar)) DEALLOCATE(maskvar)
    DEALLOCATE(sub_area)
    DEALLOCATE(sub_index)
    IF (ALLOCATED(resol_in)) DEALLOCATE(resol_in)

    RETURN

  END SUBROUTINE interpweight_4Dcont

!! ================================================================================================================================
!! SUBROUTINE   : interpweight_calc_resolution_in
!!
!>\BRIEF        ! Subroutine to compute the resolution of the input data
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): resolin
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE interpweight_calc_resolution_in(lons, lats, dx, dy, mcos, REarth, piv, resolin)

    USE constantes_var

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !
    !! 0.1 Input variables
    !
    INTEGER(i_std), INTENT(in)                           :: dx, dy        !! Resolution input data
    REAL(r_std), INTENT(in)                              :: mcos          !! cosinus
    REAL(r_std), INTENT(in)                              :: REarth        !! Earth Radius
    REAL(r_std), INTENT(in)                              :: piv           !! pi
    REAL(r_std), DIMENSION(dx,dy), INTENT(in)            :: lons, lats    !! longitude and latitudes (degrees)
    !! 0.2 Modified variables
    !
    !! 0.3 Output variables
    !
    REAL(r_std), DIMENSION(dx,dy,2), INTENT(out)         :: resolin       !! distance between grid points (km)
    !
    !! 0.4 Local variables
    !
    INTEGER                                              :: ip,jp
    REAL(r_std)                                          :: coslat

    DO ip=1,dx
       DO jp=1,dy
          !
          ! Resolution in longitude
          !
          coslat = MAX( COS( lats(ip,jp) * piv/180. ), mcos )     
          IF ( ip .EQ. 1 ) THEN
             resolin(ip,jp,1) = ABS( lons(ip+1,jp) - lons(ip,jp) ) * piv/180. * REarth * coslat
          ELSEIF ( ip .EQ. dx ) THEN
             resolin(ip,jp,1) = ABS( lons(ip,jp) - lons(ip-1,jp) ) * piv/180. * REarth * coslat
          ELSE
             resolin(ip,jp,1) = ABS( lons(ip+1,jp) - lons(ip-1,jp) )/2. * piv/180. * REarth * coslat
          ENDIF
          !
          ! Resolution in latitude
          !
          IF ( jp .EQ. 1 ) THEN
             resolin(ip,jp,2) = ABS( lats(ip,jp) - lats(ip,jp+1) ) * piv/180. * REarth
          ELSEIF ( jp .EQ. dy ) THEN
             resolin(ip,jp,2) = ABS( lats(ip,jp-1) - lats(ip,jp) ) * piv/180. * REarth
          ELSE
             resolin(ip,jp,2) =  ABS( lats(ip,jp-1) - lats(ip,jp+1) )/2. * piv/180. * REarth
          ENDIF
          !
       ENDDO
    ENDDO

  END SUBROUTINE interpweight_calc_resolution_in

!! ================================================================================================================================
!! SUBROUTINE   : interpweight_modifying_input1D
!!
!>\BRIEF        ! modify the initial 1D values from the given file
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ivar
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE interpweight_modifying_input1D(dim1, dim2, dim3, dim4, stime, etime, zeroval, ivar)

    USE constantes_var
    USE ioipsl_para

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !! 0.1 Input variables
    INTEGER, INTENT(in)                                  :: dim1, dim2, dim3, dim4  !! size of the dimensions to get
    INTEGER, INTENT(in)                                  :: stime, etime  !! starting/ending time step to get
    REAL(r_std), INTENT(in)                              :: zeroval       !! zero value
    !! 0.2 Modified variables
    REAL(r_std), DIMENSION(dim1), INTENT(inout)          :: ivar          !! modified input data


    WHERE (ivar(:) < zeroval )
      ivar(:) = zeroval
    END WHERE

  END SUBROUTINE interpweight_modifying_input1D

!! ================================================================================================================================
!! SUBROUTINE   : interpweight_modifying_input2D
!!
!>\BRIEF        ! modify the initial 2D values from the given file
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ivar
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE interpweight_modifying_input2D(dim1, dim2, dim3, dim4, stime, etime, zeroval, ivar)

    USE constantes_var
    USE ioipsl_para

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !! 0.1 Input variables
    INTEGER, INTENT(in)                                  :: dim1, dim2, & !! size of the dimensions to get
      dim3, dim4
    INTEGER, INTENT(in)                                  :: stime, etime  !! starting/ending time step to get
    REAL(r_std), INTENT(in)                              :: zeroval       !! zero value
    !! 0.2 Modified variables
    REAL(r_std), DIMENSION(dim1,dim2), INTENT(inout)     :: ivar          !! modified input data

    WHERE (ivar(:,:) < zeroval )
      ivar(:,:) = zeroval
    END WHERE

  END SUBROUTINE interpweight_modifying_input2D

!! ================================================================================================================================
!! SUBROUTINE   : interpweight_modifying_input3D
!!
!>\BRIEF        ! modify the initial 3D values from the given file
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ivar
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE interpweight_modifying_input3D(dim1, dim2, dim3, dim4, stime, etime, zeroval, ivar)

    USE constantes_var
    USE ioipsl_para

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !! 0.1 Input variables
    INTEGER, INTENT(in)                                  :: dim1, dim2, & !! size of the dimensions to get
      dim3, dim4
    INTEGER, INTENT(in)                                  :: stime, etime  !! starting/ending time step to get
    REAL(r_std), INTENT(in)                              :: zeroval       !! zero value
    !! 0.2 Modified variables
    REAL(r_std), DIMENSION(dim1,dim2,dim3), INTENT(inout):: ivar          !! modified input data

    WHERE (ivar(:,:,:) < zeroval )
      ivar(:,:,:) = zeroval
    END WHERE

  END SUBROUTINE interpweight_modifying_input3D

!! ================================================================================================================================
!! SUBROUTINE   : interpweight_modifying_input4D
!!
!>\BRIEF        ! modify the initial 4D values from the given file
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ivar
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE interpweight_modifying_input4D(dim1, dim2, dim3, dim4, stime, etime, zeroval, ivar)

    USE constantes_var
    USE ioipsl_para

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !! 0.1 Input variables
    INTEGER, INTENT(in)                                  :: dim1, dim2, & !! size of the dimensions to get
      dim3, dim4
    INTEGER, INTENT(in)                                  :: stime, etime  !! starting/ending time step to get
    REAL(r_std), INTENT(in)                              :: zeroval       !! zero value
    !! 0.2 Modified variables
    REAL(r_std), DIMENSION(dim1,dim2,dim3,dim4), INTENT(inout) :: ivar          !! modified input data

    WHERE (ivar(:,:,:,:) < zeroval )
      ivar(:,:,:,:) = zeroval
    END WHERE

  END SUBROUTINE interpweight_modifying_input4D

!! ================================================================================================================================
!! SUBROUTINE   : interpweight_interpweight_1D
!!
!>\BRIEF        ! mask 1D input values
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): msk
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE interpweight_masking_input1D(dx, dy, d1, d2, d3, d4, mtype, mvalues, var,                &
    msk)

    USE constantes_var

    IMPLICIT NONE
  
    !! 0. Variables and parameter declaration
    !! 0.1 Input variables
    INTEGER, INTENT(in)                                  :: dx, dy, d1, & !! shapes of the input variables
      d2, d3, d4
    CHARACTER(LEN=50), INTENT(in)                        :: mtype         !! Type of masking
                                                                          !!   'nomask': no-mask is applied
                                                                          !!   'mbelow': take values below values(1)
                                                                          !!   'mabove': take values above values(1)
                                                                          !!   'msumrange': take values within 2 ranges;
                                                                          !!      values(2) <= SUM(vals(k)) <= values(1)
                                                                          !!      values(1) < SUM(vals(k)) <= values(3): 
                                                                          !!        renormalize
                                                                          !!   'var': mask values are taken from a 
                                                                          !!     variable (>0)
    REAL(r_std), DIMENSION(3), INTENT(in)                :: mvalues       !! values to use to mask (according to `mtype') 
    !! 0.2 Modified variables
    REAL(r_std), DIMENSION(dx), INTENT(inout)            :: var           !! variable from which provide the mask
    INTEGER(i_std), DIMENSION(dx), INTENT(inout)         :: msk           !! mask to retrieve
    !! 0.3 Output variables
    !! 0.4 Local variables
    INTEGER                                              :: i,j,k,l

    SELECT CASE (mtype)
      CASE('nomask')
        msk=un
      CASE('mbelow')
        ! e.g.: 
        ! Exclude the points where there is never a LAI value. It is probably 
        ! an ocean point.
        !
        DO i=1,dx
          IF ( var(i) < mvalues(1) ) msk(i) = un
        END DO
      CASE('mabove')
        DO i=1,dx
          IF ( var(i) > mvalues(1) ) msk(i) = un 
        ENDDO
      CASE('msumrange')
        CALL ipslerr_p(3, 'interpweight_masking_input1D',"'msumrange' no sens on a 1D variable",'','')
      CASE DEFAULT
        CALL ipslerr_p(3, 'interpweight_masking_input1D',"mask with '" // TRIM(mtype) // "' not ready !!",'','')
    END SELECT

  END SUBROUTINE interpweight_masking_input1D

!! ================================================================================================================================
!! SUBROUTINE   : interpweight_masking_input2D
!!
!>\BRIEF        ! mask 2D input values
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): msk
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE interpweight_masking_input2D(dx, dy, d1, d2, d3, d4, mtype, mvalues, var,                &
    msk)

    USE constantes_var

    IMPLICIT NONE
  
    !! 0. Variables and parameter declaration
    !! 0.1 Input variables
    INTEGER, INTENT(in)                                  :: dx, dy, d1, & !! shapes of the input variables
      d2, d3, d4
    CHARACTER(LEN=50), INTENT(in)                        :: mtype         !! Type of masking
                                                                          !!   'nomask': no-mask is applied
                                                                          !!   'mbelow': take values below values(1)
                                                                          !!   'mabove': take values above values(1)
                                                                          !!   'msumrange': take values within 2 ranges;
                                                                          !!      values(2) <= SUM(vals(k)) <= values(1)
                                                                          !!      values(1) < SUM(vals(k)) <= values(3): 
                                                                          !!        renormalize
                                                                          !!   'var': mask values are taken from a 
                                                                          !!     variable (>0)
    REAL(r_std), DIMENSION(3), INTENT(in)                :: mvalues       !! values to use to mask (according to `mtype') 
    !! 0.2 Modified variables
    REAL(r_std), DIMENSION(dx,dy), INTENT(inout)         :: var           !! variable from which provide the mask
    INTEGER(i_std), DIMENSION(dx,dy), INTENT(inout)      :: msk           !! mask to retrieve
    !! 0.3 Output variables
    !! 0.4 Local variables
    INTEGER                                              :: i,j,k,l
    REAL(r_std)                                          :: sumval


    SELECT CASE (mtype)
      CASE('nomask')
        msk=un
      CASE('mbelow')
        ! e.g.: 
        ! Exclude the points where there is never a LAI value. It is probably 
        ! an ocean point.
        !
        DO i=1,dx
          DO j=1,dy
            IF ( var(i,j) < mvalues(1) ) msk(i,j) = un
          END DO
        END DO
      CASE('mabove')
        DO i=1,dx
          DO j=1,dy
            IF ( var(i,j) > mvalues(1) ) msk(i,j) = un 
          ENDDO
        ENDDO
      CASE('msumrange')
        IF (printlev>=3) WRITE(numout,*) &
             'interpweight_masking_input2D: masking with mask sum range. Interval:', mvalues
         DO i=1,dx
           DO j=1,dy
              sumval=SUM(var(i,:))
              IF ( sumval .GE. mvalues(2) .AND. sumval .LE. mvalues(1) ) THEN
                 msk(i,j) = un 
              ELSEIF ( sumval .GT. mvalues(1) .AND. sumval .LE. mvalues(3) ) THEN
                 ! normalization
                 var(i,j) = var(i,j) / sumval
                 msk(i,j) = un 
              ENDIF
           ENDDO
        ENDDO
      CASE DEFAULT
        CALL ipslerr_p(3, 'interpweight_masking_input2D',"mask with '" // TRIM(mtype) // "' not ready !!",'','')
    END SELECT

  END SUBROUTINE interpweight_masking_input2D

!! ================================================================================================================================
!! SUBROUTINE   : interpweight_masking_input3D
!!
!>\BRIEF        ! mask 3D input values
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): msk
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE interpweight_masking_input3D(dx, dy, d1, d2, d3, d4, mtype, mvalues, var,                &
    msk)

    USE constantes_var

    IMPLICIT NONE
  
    !! 0. Variables and parameter declaration
    !! 0.1 Input variables
    INTEGER, INTENT(in)                                  :: dx, dy, d1, & !! shapes of the input variables
      d2, d3, d4
    CHARACTER(LEN=50), INTENT(in)                        :: mtype         !! Type of masking
                                                                          !!   'nomask': no-mask is applied
                                                                          !!   'mbelow': take values below values(1)
                                                                          !!   'mabove': take values above values(1)
                                                                          !!   'msumrange': take values within 2 ranges;
                                                                          !!      values(2) <= SUM(vals(k)) <= values(1)
                                                                          !!      values(1) < SUM(vals(k)) <= values(3): 
                                                                          !!        renormalize
                                                                          !!   'var': mask values are taken from a 
                                                                          !!     variable (>0)
    REAL(r_std), DIMENSION(3), INTENT(in)                :: mvalues       !! values to use to mask (according to `mtype') 
    !! 0.2 Modified variables
    REAL(r_std), DIMENSION(dx,dy,d3), INTENT(inout)      :: var           !! variable from which provide the mask
    INTEGER(i_std), DIMENSION(dx,dy), INTENT(inout)      :: msk           !! mask to retrieve
    !! 0.3 Output variables
    !! 0.4 Local variables
    INTEGER                                              :: i,j,k,l
    REAL(r_std)                                          :: sumval


    SELECT CASE (mtype)
      CASE('nomask')
        msk=un
      CASE('mbelow')
        ! e.g.: 
        ! Exclude the points where there is never a LAI value. It is probably 
        ! an ocean point.
        !
        DO i=1,dx
          DO j=1,dy
            IF ( ANY(var(i,j,:) < mvalues(1)) ) msk(i,j) = un
          END DO
        END DO
      CASE('mabove')
        DO i=1,dx
          DO j=1,dy
            IF ( ANY(var(i,j,:) > mvalues(1)) ) msk(i,j) = un 
          ENDDO
        ENDDO
      CASE('msumrange')
         DO i=1,dx
           DO j=1,dy
             sumval=SUM(var(i,j,:))
             IF ( sumval .GE. mvalues(2) .AND. sumval .LE. mvalues(1)) THEN
                msk(i,j) = un 
             ELSEIF ( sumval .GT. mvalues(1) .AND. sumval .LE. mvalues(3)) THEN
                ! normalization
                var(i,j,:) = var(i,j,:) / sumval
                msk(i,j) = un 
             ENDIF
           ENDDO
        ENDDO
      CASE DEFAULT
        CALL ipslerr_p(3, 'interpweight_masking_input3D',"mask with '" // TRIM(mtype) // "' not ready !!",'','')
    END SELECT
  END SUBROUTINE interpweight_masking_input3D

!! ================================================================================================================================
!! SUBROUTINE   : interpweight_masking_input4D
!!
!>\BRIEF        ! mask 4D input values
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): msk
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE interpweight_masking_input4D(dx, dy, d1, d2, d3, d4, mtype, mvalues, var,                &
    msk)

    USE constantes_var

    IMPLICIT NONE
  
    !! 0. Variables and parameter declaration
    !! 0.1 Input variables
    INTEGER, INTENT(in)                                  :: dx, dy, d1, & !! shapes of the input variables
      d2, d3, d4
    CHARACTER(LEN=50), INTENT(in)                        :: mtype         !! Type of masking
                                                                          !!   'nomask': no-mask is applied
                                                                          !!   'mbelow': take values below values(1)
                                                                          !!   'mabove': take values above values(1)
                                                                          !!   'msumrange': take values within 2 ranges;
                                                                          !!      values(2) <= SUM(vals(k)) <= values(1)
                                                                          !!      values(1) < SUM(vals(k)) <= values(3): 
                                                                          !!        renormalize
                                                                          !!   'var': mask values are taken from a 
                                                                          !!     variable (>0)
    REAL(r_std), DIMENSION(3), INTENT(in)                :: mvalues       !! values to use to mask (according to `mtype') 
    !! 0.2 Modified variables
    REAL(r_std), DIMENSION(dx,dy,d3,d4), INTENT(inout)   :: var           !! variable from which provide the mask
    INTEGER(i_std), DIMENSION(dx,dy), INTENT(inout)      :: msk           !! mask to retrieve
    !! 0.3 Output variables
    !! 0.4 Local variables
    INTEGER                                              :: i,j,k,l
    REAL(r_std)                                          :: sumval


    SELECT CASE (mtype)
      CASE('nomask')
        msk=un
      CASE('mbelow')
        ! e.g.: 
        ! Exclude the points where there is never a LAI value. It is probably 
        ! an ocean point.
        !
        DO i=1,dx
          DO j=1,dy
            IF ( ANY(var(i,j,:,:) < mvalues(1)) ) msk(i,j) = un
          END DO
        END DO
      CASE('mabove')
        DO i=1,dx
          DO j=1,dy
            IF ( ANY(var(i,j,:,:) > mvalues(1)) ) msk(i,j) = un 
          ENDDO
        ENDDO
      CASE('msumrange')
        IF (printlev>=3) WRITE(numout,*) &
             'interpweight_masking_input4D: masking with mask sum range. Interval:', mvalues
         DO i=1,dx
           DO j=1,dy
             DO l=1,d4
               sumval=SUM(var(i,j,:,l))
               IF ( sumval .GE. mvalues(2) .AND. sumval .LE. mvalues(1) .AND. msk(i,j) /= un) THEN
                  msk(i,j) = un 
               ELSEIF ( sumval .GT. mvalues(1) .AND. sumval .LE. mvalues(3)) THEN
                  ! normalization
                  var(i,j,:,l) = var(i,j,:,l) / sumval
                  msk(i,j) = un 
               ENDIF
             ENDDO
           ENDDO
        ENDDO
      CASE DEFAULT
        CALL ipslerr_p(3, 'interpweight_masking_input4D',"mask with '" // TRIM(mtype) // "' not ready !!",'','')
    END SELECT

  END SUBROUTINE interpweight_masking_input4D

!! ================================================================================================================================
!! SUBROUTINE   : interpweight_provide_fractions1D
!!
!>\BRIEF        ! provide fractions from a 1D incoming variable
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ovar, aovar
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE interpweight_provide_fractions1D(tint, npt, Ntypes, valstypes,                           &
    dx, dy, d3, d4, nbmax, zeroval, ivar1D, sarea, sindex, vmax, vmin,                                &
    two, latlon, ovar, aovar)

    USE constantes_var

    IMPLICIT NONE
  
    !! 0. Variables and parameter declaration
    !! 0.1 Input variables
    CHARACTER(LEN=50), INTENT(in)                        :: tint          !! tint: type of interpolation
                                                                          !!   'XYKindTime': Input values are kinds 
                                                                          !!     of something with a temporal 
                                                                          !!     evolution on the dx*dy matrix
    INTEGER, INTENT(in)                                  :: npt           !! total number of grid points
    INTEGER, INTENT(in)                                  :: dx,dy,d3,d4   !! shape of the input values
    INTEGER, INTENT(in)                                  :: nbmax         !! maximum input grid cells to use 
                                                                          !!   for in 1 grid point
    REAL(r_std), INTENT(in)                              :: zeroval, two  !! zero and 2. real values
    INTEGER, INTENT(in)                                  :: Ntypes        !! number of types
    REAL(r_std), DIMENSION(Ntypes), INTENT(in)           :: valstypes     !! value for each type 
    REAL(r_std), DIMENSION(dx), INTENT(in)               :: ivar1D        !! values of the input data
    REAL(r_std), DIMENSION(npt,nbmax), INTENT(in)        :: sarea         !! overlaping area of a given input 
                                                                          !!   grid point
    INTEGER(i_std), DIMENSION(npt,nbmax), INTENT(in)     :: sindex        !! index of the overlaping input 
                                                                          !!   grid point
    REAL(r_std), DIMENSION(npt,2), INTENT(in)            :: latlon        !! longitude and latitude of the destiny
                                                                          !!   grid points
    !! 0.2 Modified variables
    REAL(r_std), INTENT(inout)                           :: vmax, vmin    !! thresholds of the values
    !! 0.3 Output variables
    REAL(r_std), DIMENSION(npt,Ntypes), INTENT(out)      :: ovar          !! output variable (on the nbpt space)
    REAL(r_std), DIMENSION(npt), INTENT(out)             :: aovar         !! availability of input data to 
                                                                          !!   interpolate output variable (on 
                                                                          !!   the nbpt space)
                                                                          !!     >0.: As the total sum of the 
                                                                          !!       space fraction used from the 
                                                                          !!       source
                                                                          !!     -1: Did not find any required values 
                                                                          !!       from input data
    !! 0.4 Local variables
    INTEGER                                              :: i,j,k,ij,ib,idi,jj,ip,jp,jv,it
    CHARACTER(LEN=250)                                   :: msg


! Initialization
!!
    ovar = zeroval

    SELECT CASE (tint)
      CASE ('default')
        IF (printlev>=4) WRITE(numout,*)'  interpweight_provide_fractions1D Fractions following default1D method'
        IF (ALL(ivar1D == zeroval)) THEN
          msg = "' wrong 1D input variable"
          CALL ipslerr_p(3,'interpweight_provide_fractions1D',"'" // TRIM(tint) // TRIM(msg),'','')
        END IF 
        DO ib=1,npt
          aovar(ib) = zeroval
          idi = COUNT(sarea(ib,:) > zeroval)
          IF ( idi > 0 ) THEN
            DO jj=1,idi
               ip = sindex(ib,jj)
               DO it = 1, Ntypes
                 IF (ivar1d(ip) == valstypes(it)*un) THEN
                   ovar(ib,it) = ovar(ib,it) + sarea(ib,jj)
                 END IF
               END DO
               aovar(ib) = aovar(ib) + sarea(ib,jj)
            ENDDO
            !
            ! Normalize
            !
            ovar(ib,:) = ovar(ib,:)/aovar(ib)
            !
            ! Quality variable
          ELSE
             ! No points found.
            aovar(ib) = -un
            ovar(ib,INT((vmax+vmin)/two)) = un
          ENDIF
          IF (SUM(ovar(ib,:)) > 1.0000001) THEN
            IF (printlev>=3) WRITE(numout,*) '  interpweight_provide_fractions1D On point ', ib, &
                 ' location: ', latlon(ib,2), latlon(ib,1),' total fractions=', SUM(ovar(ib,:))
            IF (printlev>=3) THEN
               WRITE(numout,*) '  type fraction: '
               DO i=1, Ntypes
                  WRITE(numout,*) i, ovar(ib,i)
               END DO
            END IF
            msg='total of fractions above 1.!'
            CALL ipslerr_p(3,'interpweight_provide_fractions1D',TRIM(msg),'','')
          END IF
        END DO

        IF (printlev >= 3 .AND. ANY(aovar == -un)) THEN
           WRITE(numout,*) 'interpweight_provide_fractions1D: '
           WRITE(numout,*) 'Some points on the model grid did not have any points on the source grid for interpolation.'
           WRITE(numout,*) 'These points are initialized by putting the average VAR ', INT((vmax+vmin)/two), ' to 1.'
           WRITE(numout,*) 'The points are (ib, lon, lat) :'
           DO ib=1, npt
              IF (aovar(ib) == -un) WRITE(numout,*) ib, latlon(ib,2), latlon(ib,1)
           END DO
        END IF
      CASE DEFAULT
        CALL ipslerr_p(3,'interpweight_provide_fractions1D',"Inerpolation type '" // TRIM(tint) // "' not ready!",'','')
    END SELECT

  END SUBROUTINE interpweight_provide_fractions1D

!! ================================================================================================================================
!! SUBROUTINE   : interpweight_provide_fractions2D
!!
!>\BRIEF        ! provide fractions from a 2D incoming variable
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ovar, aovar
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE interpweight_provide_fractions2D(tint, npt, Ntypes, valstypes,                           &
    dx, dy, d3, d4, nbmax, zeroval, ivar2D, sarea, sindex, vmax, vmin,                                &
    two, latlon, ovar, aovar)

    USE constantes_var

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !! 0.1 Input variables
    CHARACTER(LEN=50), INTENT(in)                        :: tint          !! tint: type of interpolation
                                                                          !!   'XYKindTime': Input values are kinds 
                                                                          !!     of something with a temporal 
                                                                          !!     evolution on the dx*dy matrix
    INTEGER, INTENT(in)                                  :: npt           !! total number of grid points
    INTEGER, INTENT(in)                                  :: dx,dy,d3,d4   !! shape of the input values
    INTEGER, INTENT(in)                                  :: nbmax         !! maximum input grid cells to use 
                                                                          !!   for in 1 grid point
    REAL(r_std), INTENT(in)                              :: zeroval, two  !! zero and 2. real values
    INTEGER, INTENT(in)                                  :: Ntypes        !! number of types
    REAL(r_std), DIMENSION(Ntypes), INTENT(in)           :: valstypes     !! value for each type 
    REAL(r_std), DIMENSION(dx,dy), INTENT(in)            :: ivar2D        !! values of the input data
    REAL(r_std), DIMENSION(npt,nbmax), INTENT(in)        :: sarea         !! overlaping area of a given input 
                                                                          !!   grid point
    INTEGER(i_std), DIMENSION(npt,nbmax,2), INTENT(in)   :: sindex        !! index of the overlaping input 
                                                                          !!   grid point
    REAL(r_std), DIMENSION(npt,2), INTENT(in)            :: latlon        !! longitude and latitude of the destiny
                                                                          !!   grid points
    !! 0.2 Modified variables
    REAL(r_std), INTENT(inout)                           :: vmax, vmin    !! thresholds of the values
    !! 0.3 Output variables
    REAL(r_std), DIMENSION(npt,Ntypes), INTENT(out)      :: ovar          !! output variable (on the nbpt space)
    REAL(r_std), DIMENSION(npt), INTENT(out)             :: aovar         !! availability of input data to 
                                                                          !!   interpolate output variable (on 
                                                                          !!   the nbpt space)
                                                                          !!     >0.: As the total sum of the 
                                                                          !!       space fraction used from the 
                                                                          !!       source
                                                                          !!     -1: Found any required values 
                                                                          !!       from input data
    !! 0.4 Local variables
    INTEGER                                              :: i,j,k,ij,ib,idi,jj,ip,jp,jv,it
    CHARACTER(LEN=250)                                   :: msg


! Initialization
!!
    ovar = zeroval

    SELECT CASE (tint)
      CASE ('default')
        IF (printlev>=3) WRITE(numout,*)'  interpweight_provide_fractions2D Fractions following default2D method'
        IF (ALL(ivar2D == zeroval)) THEN
          msg = "' wrong 2D input variable"
          CALL ipslerr_p(3,'interpweight_provide_fractions2D',"'" // TRIM(tint) // TRIM(msg),'','')
        END IF 
        DO ib=1,npt
          aovar(ib) = zeroval
          idi = COUNT(sarea(ib,:) > zeroval)
          IF ( idi > 0 ) THEN
            DO jj=1,idi
               ip = sindex(ib,jj,1)
               jp = sindex(ib,jj,2)
               DO it = 1, Ntypes
                 IF (ivar2d(ip,jp) == valstypes(it)*un) THEN
                   ovar(ib,it) = ovar(ib,it) + sarea(ib,jj)
                 END IF
               END DO
               aovar(ib) = aovar(ib) + sarea(ib,jj)
            ENDDO
            !
            ! Normalize
            !
            ovar(ib,:) = ovar(ib,:)/aovar(ib)
            !
            ! Quality variable
          ELSE
            ! No points on the source grid were found for the current point. 
            ! We provide here a default value.
            aovar(ib) = -un
            ovar(ib,INT((vmax+vmin)/two)) = un
          ENDIF

          IF (SUM(ovar(ib,:)) > 1.0000001) THEN
            WRITE(numout,*) '  interpweight_provide_fractions2D On point ', ib, ' location: ', latlon(ib,2),     &
              latlon(ib,1),' total fractions=', SUM(ovar(ib,:))
            WRITE(numout,*) '  type fraction _______'
            DO i=1, Ntypes
              WRITE(numout,*) i, ovar(ib,i)
            END DO
            msg='total of fractions above 1.!'
            CALL ipslerr_p(3,'interpweight_provide_fractions2D',TRIM(msg),'','')
          END IF
        END DO

        IF (printlev >= 3 .AND. ANY(aovar == -un)) THEN
           WRITE(numout,*) 'interpweight_provide_fractions2D: '
           WRITE(numout,*) 'Some points on the model grid did not have any points on the source grid for interpolation.'
           WRITE(numout,*) 'These points are initialized by putting the average VAR ', INT((vmax+vmin)/two), ' to 1.'
           WRITE(numout,*) 'The points are (ib, lon, lat) :'
           DO ib=1, npt
              IF (aovar(ib) == -un) WRITE(numout,*) ib, latlon(ib,2), latlon(ib,1)
           END DO
        END IF
 
      CASE DEFAULT
        CALL ipslerr_p(3,'interpweight_provide_fractions2D',"Inerpolation type '" // TRIM(tint) // "' not ready!",'','')
    END SELECT

  END SUBROUTINE interpweight_provide_fractions2D

!! ================================================================================================================================
!! SUBROUTINE   : interpweight_provide_fractions3D
!!
!>\BRIEF        ! provide fractions from a 3D incoming variable
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ovar, aovar
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE interpweight_provide_fractions3D(tint, npt, Ntypes, valstypes,                           &
    dx, dy, d3, d4, nbmax, zeroval, ivar3D, sarea, sindex, vmax, vmin,                                &
    two, latlon, ovar, aovar)

    USE constantes_var

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !! 0.1 Input variables
    CHARACTER(LEN=50), INTENT(in)                        :: tint          !! tint: type of interpolation
                                                                          !!   'XYKindTime': Input values are kinds 
                                                                          !!     of something with a temporal 
                                                                          !!     evolution on the dx*dy matrix
    INTEGER, INTENT(in)                                  :: npt           !! total number of grid points
    INTEGER, INTENT(in)                                  :: dx,dy,d3,d4   !! shape of the input values
    INTEGER, INTENT(in)                                  :: nbmax         !! maximum input grid cells to use 
                                                                          !!   for in 1 grid point
    REAL(r_std), INTENT(in)                              :: zeroval, two  !! zero and 2. real values
    INTEGER, INTENT(in)                                  :: Ntypes        !! number of types
    REAL(r_std), DIMENSION(Ntypes), INTENT(in)           :: valstypes     !! value for each type 
    REAL(r_std), DIMENSION(dx,dy,d3), INTENT(in)         :: ivar3D        !! values of the input data
    REAL(r_std), DIMENSION(npt,nbmax), INTENT(in)        :: sarea         !! overlaping area of a given input 
                                                                          !!   grid point
    INTEGER(i_std), DIMENSION(npt,nbmax,2), INTENT(in)   :: sindex        !! index of the overlaping input 
                                                                          !!   grid point
    REAL(r_std), DIMENSION(npt,2), INTENT(in)            :: latlon        !! longitude and latitude of the destiny
                                                                          !!   grid points
    !! 0.2 Modified variables
    REAL(r_std), DIMENSION(Ntypes), INTENT(inout)        :: vmax, vmin    !! thresholds of the values
    !! 0.3 Output variables
    REAL(r_std), DIMENSION(npt,Ntypes), INTENT(out)      :: ovar          !! output variable (on the nbpt space)
    REAL(r_std), DIMENSION(npt), INTENT(out)             :: aovar         !! availability of input data to 
                                                                          !!   interpolate output variable (on 
                                                                          !!   the nbpt space)
                                                                          !!     >0.: As the total sum of the 
                                                                          !!       space fraction used from the 
                                                                          !!       source
                                                                          !!     -1: Found any required values 
                                                                          !!       from input data
    !! 0.4 Local variables
    INTEGER                                              :: i,j,k,ij,ib,idi,jj,ip,jp,jv,it
    CHARACTER(LEN=250)                                   :: msg
    INTEGER                                              :: Nzeros
    INTEGER                                              :: ALLOC_ERR 
    REAL(r_std), ALLOCATABLE, DIMENSION(:)               :: subfrac
    REAL(r_std)                                          :: sumpt


! Initialization
!!
    ovar = zeroval

    IF (Ntypes /= d3) THEN
      msg = 'Different number of types than third dimension in the input data!!'
      CALL ipslerr_p(3,'interpweight_provide_fractions3D',TRIM(msg),'','')
    END IF

    Nzeros = 0
    SELECT CASE (tint)
      CASE ('default')
        IF (printlev>=3) WRITE(numout,*)'Fractions following default3D method'
        IF (ALL(ivar3D == zeroval)) THEN
          msg = "' wrong 3D input variable"
          CALL ipslerr_p(3,'interpweight_provide_fractions3D',"'" // TRIM(tint) // TRIM(msg),'','')
        END IF
        DO ib=1,npt
          aovar(ib) = zeroval
          idi = COUNT(sarea(ib,:) > zeroval)
          IF ( idi > 0 ) THEN
            sumpt = zeroval
            DO jj=1,idi
               ip = sindex(ib,jj,1)
               jp = sindex(ib,jj,2)
               DO it = 1, Ntypes
! Do not get that areas with zero input value and with missing value
                 IF (ivar3D(ip,jp,it) > zeroval .AND. ivar3D(ip,jp,it) < two) THEN
                   ovar(ib,it) = ovar(ib,it) + ivar3D(ip,jp,it)*sarea(ib,jj)
                   sumpt = sumpt + ivar3D(ip,jp,it)*sarea(ib,jj)
                 END IF
               END DO
               aovar(ib) = aovar(ib) + sarea(ib,jj)
            ENDDO
            !
            ! Normalize
            !
            ovar(ib,:) = ovar(ib,:)/aovar(ib)
            !
            ! Quality variable
          ELSE
            Nzeros = Nzeros + 1
            aovar(ib) = -un

            DO it = 1, d3
              ovar(ib,INT((vmax(it)+vmin(it))/two)) = un
            END DO
          ENDIF
! Tests
!
          IF (SUM(ovar(ib,:)) > 1.0000001) THEN
            idi = COUNT(sarea(ib,:) > zeroval)
            ALLOCATE(subfrac(idi), STAT=ALLOC_ERR)
            IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'interpweight_provide_fractions3D','Problem in allocation of subfrac','','')
            IF (printlev>=3) THEN
               WRITE(numout,*) 'interpweight_provide_fractions3D On point ', ib, &
                    ' location: ', latlon(ib,2), latlon(ib,1),' total fractions=', SUM(ovar(ib,:)), &
                    ' total points from invar: ', idi, ' total invar area: ', aovar(ib)
               WRITE(numout,*) '  type ivar3D jj sarea fraction: '

               DO it = 1, Ntypes
                  DO jj=1,idi
                     ip = sindex(ib,jj,1)
                     jp = sindex(ib,jj,2)
                     IF (ivar3D(ip,jp,it) > zeroval) THEN
                        WRITE(numout,*) it, ivar3D(ip,jp,it), jj, sarea(ib,jj), ovar(ib,it)
                     END IF
                  END DO
               END DO
            END IF
            msg='total of fractions above 1.!'
            CALL ipslerr_p(3,'interpweight_provide_fractions3D',TRIM(msg),'','')
         END IF

        END DO

        IF (printlev >= 3 .AND. ANY(aovar == -un)) THEN
           WRITE(numout,*) 'interpweight_provide_fractions3D: '
           WRITE(numout,*) 'Some points on the model grid did not have any points on the source grid for interpolation.'
           WRITE(numout,*) 'These points are initialized by putting the average VAR ', INT((vmax+vmin)/two), ' to 1.'
           WRITE(numout,*) 'The points are (ib, lon, lat) :'
           DO ib=1, npt
              IF (aovar(ib) == -un) WRITE(numout,*) ib, latlon(ib,2), latlon(ib,1)
           END DO
        END IF

        IF (COUNT(aovar == -1) == npt) THEN
          msg='No grid points could be interpolated from the source grid for the current process'
          CALL ipslerr_p(2,'interpweight_provide_fractions3D',TRIM(msg),'','')
        END IF

        IF (printlev >=3) WRITE(numout,*) 'interpweight_provide_fractions3D: Nzeros=', &
             Nzeros,' aovar=', COUNT(aovar == -1), ' nbpt=', npt

        IF (COUNT(aovar == -1) /= Nzeros) THEN
          msg='Something went wrong COUNT(aovar == -1) /= Nzeros'
          CALL ipslerr_p(3,'interpweight_provide_fractions3D',TRIM(msg),'','')
        END IF
 
      CASE DEFAULT
        CALL ipslerr_p(3,'interpweight_provide_fractions3D',"Inerpolation type '" // TRIM(tint) // "' not ready!",'','')
    END SELECT

  END SUBROUTINE interpweight_provide_fractions3D

!! ================================================================================================================================
!! SUBROUTINE   : interpweight_provide_fractions4D
!!
!>\BRIEF        ! provide fractions from a 4D incoming variable
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ovar, aovar
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE interpweight_provide_fractions4D(tint, npt, Ntypes, valstypes,                           &
    dx, dy, d3, d4, nbmax, zeroval, ivar4D, sarea, sindex, vmax, vmin,                                &
    two, latlon, ovar, aovar)

    USE constantes_var

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !! 0.1 Input variables
    CHARACTER(LEN=50), INTENT(in)                        :: tint          !! tint: type of interpolation
                                                                          !!   'XYKindTime': Input values are kinds 
                                                                          !!     of something with a temporal 
                                                                          !!     evolution on the dx*dy matrix
    INTEGER, INTENT(in)                                  :: npt           !! total number of grid points
    INTEGER, INTENT(in)                                  :: dx,dy,d3,d4   !! shape of the input values
    INTEGER, INTENT(in)                                  :: nbmax         !! maximum input grid cells to use 
                                                                          !!   for in 1 grid point
    REAL(r_std), INTENT(in)                              :: zeroval, two  !! zero and 2. real values
    INTEGER, INTENT(in)                                  :: Ntypes        !! number of types
    REAL(r_std), DIMENSION(Ntypes), INTENT(in)           :: valstypes     !! value for each type 
    REAL(r_std), DIMENSION(dx,dy,d3,d4), INTENT(in)      :: ivar4D        !! values of the input data
    REAL(r_std), DIMENSION(npt,nbmax), INTENT(in)        :: sarea         !! overlaping area of a given input 
                                                                          !!   grid point
    INTEGER(i_std), DIMENSION(npt,nbmax,2), INTENT(in)   :: sindex        !! index of the overlaping input 
                                                                          !!   grid point
    REAL(r_std), DIMENSION(npt,2), INTENT(in)            :: latlon        !! longitude and latitude of the destiny
                                                                          !!   grid points
    !! 0.2 Modified variables
    REAL(r_std), DIMENSION(Ntypes), INTENT(inout)        :: vmax, vmin    !! thresholds of the values
    !! 0.3 Output variables
    REAL(r_std), DIMENSION(npt,Ntypes,d4), INTENT(out)   :: ovar          !! output variable (on the nbpt space)
    REAL(r_std), DIMENSION(npt), INTENT(out)             :: aovar         !! availability of input data to 
                                                                          !!   interpolate output variable (on 
                                                                          !!   the nbpt space)
                                                                          !!     >0.: As the total sum of the 
                                                                          !!       space fraction used from the 
                                                                          !!       source
                                                                          !!     -1: Found any required values 
                                                                          !!       from input data
    !! 0.4 Local variables
    INTEGER                                              :: i,j,k,ij,ib,idi,jj,ip,jp,jv,it,itt
    CHARACTER(LEN=250)                                   :: msg


! Initialization
!!
    ovar = zeroval

    IF (Ntypes /= d3) THEN
      msg = 'Different number of types than third dimension in the input data!!'
      CALL ipslerr_p(3,'interpweight_provide_fractions4D',TRIM(msg),'','')
    END IF

    SELECT CASE (tint)
      CASE ('default')
        IF (printlev>=3) WRITE(numout,*)'Fractions following default4D method'
        IF (ALL(ivar4D == zeroval)) THEN
          msg = "' wrong 4D input variable"
          CALL ipslerr_p(3,'interpweight_provide_fractions4D',"'" // TRIM(tint) // TRIM(msg),'','')
        END IF
        DO ib=1,npt
          aovar(ib) = zeroval
          idi = COUNT(sarea(ib,:) > zeroval)
          IF ( idi > 0 ) THEN
            DO jj=1,idi
               ip = sindex(ib,jj,1)
               jp = sindex(ib,jj,2)
               DO it = 1, Ntypes
                 DO itt = 1, d4
                   ! Do not get that areas with zero input value and with missing value
                   IF (ivar4D(ip,jp,it,itt) > zeroval .AND. ivar4D(ip,jp,it,itt) < 20.) THEN
                     ovar(ib,it,itt) = ovar(ib,it,itt) + ivar4D(ip,jp,it,itt)*sarea(ib,jj)
                   END IF
                 END DO
               END DO
               aovar(ib) = aovar(ib) + sarea(ib,jj)
            ENDDO
            !
            ! Normalize
            !
            ovar(ib,:,:) = ovar(ib,:,:)/aovar(ib)
            !
            ! Quality variable
          ELSE
            ! No points found. Set the average of input threashold (vmax, vmin) as value. 
            aovar(ib) = -un
            IF (printlev >= 3) WRITE(numout,*) 'interpweight_provide_fractions4D: ',&
                 'No points found for interpolating of point ib=',ib,&
                 ' with location lon, lat = ',latlon(ib,2), latlon(ib,1)
            DO it=1, Ntypes
              DO itt=1, d4
                ovar(ib,INT((vmax(it)+vmin(it))/two),itt) = un
              END DO
            END DO
          ENDIF
!          DO itt = 1, d4
!             IF (SUM(ovar(ib,:,itt)) > 1.0000001) THEN
!                IF (printlev>=3) THEN
!                   WRITE(numout,*) '  interpweight_provide_fractions4D On point ', ib, &
!                        ' location: ', latlon(ib,2), latlon(ib,1),' itt=', itt, &
!                        ' total fractions=', SUM(ovar(ib,:,itt))
!                   WRITE(numout,*) '  type fraction:'
!                   DO i=1, Ntypes
!                      WRITE(numout,*) i, ovar(ib,:,itt)
!                   END DO
!                END IF
!                msg='total of fractions above 1.!'
!                CALL ipslerr_p(3,'interpweight_provide_fractions4D',TRIM(msg),'','')
!             END IF
!          END DO
       END DO
       
      CASE DEFAULT
        CALL ipslerr_p(3,'interpweight_provide_fractions4D',"Inerpolation type '" // TRIM(tint) // "' not ready!",'','')
    END SELECT

  END SUBROUTINE interpweight_provide_fractions4D

!! ================================================================================================================================
!! SUBROUTINE   : interpweight_provide_interpolation2D
!!
!>\BRIEF        ! perform the interpolation for a continuos 2D field
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ovar, aovar
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE interpweight_provide_interpolation2D(tint, npt, d1, d2, dx, dy, d3, d4,                  &
    nbmax, zeroval, ivar2D, sarea, sindex, vmax, vmin, one, two, latlon,                              &
    defaultval, defaultNOval, ovar, aovar)

    USE constantes_var

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !! 0.1 Input variables
    CHARACTER(LEN=50), INTENT(in)                        :: tint          !! tint: type of interpolation
                                                                          !!   'XYKindTime': Input values are kinds 
                                                                          !!     of something with a temporal 
                                                                          !!     evolution on the dx*dy matrix
    INTEGER, INTENT(in)                                  :: npt           !! total number of grid points
    INTEGER, INTENT(in)                                  :: d1, d2, d3, & !! shape of the input values
      d4, dx, dy
    INTEGER, INTENT(in)                                  :: nbmax         !! maximum input grid cells to use 
                                                                          !!   for in 1 grid point
    REAL(r_std), INTENT(in)                              :: zeroval,one,two !! zero 1. and 2. real values
    REAL(r_std), DIMENSION(dx,dy), INTENT(in)            :: ivar2D        !! values of the input data
    REAL(r_std), DIMENSION(npt,nbmax), INTENT(in)        :: sarea         !! overlaping area of a given input 
                                                                          !!   grid point
    INTEGER(i_std), DIMENSION(npt,nbmax,2), INTENT(in)   :: sindex        !! index of the overlaping input 
                                                                          !!   grid point
    REAL(r_std), DIMENSION(npt,2), INTENT(in)            :: latlon        !! longitude and latitude of the destiny
                                                                          !!   grid points
    REAL(r_std), INTENT(in)                              :: defaultval    !! default interpolated value, only usef if tint='default'
    REAL(r_std), INTENT(in)                              :: defaultNOval  !! default interpolated NO value, only used if tint='slopecalc'
    !! 0.2 Modified variables
    REAL(r_std), INTENT(in)                              :: vmax, vmin    !! thresholds of the values
    !! 0.3 Output variables
    REAL(r_std), DIMENSION(npt), INTENT(out)             :: ovar          !! output variable (on the nbpt space)
    REAL(r_std), DIMENSION(npt), INTENT(out)             :: aovar         !! availability of input data to 
                                                                          !!   interpolate output variable (on 
                                                                          !!   the nbpt space)
                                                                          !!     >0.: As the total sum of the 
                                                                          !!       space fraction used from the 
                                                                          !!       source
                                                                          !!     -1: Found any required values 
                                                                          !!       from input data
    !! 0.4 Local variables
    INTEGER                                              :: i,j,k,ij,ib,idi,jj,ip,jp,jv,it
    INTEGER                                              :: idi_last
    REAL(r_std)                                          :: int1D
    CHARACTER(LEN=250)                                   :: msg


! Initialization
!!
    ovar = zeroval

    SELECT CASE (tint)
      CASE ('default')
        IF (printlev>=3) WRITE(numout,*)'Fractions following default method'
        IF (ALL(ivar2D == zeroval)) THEN
          msg = "' wrong 2D input variable"
          CALL ipslerr_p(3,'interpweight_provide_interpolation2D',"'" // TRIM(tint) // TRIM(msg),'','')
        END IF

        DO ib = 1, npt
          !-
          !- Reinfiltration coefficient due to the slope: Calculation with parameteres maxlope_ro 
          !-
          int1D = zeroval

          ! Initialize last index to the highest possible 
          idi_last=nbmax
          aovar(ib) = zero
          DO idi=1, nbmax
             ! Leave the do loop if all sub areas are treated, sub_area <= 0
             IF ( sarea(ib,idi) <= zeroval ) THEN
                ! Set last index to the last one used
                idi_last=idi-1
                ! Exit do loop
                EXIT
             END IF

             ip = sindex(ib,idi,1)
             jp = sindex(ib,idi,2)

             int1D = int1D + ivar2D(ip,jp) * sarea(ib,idi)
             aovar(ib) = aovar(ib) + sarea(ib,idi)
          ENDDO

          IF ( idi_last >= 1 ) THEN
             ovar(ib) = int1D / aovar(ib)
          ELSE
             ovar(ib) = defaultval
            !
            ! Quality variable
            aovar(ib) = -1.
          ENDIF
        ENDDO
      CASE ('slopecalc')
        IF (printlev>=3) WRITE(numout,*)'Fractions following slopecalc method'
        IF (ALL(ivar2D == zeroval)) THEN
          msg = "' requires a 2D input variable"
          CALL ipslerr_p(3,'interpweight_provide_interpolation2D',"'" // TRIM(tint) // TRIM(msg),'','')
        END IF

        DO ib = 1, npt
          !-
          !- Reinfiltration coefficient due to the slope: Calculation with parameteres maxlope_ro 
          !-
          int1D = zeroval

          ! Initialize last index to the highest possible 
          idi_last=nbmax
          aovar(ib) = zero
          DO idi=1, nbmax
             ! Leave the do loop if all sub areas are treated, sub_area <= 0
             IF ( sarea(ib,idi) <= zeroval ) THEN
                ! Set last index to the last one used
                idi_last=idi-1
                ! Exit do loop
                EXIT
             END IF

             ip = sindex(ib,idi,1)
             jp = sindex(ib,idi,2)

             int1D = int1D + MIN(ivar2D(ip,jp)/defaultNOval,un) * sarea(ib,idi)
             aovar(ib) = aovar(ib) + sarea(ib,idi)
          ENDDO

          IF ( idi_last >= 1 ) THEN
             ovar(ib) = un - int1D / aovar(ib)
          ELSE
             ovar(ib) = defaultval
            !
            ! Quality variable
            aovar(ib) = -1.
          ENDIF
        ENDDO

      CASE DEFAULT
        CALL ipslerr_p(3,'interpweight_provide_interpolation2D',"Inerpolation type '" // TRIM(tint) // "' not ready!",'','')
    END SELECT

  END SUBROUTINE interpweight_provide_interpolation2D

!! ================================================================================================================================
!! SUBROUTINE   : interpweight_provide_interpolation3D
!!
!>\BRIEF        ! perform the interpolation for a continuos 3D field
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ovar, aovar
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE interpweight_provide_interpolation3D(tint, npt, d1, d2, dx, dy, d3, d4,                  &
    nbmax, zeroval, ivar3D, sarea, sindex, vmax, vmin, one, two, latlon,                              &
    defaultval, defaultNOval, ovar, aovar)

    USE constantes_var

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !! 0.1 Input variables
    CHARACTER(LEN=50), INTENT(in)                        :: tint          !! tint: type of interpolation
                                                                          !!   'XYKindTime': Input values are kinds 
                                                                          !!     of something with a temporal 
                                                                          !!     evolution on the dx*dy matrix
    INTEGER, INTENT(in)                                  :: npt           !! total number of grid points
    INTEGER, INTENT(in)                                  :: d1, d2, d3, & !! shape of the input values
      d4, dx, dy
    INTEGER, INTENT(in)                                  :: nbmax         !! maximum input grid cells to use 
                                                                          !!   for in 1 grid point
    REAL(r_std), INTENT(in)                              :: zeroval,two,one  !! zero 1. and 2. real values
    REAL(r_std), DIMENSION(dx,dy,d3), INTENT(in)      :: ivar3D        !! values of the input data
    REAL(r_std), DIMENSION(npt,nbmax), INTENT(in)        :: sarea         !! overlaping area of a given input 
                                                                          !!   grid point
    INTEGER(i_std), DIMENSION(npt,nbmax,2), INTENT(in)   :: sindex        !! index of the overlaping input 
                                                                          !!   grid point
    REAL(r_std), DIMENSION(npt,2), INTENT(in)            :: latlon        !! longitude and latitude of the destiny
                                                                          !!   grid points
    REAL(r_std), DIMENSION(d1), INTENT(in)            :: defaultval    !! default interpolated value
    REAL(r_std), INTENT(in)                              :: defaultNOval  !! default interpolated NO value
    !! 0.2 Modified variables
    REAL(r_std), INTENT(in)                              :: vmax, vmin    !! thresholds of the values
    !! 0.3 Output variables
    REAL(r_std), DIMENSION(npt,d1), INTENT(out)       :: ovar          !! output variable (on the nbpt space)
    REAL(r_std), DIMENSION(npt), INTENT(out)             :: aovar         !! availability of input data to 
                                                                          !!   interpolate output variable (on 
                                                                          !!   the nbpt space)
                                                                          !!     >0.: As the total sum of the 
                                                                          !!       space fraction used from the 
                                                                          !!       source
                                                                          !!     -1: Found any required values 
                                                                          !!       from input data
    !! 0.4 Local variables
    INTEGER                                              :: i,j,k,ij,ib,idi,jj,ip,jp,jv,it
    INTEGER                                              :: ilf, fopt
    REAL(r_std)                                          :: totarea
    CHARACTER(LEN=250)                                   :: msg
! Debug
    INTEGER                                              :: nbptlL
    REAL(r_std), DIMENSION(npt,d1,d2)                    :: ovartst


! Initialization
!!
    ovar = zero
    ovartst = zero

    SELECT CASE (tint)
      CASE ('default')
        IF (printlev>=3) WRITE(numout,*)'Fractions following default method'
        IF ( ALL(ivar3D == zeroval) ) THEN
          msg = "' wrong 3D input variable"
          CALL ipslerr_p(3,'provide_interpolation3D',"'" // TRIM(tint) // TRIM(msg),'','')
        END IF

        DO ib = 1, npt

          fopt = COUNT(sarea(ib,:) > zero)
          aovar(ib) = zero

          IF ( fopt > 0 ) THEN
            !-
            !- Reinfiltration coefficient
            !-
            totarea = zeroval
 
            DO ilf=1, fopt
               ! Leave the do loop if all sub areas are treated, sub_area <= 0
               ip = sindex(ib,ilf,1)
               jp = sindex(ib,ilf,2)

               DO jv=1,d1
                 ovar(ib,jv) = ovar(ib,jv) + ivar3D(ip,jp,jv)*sarea(ib,ilf)
               END DO
               totarea = totarea + sarea(ib,ilf)
            ENDDO
            ! Normalize
            ovar(ib,:) = ovar(ib,:)/totarea
            aovar(ib) = totarea
          ELSE
        ! Set defalut value for points where the interpolation fail
            IF (printlev>=3) THEN
               WRITE(numout,*) 'provide_interpolation3D: no point found !!'
               WRITE(numout,*) 'On point ', ib, ' no points were found for interpolation data. ' //      &
                    'Default values are used.'
               WRITE(numout,*) 'Location : ', latlon(ib,2), latlon(ib,1)
            END IF
            ovar(ib,ivis) = defaultval(ivis)
            ovar(ib,inir) = defaultval(inir) 
            aovar(ib) = -un
          END IF
        ENDDO

      CASE DEFAULT
        CALL ipslerr_p(3,'provide_interpolation3D',"Inerpolation type '" // TRIM(tint) // "' not ready!",'','')
    END SELECT

  END SUBROUTINE interpweight_provide_interpolation3D

!! ================================================================================================================================
!! SUBROUTINE   : interpweight_provide_interpolation2D_tsteps
!!
!>\BRIEF        ! perform the interpolation for a continuos 2D field with time steps
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ovar, aovar
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE interpweight_provide_interpolation2D_tsteps(tint, npt, d1, d2, dx, dy, d3, tsteps,       &
    nbmax, zeroval, ivar3D, sarea, sindex, vmax, vmin, one, two, latlon,                              &
    defaultval, defaultNOval, ovar, aovar)

    USE constantes_var

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !! 0.1 Input variables
    CHARACTER(LEN=50), INTENT(in)                        :: tint          !! tint: type of interpolation
                                                                          !!   'XYKindTime': Input values are kinds 
                                                                          !!     of something with a temporal 
                                                                          !!     evolution on the dx*dy matrix
    INTEGER, INTENT(in)                                  :: npt           !! total number of grid points
    INTEGER, INTENT(in)                                  :: d1, d2, d3, & !! shape of the input values
      tsteps, dx, dy
    INTEGER, INTENT(in)                                  :: nbmax         !! maximum input grid cells to use 
                                                                          !!   for in 1 grid point
    REAL(r_std), INTENT(in)                              :: zeroval,one,two !! zero 1. and 2. real values
    REAL(r_std), DIMENSION(dx,dy,tsteps), INTENT(in)      :: ivar3D        !! values of the input data
    REAL(r_std), DIMENSION(npt,nbmax), INTENT(in)        :: sarea         !! overlaping area of a given input 
                                                                          !!   grid point
    INTEGER(i_std), DIMENSION(npt,nbmax,2), INTENT(in)   :: sindex        !! index of the overlaping input 
                                                                          !!   grid point
    REAL(r_std), DIMENSION(npt,2), INTENT(in)            :: latlon        !! longitude and latitude of the destiny
                                                                          !!   grid points
    REAL(r_std), INTENT(in)                              :: defaultval    !! default interpolated value, only usef if tint='default'
    REAL(r_std), INTENT(in)                              :: defaultNOval  !! default interpolated NO value, only used if tint='slopecalc'
    !! 0.2 Modified variables
    REAL(r_std), INTENT(in)                              :: vmax, vmin    !! thresholds of the values
    !! 0.3 Output variables
    REAL(r_std), DIMENSION(npt,tsteps), INTENT(out)             :: ovar          !! output variable (on the nbpt space)
    REAL(r_std), DIMENSION(npt), INTENT(out)             :: aovar         !! availability of input data to 
                                                                          !!   interpolate output variable (on 
                                                                          !!   the nbpt space)
                                                                          !!     >0.: As the total sum of the 
                                                                          !!       space fraction used from the 
                                                                          !!       source
                                                                          !!     -1: Found any required values 
                                                                          !!       from input data
    !! 0.4 Local variables
    INTEGER                                              :: i,j,k,ij,ib,idi,jj,ip,jp,jv,it
    INTEGER                                              :: idi_last
    REAL(r_std), DIMENSION(tsteps)                        :: int1D
    CHARACTER(LEN=250)                                   :: msg


! Initialization
!!
    ovar = zeroval

    SELECT CASE (tint)
      CASE ('default')
        IF (printlev>=3) WRITE(numout,*)'Fractions following default method'
        IF (ALL(ivar3D == zeroval)) THEN
          msg = "' wrong 3D input variable"
          CALL ipslerr_p(3,'interpweight_provide_interpolation2D_tsteps',"'" // TRIM(tint) // TRIM(msg),'','')
        END IF

        DO ib = 1, npt
          !-
          !- Reinfiltration coefficient due to the slope: Calculation with parameteres maxlope_ro 
          !-
          int1D = zeroval

          ! Initialize last index to the highest possible 
          idi_last=nbmax
          aovar(ib) = zero
          DO idi=1, nbmax
             ! Leave the do loop if all sub areas are treated, sub_area <= 0
             IF ( sarea(ib,idi) <= zeroval ) THEN
                ! Set last index to the last one used
                idi_last=idi-1
                ! Exit do loop
                EXIT
             END IF

             ip = sindex(ib,idi,1)
             jp = sindex(ib,idi,2)

!            DO tstep=1, tml ! original implementation
!               data_out(ib,tstep) = data_out(ib,tstep) + data_out_file(ip,jp,tstep) * sub_area(ib,ilf) ! avg
!            ENDDO
             int1D(:) = int1D(:) + ivar3D(ip,jp,:) * sarea(ib,idi)
             aovar(ib) = aovar(ib) + sarea(ib,idi)
          ENDDO

          IF ( idi_last >= 1 ) THEN
             ovar(ib,:) = int1D(:) / aovar(ib)
          ELSE
             ovar(ib,:) = defaultval
            !
            ! Quality variable
            aovar(ib) = -1.
          ENDIF
        ENDDO
      CASE ('slopecalc')
        IF (printlev>=3) WRITE(numout,*)'Fractions following slopecalc method'
        IF (ALL(ivar3D == zeroval)) THEN
          msg = "' requires a 3D input variable"
          CALL ipslerr_p(3,'interpweight_provide_interpolation2D_tsteps',"'" // TRIM(tint) // TRIM(msg),'','')
        END IF

        DO ib = 1, npt
          !-
          !- Reinfiltration coefficient due to the slope: Calculation with parameteres maxlope_ro 
          !-
          int1D = zeroval

          ! Initialize last index to the highest possible 
          idi_last=nbmax
          aovar(ib) = zero
          DO idi=1, nbmax
             ! Leave the do loop if all sub areas are treated, sub_area <= 0
             IF ( sarea(ib,idi) <= zeroval ) THEN
                ! Set last index to the last one used
                idi_last=idi-1
                ! Exit do loop
                EXIT
             END IF

             ip = sindex(ib,idi,1)
             jp = sindex(ib,idi,2)

             int1D(:) = int1D(:) + MIN(ivar3D(ip,jp,:)/defaultNOval,un) * sarea(ib,idi)
             aovar(ib) = aovar(ib) + sarea(ib,idi)
          ENDDO

          IF ( idi_last >= 1 ) THEN
             ovar(ib,:) = un - int1D(:) / aovar(ib)
          ELSE
             ovar(ib,:) = defaultval
            !
            ! Quality variable
            aovar(ib) = -1.
          ENDIF
        ENDDO

      CASE DEFAULT
        CALL ipslerr_p(3,'interpweight_provide_interpolation2D',"Inerpolation type '" // TRIM(tint) // "' not ready!",'','')
    END SELECT

  END SUBROUTINE interpweight_provide_interpolation2D_tsteps

!! ================================================================================================================================
!! SUBROUTINE   : interpweight_provide_interpolation3D_tsteps
!!
!>\BRIEF        ! perform the interpolation for a continuos 3D field with time steps
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ovar, aovar
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE interpweight_provide_interpolation3D_tsteps(tint, npt, d1, d2, dx, dy, d3, tsteps,       &
    nbmax, zeroval, ivar4D, sarea, sindex, vmax, vmin, one, two, latlon,                              &
    defaultval, defaultNOval, ovar, aovar)

    USE constantes_var

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !! 0.1 Input variables
    CHARACTER(LEN=50), INTENT(in)                        :: tint          !! tint: type of interpolation
                                                                          !!   'XYKindTime': Input values are kinds 
                                                                          !!     of something with a temporal 
                                                                          !!     evolution on the dx*dy matrix
    INTEGER, INTENT(in)                                  :: npt           !! total number of grid points
    INTEGER, INTENT(in)                                  :: d1, d2, d3, & !! shape of the input values
      tsteps, dx, dy
    INTEGER, INTENT(in)                                  :: nbmax         !! maximum input grid cells to use 
                                                                          !!   for in 1 grid point
    REAL(r_std), INTENT(in)                              :: zeroval,one,two !! zero 1. and 2. real values
    REAL(r_std), DIMENSION(dx,dy,d3,tsteps), INTENT(in)  :: ivar4D        !! values of the input data
    REAL(r_std), DIMENSION(npt,nbmax), INTENT(in)        :: sarea         !! overlaping area of a given input 
                                                                          !!   grid point
    INTEGER(i_std), DIMENSION(npt,nbmax,2), INTENT(in)   :: sindex        !! index of the overlaping input 
                                                                          !!   grid point
    REAL(r_std), DIMENSION(npt,2), INTENT(in)            :: latlon        !! longitude and latitude of the destiny
                                                                          !!   grid points
    REAL(r_std), INTENT(in)                              :: defaultval    !! default interpolated value, only usef if tint='default'
    REAL(r_std), INTENT(in)                              :: defaultNOval  !! default interpolated NO value, only used if tint='slopecalc'
    !! 0.2 Modified variables
    REAL(r_std), INTENT(in)                              :: vmax, vmin    !! thresholds of the values
    !! 0.3 Output variables
    REAL(r_std), DIMENSION(npt,d3,tsteps), INTENT(out)   :: ovar          !! output variable (on the nbpt space)
    REAL(r_std), DIMENSION(npt), INTENT(out)             :: aovar         !! availability of input data to 
                                                                          !!   interpolate output variable (on 
                                                                          !!   the nbpt space)
                                                                          !!     >0.: As the total sum of the 
                                                                          !!       space fraction used from the 
                                                                          !!       source
                                                                          !!     -1: Found any required values 
                                                                          !!       from input data
    !! 0.4 Local variables
    INTEGER                                              :: i,j,k,ij,ib,idi,jj,ip,jp,jv,it
    INTEGER                                              :: idi_last
    REAL(r_std), DIMENSION(d3,tsteps)                        :: int2D
    CHARACTER(LEN=250)                                   :: msg


! Initialization
!!
    ovar = zeroval

    SELECT CASE (tint)
      CASE ('default')
        IF (printlev>=3) WRITE(numout,*)'Fractions following default method'
        IF (ALL(ivar4D == zeroval)) THEN
          msg = "' wrong 4D input variable"
          CALL ipslerr_p(3,'interpweight_provide_interpolation3D_tsteps',"'" // TRIM(tint) // TRIM(msg),'','')
        END IF

        DO ib = 1, npt
          !-
          !- Reinfiltration coefficient due to the slope: Calculation with parameteres maxlope_ro 
          !-
          int2D = zeroval

          ! Initialize last index to the highest possible 
          idi_last=nbmax
          aovar(ib) = zero
          DO idi=1, nbmax
             ! Leave the do loop if all sub areas are treated, sub_area <= 0
             IF ( sarea(ib,idi) <= zeroval ) THEN
                ! Set last index to the last one used
                idi_last=idi-1
                ! Exit do loop
                EXIT
             END IF

             ip = sindex(ib,idi,1)
             jp = sindex(ib,idi,2)

!            DO tstep=1, tml ! original implementation
!               data_out(ib,tstep) = data_out(ib,tstep) + data_out_file(ip,jp,tstep) * sub_area(ib,ilf) ! avg
!            ENDDO
             int2D(:,:) = int2D(:,:) + ivar4D(ip,jp,:,:) * sarea(ib,idi)
             aovar(ib) = aovar(ib) + sarea(ib,idi)
          ENDDO

          IF ( idi_last >= 1 ) THEN
             ovar(ib,:,:) = int2D(:,:) / aovar(ib)
          ELSE
             ovar(ib,:,:) = defaultval
            !
            ! Quality variable
            aovar(ib) = -1.
          ENDIF
        ENDDO
      CASE ('slopecalc')
        IF (printlev>=3) WRITE(numout,*)'Fractions following slopecalc method'
        IF (ALL(ivar4D == zeroval)) THEN
          msg = "' requires a 4D input variable"
          CALL ipslerr_p(3,'interpweight_provide_interpolation4D_tsteps',"'" // TRIM(tint) // TRIM(msg),'','')
        END IF

        DO ib = 1, npt
          !-
          !- Reinfiltration coefficient due to the slope: Calculation with parameteres maxlope_ro 
          !-
          int2D = zeroval

          ! Initialize last index to the highest possible 
          idi_last=nbmax
          aovar(ib) = zero
          DO idi=1, nbmax
             ! Leave the do loop if all sub areas are treated, sub_area <= 0
             IF ( sarea(ib,idi) <= zeroval ) THEN
                ! Set last index to the last one used
                idi_last=idi-1
                ! Exit do loop
                EXIT
             END IF

             ip = sindex(ib,idi,1)
             jp = sindex(ib,idi,2)

             int2D(:,:) = int2D(:,:) + MIN(ivar4D(ip,jp,:,:)/defaultNOval,un) * sarea(ib,idi)
             aovar(ib) = aovar(ib) + sarea(ib,idi)
          ENDDO

          IF ( idi_last >= 1 ) THEN
             ovar(ib,:,:) = un - int2D(:,:) / aovar(ib)
          ELSE
             ovar(ib,:,:) = defaultval
            !
            ! Quality variable
            aovar(ib) = -1.
          ENDIF
        ENDDO

      CASE DEFAULT
        CALL ipslerr_p(3,'interpweight_provide_interpolation3D_tsteps',"Inerpolation type '" // TRIM(tint) // "' not ready!",'','')
    END SELECT

  END SUBROUTINE interpweight_provide_interpolation3D_tsteps

!! ================================================================================================================================
!! SUBROUTINE   : interpweight_provide_interpolation4D
!!
!>\BRIEF        ! perform the interpolation for a continuos 4D field
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ovar, aovar
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE interpweight_provide_interpolation4D(tint, npt, d1, d2, dx, dy, d3, d4,                  &
    nbmax, zeroval, ivar4D, sarea, sindex, vmax, vmin, one, two, latlon,                              &
    defaultval, defaultNOval, ovar, aovar)

    USE constantes_var

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !! 0.1 Input variables
    CHARACTER(LEN=50), INTENT(in)                        :: tint          !! tint: type of interpolation
                                                                          !!   'XYKindTime': Input values are kinds 
                                                                          !!     of something with a temporal 
                                                                          !!     evolution on the dx*dy matrix
    INTEGER, INTENT(in)                                  :: npt           !! total number of grid points
    INTEGER, INTENT(in)                                  :: d1, d2, d3, & !! shape of the input values
      d4, dx, dy
    INTEGER, INTENT(in)                                  :: nbmax         !! maximum input grid cells to use 
                                                                          !!   for in 1 grid point
    REAL(r_std), INTENT(in)                              :: zeroval,two,one  !! zero 1. and 2. real values
    REAL(r_std), DIMENSION(dx,dy,d3,d4), INTENT(in)      :: ivar4D        !! values of the input data
    REAL(r_std), DIMENSION(npt,nbmax), INTENT(in)        :: sarea         !! overlaping area of a given input 
                                                                          !!   grid point
    INTEGER(i_std), DIMENSION(npt,nbmax,2), INTENT(in)   :: sindex        !! index of the overlaping input 
                                                                          !!   grid point
    REAL(r_std), DIMENSION(npt,2), INTENT(in)            :: latlon        !! longitude and latitude of the destiny
                                                                          !!   grid points
    REAL(r_std), DIMENSION(d1,d2), INTENT(in)            :: defaultval    !! default interpolated value
    REAL(r_std), INTENT(in)                              :: defaultNOval  !! default interpolated NO value
    !! 0.2 Modified variables
    REAL(r_std), INTENT(in)                              :: vmax, vmin    !! thresholds of the values
    !! 0.3 Output variables
    REAL(r_std), DIMENSION(npt,d1,d2), INTENT(out)       :: ovar          !! output variable (on the nbpt space)
    REAL(r_std), DIMENSION(npt), INTENT(out)             :: aovar         !! availability of input data to 
                                                                          !!   interpolate output variable (on 
                                                                          !!   the nbpt space)
                                                                          !!     >0.: As the total sum of the 
                                                                          !!       space fraction used from the 
                                                                          !!       source
                                                                          !!     -1: Found any required values 
                                                                          !!       from input data
    !! 0.4 Local variables
    INTEGER                                              :: i,j,k,ij,ib,idi,jj,ip,jp,jv,it
    INTEGER                                              :: ilf, fopt
    REAL(r_std)                                          :: totarea
    CHARACTER(LEN=250)                                   :: msg
! Debug
    INTEGER                                              :: nbptlL
    REAL(r_std), DIMENSION(npt,d1,d2)                    :: ovartst


! Initialization
!!
    ovar = zero
    ovartst = zero

    SELECT CASE (tint)
      CASE ('default')
        IF (printlev>=3) WRITE(numout,*)'Fractions following default method'
        IF ( ALL(ivar4D == zeroval) ) THEN
          msg = "' wrong 4D input variable"
          CALL ipslerr_p(3,'provide_interpolation4D',"'" // TRIM(tint) // TRIM(msg),'','')
        END IF

        DO ib = 1, npt

          fopt = COUNT(sarea(ib,:) > zero)
          aovar(ib) = zero

          IF ( fopt > 0 ) THEN
            !-
            !- Reinfiltration coefficient
            !-
            totarea = zeroval
 
            DO ilf=1, fopt
               ! Leave the do loop if all sub areas are treated, sub_area <= 0
               ip = sindex(ib,ilf,1)
               jp = sindex(ib,ilf,2)

               DO jv=1,d1
                 DO it=1,d2
                   ovar(ib,jv,it) = ovar(ib,jv,it) + ivar4D(ip,jp,jv,it)*sarea(ib,ilf)
                 END DO
               END DO
               totarea = totarea + sarea(ib,ilf)
            ENDDO
            ! Normalize
            ovar(ib,:,:) = ovar(ib,:,:)/totarea
            aovar(ib) = totarea
          ELSE
        ! Set defalut value for points where the interpolation fail
            IF (printlev>=3) THEN
               WRITE(numout,*) 'provide_interpolation4D: no point found !!'
               WRITE(numout,*) 'On point ', ib, ' no points were found for interpolation data. ' //      &
                    'Default values are used.'
               WRITE(numout,*) 'Location : ', latlon(ib,2), latlon(ib,1)
            END IF
            ovar(ib,ivis,:) = defaultval(ivis,:)
            ovar(ib,inir,:) = defaultval(inir,:) 
            aovar(ib) = -un
          END IF
        ENDDO

      CASE DEFAULT
        CALL ipslerr_p(3,'provide_interpolation4D',"Inerpolation type '" // TRIM(tint) // "' not ready!",'','')
    END SELECT

  END SUBROUTINE interpweight_provide_interpolation4D


!!================================================================================================================================
!! FUNCTION     : interpweight_get_var2dims_file
!!
!>\BRIEF        ! Function to get the dimensions of a given 2D variable inside a file
!!
!! DESCRIPTION : None
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : interpweight_get_var2dims_file
!!
!! REFERENCE(S) : See above, module description.
!!
!! FLOWCHART    : None
!! \n
!_================================================================================================================================
  FUNCTION interpweight_get_var2dims_file(filename, varname)

    USE netcdf

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !! 0.1 Input variables
    CHARACTER(LEN=*), INTENT(in)                         :: filename      !! filename: name of the file to open
    CHARACTER(LEN=*), INTENT(in)                         :: varname       !! varname: name of the variable 
    !! 0.2 Modified variables
    !! 0.3 Output variables
! Following: http://stackoverflow.com/questions/3828094/function-returning-an-array-in-fortran
    INTEGER, DIMENSION(2)                                :: interpweight_get_var2dims_file
    !! 0.4 Local variables
    INTEGER                                              :: nid, vid, Ndims
    INTEGER                                              :: rcode
    INTEGER, DIMENSION(2)                                :: dimsid
    CHARACTER(LEN=250)                                   :: msg


    IF (LEN_TRIM(varname) < 1) THEN
      msg = "  interpweight_get_var2dims_file: any variable name '" // TRIM(varname) // "' was provided !!"
      CALL ipslerr_p(3,'interpweight_get_var2dims_file',TRIM(msg), '', '')
    END IF

    IF (is_root_prc) THEN
       ! Open file for read only mode
       rcode = nf90_open(TRIM(filename), NF90_NOWRITE, nid)
       IF (rcode /= NF90_NOERR) CALL ipslerr_p(3, 'interpweight_get_var2dims_file', 'Error in ng90_open', &
            TRIM(nf90_strerror(rcode)),'')

       rcode = nf90_inq_varid(nid, varname, vid)
       IF (rcode /= NF90_NOERR) CALL ipslerr_p(3, 'interpweight_get_var2dims_file', 'Error in ng90_inq_varid', &
            TRIM(nf90_strerror(rcode)),'')

       rcode = nf90_inquire_variable(nid, vid, NDIMS = Ndims)
       IF (rcode /= NF90_NOERR) CALL ipslerr_p(3, 'interpweight_get_var2dims_file', 'Error in ng90_inquire_variable 1', &
            TRIM(nf90_strerror(rcode)),'')
       
       IF (Ndims /= 2) THEN
          msg = "variable '" // TRIM(varname) // "' has not 2 dimensions!!"
          CALL ipslerr_p(3,'interpweight_get_var2dims_file',TRIM(msg), '', '')
       END IF
       
       rcode = nf90_inquire_variable(nid, vid, DIMIDS = dimsid)
       IF (rcode /= NF90_NOERR) CALL ipslerr_p(3, 'interpweight_get_var2dims_file', 'Error in ng90_inquire_variable 2', &
            TRIM(nf90_strerror(rcode)),'')
       
       rcode = nf90_inquire_dimension(nid, dimsid(1), LEN = interpweight_get_var2dims_file(1))
       IF (rcode /= NF90_NOERR) CALL ipslerr_p(3, 'interpweight_get_var2dims_file', 'Error in ng90_inquire_dimension 1', &
            TRIM(nf90_strerror(rcode)),'')
       
       rcode = nf90_inquire_dimension(nid, dimsid(2), LEN = interpweight_get_var2dims_file(2))
       IF (rcode /= NF90_NOERR) CALL ipslerr_p(3, 'interpweight_get_var2dims_file', 'Error in ng90_inquire_dimension 2', &
            TRIM(nf90_strerror(rcode)),'')
       
       rcode = NF90_CLOSE(nid)
       IF (rcode /= NF90_NOERR) CALL ipslerr_p(3, 'interpweight_get_var2dims_file', 'Error in ng90_close', &
            TRIM(nf90_strerror(rcode)),'')
    END IF

    CALL bcast(interpweight_get_var2dims_file)

  END FUNCTION interpweight_get_var2dims_file

!!================================================================================================================================
!! FUNCTION     : interpweight_get_var3dims_file
!!
!>\BRIEF        ! Function to get the dimensions of a given 3D variable inside a file
!!
!! DESCRIPTION : None
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : interpweight_get_var3dims_file
!!
!! REFERENCE(S) : See above, module description.
!!
!! FLOWCHART    : None
!! \n
!_================================================================================================================================
  FUNCTION interpweight_get_var3dims_file(filename, varname)

    USE netcdf

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !! 0.1 Input variables
    CHARACTER(LEN=*), INTENT(in)                         :: filename      !! filename: name of the file to open
    CHARACTER(LEN=*), INTENT(in)                         :: varname       !! varname: name of the variable 
    !! 0.2 Modified variables
    !! 0.3 Output variables
    INTEGER, DIMENSION(3)                                :: interpweight_get_var3dims_file
    !! 0.4 Local variables
    INTEGER                                              :: nid, vid, Ndims
    INTEGER                                              :: rcode
    INTEGER, DIMENSION(3)                                :: dimsid
    CHARACTER(LEN=250)                                   :: msg


    IF (LEN_TRIM(varname) < 1) THEN
      msg = "  interpweight_get_var3dims_file: any variable name '" // TRIM(varname) // "' was provided !!"
      CALL ipslerr_p(3,'interpweight_get_var3dims_file',TRIM(msg), '', '')
    END IF

    IF (is_root_prc) THEN
       ! Open file for read only mode
       rcode = nf90_open(TRIM(filename), NF90_NOWRITE, nid)
       IF (rcode /= NF90_NOERR) CALL ipslerr_p(3,'interpweight_get_var3dims_file', 'Error in ng90_open', &
            TRIM(nf90_strerror(rcode)),'')
       
       rcode = nf90_inq_varid(nid, varname, vid)
       IF (rcode /= NF90_NOERR) CALL ipslerr_p(3,'interpweight_get_var3dims_file', 'Error in ng90_inq_varid', &
            TRIM(nf90_strerror(rcode)),'')

       rcode = nf90_inquire_variable(nid, vid, NDIMS = Ndims)
       IF (rcode /= NF90_NOERR) CALL ipslerr_p(3,'interpweight_get_var3dims_file', 'Error in ng90_inquire_variable 1', &
            TRIM(nf90_strerror(rcode)),'')
       
       IF (Ndims /= 3) THEN
          msg = "variable '" // TRIM(varname) // "' has not 3 dimensions!!"
          CALL ipslerr_p(3,'interpweight_get_var3dims_file',TRIM(msg), '', '')
       END IF
       
       rcode = nf90_inquire_variable(nid, vid, DIMIDS = dimsid)
       IF (rcode /= NF90_NOERR) CALL ipslerr_p(3,'interpweight_get_var3dims_file', 'Error in ng90_inquire_variable 2', &
            TRIM(nf90_strerror(rcode)),'')
       
       rcode = nf90_inquire_dimension(nid, dimsid(1), LEN = interpweight_get_var3dims_file(1))
       IF (rcode /= NF90_NOERR) CALL ipslerr_p(3,'interpweight_get_var3dims_file', 'Error in ng90_dimesion 1', &
            TRIM(nf90_strerror(rcode)),'')
       
       rcode = nf90_inquire_dimension(nid, dimsid(2), LEN = interpweight_get_var3dims_file(2))
       IF (rcode /= NF90_NOERR) CALL ipslerr_p(3,'interpweight_get_var3dims_file', 'Error in ng90_dimension 2', &
            TRIM(nf90_strerror(rcode)),'')
       
       rcode = nf90_inquire_dimension(nid, dimsid(3), LEN = interpweight_get_var3dims_file(3))
       IF (rcode /= NF90_NOERR) CALL ipslerr_p(3,'interpweight_get_var3dims_file', 'Error in ng90_dimension 3', &
            TRIM(nf90_strerror(rcode)),'')
       
       rcode = NF90_CLOSE(nid)
       IF (rcode /= NF90_NOERR) CALL ipslerr_p(3,'interpweight_get_var3dims_file', 'Error in ng90_close', &
            TRIM(nf90_strerror(rcode)),'')
    END IF
    CALL bcast(interpweight_get_var3dims_file)

  END FUNCTION interpweight_get_var3dims_file

!!================================================================================================================================
!! FUNCTION     : interpweight_get_var4dims_file
!!
!>\BRIEF        ! Function to get the dimensions of a given 4D variable inside a file
!!
!! DESCRIPTION : None
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : interpweight_get_var4dims_file
!!
!! REFERENCE(S) : See above, module description.
!!
!! FLOWCHART    : None
!! \n
!_================================================================================================================================
  FUNCTION interpweight_get_var4dims_file(filename, varname)

    USE netcdf

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !! 0.1 Input variables
    CHARACTER(LEN=*), INTENT(in)                         :: filename      !! filename: name of the file to open
    CHARACTER(LEN=*), INTENT(in)                         :: varname       !! varname: name of the variable 
    !! 0.2 Modified variables
    !! 0.3 Output variables
    INTEGER, DIMENSION(4)                                :: interpweight_get_var4dims_file
    !! 0.4 Local variables
    INTEGER                                              :: nid, vid, Ndims
    INTEGER                                              :: rcode
    INTEGER, DIMENSION(4)                                :: dimsid
    CHARACTER(LEN=250)                                   :: msg

    IF (LEN_TRIM(varname) < 1) THEN
      msg = "  interpweight_get_var4dims_file: any variable name '" // TRIM(varname) // "' was provided !!"
      CALL ipslerr_p(3,'interpweight_get_var4dims_file',TRIM(msg), '', '')
    END IF

    IF (is_root_prc) THEN
       ! Open file for read only mode
       rcode = nf90_open(TRIM(filename), NF90_NOWRITE, nid)
       IF (rcode /= NF90_NOERR) CALL ipslerr_p(3,'interpweight_get_var4dims_file','Problem with nf90_open', &
            TRIM(nf90_strerror(rcode)),'')
       
       rcode = nf90_inq_varid(nid, varname, vid)
       IF (rcode /= NF90_NOERR) CALL ipslerr_p(3,'interpweight_get_var4dims_file','Problem with nf90_inq_varid', &
            TRIM(nf90_strerror(rcode)),'')

       rcode = nf90_inquire_variable(nid, vid, NDIMS = Ndims)
       IF (rcode /= NF90_NOERR) CALL ipslerr_p(3,'interpweight_get_var4dims_file','Problem with nf90_inquire_variable 1', &
            TRIM(nf90_strerror(rcode)),'')

       IF (Ndims /= 4) THEN
          msg = "variable '" // TRIM(varname) // "' has not 4 dimensions!!"
          CALL ipslerr_p(3,'interpweight_get_var4dims_file',TRIM(msg), '', '')
       END IF

       rcode = nf90_inquire_variable(nid, vid, DIMIDS = dimsid)
       IF (rcode /= NF90_NOERR) CALL ipslerr_p(3,'interpweight_get_var4dims_file','Problem with nf90_inquire_variable 2', &
            TRIM(nf90_strerror(rcode)),'')
       
       rcode = nf90_inquire_dimension(nid, dimsid(1), LEN = interpweight_get_var4dims_file(1))
       IF (rcode /= NF90_NOERR) CALL ipslerr_p(3,'interpweight_get_var4dims_file','Problem with nf90_inquire_dimension 1', &
            TRIM(nf90_strerror(rcode)),'')
       
       rcode = nf90_inquire_dimension(nid, dimsid(2), LEN = interpweight_get_var4dims_file(2))
       IF (rcode /= NF90_NOERR) CALL ipslerr_p(3,'interpweight_get_var4dims_file','Problem with nf90_inquire_dimension 2', &
            TRIM(nf90_strerror(rcode)),'')

       rcode = nf90_inquire_dimension(nid, dimsid(3), LEN = interpweight_get_var4dims_file(3))
       IF (rcode /= NF90_NOERR) CALL ipslerr_p(3,'interpweight_get_var4dims_file','Problem with nf90_inquire_dimension 3', &
            TRIM(nf90_strerror(rcode)),'')
       
       rcode = nf90_inquire_dimension(nid, dimsid(4), LEN = interpweight_get_var4dims_file(4))
       IF (rcode /= NF90_NOERR) CALL ipslerr_p(3,'interpweight_get_var4dims_file','Problem with nf90_inquire_dimension 4', &
            TRIM(nf90_strerror(rcode)),'')

       rcode = NF90_CLOSE(nid)
       IF (rcode /= NF90_NOERR) CALL ipslerr_p(3,'interpweight_get_var4dims_file','Problem with nf90_close', &
            TRIM(nf90_strerror(rcode)),'')
    END iF
    CALL bcast(interpweight_get_var4dims_file)

  END FUNCTION interpweight_get_var4dims_file

!!================================================================================================================================
!! FUNCTION     : interpweight_get_varNdims_file
!!
!>\BRIEF        ! Function to get the rank(number of dimensions) of a given variable inside a file
!!
!! DESCRIPTION : None
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : interpweight_get_varNdims_file
!!
!! REFERENCE(S) : See above, module description.
!!
!! FLOWCHART    : None
!! \n
!_================================================================================================================================
  INTEGER FUNCTION interpweight_get_varNdims_file(filename, varname)

    USE netcdf

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !! 0.1 Input variables
    CHARACTER(LEN=*), INTENT(in)                         :: filename      !! filename: name of the file to open
    CHARACTER(LEN=*), INTENT(in)                         :: varname       !! varname: name of the variable 
    !! 0.2 Modified variables
    !! 0.3 Output variables
    !! 0.4 Local variables
    INTEGER                                              :: nid, vid
    INTEGER                                              :: rcode
    CHARACTER(LEN=256)                                   :: msg


    IF (LEN_TRIM(varname) < 1) THEN
      msg = "  interpweight_get_varNdims_file: any variable name '" // TRIM(varname) // "' was provided !!"
      CALL ipslerr_p(3,'interpweight_get_varNdims_file',TRIM(msg), '', '')
    END IF

    IF (is_root_prc) THEN
       
       ! Open file for read only mode
       rcode = nf90_open(TRIM(filename), NF90_NOWRITE, nid)
       IF (rcode /= NF90_NOERR) CALL ipslerr_p(3,'interpweight_get_varNdims_file','Problem with nf90_open', &
            TRIM(nf90_strerror(rcode)), '')
       
       rcode = nf90_inq_varid(nid, varname, vid)
       IF (rcode /= NF90_NOERR) CALL ipslerr_p(3,'interpweight_get_varNdims_file','Problem with nf90_inq_varid', &
            TRIM(nf90_strerror(rcode)), '')
       
       rcode = nf90_inquire_variable(nid, vid, NDIMS = interpweight_get_varNdims_file)
       IF (rcode /= NF90_NOERR) CALL ipslerr_p(3,'interpweight_get_varNdims_file','Problem with nf90_inquire_variable', &
            TRIM(nf90_strerror(rcode)), '')

       rcode = NF90_CLOSE(nid)
       IF (rcode /= NF90_NOERR) CALL ipslerr_p(3,'interpweight_get_varNdims_file','Problem with nf90_close', &
            TRIM(nf90_strerror(rcode)), '')
    END IF
    CALL bcast(interpweight_get_varNdims_file)

  END FUNCTION interpweight_get_varNdims_file

!!================================================================================================================================
!! FUNCTION     : isin_file
!!
!>\BRIEF        ! Function to tell if a given variable is inside a file
!!
!! DESCRIPTION : None
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : isin_file
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_================================================================================================================================
LOGICAL FUNCTION isin_file(filename, varname)

    USE netcdf

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !! 0.1 Input variables
    CHARACTER(LEN=*), INTENT(in)                         :: filename      !! filename: name of the file to open
    CHARACTER(LEN=*), INTENT(in)                         :: varname       !! varname: name of the variable 
    !! 0.2 Modified variables
    !! 0.3 Output variables
    !! 0.4 Local variables
    INTEGER                                              :: nid, vid, Ndims, Nvars
    INTEGER                                              :: iv, rcode
    CHARACTER(LEN=1000)                                  :: varinfile
    CHARACTER(LEN=250)                                   :: msg


    IF (is_root_prc) THEN
       rcode = nf90_open(TRIM(filename), NF90_NOWRITE, nid)
       IF (rcode /= NF90_NOERR) CALL ipslerr_p(3,'isin_file','Problem with nf90_open',TRIM(nf90_strerror(rcode)),'')

       rcode = nf90_inquire(nid, Ndims, Nvars)
       IF (rcode /= NF90_NOERR) CALL ipslerr_p(3,'isin_file','Problem with nf90_inquire',TRIM(nf90_strerror(rcode)),'')

       DO iv=1, Nvars
          rcode = nf90_inquire_variable(nid, iv, name=varinfile)
          IF (rcode /= NF90_NOERR) CALL ipslerr_p(3,'isin_file','Problem with nf90_inquire_variable',TRIM(nf90_strerror(rcode)),'')
          
          IF (TRIM(varinfile) == TRIM(varname)) THEN
             isin_file = .TRUE.
             EXIT
          ELSE
             isin_file = .FALSE.
          END IF
       END DO

       rcode = NF90_CLOSE(nid)
       IF (rcode /= NF90_NOERR) CALL ipslerr_p(3,'isin_file','Problem with nf90_close',TRIM(nf90_strerror(rcode)),'')
    END IF

    CALL bcast(isin_file)

  END FUNCTION isin_file

!!!!!!! !!!!!! !!!!! !!!! !!! !! !
! Generic functions no netCDF

!!================================================================================================================================
!! FUNCTION     : interpweight_Index1DArrayI
!!
!>\BRIEF        ! Function to provide the first index of a given value inside a 1D intger array
!!
!! DESCRIPTION : None
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : interpweight_Index1DArrayI
!!
!! REFERENCE(S) : See above, module description.
!!
!! FLOWCHART    : None
!! \n
!_================================================================================================================================
  INTEGER FUNCTION interpweight_Index1DArrayI(array1D, d1, val)

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !! 0.1 Input variables
    INTEGER, INTENT(in)                                  :: d1            !! size of the array
    INTEGER, INTENT(in)                                  :: val           !! value to search
    INTEGER, DIMENSION(d1), INTENT(in)                   :: array1D       !! values of the array
    !! 0.2 Modified variables
    !! 0.3 Output variables
    !! 0.4 Local variables
    INTEGER                                              :: i


    interpweight_Index1DArrayI = -1

    DO i=1,d1
      IF (array1d(i) == val) THEN
        interpweight_Index1DArrayI = i
        EXIT
      END IF
    END DO

  END FUNCTION interpweight_Index1DArrayI

!!================================================================================================================================
!! FUNCTION     : interpweight_Index1DArrayR
!!
!>\BRIEF        ! Function to provide the first index of a given value inside a 1D real array
!!
!! DESCRIPTION : None
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : interpweight_Index1DArrayR
!!
!! REFERENCE(S) : See above, module description.
!!
!! FLOWCHART    : None
!! \n
!_================================================================================================================================
  INTEGER FUNCTION interpweight_Index1DArrayR(array1D, d1, val)

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !! 0.1 Input variables
    INTEGER, INTENT(in)                                  :: d1            !! size of the array 
    REAL, INTENT(in)                                     :: val           !! value to search
    REAL, DIMENSION(d1), INTENT(in)                      :: array1D       !! values of the array
    !! 0.2 Modified variables
    !! 0.3 Output variables
    !! 0.4 Local variables
    INTEGER                                              :: i

    interpweight_Index1DArrayR = -1

    DO i=1,d1
      IF (array1d(i) == val) THEN
        interpweight_Index1DArrayR = i
        EXIT
      END IF
    END DO

  END FUNCTION interpweight_Index1DArrayR

!!================================================================================================================================
!! FUNCTION     : interpweight_Index2DArrayI
!!
!>\BRIEF        ! Function to provide the first index of a given value inside a 2D integer array
!!
!! DESCRIPTION : None
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : interpweight_Index2DArrayI
!!
!! REFERENCE(S) : See above, module description.
!!
!! FLOWCHART    : None
!! \n
!_================================================================================================================================
  FUNCTION interpweight_Index2DArrayI(array2D, d1, d2, val)

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !! 0.1 Input variables
    INTEGER, INTENT(in)                                  :: d1, d2        !! size of the array 
    INTEGER, INTENT(in)                                  :: val           !! value to search
    INTEGER, DIMENSION(d1,d2), INTENT(in)                :: array2D       !! values of the array
    !! 0.2 Modified variables
    !! 0.3 Output variables
    INTEGER, DIMENSION(2)                                :: interpweight_Index2DArrayI
    !! 0.4 Local variables
    INTEGER                                              :: i, j


    interpweight_Index2DArrayI = -1

    DO i=1,d1
      DO j=1,d2
        IF (array2d(i,j) == val) THEN
          interpweight_Index2DArrayI(1) = i
          interpweight_Index2DArrayI(2) = j
          EXIT
        END IF
      END DO
    END DO

  END FUNCTION interpweight_Index2DArrayI

!!================================================================================================================================
!! FUNCTION     : interpweight_Index2DArrayR
!!
!>\BRIEF        ! Function to provide the first index of a given value inside a 2D real array
!!
!! DESCRIPTION : None
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : bm_vol
!!
!! REFERENCE(S) : See above, module description.
!!
!! FLOWCHART    : None
!! \n
!_================================================================================================================================
  FUNCTION interpweight_Index2DArrayR(array2D, d1, d2, val)

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !! 0.1 Input variables
    INTEGER, INTENT(in)                                  :: d1, d2        !! size of the array 
    REAL, INTENT(in)                                     :: val           !! value to search
    REAL, DIMENSION(d1,d2), INTENT(in)                   :: array2D       !! values of the array
    !! 0.2 Modified variables
    !! 0.3 Output variables
    INTEGER, DIMENSION(2)                                :: interpweight_Index2DArrayR
    !! 0.4 Local variables
    INTEGER                                              :: i, j

    interpweight_Index2DArrayR = -1

    DO i=1,d1
      DO j=1,d2
        IF (array2d(i,j) == val) THEN
          interpweight_Index2DArrayR(1) = i
          interpweight_Index2DArrayR(2) = j
          EXIT
        END IF
      END DO
    END DO

  END FUNCTION interpweight_Index2DArrayR

!!================================================================================================================================
!! FUNCTION     : interpweight_Index1DLonLat
!!
!>\BRIEF        ! Function to provide the first index of a given pair of 
!! longitude and latitude from a set of lon, lat 1D arrays
!!
!! DESCRIPTION : None
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : interpweight_Index1DLonLat
!!
!! REFERENCE(S) : See above, module description.
!!
!! FLOWCHART    : None
!! \n
!_================================================================================================================================
  INTEGER FUNCTION interpweight_Index1DLonLat(lon, lat, d1, lonval, latval)

    USE constantes_var

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !! 0.1 Input variables
    INTEGER, INTENT(in)                                  :: d1            !! size of the array 
    REAL(r_std), INTENT(in)                              :: lonval,latval !! longitude and latitude value
    REAL(r_std), DIMENSION(d1), INTENT(in)               :: lon, lat      !! longitude and latitudes to search in
    !! 0.2 Modified variables
    !! 0.3 Output variables
    !! 0.4 Local variables
    INTEGER                                              :: i
    REAL(r_std)                                          :: mindist
    REAL(r_std), DIMENSION(d1)                           :: dist


    interpweight_Index1DLonLat = -1

    dist = SQRT((lon - lonval)**2. + (lat - latval)**2.)
!    mindist = MINVAL(ABS(dist))
    mindist = 0.

    DO i=1,d1
      IF (dist(i) == mindist) THEN
        interpweight_Index1DLonLat = i
        EXIT
      END IF
    END DO

  END FUNCTION interpweight_Index1DLonLat

!!================================================================================================================================
!! FUNCTION     : interpweight_Index2DLonLat
!!
!>\BRIEF          Function to provide the first index of a given pair of longitude and latitude
!!                from a set of lon, lat 2D arrays
!!
!!
!! DESCRIPTION : None
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : interpweight_Index2DLonLat
!!
!! REFERENCE(S) : See above, module description.
!!
!! FLOWCHART    : None
!! \n
!_================================================================================================================================
  FUNCTION interpweight_Index2DLonLat(lon, lat, d1, d2, lonval, latval)

    USE constantes_var

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !! 0.1 Input variables
    INTEGER, INTENT(in)                                  :: d1, d2        !! size of the array 
    REAL(r_std), INTENT(in)                              :: lonval,latval !! longitude and latitude value
    REAL(r_std), DIMENSION(d1,d2), INTENT(in)            :: lon, lat      !! longitude and latitudes to search in
    !! 0.2 Modified variables
    !! 0.3 Output variables
    INTEGER, DIMENSION(2)                                :: interpweight_Index2DLonLat
    !! 0.4 Local variables
    INTEGER                                              :: i, j
    REAL(r_std)                                          :: mindist
    REAL(r_std), DIMENSION(d1,d2)                        :: dist


    interpweight_Index2DLonLat = -1

    dist = SQRT((lon - lonval)**2. + (lat - latval)**2.)
    mindist = 0.

    DO i=1,d1
      DO j=1,d2
        IF (dist(i,j) == mindist) THEN
          interpweight_Index2DLonLat(1) = i
          interpweight_Index2DLonLat(2) = j
          EXIT
        END IF
      END DO
    END DO

  END FUNCTION interpweight_Index2DLonLat

!!================================================================================================================================
!! FUNCTION     : interpweight_RangeI
!!
!>\BRIEF        ! Function to provide a range of d1 values from 'iniv' to 
!!  'endv', of integer values in a vector
!!
!! DESCRIPTION : None
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : interpweight_RangeI
!!
!! REFERENCE(S) : See above, module description.
!!
!! FLOWCHART    : None
!! \n
!_================================================================================================================================
  FUNCTION interpweight_RangeI(d1, iniv, endv)

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !! 0.1 Input variables
    INTEGER, INTENT(in)                                  :: d1            !! number of values
    INTEGER, INTENT(in)                                  :: iniv          !! initial value
    INTEGER, INTENT(in)                                  :: endv          !! end value
    !! 0.2 Modified variables
    !! 0.3 Output variables
    INTEGER, DIMENSION(d1)                               :: interpweight_RangeI
    !! 0.4 Local variables
    INTEGER                                              :: i, intv

    intv = (endv - iniv) / (d1*1 - 1)

    interpweight_RangeI(1) = iniv
    DO i=2,d1
      interpweight_RangeI(i) = interpweight_RangeI(i-1) + intv
    END DO

  END FUNCTION interpweight_RangeI 

!!================================================================================================================================
!! FUNCTION     : interpweight_RangeR
!!
!>\BRIEF        ! Function to provide a range of d1 values from 'iniv' to 
!!  'endv', of real values in a vector
!!
!! DESCRIPTION : None
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : interpweight_RangeR
!!
!! REFERENCE(S) : See above, module description.
!!
!! FLOWCHART    : None
!! \n
!_================================================================================================================================
  FUNCTION interpweight_RangeR(d1, iniv, endv)

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !! 0.1 Input variables
    INTEGER, INTENT(in)                                  :: d1            !! number of values
    REAL, INTENT(in)                                     :: iniv          !! initial value
    REAL, INTENT(in)                                     :: endv          !! end value
    !! 0.2 Modified variables
    !! 0.3 Output variables
    REAL, DIMENSION(d1)                                  :: interpweight_RangeR
    !! 0.4 Local variables
    INTEGER                                              :: i
    REAL                                                 :: intv

    intv = (endv - iniv) / (d1*1. - 1.)

    interpweight_RangeR(1) = iniv
    DO i=2,d1
      interpweight_RangeR(i) = interpweight_RangeR(i-1) + intv
    END DO

  END FUNCTION interpweight_RangeR

!!================================================================================================================================
!! FUNCTION     : interpweight_ValVecI
!!
!>\BRIEF        ! Function to provide the number of times and where that a 
!!  given value 'oper' on a vector of integers
!!
!! DESCRIPTION : None
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : interpweight_ValVecI
!!
!! REFERENCE(S) : See above, module description.
!!
!! FLOWCHART    : None
!! \n
!_================================================================================================================================
  FUNCTION interpweight_ValVecI(vec, d1, val, oper)

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !! 0.1 Input variables
    INTEGER, INTENT(in)                                  :: d1            !! Number of values
    INTEGER, DIMENSION(d1), INTENT(in)                   :: vec           !! vector of integer values
    INTEGER, INTENT(in)                                  :: val           !! value to search
    CHARACTER(LEN=*), INTENT(in)                         :: oper          !! operation:
                                                                          !!   'eq': equal
                                                                          !!   'ge': greater or equal
                                                                          !!   'le': less or equal
                                                                          !!   'ne': not equal
    !! 0.2 Modified variables
    !! 0.3 Output variables
    INTEGER, DIMENSION(d1)                               :: interpweight_ValVecI !! Vector with positions 
                                                                                 !! (NtimesVal, pos1, pos2, ..., posNtimesVal, ...)
                                                                                 !! NtimesVal=number of correspondencies with 'oper'
    !! 0.4 Local variables
    INTEGER                                              :: i, NtimesVal
    CHARACTER(LEN=50)                                    :: errormsg

    errormsg = 'ERROR -- error -- ERROR -- error'

    NtimesVal = 0
    DO i=1,d1
      SELECT CASE(oper)
        CASE ('eq')
          IF (vec(i) == val) THEN
            NtimesVal = NtimesVal + 1
            IF (NtimesVal < d1) interpweight_ValVecI(NtimesVal + 1) = i
          END IF
        CASE ('ge')
          IF (vec(i) >= val) THEN
            NtimesVal = NtimesVal + 1
            IF (NtimesVal < d1) interpweight_ValVecI(NtimesVal + 1) = i
          END IF
        CASE ('le')
          IF (vec(i) <= val) THEN
            NtimesVal = NtimesVal + 1
            IF (NtimesVal < d1) interpweight_ValVecI(NtimesVal + 1) = i
          END IF
        CASE ('neq')
          IF (vec(i) /= val) THEN
            NtimesVal = NtimesVal + 1
            IF (NtimesVal < d1) interpweight_ValVecI(NtimesVal + 1) = i
          END IF
        CASE DEFAULT
          WRITE(numout,*)TRIM(errormsg)
          WRITE(numout,*)"  interpweight_ValVecI: operation '" // TRIM(oper) // "' not ready!!"
          WRITE(numout,*)"    only available: 'eq', 'ge', 'le', 'neq'"
          CALL ipslerr_p (3, 'interpweight_ValVecI', 'wrong operation', 'provide another one', '')
      END SELECT
    END DO
    interpweight_ValVecI(1) = NtimesVal


  END FUNCTION interpweight_ValVecI

!!================================================================================================================================
!! FUNCTION     : interpweight_ValVecR
!!
!>\BRIEF         Function to provide the number of times and where that a 
!!  given value 'oper' on a vector of reals
!!
!! DESCRIPTION : None
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : interpweight_ValVecR
!!
!! REFERENCE(S) : See above, module description.
!!
!! FLOWCHART    : None
!! \n
!_================================================================================================================================
  FUNCTION interpweight_ValVecR(vec, d1, val, oper)

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !! 0.1 Input variables
    INTEGER, INTENT(in)                                  :: d1            !! Number of values
    REAL, DIMENSION(:), INTENT(in)                       :: vec           !! vector of integer values
    REAL, INTENT(in)                                     :: val           !! value to search
    CHARACTER(LEN=*), INTENT(in)                         :: oper          !! operation:
                                                                          !!   'eq': equal
                                                                          !!   'ge': greater or equal
                                                                          !!   'le': less or equal
                                                                          !!   'neq': not equal
    !! 0.2 Modified variables
    !! 0.3 Output variables
    INTEGER, DIMENSION(d1)                               :: interpweight_ValVecR !! Vector with positions 
                                                                                 !! (NtimesVal, pos1, pos2, ..., posNtimesVal, ...)
                                                                                 !! NtimesVal=number of correspondencies with 'oper'
    !! 0.4 Local variables
    INTEGER                                              :: i, NtimesVal
    CHARACTER(LEN=50)                                    :: errormsg

    errormsg = 'ERROR -- error -- ERROR -- error'

    NtimesVal = 0
    DO i=1,d1
      SELECT CASE(oper)
        CASE ('eq')
          IF (vec(i) == val) THEN
            NtimesVal = NtimesVal + 1
            IF (NtimesVal < d1) interpweight_ValVecR(NtimesVal + 1) = i
          END IF
        CASE ('ge')
          IF (vec(i) >= val) THEN
            NtimesVal = NtimesVal + 1
            IF (NtimesVal < d1) interpweight_ValVecR(NtimesVal + 1) = i
          END IF
        CASE ('le')
          IF (vec(i) <= val) THEN
            NtimesVal = NtimesVal + 1
            IF (NtimesVal < d1) interpweight_ValVecR(NtimesVal + 1) = i
          END IF
        CASE ('neq')
          IF (vec(i) /= val) THEN
            NtimesVal = NtimesVal + 1
            IF (NtimesVal < d1) interpweight_ValVecR(NtimesVal + 1) = i
          END IF
        CASE DEFAULT
          WRITE(numout,*)TRIM(errormsg)
          WRITE(numout,*)"  interpweight_ValVecR: operation '" // TRIM(oper) // "' not ready!!"
          WRITE(numout,*)"    only available: 'eq', 'ge', 'le', 'neq'"
          CALL ipslerr_p (3, 'interpweight_ValVecR', 'wrong operation', 'provide another one', '')
      END SELECT
    END DO
    interpweight_ValVecR(1) = NtimesVal


  END FUNCTION interpweight_ValVecR

!!================================================================================================================================
!! FUNCTION     : interpweight_From2DmatTo1Dvec
!!
!>\BRIEF        ! Suroutine to provide value on a transformation from a 2D 
!!  matrix to a 1D vector
!!    e.g.: from T2(i,j) to T2M(ijklon)
!!
!! DESCRIPTION : None
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : interpweight_From2DmatTo1Dvec
!!
!! REFERENCE(S) : See above, module description.
!!
!! FLOWCHART    : None
!! \n
!_================================================================================================================================
  INTEGER FUNCTION interpweight_From2DmatTo1Dvec(ix,iy,dx,dy)

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !! 0.1 Input variables
    INTEGER, INTENT(in)                                  :: ix,iy         !! location within the 2D matrix
    INTEGER, INTENT(in)                                  :: dx,dy         !! dimension of the 2D matrix

    IF (ix > dx) THEN
      WRITE(numout,*) "Error: "
      WRITE(numout,*)"  interpweight_From2DmatTo1Dvec: 'ix' too large ix > dx ", ix, '>', dx, '!!'
      CALL ipslerr_p (3, 'interpweight_From2DmatTo1Dvec', 'ix value too big', 'out of bounds', '')
    END IF

    IF (iy > dy) THEN
      WRITE(numout,*) "Error: "
      WRITE(numout,*)"  interpweight_From2DmatTo1Dvec: 'iy' too large iy > dy ", iy, '>', dy, '!!'
      CALL ipslerr_p (3, 'interpweight_From2DmatTo1Dvec', 'iy value too big', 'out of bounds', '')
    END IF

    interpweight_From2DmatTo1Dvec = (ix-1)*dy + iy

  END FUNCTION interpweight_From2DmatTo1Dvec

!!================================================================================================================================
!! FUNCTION     : interpweight_From1DvecTo2Dmat
!!
!>\BRIEF        ! Suroutine to provide value on a transformation from a 1D 
!!  vector to a 2D matrix 
!!     e.g.: from T2M(ijklon) to T2(i,j)
!!
!! DESCRIPTION : None
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : interpweight_From1DvecTo2Dmat
!!
!! REFERENCE(S) : See above, module description.
!!
!! FLOWCHART    : None
!! \n
!_================================================================================================================================
  FUNCTION interpweight_From1DvecTo2Dmat(ivec,dx,dy)

    IMPLICIT NONE

    !! 0. Variables and parameter declaration
    !! 0.1 Input variables
    INTEGER, INTENT(in)                                  :: ivec              !! position within the 1D vec
    INTEGER, INTENT(in)                                  :: dx,dy             !! dimension of the 2D matrix
    !! 0.2 Modified variables
    !! 0.3 Output variables
    INTEGER, DIMENSION(2)                                :: interpweight_From1DvecTo2Dmat


    IF (ivec > dx*dy) THEN
      WRITE(numout,*) "Error: "
      WRITE(numout,*)"  interpweight_From1DvecTo2Dmat: 'ivec' too large ivec > dx*dy ", ivec, '>', dx*dy, '!!'
      CALL ipslerr_p (3, 'interpweight_From1DvecTo2Dmat', "'ivec' too large", 'out of bounds!!', '')
    END IF

    interpweight_From1DvecTo2Dmat(1) = INT((ivec-1)/dy) + 1
    interpweight_From1DvecTo2Dmat(2) = ivec - (interpweight_From1DvecTo2Dmat(1)-1)*dy

  END FUNCTION interpweight_From1DvecTo2Dmat

END MODULE interpweight
