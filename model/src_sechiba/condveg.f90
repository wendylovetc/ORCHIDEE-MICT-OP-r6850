! ===============================================================================================================================
! MODULE       : condveg
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        Initialise, compute and update the surface parameters emissivity,
!! roughness and albedo. 
!!
!! \n DESCRIPTION : The module uses 3 settings to control its flow:\n
!! 1. :: rough_dyn to choose between two methods to calculate 
!!    the roughness height. If set to false: the roughness height is calculated by the old formulation 
!!    which does not distinguish between z0m and z0h and which does not vary with LAI
!!    If set to true: the grid average is calculated by the formulation proposed by Su et al. (2001)
!! 2. :: impaze for choosing surface parameters. If set to false, the values for the 
!!    soil albedo, emissivity and roughness height are set to default values which are read from 
!!    the run.def. If set to true, the user imposes its own values, fixed for the grid point. This is useful if 
!!    one performs site simulations, however, 
!!    it is not recommended to do so for spatialized simulations.
!!     roughheight_scal imposes the roughness height in (m) , 
!!	same for emis_scal (in %), albedo_scal (in %), zo_scal (in m)                       
!!     Note that these values are only used if 'impaze' is true.\n
!! 3. :: alb_bare_model for choosing values of bare soil albedo. If set to TRUE bare 
!!    soil albedo depends on soil wetness. If set to FALSE bare soil albedo is the mean 
!!    values of wet and dry soil albedos.\n
!!   The surface fluxes are calculated between two levels: the atmospheric level reference and the effective roughness height 
!! defined as the difference between the mean height of the vegetation and the displacement height (zero wind 
!!    level). Over bare soils, the zero wind level is equal to the soil roughness. Over vegetation, the zero wind level
!!    is increased by the displacement height
!!    which depends on the height of the vegetation. For a grid point composed of different types of vegetation, 
!! an effective surface roughness has to be calculated
!!
!! RECENT CHANGE(S): Added option rough_dyn and subroutine condveg_z0cdrag_dyn. Removed subroutine condveg_z0logz. June 2016.
!!
!! REFERENCES(S)    : None
!!
!! SVN              :
!! $HeadURL$
!! $Date$
!! $Revision$
!! \n
!_ ================================================================================================================================

MODULE condveg
  
  USE ioipsl
  USE xios_orchidee
  !
  ! modules used :
  USE constantes
  USE constantes_soil
  USE pft_parameters
  USE qsat_moisture
  USE interpol_help
  USE mod_orchidee_para
  USE ioipsl_para
  USE sechiba_io_p
  USE grid

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: condveg_main, condveg_initialize, condveg_finalize, condveg_clear 

  !
  ! variables used inside condveg module
  !
  LOGICAL, SAVE                     :: l_first_condveg=.TRUE.           !! To keep first call's trace
!$OMP THREADPRIVATE(l_first_condveg)
  REAL(r_std), ALLOCATABLE, SAVE    :: soilalb_dry(:,:)                 !! Albedo values for the dry bare soil (unitless)
!$OMP THREADPRIVATE(soilalb_dry)
  REAL(r_std), ALLOCATABLE, SAVE    :: soilalb_wet(:,:)                 !! Albedo values for the wet bare soil (unitless)
!$OMP THREADPRIVATE(soilalb_wet)
  REAL(r_std), ALLOCATABLE, SAVE    :: soilalb_moy(:,:)                 !! Albedo values for the mean bare soil (unitless)
!$OMP THREADPRIVATE(soilalb_moy)
  REAL(r_std), ALLOCATABLE, SAVE    :: soilalb_bg(:,:)                  !! Albedo values for the background bare soil (unitless)
!$OMP THREADPRIVATE(soilalb_bg)
  INTEGER, SAVE                     :: printlev_loc                     !! Output debug level
!$OMP THREADPRIVATE(printlev_loc)
  
CONTAINS

  !!  =============================================================================================================================
  !! SUBROUTINE		 		    : condveg_initialize
  !!
  !>\BRIEF			            Allocate module variables, read from restart file or initialize with default values
  !!
  !! DESCRIPTION			    : Allocate module variables, read from restart file or initialize with default values.
  !!                                          condveg_snow is called to initialize corresponding variables.
  !!
  !! RECENT CHANGE(S)			    : None
  !!
  !! MAIN OUTPUT VARIABLE(S)
  !!
  !! REFERENCE(S)			    : None
  !! 
  !! FLOWCHART                              : None
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE condveg_initialize (kjit, kjpindex, index, veget, &
               veget_max, frac_nobio, totfrac_nobio, &
               lalo, neighbours, resolution, contfrac, rest_id, &
               zlev, drysoil_frac, height, snowdz, snowrho, tot_bare_soil, &
               snow, snow_age, snow_nobio, snow_nobio_age, &
               temp_air, pb, u, v, lai, &
               emis, albedo,   z0m, z0h, roughheight, roughheight_pft, &
              frac_snow_veg, frac_snow_nobio)

    !! 0. Variable and parameter declaration

    !! 0.1. Input variables  

    INTEGER(i_std), INTENT(in)                       :: kjit             !! Time step number (unitless)          
    INTEGER(i_std), INTENT(in)                       :: kjpindex         !! Domain size - Number of land pixels  (unitless)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in) :: index            !! Index for the points on the map (unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in):: veget            !! PFT coverage fraction of a PFT (= ind*cn_ind) 
                                                                         !! (m^2 m^{-2}) 
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (in)  :: lalo             !! Geographical coordinates (degree)
    INTEGER(i_std),DIMENSION (kjpindex,NbNeighb), INTENT(in):: neighbours       !! Neighbouring land grid cell
    REAL(r_std), DIMENSION (kjpindex,2), INTENT(in)  :: resolution       !! Size of grid in x and y direction (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)     :: contfrac         !! Fraction of land in each grid box 
    INTEGER(i_std), INTENT(in)                       :: rest_id          !! Restart file identifier 

    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(in) :: veget_max        !! Fraction of vegetation type
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(in) :: frac_nobio    !! Fraction of continental ice, lakes, ...
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)     :: totfrac_nobio    !! total fraction of continental ice+lakes+...
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)    :: zlev             !! Height of first layer
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)     :: drysoil_frac     !! Fraction of visibly Dry soil(between 0 and 1)
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(in):: snowdz           !! Snow depth at each snow layer
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(in):: snowrho          !! Snow density at each snow layer
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(in) :: height           !! Vegetation Height (m)

    REAL(r_std),DIMENSION (kjpindex), INTENT(in)     :: snow             !! Snow mass [Kg/m^2]
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)     :: snow_age         !! Snow age
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(in) :: snow_nobio    !! Snow mass [Kg/m^2] on ice, lakes, ...
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)     :: snow_nobio_age   !! Snow age on ice, lakes, ...
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)    :: tot_bare_soil    !! Total evaporating bare soil fraction 
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)       :: temp_air         !! Air temperature
    REAL(r_std), DIMENSION(kjpindex), INTENT(in)     :: pb               !! Surface pressure (hPa)
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)       :: u                !! Horizontal wind speed, u direction 
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)       :: v                !! Horizontal wind speed, v direction
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(in) :: lai              !! Leaf area index (m2[leaf]/m2[ground])

    !! 0.2. Output variables
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)   :: emis             !! Emissivity
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (out) :: albedo           !! Albedo, vis(1) and nir(2)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)   :: z0m              !! Roughness for momentum (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)   :: z0h              !! Roughness for heat (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)   :: roughheight      !! Effective height for roughness
    REAL(r_std),DIMENSION (kjpindex, nvm), INTENT (out)   :: roughheight_pft 
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)    :: frac_snow_veg    !! Snow cover fraction on vegeted area
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(out):: frac_snow_nobio  !! Snow cover fraction on non-vegeted area
  
    !! 0.3 Modified variables          
  
    !! 0.4 Local variables
    INTEGER(i_std)                                  :: ji, jv           !! Index
    INTEGER(i_std)                                  :: ier              !! Check errors in memory allocation
    CHARACTER(LEN=80)                               :: var_name         !! To store variables names for I/O

    REAL(r_std), DIMENSION(kjpindex,2)               :: albedo_snow      !! Snow albedo for visible and near-infrared range(unitless)
    REAL(r_std), DIMENSION(kjpindex,2)               :: alb_bare         !! Mean bare soil albedo for visible and near-infrared 
                                                                         !! range (unitless) 
    REAL(r_std), DIMENSION(kjpindex,2)               :: alb_veget        !! Mean vegetation albedo for visible and near-infrared 
!                                                                        !! range (unitless) 
!_ ================================================================================================================================
  
    !! 1. Choice of calculation of snow albedo and soil albedo
    IF (.NOT. l_first_condveg) CALL ipslerr_p(3,'condveg_initialize','Error: initialization already done','','')
    l_first_condveg=.FALSE.

    !! Initialize local printlev
    printlev_loc=get_printlev('condveg')    

    IF (printlev>=3) WRITE (numout,*) 'Start condveg_initialize'

    !! 1. Allocate module variables and read from restart or initialize
    IF (alb_bg_modis) THEN
       ! Allocate background soil albedo
       ALLOCATE (soilalb_bg(kjpindex,2),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'condveg_initialize','Pb in allocation for soilalb_bg','','')
       
       ! Read background albedo from restart file
       CALL ioconf_setatt_p('UNITS', '-')
       CALL ioconf_setatt_p('LONG_NAME','Background soil albedo for visible and near-infrared range')
       CALL restget_p (rest_id, 'soilalbedo_bg', nbp_glo, 2, 1, kjit, .TRUE., soilalb_bg, "gather", nbp_glo, index_g)

       ! Initialize by interpolating from file if the variable was not in restart file
       IF ( ALL(soilalb_bg(:,:) == val_exp) ) THEN
          CALL condveg_background_soilalb(kjpindex, lalo, neighbours, resolution, contfrac)
       END IF
       CALL xios_orchidee_send_field("soilalb_bg",soilalb_bg)

    ELSE
       ! Allocate
       ! Dry soil albedo
       ALLOCATE (soilalb_dry(kjpindex,2),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'condveg_initialize','Pb in allocation for soilalb_dry','','')
       
       ! Wet soil albedo
       ALLOCATE (soilalb_wet(kjpindex,2),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'condveg_initialize','Pb in allocation for soilalb_wet','','')
       
       ! Mean soil albedo
       ALLOCATE (soilalb_moy(kjpindex,2),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'condveg_initialize','Pb in allocation for soilalb_moy','','')
       
       ! Read variables from restart file
       ! dry soil albedo
       CALL ioconf_setatt_p('UNITS', '-')
       CALL ioconf_setatt_p('LONG_NAME','Dry bare soil albedo')
       CALL restget_p (rest_id,'soilalbedo_dry' , nbp_glo, 2, 1, kjit, .TRUE., soilalb_dry, "gather", nbp_glo, index_g)
       
       ! wet soil albedo
       CALL ioconf_setatt_p('UNITS', '-')
       CALL ioconf_setatt_p('LONG_NAME','Wet bare soil albedo')
       CALL restget_p (rest_id, 'soilalbedo_wet', nbp_glo, 2, 1, kjit, .TRUE., soilalb_wet, "gather", nbp_glo, index_g)
       
       ! mean soil aledo
       CALL ioconf_setatt_p('UNITS', '-')
       CALL ioconf_setatt_p('LONG_NAME','Mean bare soil albedo')
       CALL restget_p (rest_id, 'soilalbedo_moy', nbp_glo, 2, 1, kjit, .TRUE., soilalb_moy, "gather", nbp_glo, index_g)


       ! Initialize the variables if not found in restart file
       IF ( ALL(soilalb_wet(:,:) == val_exp) .OR. &
            ALL(soilalb_dry(:,:) == val_exp) .OR. &
            ALL(soilalb_moy(:,:) == val_exp)) THEN
          ! One or more of the variables were not in the restart file. 
          ! Call routine condveg_soilalb to calculate them.
          CALL condveg_soilalb(kjpindex, lalo, neighbours, resolution, contfrac)
          WRITE(numout,*) '---> val_exp ', val_exp
          WRITE(numout,*) '---> ALBEDO_wet VIS:', MINVAL(soilalb_wet(:,ivis)), MAXVAL(soilalb_wet(:,ivis))
          WRITE(numout,*) '---> ALBEDO_wet NIR:', MINVAL(soilalb_wet(:,inir)), MAXVAL(soilalb_wet(:,inir))
          WRITE(numout,*) '---> ALBEDO_dry VIS:', MINVAL(soilalb_dry(:,ivis)), MAXVAL(soilalb_dry(:,ivis))
          WRITE(numout,*) '---> ALBEDO_dry NIR:', MINVAL(soilalb_dry(:,inir)), MAXVAL(soilalb_dry(:,inir))
          WRITE(numout,*) '---> ALBEDO_moy VIS:', MINVAL(soilalb_moy(:,ivis)), MAXVAL(soilalb_moy(:,ivis))
          WRITE(numout,*) '---> ALBEDO_moy NIR:', MINVAL(soilalb_moy(:,inir)), MAXVAL(soilalb_moy(:,inir))
       ENDIF
    END IF

    ! z0m
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Roughness for momentum')
    CALL restget_p (rest_id, 'z0m', nbp_glo, 1, 1, kjit, .TRUE., z0m, "gather", nbp_glo, index_g)

    ! z0h
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Roughness for heat')
    CALL restget_p (rest_id, 'z0h', nbp_glo, 1, 1, kjit, .TRUE., z0h, "gather", nbp_glo, index_g)

    ! roughness height
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Roughness height')
    CALL restget_p (rest_id, 'roughheight', nbp_glo, 1, 1, kjit, .TRUE., roughheight, "gather", nbp_glo, index_g)
    CALL restget_p (rest_id, 'roughheight_pft', nbp_glo, nvm, 1, kjit, .TRUE., roughheight_pft, "gather", nbp_glo, index_g)
       

   !! 2. Calculate emissivity
    ! If TRUE read in default values for emissivity
    IF ( impaze ) THEN
       !
       emis(:) = emis_scal
       !
    ! If FALSE set emissivity to 1.
    ELSE
       emis_scal = un
       emis(:) = emis_scal
    ENDIF


    !! 3. Calculate the fraction of snow on vegetation and nobio
    CALL condveg_frac_snow(kjpindex, snow, snow_nobio, snowrho, snowdz, &
                           frac_snow_veg, frac_snow_nobio)

    !! 4. Calculate roughness height if it was not found in the restart file
    IF ( ALL(z0m(:) == val_exp) .OR. ALL(z0h(:) == val_exp) .OR. ALL(roughheight(:) == val_exp)) THEN
       !! Calculate roughness height
       ! Chooses between two methods to calculate the grid average of the roughness.
       ! If impaze set to true:  The grid average is calculated by averaging the drag coefficients over PFT.
       ! If impaze set to false: The grid average is calculated by averaging the logarithm of the roughness length per PFT.
       IF ( impaze ) THEN
          ! Use parameter CONDVEG_Z0 and ROUGHHEIGHT from run.def
          z0m(:) = z0_scal
          z0h(:) = z0_scal
          roughheight(:) = roughheight_scal
          roughheight_pft(:,:) = roughheight_scal
       ELSE
          ! Caluculate roughness height
          IF( rough_dyn ) THEN
             CALL condveg_z0cdrag_dyn(kjpindex, veget, veget_max, frac_nobio, totfrac_nobio, zlev, &
                  &               height, temp_air, pb, u, v, lai, frac_snow_veg, z0m, z0h, roughheight, roughheight_pft)
          ELSE
             CALL condveg_z0cdrag(kjpindex, veget, veget_max, frac_nobio, totfrac_nobio, zlev, &
                  height, tot_bare_soil, frac_snow_veg, z0m, z0h, roughheight, roughheight_pft)
          ENDIF
       END IF
    END IF

    !! 5. Calculate albedo
    CALL condveg_albedo (kjpindex,       veget,         veget_max,     drysoil_frac, frac_nobio,     &
                         totfrac_nobio,  snow,          snow_age,      snow_nobio,     &
                         snow_nobio_age, snowdz,        snowrho,                       &
                         tot_bare_soil,  frac_snow_veg, frac_snow_nobio,               &
                         albedo,         albedo_snow,   alb_bare,     alb_veget)
         
    IF (printlev>=3) WRITE (numout,*) 'condveg_initialize done ' 

  END SUBROUTINE condveg_initialize

!! ==============================================================================================================================
!! SUBROUTINE   : condveg_main
!!
!>\BRIEF        Calls the subroutines to initialise the variables, update the variables
!! and write out data and restart files. 
!!
!!
!! MAIN OUTPUT VARIABLE(S):  emis (emissivity), albedo (albedo of 
!! vegetative PFTs in visible and near-infrared range), z0 (surface roughness height),
!! roughheight (grid effective roughness height), soil type (fraction of soil types) 
!! 
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!!
!! REVISION(S)  : None
!!
!_ ================================================================================================================================

  SUBROUTINE condveg_main (kjit, kjpindex, index,&
       & lalo, neighbours, resolution, contfrac, veget, veget_max, frac_nobio, totfrac_nobio, &
       & zlev, snow, snow_age, snow_nobio, snow_nobio_age, tot_bare_soil, &
         temp_air, pb, u, v, lai, &
       & drysoil_frac, height, snowdz, snowrho, emis, albedo, &
       & frac_snow_veg, frac_snow_nobio, &
       & z0m, z0h, roughheight, roughheight_pft, rest_id, hist_id, hist2_id) 

     !! 0. Variable and parameter declaration

    !! 0.1 Input variables  

    INTEGER(i_std), INTENT(in)                       :: kjit             !! Time step number 
    INTEGER(i_std), INTENT(in)                       :: kjpindex         !! Domain size
    INTEGER(i_std),INTENT (in)                       :: rest_id          !! _Restart_ file identifier
    INTEGER(i_std),INTENT (in)                       :: hist_id          !! _History_ file identifier
    INTEGER(i_std), OPTIONAL, INTENT (in)            :: hist2_id          !! _History_ file 2 identifier
    ! input fields
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in) :: index            !! Indeces of the points on the map
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (in)  :: lalo             !! Geographical coordinates
    INTEGER(i_std),DIMENSION (kjpindex,NbNeighb), INTENT(in):: neighbours!! neighoring grid points if land
    REAL(r_std), DIMENSION (kjpindex,2), INTENT(in)  :: resolution       !! size in x an y of the grid (m)
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)    :: contfrac         ! Fraction of land in each grid box.
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(in) :: veget            !! Fraction of vegetation types
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(in) :: veget_max        !! Fraction of vegetation type
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(in) :: frac_nobio    !! Fraction of continental ice, lakes, ...
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)     :: totfrac_nobio    !! total fraction of continental ice+lakes+...
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)    :: zlev             !! Height of first layer
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)     :: snow             !! Snow mass [Kg/m^2]
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)     :: snow_age         !! Snow age
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(in) :: snow_nobio    !! Snow mass [Kg/m^2] on ice, lakes, ...
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)     :: snow_nobio_age   !! Snow age on ice, lakes, ...
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)     :: drysoil_frac     !! Fraction of visibly Dry soil(between 0 and 1)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(in) :: height           !! Vegetation Height (m)
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(in):: snowdz           !! Snow depth at each snow layer
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(in):: snowrho          !! Snow density at each snow layer
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)    :: tot_bare_soil    !! Total evaporating bare soil fraction 
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)    :: frac_snow_veg    !! Snow cover fraction on vegeted area
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(out):: frac_snow_nobio  !! Snow cover fraction on non-vegeted area
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)       :: temp_air         !! Air temperature
    REAL(r_std), DIMENSION(kjpindex), INTENT(in)     :: pb               !! Surface pressure (hPa)
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)       :: u                !! Horizontal wind speed, u direction 
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)       :: v                !! Horizontal wind speed, v direction
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(in) :: lai              !! Leaf area index (m2[leaf]/m2[ground])

    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (out)   :: emis             !! Emissivity
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (out) :: albedo           !! Albedo, vis(1) and nir(2)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)   :: z0m              !! Roughness for momentum (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)   :: z0h              !! Roughness for heat (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)   :: roughheight      !! Effective height for roughness
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)   :: roughheight_pft    

    !! 0.3 Modified variables

    !! 0.4 Local variables   
    REAL(r_std), DIMENSION(kjpindex,2)               :: albedo_snow      !! Snow albedo (unitless ratio)     
    REAL(r_std), DIMENSION(kjpindex)                 :: albedo_snow_mean !! Mean snow albedo over all wave length, for diag (unitless ratio)     
    REAL(r_std), DIMENSION(kjpindex,2)               :: alb_bare         !! Mean bare soil albedo for visible and near-infrared 
                                                                         !! range (unitless) 
    REAL(r_std), DIMENSION(kjpindex,2)               :: alb_veget        !! Mean vegetation albedo for visible and near-infrared 
                                                                         !! range (unitless) 
    INTEGER(i_std)                                   :: ji
!_ ================================================================================================================================
  
    !! 1. Calculate the fraction of snow on vegetation and nobio
    CALL condveg_frac_snow(kjpindex, snow, snow_nobio, snowrho, snowdz, &
                           frac_snow_veg, frac_snow_nobio)

    !! 2. Calculate emissivity
    emis(:) = emis_scal
    
    !! 3. Calculate roughness height
    
    ! If TRUE read in prescribed values for roughness height
    IF ( impaze ) THEN

       DO ji = 1, kjpindex
         z0m(ji) = z0_scal
         z0h(ji) = z0_scal
         roughheight(ji) = roughheight_scal
         roughheight_pft(ji,:) = roughheight_scal
      ENDDO

    ! Calculate roughness height
    ELSE
     
       IF ( rough_dyn ) THEN
          CALL condveg_z0cdrag_dyn (kjpindex, veget, veget_max, frac_nobio, totfrac_nobio, zlev, height, temp_air, pb, u, v, lai, &
               frac_snow_veg, z0m, z0h, roughheight, roughheight_pft)
       ELSE
          CALL condveg_z0cdrag (kjpindex, veget, veget_max, frac_nobio, totfrac_nobio, zlev, height, tot_bare_soil, &
               frac_snow_veg, z0m, z0h, roughheight,  roughheight_pft)
       ENDIF
     
    ENDIF

    !! 4. Calculate albedo
    CALL condveg_albedo (kjpindex,       veget,         veget_max,     drysoil_frac, frac_nobio,     &
                         totfrac_nobio,  snow,          snow_age,      snow_nobio,     &
                         snow_nobio_age, snowdz,        snowrho,                       &
                         tot_bare_soil,  frac_snow_veg, frac_snow_nobio,               &
                         albedo,         albedo_snow,   alb_bare,      alb_veget)
         


    !! 5. Output diagnostics
    CALL xios_orchidee_send_field("soilalb_vis",alb_bare(:,1))
    CALL xios_orchidee_send_field("soilalb_nir",alb_bare(:,2))
    CALL xios_orchidee_send_field("vegalb_vis",alb_veget(:,1))
    CALL xios_orchidee_send_field("vegalb_nir",alb_veget(:,2))
    CALL xios_orchidee_send_field("albedo_vis",albedo(:,1))
    CALL xios_orchidee_send_field("albedo_nir",albedo(:,2))

    ! Calculcate albedo_snow mean over wave length, setting xios_default_val when there is no snow    
    DO ji=1,kjpindex
       IF (snow(ji) > 0) THEN
          albedo_snow_mean(ji) = (albedo_snow(ji,1) + albedo_snow(ji,2))/2
       ELSE
          albedo_snow_mean(ji) = xios_default_val
       END IF
    END DO
    CALL xios_orchidee_send_field("albedo_snow", albedo_snow_mean)
    
    IF ( almaoutput ) THEN
       CALL histwrite_p(hist_id, 'Albedo', kjit, (albedo(:,1) + albedo(:,2))/2, kjpindex, index)
       CALL histwrite_p(hist_id, 'SAlbedo', kjit, (albedo_snow(:,1) + albedo_snow(:,2))/2, kjpindex, index)
       IF ( hist2_id > 0 ) THEN
          CALL histwrite_p(hist2_id, 'Albedo', kjit, (albedo(:,1) + albedo(:,2))/2, kjpindex, index)
          CALL histwrite_p(hist2_id, 'SAlbedo', kjit, (albedo_snow(:,1) + albedo_snow(:,2))/2, kjpindex, index)
       ENDIF
    ELSE
       CALL histwrite_p(hist_id, 'soilalb_vis', kjit, alb_bare(:,1), kjpindex, index)
       CALL histwrite_p(hist_id, 'soilalb_nir', kjit, alb_bare(:,2), kjpindex, index)
       CALL histwrite_p(hist_id, 'vegalb_vis',  kjit, alb_veget(:,1), kjpindex, index)
       CALL histwrite_p(hist_id, 'vegalb_nir',  kjit, alb_veget(:,2), kjpindex, index)
       IF ( hist2_id > 0 ) THEN
          CALL histwrite_p(hist2_id, 'soilalb_vis', kjit, alb_bare(:,1), kjpindex, index)
          CALL histwrite_p(hist2_id, 'soilalb_nir', kjit, alb_bare(:,2), kjpindex, index)
          CALL histwrite_p(hist2_id, 'vegalb_vis',  kjit, alb_veget(:,1), kjpindex, index)
          CALL histwrite_p(hist2_id, 'vegalb_nir',  kjit, alb_veget(:,2), kjpindex, index)
       ENDIF
    ENDIF

    IF (printlev>=3) WRITE (numout,*)' condveg_main done '

  END SUBROUTINE condveg_main

!!
!=============================================================================================================================
!! SUBROUTINE                             : condveg_finalize
!!
!>\BRIEF                                    Write to restart file
!!
!! DESCRIPTION                            : This subroutine writes the module
!variables and variables calculated in condveg
!!                                          to restart file
!!
!! RECENT CHANGE(S)                       : None
!!
!! REFERENCE(S)                           : None
!!
!! FLOWCHART                              : None
!! \n
!_
!==============================================================================================================================
  SUBROUTINE condveg_finalize (kjit, kjpindex, rest_id, z0m, z0h, roughheight, roughheight_pft)

    !! 0. Variable and parameter declaration
    !! 0.1 Input variables 
    INTEGER(i_std), INTENT(in)                  :: kjit             !! Time step number
    INTEGER(i_std), INTENT(in)                  :: kjpindex         !! Domain size
    INTEGER(i_std),INTENT (in)                  :: rest_id          !! Restart file identifier
    REAL(r_std),DIMENSION(kjpindex), INTENT(in) :: z0m              !! Roughness for momentum
    REAL(r_std),DIMENSION(kjpindex), INTENT(in) :: z0h              !! Roughness for heat
    REAL(r_std),DIMENSION(kjpindex), INTENT(in) :: roughheight      !! Grid effective roughness height (m)     
    REAL(r_std),DIMENSION(kjpindex,nvm), INTENT(in) :: roughheight_pft    
    
    !_ ================================================================================================================================
    
    CALL restput_p (rest_id, 'z0m', nbp_glo, 1, 1, kjit, z0m, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'z0h', nbp_glo, 1, 1, kjit, z0h, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'roughheight', nbp_glo, 1, 1, kjit, roughheight, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'roughheight_pft', nbp_glo, nvm, 1, kjit, roughheight_pft, 'scatter', nbp_glo, index_g)
    
    IF ( alb_bg_modis ) THEN
       CALL restput_p (rest_id, 'soilalbedo_bg', nbp_glo, 2, 1, kjit, soilalb_bg, 'scatter',  nbp_glo, index_g)
    ELSE
       CALL restput_p (rest_id, 'soilalbedo_dry', nbp_glo, 2, 1, kjit, soilalb_dry, 'scatter',  nbp_glo, index_g)
       CALL restput_p (rest_id, 'soilalbedo_wet', nbp_glo, 2, 1, kjit, soilalb_wet, 'scatter',  nbp_glo, index_g)
       CALL restput_p (rest_id, 'soilalbedo_moy', nbp_glo, 2, 1, kjit, soilalb_moy, 'scatter',  nbp_glo, index_g)
    END IF
  END SUBROUTINE condveg_finalize

!! ==============================================================================================================================
!! SUBROUTINE 	: condveg_clear
!!
!>\BRIEF        Deallocate albedo variables
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): None 
!!
!! REFERENCES	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE condveg_clear  ()

      l_first_condveg=.TRUE.
       
      ! Dry soil albedo
       IF (ALLOCATED (soilalb_dry)) DEALLOCATE (soilalb_dry)
       ! Wet soil albedo
       IF (ALLOCATED(soilalb_wet))  DEALLOCATE (soilalb_wet)
       ! Mean soil albedo
       IF (ALLOCATED(soilalb_moy))  DEALLOCATE (soilalb_moy)
       ! BG soil albedo
       IF (ALLOCATED(soilalb_bg))  DEALLOCATE (soilalb_bg)

  END SUBROUTINE condveg_clear

!! ==============================================================================================================================\n
!! SUBROUTINE   : condveg_albedo
!!
!>\BRIEF        Calculate albedo
!!
!! DESCRIPTION  : The albedo is calculated for both the visible and near-infrared 
!! domain. First the mean albedo of the bare soil is calculated. Two options exist: 
!! either the soil albedo depends on soil wetness (drysoil_frac variable), or the soil albedo 
!! is set to a mean soil albedo value.
!! The snow albedo scheme presented below belongs to prognostic albedo 
!! category, i.e. the snow albedo value at a time step depends on the snow albedo value 
!! at the previous time step.
!!
!! First, the following formula (described in Chalita and Treut 1994) is used to describe 
!! the change in snow albedo with snow age on each PFT and each non-vegetative surfaces, 
!! i.e. continental ice, lakes, etc.: \n 
!! \latexonly 
!! \input{SnowAlbedo.tex}
!! \endlatexonly
!! \n
!! Where snowAge is snow age, tcstSnowa is a critical aging time (tcstSnowa=5 days)
!! snowaIni and snowaIni+snowaDec corresponds to albedos measured for aged and
!! fresh snow respectively, and their values for each PFT and each non-vegetative surfaces 
!! is precribed in in constantes_veg.f90.\n
!! In order to estimate gridbox snow albedo, snow albedo values for each PFT and 
!! each  non-vegetative surfaces with a grid box are weightedly summed up by their 
!! respective fractions.\n
!! Secondly, the snow cover fraction is computed as:
!! \latexonly 
!! \input{SnowFraction.tex}
!! \endlatexonly
!! \n
!! Where fracSnow is the fraction of snow on total vegetative or total non-vegetative
!! surfaces, snow is snow mass (kg/m^2) on total vegetated or total nobio surfaces.\n 
!! Finally, the surface albedo is then updated as the weighted sum of fracSnow, total
!! vegetated fraction, total nobio fraction, gridbox snow albedo, and previous
!! time step surface albedo.
!!
!! RECENT CHANGE(S): These calculations were previously done in condveg_albcalc and condveg_snow
!!
!! MAIN OUTPUT VARIABLE(S): :: albedo; surface albedo. :: albedo_snow; snow
!! albedo
!!
!! REFERENCE(S) :  
!! Chalita, S. and H Le Treut (1994), The albedo of temperate and boreal forest and 
!!  the Northern Hemisphere climate: a sensitivity experiment using the LMD GCM,
!!  Climate Dynamics, 10 231-240.
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE condveg_albedo  (kjpindex,       veget,         veget_max,     drysoil_frac, frac_nobio,       &
                              totfrac_nobio,  snow,          snow_age,      snow_nobio,       &
                              snow_nobio_age, snowdz,        snowrho,                         &
                              tot_bare_soil,  frac_snow_veg, frac_snow_nobio,                 &
                              albedo,         albedo_snow,   alb_bare,   alb_veget)

    !! 0. Variable and parameter declarations

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                          :: kjpindex        !! Domain size - Number of land pixels  (unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)   :: veget           !! PFT coverage fraction of a PFT (= ind*cn_ind) 
                                                                           !! (m^2 m^{-2})   
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(in)    :: veget_max
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)        :: drysoil_frac    !! Fraction of visibly Dry soil(between 0 and 1)
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(in) :: frac_nobio      !! Fraction of non-vegetative surfaces, i.e. 
                                                                           !! continental ice, lakes, etc. (unitless)     
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)        :: totfrac_nobio   !! Total fraction of non-vegetative surfaces, i.e. 
                                                                           !! continental ice, lakes, etc. (unitless)   
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)        :: snow            !! Snow mass in vegetation (kg m^{-2})           
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(in) :: snow_nobio      !! Snow mass on continental ice, lakes, etc. (kg m^{-2})      
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)        :: snow_age        !! Snow age (days)        
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(in) :: snow_nobio_age  !! Snow age on continental ice, lakes, etc. (days)    
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(in)   :: snowdz          !! Snow depth at each snow layer
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(in)   :: snowrho         !! Snow density at each snow layer
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)       :: tot_bare_soil   !! Total evaporating bare soil fraction 
    REAL(r_std), DIMENSION(kjpindex), INTENT(in)        :: frac_snow_veg   !! Fraction of snow on vegetation (unitless ratio)
    REAL(r_std), DIMENSION(kjpindex,nnobio), INTENT(in) :: frac_snow_nobio !! Fraction of snow on continental ice, lakes, etc. (unitless ratio) 

    !! 0.2 Output variables
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (out)    :: albedo          !! Albedo (unitless ratio)          
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (out)    :: albedo_snow     !! Snow albedo (unitless ratio)     
    REAL(r_std), DIMENSION(kjpindex,2), INTENT(out)     :: alb_bare        !! Mean bare soil albedo for visible and near-infrared 
                                                                           !! range (unitless). Only calculated for .NOT. impaze
    REAL(r_std), DIMENSION(kjpindex,2), INTENT(out)     :: alb_veget       !! Mean vegetation albedo for visible and near-infrared 
                                                                           !! range (unitless). Only calculated for .NOT. impaze

    !! 0.3 Local variables
    INTEGER(i_std)                                      :: ji, jv, jb,ks   !! indices (unitless)
    REAL(r_std), DIMENSION(kjpindex,2)                  :: snowa_veg       !! Albedo of snow covered area on vegetation
                                                                           !! (unitless ratio)
    REAL(r_std), DIMENSION(kjpindex,nnobio,2)           :: snowa_nobio     !! Albedo of snow covered area on continental ice, 
                                                                           !! lakes, etc. (unitless ratio)     
    REAL(r_std), DIMENSION(kjpindex)                    :: fraction_veg    !! Total vegetation fraction (unitless ratio)
    REAL(r_std), DIMENSION(kjpindex)                    :: agefunc_veg     !! Age dependency of snow albedo on vegetation 
                                                                           !! (unitless)
    REAL(r_std), DIMENSION(kjpindex,nnobio)             :: agefunc_nobio   !! Age dependency of snow albedo on ice, 
                                                                           !! lakes, .. (unitless)
    REAL(r_std)                                         :: alb_nobio       !! Albedo of continental ice, lakes, etc. 
                                                                           !!(unitless ratio)
    REAL(r_std),DIMENSION (nvm,2)                       :: alb_leaf_tmp    !! Variables for albedo values for all PFTs and 
    REAL(r_std),DIMENSION (nvm,2)                       :: snowa_aged_tmp  !! spectral domains (unitless) 
    REAL(r_std),DIMENSION (nvm,2)                       :: snowa_dec_tmp
!_ ================================================================================================================================


    snowa_aged_tmp(:,ivis) = snowa_aged_vis(:)
    snowa_aged_tmp(:,inir) = snowa_aged_nir(:)
    snowa_dec_tmp(:,ivis) = snowa_dec_vis(:)
    snowa_dec_tmp(:,inir) = snowa_dec_nir(:)

    !! 1. Preliminary calculation without considering snow
    snowa_aged_tmp(:,ivis) = snowa_aged_vis(:)
    snowa_aged_tmp(:,inir) = snowa_aged_nir(:)
    snowa_dec_tmp(:,ivis) = snowa_dec_vis(:)
    snowa_dec_tmp(:,inir) = snowa_dec_nir(:)

    IF ( impaze ) THEN
       !! No caluculation, set default value
       albedo(:,ivis) = albedo_scal(ivis)
       albedo(:,inir) = albedo_scal(inir)
       !! need these variables for snow albedo and for diagnostic output
       alb_veget(:,ivis) = albedo_scal(ivis)
       alb_veget(:,inir) = albedo_scal(inir)
       alb_bare(:,ivis) = albedo_scal(ivis)
       alb_bare(:,inir) = albedo_scal(inir)

       ! These variables are needed for snow albedo and for diagnostic output
       alb_veget(:,ivis) = albedo_scal(ivis)
       alb_veget(:,inir) = albedo_scal(inir)
       alb_bare(:,ivis) = albedo_scal(ivis)
       alb_bare(:,inir) = albedo_scal(inir)
    ELSE
       !! Preliminary calculation without considering snow (previously done in condveg_albcalc)    
       ! Assign values of leaf and snow albedo for visible and near-infrared range
       ! to local variable (constantes_veg.f90)
       alb_leaf_tmp(:,ivis) = alb_leaf_vis(:)
       alb_leaf_tmp(:,inir) = alb_leaf_nir(:)
       
       !! 1.1 Calculation and assignment of soil albedo
       
       DO ks = 1, 2! Loop over # of spectra
          
          ! If alb_bg_modis=TRUE, the background soil albedo map for the current simulated month is used
          ! If alb_bg_modis=FALSE and alb_bare_model=TRUE, the soil albedo calculation depends on soil moisture
          ! If alb_bg_modis=FALSE and alb_bare_model=FALSE, the mean soil albedo is used without the dependance on soil moisture
          ! see subroutines 'condveg_soilalb' and 'condveg_background_soilalb'
          IF ( alb_bg_modis ) THEN
             alb_bare(:,ks) = soilalb_bg(:,ks)
          ELSE
             IF ( alb_bare_model ) THEN
                alb_bare(:,ks) = soilalb_wet(:,ks) + drysoil_frac(:) * (soilalb_dry(:,ks) -  soilalb_wet(:,ks))
             ELSE
                alb_bare(:,ks) = soilalb_moy(:,ks)
             ENDIF
          ENDIF
          
          ! Soil albedo is weighed by fraction of bare soil          
          albedo(:,ks) = tot_bare_soil(:) * alb_bare(:,ks)
          
          !! 1.2 Calculation of mean albedo of over the grid cell 
          
          ! Calculation of mean albedo of over the grid cell and
          !    mean albedo of only vegetative PFTs over the grid cell
          alb_veget(:,ks) = zero
          
          DO jv = 2, nvm  ! Loop over # of PFTs
             
             ! Mean albedo of grid cell for visible and near-infrared range
             albedo(:,ks) = albedo(:,ks) + veget(:,jv)*alb_leaf_tmp(jv,ks)
             
             ! Mean albedo of vegetation for visible and near-infrared range
             alb_veget(:,ks) = alb_veget(:,ks) + veget(:,jv)*alb_leaf_tmp(jv,ks)
          ENDDO ! Loop over # of PFTs
          
       ENDDO
    END IF


    !! 2. Calculate snow albedos on both total vegetated and total nobio surfaces
 
    ! The snow albedo could be either prescribed (in condveg_init.f90) or 
    !  calculated following Chalita and Treut (1994).
    ! Check if the precribed value fixed_snow_albedo exists
    IF (ABS(fixed_snow_albedo - undef_sechiba) .GT. EPSILON(undef_sechiba)) THEN
       snowa_veg(:,:) = fixed_snow_albedo
       snowa_nobio(:,:,:) = fixed_snow_albedo
       fraction_veg(:) = un - totfrac_nobio(:)
    ELSE ! calculated following Chalita and Treut (1994)
       
       !! 2.1 Calculate age dependence
       
       ! On vegetated surfaces
       DO ji = 1, kjpindex
          agefunc_veg(ji) = EXP(-snow_age(ji)/tcst_snowa)
       ENDDO
       
       ! On non-vegtative surfaces
       DO jv = 1, nnobio ! Loop over # nobio types
          DO ji = 1, kjpindex
             agefunc_nobio(ji,jv) = EXP(-snow_nobio_age(ji,jv)/tcst_snowa)
          ENDDO
       ENDDO
       
       !! 2.1 Calculate snow albedo 
       ! For vegetated surfaces
       fraction_veg(:) = un - totfrac_nobio(:)
       snowa_veg(:,:) = zero
!
!      Alternative formulation based on veget and not veget_max that needs to be tested 
!      See ticket 223 
!
!       IF (ok_dgvm) THEN
!          DO jb = 1, 2
!             DO ji = 1, kjpindex
!                IF ( fraction_veg(ji) .GT. min_sechiba ) THEN
!                   snowa_veg(ji,jb) = snowa_veg(ji,jb) + &
!                        tot_bare_soil(ji)/fraction_veg(ji) * ( snowa_aged_tmp(1,jb)+snowa_dec_tmp(1,jb)*agefunc_veg(ji) )
!                END IF
!             END DO
!          END DO
!          
!          DO jb = 1, 2
!             DO jv = 2, nvm
!                DO ji = 1, kjpindex
!                   IF ( fraction_veg(ji) .GT. min_sechiba ) THEN
!                      snowa_veg(ji,jb) = snowa_veg(ji,jb) + &
!                           veget(ji,jv)/fraction_veg(ji) * ( snowa_aged_tmp(jv,jb)+snowa_dec_tmp(jv,jb)*agefunc_veg(ji) )
!                   ENDIF
!                ENDDO
!             ENDDO
!          ENDDO
!       ELSE
          DO jb = 1, 2
             DO jv = 1, nvm
                DO ji = 1, kjpindex
                   IF ( fraction_veg(ji) .GT. min_sechiba ) THEN
                      snowa_veg(ji,jb) = snowa_veg(ji,jb) + &
                           veget_max(ji,jv)/fraction_veg(ji) * ( snowa_aged_tmp(jv,jb)+snowa_dec_tmp(jv,jb)*agefunc_veg(ji) )
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
!       ENDIF
       !
       ! snow albedo on other surfaces
       !
       DO jb = 1, 2
          DO jv = 1, nnobio
             DO ji = 1, kjpindex
                snowa_nobio(ji,jv,jb) = ( snowa_aged_tmp(1,jb) + snowa_dec_tmp(1,jb) * agefunc_nobio(ji,jv) ) 
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    
    !! 3. Update surface albedo
    
    ! Update surface albedo using the weighted sum of previous time step surface albedo,
    ! total vegetated fraction, total nobio fraction, snow cover fraction (both vegetated and 
    ! non-vegetative surfaces), and snow albedo (both vegetated and non-vegetative surfaces). 
    ! Although both visible and near-infrared surface albedo are presented, their calculations 
    ! are the same.
    DO jb = 1, 2

       albedo(:,jb) = ( fraction_veg(:) ) * &
            ( (un-frac_snow_veg(:)) * albedo(:,jb) + &
            ( frac_snow_veg(:)  ) * snowa_veg(:,jb)    )
       DO jv = 1, nnobio ! Loop over # nobio surfaces
          
          IF ( jv .EQ. iice ) THEN
             alb_nobio = alb_ice(jb)
          ELSE
             WRITE(numout,*) 'jv=',jv
             WRITE(numout,*) 'DO NOT KNOW ALBEDO OF THIS SURFACE TYPE'
             CALL ipslerr_p(3,'condveg_snow','DO NOT KNOW ALBEDO OF THIS SURFACE TYPE','','')
          ENDIF
          
          albedo(:,jb) = albedo(:,jb) + &
               ( frac_nobio(:,jv) ) * &
               ( (un-frac_snow_nobio(:,jv)) * alb_nobio + &
               ( frac_snow_nobio(:,jv)  ) * snowa_nobio(:,jv,jb)   )
       ENDDO
       
    END DO
    
    ! Calculate snow albedo
    DO jb = 1, 2
       albedo_snow(:,jb) =  fraction_veg(:) * frac_snow_veg(:) * snowa_veg(:,jb)
       DO jv = 1, nnobio 
          albedo_snow(:,jb) = albedo_snow(:,jb) + &
               frac_nobio(:,jv) * frac_snow_nobio(:,jv) * snowa_nobio(:,jv,jb)
       ENDDO
    ENDDO
    
    IF (printlev>=3) WRITE (numout,*) ' condveg_albedo done '
    
  END SUBROUTINE condveg_albedo


  
!! ==============================================================================================================================
!! SUBROUTINE   : condveg_frac_snow
!!
!>\BRIEF        This subroutine calculates the fraction of snow on vegetation and nobio
!!
!! DESCRIPTION  
!!
!! RECENT CHANGE(S): These calculations were previously done in condveg_snow.
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  
  SUBROUTINE condveg_frac_snow(kjpindex, snow, snow_nobio, snowrho, snowdz, &
                               frac_snow_veg, frac_snow_nobio)
    !! 0. Variable and parameter declaration
    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                          :: kjpindex        !! Domain size
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)        :: snow            !! Snow mass in vegetation (kg m^{-2})
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(in) :: snow_nobio      !! Snow mass on continental ice, lakes, etc. (kg m^{-2})
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(in)   :: snowrho         !! Snow density at each snow layer
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(in)   :: snowdz          !! Snow depth at each snow layer

    !! 0.2 Output variables
    REAL(r_std), DIMENSION(kjpindex), INTENT(out)       :: frac_snow_veg   !! Fraction of snow on vegetation (unitless ratio)
    REAL(r_std), DIMENSION(kjpindex,nnobio), INTENT(out):: frac_snow_nobio !! Fraction of snow on continental ice, lakes, etc. 

    !! 0.3 Local variables
    REAL(r_std), DIMENSION(kjpindex)                    :: snowrho_ave     !! Average snow density
    REAL(r_std), DIMENSION(kjpindex)                    :: snowdepth       !! Snow depth
    REAL(r_std), DIMENSION(kjpindex)                    :: snowrho_snowdz       !! Snow rho time snowdz
    INTEGER(i_std)                                      :: jv
   
    !! Calculate snow cover fraction for both total vegetated and total non-vegetative surfaces.
    IF (ok_explicitsnow) THEN
       snowdepth=sum(snowdz,2)
       snowrho_snowdz=sum(snowrho*snowdz,2)
       WHERE(snowdepth(:) .LT. min_sechiba)
          frac_snow_veg(:) = 0.
       ELSEWHERE
          snowrho_ave(:)=snowrho_snowdz(:)/snowdepth(:)
          frac_snow_veg(:) = tanh(snowdepth(:)/(0.025*(snowrho_ave(:)/50.)))
       END WHERE
    ELSE
       frac_snow_veg(:) = MIN(MAX(snow(:),zero)/(MAX(snow(:),zero)+snowcri_alb*sn_dens/100.0),un)
    END IF
    
    DO jv = 1, nnobio
      frac_snow_nobio(:,jv) = MIN(MAX(snow_nobio(:,jv),zero)/(MAX(snow_nobio(:,jv),zero)+snowcri_alb),un)
    ENDDO

    IF (printlev>=3) WRITE (numout,*) ' condveg_frac_snow done '
    
  END SUBROUTINE condveg_frac_snow

  
!! ==============================================================================================================================
!! SUBROUTINE   : condveg_soilalb
!!
!>\BRIEF        This subroutine calculates the albedo of soil (without snow).
!!
!! DESCRIPTION  This subroutine reads the soil colour maps in 1 x 1 deg resolution 
!! from the Henderson-Sellers & Wilson database. These values are interpolated to 
!! the model's resolution and transformed into 
!! dry and wet albedos.\n
!!
!! If the soil albedo is calculated without the dependence of soil moisture, the
!! soil colour values are transformed into mean soil albedo values.\n 
!!
!! The calculations follow the assumption that the grid of the data is regular and 
!! it covers the globe. The calculation for the model grid are based on the borders 
!! of the grid of the resolution.
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S):    soilalb_dry for visible and near-infrared range,
!!                             soilalb_wet for visible and near-infrared range, 
!!                             soilalb_moy for visible and near-infrared range 
!!
!! REFERENCE(S) : 
!! -Wilson, M.F., and A. Henderson-Sellers, 1985: A global archive of land cover and
!!  soils data for use in general circulation climate models. J. Clim., 5, 119-143.
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  
  SUBROUTINE condveg_soilalb(nbpt, lalo, neighbours, resolution, contfrac)
 
    USE interpweight

    IMPLICIT NONE

 
    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                    :: nbpt                  !! Number of points for which the data needs to be 
                                                                           !! interpolated (unitless)             
    REAL(r_std), INTENT(in)                       :: lalo(nbpt,2)          !! Vector of latitude and longitudes (degree)        
    INTEGER(i_std), INTENT(in)                    :: neighbours(nbpt,NbNeighb)!! Vector of neighbours for each grid point 
                                                                           !! (1=N, 2=E, 3=S, 4=W)  
    REAL(r_std), INTENT(in)                       :: resolution(nbpt,2)    !! The size of each grid cell in X and Y (km)
    REAL(r_std), INTENT(in)                       :: contfrac(nbpt)        !! Fraction of land in each grid cell (unitless)   

    !! 0.4 Local variables

    CHARACTER(LEN=80)                             :: filename              !! Filename of soil colour map
    INTEGER(i_std)                                :: i, ib, ip, nbexp      !! Indices
    INTEGER                                       :: ALLOC_ERR             !! Help varialbe to count allocation error
    REAL(r_std), DIMENSION(nbpt)                  :: asoilcol              !! Availability of the soilcol interpolation
    REAL(r_std), DIMENSION(:), ALLOCATABLE        :: variabletypevals      !! Values for all the types of the variable
                                                                           !!   (variabletypevals(1) = -un, not used)
    REAL(r_std), DIMENSION(:,:), ALLOCATABLE      :: soilcolrefrac         !! soilcol fractions re-dimensioned
    REAL(r_std)                                   :: vmin, vmax            !! min/max values to use for the 
                                                                           !!   renormalization
    CHARACTER(LEN=80)                             :: variablename          !! Variable to interpolate
    CHARACTER(LEN=80)                             :: lonname, latname      !! lon, lat names in input file
    CHARACTER(LEN=50)                             :: fractype              !! method of calculation of fraction
                                                                           !!   'XYKindTime': Input values are kinds 
                                                                           !!     of something with a temporal 
                                                                           !!     evolution on the dx*dy matrix'
    LOGICAL                                       :: nonegative            !! whether negative values should be removed
    CHARACTER(LEN=50)                             :: maskingtype           !! Type of masking
                                                                           !!   'nomask': no-mask is applied
                                                                           !!   'mbelow': take values below maskvals(1)
                                                                           !!   'mabove': take values above maskvals(1)
                                                                           !!   'msumrange': take values within 2 ranges;
                                                                           !!      maskvals(2) <= SUM(vals(k)) <= maskvals(1)
                                                                           !!      maskvals(1) < SUM(vals(k)) <= maskvals(3)
                                                                           !!        (normalized by maskvals(3))
                                                                           !!   'var': mask values are taken from a 
                                                                           !!     variable inside the file (>0)
    REAL(r_std), DIMENSION(3)                     :: maskvals              !! values to use to mask (according to 
                                                                           !!   `maskingtype') 
    CHARACTER(LEN=250)                            :: namemaskvar           !! name of the variable to use to mask 
    CHARACTER(LEN=250)                            :: msg
    INTEGER                                       :: fopt
    INTEGER(i_std), DIMENSION(:), ALLOCATABLE     :: vecpos
    INTEGER(i_std), DIMENSION(:), ALLOCATABLE     :: solt

!_ ================================================================================================================================
  !! 1. Open file and allocate memory

  ! Open file with soil colours 

  !Config Key   = SOILALB_FILE
  !Config Desc  = Name of file from which the bare soil albedo
  !Config Def   = soils_param.nc
  !Config If    = NOT(IMPOSE_AZE)
  !Config Help  = The name of the file to be opened to read the soil types from 
  !Config         which we derive then the bare soil albedos. This file is 1x1 
  !Config         deg and based on the soil colors defined by Wilson and Henderson-Seller.
  !Config Units = [FILE]
  !
  filename = 'soils_param.nc'
  CALL getin_p('SOILALB_FILE',filename)


  ALLOCATE(soilcolrefrac(nbpt, classnb), STAT=ALLOC_ERR) 
  IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable soilcolrefrac','','')
  ALLOCATE(vecpos(classnb), STAT=ALLOC_ERR) 
  IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable vecpos','','')
  ALLOCATE(solt(classnb), STAT=ALLOC_ERR)
  IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable solt','','')

! Assigning values to vmin, vmax
  vmin = 1.0
  vmax = classnb

  ALLOCATE(variabletypevals(classnb),STAT=ALLOC_ERR)
  IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variabletypevals','','')
  variabletypevals = -un

  !! Variables for interpweight
  ! Type of calculation of cell fractions
  fractype = 'default'
  ! Name of the longitude and latitude in the input file
  lonname = 'nav_lon'
  latname = 'nav_lat'
  ! Should negative values be set to zero from input file?
  nonegative = .FALSE.
  ! Type of mask to apply to the input data (see header for more details)
  maskingtype = 'mabove'
  ! Values to use for the masking
  maskvals = (/ min_sechiba, undef_sechiba, undef_sechiba /)
  ! Name of the variable with the values for the mask in the input file (only if maskkingtype='var') (here not used)
  namemaskvar = ''

  ! Interpolate variable soilcolor
  variablename = 'soilcolor'
  IF (printlev_loc >= 1) WRITE(numout,*) "condveg_soilalb: Read and interpolate " &
       // TRIM(filename) // " for variable " // TRIM(variablename)
  CALL interpweight_2D(nbpt, classnb, variabletypevals, lalo, resolution, neighbours,          &
    contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,          &
    maskvals, namemaskvar, 0, 0, -1, fractype,                                                        &
    -1., -1., soilcolrefrac, asoilcol)
  IF (printlev_loc >= 5) WRITE(numout,*)'  condveg_soilalb after interpweight_2D'
      
  ! Check how many points with soil information are found
  nbexp = 0

  soilalb_dry(:,:) = zero
  soilalb_wet(:,:) = zero
  soilalb_moy(:,:) = zero
  IF (printlev_loc >= 5) THEN
    WRITE(numout,*)'  condveg_soilalb before starting loop nbpt:', nbpt
    WRITE(numout,*)'  condveg_soilalb initial values classnb: ',classnb
    WRITE(numout,*)'  condveg_soilalb vis_dry. SUM:',SUM(vis_dry),' vis_dry= ',vis_dry
    WRITE(numout,*)'  condveg_soilalb nir_dry. SUM:',SUM(nir_dry),' nir_dry= ',nir_dry
    WRITE(numout,*)'  condveg_soilalb vis_wet. SUM:',SUM(vis_wet),' vis_wet= ',vis_wet
    WRITE(numout,*)'  condveg_soilalb nir_wet. SUM:',SUM(nir_wet),' nir_wet= ',nir_wet
  END IF

  DO ib=1,nbpt ! Loop over domain size

     ! vecpos: List of positions where textures were not zero
     !   vecpos(1): number of not null textures found
     vecpos = interpweight_ValVecR(soilcolrefrac(ib,:),classnb,zero,'neq')
     fopt = vecpos(1)
     IF (fopt == classnb) THEN
        ! All textures are not zero
        solt(:) = (/(i,i=1,classnb)/)
     ELSE IF (fopt == 0) THEN
        WRITE(numout,*)'  condveg_soilalb: for point=', ib, ' no soil class!'
     ELSE
        DO ip = 1,fopt
           solt(ip) = vecpos(ip+1)
        END DO
     END IF
        
     !! 3. Compute the average bare soil albedo parameters
     
     IF ( (fopt .EQ. 0) .OR. (asoilcol(ib) .LT. min_sechiba)) THEN
        ! Initialize with mean value if no points were interpolated or if no data was found
        nbexp = nbexp + 1
        soilalb_dry(ib,ivis) = (SUM(vis_dry)/classnb + SUM(vis_wet)/classnb)/deux
        soilalb_dry(ib,inir) = (SUM(nir_dry)/classnb + SUM(nir_wet)/classnb)/deux
        soilalb_wet(ib,ivis) = (SUM(vis_dry)/classnb + SUM(vis_wet)/classnb)/deux
        soilalb_wet(ib,inir) = (SUM(nir_dry)/classnb + SUM(nir_wet)/classnb)/deux
        soilalb_moy(ib,ivis) = SUM(albsoil_vis)/classnb
        soilalb_moy(ib,inir) = SUM(albsoil_nir)/classnb
     ELSE          
        ! If points were interpolated
        DO ip=1, fopt
           IF ( solt(ip) .LE. classnb) THEN
              ! Set to zero if the value is below min_sechiba
              IF (soilcolrefrac(ib,solt(ip)) < min_sechiba) soilcolrefrac(ib,solt(ip)) = zero

              soilalb_dry(ib,ivis) = soilalb_dry(ib,ivis) + vis_dry(solt(ip))*soilcolrefrac(ib,solt(ip))
              soilalb_dry(ib,inir) = soilalb_dry(ib,inir) + nir_dry(solt(ip))*soilcolrefrac(ib,solt(ip))
              soilalb_wet(ib,ivis) = soilalb_wet(ib,ivis) + vis_wet(solt(ip))*soilcolrefrac(ib,solt(ip))
              soilalb_wet(ib,inir) = soilalb_wet(ib,inir) + nir_wet(solt(ip))*soilcolrefrac(ib,solt(ip))
              soilalb_moy(ib,ivis) = soilalb_moy(ib,ivis) + albsoil_vis(solt(ip))*                    &
                soilcolrefrac(ib,solt(ip))
              soilalb_moy(ib,inir) = soilalb_moy(ib,inir) + albsoil_nir(solt(ip))*                    &
                soilcolrefrac(ib,solt(ip))
           ELSE
              msg = 'The file contains a soil color class which is incompatible with this program'
              CALL ipslerr_p(3,'condveg_soilalb',TRIM(msg),'','')
           ENDIF
        ENDDO
     ENDIF

  ENDDO

  IF ( nbexp .GT. 0 ) THEN
     WRITE(numout,*) 'condveg_soilalb _______'
     WRITE(numout,*) 'condveg_soilalb: The interpolation of the bare soil albedo had ', nbexp
     WRITE(numout,*) 'condveg_soilalb: points without data. This are either coastal points or'
     WRITE(numout,*) 'condveg_soilalb: ice covered land.'
     WRITE(numout,*) 'condveg_soilalb: The problem was solved by using the average of all soils'
     WRITE(numout,*) 'condveg_soilalb: in dry and wet conditions'
     WRITE(numout,*) 'condveg_soilalb: Use the diagnostic output field asoilcol to see location of these points'
  ENDIF

  DEALLOCATE (soilcolrefrac)
  DEALLOCATE (variabletypevals)

  ! Write diagnostics
  CALL xios_orchidee_send_field("asoilcol",asoilcol)


  IF (printlev_loc >= 3) WRITE(numout,*)'  condveg_soilalb ended'

  END SUBROUTINE condveg_soilalb


!! ==============================================================================================================================
!! SUBROUTINE 	: condveg_background_soilalb
!!
!>\BRIEF        This subroutine reads the albedo of bare soil
!!
!! DESCRIPTION  This subroutine reads the background albedo map in 0.5 x 0.5 deg resolution 
!! derived from JRCTIP product to be used as bare soil albedo. These values are then interpolated
!! to the model's resolution.\n
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): soilalb_bg for visible and near-infrared range 
!!
!! REFERENCES	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  
  SUBROUTINE condveg_background_soilalb(nbpt, lalo, neighbours, resolution, contfrac)

    USE interpweight

    IMPLICIT NONE

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                    :: nbpt                  !! Number of points for which the data needs to be 
                                                                           !! interpolated (unitless)             
    REAL(r_std), INTENT(in)                       :: lalo(nbpt,2)          !! Vector of latitude and longitudes (degree)        
    INTEGER(i_std), INTENT(in)                    :: neighbours(nbpt,NbNeighb)!! Vector of neighbours for each grid point 
                                                                           !! (1=N, 2=E, 3=S, 4=W)  
    REAL(r_std), INTENT(in)                       :: resolution(nbpt,2)    !! The size of each grid cell in X and Y (km)
    REAL(r_std), INTENT(in)                       :: contfrac(nbpt)        !! Fraction of land in each grid cell (unitless)   

    !! 0.4 Local variables

    CHARACTER(LEN=80)                             :: filename              !! Filename of background albedo
    REAL(r_std), DIMENSION(nbpt)                  :: aalb_bg               !! Availability of the interpolation
    REAL(r_std), ALLOCATABLE, DIMENSION(:)        :: lat_lu, lon_lu        !! Latitudes and longitudes read from input file
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)      :: lat_rel, lon_rel      !! Help variable to read file data and allocate memory
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)      :: mask_lu               !! Help variable to read file data and allocate memory
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:)   :: mask
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)      :: soilalbedo_bg         !! Help variable to read file data and allocate memory
    INTEGER                                       :: ALLOC_ERR             !! Help varialbe to count allocation error
    REAL(r_std)                                   :: vmin, vmax            !! min/max values to use for the 
                                                                           !!   renormalization
    CHARACTER(LEN=80)                             :: variablename          !! Variable to interpolate
    CHARACTER(LEN=250)                            :: maskvname             !! Variable to read the mask from 
                                                                           !! the file
    CHARACTER(LEN=80)                             :: lonname, latname      !! lon, lat names in input file
    CHARACTER(LEN=50)                             :: fractype              !! method of calculation of fraction
                                                                           !!   'XYKindTime': Input values are kinds 
                                                                           !!     of something with a temporal 
                                                                           !!     evolution on the dx*dy matrix'
    LOGICAL                                       :: nonegative            !! whether negative values should be removed
    CHARACTER(LEN=50)                             :: maskingtype           !! Type of masking
                                                                           !!   'nomask': no-mask is applied
                                                                           !!   'mbelow': take values below maskvals(1)
                                                                           !!   'mabove': take values above maskvals(1)
                                                                           !!   'msumrange': take values within 2 ranges;
                                                                           !!      maskvals(2) <= SUM(vals(k)) <= maskvals(1)
                                                                           !!      maskvals(1) < SUM(vals(k)) <= maskvals(3)
                                                                           !!        (normalized by maskedvals(3))
                                                                           !!   'var': mask values are taken from a 
                                                                           !!     variable inside the file (>0)
    REAL(r_std), DIMENSION(3)                     :: maskvals              !! values to use to mask (according to 
                                                                           !!   `maskingtype') 
    CHARACTER(LEN=250)                            :: namemaskvar           !! name of the variable to use to mask 
    REAL(r_std)                                   :: albbg_norefinf        !! No value
    REAL(r_std), ALLOCATABLE, DIMENSION(:)        :: albbg_default         !! Default value

!_ ================================================================================================================================
  
  !! 1. Open file and allocate memory

  ! Open file with background albedo

  !Config Key   = ALB_BG_FILE
  !Config Desc  = Name of file from which the background albedo is read 
  !Config Def   = alb_bg.nc
  !Config If    = 
  !Config Help  = The name of the file to be opened to read background albedo 
  !Config Units = [FILE]
  !
  filename = 'alb_bg.nc'
  CALL getin_p('ALB_BG_FILE',filename)

 
  ALLOCATE(albbg_default(2), STAT=ALLOC_ERR)
  IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'condveg_background_soilalb','Pb in allocation for albbg_default','','')

! For this case there are not types/categories. We have 'only' a continuos field
! Assigning values to vmin, vmax

  vmin = 0.
  vmax = 9999.

  !! Variables for interpweight
  ! Type of calculation of cell fractions (not used here)
  fractype = 'default'
  ! Name of the longitude and latitude in the input file
  lonname = 'longitude'
  latname = 'latitude'
  ! Default value when no value is get from input file
  albbg_default(ivis) = 0.129
  albbg_default(inir) = 0.247
  ! Reference value when no value is get from input file (not used here)
  albbg_norefinf = undef_sechiba
  ! Should negative values be set to zero from input file?
  nonegative = .FALSE.
  ! Type of mask to apply to the input data (see header for more details)
  maskingtype = 'var'
  ! Values to use for the masking (here not used)
  maskvals = (/ undef_sechiba, undef_sechiba, undef_sechiba /)
  ! Name of the variable with the values for the mask in the input file (only if maskkingtype='var')
  namemaskvar = 'mask'

  ! There is a variable for each chanel 'infrared' and 'visible'
  ! Interpolate variable bg_alb_vis
  variablename = 'bg_alb_vis'
  IF (printlev_loc >= 1) WRITE(numout,*) "condveg_background_soilalb: Read and interpolate " &
       // TRIM(filename) // " for variable " // TRIM(variablename)
  CALL interpweight_2Dcont(nbpt, 0, 0, lalo, resolution, neighbours,                                  &
    contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,          &
    maskvals, namemaskvar, -1, fractype, albbg_default(ivis), albbg_norefinf,                         &
    soilalb_bg(:,ivis), aalb_bg)
  IF (printlev_loc >= 5) WRITE(numout,*)"  condveg_background_soilalb after InterpWeight2Dcont for '" //   &
    TRIM(variablename) // "'"

  ! Interpolate variable bg_alb_nir in the same file
  variablename = 'bg_alb_nir'
  IF (printlev_loc >= 1) WRITE(numout,*) "condveg_background_soilalb: Read and interpolate " &
       // TRIM(filename) // " for variable " // TRIM(variablename)
  CALL interpweight_2Dcont(nbpt, 0, 0, lalo, resolution, neighbours,                                  &
    contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,          &
    maskvals, namemaskvar, -1, fractype, albbg_default(inir), albbg_norefinf,                         &
    soilalb_bg(:,inir), aalb_bg)
  IF (printlev_loc >= 5) WRITE(numout,*)"  condveg_background_soilalb after InterpWeight2Dcont for '" //   &
    TRIM(variablename) // "'"

  IF (ALLOCATED(albbg_default)) DEALLOCATE(albbg_default)

  IF (printlev_loc >= 3) WRITE(numout,*)'  condveg_background_soilalb ended'

  ! Write diagnostics
  CALL xios_orchidee_send_field("aalb_bg",aalb_bg)

  END SUBROUTINE condveg_background_soilalb


!! ==============================================================================================================================
!! SUBROUTINE   : condveg_z0cdrag
!!
!>\BRIEF        Computation of grid average of roughness length by calculating 
!! the drag coefficient.
!!
!! DESCRIPTION  : This routine calculates the mean roughness height and mean 
!! effective roughness height over the grid cell. The mean roughness height (z0) 
!! is computed by averaging the drag coefficients  \n
!!
!! \latexonly 
!! \input{z0cdrag1.tex}
!! \endlatexonly
!! \n 
!!
!! where C is the drag coefficient at the height of the vegetation, kappa is the 
!! von Karman constant, z (Ztmp) is the height at which the fluxes are estimated and z0 the roughness height. 
!! The reference level for z needs to be high enough above the canopy to avoid 
!! singularities of the LOG. This height is set to  minimum 10m above ground. 
!! The drag coefficient increases with roughness height to represent the greater 
!! turbulence generated by rougher surfaces. 
!! The roughenss height is obtained by the inversion of the drag coefficient equation.\n
!!
!! The roughness height for the non-vegetative surfaces is calculated in a second step. 
!! In order to calculate the transfer coefficients the 
!! effective roughness height is calculated. This effective value is the difference
!! between the height of the vegetation and the zero plane displacement height.\nn
!!
!! RECENT CHANGE(S): None
!! 
!! MAIN OUTPUT VARIABLE(S):  :: roughness height(z0) and grid effective roughness height(roughheight)
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE condveg_z0cdrag (kjpindex,veget,veget_max,frac_nobio,totfrac_nobio,zlev, height, tot_bare_soil, frac_snow_veg, &
       &                      z0m, z0h, roughheight, roughheight_pft)

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables
   
    INTEGER(i_std), INTENT(in)                          :: kjpindex      !! Domain size - Number of land pixels  (unitless)
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(in)    :: veget         !! PFT coverage fraction of a PFT (= ind*cn_ind) 
                                                                         !! (m^2 m^{-2})
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(in)    :: veget_max     !! PFT "Maximal" coverage fraction of a PFT 
                                                                         !! (= ind*cn_ind) (m^2 m^{-2})
    REAL(r_std), DIMENSION(kjpindex,nnobio), INTENT(in) :: frac_nobio    !! Fraction of non-vegetative surfaces, 
                                                                         !! i.e. continental ice, lakes, etc. (unitless)
    REAL(r_std), DIMENSION(kjpindex), INTENT(in)        :: totfrac_nobio !! Total fraction of non-vegetative surfaces, 
                                                                         !! i.e. continental ice, lakes, etc. (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: zlev          !! Height of first layer (m)           
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(in)    :: height        !! Vegetation height (m)
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)       :: tot_bare_soil !! Total evaporating bare soil fraction 
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)        :: frac_snow_veg !! Snow cover fraction on vegeted area

    !! 0.2 Output variables

    REAL(r_std), DIMENSION(kjpindex), INTENT(out)       :: z0m           !! Roughness height for momentum (m)
    REAL(r_std), DIMENSION(kjpindex), INTENT(out)       :: z0h           !! Roughness height for heat (m) 
    REAL(r_std), DIMENSION(kjpindex), INTENT(out)       :: roughheight   !! Grid effective roughness height (m) 
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(out)   :: roughheight_pft 
    
    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                      :: jv            !! Loop index over PFTs (unitless)
    REAL(r_std), DIMENSION(kjpindex)                    :: sumveg        !! Fraction of bare soil (unitless)
    REAL(r_std), DIMENSION(kjpindex)                    :: ztmp          !! Max height of the atmospheric level (m)
    REAL(r_std), DIMENSION(kjpindex)                    :: ave_height    !! Average vegetation height (m)
    REAL(r_std), DIMENSION(kjpindex)                    :: d_veg         !! PFT coverage of vegetative PFTs 
                                                                         !! (= ind*cn_ind) (m^2 m^{-2})
    REAL(r_std), DIMENSION(kjpindex)                    :: zhdispl       !! Zero plane displacement height (m)
    REAL(r_std), DIMENSION(kjpindex,nvm)                :: zhdispl_pft 
    REAL(r_std)                                         :: z0_nobio      !! Roughness height of non-vegetative fraction (m),  
                                                                         !! i.e. continental ice, lakes, etc. 
    REAL(r_std), DIMENSION(kjpindex)                    :: dragm         !! Drag coefficient for momentum
    REAL(r_std), DIMENSION(kjpindex)                    :: dragh         !! Drag coefficient for heat
    REAL(r_std), DIMENSION(kjpindex)                    :: z0_ground     !! z0m value used for ground surface
!_ ================================================================================================================================
    
    !! 1. Preliminary calculation

    ! Set maximal height of first layer
    ztmp(:) = MAX(10., zlev(:))

    z0_ground(:) = (1.-frac_snow_veg(:))*z0_bare + frac_snow_veg(:)*z0_bare/10.

    ! Calculate roughness for non-vegetative surfaces
    ! with the von Karman constant 
    dragm(:) = tot_bare_soil(:) * (ct_karman/LOG(ztmp(:)/z0_ground))**2
    dragh(:) = tot_bare_soil(:) * (ct_karman/LOG(ztmp(:)/(z0_ground/ratio_z0m_z0h(1))))*(ct_karman/LOG(ztmp(:)/z0_ground))
    ! Fraction of bare soil
    sumveg(:) = tot_bare_soil(:)

    ! Set average vegetation height to zero
    ave_height(:) = zero
    
    !! 2. Calculate the mean roughness height 
    
    ! Calculate the mean roughness height of
    ! vegetative PFTs over the grid cell
    DO jv = 2, nvm

       ! In the case of forest, use parameter veget_max because 
       ! tree trunks influence the roughness even when there are no leaves
       IF ( is_tree(jv) ) THEN
          ! In the case of grass, use parameter veget because grasses 
          ! only influence the roughness during the growing season
          d_veg(:) = veget_max(:,jv)
       ELSE
          ! grasses only have an influence if they are really there!
          d_veg(:) = veget(:,jv)
       ENDIF
       
       ! Calculate the average roughness over the grid cell:
       ! The unitless drag coefficient is per vegetative PFT
       ! calculated by use of the von Karman constant, the height 
       ! of the first layer and the roughness. The roughness
       ! is calculated as the vegetation height  per PFT 
       ! multiplied by the roughness  parameter 'z0_over_height= 1/16'. 
       ! If this scaled value is lower than 0.01 then the value for 
       ! the roughness of bare soil (0.01) is used. 
       ! The sum over all PFTs gives the average roughness 
       ! per grid cell for the vegetative PFTs.
       dragm(:) = dragm(:) + d_veg(:) * (ct_karman/LOG(ztmp(:)/MAX(height(:,jv)*z0_over_height(jv),z0_ground)))**2
       dragh(:) = dragh(:) + d_veg(:) * (ct_karman/LOG(ztmp(:)/(MAX(height(:,jv)*z0_over_height(jv),z0_ground) / &
            ratio_z0m_z0h(jv)))) * (ct_karman/LOG(ztmp(:)/MAX(height(:,jv)*z0_over_height(jv),z0_ground)))

       ! Sum of bare soil and fraction vegetated fraction
       sumveg(:) = sumveg(:) + d_veg(:)
       
       ! Weigh height of vegetation with maximal cover fraction
       ave_height(:) = ave_height(:) + veget_max(:,jv)*height(:,jv)
       
    ENDDO
    
    !! 3. Calculate the mean roughness height of vegetative PFTs over the grid cell
    
    !  Search for pixels with vegetated part to normalise 
    !  roughness height
    WHERE ( sumveg(:) .GT. min_sechiba ) 
       dragm(:) = dragm(:) / sumveg(:)
       dragh(:) = dragh(:) / sumveg(:)
    ENDWHERE
    ! Calculate fraction of roughness for vegetated part 
    dragm(:) = (un - totfrac_nobio(:)) * dragm(:)
    dragh(:) = (un - totfrac_nobio(:)) * dragh(:)

    DO jv = 1, nnobio ! Loop over # of non-vegative surfaces

       ! Set rougness for ice
       IF ( jv .EQ. iice ) THEN
          z0_nobio = z0_ice
       ELSE
          WRITE(numout,*) 'jv=',jv
          WRITE(numout,*) 'DO NOT KNOW ROUGHNESS OF THIS SURFACE TYPE'
          CALL ipslerr_p(3,'condveg_z0cdrag','DO NOT KNOW ROUGHNESS OF THIS SURFACE TYPE','','')
       ENDIF
       
       ! Sum of vegetative roughness length and non-vegetative
       ! roughness length
       dragm(:) = dragm(:) + frac_nobio(:,jv) * (ct_karman/LOG(ztmp(:)/z0_nobio))**2
       dragh(:) = dragh(:) + frac_nobio(:,jv) * (ct_karman/LOG(ztmp(:)/z0_nobio/ratio_z0m_z0h(1)))*(ct_karman/LOG(ztmp(:)/z0_nobio))

    ENDDO ! Loop over # of non-vegative surfaces
    
    !! 4. Calculate the zero plane displacement height and effective roughness length

    !  Take the exponential of the roughness 
    z0m(:) = ztmp(:) / EXP(ct_karman/SQRT(dragm(:)))
    z0h(:) = ztmp(:) / EXP((ct_karman**2.)/(dragh(:)*LOG(ztmp(:)/z0m(:))))

    ! Compute the zero plane displacement height which
    ! is an equivalent height for the absorption of momentum
    zhdispl(:) = ave_height(:) * height_displacement
    DO jv = 2,nvm
        zhdispl_pft(:,jv) = height(:,jv) * height_displacement
    ENDDO

    ! In order to calculate the fluxes we compute what we call the grid effective roughness height.
    ! This is the height over which the roughness acts. It combines the
    ! zero plane displacement height and the vegetation height.
    roughheight(:) = ave_height(:) - zhdispl(:)
    DO jv = 2,nvm
        roughheight_pft(:,jv) = height(:,jv) - zhdispl_pft(:,jv)
        !!! note that we only use roughheight_pft for croplands
    ENDDO

  END SUBROUTINE condveg_z0cdrag


!! ==============================================================================================================================
!! SUBROUTINE   : condveg_z0cdrag_dyn
!!
!>\BRIEF        Computation of grid average of roughness length by calculating 
!! the drag coefficient based on formulation proposed by Su et al. (2001). 
!!
!! DESCRIPTION  : This routine calculates the mean roughness height and mean 
!! effective roughness height over the grid cell. The mean roughness height (z0) 
!! is computed by averaging the drag coefficients  \n
!!
!! \latexonly 
!! \input{z0cdrag1.tex}
!! \endlatexonly
!! \n 
!!
!! where C is the drag coefficient at the height of the vegetation, kappa is the 
!! von Karman constant, z (Ztmp) is the height at which the fluxes are estimated and z0 the roughness height. 
!! The reference level for z needs to be high enough above the canopy to avoid 
!! singularities of the LOG. This height is set to  minimum 10m above ground. 
!! The drag coefficient increases with roughness height to represent the greater 
!! turbulence generated by rougher surfaces. 
!! The roughenss height is obtained by the inversion of the drag coefficient equation.\n
!! In the formulation of Su et al. (2001), one distinguishes the roughness height for
!! momentum (z0m) and the one for heat (z0h). 
!! z0m is computed as a function of LAI (z0m increases with LAI) and z0h is computed  
!! with a so-called kB-1 term (z0m/z0h=exp(kB-1))
!!
!! RECENT CHANGE(S): Written by N. Vuichard (2016)
!! 
!! MAIN OUTPUT VARIABLE(S):  :: roughness height(z0) and grid effective roughness height(roughheight)
!!
!! REFERENCE(S) : 
!! - Su, Z., Schmugge, T., Kustas, W.P., Massman, W.J., 2001. An Evaluation of Two Models for 
!! Estimation of the Roughness Height for Heat Transfer between the Land Surface and the Atmosphere. J. Appl. 
!! Meteorol. 40, 1933–1951. doi:10.1175/1520-0450(2001)
!! - Ershadi, A., McCabe, M.F., Evans, J.P., Wood, E.F., 2015. Impact of model structure and parameterization 
!! on Penman-Monteith type evaporation models. J. Hydrol. 525, 521–535. doi:10.1016/j.jhydrol.2015.04.008
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE condveg_z0cdrag_dyn (kjpindex,veget,veget_max,frac_nobio,totfrac_nobio,zlev, height, &
       &                      temp_air, pb, u, v, lai, frac_snow_veg, z0m, z0h, roughheight, roughheight_pft)

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables
   
    INTEGER(i_std), INTENT(in)                          :: kjpindex      !! Domain size - Number of land pixels  (unitless)
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(in)    :: veget         !! PFT coverage fraction of a PFT (= ind*cn_ind) 
                                                                         !! (m^2 m^{-2})
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(in)    :: veget_max     !! PFT "Maximal" coverage fraction of a PFT 
                                                                         !! (= ind*cn_ind) (m^2 m^{-2})
    REAL(r_std), DIMENSION(kjpindex,nnobio), INTENT(in) :: frac_nobio    !! Fraction of non-vegetative surfaces, 
                                                                         !! i.e. continental ice, lakes, etc. (unitless)
    REAL(r_std), DIMENSION(kjpindex), INTENT(in)        :: totfrac_nobio !! Total fraction of non-vegetative surfaces, 
                                                                         !! i.e. continental ice, lakes, etc. (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: zlev          !! Height of first layer (m)           
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(in)    :: height        !! Vegetation height (m)    
    REAL(r_std), DIMENSION(kjpindex), INTENT(in)        :: temp_air      !! 2m air temperature (K)
    REAL(r_std), DIMENSION(kjpindex), INTENT(in)        :: pb            !! Surface pressure (hPa)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: u             !! Lowest level wind speed in direction u 
                                                                         !! @tex $(m.s^{-1})$ @endtex 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: v             !! Lowest level wind speed in direction v 
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(in)    :: lai           !! Leaf area index (m2[leaf]/m2[ground])
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)        :: frac_snow_veg    !! Snow cover fraction on vegeted area
    !! 0.2 Output variables

    REAL(r_std), DIMENSION(kjpindex), INTENT(out)       :: z0m           !! Roughness height for momentum (m)
    REAL(r_std), DIMENSION(kjpindex), INTENT(out)       :: z0h           !! Roughness height for heat (m)
    REAL(r_std), DIMENSION(kjpindex), INTENT(out)       :: roughheight   !! Grid effective roughness height (m) 
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(out)   :: roughheight_pft  
    
    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                      :: jv            !! Loop index over PFTs (unitless)
    REAL(r_std), DIMENSION(kjpindex)                    :: sumveg        !! Fraction of bare soil (unitless)
    REAL(r_std), DIMENSION(kjpindex)                    :: ztmp          !! Max height of the atmospheric level (m)
    REAL(r_std), DIMENSION(kjpindex)                    :: ave_height    !! Average vegetation height (m)
    REAL(r_std), DIMENSION(kjpindex)                    :: zhdispl       !! Zero plane displacement height (m)
    REAL(r_std), DIMENSION(kjpindex,nvm)                :: zhdispl_pft    
    REAL(r_std)                                         :: z0_nobio      !! Roughness height of non-vegetative fraction (m),  
                                                                         !! i.e. continental ice, lakes, etc. 
    REAL(r_std), DIMENSION(kjpindex)                    :: z0m_pft       !! Roughness height for momentum for a specific PFT
    REAL(r_std), DIMENSION(kjpindex)                    :: z0h_pft       !! Roughness height for heat for a specific PFT
    REAL(r_std), DIMENSION(kjpindex)                    :: dragm         !! Drag coefficient for momentum
    REAL(r_std), DIMENSION(kjpindex)                    :: dragh         !! Drag coefficient for heat
    REAL(r_std), DIMENSION(kjpindex)                    :: eta           !! Ratio of friction velocity to the wind speed at the canopy top - See Ershadi et al. (2015)
    REAL(r_std), DIMENSION(kjpindex)                    :: eta_ec        !! Within-canopy wind speed profile estimation coefficient - See Ershadi et al. (2015)
    REAL(r_std), DIMENSION(kjpindex)                    :: Ct_star       !! Heat transfer coefficient of the soil - see Su et al. (2001)
    REAL(r_std), DIMENSION(kjpindex)                    :: kBs_m1        !! Canopy model of Brutsaert (1982) for a bare soil surface - used in the calculation of kB_m1 (see Ershadi et al. (2015))
    REAL(r_std), DIMENSION(kjpindex)                    :: kB_m1         !! kB**-1: Term used in the calculation of z0h where B-1 is the inverse Stanton number (see Ershadi et al. (2015))
    REAL(r_std), DIMENSION(kjpindex)                    :: fc            !! fractional canopy coverage
    REAL(r_std), DIMENSION(kjpindex)                    :: fs            !! fractional soil coverage
    REAL(r_std), DIMENSION(kjpindex)                    :: Reynolds      !! Reynolds number
    REAL(r_std), DIMENSION(kjpindex)                    :: wind          !! wind Speed (m)
    REAL(r_std), DIMENSION(kjpindex)                    :: u_star        !! friction velocity
    REAL(r_std), DIMENSION(kjpindex)                    :: z0_ground     !! z0m value used for ground surface
!_ ================================================================================================================================
    
    !! 1. Preliminary calculation

    ! Set maximal height of first layer
    ztmp(:) = MAX(10., zlev(:))
    
    z0_ground(:) = (1.-frac_snow_veg(:))*z0_bare + frac_snow_veg(:)*z0_bare/10.

    ! Calculate roughness for non-vegetative surfaces
    ! with the von Karman constant 
    dragm(:) = veget_max(:,1) * (ct_karman/LOG(ztmp(:)/z0_ground(:)))**2

    wind(:) = SQRT(u(:)*u(:)+v(:)*v(:))
    u_star(:)= ct_karman * MAX(min_wind,wind(:)) / LOG(zlev(:)/z0_ground(:))
    Reynolds(:) = z0_ground(:) * u_star(:) &
         / (1.327*1e-5 * (pb_std/pb(:)) * (temp_air(:)/ZeroCelsius)**(1.81))
    
    kBs_m1(:) = 2.46 * reynolds**(1./4.) - LOG(7.4)

    dragh(:) = veget_max(:,1) * (ct_karman/LOG(ztmp(:)/z0_ground(:)))*(ct_karman/LOG(ztmp(:)/(z0_ground(:)/ exp(kBs_m1(:))) ))

    ! Fraction of bare soil
    sumveg(:) = veget_max(:,1)

    ! Set average vegetation height to zero
    ave_height(:) = zero
    
    !! 2. Calculate the mean roughness height 
    
    ! Calculate the mean roughness height of
    ! vegetative PFTs over the grid cell
    DO jv = 2, nvm
       
       WHERE(veget_max(:,jv) .GT. zero)       
          ! Calculate the average roughness over the grid cell:
          ! The unitless drag coefficient is per vegetative PFT
          ! calculated by use of the von Karman constant, the height 
          ! of the first layer and the roughness. The roughness
          ! is calculated as the vegetation height  per PFT 
          ! multiplied by the roughness  parameter 'z0_over_height= 1/16'. 
          ! If this scaled value is lower than 0.01 then the value for 
          ! the roughness of bare soil (0.01) is used. 
          ! The sum over all PFTs gives the average roughness 
          ! per grid cell for the vegetative PFTs.
          eta(:) = c1 - c2 * exp(-c3 * Cdrag_foliage * lai(:,jv))
          
          z0m_pft(:) = (height(:,jv)*(1-height_displacement)*(exp(-ct_karman/eta(:))-exp(-ct_karman/(c1-c2)))) &
               + z0_ground(:)
   
          dragm(:) = dragm(:) + veget_max(:,jv) * (ct_karman/LOG(ztmp(:)/z0m_pft(:)))**2
   
          fc(:) = veget(:,jv)/veget_max(:,jv)
          fs(:) = 1. - fc(:)

          eta_ec(:) = ( Cdrag_foliage * lai(:,jv)) / (2 * eta(:)*eta(:))
          wind(:) = SQRT(u(:)*u(:)+v(:)*v(:))
          u_star(:)= ct_karman * MAX(min_wind,wind(:)) / LOG((zlev(:)+(height(:,jv)*(1-height_displacement)))/z0m_pft(:))
          Reynolds(:) = z0_ground(:) * u_star(:) &
               / (1.327*1e-5 * (pb_std/pb(:)) * (temp_air(:)/ZeroCelsius)**(1.81))
                 
          kBs_m1(:) = 2.46 * reynolds**(1./4.) - LOG(7.4)
          Ct_star(:) = Prandtl**(-2./3.) * SQRT(1./Reynolds(:))
   
          WHERE(lai(:,jv) .GT. min_sechiba)
             kB_m1(:) = (ct_karman * Cdrag_foliage) / (4 * Ct * eta(:) * (1 - exp(-eta_ec(:)/2.))) * fc(:)**2. &
                  + 2*fc(:)*fs(:) * (ct_karman * eta(:) * z0m_pft(:) / height(:,jv)) / Ct_star(:) &
                  + kBs_m1(:) * fs(:)**2. 
          ELSEWHERE
             kB_m1(:) = kBs_m1(:) * fs(:)**2. 
          ENDWHERE
   
          z0h_pft(:) = z0m_pft(:) / exp(kB_m1(:))
   
          dragh(:) = dragh(:) + veget_max(:,jv) * (ct_karman/LOG(ztmp(:)/z0m_pft(:)))*(ct_karman/LOG(ztmp(:)/z0h_pft(:)))
   
          ! Sum of bare soil and fraction vegetated fraction
          sumveg(:) = sumveg(:) + veget_max(:,jv)

          ! Weigh height of vegetation with maximal cover fraction
          ave_height(:) = ave_height(:) + veget_max(:,jv)*height(:,jv)

       ENDWHERE
    ENDDO
    
    !! 3. Calculate the mean roughness height of vegetative PFTs over the grid cell
    
    !  Search for pixels with vegetated part to normalise 
    !  roughness height
    WHERE ( sumveg(:) .GT. min_sechiba ) 
       dragh(:) = dragh(:) / sumveg(:)
       dragm(:) = dragm(:) / sumveg(:)
    ENDWHERE

    ! Calculate fraction of roughness for vegetated part 
    dragh(:) = (un - totfrac_nobio(:)) * dragh(:)
    dragm(:) = (un - totfrac_nobio(:)) * dragm(:)

    DO jv = 1, nnobio ! Loop over # of non-vegative surfaces

       ! Set rougness for ice
       IF ( jv .EQ. iice ) THEN
          z0_nobio = z0_ice
       ELSE
          WRITE(numout,*) 'jv=',jv
          WRITE(numout,*) 'DO NOT KNOW ROUGHNESS OF THIS SURFACE TYPE'
          CALL ipslerr_p(3,'condveg_z0cdrag_dyn','DO NOT KNOW ROUGHNESS OF THIS SURFACE TYPE','','')
       ENDIF
       
       ! Sum of vegetative roughness length and non-vegetative roughness length
       ! Note that z0m could be made dependent of frac_snow_nobio
       dragm(:) = dragm(:) + frac_nobio(:,jv) * (ct_karman/LOG(ztmp(:)/z0_nobio))**2
       
       u_star(:)= ct_karman * MAX(min_wind,wind(:)) / LOG(zlev(:)/z0_nobio)
       Reynolds(:) = z0_nobio * u_star(:) &
            / (1.327*1e-5 * (pb_std/pb(:)) * (temp_air(:)/ZeroCelsius)**(1.81))
       
       kBs_m1(:) = 2.46 * reynolds**(1./4.) - LOG(7.4)
   
       dragh(:) = dragh(:) + frac_nobio(:,jv) *  (ct_karman/LOG(ztmp(:)/z0_nobio)) * &
            (ct_karman/LOG(ztmp(:)/(z0_nobio/ exp(kBs_m1(:))) ))
    ENDDO ! Loop over # of non-vegative surfaces
    
    !! 4. Calculate the zero plane displacement height and effective roughness length
    !  Take the exponential of the roughness 
    z0m(:) = ztmp(:) / EXP(ct_karman/SQRT(dragm(:)))
    z0h(:) = ztmp(:) / EXP((ct_karman**2.)/(dragh(:)*LOG(ztmp(:)/z0m(:))))

    ! Compute the zero plane displacement height which
    ! is an equivalent height for the absorption of momentum
    zhdispl(:) = ave_height(:) * height_displacement
    DO jv = 2,nvm
        zhdispl_pft(:,jv) = height(:,jv) * height_displacement
    ENDDO

    ! In order to calculate the fluxes we compute what we call the grid effective roughness height.
    ! This is the height over which the roughness acts. It combines the
    ! zero plane displacement height and the vegetation height.
    roughheight(:) = ave_height(:) - zhdispl(:)
    DO jv = 2,nvm
        roughheight_pft(:,jv) = height(:,jv) - zhdispl_pft(:,jv)
        !!!! roughheight_pft only used for croplands, xuhui
    ENDDO

  END SUBROUTINE condveg_z0cdrag_dyn


END MODULE condveg
