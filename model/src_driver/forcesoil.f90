! =================================================================================================================================
! PROGRAM       : forcesoil
!
! CONTACT	: orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      	: IPSL (2006). This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        This subroutine runs the soilcarbon submodel using specific initial conditions 
!! and driving variables in order to obtain soil carbon stocks closed to the steady-state values 
!! quicker than when using the ''full'' ORCHIDEE.  
!!	
!!\n DESCRIPTION: None
!! This subroutine computes the soil carbon stocks by calling the soilcarbon routine at each time step. \n
!! The aim is to obtain soil carbon stocks closed to the steady-state values and ultimately to create  
!! an updated stomate restart file for the stomate component. The state variables of the subsystem are the clay content 
!! (fixed value) and the soil carbon stocks. Initial conditions for the state variables are read in an  
!! input stomate restart file. Driving variables are Soil carbon input, Water and Temperature stresses on 
!! Organic Matter decomposition. Driving variables are read from a specific forcing file produced by a former run of ORCHIDEE
!! (SECHIBA+STOMATE). \n 
!! The FORCESOIL program first consists in reading a set of input files, allocating variables and 
!! preparing output stomate restart file. \n                                                             
!! Then, a loop over time is performed in which the soilcarbon routine is called at each time step. \n
!! Last, final values of the soil carbon stocks are written into the output stomate restart file. \n
!! No flag is associated with the use of the FORCESOIL program. \n
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S) : None    
!!
!! FLOWCHART    : None
!!
!! SVN		:
!! $HeadURL: $ 
!! $Date: $
!! $Revision: $
!! \n
!_ =================================================================================================================================

PROGRAM forcesoil
 
  USE netcdf
  !-
  USE utils
  USE defprec
  USE constantes
  USE constantes_var
  USE constantes_mtc
  USE constantes_soil
  USE pft_parameters 
  USE stomate_data
  USE ioipsl_para
  USE mod_orchidee_para
  USE stomate_soilcarbon
  USE stomate_permafrost_soilcarbon
  USE constantes_soil_var
  USE stomate
  USE stomate_io_carbon_permafrost
#ifdef CPP_PARA
  USE mpi
#endif
  !-
  IMPLICIT NONE
  !-
  !-
  CHARACTER(LEN=80)                          :: sto_restname_in,sto_restname_out
  INTEGER(i_std)                             :: iim,jjm                !! Indices (unitless)

  INTEGER(i_std),PARAMETER                   :: llm = 1                !! Vertical Layers (requested by restini routine) (unitless)
  INTEGER(i_std)                             :: kjpindex               !! Domain size (unitless)

  INTEGER(i_std)                             :: itau_dep,itau_len      !! Time step read in the restart file (?) 
                                                                       !! and number of time steps of the simulation (unitless) 
  CHARACTER(LEN=30)                          :: time_str               !! Length of the simulation (year)
  REAL(r_std)                                :: dt_files               !! time step between two successive itaus (?) 
                                                                       !! (requested by restini routine) (seconds)
  REAL(r_std)                                :: date0                  !! Time at which itau = 0 (requested by restini routine) (?)
  INTEGER(i_std)                             :: rest_id_sto            !! ID of the input restart file (unitless)
  CHARACTER(LEN=20), SAVE                    :: thecalendar = 'noleap' !! Type of calendar defined in the input restart file 
                                                                       !! (unitless)
  !-
  CHARACTER(LEN=100)                         :: Cforcing_name          !! Name of the forcing file (unitless)
  INTEGER                                    :: Cforcing_id            !! ID of the forcing file (unitless)
  INTEGER                                    :: v_id                   !! ID of the variable 'Index' stored in the forcing file 
                                                                       !! (unitless)
  REAL(r_std)                                :: dt_forcesoil           !! Time step at which soilcarbon routine is called (days) 
  INTEGER                                    :: nparan                 !! Number of values stored per year in the forcing file 
                                                                       !! (unitless)
  INTEGER                                    :: nbyear
  INTEGER(i_std),DIMENSION(:),ALLOCATABLE    :: indices                !! Grid Point Index used per processor (unitless)
  INTEGER(i_std),DIMENSION(:),ALLOCATABLE    :: indices_g              !! Grid Point Index for all processor (unitless)
  REAL(r_std),DIMENSION(:),ALLOCATABLE       :: x_indices_g            !! Grid Point Index for all processor (unitless)
  REAL(r_std),DIMENSION(:,:),ALLOCATABLE     :: lon, lat               !! Longitude and Latitude of each grid point defined 
                                                                       !! in lat/lon (2D) (degrees)
  REAL(r_std),DIMENSION(llm)                 :: lev                    !! Number of level (requested by restini routine) (unitless)


  INTEGER                                    :: i,m,iatt,iv,iyear      !! counters (unitless)
                                                                       
  CHARACTER(LEN=100)                          :: var_name              
  CHARACTER(LEN=8000)                        :: taboo_vars             !! string used for storing the name of the variables 
                                                                       !! of the stomate restart file that are not automatically 
                                                                       !! duplicated from input to output restart file (unitless)
  REAL(r_std),DIMENSION(1)                   :: xtmp                   !! scalar read/written in restget/restput routines (unitless)
  INTEGER(i_std),PARAMETER                   :: nbvarmax=1000           !! maximum # of variables assumed in the stomate restart file 
                                                                       !! (unitless)
  INTEGER(i_std)                             :: nbvar                  !! # of variables effectively present 
                                                                       !! in the stomate restart file (unitless)
  CHARACTER(LEN=1000),DIMENSION(nbvarmax)      :: varnames              !! list of the names of the variables stored 
                                                                       !! in the stomate restart file (unitless)
  INTEGER(i_std)                             :: varnbdim               !! # of dimensions of a given variable 
                                                                       !! of the stomate restart file
  INTEGER(i_std),PARAMETER                   :: varnbdim_max=20        !! maximal # of dimensions assumed for any variable 
                                                                       !! of the stomate restart file 
  INTEGER,DIMENSION(varnbdim_max)            :: vardims                !! length of each dimension of a given variable 
                                                                       !! of the stomate restart file
  LOGICAL                                    :: l1d                    !! boolean : TRUE if all dimensions of a given variable 
                                                                       !! of the stomate restart file are of length 1 (ie scalar) 
                                                                       !! (unitless)
  REAL(r_std),DIMENSION(:),ALLOCATABLE         :: var_2d               !! matrix read/written in restget/restput routines (unitless)
  REAL(r_std),DIMENSION(:,:),ALLOCATABLE       :: var_3d               !! matrix read/written in restget/restput routines (unitless)
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE     :: var_4d               !! matrix read/written in restget/restput routines (unitless)
  REAL(r_std),DIMENSION(:,:,:,:),ALLOCATABLE   :: var_5d               !! matrix read/written in restget/restput routines (unitless)
  REAL(r_std),DIMENSION(:,:,:,:,:),ALLOCATABLE :: var_6d               !! matrix read/written in restget/restput routines (unitless)
  REAL(r_std)                                :: x_tmp                  !! temporary variable used to store return value 
  INTEGER(i_std)                               :: orch_vardims         !! Orchidee dimensions (different to IOIPSL -exclude time dim-)
                                                                       !! from nf90_get_att (unitless)
  CHARACTER(LEN=10)  :: part_str                                       !! string suffix indicating the index of a PFT 
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)  :: clay_g                 !! clay fraction (nbpglo) (unitless)
                                                                       !! (\f$gC m^{-2} dt_forcesoil^{-1}\f$) 
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE   :: control_temp_g         !! Temperature control (nbp_glo,above/below,time) on OM decomposition 
                                                                       !! (unitless)
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE   :: control_moist_g        !! Moisture control (nbp_glo,abo/below,time) on OM decomposition 
                                                                       !! ?? Should be defined per PFT as well (unitless)
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE   :: carbon_g               !! Soil carbon stocks (nbp_glo,ncarb,nvm) (\f$gC m^{-2}\f$)
                                                                       
  REAL(r_std),ALLOCATABLE :: clay(:)                                   !! clay fraction (nbp_loc) (unitless)
  REAL(r_std),ALLOCATABLE :: soilcarbon_input(:,:,:,:)                 !! soil carbon input (nbp_loc,ncarb,nvm,time) 
                                                                       !! (\f$gC m^{-2} dt_forcesoil^{-1}\f$) 
  REAL(r_std),ALLOCATABLE :: control_temp(:,:,:)                       !! Temperature control (nbp_loc,above/below,time) on OM decomposition 
                                                                       !! (unitless)
  REAL(r_std),ALLOCATABLE :: control_moist(:,:,:)                      !! Moisture control (nbp_loc,abo/below,time) on OM decomposition 
                                                                       !! ?? Should be defined per PFT as well (unitless)
  REAL(r_std),ALLOCATABLE :: carbon(:,:,:)                             !! Soil carbon stocks (nbp_loc,ncarb,nvm) (\f$gC m^{-2}\f$)
  REAL(r_std),DIMENSION(:,:),ALLOCATABLE     :: resp_hetero_soil       !! Heterotrophic respiration (\f$gC m^{-2} dt_forcesoil^{-1}\f$) 
                                                                       !! (requested by soilcarbon routine but not used here) 

  INTEGER(i_std)                             :: printlev_loc           !! Local write level                                                                     
  INTEGER(i_std)                             :: ier,iret               !! Used for hangling errors 
                                                                       
  CHARACTER(LEN=50) :: temp_name                                       
  CHARACTER(LEN=100) :: msg3                                       
  LOGICAL :: debug                                                     !! boolean used for printing messages
  LOGICAL :: l_error                                                   !! boolean for memory allocation
  REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:,:) :: matrixA              !! Carbon fluxes matrix 
  ! allocateable arrays needed for permafrost carbon 
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE  :: deepC_a_g               !! active carbon concentration
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE  :: deepC_s_g               !! slow carbon concentration  
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE  :: deepC_p_g               !! passive carbon concentration
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE  :: O2_soil_g               !! oxygen in the soil
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE  :: CH4_soil_g              !! methane in the soil
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE  :: O2_snow_g               !! oxygen in the snow
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE  :: CH4_snow_g              !! methane in the snow 
  REAL(r_std),DIMENSION(:),ALLOCATABLE  :: z_organic_g                 !! organic carbon depth
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE :: snowdz_g                 !! snow depth at each layer
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE :: snowrho_g                !! snow density at each layer
  REAL(r_std),DIMENSION(:),ALLOCATABLE  :: zz_deep                     !! deep vertical profile 
  REAL(r_std),DIMENSION(:),ALLOCATABLE  :: zz_coef_deep                !! deep vertical profile
  REAL(r_std),DIMENSION(:),ALLOCATABLE  :: z_organic                   !! organic carbon depth
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE  :: deepC_a                 !! active carbon concentration
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE  :: deepC_s                 !! slow carbon concentration
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE  :: deepC_p                 !! passive carbon concentration
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE  :: O2_soil                 !! oxygen in the soil
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE  :: CH4_soil                !! methane in the soil
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE  :: O2_snow                 !! oxygen in the snow
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE  :: CH4_snow                !! methane in the snow 
  REAL(r_std),DIMENSION(:,:),ALLOCATABLE  :: pb                        !! surface pressure
  REAL(r_std),DIMENSION(:,:),ALLOCATABLE  :: snow                      !! snow mass
  REAL(r_std),DIMENSION(:,:,:,:),ALLOCATABLE  :: tprof                 !! deep soil temperature profile
  REAL(r_std),DIMENSION(:,:,:,:),ALLOCATABLE  :: fbact                 !! factor for soil carbon decomposition
  REAL(r_std),DIMENSION(:,:,:,:),ALLOCATABLE  :: hslong                !! deep soil humidity
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE  :: veget_max               !! maximum vegetation fraction
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE  :: rprof                   !! PFT rooting depth
  REAL(r_std),DIMENSION(:,:),ALLOCATABLE  :: tsurf                     !! surface temperature
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE :: snowrho                  !! snow density
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE :: snowdz                   !! snow depth
  REAL(r_std),DIMENSION(:,:),ALLOCATABLE      :: lalo                  !! Geogr. coordinates (latitude,longitude) (degrees)    
  REAL(r_std), DIMENSION(:,:,:), ALLOCATABLE    :: heat_Zimov                 !! heating associated with decomposition  [W/m**3 soil]
  REAL(R_STD), DIMENSION(:), ALLOCATABLE      :: sfluxCH4_deep, sfluxCO2_deep !! [g / m**2]
  REAL(R_STD), DIMENSION(:,:), ALLOCATABLE      :: altmax                     !! active layer thickness (m)
  REAL(R_STD), DIMENSION(:,:), ALLOCATABLE      :: altmax_g                   !! global active layer thickness (m)
  REAL(r_std), DIMENSION(:,:,:),  ALLOCATABLE  :: carbon_surf_g
  REAL(r_std), DIMENSION(:,:,:),  ALLOCATABLE  :: carbon_surf                 !! vertically-integrated (diagnostic) soil carbon pool: active, slow, or passive, (gC/(m**2 of ground))
  REAL(R_STD), ALLOCATABLE, DIMENSION(:,:)      :: fixed_cryoturbation_depth  !! depth to hold cryoturbation to for fixed runs
  LOGICAL, SAVE                             :: satsoil = .FALSE.
  LOGICAL                                   :: reset_soilc = .false.

  INTEGER(i_std)                            :: start_2d(2), count_2d(2) 
  INTEGER(i_std)                            :: start_4d(4), count_4d(4), start_3d(3), count_3d(3)

  INTEGER(i_std)                            :: nc_opts
!_ =================================================================================================================================
 
  CALL Init_orchidee_para
  CALL init_timer

! Set specific write level to forcesoil using WRITELEVEL_forcesoil=[0-4] in run.def. 
! The global printlev is used as default value. 
  printlev_loc=get_printlev('forcesoil')

!-
! Configure the number of PFTS 
!-
  ok_pc=.FALSE. 
  CALL getin_p('OK_PC',ok_pc)
  ! 1. Read the number of PFTs
  !
  !Config Key   = NVM
  !Config Desc  = number of PFTs  
  !Config If    = OK_SECHIBA or OK_STOMATE
  !Config Def   = 13
  !Config Help  = The number of vegetation types define by the user
  !Config Units = [-]
  CALL getin_p('NVM',nvm)

  ! 2. Allocation
  ALLOCATE(pft_to_mtc(nvm),stat=ier)
  IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'pft_to_mtc : error in memory allocation', '', '')

  ! 3. Initialisation of the correspondance table
  pft_to_mtc(:) = undef_int
  
  ! 4.Reading of the conrrespondance table in the .def file
  !
  !Config Key   = PFT_TO_MTC
  !Config Desc  = correspondance array linking a PFT to MTC
  !Config if    = OK_SECHIBA or OK_STOMATE
  !Config Def   = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13
  !Config Help  =
  !Config Units = [-]
  CALL getin_p('PFT_TO_MTC',pft_to_mtc)

  ! 4.1 if nothing is found, we use the standard configuration
  IF(nvm <= nvmc ) THEN
     IF(pft_to_mtc(1) == undef_int) THEN
        WRITE(numout,*) 'Note to the user : we will use ORCHIDEE to its standard configuration'
        pft_to_mtc(:) = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 /)
     ENDIF
  ELSE   
     IF(pft_to_mtc(1) == undef_int) THEN
        WRITE(numout,*)' The array PFT_TO_MTC is empty : we stop'
     ENDIF
  ENDIF
  
  ! 4.2 What happened if pft_to_mtc(j) > nvmc (if the mtc doesn't exist)?
  DO i = 1, nvm
     IF(pft_to_mtc(i) > nvmc) THEN
        CALL ipslerr_p(3, 'forcesoil', 'the MTC you chose doesnt exist', 'we stop reading pft_to_mtc', '')
     ENDIF
  ENDDO
  
  ! 4.3 Check if pft_to_mtc(1) = 1 
  IF(pft_to_mtc(1) /= 1) THEN
     CALL ipslerr_p(3, 'forcesoil', 'the first pft has to be the bare soil', 'we stop reading next values of pft_to_mtc', '')
  ENDIF

  DO i = 2,nvm
     IF(pft_to_mtc(i) == 1) THEN
        CALL ipslerr_p(3, 'forcesoil', 'only pft_to_mtc(1) has to be the bare soil', 'we stop reading next values of pft_to_mtc', '')
     ENDIF
  ENDDO
  
  ! 5. Allocate and initialize natural and is_c4
  
  ! 5.1 Memory allocation
  ALLOCATE(natural(nvm),stat=ier)
  IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'natural : error in memory allocation', '', '')

  ALLOCATE(is_c4(nvm),stat=ier)
  IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'is_c4 : error in memory allocation', '', '')

  ALLOCATE(permafrost_veg_exists(nvm),stat=ier)
  IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'permafrost_veg_exists : error in memory allocation', '', '')

  ! 5.2 Initialisation
  DO i = 1, nvm
     natural(i) = natural_mtc(pft_to_mtc(i))
     is_c4(i) = is_c4_mtc(pft_to_mtc(i))
  ENDDO

  DO i = 1, nvm
     permafrost_veg_exists(i) = permafrost_veg_exists_mtc(pft_to_mtc(i))
  ENDDO

!!!! yidi  Allocate and initialize is_oilpalm/is_oilpalm_ffbharvest
  ALLOCATE(is_oilpalm(nvm),stat=ier)
  IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'is_oilpalm : error in memory allocation', '', '')

  DO i = 1, nvm
     is_oilpalm(i) = is_oilpalm_mtc(pft_to_mtc(i))
  ENDDO

  ALLOCATE(is_oilpalm_ffbharvest(nvm),stat=ier)
  IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'is_oilpalm_ffbharvest : error in memory allocation', '', '')

  DO i = 1, nvm
     is_oilpalm_ffbharvest(i) = is_oilpalm_ffbharvest_mtc(pft_to_mtc(i))
  ENDDO
!!!! yidi

  !!- 
  !! 1. Initialisation stage
  !! Reading a set of input files, allocating variables and preparing output restart file.     
  !!-
  ! Define restart file name
  ! for reading initial conditions (sto_restname_in var) and for writting final conditions (sto_restname_out var). 
  ! User values are used if present in the .def file.
  ! If not present, default values (stomate_start.nc and stomate_rest_out.c) are used.
  !-
  IF (is_root_prc) THEN
     sto_restname_in = 'stomate_start.nc'
     CALL getin ('STOMATE_RESTART_FILEIN',sto_restname_in)
     WRITE(numout,*) 'STOMATE INPUT RESTART_FILE: ',TRIM(sto_restname_in)
     sto_restname_out = 'stomate_rest_out.nc'
     CALL getin ('STOMATE_RESTART_FILEOUT',sto_restname_out)
     WRITE(numout,*) 'STOMATE OUTPUT RESTART_FILE: ',TRIM(sto_restname_out)
     IF (ok_pc) CALL getin ('satsoil', satsoil)
     !-
     ! Open the input file and Get some Dimension and Attributes ID's 
     !-
     CALL nccheck( NF90_OPEN (sto_restname_in, NF90_NOWRITE, rest_id_sto))
     CALL nccheck( NF90_INQUIRE_DIMENSION (rest_id_sto,1,len=iim_g))
     CALL nccheck( NF90_INQUIRE_DIMENSION (rest_id_sto,2,len=jjm_g))
     CALL nccheck( NF90_INQ_VARID (rest_id_sto, "time", iv))
     CALL nccheck( NF90_GET_ATT (rest_id_sto, iv, 'calendar',thecalendar))
     CALL nccheck( NF90_CLOSE (rest_id_sto))
     i=INDEX(thecalendar,ACHAR(0))
     IF ( i > 0 ) THEN
        thecalendar(i:20)=' '
     ENDIF
     !-
     ! Allocate longitudes and latitudes
     !-
     ALLOCATE (lon(iim_g,jjm_g), stat=ier)
     IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'lon : error in memory allocation', '', '')
     ALLOCATE (lat(iim_g,jjm_g), stat=ier)
     IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'lat : error in memory allocation', '', '')
     lon(:,:) = zero
     lat(:,:) = zero
     lev(1)   = zero
     !-
     CALL restini &
          & (sto_restname_in, iim_g, jjm_g, lon, lat, llm, lev, &
          &  sto_restname_out, itau_dep, date0, dt_files, rest_id_sto, &
          &  use_compression=NC_COMPRESSION_ENABLE )
  ENDIF

  CALL bcast(date0)
  CALL bcast(thecalendar)
  WRITE(numout,*) "calendar = ",thecalendar
  !-
  ! calendar
  !-
  CALL ioconf_calendar (thecalendar)
  CALL ioget_calendar  (one_year,one_day)
  CALL ioconf_startdate(date0)
  !
  IF (ok_pc) THEN
         !- Permafrost variables (zz_deep and zz_coef_deep are constants)
         ALLOCATE (zz_deep(ndeep), stat=ier)
         IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'zz_deep : error in memory allocation', '', '')
         ALLOCATE (zz_coef_deep(ndeep), stat=ier)
         IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'zz_coef_deep : error in memory allocation', '', '')
  ENDIF
  !-
  ! define forcing file's name (Cforcing_name var)
  ! User value is used if present in the .def file
  ! If not, default (NONE) is used
  !-
  IF (.NOT.ok_pc) THEN
      Cforcing_name = 'NONE'
      CALL getin ('STOMATE_CFORCING_NAME',Cforcing_name)
  ELSE
      Cforcing_name = 'stomate_Cforcing_permafrost.nc'
      CALL getin ('STOMATE_CFORCING_PF_NM',Cforcing_name)
  ENDIF
  !
  !! For master process only
  !
  IF (is_root_prc) THEN
     !-
     ! Open FORCESOIL's forcing file to read some basic info (dimensions, variable ID's)
     ! and allocate variables.
     !-
     CALL nccheck( NF90_OPEN (TRIM(Cforcing_name),IOR(NF90_NOWRITE,NF90_NETCDF4),Cforcing_id))
     !-
     ! Total Domain size is stored in nbp_glo variable
     !-
     CALL nccheck( NF90_GET_ATT (Cforcing_id,NF90_GLOBAL,'kjpindex',x_tmp))
     nbp_glo = NINT(x_tmp)
     !-
     ! Number of values stored per year in the forcing file is stored in nparan var.
     !-
     CALL nccheck( NF90_GET_ATT (Cforcing_id,NF90_GLOBAL,'nparan',x_tmp))
     nparan = NINT(x_tmp)
     CALL nccheck( NF90_GET_ATT (Cforcing_id,NF90_GLOBAL,'nbyear',x_tmp))
     nbyear = NINT(x_tmp)
     !-
     ALLOCATE (indices_g(nbp_glo), stat=ier)
     IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'indices_g : error in memory allocation', '', '')
     ALLOCATE (clay_g(nbp_glo), stat=ier)
     IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'clay_g : error in memory allocation', '', '')
     !-
     ALLOCATE (x_indices_g(nbp_glo),stat=ier)
     IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'x_indices_g : error in memory allocation', '', '')
     CALL nccheck( NF90_INQ_VARID (Cforcing_id,'index',v_id))
     CALL nccheck( NF90_GET_VAR   (Cforcing_id,v_id,x_indices_g))
     indices_g(:) = NINT(x_indices_g(:))
     WRITE(numout,*) mpi_rank,"indices globaux : ",indices_g
     DEALLOCATE (x_indices_g)
     !-
     CALL nccheck( NF90_INQ_VARID (Cforcing_id,'clay',v_id))
     CALL nccheck( NF90_GET_VAR   (Cforcing_id,v_id,clay_g))
     !-
     IF (ok_pc) THEN
         !- Permafrost variables (zz_deep and zz_coef_deep are constants)
!         ALLOCATE (zz_deep(ndeep))
!         ALLOCATE (zz_coef_deep(ndeep))
         ALLOCATE (z_organic_g(nbp_glo), stat=ier)
         IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'z_organic_g : error in memory allocation', '', '')
         CALL nccheck( NF90_INQ_VARID (Cforcing_id,'zz_deep',v_id))
         CALL nccheck( NF90_GET_VAR (Cforcing_id,v_id,zz_deep))
         CALL nccheck( NF90_INQ_VARID (Cforcing_id,'zz_coef_deep',v_id))
         CALL nccheck( NF90_GET_VAR (Cforcing_id,v_id,zz_coef_deep))
         CALL nccheck( NF90_INQ_VARID (Cforcing_id,'z_organic',v_id))
         CALL nccheck( NF90_GET_VAR (Cforcing_id,v_id,z_organic_g))
     ENDIF
     CALL nccheck( NF90_CLOSE(Cforcing_id) )
     ! time step of forcesoil program (in days)
     !-
     dt_forcesoil = one_year / FLOAT(nparan)
     WRITE(numout,*) 'time step (d): ',dt_forcesoil
     WRITE(numout,*) 'nparan: ',nparan
     WRITE(numout,*) 'nbyear: ',nbyear    
     !-
     ! read and write the variables in the output restart file we do not modify within the Forcesoil program
     ! ie all variables stored in the input restart file except those stored in taboo_vars
     !-
     IF (.NOT. ok_pc) THEN
         !-
         taboo_vars ='$lon$ $lat$ $lev$ $nav_lon$ $nav_lat$ $nav_lev$ $time$ $time_steps$ '// &
              &             '$day_counter$ $dt_days$ $date$ $carbon$ '
         !-
     ELSE
         !-
         taboo_vars = '$nav_lon$ $nav_lat$ $nav_lev$ $time$ $time_steps$ '// &
         &            '$day_counter$ $dt_days$ $date$ $deepC_a$ $deepC_s$ '// &
         &            '$deepC_p$ $O2_soil$ $CH4_soil$ $O2_snow$ $CH4_snow$ '// &
         &            '$altmax$ '  
         !-
     ENDIF ! ok_pc
     !-
     CALL ioget_vname(rest_id_sto, nbvar, varnames)
     !-
     ! read and write some special variables (1D or variables that we need)
     !-
     var_name = 'day_counter'
     CALL restget (rest_id_sto, var_name, 1, 1, 1, itau_dep, .TRUE., xtmp)
     CALL restput (rest_id_sto, var_name, 1, 1, 1, itau_dep, xtmp)
     !-
     var_name = 'dt_days'
     CALL restget (rest_id_sto, var_name, 1, 1, 1, itau_dep, .TRUE., xtmp)
     CALL restput (rest_id_sto, var_name, 1, 1, 1, itau_dep, xtmp)
     !-
     var_name = 'date'
     CALL restget (rest_id_sto, var_name, 1, 1, 1, itau_dep, .TRUE., xtmp)
     CALL restput (rest_id_sto, var_name, 1, 1, 1, itau_dep, xtmp)
     !-
     DO iv=1,nbvar
        !-- check if the variable is to be written here
        IF (INDEX(taboo_vars,'$'//TRIM(varnames(iv))//'$') == 0 ) THEN
           !---- get variable dimensions, especially 3rd dimension
           CALL ioget_vdim &
                &      (rest_id_sto, varnames(iv), varnbdim_max, varnbdim, vardims)
           l1d = ALL(vardims(1:varnbdim) == 1)
           WRITE(*,*) TRIM(varnames(iv)),": ", varnbdim, "-", vardims(1:varnbdim)," l1d=" ,l1d
           !---- read it
           IF (l1d) THEN
              CALL restget &
                   &        (rest_id_sto, TRIM(varnames(iv)), 1, 1, &
                   &         1, itau_dep, .TRUE., xtmp)
              CALL restput &
                   &        (rest_id_sto, TRIM(varnames(iv)), 1, 1, &
                   &         1, itau_dep, xtmp)
           ELSE
              orch_vardims = varnbdim - 1 ! exclude time dimension introduced by IOIPSL
              ! Deal for different number of dimension
              IF (orch_vardims == 2) THEN
                 !----
                 ! vardims X(1), Y(2), remaining variables, time(last position)
                 ALLOCATE( var_2d(nbp_glo), stat=ier)
                 IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'var_4d : error in memory allocation', '', '')
                 !----
                 CALL restget &
                      &        (rest_id_sto, TRIM(varnames(iv)), nbp_glo, 1, &
                      &         1, itau_dep, .TRUE., var_2d, "gather", nbp_glo, indices_g)
                 CALL restput &
                      &        (rest_id_sto, TRIM(varnames(iv)), nbp_glo, 1, &
                      &         1, itau_dep, var_2d, 'scatter',  nbp_glo, indices_g)
                 !----
                 DEALLOCATE(var_2d)
              ELSE IF (orch_vardims == 3) THEN
                 !----
                 ! vardims X(1), Y(2), remaining variables, time(last position)
                 ALLOCATE( var_3d(nbp_glo,vardims(3)), stat=ier)
                 IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'var_3d : error in memory allocation', '', '')
                 !----
                 CALL restget &
                      &        (rest_id_sto, TRIM(varnames(iv)), nbp_glo, vardims(3), &
                      &         1, itau_dep, .TRUE., var_3d, "gather", nbp_glo, indices_g)
                 CALL restput &
                      &        (rest_id_sto, TRIM(varnames(iv)), nbp_glo, vardims(3), &
                      &         1, itau_dep, var_3d, 'scatter',  nbp_glo, indices_g)
                 !----
                 DEALLOCATE(var_3d)
              ELSE IF (orch_vardims == 4) THEN
                 !----
                 ! vardims X(1), Y(2), remaining variables, time(last position)
                 ALLOCATE( var_4d(nbp_glo,vardims(3),vardims(4)), stat=ier)
                 IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'var_4d : error in memory allocation', '', '')
                 !----
                 CALL restget &
                      &        (rest_id_sto, TRIM(varnames(iv)), nbp_glo, vardims(3), &
                      &         vardims(4), itau_dep, .TRUE., var_4d, "gather", nbp_glo, indices_g)
                 CALL restput &
                      &        (rest_id_sto, TRIM(varnames(iv)), nbp_glo, vardims(3), &
                      &         vardims(4), itau_dep, var_4d, 'scatter',  nbp_glo, indices_g)
                 !----
                 DEALLOCATE(var_4d)
              ELSE IF (orch_vardims == 5) THEN
                 !----
                 ! vardims X(1), Y(2), remaining variables, time(last position)
                 ALLOCATE( var_5d(nbp_glo,vardims(3),vardims(4),vardims(5)), stat=ier)
                 IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'var_5d : error in memory allocation', '', '')
                 !----
                 CALL restget &
                      &        (rest_id_sto, TRIM(varnames(iv)), nbp_glo, vardims(3), &
                      &         vardims(4), vardims(5), itau_dep, .TRUE., var_5d, "gather", nbp_glo, indices_g)
                 CALL restput &
                      &        (rest_id_sto, TRIM(varnames(iv)), nbp_glo, vardims(3), &
                      &         vardims(4), vardims(5), itau_dep, var_5d, 'scatter',  nbp_glo, indices_g)
                 !----
                 DEALLOCATE(var_5d)
              ELSE IF (orch_vardims == 6) THEN
                 !----
                 ! vardims X(1), Y(2), remaining variables, time(last position)
                 ALLOCATE( var_6d(nbp_glo,vardims(3),vardims(4),vardims(5),vardims(6)), stat=ier)
                 IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'var_6d : error in memory allocation', '', '')
                 !----
                 CALL restget &
                      &        (rest_id_sto, TRIM(varnames(iv)), nbp_glo, vardims(3), &
                      &         vardims(4), vardims(5), vardims(6), itau_dep, .TRUE., var_6d, "gather", nbp_glo, indices_g)
                 CALL restput &
                      &        (rest_id_sto, TRIM(varnames(iv)), nbp_glo, vardims(3), &
                      &         vardims(4), vardims(5), vardims(6), itau_dep, var_6d, 'scatter',  nbp_glo, indices_g)
                 !----
                 DEALLOCATE(var_6d)
              ELSE
                 WRITE( msg3, '(i5)' ) orch_vardims
                 CALL ipslerr(3, 'forcesoil', 'Varialbes(1) Restart Read/Write not implement for N dimensions(2)', &
                            & TRIM(varnames(iv)), TRIM(msg3))
              ENDIF
           ENDIF
        ENDIF
     ENDDO
     ! Length of the run (in Years)
     ! User value is used if present in the .def file
     ! If not, default value (10000 Years) is used
     !-
     WRITE(time_str,'(a)') '10000Y'
     CALL getin('TIME_LENGTH', time_str)
     write(numout,*) 'Number of years for carbon spinup : ',time_str
     ! transform into itau
     CALL tlen2itau(time_str, dt_forcesoil*one_day, date0, itau_len)
     write(numout,*) 'Number of time steps to do: ',itau_len

     ! read soil carbon stocks values stored in the input restart file
     !-
     IF (.NOT. ok_pc) THEN
           ALLOCATE(carbon_g(nbp_glo,ncarb,nvm), stat=ier)
           IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'carbon_g : error in memory allocation', '', '')
           carbon_g(:,:,:) = val_exp
           CALL restget &
                &    (rest_id_sto, 'carbon', nbp_glo, ncarb , nvm, itau_dep, &
                &     .TRUE., carbon_g, 'gather', nbp_glo, indices_g)
           IF (ALL(carbon_g == val_exp)) carbon_g = zero
           WRITE(numout,*) "date0 : ",date0, itau_dep
     ELSE
           !-
           ! Permafrost carbon
           !-
           ALLOCATE(carbon_g(nbp_glo,ncarb,nvm), stat=ier)
           IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'carbon_g : error in memory allocation', '', '')
           carbon_g(:,:,:) = 0.
           ALLOCATE(carbon_surf_g(nbp_glo,ncarb,nvm), stat=ier)
           IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'carbon_surf_g : error in memory allocation', '', '')
           carbon_surf_g(:,:,:) = 0.
           ALLOCATE(deepC_a_g(nbp_glo,ndeep,nvm), stat=ier)
           IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'deepC_a_g : error in memory allocation', '', '')
           ALLOCATE(deepC_s_g(nbp_glo,ndeep,nvm), stat=ier)
           IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'deepC_s_g : error in memory allocation', '', '')
           ALLOCATE(deepC_p_g(nbp_glo,ndeep,nvm), stat=ier)
           IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'deepC_p_g : error in memory allocation', '', '')
           ALLOCATE(O2_soil_g(nbp_glo,ndeep,nvm), stat=ier)
           IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'O2_soil_g : error in memory allocation', '', '')
           ALLOCATE(CH4_soil_g(nbp_glo,ndeep,nvm), stat=ier)
           IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'CH4_soil_g : error in memory allocation', '', '')
           ALLOCATE(O2_snow_g(nbp_glo,nsnow,nvm), stat=ier)
           IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'O2_snow_g : error in memory allocation', '', '')
           ALLOCATE(CH4_snow_g(nbp_glo,nsnow,nvm), stat=ier)
           IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'CH4_snow_g : error in memory allocation', '', '')
           ALLOCATE(altmax_g(nbp_glo,nvm), stat=ier)
           IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'altmax_g : error in memory allocation', '', '')

           deepC_a_g(:,:,:) = val_exp
           CALL restget (rest_id_sto, 'deepC_a', nbp_glo, ndeep, nvm, itau_dep, &
                &               .TRUE., deepC_a_g, 'gather', nbp_glo, indices_g)
           IF (ALL(deepC_a_g == val_exp)) deepC_a_g = zero
    
           deepC_s_g(:,:,:) = val_exp
           CALL restget (rest_id_sto, 'deepC_s', nbp_glo, ndeep, nvm, itau_dep, &
                &               .TRUE., deepC_s_g, 'gather', nbp_glo, indices_g)
           IF (ALL(deepC_s_g == val_exp)) deepC_s_g = zero
         
           deepC_p_g(:,:,:) = val_exp
           CALL restget (rest_id_sto, 'deepC_p', nbp_glo, ndeep, nvm, itau_dep, &
                &               .TRUE., deepC_p_g, 'gather', nbp_glo, indices_g)
           IF (ALL(deepC_p_g == val_exp)) deepC_p_g = zero

           var_name= 'altmax'
           altmax_g(:,:) = val_exp
           CALL restget (rest_id_sto, var_name, nbp_glo, nvm, 1, itau_dep, .TRUE., altmax_g, "gather", nbp_glo, indices_g)
           IF ( ALL( altmax_g(:,:) .EQ. val_exp ) ) THEN
               CALL ipslerr(3, 'forcesoil', 'altmax is not found in stomate restart file', '', '')
           END IF

           CALL getin('reset_soilc', reset_soilc)
           IF (reset_soilc) THEN
              CALL ipslerr(1, 'forcesoil', 'deepC_a, deepC_s and deeC_p',  & 
                            'are ignored and set to zero value due to', 'reset_soilc option')
              deepC_a_g(:,:,:) = zero
              deepC_s_g(:,:,:) = zero
              deepC_p_g(:,:,:) = zero
           ENDIF
         
           O2_soil_g(:,:,:) = val_exp
           CALL restget (rest_id_sto, 'O2_soil', nbp_glo, ndeep, nvm, itau_dep, &
                &               .TRUE., O2_soil_g, 'gather', nbp_glo, indices_g)
           IF (ALL(O2_soil_g == val_exp)) O2_soil_g = O2_init_conc
    
           CH4_soil_g(:,:,:) = val_exp
           CALL restget (rest_id_sto, 'CH4_soil', nbp_glo, ndeep, nvm, itau_dep, &
               &               .TRUE., CH4_soil_g, 'gather', nbp_glo, indices_g)
           IF (ALL(CH4_soil_g == val_exp)) CH4_soil_g =  CH4_init_conc
    
           O2_snow_g(:,:,:) = val_exp
           CALL restget (rest_id_sto, 'O2_snow', nbp_glo, nsnow, nvm, itau_dep, &
                &               .TRUE., O2_snow_g, 'gather', nbp_glo, indices_g)
           IF (ALL(O2_snow_g == val_exp)) O2_snow_g =  O2_init_conc
    
           CH4_snow_g(:,:,:) = val_exp
           CALL restget (rest_id_sto, 'CH4_snow', nbp_glo, nsnow, nvm, itau_dep, &
               &               .TRUE., CH4_snow_g, 'gather', nbp_glo, indices_g)
           IF (ALL(CH4_snow_g == val_exp)) CH4_snow_g = CH4_init_conc
     ENDIF ! ok_pc
  ENDIF ! is_root_prc
  !
  CALL bcast(nbp_glo)
  CALL bcast(iim_g)
  CALL bcast(jjm_g)
  IF (.NOT. ALLOCATED(indices_g)) ALLOCATE (indices_g(nbp_glo))
  CALL bcast(indices_g)
  CALL bcast(nparan)
  CALL bcast(nbyear)
  CALL bcast(dt_forcesoil)
  CALL bcast(itau_dep)
  CALL bcast(itau_len)
  IF (ok_pc) THEN
      CALL bcast(zz_deep)
      CALL bcast(zz_coef_deep) 
  ENDIF
  !
  ! we must initialize data_para :
  CALL init_orchidee_data_para_driver(nbp_glo,indices_g)

  kjpindex=nbp_loc
  jjm=jj_nb
  iim=iim_g
  IF (printlev_loc>=3) WRITE(numout,*) "Local grid : ",kjpindex,iim,jjm
  !-
  ! Analytical spinup is set to false
  !
  spinup_analytic = .FALSE.
  !-
  ! read soil carbon inputs, water and temperature stresses on OM
  ! decomposition 
  ! into the forcing file - We read an average year.
  !-
  IF (.NOT. ok_pc) THEN
     !-
     ! Open FORCESOIL's forcing file to read some basic info (dimensions, variable ID's)
     ! and allocate variables.
     !-
#ifdef CPP_PARA
     nc_opts = IOR(NF90_NOWRITE, NF90_MPIIO)
     CALL nccheck( NF90_OPEN (TRIM(Cforcing_name),IOR(nc_opts, NF90_NETCDF4),Cforcing_id, &
                 & comm = MPI_COMM_ORCH, info = MPI_INFO_NULL ))
#else
     CALL nccheck( NF90_OPEN (TRIM(Cforcing_name), IOR(NF90_NOWRITE, NF90_NETCDF4),Cforcing_id) )
#endif
     ALLOCATE(control_temp(kjpindex,nlevs,nparan*nbyear), stat=ier)
     IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'control_temp : error in memory allocation', '', '')
     ALLOCATE(control_moist(kjpindex,nlevs,nparan*nbyear), stat=ier)
     IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'control_moist : error in memory allocation', '', '')
     ALLOCATE(soilcarbon_input(kjpindex,ncarb,nvm,nparan*nbyear), stat=ier)
     IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'soilcarbon_input : error in memory allocation', '', '')
     ALLOCATE(veget_max(kjpindex,nvm,nparan*nbyear), stat=ier)
     IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'veget_max : error in memory allocation', 'No ok_pc', '')
     !
     CALL nccheck( NF90_GET_ATT (Cforcing_id,NF90_GLOBAL,'dt_sechiba',dt_sechiba))
     !-
     start_4d = (/ nbp_mpi_para_begin(mpi_rank), 1, 1, 1 /)
     count_4d = (/ nbp_mpi_para(mpi_rank), ncarb, nvm, nparan*nbyear /)
     CALL nccheck( NF90_INQ_VARID (Cforcing_id,'soilcarbon_input',v_id))
     CALL nccheck( NF90_GET_VAR   (Cforcing_id,v_id,soilcarbon_input,  &
                    &  start = start_4d, count = count_4d ))
     !
     start_3d = (/ nbp_mpi_para_begin(mpi_rank), 1, 1 /)
     count_3d = (/ nbp_mpi_para(mpi_rank), nlevs, nparan*nbyear /)
     CALL nccheck( NF90_INQ_VARID (Cforcing_id,   'control_moist',v_id))
     CALL nccheck( NF90_GET_VAR   (Cforcing_id,v_id,control_moist, &
                       & start = start_3d, count = count_3d))
     CALL nccheck( NF90_INQ_VARID (Cforcing_id,    'control_temp',v_id))
     CALL nccheck( NF90_GET_VAR   (Cforcing_id,v_id,control_temp, &
                       & start = start_3d, count = count_3d))
     count_3d = (/ nbp_mpi_para(mpi_rank), nvm, nparan*nbyear /)
     CALL nccheck( NF90_INQ_VARID (Cforcing_id,    'veget_max',v_id))
     CALL nccheck( NF90_GET_VAR   (Cforcing_id,v_id,veget_max, &
                       & start = start_3d, count = count_3d))
     !- Close Netcdf carbon permafrost file reference
     CALL nccheck( NF90_CLOSE (Cforcing_id))
  ELSE
     !-
     ! Read permafrost-related soil carbon from Cforcing_id
     !-
     ALLOCATE(pb(kjpindex,nparan*nbyear), stat=ier)
     IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'pb : error in memory allocation', '', '')
     ALLOCATE(snow(kjpindex,nparan*nbyear), stat=ier)
     IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'snow : error in memory allocation', '', '')
     ALLOCATE(tprof(kjpindex,ndeep,nvm,nparan*nbyear), stat=ier)
     IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'tprof : error in memory allocation', '', '')
     ALLOCATE(fbact(kjpindex,ndeep,nvm,nparan*nbyear), stat=ier)
     IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'fbact : error in memory allocation', '', '')
     ALLOCATE(hslong(kjpindex,ndeep,nvm,nparan*nbyear), stat=ier)
     IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'hslong : error in memory allocation', '', '')
     ALLOCATE(veget_max(kjpindex,nvm,nparan*nbyear), stat=ier)
     IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'veget_max : error in memory allocation', '', '')
     ALLOCATE(rprof(kjpindex,nvm,nparan*nbyear), stat=ier)
     IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'rprof : error in memory allocation', '', '')
     ALLOCATE(tsurf(kjpindex,nparan*nbyear), stat=ier)
     IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'tsurf : error in memory allocation', '', '')
     ALLOCATE(lalo(kjpindex,2), stat=ier)
     IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'lalo : error in memory allocation', '', '')
     ALLOCATE(snowdz(kjpindex,nsnow,nparan*nbyear), stat=ier)
     IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'snowdz_ : error in memory allocation', '', '')
     ALLOCATE(snowrho(kjpindex,nsnow,nparan*nbyear), stat=ier)
     IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'snowrho : error in memory allocation', '', '')
     ALLOCATE(soilcarbon_input(kjpindex,ncarb,nvm,nparan*nbyear), stat=ier)
     IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'soilcarbon_input : error in memory allocation', '', '')
     !-
     CALL stomate_io_carbon_permafrost_read(Cforcing_name,  nparan,      nbyear,&
                nbp_mpi_para_begin(mpi_rank),   nbp_mpi_para(mpi_rank),         &
                soilcarbon_input,               pb,         snow,       tsurf,  &
                tprof,                          fbact,      hslong,     rprof,  &
                lalo,                           snowdz,     snowrho,    veget_max )
  ENDIF ! ok_pc
  !---
  !--- Create the index table
  !---
  !--- This job returns a LOCAL kindex.
  !---
  ALLOCATE (indices(kjpindex),stat=ier)
  IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'indices : error in memory allocation', '', '')
  !
  !! scattering to all processes in parallel mode
  !
  CALL scatter(indices_g,indices)
  indices(1:kjpindex)=indices(1:kjpindex)-(jj_begin-1)*iim_g
  IF (printlev_loc>=3) WRITE(numout,*) mpi_rank,"indices locaux = ",indices(1:kjpindex)
  !
  ! Initialize _index variables
  CALL stomate_init_index(nbp_glo, kjpindex, indices)
  !-
  ! Allocation of the variables for a processor
  !-
  ALLOCATE(clay(kjpindex), stat=ier)
  IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'clay : error in memory allocation', '', '')
  ALLOCATE(carbon(kjpindex,ncarb,nvm), stat=ier)
  IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'indices : error in memory allocation', '', '')
  !-
  IF (.NOT. ok_pc) THEN
     ALLOCATE(resp_hetero_soil(kjpindex,nvm), stat=ier)
     IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'resp_hetero_soil : error in memory allocation', '', '')
     ALLOCATE(matrixA(kjpindex,nvm,nbpools,nbpools), stat=ier)
     IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'matrixA : error in memory allocation', '', '')
     DO i = 1,nbpools
        matrixA(:,:,i,i) = un
     ENDDO
  ELSE
     ALLOCATE(carbon_surf(kjpindex,ncarb,nvm), stat=ier)
     IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'carbon_surf : error in memory allocation', '', '')
     ALLOCATE(deepC_a(kjpindex,ndeep,nvm), stat=ier)
     IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'deepC_a : error in memory allocation', '', '')
     ALLOCATE(deepC_s(kjpindex,ndeep,nvm), stat=ier)
     IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'deepC_s : error in memory allocation', '', '')
     ALLOCATE(deepC_p(kjpindex,ndeep,nvm), stat=ier)
     IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'deepC_p : error in memory allocation', '', '')
     ALLOCATE(O2_soil(kjpindex,ndeep,nvm), stat=ier)
     IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'O2_soil : error in memory allocation', '', '')
     ALLOCATE(CH4_soil(kjpindex,ndeep,nvm), stat=ier)
     IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'CH4_soil : error in memory allocation', '', '')
     ALLOCATE(O2_snow(kjpindex,nsnow,nvm), stat=ier)
     IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'O2_snow : error in memory allocation', '', '')
     ALLOCATE(CH4_snow(kjpindex,nsnow,nvm), stat=ier)
     IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'CH4_snow : error in memory allocation', '', '')
     ALLOCATE(altmax(kjpindex,nvm), stat=ier)
     IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'altmax : error in memory allocation', '', '')
     ALLOCATE(z_organic(kjpindex), stat=ier)
     IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'z_organic : error in memory allocation', '', '')
  ENDIF
  !-
  ! Initialization of the variables for a processor
  !-
  CALL Scatter(clay_g,clay)
  CALL Scatter(carbon_g,carbon)
  IF (ok_pc) THEN
     CALL Scatter(z_organic_g,z_organic)
     CALL Scatter(carbon_surf_g,carbon_surf)
     CALL Scatter(deepC_a_g,deepC_a)
     CALL Scatter(deepC_s_g,deepC_s)
     CALL Scatter(deepC_p_g,deepC_p)
     CALL Scatter(O2_soil_g,O2_soil)
     CALL Scatter(CH4_soil_g,CH4_soil)
     CALL Scatter(O2_snow_g,O2_snow)
     CALL Scatter(CH4_snow_g,CH4_snow)
     CALL Scatter(altmax_g,altmax)
   ENDIF
!-
! Configuration of the parameters
!-
  !Config Key   = FRAC_CARB_AP
  !Config Desc  = frac carb coefficients from active pool: depends on clay
  !content
  !Config if    = OK_STOMATE 
  !Config Def   = 0.004
  !Config Help  = fraction of the active pool going to the passive pool
  !Config Units = [-]
  CALL getin_p('FRAC_CARB_AP',frac_carb_ap)
  !
  !Config Key   = FRAC_CARB_SA
  !Config Desc  = frac_carb_coefficients from slow pool
  !Config if    = OK_STOMATE 
  !Config Def   = 0.42
  !Config Help  = fraction of the slow pool going to the active pool
  !Config Units = [-] 
  CALL getin_p('FRAC_CARB_SA',frac_carb_sa)
  !
  !Config Key   = FRAC_CARB_SP
  !Config Desc  = frac_carb_coefficients from slow pool
  !Config if    = OK_STOMATE 
  !Config Def   = 0.03
  !Config Help  = fraction of the slow pool going to the passive pool
  !Config Units = [-] 
  CALL getin_p('FRAC_CARB_SP',frac_carb_sp)
  !
  !Config Key   = FRAC_CARB_PA
  !Config Desc  = frac_carb_coefficients from passive pool
  !Config if    = OK_STOMATE 
  !Config Def   = 0.45
  !Config Help  = fraction of the passive pool going to the passive pool
  !Config Units = [-] 
  CALL getin_p('FRAC_CARB_PA',frac_carb_pa)
  !
  !Config Key   = FRAC_CARB_PS
  !Config Desc  = frac_carb_coefficients from passive pool
  !Config if    = OK_STOMATE 
  !Config Def   = 0.0
  !Config Help  = fraction of the passive pool going to the passive pool
  !Config Units = [-]
  CALL getin_p('FRAC_CARB_PS',frac_carb_ps)
  !
  !Config Key   = ACTIVE_TO_PASS_CLAY_FRAC
  !Config Desc  = 
  !Config if    = OK_STOMATE 
  !Config Def   =  .68  
  !Config Help  =
  !Config Units = [-]
  CALL getin_p('ACTIVE_TO_PASS_CLAY_FRAC',active_to_pass_clay_frac)
  !
  !Config Key   = CARBON_TAU_IACTIVE
  !Config Desc  = residence times in carbon pools
  !Config if    = OK_STOMATE 
  !Config Def   = 0.149
  !Config Help  =
  !Config Units = [days] 
  CALL getin_p('CARBON_TAU_IACTIVE',carbon_tau_iactive)
  !
  !Config Key   = CARBON_TAU_ISLOW
  !Config Desc  = residence times in carbon pools
  !Config if    = OK_STOMATE 
  !Config Def   = 5.48
  !Config Help  =
  !Config Units = [days]
  CALL getin_p('CARBON_TAU_ISLOW',carbon_tau_islow)
  !
  !Config Key   = CARBON_TAU_IPASSIVE
  !Config Desc  = residence times in carbon pools
  !Config if    = OK_STOMATE 
  !Config Def   = 241.
  !Config Help  =
  !Config Units = [days]
  CALL getin_p('CARBON_TAU_IPASSIVE',carbon_tau_ipassive)
  !
  !Config Key   = FLUX_TOT_COEFF
  !Config Desc  =
  !Config if    = OK_STOMATE 
  !Config Def   = 1.2, 1.4,.75
  !Config Help  =
  !Config Units = [days]
  CALL getin_p('FLUX_TOT_COEFF',flux_tot_coeff)

  !
  !! 2. Computational step
  !! Loop over time - Call of soilcarbon routine at each time step 
  !! Updated soil carbon stocks are stored into carbon variable
  !! We only keep the last value of carbon variable (no time dimension).
  !!-
  IF (.NOT.ok_pc) THEN
      iyear=1
      iatt = 0
      DO i=1,itau_len
         iatt = iatt+1
         IF (iatt > nparan*nbyear) THEN
            IF (printlev>=3) WRITE(numout,*) iyear
            iatt = 1
            iyear=iyear+1
         ENDIF
         CALL soilcarbon &
              &    (kjpindex, dt_forcesoil, clay, &
              &     soilcarbon_input(:,:,:,iatt), &
              &     control_temp(:,:,iatt), control_moist(:,:,iatt), veget_max(:,:,iatt), &
              &     carbon, resp_hetero_soil, &
              &     matrixA)
      ENDDO
      WRITE(numout,*) "End of soilcarbon LOOP."
  ELSE
      IF ( satsoil )  hslong(:,:,:,:) = 1.
      !these variables are only ouputs from deep_carbcycle (thus not necessary for
      !Gather and Scatter)
      ALLOCATE(heat_Zimov(kjpindex,ndeep,nvm), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'heat_Zimov : error in memory allocation', '', '')
      ALLOCATE(sfluxCH4_deep(kjpindex), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'sfluxCH4_deep : error in memory allocation', '', '')
      ALLOCATE(sfluxCO2_deep(kjpindex), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'sfluxCO2_deep : error in memory allocation', '', '')
      ALLOCATE(fixed_cryoturbation_depth(kjpindex,nvm), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'fixed_cryoturbation_depth : error in memory allocation', '', '')
      ALLOCATE(resp_hetero_soil(kjpindex, nvm), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'resp_hetero_soil : error in memory allocation', '', '')
    
      iatt = 0
      iyear=1
      DO i=1,itau_len
        iatt = iatt+1
        IF (iatt > nparan*nbyear) THEN
            IF (printlev>=3) WRITE(numout,*) iyear
            iatt = 1
            iyear=iyear+1
        ENDIF
        WRITE(numout, *) "Forcesoil:: deep_carbcycle, iyear=", iyear
        CALL deep_carbcycle(kjpindex, indices, iatt, dt_forcesoil*one_day, lalo, clay, &
             tsurf(:,iatt), tprof(:,:,:,iatt), hslong(:,:,:,iatt), snow(:,iatt), heat_Zimov, pb(:,iatt), &
             sfluxCH4_deep, sfluxCO2_deep,  &
             deepC_a, deepC_s, deepC_p, O2_soil, CH4_soil, O2_snow, CH4_snow, &
             zz_deep, zz_coef_deep, z_organic, soilcarbon_input(:,:,:,iatt), &
             veget_max(:,:,iatt), rprof(:,:,iatt), altmax,  carbon, carbon_surf, resp_hetero_soil, &
             fbact(:,:,:,iatt), fixed_cryoturbation_depth, snowdz(:,:,iatt), snowrho(:,:,iatt))
      ENDDO
    
  ENDIF
  !!-
  !! 3. write new carbon stocks into the ouput restart file
  !!-
  CALL restput_p (rest_id_sto, 'carbon', nbp_glo, ncarb , nvm, itau_dep, &
         &     carbon, 'scatter', nbp_glo, indices_g)

  IF (ok_pc) THEN
     CALL restput_p (rest_id_sto, 'deepC_a', nbp_glo, ndeep, nvm, itau_dep, &
            &               deepC_a, 'scatter', nbp_glo, indices_g)
     CALL restput_p (rest_id_sto, 'deepC_s', nbp_glo, ndeep, nvm, itau_dep, &
            &               deepC_s, 'scatter', nbp_glo, indices_g)
     CALL restput_p (rest_id_sto, 'deepC_p', nbp_glo, ndeep, nvm, itau_dep, &
            &               deepC_p, 'scatter', nbp_glo, indices_g)
     CALL restput_p (rest_id_sto, 'O2_soil', nbp_glo, ndeep, nvm, itau_dep, &
            &               O2_soil, 'scatter', nbp_glo, indices_g)
     CALL restput_p (rest_id_sto, 'CH4_soil', nbp_glo, ndeep, nvm, itau_dep, &
            &               CH4_soil, 'scatter', nbp_glo, indices_g)
     CALL restput_p (rest_id_sto, 'O2_snow', nbp_glo, nsnow, nvm, itau_dep, &
            &               O2_snow, 'scatter', nbp_glo, indices_g)
     CALL restput_p (rest_id_sto, 'CH4_snow', nbp_glo, nsnow, nvm, itau_dep, &
            &               CH4_snow, 'scatter', nbp_glo, indices_g)
     CALL restput_p (rest_id_sto, 'altmax', nbp_glo, nvm, 1, itau_dep,     &
            &               altmax, 'scatter',  nbp_glo, indices_g)
  ENDIF
  !-
  IF (is_root_prc) THEN
        !- Close restart files
        CALL getin_dump
        CALL restclo
  ENDIF
  !-
#ifdef CPP_PARA
  CALL MPI_FINALIZE(ier)
#endif
  WRITE(numout,*) "End of forcesoil."
  !--------------------
END PROGRAM forcesoil
