PROGRAM driver2oasis
  !---------------------------------------------------------------------
  !-
  !- Reads the forcing file in the ALMA format and feeds the data to the
  !- OASIS coupler.
  !- 
  !---------------------------------------------------------------------
  USE defprec
  !
  USE netcdf
  !
  USE constantes
  !
  USE ioipsl_para
  USE mod_orchidee_para
  !
  USE grid
  USE globgrd
  USE forcing_tools
  !
  USE mod_oasis
  USE timer
  USE mpi
  !-
  IMPLICIT NONE
  !-
  CHARACTER(LEN=80) :: gridfilename
  CHARACTER(LEN=80), DIMENSION(100) :: forfilename
  CHARACTER(LEN=6)   :: comp_name = 'driver'
  !
  !
  CHARACTER(LEN=8)  :: model_guess
  INTEGER(i_std)    :: iim_glo, jjm_glo, file_id
  !- 
  REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: lon_glo, lat_glo, area_glo
  REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: mask_glo
  INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:) :: maskinv_glo
  REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:,:) :: corners_glo
  REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:) :: corners_lon, corners_lat
  INTEGER(i_std) :: nbindex_g, kjpindex
  INTEGER(i_std), ALLOCATABLE, DIMENSION(:) :: kindex_g
  REAL(r_std), DIMENSION(2) :: zoom_lon, zoom_lat
  CHARACTER(LEN=20) :: calendar
  !-
  !- Variables local to each processors.
  !-
  INTEGER(i_std) :: i, j, ik, nbdt, first_point
  INTEGER(i_std) :: nb_forcefile
  INTEGER(i_std) :: itau, itau_offset, itau_sechiba
  REAL(r_std)    :: date_start, dt
  REAL(r_std)    :: timestep_interval(2), timestep_int_next(2), julian
  INTEGER(i_std) :: rest_id, rest_id_stom
  INTEGER(i_std) :: hist_id, hist2_id, hist_id_stom, hist_id_stom_IPCC
!-
!- input fields
!-
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: u             !! Lowest level wind speed
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: v             !! Lowest level wind speed 
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: pb            !! Surface pressure
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: zlev_tq       !! Height of layer for T and q
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: zlev_uv       !! Height of layer for u and v
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: temp_air      !! Air temperature in Kelvin
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: qair          !! Lowest level specific humidity
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: ccanopy       !! CO2 concentration in the canopy
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: cdrag         !! Cdrag
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: precip_rain   !! Rain precipitation
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: precip_snow   !! Snow precipitation
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: lwdown        !! Down-welling long-wave flux 
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: swdown        !! Downwelling surface short-wave flux
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: swnet         !! Net surface short-wave flux
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: solarang      !! Cosine of solar zenith angle
  !-
  !- output fields
  !-
  REAL(r_std), ALLOCATABLE, DIMENSION (:)            :: z0            !! Surface roughness
  REAL(r_std), ALLOCATABLE, DIMENSION (:)            :: coastalflow   !! Diffuse flow of water into the ocean (m^3/dt)
  REAL(r_std), ALLOCATABLE, DIMENSION (:)            :: riverflow     !! Largest rivers flowing into the ocean (m^3/dt)
  REAL(r_std), ALLOCATABLE, DIMENSION (:)            :: tsol_rad      !! Radiative surface temperature
  REAL(r_std), ALLOCATABLE, DIMENSION (:)            :: vevapp        !! Total of evaporation
  REAL(r_std), ALLOCATABLE, DIMENSION (:)            :: temp_sol_new  !! New soil temperature
  REAL(r_std), ALLOCATABLE, DIMENSION (:)            :: qsurf         !! Surface specific humidity
  REAL(r_std), ALLOCATABLE, DIMENSION (:,:)          :: albedo        !! Albedo
  REAL(r_std), ALLOCATABLE, DIMENSION (:)            :: fluxsens      !! Sensible chaleur flux
  REAL(r_std), ALLOCATABLE, DIMENSION (:)            :: fluxlat       !! Latent chaleur flux
  REAL(r_std), ALLOCATABLE, DIMENSION (:)            :: emis          !! Emissivity
  REAL(r_std), ALLOCATABLE, DIMENSION (:)            :: netco2        !! netco2flux
  REAL(r_std), ALLOCATABLE, DIMENSION (:)            :: carblu        !! fco2_land_use
  !-
  !- Declarations for OASIS
  !-
  INTEGER(i_std) :: glo_rank, glo_size ! rank and  number of pe
  INTEGER(i_std) :: loc_rank, loc_size ! rank and  number of pe
  INTEGER(i_std) :: LOCAL_OASIS_COMM  ! local MPI communicator and Initialized
  INTEGER(i_std) :: comp_id    ! component identification
  INTEGER(i_std) :: ierror, flag, oasis_info
  CHARACTER(LEN=4) :: drv_gridname
  CHARACTER(LEN=8) :: varname
  INTEGER(i_std), DIMENSION(3) :: ig_paral
  ! OASIS Output
  INTEGER(i_std) :: il_part_id, tair_id, qair_id, zlevtq_id, zlevuv_id
  INTEGER(i_std) :: rainf_id, snowf_id, swnet_id, lwdown_id, solarang_id
  INTEGER(i_std) :: u_id, v_id, ps_id, cdrag_id
  ! OASIS Input
  INTEGER(i_std) :: vevapp_id, fluxsens_id, fluxlat_id, coastal_id, river_id
  INTEGER(i_std) :: netco2_id, carblu_id, tsolrad_id, tsolnew_id, qsurf_id
  INTEGER(i_std) :: albnir_id, albvis_id, emis_id, z0_id
  !
  INTEGER(i_std), DIMENSION(2) :: var_nodims, var_shape
  INTEGER(i_std) :: nbarg, iret, helpmsg = 0
  CHARACTER(LEN=10) :: arg
  LOGICAL :: initmode = .FALSE.
  !-
  INTEGER(i_std) :: freq_diag=10, debug_lev=1
  INTEGER(i_std) :: w_unit=737
  !
  ! Timer variables
  !
  LOGICAL, PARAMETER :: timemeasure=.TRUE.
  REAL(r_std) :: waitput_cputime=0.0, waitget_cputime=0.0, preparation_cputime=0.0
  REAL(r_std) :: waitput_walltime=0.0, waitget_walltime=0.0, preparation_walltime=0.0
  !
  ! Print point
  !
!!  REAL(r_std), DIMENSION(2) :: testpt=(/44.8,-25.3/)
!!  REAL(r_std), DIMENSION(2) :: testpt=(/44.8,-18.3/)
!!  REAL(r_std), DIMENSION(2) :: testpt=(/-60.25,-5.25/)
!!  REAL(r_std), DIMENSION(2) :: testpt=(/-5.25,41.25/)
  REAL(r_std), DIMENSION(2) :: testpt=(/9999.99,9999.99/)
!!  REAL(r_std), DIMENSION(2) :: testpt=(/46.7,10.3/)
!!  REAL(r_std), DIMENSION(2) :: testpt=(/0.25,49.25/)
  !
  INTEGER iargc, getarg 
  EXTERNAL iargc, getarg 
  !-
  !---------------------------------------------------------------------------------------
  !-
  !- The code has 2 execution mode :
  !-          1) initialisation of the grid file based on the forcing file
  !-          2) Reading the forcing file and sending the data out with OASIS
  !-
  !-
  nbarg = iargc()
  IF ( nbarg > 1 ) THEN
     helpmsg = 1
  ELSE IF ( nbarg == 1 ) THEN
     iret = getarg(1,arg)
     SELECT CASE(arg)
        !
     CASE('-h')
        helpmsg = 1
     CASE('-init')
        initmode = .TRUE.
     CASE DEFAULT
        helpmsg = 1
        !
     END SELECT
  ELSE
     initmode = .FALSE.
  ENDIF
  !
  ! Does the user just want help ?
  !
  IF ( helpmsg > 0 ) THEN
     WRITE(*,*) "USAGE : driver2oasis [-init] " 
     WRITE(*,*) "             The program will read the forcing file provided by variable" 
     WRITE(*,*) "             FORCING_FILE in the run.def file and do one of 2 things :"
     WRITE(*,*) "        "
     WRITE(*,*) "      -init  a grid description file will be generated for the specified "
     WRITE(*,*) "             forcing file. This grid description file will be written into"
     WRITE(*,*) "             the file name provided by the variable GRID_FILE "
     WRITE(*,*) "             of the run.def. "
     WRITE(*,*) "       "
     WRITE(*,*) "             If no arguments are provided driver2oasis will initiate and "
     WRITE(*,*) "             send the forcing data out via OASIS."
     STOP "HELP from driver2oasis"
  ENDIF
  !-
  !- Open output file for driver
  !-
  OPEN(UNIT=w_unit, FILE="out_driver_mono", FORM="formatted")
  !---------------------------------------------------------------------------------------
  !-
  !- Get the general information we need 
  !-
  !---------------------------------------------------------------------------------------
  !-
!!  CALL getin_name("run.def")
  !
  !Config Key   = FORCING_FILE
  !Config Desc  = Name of file containing the forcing data
  !Config If    = [-]
  !Config Def   = forcing_file.nc
  !Config Help  = This is the name of the file which should be opened
  !Config         for reading the forcing data of the dim0 model.
  !Config         The format of the file has to be netCDF and COADS
  !Config         compliant.
  !Config Units = [FILE] 
  !-
  forfilename(:)=" "
  forfilename(1)='forcing_file.nc'
  CALL getin('FORCING_FILE', forfilename)
  !
  !Config Key   = GRID_FILE
  !Config Desc  = Name of file containing the forcing data
  !Config If    = [-]
  !Config Def   = grid_file.nc
  !Config Help  = This is the name of the file from which we will read
  !Config         or write into it the description of the grid from 
  !Config         the forcing file.
  !Config         compliant.
  !Config Units = [FILE] 
  !- 
  gridfilename='grid_file.nc'
  CALL getin('GRID_FILE', gridfilename)
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
  !-
  debug_lev=1
  CALL getin('BAVARD', debug_lev)
  IF ( debug_lev .EQ. 4 ) THEN
     freq_diag=1
  ELSE IF ( debug_lev .EQ. 3 ) THEN
     freq_diag=10
  ELSE IF ( debug_lev .EQ. 2 ) THEN
     freq_diag=24
  ELSE
     freq_diag=96
  ENDIF
  WRITE(w_unit,*) "Debug_lev, freq_diag = ", debug_lev, freq_diag
  !-
  !---------------------------------------------------------------------------------------
  !-
  !- We go into the mode which initialises the grid of the forcing file and writes it
  !- for future usage by the driver and ORCHIDEE.
  !-
  !---------------------------------------------------------------------------------------
  !-
  IF ( initmode ) THEN

     WRITE(w_unit,*) 'Forcing files : ', forfilename(:)
     WRITE(w_unit,*) 'Grid files : ', gridfilename

     nb_forcefile = 0
     DO ik=1,100
        IF ( INDEX(forfilename(ik), '.nc') > 0 ) nb_forcefile = nb_forcefile+1
     ENDDO
     !
     ! This mode of driver2oasis is monoproc and thus we need to force is_root_prc
     !
     is_root_prc = .TRUE.
     !
     CALL forcing_getglogrid (nb_forcefile, forfilename, iim_glo, jjm_glo, nbindex_g, .TRUE.)

     CALL forcing_zoomgrid (zoom_lon, zoom_lat, forfilename(1), .TRUE.)
 
     CALL globgrd_writegrid (gridfilename)

     STOP "Grid file sucessfully generated"
  ENDIF
  !-
  !---------------------------------------------------------------------------------------
  !-
  !- This mode will read the grid file and the forcing file and produce the
  !- data to be handed over to OASIS.
  !-
  !---------------------------------------------------------------------------------------
  !---------------------------------------------------------------------------------------
  !-
  !- Define MPI communicator and set-up OASIS
  !-
  CALL oasis_init_comp(comp_id, comp_name, ierror)
  CALL oasis_get_localcomm(LOCAL_OASIS_COMM, ierror)
  !
  CALL Init_orchidee_para(LOCAL_OASIS_COMM)
  !-
  !-
  !---------------------------------------------------------------------------------------
  !-
  !- Get the grid associated to the forcing file ... it should be generated by
  !- this code with a special option before the OASIS coupled case is run.
  !-
  !---------------------------------------------------------------------------------------
  !-
  IF ( timemeasure ) THEN
     CALL init_timer
     CALL start_timer(timer_global)
  ENDIF
  !
  CALL globgrd_getdomsz(gridfilename, iim_glo, jjm_glo, nbindex_g, model_guess, file_id)
  !-
  !- Allocation of memory
  !- variables over the entire grid (thus in x,y)
  ALLOCATE(lon_glo(iim_glo, jjm_glo))
  ALLOCATE(lat_glo(iim_glo, jjm_glo))
  ALLOCATE(mask_glo(iim_glo, jjm_glo))
  ALLOCATE(area_glo(iim_glo, jjm_glo))
  ALLOCATE(corners_glo(iim_glo, jjm_glo, 4, 2))
  !
  ! Gathered variables
  ALLOCATE(kindex_g(nbindex_g))
  ALLOCATE(contfrac(nbindex_g))
  !-
  !-
  !-
  CALL globgrd_getgrid(file_id, iim_glo, jjm_glo, nbindex_g, model_guess, &
       &               lon_glo, lat_glo, mask_glo, area_glo, corners_glo,&
       &               kindex_g, contfrac, calendar)
  !-
  !- Set the calendar and get some information
  !-
  CALL ioconf_calendar(calendar)
  CALL ioget_calendar(one_year, one_day)
  !-
  !- lalo needs to be created before going into the parallel region
  !-
  ALLOCATE(lalo(nbindex_g,2))
  DO ik=1,nbindex_g
     !
     j = ((kindex_g(ik)-1)/iim_glo)+1
     i = (kindex_g(ik)-(j-1)*iim_glo)
     !
     IF ( i > iim_glo .OR. j > jjm_glo ) THEN
        WRITE(w_unit,*) "Error in the indexing (ik, kindex_g, i, j) : ", ik, kindex_g(ik), i, j
        STOP "ERROR in driver2oasis"
     ENDIF
     !
     lalo(ik,1) = lat_glo(i,j)
     lalo(ik,2) = lon_glo(i,j)
     !
  ENDDO
  !
  WRITE(w_unit,*) "Rank", mpi_rank, " Before opening forcingfile. All land points : ",  nbindex_g
  WRITE(w_unit,*) "Rank", mpi_rank, " from ", iim_glo, " point in Lon. and ", jjm_glo, "in Lat."
  !
  ! Set-up the paralelisation so that all gather and scatters work properly on this monoproc task.
  !
  CALL grid_set_glo(iim_glo, jjm_glo, nbindex_g)
  CALL grid_allocate_glo(4)
  CALL bcast(nbindex_g)
  CALL bcast(kindex_g)
  !
  WRITE(numout,*) "Rank", mpi_rank, "Into Init_orchidee_data_para_driver with ", nbindex_g,index_g(1)
  !
  CALL Init_orchidee_data_para_driver(nbindex_g,kindex_g)
  CALL init_ioipsl_para 
  !-
  !---------------------------------------------------------------------------------------
  !-
  !- Open the forcing file and get the time information.
  !-
  !---------------------------------------------------------------------------------------
  !-
  ! Add w_unit as last arguments in order to get some print_out from the forcing_open routine.
  !
  CALL forcing_open(forfilename, iim_glo,  jjm_glo, lon_glo, lat_glo, nbindex_g, zoom_lon, zoom_lat, &
       &            kindex_g, nbindex_g, w_unit)
  CALL forcing_integration_time(date_start, dt, nbdt)
  !
  !
  ALLOCATE(zlev_tq(nbindex_g), zlev_uv(nbindex_g))
  ALLOCATE(u(nbindex_g), v(nbindex_g), pb(nbindex_g))
  ALLOCATE(temp_air(nbindex_g))
  ALLOCATE(qair(nbindex_g))
  ALLOCATE(ccanopy(nbindex_g))
  ALLOCATE(cdrag(nbindex_g))
  ALLOCATE(precip_rain(nbindex_g))
  ALLOCATE(precip_snow(nbindex_g))
  ALLOCATE(swdown(nbindex_g))
  ALLOCATE(swnet(nbindex_g))
  ALLOCATE(lwdown(nbindex_g))
  ALLOCATE(solarang(nbindex_g))
  ALLOCATE(vevapp(nbindex_g))
  ALLOCATE(fluxsens(nbindex_g))
  ALLOCATE(fluxlat(nbindex_g))
  ALLOCATE(coastalflow(nbindex_g))
  ALLOCATE(riverflow(nbindex_g))
  ALLOCATE(netco2(nbindex_g))
  ALLOCATE(carblu(nbindex_g))
  ALLOCATE(tsol_rad(nbindex_g))
  ALLOCATE(temp_sol_new(nbindex_g))
  ALLOCATE(qsurf(nbindex_g))
  ALLOCATE(albedo(nbindex_g,2))
  ALLOCATE(emis(nbindex_g))
  ALLOCATE(z0(nbindex_g))
  !
  !-
  !---------------------------------------------------------------------------------------
  !-
  !- OASIS Diagnostics
  !-
  !---------------------------------------------------------------------------------------
  !-
  ! Unit for output messages : one file for each process
  !-
  CALL MPI_Comm_Size (LOCAL_OASIS_COMM, loc_size, ierror )
  CALL MPI_Comm_Size (MPI_COMM_WORLD, glo_size, ierror )
  IF (ierror /= 0) THEN
     WRITE(w_unit,*) 'MPI_comm_size abort by model1 compid ',comp_id
  ENDIF
  CALL MPI_Comm_Rank (LOCAL_OASIS_COMM, loc_rank, ierror )
  CALL MPI_Comm_Rank (MPI_COMM_WORLD, glo_rank, ierror )
  IF (ierror /= 0) THEN
     WRITE(0,*) 'MPI_Comm_Rank abort by orchidee comp_id ',comp_id
  ENDIF
  !
  WRITE (w_unit,*) '-----------------------------------------------------------'
  WRITE (w_unit,*) TRIM(comp_name), ' Running with reals compiled as kind =',r_std
  WRITE (w_unit,*) 'I am component ', TRIM(comp_name), ' Local rank :',loc_rank, &
       &           " Global rank : ", glo_rank
  WRITE (w_unit,*) '----------------------------------------------------------'
  CALL flush(w_unit)
  !
  !
  WRITE (w_unit,*) 'DD I am the ', TRIM(comp_name), 'local rank', loc_rank, " with ", iim_glo*jjm_glo, "points."
  WRITE (w_unit,*) 'DD Local number of processors :', loc_size, " Global value :", glo_size
  WRITE (w_unit,*) 'DD Local MPI communicator is :', LOCAL_OASIS_COMM,  MPI_COMM_WORLD
  CALL flush(w_unit)
  !-
  !-
  !---------------------------------------------------------------------------------------
  !-
  !- Send the grid to OASIS ... only on the root processor
  !-
  !---------------------------------------------------------------------------------------
  !-
  IF (loc_rank == 0) THEN
     !
     ! TOCOMPLETE - Put here OASIS grid, corner, areas and mask writing calls !
     !
     drv_gridname = "DRIV"
     !
     CALL oasis_start_grids_writing(flag)
     !
     CALL oasis_write_grid(drv_gridname, iim_glo, jjm_glo, lon_glo, lat_glo)
     !
     ALLOCATE(corners_lon(iim_glo, jjm_glo, 4), corners_lat(iim_glo, jjm_glo, 4))
     corners_lon(:,:,:) = corners_glo(:,:,:,1)
     corners_lat(:,:,:) = corners_glo(:,:,:,2)
     CALL oasis_write_corner(drv_gridname, iim_glo, jjm_glo, 4, corners_lon, corners_lat)
     DEALLOCATE(corners_lon, corners_lat)
     !
     ALLOCATE(maskinv_glo(iim_glo, jjm_glo))
     DO i=1,iim_glo
        DO j=1,jjm_glo
           IF (mask_glo(i,j) == 0) THEN
              maskinv_glo(i,j) = 1
           ELSE
              maskinv_glo(i,j) = 0
           ENDIF
        ENDDO
     ENDDO
     CALL oasis_write_mask(drv_gridname, iim_glo, jjm_glo, maskinv_glo)
     !
     CALL oasis_write_area(drv_gridname, iim_glo, jjm_glo, area_glo)
     !
     CALL oasis_terminate_grids_writing()
  ENDIF
  !
  WRITE(w_unit,*) 'After grids writing'
  call flush(w_unit)
  !-
  !---------------------------------------------------------------------------------------
  !-
  !- Declare the variables
  !-
  !---------------------------------------------------------------------------------------
  !-
  ig_paral(1) = 0
  ig_paral(2) = 0
  ig_paral(3) = nbindex_g

  CALL oasis_def_partition (il_part_id, ig_paral, ierror)

  var_nodims(1) = 1
  var_nodims(2) = 1
  var_shape(1) = 1
  var_shape(1) = nbindex_g
  !
  ! Variables SENT to ORCHIDEE
  ! ==========================
  !
  ! Define levels
  !
  varname="ZLTQDRIV"
  CALL oasis_def_var(zlevtq_id, varname, il_part_id, var_nodims, OASIS_Out, var_shape, OASIS_Real, ierror)
  varname="ZLUVDRIV"
  CALL oasis_def_var(zlevuv_id, varname, il_part_id, var_nodims, OASIS_Out, var_shape, OASIS_Real, ierror)
  !
  ! Define scalar atmospheric variables
  !
  varname="TAIRDRIV"
  CALL oasis_def_var(tair_id, varname, il_part_id, var_nodims, OASIS_Out, var_shape, OASIS_Real, ierror)
  varname="QAIRDRIV"
  CALL oasis_def_var(qair_id, varname, il_part_id, var_nodims, OASIS_Out, var_shape, OASIS_Real, ierror)
  !
  ! Define precipitation variables
  !
  varname="RAINDRIV"
  CALL oasis_def_var(rainf_id, varname, il_part_id, var_nodims, OASIS_Out, var_shape, OASIS_Real, ierror)
  varname="SNOWDRIV"
  CALL oasis_def_var(snowf_id, varname, il_part_id, var_nodims, OASIS_Out, var_shape, OASIS_Real, ierror)
  !
  ! Define radiation variables
  !
  varname="SWNEDRIV"
  CALL oasis_def_var(swnet_id, varname, il_part_id, var_nodims, OASIS_Out, var_shape, OASIS_Real, ierror)
  varname="LWDODRIV"
  CALL oasis_def_var(lwdown_id, varname, il_part_id, var_nodims, OASIS_Out, var_shape, OASIS_Real, ierror)
  varname="SOLADRIV"
  CALL oasis_def_var(solarang_id, varname, il_part_id, var_nodims, OASIS_Out, var_shape, OASIS_Real, ierror)
  !
  ! Finaly pressure and wind
  !
  varname="UWINDRIV"
  CALL oasis_def_var(u_id, varname, il_part_id, var_nodims, OASIS_Out, var_shape, OASIS_Real, ierror)
  varname="VWINDRIV"
  CALL oasis_def_var(v_id, varname, il_part_id, var_nodims, OASIS_Out, var_shape, OASIS_Real, ierror)
  varname="PRESDRIV"
  CALL oasis_def_var(ps_id, varname, il_part_id, var_nodims, OASIS_Out, var_shape, OASIS_Real, ierror)
  varname="DRAGDRIV"
  CALL oasis_def_var(cdrag_id, varname, il_part_id, var_nodims, OASIS_Out, var_shape, OASIS_Real, ierror)
  !
  !
  ! Variables RECEIVED from ORCHIDEE
  ! ================================
  !
  !
  ! Turbulent fluxes
  !
  varname="EVAPDRIV"
  CALL oasis_def_var(vevapp_id, varname, il_part_id, var_nodims, OASIS_In, var_shape, OASIS_Real, ierror)
  !
  varname="SENSDRIV"
  CALL oasis_def_var(fluxsens_id, varname, il_part_id, var_nodims, OASIS_In, var_shape, OASIS_Real, ierror)
  !
  varname="LATEDRIV"
  CALL oasis_def_var(fluxlat_id, varname, il_part_id, var_nodims, OASIS_In, var_shape, OASIS_Real, ierror)
  !
  ! Discharge to the oceans
  !
  varname="COASDRIV"
  CALL oasis_def_var(coastal_id, varname, il_part_id, var_nodims, OASIS_In, var_shape, OASIS_Real, ierror)
  !
  varname="RIVEDRIV"
  CALL oasis_def_var(river_id, varname, il_part_id, var_nodims, OASIS_In, var_shape, OASIS_Real, ierror)
  !
  ! Carbon fluxes
  !
  varname="NECODRIV"
  CALL oasis_def_var(netco2_id, varname, il_part_id, var_nodims, OASIS_In, var_shape, OASIS_Real, ierror)
  !
  varname="LUCODRIV"
  CALL oasis_def_var(carblu_id, varname, il_part_id, var_nodims, OASIS_In, var_shape, OASIS_Real, ierror)
  !
  ! Surface states
  !
  varname="TRADDRIV"
  CALL oasis_def_var(tsolrad_id, varname, il_part_id, var_nodims, OASIS_In, var_shape, OASIS_Real, ierror)
  !
  varname="TNEWDRIV"
  CALL oasis_def_var(tsolnew_id, varname, il_part_id, var_nodims, OASIS_In, var_shape, OASIS_Real, ierror)
  !
  varname="QSURDRIV"
  CALL oasis_def_var(qsurf_id, varname, il_part_id, var_nodims, OASIS_In, var_shape, OASIS_Real, ierror)
  !
  varname="ANIRDRIV"
  CALL oasis_def_var(albnir_id, varname, il_part_id, var_nodims, OASIS_In, var_shape, OASIS_Real, ierror)
  !
  varname="AVISDRIV"
  CALL oasis_def_var(albvis_id, varname, il_part_id, var_nodims, OASIS_In, var_shape, OASIS_Real, ierror)
  !
  varname="EMISDRIV"
  CALL oasis_def_var(emis_id, varname, il_part_id, var_nodims, OASIS_In, var_shape, OASIS_Real, ierror)
  ! 
  varname="ROUGDRIV"
  CALL oasis_def_var(z0_id, varname, il_part_id, var_nodims, OASIS_In, var_shape, OASIS_Real, ierror)
  !
  CALL oasis_enddef (ierror)
  !-
  !---------------------------------------------------------------------------------------
  !-
  !- Get a first set of forcing data
  !-
  !---------------------------------------------------------------------------------------
  !-
  timestep_interval(1) = date_start
  timestep_interval(2) = date_start + (dt/one_day)
  CALL forcing_getvalues(timestep_interval, dt, zlev_tq, zlev_uv, temp_air, qair, &
       &                 precip_rain, precip_snow, swdown, lwdown, solarang, u, v, pb)
  ! An atmsopheric model will provide this information. If zero ORCHIDEE will compute it.
  cdrag(:) = 0.0
  ! Transform the swdown into a swnet
  swnet(:) = 0.13*swdown(:)
  !-
  !---------------------------------------------------------------------------------------
  !-
  !- Go into the time loop
  !-
  !---------------------------------------------------------------------------------------
  !-
  IF ( timemeasure ) THEN
     WRITE(w_unit,*) '------> CPU Time for start-up of driver : ',Get_cpu_Time(timer_global)
     WRITE(w_unit,*) '------> Real Time for start-up of driver : ',Get_real_Time(timer_global)
     CALL stop_timer(timer_global)
     CALL start_timer(timer_global)
  ENDIF
  !-
  DO itau = 0,nbdt-1
     !
     timestep_interval(1) = date_start + itau*(dt/one_day)
     timestep_interval(2) = date_start + (itau+1)*(dt/one_day)
     julian = date_start + (itau+0.5)*(dt/one_day)
     !
     ! Write some diagnostics to look at the shape of the maps
     !
     IF ( itau == 0 ) THEN
       CALL globgrd_writevar(iim_glo, jjm_glo, lon_glo, lat_glo, &
          &                  nbindex_g, lalo, temp_air, "TAIR", "Tair_Driver_0.nc")
     ENDIF
     !
     ! Put first the height of the atmospheric levels
     !
     CALL forcing_printpoint(julian, testpt(1), testpt(2), zlev_tq, "TQ height")
     CALL oasis_put(zlevtq_id, NINT(itau*dt), zlev_tq, oasis_info)
     CALL oasis_put(zlevuv_id, NINT(itau*dt), zlev_uv, oasis_info)
     !
     !
     CALL forcing_printpoint(julian, testpt(1), testpt(2), temp_air, "Air temperature")
     CALL forcing_printpoint(julian, testpt(1), testpt(2), qair, "Air humidity")
     !
     CALL oasis_put(tair_id, NINT(itau*dt), temp_air, oasis_info)
     CALL oasis_put(qair_id, NINT(itau*dt), qair, oasis_info)
     !
     ! Precipitation fluxes
     !
     CALL forcing_printpoint(julian, testpt(1), testpt(2), precip_rain*one_day, "Rainfall")
     CALL forcing_printpoint(julian, testpt(1), testpt(2), precip_snow*one_day, "Snowfall")
     CALL oasis_put(rainf_id, NINT(itau*dt), precip_rain, oasis_info)
     CALL oasis_put(snowf_id, NINT(itau*dt), precip_snow, oasis_info)
     !
     ! Radiation fluxes
     !
     CALL forcing_printpoint(julian, testpt(1), testpt(2), swnet, "Net solar")
     CALL oasis_put(swnet_id, NINT(itau*dt), swnet, oasis_info)
     CALL oasis_put(lwdown_id, NINT(itau*dt), lwdown, oasis_info)
     CALL forcing_printpoint(julian, testpt(1), testpt(2), lwdown, "Downward Longwave") 
     CALL oasis_put(solarang_id, NINT(itau*dt), solarang, oasis_info)
     !
     ! Dynamical fields
     !
     CALL forcing_printpoint(julian, testpt(1), testpt(2), u, "East-ward wind")
     CALL forcing_printpoint(julian, testpt(1), testpt(2), pb, "Surface Pressure")
     CALL oasis_put(u_id, NINT(itau*dt), u, oasis_info)
     CALL oasis_put(v_id, NINT(itau*dt), v, oasis_info)
     CALL oasis_put(ps_id, NINT(itau*dt), pb, oasis_info)
     CALL oasis_put(cdrag_id, NINT(itau*dt), cdrag, oasis_info)
     !
     IF ( timemeasure ) THEN
        waitput_cputime = waitput_cputime + Get_cpu_Time(timer_global)
        waitput_walltime = waitput_walltime + Get_real_Time(timer_global)
        CALL stop_timer(timer_global)
        CALL start_timer(timer_global)
     ENDIF
     !
     ! Get the forcing data for the next time step before we go into the blocking gets
     !
     IF ( itau < (nbdt-1) ) THEN
        timestep_int_next(1) = date_start + (itau+1)*(dt/one_day)
        timestep_int_next(2) = date_start + (itau+2)*(dt/one_day)
        CALL forcing_getvalues(timestep_int_next, dt, zlev_tq, zlev_uv, temp_air, qair, &
             &                 precip_rain, precip_snow, swdown, lwdown, solarang, u, v, pb)
        ! An atmsopheric model will provide this information. If zero ORCHIDEE will compute it.
        cdrag(:) = 0.0
        ! Compute swnet with the albedo computed in the previous call to ORCHIDEE
        swnet(:) = (1.-(albedo(:,1)+albedo(:,2))/2.)*swdown(:)
        !
     ENDIF
     !
     ! Timing of computations for next forcing data
     !
     IF ( timemeasure ) THEN
        preparation_cputime = preparation_cputime + Get_cpu_Time(timer_global)
        preparation_walltime = preparation_walltime + Get_real_Time(timer_global)
        CALL stop_timer(timer_global)
        CALL start_timer(timer_global)
     ENDIF
     !
     ! Go into the blocking gets from OASIS
     !
     CALL oasis_get(vevapp_id, NINT(itau*dt), vevapp, oasis_info)
     CALL oasis_get(fluxsens_id, NINT(itau*dt), fluxsens, oasis_info)
     CALL oasis_get(fluxlat_id, NINT(itau*dt), fluxlat, oasis_info)
     !
     CALL oasis_get(coastal_id, NINT(itau*dt), coastalflow, oasis_info)
     CALL oasis_get(river_id, NINT(itau*dt), riverflow, oasis_info)
     !
     CALL oasis_get(netco2_id, NINT(itau*dt), netco2, oasis_info)
     CALL oasis_get(carblu_id, NINT(itau*dt), carblu, oasis_info)
     !
     CALL oasis_get(tsolrad_id, NINT(itau*dt), tsol_rad, oasis_info)
     CALL oasis_get(tsolnew_id, NINT(itau*dt), temp_sol_new, oasis_info)
     CALL oasis_get(qsurf_id, NINT(itau*dt), qsurf, oasis_info)
     !
     CALL oasis_get(albvis_id, NINT(itau*dt), albedo(:,1), oasis_info)
     CALL oasis_get(albnir_id, NINT(itau*dt), albedo(:,2), oasis_info)
     CALL oasis_get(emis_id, NINT(itau*dt), emis, oasis_info)
     CALL oasis_get(z0_id, NINT(itau*dt), z0, oasis_info)
     CALL forcing_printpoint(julian, testpt(1), testpt(2), vevapp, "Recived evap")
     !
     ! Timing of waits
     !
     IF ( timemeasure ) THEN
        waitget_cputime = waitget_cputime + Get_cpu_Time(timer_global)
        waitget_walltime = waitget_walltime + Get_real_Time(timer_global)
        CALL stop_timer(timer_global)
        CALL start_timer(timer_global)
     ENDIF
     !
     !---------------------------------------------------------------------------------------
     ! Some dianostics
     !---------------------------------------------------------------------------------------
     !
     IF ( MOD(itau, freq_diag) == 0 ) THEN
        !
        CALL forcing_printdate(timestep_interval(1), "Diagnostics START", w_unit)
        !
        WRITE(w_unit,*) MINVAL(temp_sol_new), " << TEMP_SOL_NEW << ", MAXVAL(temp_sol_new)
        WRITE(w_unit,*) MINVAL(vevapp), " << VEVAPP << ", MAXVAL(vevapp)
        WRITE(w_unit,*) MINVAL(fluxsens), " << FLUXSENS << ", MAXVAL(fluxsens)
        WRITE(w_unit,*) MINVAL(fluxlat), " << FLUXLAT << ", MAXVAL(fluxlat)
        WRITE(w_unit,*) MINVAL(coastalflow), " <<  COASTALFLOW << ", MAXVAL(coastalflow)
        WRITE(w_unit,*) MINVAL(riverflow), " << RIVERFLOW << ", MAXVAL(riverflow)
        WRITE(w_unit,*) MINVAL(netco2), " << NETCO2 << ", MAXVAL(netco2)
        WRITE(w_unit,*) MINVAL(carblu), " << CARBLU << ", MAXVAL(carblu)
        WRITE(w_unit,*) MINVAL(tsol_rad), " << TSOL_RAD, << ", MAXVAL(tsol_rad)
        WRITE(w_unit,*) MINVAL(temp_sol_new), " << TEMP_SOL_NEW << ", MAXVAL(temp_sol_new)
        WRITE(w_unit,*) MINVAL(qsurf), " << QSURF << ", MAXVAL(qsurf)
        WRITE(w_unit,*) MINVAL(albedo(:,1)), " << ALBEDO VIS << ", MAXVAL(albedo(:,1))
        WRITE(w_unit,*) MINVAL(albedo(:,2)), " << ALBEDO NIR << ", MAXVAL(albedo(:,2))
        WRITE(w_unit,*) MINVAL(emis), " << EMIS << ", MAXVAL(emis)
        WRITE(w_unit,*) MINVAL(z0), " << Z0 << ", MAXVAL(z0)
        !
        CALL forcing_printdate(timestep_interval(2), "Diagnostics END", w_unit)
        !
     ENDIF
     !
  ENDDO
  !-
  !-
  !-
  IF ( timemeasure ) THEN
     WRITE(w_unit,*) '------> Total CPU Time waiting for put to ORCH : ',waitput_cputime
     WRITE(w_unit,*) '------> Total Real Time waiting for put to ORCH : ',waitput_walltime
     WRITE(w_unit,*) '------> Total CPU Time for preparing forcing : ', preparation_cputime
     WRITE(w_unit,*) '------> Total Real Time for preparing forcing : ', preparation_walltime
     WRITE(w_unit,*) '------> Total CPU Time waiting for get from ORCH : ',waitget_cputime
     WRITE(w_unit,*) '------> Total Real Time waiting for get from ORCH : ',waitget_walltime
     CALL stop_timer(timer_global)
  ENDIF
  !-
  !---------------------------------------------------------------------------------------
  !-
  !- Close OASIS and MPI
  !-
  !---------------------------------------------------------------------------------------
  !-
  CALL oasis_terminate(ierror)
  !
  CALL forcing_close()

  ! Deallocate all variables and reset init flags
  CALL forcing_tools_clear()

  CLOSE(UNIT=w_unit)
  !-
END PROGRAM driver2oasis
