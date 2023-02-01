! =================================================================================================================================
! PROGRAM       : orchideedriveroasis
!
! CONTACT       : jan.polcher@lmd.jussieu.fr
!
! LICENCE      : IPSL (2016)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF      This is the main program for the driver which gets the forcing data from OASIS. It is used as an interface
!!            for WRF for instance. It allows to have another domain decomposition than the atmosphere. Particularly useful
!!            for coupling to atmospheric models not decomposed by the "apple method".
!!
!!\n DESCRIPTION: This program only organises the data and calls sechiba_main after having received the data from OASIS.
!!                The main work for the OASIS interface is done in orchoasis_tools.f90 module.
!!                Call the various modules to get the forcing data and provide it to SECHIBA. The only complexity
!!                is setting-up the domain decomposition and distributing the grid information.
!!                The code is parallel from tip to toe using the domain decomposition inherited from LMDZ.
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S)	:
!!
!! SVN          :
!! $HeadURL:  $
!! $Date:  $
!! $Revision: $
!! \n
!_ ================================================================================================================================
PROGRAM orchideeoasis
  !---------------------------------------------------------------------
  !-
  !- 
  !---------------------------------------------------------------------
  USE defprec
  USE netcdf
  !
  !
  USE ioipsl_para
  USE mod_orchidee_para
  !
  USE grid
  USE timer

  USE globgrd
  USE orchoasis_tools
  !
  USE sechiba
  USE control
  USE ioipslctrl
  !
  USE mod_oasis
  USE mpi
  !-
  IMPLICIT NONE
  !-
  CHARACTER(LEN=80) :: gridfilename
  CHARACTER(LEN=8)  :: model_guess
  INTEGER(i_std)    :: iim_glo, jjm_glo, file_id
  !- 
  INTEGER(i_std)    :: nbseg
  REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: lon_glo, lat_glo, area_glo
  REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: mask_glo
  INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:) :: maskinv_glo
  REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:,:) :: corners_glo
  REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:) :: corners_lon, corners_lat
  INTEGER(i_std) :: nbindex_g, kjpindex
  INTEGER(i_std), ALLOCATABLE, DIMENSION(:) :: kindex, kindex_g
  !
  ! Variables for the global grid available on all procs and used
  ! to fill the ORCHIDEE variable on the root_proc
  !
  REAL(r_std), ALLOCATABLE, DIMENSION(:,:)    :: lalo_glo
  REAL(r_std), ALLOCATABLE, DIMENSION(:)      :: contfrac_glo
  CHARACTER(LEN=20)                           :: calendar
  !-
  !- Variables local to each processors.
  !-
  INTEGER(i_std) :: i, j, ik, nbdt, first_point
  INTEGER(i_std) :: itau, itau_offset, itau_sechiba
  REAL(r_std)    :: date0, date0_shifted, dt, julian, julian0
  INTEGER(i_std) :: rest_id, rest_id_stom
  INTEGER(i_std) ::  hist_id, hist2_id, hist_id_stom, hist_id_stom_IPCC
  REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: lalo_loc
  INTEGER(i_std) :: iim, jjm
  REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: lon, lat
  !-
  !- input fields
  !-
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: u             !! Lowest level wind speed
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: v             !! Lowest level wind speed 
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: zlev_uv       !! Height of first layer
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: zlev_tq       !! Height of first layer
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: qair          !! Lowest level specific humidity
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: precip_rain   !! Rain precipitation
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: precip_snow   !! Snow precipitation
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: lwdown        !! Down-welling long-wave flux 
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: swdown        !! Downwelling surface short-wave flux
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: sinang        !! cosine of solar zenith angle
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: temp_air      !! Air temperature in Kelvin
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: epot_air      !! Air potential energy
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: ccanopy       !! CO2 concentration in the canopy
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: petAcoef      !! Coeficients A from the PBL resolution
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: peqAcoef      !! One for T and another for q
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: petBcoef      !! Coeficients B from the PBL resolution
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: peqBcoef      !! One for T and another for q
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: cdrag         !! Cdrag
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: pb            !! Lowest level pressure
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
  REAL(r_std), ALLOCATABLE, DIMENSION (:,:)          :: veget         !! Fraction of vegetation type (unitless, 0-1)
  REAL(r_std), ALLOCATABLE, DIMENSION (:,:)          :: lai           !! Leaf area index (m^2 m^{-2}
  REAL(r_std), ALLOCATABLE, DIMENSION (:,:)          :: height        !! Vegetation Height (m)
  !-
  !- Declarations for OASIS
  !-
  CHARACTER(LEN=6)   :: comp_name = 'orchid'
  INTEGER(i_std) :: loc_size, glo_size ! number of pes used
  INTEGER(i_std) :: loc_rank, glo_rank ! rank of current pe 
  INTEGER(i_std) :: localComm, LOCAL_OASIS_COMM  ! local MPI communicator and Initialized
  INTEGER(i_std) :: comp_id    ! component identification
  INTEGER(i_std) :: ierror, flag, oasis_info
  CHARACTER(LEN=3) :: chout
  CHARACTER(LEN=4) :: oasis_gridname
  CHARACTER(LEN=30) :: diagfilename
  !-
  !-
  !-
  REAL(r_std) :: atmco2
  REAL(r_std), ALLOCATABLE, DIMENSION (:)  :: u_tq, v_tq, swnet
  LOGICAL :: lrestart_read = .TRUE. !! Logical for _restart_ file to read
  LOGICAL :: lrestart_write = .FALSE. !! Logical for _restart_ file to write'
  !
  ! Time  variables
  !
  LOGICAL, PARAMETER :: timemeasure=.TRUE.
  REAL(r_std) :: waitput_cputime=0.0, waitget_cputime=0.0, orchidee_cputime=0.0
  REAL(r_std) :: waitput_walltime=0.0, waitget_walltime=0.0, orchidee_walltime=0.0
  !
  ! Print point
  !
!!  REAL(r_std), DIMENSION(2) :: testpt=(/44.8,-25.3/)
!!  REAL(r_std), DIMENSION(2) :: testpt=(/44.8,-18.3/)
!!  REAL(r_std), DIMENSION(2) :: testpt=(/-60.25,-5.25/)
!!  REAL(r_std), DIMENSION(2) :: testpt=(/46.7,10.3/)
!!  REAL(r_std), DIMENSION(2) :: testpt=(/0.25,49.25/)
  ! Case when no ouput is desired.
  REAL(r_std), DIMENSION(2) :: testpt=(/9999.99,9999.99/)
  INTEGER(i_std) :: ktest
  !
  !-
  !---------------------------------------------------------------------------------------
  !- 
  !- Define MPI communicator and set-up OASIS
  !- 
  !---------------------------------------------------------------------------------------
  !-
  !
  CALL oasis_init_comp(comp_id, comp_name, ierror)
  CALL oasis_get_localcomm (LOCAL_OASIS_COMM, ierror)
  !
  ! Set parallel processing in ORCHIDEE
  !
  CALL Init_orchidee_para(LOCAL_OASIS_COMM) 
  !====================================================================================
  !
  ! Start timer now that the paralelisation is started.
  !
  IF ( timemeasure ) THEN
     CALL init_timer
     CALL start_timer(timer_global)
     CALL start_timer(timer_mpi)
  ENDIF
  !-
  !
  !---------------------------------------------------------------------------------------
  !-
  !-  Print some diagnostics for this main of ORCHIDEE
  !-
  !---------------------------------------------------------------------------------------
  !-
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, glo_rank, ierror)
  CALL MPI_COMM_RANK(LOCAL_OASIS_COMM, loc_rank, ierror)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, glo_size, ierror)
  CALL MPI_COMM_SIZE(LOCAL_OASIS_COMM, loc_size, ierror)
  !-
  !
  WRITE (numout,*) '-----------------------------------------------------------'
  WRITE (numout,*) '-------------------- comp_name= ',comp_name,'  ---------------'
  WRITE (numout,*) '-----------------------------------------------------------'
  WRITE (numout,*) TRIM(comp_name), ' Running with reals compiled as kind =',r_std
  WRITE (numout,*) 'I am component ', TRIM(comp_name), ' local rank :', loc_rank, &
       &           " Global rank :", glo_rank
  WRITE (numout,*) 'CPUs we are using for the main :', loc_size, ' Global number :', glo_size
  WRITE (numout,*) 'Local and global MPI communicators :', LOCAL_OASIS_COMM,  MPI_COMM_WORLD
  WRITE (numout,*) '----------------------------------------------------------'
  CALL flush(numout)
  !
  !---------------------------------------------------------------------------------------
  !-
  !- Start the getconf processes 
  !-
  !---------------------------------------------------------------------------------------
  !-
!!  CALL getin_name("run.def")
  !-
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
  CALL getin_p('GRID_FILE', gridfilename)
  !-
  !- Get some basic variables from the run.def
  !-
  atmco2=350.
  CALL getin_p('ATM_CO2',atmco2)
  !---------------------------------------------------------------------------------------
  !-
  !- Get the grid on all processors.
  !-
  !---------------------------------------------------------------------------------------
  !-
  !
  CALL globgrd_getdomsz(gridfilename, iim_glo, jjm_glo, nbindex_g, model_guess, file_id)
  nbseg = 4
  !
  !-
  !- Allocation of memory
  !- variables over the entire grid (thus in x,y)
  ALLOCATE(lon_glo(iim_glo, jjm_glo))
  ALLOCATE(lat_glo(iim_glo, jjm_glo))
  ALLOCATE(mask_glo(iim_glo, jjm_glo))
  ALLOCATE(area_glo(iim_glo, jjm_glo))
  ALLOCATE(corners_glo(iim_glo, jjm_glo, nbseg, 2))
  !
  ! Gathered variables
  ALLOCATE(kindex_g(nbindex_g))
  ALLOCATE(contfrac_glo(nbindex_g))
  !-
  !-
  !-
  CALL globgrd_getgrid(file_id, iim_glo, jjm_glo, nbindex_g, model_guess, &
       &               lon_glo, lat_glo, mask_glo, area_glo, corners_glo,&
       &               kindex_g, contfrac_glo, calendar)
  !-
  !- Set the calendar and get some information
  !-
  CALL ioconf_calendar(calendar)
  CALL ioget_calendar(one_year, one_day)
  !-
  !- get the time period for the run
  !-
  CALL orchoasis_time(date0, dt, nbdt)
  !-
  !- lalo needs to be created before going into the parallel region
  !-
  ALLOCATE(lalo_glo(nbindex_g,2))
  DO ik=1,nbindex_g
     !
     j = ((kindex_g(ik)-1)/iim_glo)+1
     i = (kindex_g(ik)-(j-1)*iim_glo)
     !
     IF ( i > iim_glo .OR. j > jjm_glo ) THEN
        WRITE(100+mpi_rank,*) "Error in the indexing (ik, kindex, i, j) : ", ik, kindex(ik), i, j
        STOP "ERROR in orchideeoasis"
     ENDIF
     !
     lalo_glo(ik,1) = lat_glo(i,j)
     lalo_glo(ik,2) = lon_glo(i,j)
     !
  ENDDO
  !
  WRITE(100+mpi_rank,*) "Rank", mpi_rank, " Before parallel region All land points : ",  nbindex_g
  WRITE(100+mpi_rank,*) "Rank", mpi_rank, " from ", iim_glo, " point in Lon. and ", jjm_glo, "in Lat."
  !-
  !---------------------------------------------------------------------------------------
  !-
  !- Initialise the ORCHIDEE domain on the procs (longitude, latitude, indices)
  !-
  !---------------------------------------------------------------------------------------
  !-
  !- init_data_para also transfers kindex_g to index_g (the variable used in ORCHIDEE) 
  !-
  CALL grid_set_glo(iim_glo, jjm_glo, nbindex_g)
  CALL grid_allocate_glo(nbseg)
  ! Copy the list of indexes of land points into index_g used by ORCHIDEE and then broacast to all
  ! processors
  CALL bcast(nbindex_g)
  IF ( is_root_prc) index_g = kindex_g
  CALL bcast(index_g)
  !
  WRITE(numout,*) "Rank", mpi_rank, "Into Init_orchidee_data_para_driver with ", nbindex_g,index_g(1)
  !
  CALL Init_orchidee_data_para_driver(nbindex_g,index_g)
  CALL init_ioipsl_para 
  !
  WRITE(numout,*) "Rank", mpi_rank, "After init_data_para global size : ",  nbp_glo, SIZE(index_g), iim_g, iim_glo, jjm_g, jjm_glo
  WRITE(numout,'("After init_data_para local : ij_nb, jj_nb",2I4)') iim_glo, jj_nb
  !
  ! Allocate grid on the local processor
  !
  IF ( model_guess == "regular") THEN
     CALL grid_init (nbp_loc, nbseg, "RegLonLat", "ForcingGrid")
  ELSE IF ( model_guess == "WRF") THEN
     CALL grid_init (nbp_loc, nbseg, "RegXY", "WRFGrid")
  ELSE
     CALL ipslerr(3, "orchideeoasis", "The grid found in the GRID_FILE is not supported by ORCHIDEE", "", "")
  ENDIF
  !
  ! Transfer the global grid variables to the ORCHIDEE version on the root proc
  ! *_glo -> *_g
  ! Variables *_g were allocated with the CALL init_grid
  !
  IF ( is_root_prc) THEN
     !
     lalo_g(:,:) = lalo_glo(:,:)
     lon_g(:,:) = lon_glo(:,:)
     lat_g(:,:) = lat_glo(:,:)
     !
  ENDIF
  !-
  !
  ! Set the local dimensions of the fields
  !
  iim = iim_glo
  jjm = jj_nb
  kjpindex = nbp_loc
  !
  WRITE(numout,*) mpi_rank, "DIMENSIONS of grid on processor : iim, jjm, kjpindex = ", iim, jjm, kjpindex
  !
  !  Allocate the local arrays we need :
  !
  ALLOCATE(lon(iim,jjm), lat(iim,jjm))
  ALLOCATE(kindex(kjpindex))
  !
  lon=lon_glo(:,jj_para_begin(mpi_rank):jj_para_end(mpi_rank))
  lat=lat_glo(:,jj_para_begin(mpi_rank):jj_para_end(mpi_rank))
  !
  !
  ! Redistribute the indeces on all procs (apple distribution of land points) 
  !
  CALL bcast(lon_g)
  CALL bcast(lat_g)
  CALL scatter(index_g, kindex) 
  ! 
  !
  ! Apply the offset needed so that kindex refers to the index of the land point 
  ! on the current region, i.e. the local lon lat domain.
  !
  kindex(1:kjpindex)=kindex(1:kjpindex)-(jj_begin-1)*iim_glo
  !
  ! This routine transforms the global grid into a series of polygons for all land
  ! points identified by index_g.
  !
  CALL grid_stuff(nbindex_g, iim_g, jjm_g, lon_g, lat_g, index_g, contfrac_glo)
  !
  ! Distribute the global lalo to the local processor level lalo
  !
  ALLOCATE(lalo_loc(kjpindex,2))
  CALL scatter(lalo_glo, lalo_loc)
  lalo(:,:) = lalo_loc(:,:)
  !
  !- 
  !---------------------------------------------------------------------------------------
  !-
  !- Allocate the space for the per processor coupling variables
  !-
  !---------------------------------------------------------------------------------------
  !-
  ALLOCATE(zlev_tq(kjpindex), zlev_uv(kjpindex))
  ALLOCATE(u(kjpindex), v(kjpindex), pb(kjpindex))
  ALLOCATE(temp_air(kjpindex))
  ALLOCATE(qair(kjpindex))
  ALLOCATE(precip_rain(kjpindex), precip_snow(kjpindex))
  ALLOCATE(swdown(kjpindex), lwdown(kjpindex), sinang(kjpindex))
  !
  ALLOCATE(epot_air(kjpindex), ccanopy(kjpindex), cdrag(kjpindex), swnet(kjpindex))
  !
  ALLOCATE(petAcoef(kjpindex), peqAcoef(kjpindex), petBcoef(kjpindex), peqBcoef(kjpindex))
  ALLOCATE(u_tq(kjpindex), v_tq(kjpindex))
  !
  ALLOCATE(vevapp(kjpindex), fluxsens(kjpindex), fluxlat(kjpindex), coastalflow(kjpindex), riverflow(kjpindex))
  ALLOCATE(netco2(kjpindex), carblu(kjpindex))
  ALLOCATE(tsol_rad(kjpindex), temp_sol_new(kjpindex), qsurf(kjpindex))
  ALLOCATE(albedo(kjpindex,2), emis(kjpindex), z0(kjpindex))
  ALLOCATE(veget(kjpindex,nvm))
  ALLOCATE(lai(kjpindex,nvm))
  ALLOCATE(height(kjpindex,nvm))
  !
  WRITE(numout,*) "Rank", mpi_rank, "Domain size per proc !:", iim_glo, jjm_glo, kjpindex, SIZE(kindex)
  WRITE(numout,*) "Rank", mpi_rank, "In parallel region land index starts at : ", kindex(1)
  ! 
  CALL orchoasis_time(date0, dt, nbdt)
  !
  CALL control_initialize
  dt_sechiba = dt
  ! 
  itau = 0
  !
  CALL ioipslctrl_restini(itau, date0, dt, rest_id, rest_id_stom, itau_offset, date0_shifted)
  !
  ! Get the date of the first time step
  !
  julian = date0 + 0.5*(dt/one_day)
  CALL ju2ymds (julian, year, month, day, sec)
  !
  CALL xios_orchidee_init( MPI_COMM_ORCH,                &
       date0,    year,    month,           day,          &
       lon,      lat,     soilth_lev)
  !
  itau_sechiba = itau + itau_offset
  !
  WRITE(numout,*) "orchoasis_history :", iim, jjm, SHAPE(lon)
  !- Initialize IOIPSL sechiba output files
  CALL ioipslctrl_history(iim, jjm, lon, lat,  kindex, kjpindex, itau_sechiba, &
       date0, dt, hist_id, hist2_id, hist_id_stom, hist_id_stom_IPCC)
  WRITE(numout,*) "HISTORY : Define for ", itau, date0, dt
  !
  !-
  !-
  !---------------------------------------------------------------------------------------
  !-
  !- Send the grid to OASIS ... only on the root processor
  !-
  !---------------------------------------------------------------------------------------
  !-
  IF ( is_root_prc ) THEN
     !
     ! TOCOMPLETE - Put here OASIS grid, corner, areas and mask writing calls !
     !
     oasis_gridname = model_guess(1:4)
     CALL oasis_start_grids_writing(flag)
     !
     CALL oasis_write_grid(oasis_gridname, iim_glo, jjm_glo, lon_glo, lat_glo)
     !
     ALLOCATE(corners_lon(iim_glo, jjm_glo, nbseg), corners_lat(iim_glo, jjm_glo, nbseg))
     corners_lon(:,:,:) = corners_glo(:,:,:,1)
     corners_lat(:,:,:) = corners_glo(:,:,:,2)
     CALL oasis_write_corner(oasis_gridname, iim_glo, jjm_glo, nbseg, corners_lon, corners_lat)
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
     CALL oasis_write_mask(oasis_gridname, iim_glo, jjm_glo, maskinv_glo)
     !
     CALL oasis_write_area(oasis_gridname, iim_glo, jjm_glo, area_glo)
     !
     CALL oasis_terminate_grids_writing()
     !
     IF (model_guess == "WRF") THEN
        oasis_gridname = "twrf"
        !
        CALL oasis_start_grids_writing(flag)
        !
        CALL oasis_write_grid(oasis_gridname, iim_glo, jjm_glo, lon_glo, lat_glo)
        !
        CALL oasis_write_corner(oasis_gridname, iim_glo, jjm_glo, nbseg, corners_lon, corners_lat)
        !
        CALL oasis_write_mask(oasis_gridname, iim_glo, jjm_glo, maskinv_glo)
        !
        CALL oasis_write_area(oasis_gridname, iim_glo, jjm_glo, area_glo)
        !
        CALL oasis_terminate_grids_writing()
     ENDIF
     DEALLOCATE(corners_lon, corners_lat)
     !
  ENDIF
  !
  call flush(numout)
  !-
  !---------------------------------------------------------------------------------------
  !-
  !- Declare the variables to be exchanged through OASIS
  !-
  !---------------------------------------------------------------------------------------
  !-
  CALL orchoasis_defvar(mpi_rank, kjpindex)
  !
  !-
  !---------------------------------------------------------------------------------------
  !-
  !- Get a first set of forcing data 
  !-
  !---------------------------------------------------------------------------------------
  !-
  itau = 0
  julian = date0 + (itau+0.5)*(dt/one_day)
  CALL ju2ymds (julian, year, month, day, sec)
  !-
  !---------------------------------------------------------------------------------------
  !-
  !- Going into the time loop
  !-
  !---------------------------------------------------------------------------------------
  !-
  !
  DO itau=1,nbdt
     !
     julian = date0 + (itau-0.5)*(dt/one_day)
     ! Needed for STOMATE but should be taken from the arguments of SECHIBA
     in_julian = itau2date(itau, date0, dt)
     !
     CALL ju2ymds (julian, year, month, day, sec)
     CALL ymds2ju (year,1,1,zero, julian0)
     julian_diff = in_julian-julian0
     !
     !
     IF ( itau == nbdt ) lrestart_write = .TRUE.
     !
     !---------------------------------------------------------------------------------------
     !-
     !- Get the variables from OASIS
     !-
     !---------------------------------------------------------------------------------------
     !
     ! OASIS get call are blocking so they need to be done once all is finsihed
     !
     CALL orchoasis_getvar(itau-1, dt, kjpindex, kindex - (ii_begin - 1), &
          &                zlev_tq, zlev_uv, temp_air, qair, &
          &                precip_rain, precip_snow, swnet, lwdown, sinang, u, v, pb, cdrag)
     !
     CALL orchoasis_printpoint(julian, testpt(1), testpt(2), kjpindex, lalo_loc, temp_air, "RECEIVED Air temperature")
     CALL orchoasis_printpoint(julian, testpt(1), testpt(2), kjpindex, lalo_loc, qair, "RECEIVED Air humidity")
     CALL orchoasis_printpoint(julian, testpt(1), testpt(2), kjpindex, lalo_loc, precip_rain*one_day, "RECEIVED Rainfall")
     CALL orchoasis_printpoint(julian, testpt(1), testpt(2), kjpindex, lalo_loc, precip_snow*one_day, "RECEIVED Snowfall")
     CALL orchoasis_printpoint(julian, testpt(1), testpt(2), kjpindex, lalo_loc, swnet, "RECEIVED net solar")
     CALL orchoasis_printpoint(julian, testpt(1), testpt(2), kjpindex, lalo_loc, lwdown, "RECEIVED lwdown")
     CALL orchoasis_printpoint(julian, testpt(1), testpt(2), kjpindex, lalo_loc, u, "RECEIVED East-ward wind")
     CALL orchoasis_printpoint(julian, testpt(1), testpt(2), kjpindex, lalo_loc, v, "RECEIVED North-ward wind")
     CALL orchoasis_printpoint(julian, testpt(1), testpt(2), kjpindex, lalo_loc, pb, "RECEIVED surface pressure")
     CALL orchoasis_printpoint(julian, testpt(1), testpt(2), kjpindex, lalo_loc, zlev_uv, "RECEIVED UV height")
     CALL orchoasis_printpoint(julian, testpt(1), testpt(2), kjpindex, lalo_loc, zlev_tq, "RECEIVED TQ height")
     CALL orchoasis_printpoint(julian, testpt(1), testpt(2), kjpindex, lalo_loc, sinang, "RECEIVED sinang")
     !
     !
     ! Adaptation of the forcing data to SECHIBA's needs
     !
     ! Contrary to what the documentation says, ORCHIDEE expects surface pressure in hPa.
     pb(:) = pb(:)/100.
     epot_air(:) = cp_air*temp_air(:)+cte_grav*zlev_tq(:)
     ccanopy(:) = atmco2 
     !
     petBcoef(:) = epot_air(:)
     peqBcoef(:) = qair(:)
     petAcoef(:) = zero
     peqAcoef(:) = zero
     !
     CALL orchoasis_printpoint(julian, testpt(1), testpt(2), kjpindex, lalo_loc, z0, &
          &                    "Z0 before compute", ktest)
     !
     ! Interpolate the wind (which is at hight zlev_uv) to the same height 
     ! as the temperature and humidity (at zlev_tq).
     !
     u_tq(:) = u(:)*LOG(zlev_tq(:)/z0(:))/LOG(zlev_uv(:)/z0(:))
     v_tq(:) = v(:)*LOG(zlev_tq(:)/z0(:))/LOG(zlev_uv(:)/z0(:))
     !
     CALL orchoasis_printpoint(julian, testpt(1), testpt(2), kjpindex, lalo_loc, u_tq, "USED East-ward wind")
     !
     IF ( itau .NE. 1 ) THEN
        IF ( timemeasure ) THEN
           waitget_cputime = waitget_cputime + Get_cpu_Time(timer_global)
           waitget_walltime = waitget_walltime + Get_real_Time(timer_global)
           CALL stop_timer(timer_global)
           CALL start_timer(timer_global)
        ENDIF
     ENDIF
     !
     !---------------------------------------------------------------------------------------
     !-
     !- IF first time step : Call to SECHIBA_initialize to set-up ORCHIDEE before doing an actual call
     !- which will provide the first fluxes.
     !-
     !---------------------------------------------------------------------------------------
     !
     itau_sechiba = itau+itau_offset
     !
     ! Update the calendar in xios by sending the new time step
     CALL xios_orchidee_update_calendar(itau_sechiba)
     !
     IF ( itau == 1 ) THEN
        !
        IF ( timemeasure ) THEN
           WRITE(numout,*) '------> CPU Time for start-up of main : ',Get_cpu_Time(timer_global)
           WRITE(numout,*) '------> Real Time for start-up of main : ',Get_real_Time(timer_global)
           CALL stop_timer(timer_global)
           CALL start_timer(timer_global)
        ENDIF
        !
        CALL sechiba_initialize( &
             itau_sechiba, iim*jjm,      kjpindex,      kindex,                    &
             lalo_loc,     contfrac,     neighbours,    resolution,  zlev_tq,       &
             u_tq,         v_tq,         qair,          temp_air,    temp_air,     &
             petAcoef,     peqAcoef,     petBcoef,      peqBcoef,                  &
             precip_rain,  precip_snow,  lwdown,        swnet,       swdown,       &
             pb,           rest_id,      hist_id,       hist2_id,                   &
             rest_id_stom, hist_id_stom, hist_id_stom_IPCC,                         &
             coastalflow,  riverflow,    tsol_rad,      vevapp,      qsurf,        &
             z0,           albedo,       fluxsens,      fluxlat,     emis,         &
             netco2,       carblu,       temp_sol_new,  cdrag)
        CALL orchoasis_printpoint(julian, testpt(1), testpt(2), kjpindex, lalo_loc, temp_sol_new, "Init temp_sol_new")
        !
        ! Use the obtained albedo to diagnose the Downward Solar
        !
        swdown(:) = swnet(:)/(1.-(albedo(:,1)+albedo(:,2))/2.)
        !
        lrestart_read = .FALSE.
        !
        CALL histwrite_p(hist_id, 'LandPoints',  itau+1, (/ REAL(kindex) /), kjpindex, kindex)
        CALL histwrite_p(hist_id, 'Areas',  itau+1, area, kjpindex, kindex)
        CALL histwrite_p(hist_id, 'Contfrac',  itau+1, contfrac, kjpindex, kindex)
        !
        IF ( timemeasure ) THEN
           WRITE(numout,*) '------> CPU Time for set-up of ORCHIDEE : ',Get_cpu_Time(timer_global)
           WRITE(numout,*) '------> Real Time for set-up of ORCHIDEE : ',Get_real_Time(timer_global)
           CALL stop_timer(timer_global)
           CALL start_timer(timer_global)
        ENDIF
        !
     ENDIF
     !
     !---------------------------------------------------------------------------------------
     !-
     !- Call to SECHIBA  
     !-
     !---------------------------------------------------------------------------------------
     !
     CALL sechiba_main (itau_sechiba, iim*jjm, kjpindex, kindex, &
          & lrestart_read, lrestart_write, &
          & lalo_loc, contfrac, neighbours, resolution, &
          ! First level conditions
          & zlev_tq, u_tq, v_tq, qair, qair, temp_air, temp_air, epot_air, ccanopy, &
          ! Variables for the implicit coupling
          & cdrag, petAcoef, peqAcoef, petBcoef, peqBcoef, &
          ! Rain, snow, radiation and surface pressure
          & precip_rain ,precip_snow, lwdown, swnet, swdown, sinang, pb, &
          ! Output : Fluxes
          & vevapp, fluxsens, fluxlat, coastalflow, riverflow, netco2, carblu, &
          ! Surface temperatures and surface properties
          & tsol_rad, temp_sol_new, qsurf, albedo, emis, z0, &
          ! Vegeation, lai and vegetation height
          & veget, lai, height, &
          ! File ids
          & rest_id, hist_id, hist2_id, rest_id_stom, hist_id_stom, hist_id_stom_IPCC)
     !
     ! Use the obtained albedo to diagnose the Downward Solar
     !
     swdown(:) = swnet(:)/(1.-(albedo(:,1)+albedo(:,2))/2.)
     !
     CALL orchoasis_printpoint(julian, testpt(1), testpt(2), kjpindex, lalo_loc, temp_sol_new, "Produced temp_sol_new")
     CALL orchoasis_printpoint(julian, testpt(1), testpt(2), kjpindex, lalo_loc, fluxsens, "Produced fluxsens")
     CALL orchoasis_printpoint(julian, testpt(1), testpt(2), kjpindex, lalo_loc, fluxlat, "Produced fluxlat")
     !
     IF ( timemeasure ) THEN
        orchidee_cputime = orchidee_cputime + Get_cpu_Time(timer_global)
        orchidee_walltime = orchidee_walltime + Get_real_Time(timer_global)
        CALL stop_timer(timer_global)
        CALL start_timer(timer_global)
     ENDIF
     !
     !---------------------------------------------------------------------------------------
     ! Send the ORCHIDEE output back to the driver
     !---------------------------------------------------------------------------------------
     !
     CALL orchoasis_putvar(itau, dt, kjpindex, kindex - (ii_begin - 1),&
          &                vevapp, fluxsens, fluxlat, coastalflow, riverflow, &
          &                netco2, carblu, tsol_rad, temp_sol_new, qsurf, albedo, emis, z0)
     !
     IF ( timemeasure ) THEN
        waitput_cputime = waitput_cputime + Get_cpu_Time(timer_global)
        waitput_walltime = waitput_walltime + Get_real_Time(timer_global)
        CALL stop_timer(timer_global)
        CALL start_timer(timer_global)
     ENDIF
     !
     !---------------------------------------------------------------------------------------
     !-
     !- Write diagnostics
     !-
     !---------------------------------------------------------------------------------------
     !
     !---------------------------------------------------------------------------------------
     !-
     !- Write diagnostics
     !-
     !---------------------------------------------------------------------------------------
     !
     CALL xios_orchidee_send_field("LandPoints" ,(/ ( REAL(ik), ik=1,kjpindex ) /))
     CALL xios_orchidee_send_field("Areas", area)
     CALL xios_orchidee_send_field("Contfrac",contfrac)
     CALL xios_orchidee_send_field("evap",vevapp*one_day/dt_sechiba)
     CALL xios_orchidee_send_field("evap_alma",vevapp/dt_sechiba)
     CALL xios_orchidee_send_field("coastalflow",coastalflow/dt_sechiba)
     CALL xios_orchidee_send_field("riverflow",riverflow/dt_sechiba)
     CALL xios_orchidee_send_field("temp_sol_C",temp_sol_new-ZeroCelsius)
     CALL xios_orchidee_send_field("temp_sol_K",temp_sol_new)
     CALL xios_orchidee_send_field("fluxsens",fluxsens)
     CALL xios_orchidee_send_field("fluxlat",fluxlat)
     CALL xios_orchidee_send_field("tair",temp_air)
     CALL xios_orchidee_send_field("qair",qair)
     CALL xios_orchidee_send_field("q2m",qair)
     CALL xios_orchidee_send_field("t2m",temp_air)
     CALL xios_orchidee_send_field("swnet",swnet)
     CALL xios_orchidee_send_field("swdown",swdown)
     ! zpb in hPa, output in Pa
     CALL xios_orchidee_send_field("Psurf",pb*100.)
     !
     IF ( .NOT. almaoutput ) THEN
        !
        !  ORCHIDEE INPUT variables
        !
        CALL histwrite_p (hist_id, 'swdown',   itau_sechiba, swdown,   kjpindex, kindex)
        CALL histwrite_p (hist_id, 'tair',     itau_sechiba, temp_air, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'qair',     itau_sechiba, qair, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'evap',     itau_sechiba, vevapp, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'coastalflow',itau_sechiba, coastalflow, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'riverflow',itau_sechiba, riverflow, kjpindex, kindex)
        ! 
        CALL histwrite_p (hist_id, 'temp_sol', itau_sechiba, temp_sol_new, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'tsol_max', itau_sechiba, temp_sol_new, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'tsol_min', itau_sechiba, temp_sol_new, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'fluxsens', itau_sechiba, fluxsens, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'fluxlat',  itau_sechiba, fluxlat,  kjpindex, kindex)
        CALL histwrite_p (hist_id, 'swnet',    itau_sechiba, swnet,    kjpindex, kindex)
        CALL histwrite_p (hist_id, 'alb_vis',  itau_sechiba, albedo(:,1), kjpindex, kindex)
        CALL histwrite_p (hist_id, 'alb_nir',  itau_sechiba, albedo(:,2), kjpindex, kindex)
        !
        IF ( hist2_id > 0 ) THEN
           CALL histwrite_p (hist2_id, 'swdown',   itau_sechiba, swdown, kjpindex, kindex)
           CALL histwrite_p (hist2_id, 'tair',     itau_sechiba, temp_air, kjpindex, kindex)
           CALL histwrite_p (hist2_id, 'qair',     itau_sechiba, qair, kjpindex, kindex)
           !
           CALL histwrite_p (hist2_id, 'evap',     itau_sechiba, vevapp, kjpindex, kindex)
           CALL histwrite_p (hist2_id, 'coastalflow',itau_sechiba, coastalflow, kjpindex, kindex)
           CALL histwrite_p (hist2_id, 'riverflow',itau_sechiba, riverflow, kjpindex, kindex)
           ! 
           CALL histwrite_p (hist2_id, 'temp_sol', itau_sechiba, temp_sol_new, kjpindex, kindex)
           CALL histwrite_p (hist2_id, 'tsol_max', itau_sechiba, temp_sol_new, kjpindex, kindex)
           CALL histwrite_p (hist2_id, 'tsol_min', itau_sechiba, temp_sol_new, kjpindex, kindex)
           CALL histwrite_p (hist2_id, 'fluxsens', itau_sechiba, fluxsens, kjpindex, kindex)
           CALL histwrite_p (hist2_id, 'fluxlat',  itau_sechiba, fluxlat,  kjpindex, kindex)
           CALL histwrite_p (hist2_id, 'swnet',    itau_sechiba, swnet,    kjpindex, kindex)
           !
           CALL histwrite_p (hist2_id, 'alb_vis',  itau_sechiba, albedo(:,1), kjpindex, kindex)
           CALL histwrite_p (hist2_id, 'alb_nir',  itau_sechiba, albedo(:,2), kjpindex, kindex)
        ENDIF
     ELSE
        !
        ! Input variables
        !
        CALL histwrite_p (hist_id, 'SinAng', itau_sechiba, sinang, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'LWdown', itau_sechiba, lwdown, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'SWdown', itau_sechiba, swdown, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'Tair', itau_sechiba, temp_air, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'Qair', itau_sechiba, qair, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'SurfP', itau_sechiba, pb, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'Windu', itau_sechiba, u_tq, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'Windv', itau_sechiba, v_tq, kjpindex, kindex)
        !
        CALL histwrite_p (hist_id, 'Evap', itau_sechiba, vevapp, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'SWnet',    itau_sechiba, swnet, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'Qh', itau_sechiba, fluxsens, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'Qle',  itau_sechiba, fluxlat, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'AvgSurfT', itau_sechiba, temp_sol_new, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'RadT', itau_sechiba, temp_sol_new, kjpindex, kindex)
        !
        ! There is a mess with the units passed to the coupler. To be checked with Marc
        !
        IF ( river_routing ) THEN
           CALL histwrite_p (hist_id, 'CoastalFlow',itau_sechiba, coastalflow, kjpindex, kindex)
           CALL histwrite_p (hist_id, 'RiverFlow',itau_sechiba, riverflow, kjpindex, kindex)
        ENDIF
        !
        IF ( hist2_id > 0 ) THEN
           CALL histwrite_p (hist2_id, 'Evap', itau_sechiba, vevapp, kjpindex, kindex)
           CALL histwrite_p (hist2_id, 'SWnet',    itau_sechiba, swnet, kjpindex, kindex)
           CALL histwrite_p (hist2_id, 'Qh', itau_sechiba, fluxsens, kjpindex, kindex)
           CALL histwrite_p (hist2_id, 'Qle',  itau_sechiba, fluxlat, kjpindex, kindex)
           CALL histwrite_p (hist2_id, 'AvgSurfT', itau_sechiba, temp_sol_new, kjpindex, kindex)
           CALL histwrite_p (hist2_id, 'RadT', itau_sechiba, temp_sol_new, kjpindex, kindex)
        ENDIF
     ENDIF
     !
     !
     !
  ENDDO
  !
  !-
  !---------------------------------------------------------------------------------------
  !-
  !- Get time and close IOIPSL, OASIS and MPI
  !-
  !---------------------------------------------------------------------------------------
  !-
  !
  CALL histclo
  IF(is_root_prc) THEN
     CALL restclo
     CALL getin_dump
     !-
     DEALLOCATE(maskinv_glo)
     !
  ENDIF
  !-
  !- Deallocate all variables and reset init flags
  !-
  CALL orchideeoasis_clear()
  !
  WRITE(numout,*) "MPI finalized"
  !-
  IF ( timemeasure ) THEN
     WRITE(numout,*) '------> Total CPU Time waiting to get forcing : ',waitget_cputime
     WRITE(numout,*) '------> Total Real Time waiting to get forcing : ',waitget_walltime
     WRITE(numout,*) '------> Total CPU Time for ORCHIDEE : ', orchidee_cputime
     WRITE(numout,*) '------> Total Real Time for ORCHIDEE : ', orchidee_walltime
     WRITE(numout,*) '------> Total CPU Time waiting to put fluxes : ',waitput_cputime
     WRITE(numout,*) '------> Total Real Time waiting to put fluxes : ',waitput_walltime
     WRITE(numout,*) '------> Total without MPI : CPU Time  : ', Get_cpu_Time(timer_mpi)
     WRITE(numout,*) '------> Total without MPI : Real Time : ', Get_real_Time(timer_mpi)
     CALL stop_timer(timer_global)
     CALL stop_timer(timer_mpi)
  ENDIF
  !-
  CALL oasis_terminate(ierror)
  !
  WRITE(numout,*) "OASIS terminated"
  !-
  !-
CONTAINS
  
!! ================================================================================================================================
!! SUBROUTINE   : orchideeoasis_clear
!!
!>\BRIEF         Clear orchideeoasis
!!
!! DESCRIPTION  :  Deallocate memory and reset initialization variables to there original values
!!
!_ ================================================================================================================================

  SUBROUTINE orchideeoasis_clear

    !- Deallocate all variables existing on all procs (list still incomplete)

    IF ( ALLOCATED(lon_glo) ) DEALLOCATE(lon_glo)
    IF ( ALLOCATED(lat_glo) ) DEALLOCATE(lat_glo)
    IF ( ALLOCATED(mask_glo) ) DEALLOCATE(mask_glo)
    IF ( ALLOCATED(area_glo) ) DEALLOCATE(area_glo)
    IF ( ALLOCATED(corners_glo) ) DEALLOCATE(corners_glo)
    IF ( ALLOCATED(kindex_g) ) DEALLOCATE(kindex_g)
    IF ( ALLOCATED(contfrac_glo) ) DEALLOCATE(contfrac_glo)
    IF ( ALLOCATED(lalo_glo) ) DEALLOCATE(lalo_glo)
    IF ( ALLOCATED(lon) ) DEALLOCATE(lon)
    IF ( ALLOCATED(lat) ) DEALLOCATE(lat)
    IF ( ALLOCATED(kindex) ) DEALLOCATE(kindex)
    IF ( ALLOCATED(lalo_loc) ) DEALLOCATE(lalo_loc)
    IF ( ALLOCATED(zlev_tq) ) DEALLOCATE(zlev_tq)
    IF ( ALLOCATED(zlev_uv) ) DEALLOCATE(zlev_uv)
    IF ( ALLOCATED(u) ) DEALLOCATE(u)
    IF ( ALLOCATED(v) ) DEALLOCATE(v)
    IF ( ALLOCATED(pb) ) DEALLOCATE(pb)
    IF ( ALLOCATED(temp_air) ) DEALLOCATE(temp_air)
    IF ( ALLOCATED(qair) ) DEALLOCATE(qair)
    IF ( ALLOCATED(precip_rain) ) DEALLOCATE(precip_rain)
    IF ( ALLOCATED(precip_snow) ) DEALLOCATE(precip_snow)
    IF ( ALLOCATED(swdown) ) DEALLOCATE(swdown)
    IF ( ALLOCATED(swnet) ) DEALLOCATE(swnet)
    IF ( ALLOCATED(lwdown) ) DEALLOCATE(lwdown)
    IF ( ALLOCATED(sinang) ) DEALLOCATE(sinang)
    IF ( ALLOCATED(epot_air) ) DEALLOCATE(epot_air)
    IF ( ALLOCATED(ccanopy) ) DEALLOCATE(ccanopy)
    IF ( ALLOCATED(cdrag) ) DEALLOCATE(cdrag)
    IF ( ALLOCATED(swnet) ) DEALLOCATE(swnet)
    IF ( ALLOCATED(petAcoef) ) DEALLOCATE(petAcoef)
    IF ( ALLOCATED(peqAcoef) ) DEALLOCATE(peqAcoef)
    IF ( ALLOCATED(petBcoef) ) DEALLOCATE(petBcoef)
    IF ( ALLOCATED(peqBcoef) ) DEALLOCATE(peqBcoef)
    IF ( ALLOCATED(u_tq) ) DEALLOCATE(u_tq)
    IF ( ALLOCATED(v_tq) ) DEALLOCATE(v_tq)
    IF ( ALLOCATED(vevapp) ) DEALLOCATE(vevapp)
    IF ( ALLOCATED(fluxsens) ) DEALLOCATE(fluxsens)
    IF ( ALLOCATED(fluxlat) ) DEALLOCATE(fluxlat)
    IF ( ALLOCATED(coastalflow) ) DEALLOCATE(coastalflow)
    IF ( ALLOCATED(riverflow) ) DEALLOCATE(riverflow)
    IF ( ALLOCATED(netco2) ) DEALLOCATE(netco2)
    IF ( ALLOCATED(carblu) ) DEALLOCATE(carblu)
    IF ( ALLOCATED(tsol_rad) ) DEALLOCATE(tsol_rad)
    IF ( ALLOCATED(temp_sol_new) ) DEALLOCATE(temp_sol_new)
    IF ( ALLOCATED(qsurf) ) DEALLOCATE(qsurf)
    IF ( ALLOCATED(albedo) ) DEALLOCATE(albedo)
    IF ( ALLOCATED(emis) ) DEALLOCATE(emis)
    IF ( ALLOCATED(z0) ) DEALLOCATE(z0)

    CALL sechiba_clear()
    WRITE(numout,*) "Memory deallocated"

  END SUBROUTINE orchideeoasis_clear

END PROGRAM orchideeoasis

