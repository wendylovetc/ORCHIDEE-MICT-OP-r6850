! =================================================================================================================================
! PROGRAM      : testrouting
!
! CONTACT      : jan.polcher _at_ lmd.jussieu.fr
!
! LICENCE      : :)
!
!>\BRIEF       This program tests routing scheme (from routing.f90) which routes the water over the continents
!!             into the oceans and computes the water stored in floodplains or taken for irrigation.
!!
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S) :
!!
!! SVN          :
!! $HeadURL     : svn://forge.ipsl.jussieu.fr/orchidee/trunk/ORCHIDEE/[somewhere]/testrouting.f90 $
!! $Date: 2014-10-27 10:39:00 +0200 (Mon, 27 Oct 2014) $
!! $Revision: XXXX $
!! \n
!_ ================================================================================================================================
PROGRAM testrouting
  !
  USE ioipsl_para
  USE pft_parameters
  USE mod_orchidee_para
  USE control
  USE constantes_soil_var
  USE constantes_var
  USE constantes
  USE time
  USE routing
  USE timer
  USE grid
  !
  USE getlandseamask
  !
  IMPLICIT NONE
  !
  INTEGER(i_std)                                      :: nbseg                            !! 
  !
  INTEGER(i_std)                                      :: iim                              !! Size in longitude of coarser grid
  INTEGER(i_std)                                      :: jjm                              !! Size in latitude of coarser grid
  INTEGER(i_std)                                      :: i, j                             !! Integer variable for loops 
  INTEGER(i_std)                                      :: ibegt, iendt
  INTEGER(i_std)                                      :: ni                               !! For checking nbindex
  REAL(r_std)                                         :: nbyears                          !! Lenght of simulation in years 
  INTEGER(i_std)                                      :: simlen                           !! Lenght of simulation: simlen = 365*48*nbyears
  REAL(r_std)                                         :: dx, dy                           !! Lon/Lat resolution of coarser grid
  REAL(r_std), ALLOCATABLE, DIMENSION (:,:)           :: lon, lat                         !! Lon/lat of coarser grid
  REAL(r_std), ALLOCATABLE, DIMENSION (:,:)           :: orog                             !! New orography after interpolation
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: orog_land, orog_loc
  REAL(r_std), ALLOCATABLE, DIMENSION (:,:)           :: lalo_land
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: contfrac_land, contfrac_loc
  REAL(r_std), ALLOCATABLE, DIMENSION (:,:)           :: contfrac_2d
  !
  REAL(r_std), DIMENSION (1)                          :: lev                              !! Number of level (requested by restini routine) (unitless)
  CHARACTER(LEN=80)                                   :: histname                         !! Name of history file (can not find HISTNAME?)
  INTEGER(i_std)                                      :: hori_id                          !! ID of the default horizontal longitude and latitude map.
  INTEGER(i_std)                                      :: hist_id                          !! History file identification for ???
  INTEGER(i_std)                                      :: rest_id                          !! ID of the restart file
  !
  REAL(r_std)                                         :: date0                            !! Initial date
  REAL(r_std)                                         :: date0_rest                       !! Initial date from restart file
  REAL(r_std)                                         :: date                             !! Current date
  REAL(r_std)                                         :: dt                               !! Same as dtradia ???
  REAL(r_std)                                         :: dw                               !! 86400. ???
  REAL(r_std)                                         :: one_day_loc
  !
  CHARACTER(LEN=40)                                   :: flux_op                          !! Operations to be performed on fluxes
  CHARACTER(LEN=40)                                   :: flux_scinsec                     !! Operation in seconds
  CHARACTER(LEN=40)                                   :: avescatter, once_wrt             !! The various operation to be performed
  !
  ! Input for routing_main
  !
  INTEGER(i_std)                                      :: kjit                             !! Time step number
  INTEGER(i_std)                                      :: nbindex                          !! Number of local continental points 
  REAL(r_std)                                         :: dtradia                          !! Timestep length
  !
  INTEGER(i_std), ALLOCATABLE, DIMENSION (:)          :: kindex_g                         !! Index of land point on 2D map (in local position)
  INTEGER(i_std), ALLOCATABLE, DIMENSION (:)          :: kindex                           !! index of land point per proc
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: runoff                           !! Grid-point runoff (kg/m^2/dt)
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: drainage                         !! Grid-point drainage (kg/m^2/dt)
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: humrel                           !! Soil moisture stress, root extraction potential (unitless)
  REAL(r_std), ALLOCATABLE, DIMENSION (:,:)           :: stempdiag                        !! Diagnostic soil temperature profile
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)       :: totfrac_nobio                    !! Total fraction of continental ice+lakes+cities+...
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)     :: veget_max                        !! Max. fraction of vegetation type (LAI -> infty, unitless)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)       :: floodout                         !! Flow out of floodplains from hydrol
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)     :: transpot                         !! Potential Transpiration (needed for irrigation)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)       :: precip_rain                      !! Rainfall (kg/m^2/dt)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)       :: k_litt                           !! Averaged conductivity for saturated infiltration in the 'litter' layer (kg/m^2/dt)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)       :: reinf_slope                      !! Coefficient which determines the reinfiltration ratio in the grid box due to flat areas (unitless;0-1)
  INTEGER(i_std)                                      :: hist2_id                         !! Access to history file 2 (unitless)
  !
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: evapot_corr                      !! Soil Potential Evaporation
  !
  ! Output from routing_main
  !
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: returnflow                       !! The water flow from lakes and swamps which returns to the grid box (kg/m^2/dt)
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: irrigation                       !! This is the water taken from the reservoirs and beeing put into the upper layers of the soil (kg/m^2/dt)
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: riverflow                        !! Outflow of the major rivers, will be located on the continental grid but this should be a coastal point (kg/dt)
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: coastalflow                      !! Outflow on coastal points by small basins, the water which flows in a disperse way into the ocean (kg/dt)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)       :: reinfiltration                   !! Water flow from ponds and floodplains which returns to the grid box (kg/m^2/dt)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)       :: flood_res                        !! Diagnostic of water amount in the floodplains reservoir (kg)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)       :: flood_frac                       !! Flooded fraction of the grid box (unitless;0-1)
  !
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: soil_deficit                      !! water deficit to reach IRRIG_FULFILL
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: veget          !! Fraction of vegetation type (unitless, 0-1)       
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: vegstress      !! Vegetation moisture stress (only for vegetation growth)
  !
  !
!_ ================================================================================================================================
  !
  !
  CALL Init_orchidee_para() 
  !
  CALL getlandseamask_init(iim, jjm, nbindex)
  ALLOCATE(lon(iim,jjm))
  ALLOCATE(lat(iim,jjm))
  ALLOCATE(orog(iim,jjm))
  ALLOCATE(contfrac_2d(iim,jjm))
  CALL getlandseamask_read(lon, lat, contfrac_2d, orog)
  !
  ! ALLOCATE memory needed
  !
  ALLOCATE(kindex_g(nbindex))
  ALLOCATE(lalo_land(nbindex,2))
  ALLOCATE(contfrac_land(nbindex))
  ALLOCATE(orog_land(nbindex))
  ! 
  !
  !
  ni=0
  DO j=1,jjm
     DO i=1,iim
        IF ( contfrac_2d(i,j) > 0.0 ) THEN
           ni = ni + 1
           IF ( ni .GT. nbindex ) THEN
              WRITE(*,*) "We are expecting ", nbindex, "point."
              WRITE(*,*) "We are at : ", i, j, orog(i,j)
              STOP 'Too many continental points'
           ENDIF
           kindex_g(ni) = (j-1)*iim + i
           lalo_land(ni,1) = lat(i,j)
           lalo_land(ni,2) = lon(i,j)
           contfrac_land(ni) = contfrac_2d(i,j)
           orog_land(ni) = orog(i,j)
        ENDIF
     ENDDO
  ENDDO  
  !
  !
  nbseg = 4
  !
  ! 
  CALL grid_set_glo(iim, jjm, nbindex)
  CALL grid_allocate_glo(nbseg)
  !
  CALL bcast(nbindex)
  ALLOCATE(index_g(nbindex))
  IF ( is_root_prc ) index_g(:) = kindex_g(:)
  CALL bcast(index_g)
  !
  WRITE(*,*) "GOING INTO Init_orchidee_data_para_driver", nbindex, index_g(1), SIZE(kindex_g)
  CALL Init_orchidee_data_para_driver(nbindex, index_g)
  WRITE(*,*) "OUT OF Init_orchidee_data_para_driver"
  CALL init_ioipsl_para 
  !
  WRITE(*,*) mpi_rank, "DIMENSIONS of grid on processor : iim, jjm, nbindex = ", iim, jjm, nbindex, nbp_loc
  !
  CALL grid_init (nbp_loc, nbseg, "RegLonLat", "ForcingGrid")
  !
  !==========================================================================
  !
  ! Transfer the global grid variables to the root proc
  ! *_glo -> *_g
  ! Variables *_g were allocated with the CALL init_grid
  !
  IF ( is_root_prc) THEN
     !
     lalo_g(:,:) = lalo_land(:,:)
     lon_g(:,:) = lon(:,:)
     lat_g(:,:) = lat(:,:)
     contfrac_g(:) = contfrac_land(:)
     !
  ENDIF
  !
  CALL grid_stuff(nbindex, iim, jjm, lon_g, lat_g, index_g)
  !
  !
  ! Distribute the grid to all processors
  !
  ! Redistribute the indeces on all procs (apple distribution of land points) 
  !
  ALLOCATE(kindex(nbp_loc))
  ALLOCATE(orog_loc(nbp_loc), contfrac_loc(nbp_loc))
  CALL bcast(lon_g)
  CALL bcast(lat_g)
  CALL scatter(index_g, kindex) 
  CALL scatter(lalo_land, lalo)
  CALL scatter(orog_land, orog_loc)
  CALL scatter(contfrac_land, contfrac_loc)
  ! 
  !
  ! Apply the offset needed so that kindex refers to the index of the land point 
  ! on the current region, i.e. the local lon lat domain.
  !
  kindex(1:nbp_loc)=kindex(1:nbp_loc)-(jj_begin-1)*iim
  !
  !
  !==========================================================================================
  !
  ! The grid is in place and we can start to prepare the time of integration.
  !
  CALL ioconf_calendar("gregorian")
  !
  !
  ! Determine initial step or restart one
  !
  !Config Key   = RESTART_IN
  !Config Desc  = Name of restart file to read at restart
  !Config If    = [-]
  !Config Def   = NONE
  !Config Help  = This function allows to select a restart file with which
  !               the simulation will be initialized.
  !Config Units = [-] 
  !- 
  !
  restname_in="NONE"
  CALL getin('RESTART_IN', restname_in)
  !
  !Config Key   = RESTART_OUT
  !Config Desc  = Name of restart file to be written at the end of the simulation.
  !Config If    = [-]
  !Config Def   = NONE
  !Config Help  = This function allows to select a restart file which will be written by the
  !               model and which can be used as input for a future restart.
  !Config Units = [-] 
  !- 
  !
  restname_out="restart_out.nc"
  CALL getin('RESTART_OUT', restname_out)
  !
  !Config Key   = SIMULATION_LEN
  !Config Desc  = Time step length in years for "testrouting"
  !Config If    = [-]
  !Config Def   = 1
  !Config Help  = This is time step length for testrouting
  !Config Units = [-] 
  nbyears=1.0
  CALL getin('SIMULATION_LEN', nbyears)
  !
  !Config Key   = DTRADIA
  !Config Desc  = Time step length for "testrouting"
  !Config If    = [-]
  !Config Def   = 1800.
  !Config Help  = This is time step length for testrouting
  !Config Units = [-] 
  !- 
  !DTRADIA = 1800.
  dtradia = 1800.
  CALL getin('DTRADIA', dtradia)
  !
  !- Initial date
  CALL ymds2ju (2000,1,1,0.0, date0)
  date0_rest = date0
  !
  dt = dtradia
  dw = 86400.
  !
  CALL control_initialize
  ! 
  CALL ioget_calendar(one_year,one_day_loc)
  !
  !  We have all we need and we can start to work
  !
  !
  IF (is_root_prc) THEN
     CALL restini(restname_in, iim, jjm, lon_g, lat_g, 1, lev, &
          &  restname_out, ibegt, date0_rest, dtradia, rest_id, .FALSE.)
  ELSE
     rest_id=0
  ENDIF
  CALL bcast (ibegt)
  CALL bcast (date0_rest)
  CALL bcast (dtradia)
  !
  !
  IF ( INDEX(restname_in, "NONE") > 0 ) THEN
     kjit = 1
     ibegt = 1
  ELSE
     kjit = ibegt
     date0 = date0_rest
  ENDIF
  WRITE(*,*) 'Out of restini : kjit=',kjit, " ibegt=", ibegt, " date0=", date0
  !
  !- time step length
  !
  !  Set up the history file
  !
  !Config Key   = HISTNAME
  !Config Desc  = Name of the history file
  !Config If    = [-]
  !Config Def   = out_testrouting
  !Config Help  = The name of the file which will contain all the diagnostics of routing.
  !Config Units = [-] 
  !-
  WRITE(flux_op,'("ave(scatter(X*",F8.1,"))")') one_day_loc/dt
  WRITE(flux_scinsec,'("ave(scatter(X*",F8.6,"))")') 1.0/dt
  avescatter = 'ave(scatter(X))'
  once_wrt = 'once(scatter(X))'
  ! 
  histname="out_testrouting"
  CALL getin('HISTNAME', histname)
  !
  CALL histbeg(histname, iim, lon, jjm, lat,  1, iim, 1, jjm, &
       &     kjit-1, date0, dtradia, hori_id, hist_id, domain_id=orch_domain_id)
  !
  CALL histdef(hist_id, 'Orog', 'Orography', ' ', &
       & iim, jjm, hori_id, 1,1,1, -99, 32, once_wrt, dt, dw)
  CALL histdef(hist_id, 'Contfrac', 'Fraction of continent', ' ', &
       & iim, jjm, hori_id, 1,1,1, -99, 32, once_wrt, dt, dw)
  CALL histdef(hist_id, 'Areas', 'Mesh areas', 'm2', &
       & iim,jjm, hori_id, 1,1,1, -99, 32, once_wrt, dt, dw)
  !
  CALL histdef(hist_id, 'riversret', 'Return from endorheic rivers', 'mm/d', &
       & iim,jjm, hori_id, 1,1,1, -99, 32, flux_op, dt,dw)
  CALL histdef(hist_id, 'hydrographs', 'Hydrographs of gridbox outflow', 'm^3/s', &
       & iim,jjm, hori_id, 1,1,1, -99, 32, flux_scinsec, dt,dw)
  !
  CALL histdef(hist_id, 'fastr', 'Fast flow reservoir', 'kg/m^2', &
       & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter, dt,dw)
  CALL histdef(hist_id, 'slowr', 'Slow flow reservoir', 'kg/m^2', &
       & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter, dt,dw)
  CALL histdef(hist_id, 'streamr', 'Stream flow reservoir', 'kg/m^2', &
       & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter, dt,dw)
  CALL histdef(hist_id, 'pondr', 'Volume in pond reservoir', 'kg/m^2', &
       & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter, dt,dw)
  CALL histdef(hist_id, 'lakevol', 'Volume in lake reservoir', 'kg/m^2', &
       & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter, dt,dw)
  !
  CALL histdef(hist_id, 'basinmap', 'Aproximate map of the river basins', ' ', &
       & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter, dt,dw) 
  CALL histdef(hist_id, 'nbrivers', 'Number or rivers in the outflow grid box', ' ', &
       & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter, dt,dw)
  !
  CALL histend(hist_id)
  !
  ! Put a copy of the orography into the restart
  !
  !
  CALL histwrite_p(hist_id, 'Orog', kjit+1, orog_loc, nbp_loc, kindex)
  CALL histwrite_p(hist_id, 'Contfrac', kjit+1, contfrac_loc, nbp_loc, kindex)
  CALL histwrite_p(hist_id, 'Areas',  kjit+1, area, nbp_loc, kindex)
  !
  ! Override some settings as testrouting is not as flexible as the rest of ORCHIDEE.
  !
  hist2_id=-1
  almaoutput=.FALSE.
  !================================================================================
  !
  ! Set up the routing schemes
  !
  !
  !  Allocate all the physical variables
  !
  ! Input variables
  ALLOCATE(runoff(nbp_loc), drainage(nbp_loc), humrel(nbp_loc), stempdiag(nbp_loc,nslm), transpot(nbp_loc,nvmc))
  ALLOCATE(totfrac_nobio(nbp_loc), veget_max(nbp_loc,nvmc), floodout(nbp_loc))
  ALLOCATE(precip_rain(nbp_loc), k_litt(nbp_loc),reinf_slope(nbp_loc)) 
  ALLOCATE(evapot_corr(nbp_loc))
  ! Output variables
  ALLOCATE(returnflow(nbp_loc), irrigation(nbp_loc), riverflow(nbp_loc), coastalflow(nbp_loc), reinfiltration(nbp_loc))
  ALLOCATE(flood_frac(nbp_loc), flood_res(nbp_loc))
  !
  ALLOCATE(vegstress(nbp_loc, nvmc))
  ALLOCATE(veget(nbp_loc, nvmc))
  ALLOCATE(soil_deficit(nbp_loc, nvmc))
  !
  ! Get some fake value for input arrays
  !
  runoff(:) = 1.0
  drainage(:) = 1.0
  humrel(:) = 0.75
  stempdiag(:,:) = 273.5
  transpot(:,:)=0.0
  reinfiltration(:)=0.0
  flood_frac(:)=0.0
  flood_res(:)=0.0
  totfrac_nobio = 0.1
  veget_max = 0.2
  floodout = 0.0
  precip_rain = 0.0
  k_litt = 0.0 
  reinf_slope = 0.1
  evapot_corr(:) = 10.0
  !
  vegstress = 0
  veget = 0
  soil_deficit = 0
  !
  ! 
  !
  !
  CALL routing_initialize( kjit,        nbp_loc,        kindex, &
       &                   rest_id,     hist_id,        hist2_id,   lalo, &
       &                   neighbours,  resolution,     contfrac,   stempdiag, &
       &                   returnflow,  reinfiltration, irrigation, riverflow, &
       &                   coastalflow, flood_frac,     flood_res)
  !
  ! Do loop over a number of time-steps
  !
  simlen = NINT(nbyears*365*one_day_loc/dtradia)
  ibegt=kjit
  iendt=kjit+simlen
  WRITE(*,*) "The simulation will go from ", ibegt, " to ", iendt
  !
  DO kjit = ibegt,iendt
     !
     date = date0 + (kjit-1)*(dtradia/one_day_loc)
     !
     IF ( date < date0+1 ) THEN
        ! During one day one kg/m^2d divided up in runoff and drainage
        runoff(:) = 0.5/48.
        drainage(:) = 0.5/48.
     ELSE
        runoff(:) = 0.0
        drainage(:) = 0.0
     ENDIF
     !
     vegstress = 0
     veget = 0
     evapot_corr = 0
     soil_deficit = 0
     ! 
     CALL routing_main(kjit, nbp_loc, kindex, &
           & lalo, neighbours, resolution, contfrac, totfrac_nobio, veget, veget_max, soil_deficit, floodout, runoff, &
           & drainage, transpot, evapot_corr, vegstress, precip_rain, humrel, k_litt, flood_frac, flood_res, &
           & stempdiag, reinf_slope, returnflow, reinfiltration, irrigation, riverflow, coastalflow, rest_id, hist_id, hist2_id)
     !
     WRITE(*,*) "Out of routing at time step = ",kjit,' Seconds since start', (kjit-ibegt)*dtradia
     !
  ENDDO
  !
  ! Shut everything down
  !
  CALL routing_finalize(kjit, nbp_loc, rest_id, flood_frac, flood_res)
  !
  !
  CALL histclo
  IF ( is_root_prc) THEN
     CALL restclo
  ENDIF
  !
  CALL Finalize_mpi
  !
END PROGRAM testrouting
!
