! ================================================================================================================================
!  MODULE       : intersurf
!
!  CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
!  LICENCE      : IPSL (2006)
!  This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF   Subroutines for interfacing between the driver dim2_driver and sechiba and between LMDZ and sechiba
!!
!!\n DESCRIPTION:    This module contains subroutines for the interfacing between ORCHIDEE and LMDZ and between the driver 
!!                   dim2_driver and the reste of ORCHIDEE. Follwoing subroutines exists:
!!
!!              - intersurf_initialize_2d       : Intitialization and call to sechiba_initialize, called from dim2_driver 
!!                                                before the first call to intersurf_main_2d.
!!              - intersurf_main_2d             : Main subroutine will call sechiba_main, called from dim2_driver at each time-step
!!              - init_intersurf                : Initialize grid information when coupling with LMDZ. This subroutine is called 
!!                                                from LMDZ before call to intersurf_initialize_gathered
!!              - intersurf_initialize_gathered : Intitialization and call to sechiba_initialize, called from LMDZ
!!                                                before the first call to intersurf_main_2d.
!!              - intersurf_main_gathered       : Main subroutine will call sechiba_main, called from LMDZ at each time-step
!!              - intsurf_time                  : Initialize and update time information, called in the initialization phase 
!!                                                and at each time step from the different intersurf subroutines. 
!!
!! RECENT CHANGE(S): 
!!
!! REFERENCE(S)	: None
!!
!! SVN          :
!! $HeadURL$
!! $Date$
!! $Revision$
!! \n
!_ ================================================================================================================================
MODULE intersurf

  USE IOIPSL
  USE xios_orchidee
  USE ioipsl_para 
  USE defprec
  USE sechiba
  USE constantes
  USE constantes_soil
  USE control
  USE pft_parameters
  USE mod_orchidee_para
  USE solar
  USE grid
  USE time, ONLY : one_year, one_day, dt_sechiba, time_initialize, time_nextstep, julian_diff
  USE time, ONLY : year_end, month_end, day_end, sec_end
  USE thermosoilc, ONLY : thermosoilc_levels
  USE ioipslctrl, ONLY : ioipslctrl_history, ioipslctrl_restini, ok_histsync, max_hist_level, dw

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: Init_intersurf, intersurf_main, intersurf_main_2d, intersurf_main_gathered
  PUBLIC :: intersurf_initialize_2d, intersurf_initialize_gathered, intersurf_clear


  ! Interface called from LMDZ in older verisons
  INTERFACE intersurf_main
    MODULE PROCEDURE intersurf_main_gathered
  END INTERFACE

  LOGICAL, SAVE                                      :: lstep_init_intersurf=.TRUE.!! Initialisation has to be done one time
!$OMP THREADPRIVATE(lstep_init_intersurf)
  INTEGER(i_std), SAVE                               :: printlev_loc            !! Write level to this module
!$OMP THREADPRIVATE(printlev_loc)
  INTEGER(i_std), SAVE                               :: hist_id, rest_id        !! IDs for history and restart files
!$OMP THREADPRIVATE(hist_id, rest_id)
  INTEGER(i_std), SAVE                               :: hist2_id                !! ID for the second history files (Hi-frequency ?)
!$OMP THREADPRIVATE(hist2_id)
  INTEGER(i_std), SAVE                               :: hist_id_stom, hist_id_stom_IPCC, rest_id_stom !! Dito for STOMATE
!$OMP THREADPRIVATE(hist_id_stom, hist_id_stom_IPCC, rest_id_stom)
  INTEGER(i_std), SAVE                               :: itau_offset  !! This offset is used to phase the 
                                                                     !! calendar of the GCM or the driver.
!$OMP THREADPRIVATE(itau_offset)
  REAL(r_std), SAVE                                  :: date0_shifted
!$OMP THREADPRIVATE(date0_shifted)
  REAL(r_std), SAVE :: julian0                       !! first day of this year
!$OMP THREADPRIVATE(julian0)
  LOGICAL, SAVE                                      :: ok_q2m_t2m=.TRUE. !! Flag ok if the variables are present in the call in gcm.
!$OMP THREADPRIVATE(ok_q2m_t2m)
  LOGICAL, SAVE                                      :: fatmco2           !! Flag to force the value of atmospheric CO2 for vegetation.
!$OMP THREADPRIVATE(fatmco2)
  REAL(r_std), SAVE                                  :: atmco2            !! atmospheric CO2 
!$OMP THREADPRIVATE(atmco2)
  REAL(r_std), SAVE                                  :: coeff_rel         !! Coefficient for time filter on riverflow and coastalflow
!$OMP THREADPRIVATE(coeff_rel)
  REAL(r_std), ALLOCATABLE, DIMENSION(:), SAVE       :: riverflow_cpl0    !! Value from previous time step for riverflow
!$OMP THREADPRIVATE(riverflow_cpl0)
  REAL(r_std), ALLOCATABLE, DIMENSION(:), SAVE       :: coastalflow_cpl0  !! Value from previous time step for coastalflow
!$OMP THREADPRIVATE(coastalflow_cpl0)
  INTEGER(i_std), SAVE                               :: nb_fields_in=-1   !!  Number of fields to give to the GCM
!$OMP THREADPRIVATE(nb_fields_in)  
  INTEGER(i_std), SAVE                               :: nb_fields_out=-1  !!  Number of fields to get from the GCM
!$OMP THREADPRIVATE(nb_fields_out)  

  PUBLIC lstep_init_intersurf
  
CONTAINS

!!  =============================================================================================================================
!! SUBROUTINE:    intersurf_initialize_2d
!!
!>\BRIEF	  Initialization and call to sechiba_initialize
!!
!! DESCRIPTION:	  Initialization of module variables, read options from parameter file, initialize output files and call to 
!!                sechiba_initialize.
!!
!!                This subroutine is called from dim2_driver before the first call to intersurf_main_2d.
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE intersurf_initialize_2d (kjit, iim, jjm, kjpindex, kindex, xrdt, &
       lrestart_read, lrestart_write, lon, lat, zcontfrac, zresolution, date0, &
       zlev, u, v, qair, temp_air, epot_air, ccanopy, &
       cdrag, petAcoef, peqAcoef, petBcoef, peqBcoef, &
       precip_rain, precip_snow, lwdown, swnet, swdown, pb, &
       vevapp, fluxsens, fluxlat, coastalflow, riverflow, &
       tsol_rad, temp_sol_new, qsurf, albedo, emis, z0m)

    IMPLICIT NONE

    !! 0. Variable and parameter declaration
    !! 0.1 Input variables
    INTEGER(i_std),INTENT (in)                            :: kjit            !! Time step number
    INTEGER(i_std),INTENT (in)                            :: iim, jjm        !! Dimension of input fields
    INTEGER(i_std),INTENT (in)                            :: kjpindex        !! Number of continental points
    REAL(r_std),INTENT (in)                               :: xrdt            !! Time step in seconds
    LOGICAL, INTENT (in)                                 :: lrestart_read    !! Logical for _restart_ file to read
    LOGICAL, INTENT (in)                                 :: lrestart_write   !! Logical for _restart_ file to write'
    REAL(r_std), INTENT (in)                              :: date0           !! Date at which kjit = 0
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)      :: kindex          !! Index for continental points
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: u             !! Lowest level wind speed
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: v             !! Lowest level wind speed 
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: zlev          !! Height of first layer
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: qair          !! Lowest level specific humidity
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: precip_rain   !! Rain precipitation
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: precip_snow   !! Snow precipitation
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: lwdown        !! Down-welling long-wave flux 
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: swnet         !! Net surface short-wave flux
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: swdown        !! Downwelling surface short-wave flux
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: temp_air      !! Air temperature in Kelvin
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: epot_air      !! Air potential energy
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: ccanopy       !! CO2 concentration in the canopy
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: petAcoef      !! Coeficients A from the PBL resolution
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: peqAcoef      !! One for T and another for q
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: petBcoef      !! Coeficients B from the PBL resolution
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: peqBcoef      !! One for T and another for q
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: pb            !! Surface pressure
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: lon, lat      !! Geographical coordinates
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: zcontfrac     !! Fraction of continent in the grid
    REAL(r_std),DIMENSION (iim,jjm,2), INTENT(in)           :: zresolution   !! resolution in x and y dimensions

    !! 0.2 Output variables
    REAL(r_std),DIMENSION (iim,jjm), INTENT(out)            :: z0m            !! Surface roughness
    REAL(r_std),DIMENSION (iim,jjm), INTENT(out)            :: coastalflow   !! Diffuse flow of water into the ocean (m^3/s)
    REAL(r_std),DIMENSION (iim,jjm), INTENT(out)            :: riverflow     !! Largest rivers flowing into the ocean (m^3/s)
    REAL(r_std),DIMENSION (iim,jjm), INTENT(out)            :: tsol_rad      !! Radiative surface temperature
    REAL(r_std),DIMENSION (iim,jjm), INTENT(out)            :: vevapp        !! Total of evaporation
    REAL(r_std),DIMENSION (iim,jjm), INTENT(out)            :: temp_sol_new  !! New soil temperature
    REAL(r_std),DIMENSION (iim,jjm), INTENT(out)            :: qsurf         !! Surface specific humidity
    REAL(r_std),DIMENSION (iim,jjm,2), INTENT(out)          :: albedo        !! Albedo
    REAL(r_std),DIMENSION (iim,jjm), INTENT(out)            :: fluxsens      !! Sensible chaleur flux
    REAL(r_std),DIMENSION (iim,jjm), INTENT(out)            :: fluxlat       !! Latent chaleur flux
    REAL(r_std),DIMENSION (iim,jjm), INTENT(out)            :: emis          !! Emissivity

    !! 0.3 Modified variables
    REAL(r_std),DIMENSION (iim,jjm), INTENT(inout)          :: cdrag         !! Cdrag

    !! 0.4 Local variables
    REAL(r_std),DIMENSION (kjpindex)                      :: zu            !! Work array to keep u
    REAL(r_std),DIMENSION (kjpindex)                      :: zv            !! Work array to keep v
    REAL(r_std),DIMENSION (kjpindex)                      :: zzlev         !! Work array to keep zlev
    REAL(r_std),DIMENSION (kjpindex)                      :: zqair         !! Work array to keep qair
    REAL(r_std),DIMENSION (kjpindex)                      :: zprecip_rain  !! Work array to keep precip_rain
    REAL(r_std),DIMENSION (kjpindex)                      :: zprecip_snow  !! Work array to keep precip_snow
    REAL(r_std),DIMENSION (kjpindex)                      :: zlwdown       !! Work array to keep lwdown
    REAL(r_std),DIMENSION (kjpindex)                      :: zswnet        !! Work array to keep swnet
    REAL(r_std),DIMENSION (kjpindex)                      :: zswdown       !! Work array to keep swdown
    REAL(r_std),DIMENSION (kjpindex)                      :: ztemp_air     !! Work array to keep temp_air
    REAL(r_std),DIMENSION (kjpindex)                      :: zepot_air     !! Work array to keep epot_air
    REAL(r_std),DIMENSION (kjpindex)                      :: zccanopy      !! Work array to keep ccanopy
    REAL(r_std),DIMENSION (kjpindex)                      :: zpetAcoef     !! Work array to keep petAcoef
    REAL(r_std),DIMENSION (kjpindex)                      :: zpeqAcoef     !! Work array to keep peqAcoef
    REAL(r_std),DIMENSION (kjpindex)                      :: zpetBcoef     !! Work array to keep petBcoef
    REAL(r_std),DIMENSION (kjpindex)                      :: zpeqBcoef     !! Work array to keep peqVcoef
    REAL(r_std),DIMENSION (kjpindex)                      :: zcdrag        !! Work array to keep cdrag
    REAL(r_std),DIMENSION (kjpindex)                      :: zpb           !! Work array to keep pb
    REAL(r_std),DIMENSION (kjpindex)                      :: zz0m, zz0h    !! Work array to keep zz0m, zz0h
    REAL(r_std),DIMENSION (kjpindex)                      :: zcoastal      !! Work array to keep coastalflow (m^3/dt)
    REAL(r_std),DIMENSION (kjpindex)                      :: zriver        !! Work array to keep riverflow (m^3/dt)
    REAL(r_std),DIMENSION (kjpindex)                      :: znetco2       !! Work array to keep netco2flux
    REAL(r_std),DIMENSION (kjpindex)                      :: zcarblu       !! Work array to keep fco2_land_use
    REAL(r_std),DIMENSION (kjpindex)                      :: ztsol_rad     !! Work array to keep tsol_rad
    REAL(r_std),DIMENSION (kjpindex)                      :: zvevapp       !! Work array to keep vevapp
    REAL(r_std),DIMENSION (kjpindex)                      :: ztemp_sol_new !! Work array to keep temp_sol_new
    REAL(r_std),DIMENSION (kjpindex)                      :: zqsurf        !! Work array to keep qsurf
    REAL(r_std),DIMENSION (kjpindex,2)                    :: zalbedo       !! Work array to keep albedo
    REAL(r_std),DIMENSION (kjpindex)                      :: zfluxsens     !! Work array to keep fluxsens
    REAL(r_std),DIMENSION (kjpindex)                      :: zfluxlat      !! Work array to keep fluxlat
    REAL(r_std),DIMENSION (kjpindex)                      :: zemis         !! Work array to keep emis
    REAL(r_std),ALLOCATABLE, DIMENSION (:)                :: soilth_lev    !! Vertical soil axis for thermal scheme (m)
    INTEGER(i_std)                                       :: i, j, ik
    INTEGER(i_std)                                       :: ier
    INTEGER(i_std)                                       :: itau_sechiba
    INTEGER                                              :: old_fileout   !! old Logical Int for std IO output


    IF (printlev >= 2) WRITE(numout,*) 'Start intersurf_initialize_2d'
    CALL ipslnlf_p(new_number=numout,old_number=old_fileout)

    ! Initialize variables in module time
    CALL time_initialize(kjit, date0, xrdt, "END")

    !  Configuration of SSL specific parameters
    CALL control_initialize

    ! Initialize specific write level
    printlev_loc=get_printlev('instersurf')
    
    OFF_LINE_MODE = .TRUE. 
    
    DO ik=1,kjpindex
       
       j = ((kindex(ik)-1)/iim) + 1
       i = (kindex(ik) - (j-1)*iim)
       
       !- Create the internal coordinate table
       !-
       lalo(ik,1) = lat(i,j)
       lalo(ik,2) = lon(i,j)
       !
       !- Store the fraction of the continents only once so that the user
       !- does not change them afterwards.
       !-
       contfrac(ik) = zcontfrac(i,j)
    ENDDO
    CALL gather(contfrac,contfrac_g)
    CALL gather(lalo,lalo_g)
    CALL gather2D_mpi(lon,lon_g)
    CALL gather2D_mpi(lat,lat_g)

    CALL bcast(lalo_g)
    CALL bcast(contfrac_g)
    
    CALL ioipslctrl_restini(kjit, date0, xrdt, rest_id, rest_id_stom, itau_offset, date0_shifted)
    itau_sechiba = kjit + itau_offset
    
    !!- Initialize module for output with XIOS
    !
    ! Get the vertical soil levels for the thermal scheme, to be used in xios_orchidee_init
    ALLOCATE(soilth_lev(ngrnd), stat=ier)
    IF (ier /= 0) THEN
       CALL ipslerr_p(3,'intersurf_main_2d', 'Error in allocation of soilth_lev','','')
    END IF
    IF (hydrol_cwrr) THEN
       soilth_lev(1:ngrnd) = znt(:)
    ELSE
       soilth_lev(1:ngrnd) = thermosoilc_levels()
    END IF

    CALL xios_orchidee_init( MPI_COMM_ORCH,                        &
         date0,   year_end,  month_end,     day_end,  julian_diff, &
         lon,     lat,       soilth_lev )
    
    !- Initialize IOIPSL sechiba output files
    CALL ioipslctrl_history(iim, jjm, lon, lat,  kindex, kjpindex, itau_sechiba, date0_shifted, xrdt, hist_id, &
         hist2_id, hist_id_stom, hist_id_stom_IPCC)
 
    !
    !  Shift the time step to phase the two models
    !
    itau_sechiba = kjit + itau_offset
    
    ! Update the calendar in xios by sending the new time step
    ! Special case : the model is only in initialization phase and the current itau_sechiba is not a real time step. 
    ! Therefor give itau_sechiba+1 to xios to have a correct time axis in output files. 
    CALL xios_orchidee_update_calendar(itau_sechiba+1)

    !
    ! 1. gather input fields from kindex array
    !    Warning : I'm not sure this interface with one dimension array is the good one
    !
    DO ik=1, kjpindex
      
       j = ((kindex(ik)-1)/iim) + 1
       i = (kindex(ik) - (j-1)*iim)
       
       zu(ik)           = u(i,j)
       zv(ik)           = v(i,j)
       zzlev(ik)        = zlev(i,j)
       zqair(ik)        = qair(i,j)
       zprecip_rain(ik) = precip_rain(i,j)*xrdt
       zprecip_snow(ik) = precip_snow(i,j)*xrdt
       zlwdown(ik)      = lwdown(i,j)
       zswnet(ik)       = swnet(i,j)
       zswdown(ik)      = swdown(i,j)
       ztemp_air(ik)    = temp_air(i,j)
       zepot_air(ik)    = epot_air(i,j)
       zccanopy(ik)     = ccanopy(i,j)
       zpetAcoef(ik)    = petAcoef(i,j)
       zpeqAcoef(ik)    = peqAcoef(i,j)
       zpetBcoef(ik)    = petBcoef(i,j)
       zpeqBcoef(ik)    = peqBcoef(i,j)
       zcdrag(ik)       = cdrag(i,j)
       zpb(ik)          = pb(i,j)
       
    ENDDO

    !
    ! 2. save the grid
    !
    CALL histwrite_p(hist_id, 'LandPoints',  itau_sechiba+1, (/ ( REAL(ik), ik=1,kjpindex ) /), kjpindex, kindex)
    CALL histwrite_p(hist_id, 'Areas',  itau_sechiba+1, area, kjpindex, kindex)
    IF ( ok_stomate ) THEN
       CALL histwrite_p(hist_id_stom, 'Areas',  itau_sechiba+1, area, kjpindex, kindex)
       CALL histwrite_p(hist_id_stom_IPCC, 'Areas',  itau_sechiba+1, area, kjpindex, kindex)
    ENDIF
    
    CALL histwrite_p(hist_id, 'Contfrac',  itau_sechiba+1, contfrac, kjpindex, kindex)
    IF ( is_omp_root .AND. hist_id > 0 ) THEN
       ! Always syncronize output after initialization 
       CALL histsync(hist_id)
    END IF
    
    CALL histwrite_p(hist2_id, 'LandPoints',  itau_sechiba+1, (/ ( REAL(ik), ik=1,kjpindex ) /), kjpindex, kindex)
    CALL histwrite_p(hist2_id, 'Areas',  itau_sechiba+1, area, kjpindex, kindex)
    CALL histwrite_p(hist2_id, 'Contfrac',  itau_sechiba+1, contfrac, kjpindex, kindex)
    IF ( is_omp_root .AND. hist2_id > 0 ) THEN
       ! Always syncronize output after initialization 
       CALL histsync(hist2_id)
    ENDIF
    
    !
    ! 3. call sechiba for continental points only
    !
    IF (printlev_loc >= 3) WRITE(numout,*) 'Before call to sechiba_initialize'
    
    CALL sechiba_initialize( &
         itau_sechiba, iim*jjm,      kjpindex,      kindex,      date0_shifted, &
         lalo,         contfrac,     neighbours,    resolution,  zzlev,         &
         zu,           zv,           zqair,         ztemp_air,   ztemp_air,     &
         zpetAcoef,    zpeqAcoef,    zpetBcoef,     zpeqBcoef,                  &
         zprecip_rain, zprecip_snow, zlwdown,       zswnet,      zswdown,       &
         zpb,          rest_id,      hist_id,       hist2_id,                   &
         rest_id_stom, hist_id_stom, hist_id_stom_IPCC,                         &
         zcoastal,     zriver,       ztsol_rad,     zvevapp,     zqsurf,        &
         zz0m,         zz0h,         zalbedo,      zfluxsens,     zfluxlat,     &
         zemis,        znetco2,      zcarblu,      ztemp_sol_new, zcdrag)
    
    IF (printlev_loc >= 3) WRITE(numout,*) 'After call to sechiba_initialize'
    !
    ! 5. scatter output fields
    !
    z0m(:,:)           = undef_sechiba
    coastalflow(:,:)  = undef_sechiba
    riverflow(:,:)    = undef_sechiba
    tsol_rad(:,:)     = undef_sechiba
    vevapp(:,:)       = undef_sechiba
    temp_sol_new(:,:) = undef_sechiba 
    qsurf(:,:)        = undef_sechiba 
    albedo(:,:,:)     = undef_sechiba
    fluxsens(:,:)     = undef_sechiba
    fluxlat(:,:)      = undef_sechiba
    emis(:,:)         = undef_sechiba 
    cdrag(:,:)        = undef_sechiba 
    
    DO ik=1, kjpindex
       j = ((kindex(ik)-1)/iim) + 1
       i = (kindex(ik) - (j-1)*iim)
       
       z0m(i,j)           = zz0m(ik)
       coastalflow(i,j)  = zcoastal(ik)
       riverflow(i,j)    = zriver(ik)
       tsol_rad(i,j)     = ztsol_rad(ik)
       vevapp(i,j)       = zvevapp(ik)
       temp_sol_new(i,j) = ztemp_sol_new(ik)
       qsurf(i,j)        = zqsurf(ik)
       albedo(i,j,1)     = zalbedo(ik,1)
       albedo(i,j,2)     = zalbedo(ik,2)
       fluxsens(i,j)     = zfluxsens(ik)
       fluxlat(i,j)      = zfluxlat(ik)
       emis(i,j)         = zemis(ik)
       cdrag(i,j)        = zcdrag(ik)

    ENDDO

    !
    ! 6. Transform the water fluxes into Kg/m^2s and m^3/s
    !
    DO ik=1, kjpindex
    
       j = ((kindex(ik)-1)/iim) + 1
       i = (kindex(ik) - (j-1)*iim)

       vevapp(i,j) = vevapp(i,j)/xrdt
       coastalflow(i,j) = coastalflow(i,j)/xrdt
       riverflow(i,j) = riverflow(i,j)/xrdt

    ENDDO

    IF (is_root_prc) CALL getin_dump

    lstep_init_intersurf = .FALSE.
    CALL ipslnlf_p(new_number=old_fileout)
    IF (printlev_loc >= 1) WRITE (numout,*) 'End intersurf_initialize_2d'

  END SUBROUTINE intersurf_initialize_2d


!!  =============================================================================================================================
!! SUBROUTINE:    intersurf_main_2d
!!
!>\BRIEF	  Main subroutine to call ORCHIDEE from dim2_driver using variables on a 2d grid.
!!
!! DESCRIPTION:	  This subroutine is the main interface for ORCHIDEE when it is called from the offline driver dim2_driver.
!!                The variables are all on the 2D grid including ocean points. intersurf_initialize_2d should be called before
!!                this subroutine is called. This subroutine is called at each time step.
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE intersurf_main_2d (kjit, iim, jjm, kjpindex, kindex, xrdt, &
       lrestart_read, lrestart_write, lon, lat, zcontfrac, zresolution, date0, &
       zlev, u, v, qair, temp_air, epot_air, ccanopy, &
       cdrag, petAcoef, peqAcoef, petBcoef, peqBcoef, &
       precip_rain, precip_snow, lwdown, swnet, swdown, pb, &
       vevapp, fluxsens, fluxlat, coastalflow, riverflow, &
       tsol_rad, temp_sol_new, qsurf, albedo, emis, z0m, &
       coszang)

    IMPLICIT NONE

    !! 0. Variable and parameter declaration
    !! 0.1 Input variables
    INTEGER(i_std),INTENT (in)                              :: kjit            !! Time step number
    INTEGER(i_std),INTENT (in)                              :: iim, jjm        !! Dimension of input fields
    INTEGER(i_std),INTENT (in)                              :: kjpindex        !! Number of continental points
    REAL(r_std),INTENT (in)                                 :: xrdt            !! Time step in seconds
    LOGICAL, INTENT (in)                                    :: lrestart_read   !! Logical for _restart_ file to read
    LOGICAL, INTENT (in)                                    :: lrestart_write  !! Logical for _restart_ file to write'
    REAL(r_std), INTENT (in)                                :: date0           !! Date at which kjit = 0
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)        :: kindex          !! Index for continental points
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: u             !! Lowest level wind speed
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: v             !! Lowest level wind speed 
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: zlev          !! Height of first layer
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: qair          !! Lowest level specific humidity
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: precip_rain   !! Rain precipitation
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: precip_snow   !! Snow precipitation
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: lwdown        !! Down-welling long-wave flux 
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: swnet         !! Net surface short-wave flux
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: swdown        !! Downwelling surface short-wave flux
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: coszang       !! Cosine of the solar zenith angle (unitless)
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: temp_air      !! Air temperature in Kelvin
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: epot_air      !! Air potential energy
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: ccanopy       !! CO2 concentration in the canopy
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: petAcoef      !! Coeficients A from the PBL resolution
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: peqAcoef      !! One for T and another for q
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: petBcoef      !! Coeficients B from the PBL resolution
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: peqBcoef      !! One for T and another for q
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: pb            !! Surface pressure
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: lon, lat      !! Geographical coordinates
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: zcontfrac     !! Fraction of continent in the grid
    REAL(r_std),DIMENSION (iim,jjm,2), INTENT(in)           :: zresolution   !! resolution in x and y dimensions

    !! 0.2 Output variables
    REAL(r_std),DIMENSION (iim,jjm), INTENT(out)            :: z0m            !! Surface roughness for momemtum
    REAL(r_std),DIMENSION (iim,jjm), INTENT(out)            :: coastalflow   !! Diffuse flow of water into the ocean (m^3/s)
    REAL(r_std),DIMENSION (iim,jjm), INTENT(out)            :: riverflow     !! Largest rivers flowing into the ocean (m^3/s)
    REAL(r_std),DIMENSION (iim,jjm), INTENT(out)            :: tsol_rad      !! Radiative surface temperature
    REAL(r_std),DIMENSION (iim,jjm), INTENT(out)            :: vevapp        !! Total of evaporation
    REAL(r_std),DIMENSION (iim,jjm), INTENT(out)            :: temp_sol_new  !! New soil temperature
    REAL(r_std),DIMENSION (iim,jjm), INTENT(out)            :: qsurf         !! Surface specific humidity
    REAL(r_std),DIMENSION (iim,jjm,2), INTENT(out)          :: albedo        !! Albedo
    REAL(r_std),DIMENSION (iim,jjm), INTENT(out)            :: fluxsens      !! Sensible chaleur flux
    REAL(r_std),DIMENSION (iim,jjm), INTENT(out)            :: fluxlat       !! Latent chaleur flux
    REAL(r_std),DIMENSION (iim,jjm), INTENT(out)            :: emis          !! Emissivity

    !! 0.3 Modified variables
    REAL(r_std),DIMENSION (iim,jjm), INTENT(inout)          :: cdrag         !! Cdrag

    !! 0.4 Local variables
    REAL(r_std),DIMENSION (kjpindex)                      :: zu            !! Work array to keep u
    REAL(r_std),DIMENSION (kjpindex)                      :: zv            !! Work array to keep v
    REAL(r_std),DIMENSION (kjpindex)                      :: zzlev         !! Work array to keep zlev
    REAL(r_std),DIMENSION (kjpindex)                      :: zqair         !! Work array to keep qair
    REAL(r_std),DIMENSION (kjpindex)                      :: zprecip_rain  !! Work array to keep precip_rain
    REAL(r_std),DIMENSION (kjpindex)                      :: zprecip_snow  !! Work array to keep precip_snow
    REAL(r_std),DIMENSION (kjpindex)                      :: zlwdown       !! Work array to keep lwdown
    REAL(r_std),DIMENSION (kjpindex)                      :: zswnet        !! Work array to keep swnet
    REAL(r_std),DIMENSION (kjpindex)                      :: zswdown       !! Work array to keep swdown
    REAL(r_std),DIMENSION (kjpindex)                      :: zcoszang      !! Work array to keep coszang
    REAL(r_std),DIMENSION (kjpindex)                      :: ztemp_air     !! Work array to keep temp_air
    REAL(r_std),DIMENSION (kjpindex)                      :: zepot_air     !! Work array to keep epot_air
    REAL(r_std),DIMENSION (kjpindex)                      :: zccanopy      !! Work array to keep ccanopy
    REAL(r_std),DIMENSION (kjpindex)                      :: zpetAcoef     !! Work array to keep petAcoef
    REAL(r_std),DIMENSION (kjpindex)                      :: zpeqAcoef     !! Work array to keep peqAcoef
    REAL(r_std),DIMENSION (kjpindex)                      :: zpetBcoef     !! Work array to keep petBcoef
    REAL(r_std),DIMENSION (kjpindex)                      :: zpeqBcoef     !! Work array to keep peqVcoef
    REAL(r_std),DIMENSION (kjpindex)                      :: zcdrag        !! Work array to keep cdrag
    REAL(r_std),DIMENSION (kjpindex)                      :: zpb           !! Work array to keep pb
    REAL(r_std),DIMENSION (kjpindex)                      :: zz0m, zz0h    !! Work array to keep zz0m, zz0h
    REAL(r_std),DIMENSION (kjpindex,nvm)                  :: zveget        !! Work array to keep veget
    REAL(r_std),DIMENSION (kjpindex,nvm)                  :: zlai          !! Work array to keep lai
    REAL(r_std),DIMENSION (kjpindex,nvm)                  :: zheight       !! Work array to keep height
    REAL(r_std),DIMENSION (kjpindex)                      :: zcoastal      !! Work array to keep coastalflow (m^3/dt)
    REAL(r_std),DIMENSION (kjpindex)                      :: zriver        !! Work array to keep riverflow (m^3/dt)
    REAL(r_std),DIMENSION (kjpindex)                      :: znetco2       !! Work array to keep netco2flux
    REAL(r_std),DIMENSION (kjpindex)                      :: zcarblu       !! Work array to keep fco2_land_use
    REAL(r_std),DIMENSION (kjpindex)                      :: ztsol_rad     !! Work array to keep tsol_rad
    REAL(r_std),DIMENSION (kjpindex)                      :: zvevapp       !! Work array to keep vevapp
    REAL(r_std),DIMENSION (kjpindex)                      :: ztemp_sol_new !! Work array to keep temp_sol_new
    REAL(r_std),DIMENSION (kjpindex)                      :: zqsurf        !! Work array to keep qsurf
    REAL(r_std),DIMENSION (kjpindex,2)                    :: zalbedo       !! Work array to keep albedo
    REAL(r_std),DIMENSION (kjpindex)                      :: zfluxsens     !! Work array to keep fluxsens
    REAL(r_std),DIMENSION (kjpindex)                      :: zfluxlat      !! Work array to keep fluxlat
    REAL(r_std),DIMENSION (kjpindex)                      :: zemis         !! Work array to keep emis
    REAL(r_std),ALLOCATABLE, DIMENSION (:)                :: soilth_lev    !! Vertical soil axis for thermal scheme (m)
    INTEGER(i_std)                                        :: i, j, ik
    INTEGER(i_std)                                        :: ier
    INTEGER(i_std)                                        :: itau_sechiba
    INTEGER                                               :: old_fileout   !! old Logical Int for std IO output

    IF (printlev_loc >= 3) WRITE(numout,*) 'Start intersurf_main_2d'
    CALL ipslnlf_p(new_number=numout,old_number=old_fileout)

    !
    !  Shift the time step to phase the two models
    !
    itau_sechiba = kjit + itau_offset
    !

    ! Update the calendar in xios by sending the new time step
    CALL xios_orchidee_update_calendar(itau_sechiba)

    ! Update the calendar and all time variables in module time
    CALL time_nextstep(itau_sechiba)
    !
    ! 1. gather input fields from kindex array
    !    Warning : I'm not sure this interface with one dimension array is the good one
    !
    DO ik=1, kjpindex
      
       j = ((kindex(ik)-1)/iim) + 1
       i = (kindex(ik) - (j-1)*iim)
       
       zu(ik)           = u(i,j)
       zv(ik)           = v(i,j)
       zzlev(ik)        = zlev(i,j)
       zqair(ik)        = qair(i,j)
       zprecip_rain(ik) = precip_rain(i,j)*xrdt
       zprecip_snow(ik) = precip_snow(i,j)*xrdt
       zlwdown(ik)      = lwdown(i,j)
       zswnet(ik)       = swnet(i,j)
       zswdown(ik)      = swdown(i,j)
       zcoszang(ik)     = coszang(i,j)
       ztemp_air(ik)    = temp_air(i,j)
       zepot_air(ik)    = epot_air(i,j)
       zccanopy(ik)     = ccanopy(i,j)
       zpetAcoef(ik)    = petAcoef(i,j)
       zpeqAcoef(ik)    = peqAcoef(i,j)
       zpetBcoef(ik)    = petBcoef(i,j)
       zpeqBcoef(ik)    = peqBcoef(i,j)
       zcdrag(ik)       = cdrag(i,j)
       zpb(ik)          = pb(i,j)
       
    ENDDO

    !
    ! 3. call sechiba for continental points only
    !
    IF (printlev_loc >= 3) WRITE(numout,*) 'Before call to sechiba_main from intersurf_main_2d'

    CALL sechiba_main (itau_sechiba, iim*jjm, kjpindex, kindex, date0_shifted, &
         lrestart_read, lrestart_write, &
         lalo, contfrac, neighbours, resolution, &
         zzlev, zu, zv, zqair, zqair, ztemp_air, ztemp_air, zepot_air, zccanopy, &
         zcdrag, zpetAcoef, zpeqAcoef, zpetBcoef, zpeqBcoef, &
         zprecip_rain ,zprecip_snow,  zlwdown, zswnet, zswdown, zcoszang, zpb, &
         zvevapp, zfluxsens, zfluxlat, zcoastal, zriver, znetco2, zcarblu, &
         ztsol_rad, ztemp_sol_new, zqsurf, zalbedo, zemis, zz0m, zz0h,&
         zveget, zlai, zheight, &
         rest_id, hist_id, hist2_id, rest_id_stom, hist_id_stom, hist_id_stom_IPCC ) 
    
    IF (printlev_loc >= 3) WRITE(numout,*) 'After call to sechiba_main'

    !
    ! 5. scatter output fields
    !
    z0m(:,:)          = undef_sechiba
    coastalflow(:,:)  = undef_sechiba
    riverflow(:,:)    = undef_sechiba
    tsol_rad(:,:)     = undef_sechiba
    vevapp(:,:)       = undef_sechiba
    temp_sol_new(:,:) = undef_sechiba 
    qsurf(:,:)        = undef_sechiba 
    albedo(:,:,:)     = undef_sechiba
    fluxsens(:,:)     = undef_sechiba
    fluxlat(:,:)      = undef_sechiba
    emis(:,:)         = undef_sechiba 
    cdrag(:,:)        = undef_sechiba 
    !
    DO ik=1, kjpindex
      
    
       j = ((kindex(ik)-1)/iim) + 1
       i = (kindex(ik) - (j-1)*iim)

       z0m(i,j)           = zz0m(ik)
       coastalflow(i,j)  = zcoastal(ik)
       riverflow(i,j)    = zriver(ik)
       tsol_rad(i,j)     = ztsol_rad(ik)
       vevapp(i,j)       = zvevapp(ik)
       temp_sol_new(i,j) = ztemp_sol_new(ik)
       qsurf(i,j)        = zqsurf(ik)
       albedo(i,j,1)     = zalbedo(ik,1)
       albedo(i,j,2)     = zalbedo(ik,2)
       fluxsens(i,j)     = zfluxsens(ik)
       fluxlat(i,j)      = zfluxlat(ik)
       emis(i,j)         = zemis(ik)
       cdrag(i,j)        = zcdrag(ik)

    ENDDO

    CALL xios_orchidee_send_field("LandPoints" ,(/ ( REAL(ik), ik=1,kjpindex ) /))
    CALL xios_orchidee_send_field("areas", area)
    CALL xios_orchidee_send_field("contfrac",contfrac)
    CALL xios_orchidee_send_field("temp_air",ztemp_air)
    CALL xios_orchidee_send_field("qair",zqair)
    CALL xios_orchidee_send_field("swnet",zswnet)
    CALL xios_orchidee_send_field("swdown",zswdown)
    CALL xios_orchidee_send_field("pb",zpb)
    
    IF ( .NOT. almaoutput ) THEN
       !
       !  scattered during the writing
       ! 
       CALL histwrite_p (hist_id, 'evap',     itau_sechiba, zvevapp, kjpindex, kindex)
       CALL histwrite_p (hist_id, 'coastalflow',itau_sechiba, zcoastal, kjpindex, kindex)
       CALL histwrite_p (hist_id, 'riverflow',itau_sechiba, zriver, kjpindex, kindex)
       CALL histwrite_p (hist_id, 'temp_sol', itau_sechiba, ztemp_sol_NEW, kjpindex, kindex)
       CALL histwrite_p (hist_id, 'tsol_max', itau_sechiba, ztemp_sol_NEW, kjpindex, kindex)
       CALL histwrite_p (hist_id, 'tsol_min', itau_sechiba, ztemp_sol_NEW, kjpindex, kindex)
       CALL histwrite_p (hist_id, 'fluxsens', itau_sechiba, zfluxsens, kjpindex, kindex)
       CALL histwrite_p (hist_id, 'fluxlat',  itau_sechiba, zfluxlat, kjpindex, kindex)
       CALL histwrite_p (hist_id, 'swnet',    itau_sechiba, zswnet, kjpindex, kindex)
       CALL histwrite_p (hist_id, 'swdown',   itau_sechiba, zswdown, kjpindex, kindex)
       CALL histwrite_p (hist_id, 'alb_vis',  itau_sechiba, zalbedo(:,1), kjpindex, kindex)
       CALL histwrite_p (hist_id, 'alb_nir',  itau_sechiba, zalbedo(:,2), kjpindex, kindex)
       CALL histwrite_p (hist_id, 'tair',     itau_sechiba, ztemp_air, kjpindex, kindex)
       CALL histwrite_p (hist_id, 'qair',     itau_sechiba, zqair, kjpindex, kindex)
       CALL histwrite_p (hist_id, 'q2m',     itau_sechiba, zqair, kjpindex, kindex)
       CALL histwrite_p (hist_id, 't2m',     itau_sechiba, ztemp_air, kjpindex, kindex)
       CALL histwrite_p (hist2_id, 'evap',     itau_sechiba, zvevapp, kjpindex, kindex)
       CALL histwrite_p (hist2_id, 'coastalflow',itau_sechiba, zcoastal, kjpindex, kindex)
       CALL histwrite_p (hist2_id, 'riverflow',itau_sechiba, zriver, kjpindex, kindex)
       CALL histwrite_p (hist2_id, 'temp_sol', itau_sechiba, ztemp_sol_NEW, kjpindex, kindex)
       CALL histwrite_p (hist2_id, 'tsol_max', itau_sechiba, ztemp_sol_NEW, kjpindex, kindex)
       CALL histwrite_p (hist2_id, 'tsol_min', itau_sechiba, ztemp_sol_NEW, kjpindex, kindex)
       CALL histwrite_p (hist2_id, 'fluxsens', itau_sechiba, zfluxsens, kjpindex, kindex)
       CALL histwrite_p (hist2_id, 'fluxlat',  itau_sechiba, zfluxlat, kjpindex, kindex)
       CALL histwrite_p (hist2_id, 'swnet',    itau_sechiba, zswnet, kjpindex, kindex)
       CALL histwrite_p (hist2_id, 'swdown',   itau_sechiba, zswdown, kjpindex, kindex)
       CALL histwrite_p (hist2_id, 'alb_vis',  itau_sechiba, zalbedo(:,1), kjpindex, kindex)
       CALL histwrite_p (hist2_id, 'alb_nir',  itau_sechiba, zalbedo(:,2), kjpindex, kindex)
       CALL histwrite_p (hist2_id, 'tair',     itau_sechiba, ztemp_air, kjpindex, kindex)
       CALL histwrite_p (hist2_id, 'qair',     itau_sechiba, zqair, kjpindex, kindex)
       CALL histwrite_p (hist2_id, 'q2m',     itau_sechiba, zqair, kjpindex, kindex)
       CALL histwrite_p (hist2_id, 't2m',     itau_sechiba, ztemp_air, kjpindex, kindex)
    ELSE
       CALL histwrite_p (hist_id, 'Evap', itau_sechiba, zvevapp, kjpindex, kindex)
       CALL histwrite_p (hist_id, 'SWnet',    itau_sechiba, zswnet, kjpindex, kindex)
       CALL histwrite_p (hist_id, 'Qh', itau_sechiba, zfluxsens, kjpindex, kindex)
       CALL histwrite_p (hist_id, 'Qle',  itau_sechiba, zfluxlat, kjpindex, kindex)
       CALL histwrite_p (hist_id, 'AvgSurfT', itau_sechiba, ztemp_sol_NEW, kjpindex, kindex)
       CALL histwrite_p (hist_id, 'RadT', itau_sechiba, ztemp_sol_NEW, kjpindex, kindex)
       CALL histwrite_p (hist_id, 'Tair', itau_sechiba, ztemp_air, kjpindex, kindex)
       CALL histwrite_p (hist_id, 'Qair', itau_sechiba, zqair, kjpindex, kindex)
       CALL histwrite_p (hist2_id, 'Evap', itau_sechiba, zvevapp, kjpindex, kindex)
       CALL histwrite_p (hist2_id, 'SWnet',    itau_sechiba, zswnet, kjpindex, kindex)
       CALL histwrite_p (hist2_id, 'Qh', itau_sechiba, zfluxsens, kjpindex, kindex)
       CALL histwrite_p (hist2_id, 'Qle',  itau_sechiba, zfluxlat, kjpindex, kindex)
       CALL histwrite_p (hist2_id, 'AvgSurfT', itau_sechiba, ztemp_sol_NEW, kjpindex, kindex)
       CALL histwrite_p (hist2_id, 'RadT', itau_sechiba, ztemp_sol_NEW, kjpindex, kindex)
    ENDIF
    !
    IF ( is_omp_root ) THEN
       IF ( (dw .EQ. xrdt) .AND. hist_id > 0 ) THEN
          ! Syncronize output but only if flag ok_histsync is set to true
          IF (ok_histsync) CALL histsync(hist_id)
       ENDIF
    END IF

    !
    ! 6. Transform the water fluxes into Kg/m^2s and m^3/s
    !
    DO ik=1, kjpindex
       
       j = ((kindex(ik)-1)/iim) + 1
       i = (kindex(ik) - (j-1)*iim)
       
       vevapp(i,j) = vevapp(i,j)/xrdt
       coastalflow(i,j) = coastalflow(i,j)/xrdt
       riverflow(i,j) = riverflow(i,j)/xrdt
       
    ENDDO
    
    CALL ipslnlf_p(new_number=old_fileout)
    IF (printlev_loc >= 3) WRITE (numout,*) 'End intersurf_main_2d'

  END SUBROUTINE intersurf_main_2d


!!  =============================================================================================================================
!! SUBROUTINE:    init_intersurf
!!
!>\BRIEF	  Initialize grid information
!!
!! DESCRIPTION:	  This subroutine is called from LMDZ before first call to intersurf_main_gathered or 
!!                intersurf_initialize_gathered
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE init_intersurf(nbp_l_lon,nbp_l_lat,kjpindex,kindex,orch_offset,orch_omp_size,orch_omp_rank,COMM)

    USE mod_orchidee_para
    USE timer
    IMPLICIT NONE

    INTEGER,INTENT(IN)  :: nbp_l_lon
    INTEGER,INTENT(IN)  :: nbp_l_lat
    INTEGER,INTENT(IN)  :: kjpindex
    INTEGER,INTENT(IN)  :: kindex(:)
    INTEGER,INTENT(IN)  :: orch_offset
    INTEGER,INTENT(IN)  :: COMM
    INTEGER,INTENT(IN)  :: orch_omp_size
    INTEGER,INTENT(IN)  :: orch_omp_rank

    INTEGER,DIMENSION(kjpindex)  :: kindex_offset

    IF (printlev >= 1) WRITE(*,*) 'Start ORCHIDEE'

    IF (orch_omp_rank==0) THEN
      CALL Init_timer
      CALL start_timer(timer_mpi)
      CALL grid_set_glo(nbp_l_lon,nbp_l_lat)
    ENDIF
    CALL barrier2_omp()    
    CALL init_orchidee_data_para(kjpindex,kindex,orch_offset,orch_omp_size,orch_omp_rank,COMM)
    CALL Set_stdout_file('out_orchidee')

    IF (printlev >= 1) WRITE(numout,*) 'Start ORCHIDEE intitalization phase'
    
    IF (is_omp_root) CALL grid_allocate_glo(4)
    CALL barrier2_omp()
    CALL init_ioipsl_para
          
    kindex_offset(:)=kindex(:)+offset
    CALL gather(kindex_offset,index_g)
    CALL bcast(index_g)  

    IF (printlev_loc >= 2) THEN
       WRITE(numout,*) "kjpindex = ",kjpindex
       WRITE(numout,*) "offset for OMP = ",offset_omp
       WRITE(numout,*) "Index num local for continental points = ",kindex
       WRITE(numout,*) "Index num global for continental points = ",kindex_offset
       IF (is_omp_root) THEN
          WRITE(numout,*) "ROOT OMP, Index global MPI : ",kindex_mpi(:)
       ENDIF
       IF (is_root_prc) THEN
          WRITE(numout,*) "ROOT global, Index global : ",index_g(:)
       ENDIF
    END IF

    ! Allocation of grid variables
    CALL grid_init ( kjpindex, 4, "RegLonLat", "2DGrid" )


  END SUBROUTINE init_intersurf

!!  =============================================================================================================================
!! SUBROUTINE:    intersurf_initialize_gathered
!!
!>\BRIEF	  Initialization and call to sechiba_initialize
!!
!! DESCRIPTION:	  Initialization of module variables, read options from parameter file, initialize output files and call to 
!!                sechiba_initialize.
!!
!!                This subroutine can be called directly from GCM(LMDZ). If it is not called before the first call to 
!!                intersurf_main_gathered, then it will be done from there. This possibility is done to keep backward 
!!                compatibility with LMDZ. 
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE intersurf_initialize_gathered (kjit, iim_glo, jjm_glo, kjpindex, kindex, xrdt, &
       lrestart_read, lrestart_write, latlon, zcontfrac, zneighbours, zresolution, date0, &
       zlev,  u, v, qair, temp_air, epot_air, &
       cdrag, petAcoef, peqAcoef, petBcoef, peqBcoef, &
       precip_rain, precip_snow, lwdown, swnet, swdown, pb, &
       vevapp, fluxsens, fluxlat, coastalflow_cpl, riverflow_cpl, &
       tsol_rad, temp_sol_new, qsurf, albedo, emis, z0m, lon_scat_g, lat_scat_g, &
       q2m, t2m, z0h, nvm_out, &
       field_out_names, fields_out, field_in_names, fields_in)

    USE mod_orchidee_para
    IMPLICIT NONE

    !! 0. Variable and parameter declaration
    !! 0.1 Input 
    INTEGER(i_std),INTENT (in)                             :: kjit           !! Time step number
    INTEGER(i_std),INTENT (in)                             :: iim_glo        !! Dimension of global fields
    INTEGER(i_std),INTENT (in)                             :: jjm_glo        !! Dimension of global fields
    INTEGER(i_std),INTENT (in)                             :: kjpindex       !! Number of continental points
    REAL(r_std),INTENT (in)                                :: xrdt           !! Time step in seconds
    LOGICAL, INTENT (in)                                   :: lrestart_read  !! Logical for _restart_ file to read
    LOGICAL, INTENT (in)                                   :: lrestart_write !! Logical for _restart_ file to write'
    REAL(r_std), INTENT (in)                               :: date0          !! Date at which kjit = 0
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)       :: kindex         !! Index for continental points
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)           :: u              !! Lowest level wind speed
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)           :: v              !! Lowest level wind speed 
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)           :: zlev           !! Height of first layer
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)           :: qair           !! Lowest level specific humidity
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)           :: precip_rain    !! Rain precipitation
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)           :: precip_snow    !! Snow precipitation
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)           :: lwdown         !! Down-welling long-wave flux 
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)           :: swnet          !! Net surface short-wave flux
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)           :: swdown         !! Downwelling surface short-wave flux
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)           :: temp_air       !! Air temperature in Kelvin
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)           :: epot_air       !! Air potential energy
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)           :: petAcoef       !! Coeficients A from the PBL resolution
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)           :: peqAcoef       !! One for T and another for q
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)           :: petBcoef       !! Coeficients B from the PBL resolution
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)           :: peqBcoef       !! One for T and another for q
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)           :: pb             !! Surface pressure
    REAL(r_std),DIMENSION (kjpindex,2), INTENT(in)         :: latlon         !! Geographical coordinates
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)           :: zcontfrac      !! Fraction of continent
    INTEGER(i_std),DIMENSION (kjpindex,NbNeighb), INTENT(in):: zneighbours    !! neighbours
    REAL(r_std),DIMENSION (kjpindex,2), INTENT(in)         :: zresolution    !! size of the grid box
    REAL(r_std),DIMENSION (iim_glo,jjm_glo), INTENT(IN)    :: lon_scat_g     !! Longitudes on the global 2D grid including ocean 
    REAL(r_std),DIMENSION (iim_glo,jjm_glo), INTENT(IN)    :: lat_scat_g     !! Latitudes on the global 2D grid including ocean 

    !! 0.2 Output variables
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)          :: z0m            !! Surface roughness (momentum)
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)          :: coastalflow_cpl!! Diffuse flow of water into the ocean (m^3/s)
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)          :: riverflow_cpl  !! Largest rivers flowing into the ocean (m^3/s)
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)          :: tsol_rad       !! Radiative surface temperature
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)          :: vevapp         !! Total of evaporation
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)          :: temp_sol_new   !! New soil temperature
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)          :: qsurf          !! Surface specific humidity
    REAL(r_std),DIMENSION (kjpindex,2), INTENT(out)        :: albedo         !! Albedo
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)          :: fluxsens       !! Sensible chaleur flux
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)          :: fluxlat        !! Latent chaleur flux
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)          :: emis           !! Emissivity

    !! 0.3 Modified variables
    REAL(r_std),DIMENSION (kjpindex), INTENT(inout)        :: cdrag          !! Cdrag

    !! 0.4 Optional input and output variables
    REAL(r_std),DIMENSION (kjpindex), INTENT(in), OPTIONAL :: q2m            !! Surface specific humidity
    REAL(r_std),DIMENSION (kjpindex), INTENT(in), OPTIONAL :: t2m            !! Surface air temperature
    REAL(r_std),DIMENSION (kjpindex), INTENT(out),OPTIONAL :: z0h            !! Surface roughness (heat)
    INTEGER(i_std), INTENT(out), OPTIONAL                  :: nvm_out        !! Number of vegetation types, PFTs
    CHARACTER(LEN=*),DIMENSION(:), OPTIONAL, INTENT(IN)    :: field_in_names !! Names for deposit variables to be transported
                                                                             !! from chemistry model by GCM to ORCHIDEE
    REAL(r_std),DIMENSION(:,:), OPTIONAL, INTENT(IN)       :: fields_in      !! Fields for deposit variables to be transported 
                                                                             !! from chemistry model by GCM to ORCHIDEE
    CHARACTER(LEN=*),DIMENSION(:), OPTIONAL, INTENT(IN)    :: field_out_names!! Names for emission variables to be transported 
                                                                             !! to chemistry model by GCM from ORCHIDEE
    REAL(r_std),DIMENSION(:,:), OPTIONAL, INTENT(OUT)      :: fields_out     !! Fields for emission variables to be transported
                                                                             !! to chemistry model by GCM from ORCHIDEE

    !! 0.5 Local variables
    REAL(r_std),DIMENSION (kjpindex)                      :: zprecip_rain  !! Work array to keep precip_rain
    REAL(r_std),DIMENSION (kjpindex)                      :: zprecip_snow  !! Work array to keep precip_snow
    REAL(r_std),DIMENSION (kjpindex)                      :: zz0m          !! Work array to keep zz0m
    REAL(r_std),DIMENSION (kjpindex)                      :: zz0h          !! Work array to keep zz0h
    REAL(r_std),DIMENSION (kjpindex)                      :: zcdrag        !! Work array for surface drag
    REAL(r_std),DIMENSION (kjpindex)                      :: zcoastal      !! Work array to keep coastal flow
    REAL(r_std),DIMENSION (kjpindex)                      :: zriver        !! Work array to keep river out flow
    REAL(r_std),DIMENSION (kjpindex)                      :: znetco2       !! Work array to keep netco2flux
    REAL(r_std),DIMENSION (kjpindex)                      :: zcarblu       !! Work array to keep fco2_land_use
    REAL(r_std),DIMENSION (kjpindex)                      :: ztsol_rad     !! Work array to keep tsol_rad
    REAL(r_std),DIMENSION (kjpindex)                      :: zvevapp       !! Work array to keep vevapp
    REAL(r_std),DIMENSION (kjpindex)                      :: ztemp_sol_new !! Work array to keep temp_sol_new
    REAL(r_std),DIMENSION (kjpindex)                      :: zqsurf        !! Work array to keep qsurf
    REAL(r_std),DIMENSION (kjpindex,2)                    :: zalbedo       !! Work array to keep albedo
    REAL(r_std),DIMENSION (kjpindex)                      :: zfluxsens     !! Work array to keep fluxsens
    REAL(r_std),DIMENSION (kjpindex)                      :: zfluxlat      !! Work array to keep fluxlat
    REAL(r_std),DIMENSION (kjpindex)                      :: zemis         !! Work array to keep emis
    REAL(r_std),ALLOCATABLE, DIMENSION (:)                :: soilth_lev    !! Vertical soil axis for thermal scheme (m)
    REAL(r_std),DIMENSION (kjpindex)                      :: q2m_loc       !! Work array for q2m or qair
    REAL(r_std),DIMENSION (kjpindex)                      :: t2m_loc       !! Work array for t2m or temp_air
    INTEGER(i_std)                                        :: i, j, ik
    INTEGER(i_std)                                        :: ier
    INTEGER(i_std)                                        :: itau_sechiba
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)              :: tmp_lon       !! Longitudes for local MPI process. Only available on master OMP.
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)              :: tmp_lat       !! Latitudes for local MPI process. Only available on master OMP.
    INTEGER                                               :: old_fileout   !! old Logical Int for std IO output
    REAL,ALLOCATABLE,DIMENSION(:,:)                       :: lalo_mpi
    REAL(r_std),DIMENSION (kjpindex)                      :: landpoints    !! Land point vector
    REAL(r_std)                                           :: tau_outflow   !! Number of days for temoral filter on river- and coastalflow
    
    ! Initialize specific write level
    printlev_loc=get_printlev('instersurf')

    IF (printlev_loc >= 1) WRITE(numout,*) 'Entering intersurf_initialize_gathered'
    IF (printlev_loc >= 2) WRITE(numout,*) 'using printlev_loc for intersurf:', printlev_loc

    OFF_LINE_MODE = .FALSE. 
    CALL ipslnlf_p(new_number=numout, old_number=old_fileout)

    ! Initialize variables in module time
    CALL time_initialize( kjit, date0, xrdt, "END" )

    !  Configuration of SSL specific parameters
    CALL control_initialize

    ! Set the variable ok_q2m_t2m=true if q2m and t2m are present in the call from the gcm. 
    ! If one of the variables are not present in this subroutine, set ok_q2m_t2m=.FALSE.
    ! Otherwise do not change the current value. Note that the current value for ok_q2m_t2m comes 
    ! either from default initialization (true) or from intersurf_main_gathered.
    IF (.NOT. PRESENT(q2m) .OR. .NOT. PRESENT(t2m)) THEN
       ok_q2m_t2m=.FALSE.
    END IF
    
    IF (ok_q2m_t2m) THEN
       t2m_loc=t2m
       q2m_loc=q2m
    ELSE
       t2m_loc=temp_air
       q2m_loc=qair
    END IF
    
    
    !  Create the internal coordinate table
    !
    lalo(:,:) = latlon(:,:)
    CALL gather(lalo,lalo_g)
    CALL bcast(lalo_g)
    !
    !-
    !- Store variable to help describe the grid
    !- once the points are gathered.
    !-
    neighbours(:,:) = zneighbours(:,:)
    CALL gather(neighbours,neighbours_g)
    CALL bcast(neighbours_g)
    !
    resolution(:,:) = zresolution(:,:)
    CALL gather(resolution,resolution_g)
    CALL bcast(resolution_g)
    !
    area(:) = resolution(:,1)*resolution(:,2)
    CALL gather(area,area_g)
    CALL bcast(area_g)
    !
    !- Store the fraction of the continents only once so that the user
    !- does not change them afterwards.
    !
    contfrac(:) = zcontfrac(:)
    CALL gather(contfrac,contfrac_g)
    CALL bcast(contfrac_g)
    !
    !
    !  Create the internal coordinate table
    !
    IF ( (.NOT.ALLOCATED(tmp_lon))) THEN
       ALLOCATE(tmp_lon(iim_g,jj_nb))
    ENDIF
    IF ( (.NOT.ALLOCATED(tmp_lat))) THEN
       ALLOCATE(tmp_lat(iim_g,jj_nb))
    ENDIF


    IF (is_omp_root) THEN
       ! Extract from gloabl variables the longitudes and latitudes 
       ! for the local MPI process using the indices for the latitude bands(jj_begin, jj_end).
       tmp_lon(:,:)=lon_scat_g(:,jj_begin:jj_end)
       tmp_lat(:,:)=lat_scat_g(:,jj_begin:jj_end)

       ! Save the global variables only on mpi root, to be used in other modules
       IF (is_mpi_root) THEN
          lon_g(:,:) = lon_scat_g(:,:)
          lat_g(:,:) = lat_scat_g(:,:)
       ENDIF
    ENDIF
    
    !Config Key  = FORCE_CO2_VEG
    !Config Desc = Flag to force the value of atmospheric CO2 for vegetation.
    !Config If   = Only in coupled mode
    !Config Def  = FALSE
    !Config Help = If this flag is set to true, the ATM_CO2 parameter is used
    !Config        to prescribe the atmospheric CO2.
    !Config        This Flag is only use in couple mode.
    !Config Units = [FLAG]
    fatmco2=.FALSE.
    CALL getin_p('FORCE_CO2_VEG',fatmco2)
    !
    ! Next flag is only use in couple mode with a gcm in intersurf.
    ! In forced mode, it has already been read and set in driver.
    IF ( fatmco2 ) THEN
       atmco2=350.
       CALL getin_p('ATM_CO2',atmco2)
       WRITE(numout,*) 'atmco2 ',atmco2
    ENDIF
    
    CALL ioipslctrl_restini(kjit, date0, xrdt, rest_id, rest_id_stom, itau_offset, date0_shifted)
    itau_sechiba = kjit + itau_offset
    
    !!- Initialize module for output with XIOS
    !
    ! Get the vertical soil levels for the thermal scheme, to be used in xios_orchidee_init
    ALLOCATE(soilth_lev(ngrnd), stat=ier)
    IF (ier /= 0) THEN
       CALL ipslerr_p(3,'intersurf_main_gathered', 'Error in allocation of soilth_lev','','')
    END IF
    IF (hydrol_cwrr) THEN
       soilth_lev(1:ngrnd) = znt(:)
    ELSE
       soilth_lev(1:ngrnd) = thermosoilc_levels()
    END IF

    CALL xios_orchidee_init( MPI_COMM_ORCH,                   &
         date0,    year_end, month_end,     day_end, julian_diff, &
         tmp_lon,  tmp_lat,  soilth_lev )
    
    !- Initialize IOIPSL sechiba output files
    CALL ioipslctrl_history(iim_g, jj_nb, tmp_lon, tmp_lat,  kindex, kjpindex, itau_sechiba, &
         date0_shifted, xrdt, hist_id, hist2_id, hist_id_stom, hist_id_stom_IPCC)
    
    CALL bcast_omp(hist_id)
    CALL bcast_omp(hist2_id)
    CALL bcast_omp(hist_id_stom)
    CALL bcast_omp(hist_id_stom_IPCC)
    
    ! Count number of extra output fields to the GCM if it is not already done. 
    IF (nb_fields_out == -1) THEN
       ! nb_fields_out is not yet calculated. Do it now.  
       ! This means that the call is done directly from GCM.
       IF (PRESENT(field_out_names)) THEN
          nb_fields_out=SIZE(field_out_names)
       ELSE
          nb_fields_out=0
       ENDIF
    END IF

    ! Count number of extra input fields to the GCM if it is not already done. 
    IF (nb_fields_in == -1) THEN
       ! nb_fields_in is not yet calculated. Do it now.  
       ! This means that the call is done directly from GCM.
       IF (PRESENT(field_in_names)) THEN
          nb_fields_in=SIZE(field_in_names)
       ELSE
          nb_fields_in=0
       ENDIF
    END IF


    !
    !! Change to be in the orchidee context for XIOS
    !
    CALL xios_orchidee_change_context("orchidee")
    
    !
    !  Shift the time step to phase the two models
    !
    itau_sechiba = kjit + itau_offset
    
    ! Update the calendar in xios by sending the new time step
    ! Special case : the model is only in initialization phase and the current itau_sechiba is not a real time step. 
    ! Therefor give itau_sechiba+1 to xios to have a correct time axis in output files. 
    CALL xios_orchidee_update_calendar(itau_sechiba+1)
    
    !
    ! 1. Just change the units of some input fields
    !
    DO ik=1, kjpindex
       
       zprecip_rain(ik) = precip_rain(ik)*xrdt
       zprecip_snow(ik) = precip_snow(ik)*xrdt
       zcdrag(ik)       = cdrag(ik)
       
    ENDDO
 
    ! Fields for deposit variables : to be transport from chemistry model by GCM to ORCHIDEE.
    ! There are currently no fields to be transported into ORCHIDEE in this way
    DO i = 1, nb_fields_in
       WRITE(numout,*) i," Champ = ",TRIM(field_in_names(i)) 
       SELECT CASE(TRIM(field_in_names(i)))
       CASE DEFAULT 
          CALL ipslerr_p (3,'intsurf_gathered', &
               'You ask in GCM an unknown field '//TRIM(field_in_names(i))//&
               ' to give to ORCHIDEE for this specific version.',&
               'This model won''t be able to continue.', &
               '(check your tracer parameters in GCM)')
       END SELECT
    ENDDO

    !
    ! 3. save the grid
    !
    landpoints(:)=(/ ( REAL(ik), ik=1,kjpindex ) /)
    CALL histwrite_p(hist_id, 'LandPoints',  itau_sechiba+1, landpoints, kjpindex, kindex)
    CALL histwrite_p(hist_id, 'Areas',  itau_sechiba+1, area, kjpindex, kindex)
    IF ( ok_stomate ) THEN
       CALL histwrite_p(hist_id_stom, 'Areas',  itau_sechiba+1, area, kjpindex, kindex)
       IF ( hist_id_stom_ipcc > 0 ) &
            CALL histwrite_p(hist_id_stom_IPCC, 'Areas',  itau_sechiba+1, area, kjpindex, kindex)
    ENDIF
    CALL histwrite_p(hist_id, 'Contfrac',  itau_sechiba+1, contfrac, kjpindex, kindex)
    
    ! Syncronize output but only if flag ok_histsync is set to true       
    IF (ok_histsync) THEN
       IF (is_omp_root .AND. hist_id > 0) THEN
          CALL histsync(hist_id)
       END IF
    END IF
    
    IF ( hist2_id > 0 ) THEN
       CALL histwrite_p(hist2_id, 'LandPoints',  itau_sechiba+1, landpoints, kjpindex, kindex)
       CALL histwrite_p(hist2_id, 'Areas',  itau_sechiba+1, area, kjpindex, kindex)
       CALL histwrite_p(hist2_id, 'Contfrac',  itau_sechiba+1, contfrac, kjpindex, kindex)
       
       ! Syncronize output but only if flag ok_histsync is set to true
       IF (ok_histsync .AND. is_omp_root) THEN
          CALL histsync(hist2_id)
       ENDIF
    ENDIF
    !
    
    !
    ! 4. call sechiba for continental points only
    !
    IF ( printlev_loc>=3 ) WRITE(numout,*) 'Before call to sechiba_initialize'

    CALL sechiba_initialize( &
         itau_sechiba, iim_g*jj_nb,  kjpindex,      kindex,      date0_shifted, &
         lalo,         contfrac,     neighbours,    resolution,  zlev,          &
         u,            v,            qair,          t2m_loc,     temp_air,      &
         petAcoef,     peqAcoef,     petBcoef,      peqBcoef,                   &
         zprecip_rain, zprecip_snow, lwdown,        swnet,       swdown,        &
         pb,           rest_id,      hist_id,       hist2_id,                   &
         rest_id_stom, hist_id_stom, hist_id_stom_IPCC,                         &
         zcoastal,     zriver,       ztsol_rad,     zvevapp,     zqsurf,        &
         zz0m,          zz0h,        zalbedo,      zfluxsens,     zfluxlat,    zemis,         &
         znetco2,      zcarblu,      ztemp_sol_new, zcdrag)
    
    IF ( printlev_loc>=3 ) WRITE(numout,*) 'After call to sechiba_initialize'

    !
    ! 6. scatter output fields
    !
    DO ik=1, kjpindex
       z0m(ik)          = zz0m(ik)
       tsol_rad(ik)     = ztsol_rad(ik)
       vevapp(ik)       = zvevapp(ik)/xrdt ! Transform into kg/m^2/s
       temp_sol_new(ik) = ztemp_sol_new(ik)
       qsurf(ik)        = zqsurf(ik)
       albedo(ik,1)     = zalbedo(ik,1)
       albedo(ik,2)     = zalbedo(ik,2)
       fluxsens(ik)     = zfluxsens(ik)
       fluxlat(ik)      = zfluxlat(ik)
       emis(ik)         = zemis(ik)
       cdrag(ik)        = zcdrag(ik)
    ENDDO
    ! z0h is a optional output variable. Check first if it is present in the call from LMDZ.
    IF ( PRESENT(z0h) ) z0h(:) = zz0h(:)


    !Config Key  = TAU_OUTFLOW
    !Config Desc = Number of days over which the coastal- and riverflow will be distributed
    !Config If   = Only in coupled mode
    !Config Def  = 0
    !Config Help = The default value 0 makes the distribution instanteneous
    !Config Units = [days]
    tau_outflow = 0
    CALL getin_p('TAU_OUTFLOW',tau_outflow)
    IF (tau_outflow <=xrdt/one_day) THEN
       coeff_rel = 1.0
    ELSE
       coeff_rel = (1.0 - exp(-xrdt/(tau_outflow*one_day)))
    END IF
    IF (printlev_loc >=2)  WRITE(numout,*) 'tau_outflow, coeff_rel = ', tau_outflow, coeff_rel

    ! Allocate and read riverflow_cpl0 from restart file. Initialize to 0 if the variable was not found in the restart file.
    ALLOCATE(riverflow_cpl0(kjpindex), stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'intersurf_initialize_gathered', 'Error in allocation of riverflow_cpl0','','')
    CALL restget_p (rest_id, 'riverflow_cpl0', nbp_glo, 1, 1, kjit, .TRUE., riverflow_cpl0, "gather", nbp_glo, index_g)
    IF ( ALL(riverflow_cpl0(:) == val_exp) ) riverflow_cpl0(:)=0

    ! Allocate and read coastalflow_cpl0 from restart file. Initialize to 0 if the variable was not found in the restart file.
    ALLOCATE(coastalflow_cpl0(kjpindex), stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'intersurf_initialize_gathered', 'Error in allocation of coastalflow_cpl0','','')
    CALL restget_p (rest_id, 'coastalflow_cpl0', nbp_glo, 1, 1, kjit, .TRUE., coastalflow_cpl0, "gather", nbp_glo, index_g)
    IF ( ALL(coastalflow_cpl0(:) == val_exp) ) coastalflow_cpl0(:)=0

    ! Do not applay the filter now in initialization phase. 
    ! These variables will never be used anyway in the initialization phase. 
    ! Transform into m^3/s
    riverflow_cpl = zriver/xrdt
    coastalflow_cpl = zcoastal/xrdt
    
    ! Fields for emission variables : to be transport by GCM to chemistry model.
    DO i = 1, nb_fields_out
       SELECT CASE(TRIM(field_out_names(i)))
       CASE("fCO2_land") 
          fields_out(:,i)=znetco2(:)
       CASE("fCO2_land_use")
          fields_out(:,i)=zcarblu(:)
       CASE DEFAULT 
          CALL ipslerr_p (3,'intsurf_gathered', &
               'You ask from GCM an unknown field '//TRIM(field_out_names(i))//&
               ' to ORCHIDEE for this specific version.',&
               'This model won''t be able to continue.', &
               '(check your tracer parameters in GCM)')
       END SELECT
    END  DO

    ! Copy the number of vegetation types to local output variable
    IF (PRESENT(nvm_out)) nvm_out=nvm

    IF(is_root_prc) CALL getin_dump
    lstep_init_intersurf = .FALSE.
    
    CALL ipslnlf_p(new_number=old_fileout)
    !
    !! Change back to be in the LMDZ context for XIOS
    !
    CALL xios_orchidee_change_context("LMDZ")

    IF (printlev_loc >= 2) WRITE (numout,*) 'End intersurf_initialize_gathered'
    IF (printlev_loc >= 1) WRITE (numout,*) 'Initialization phase for ORCHIDEE is finished.'

  END SUBROUTINE intersurf_initialize_gathered


!!  =============================================================================================================================
!! SUBROUTINE:    intersurf_main_gathered
!!
!>\BRIEF	  Main subroutine to call ORCHIDEE from the gcm (LMDZ) using variables on a 1D grid with only land points.
!!
!! DESCRIPTION:	  This subroutine is the main interface for ORCHIDEE when it is called from the gcm (LMDZ).
!!                The variables are all gathered before entering this subroutine on the 1D grid with only landpoints.
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE intersurf_main_gathered (kjit, iim_glo, jjm_glo, kjpindex, kindex, xrdt, &
       lrestart_read, lrestart_write, latlon, zcontfrac, zneighbours, zresolution, date0, &
       zlev,  u, v, qair, temp_air, epot_air, ccanopy, &
       cdrag, petAcoef, peqAcoef, petBcoef, peqBcoef, &
       precip_rain, precip_snow, lwdown, swnet, swdown, pb, &
       vevapp, fluxsens, fluxlat, coastalflow_cpl, riverflow_cpl, &
       tsol_rad, temp_sol_new, qsurf, albedo, emis, z0m,lon_scat_g, lat_scat_g, q2m, t2m, z0h, &
       veget, lai, height, &
       field_out_names, fields_out, field_in_names, fields_in, &
       coszang)  

    USE mod_orchidee_para
    IMPLICIT NONE

    !! 0. Variable and parameter declaration
    !! 0.1 Input variables
    INTEGER(i_std),INTENT (in)                            :: kjit            !! Time step number
    INTEGER(i_std),INTENT (in)                            :: iim_glo         !! Dimension of global fields
    INTEGER(i_std),INTENT (in)                            :: jjm_glo         !! Dimension of global fields
    INTEGER(i_std),INTENT (in)                            :: kjpindex        !! Number of continental points
    REAL(r_std),INTENT (in)                               :: xrdt            !! Time step in seconds
    LOGICAL, INTENT (in)                                  :: lrestart_read   !! Logical for _restart_ file to read
    LOGICAL, INTENT (in)                                  :: lrestart_write  !! Logical for _restart_ file to write'
    REAL(r_std), INTENT (in)                              :: date0           !! Date at which kjit = 0
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)      :: kindex          !! Index for continental points
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: u               !! Lowest level wind speed
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: v               !! Lowest level wind speed 
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: zlev            !! Height of first layer
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: qair            !! Lowest level specific humidity
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: precip_rain     !! Rain precipitation
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: precip_snow     !! Snow precipitation
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: lwdown          !! Down-welling long-wave flux 
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: swnet           !! Net surface short-wave flux
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: swdown          !! Downwelling surface short-wave flux
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: temp_air        !! Air temperature in Kelvin
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: epot_air        !! Air potential energy
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: ccanopy         !! CO2 concentration in the canopy
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: petAcoef        !! Coeficients A from the PBL resolution
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: peqAcoef        !! One for T and another for q
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: petBcoef        !! Coeficients B from the PBL resolution
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: peqBcoef        !! One for T and another for q
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: pb              !! Surface pressure
    REAL(r_std),DIMENSION (kjpindex,2), INTENT(in)        :: latlon          !! Geographical coordinates
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: zcontfrac       !! Fraction of continent
    INTEGER(i_std),DIMENSION (kjpindex,NbNeighb), INTENT(in):: zneighbours     !! neighbours
    REAL(r_std),DIMENSION (kjpindex,2), INTENT(in)        :: zresolution     !! size of the grid box
    REAL(r_std),DIMENSION (iim_glo,jjm_glo), INTENT(IN)   :: lon_scat_g      !! The scattered values for longitude 
    REAL(r_std),DIMENSION (iim_glo,jjm_glo), INTENT(IN)   :: lat_scat_g      !! The scattered values for latitude

    !! 0.2 Output variables
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)         :: z0m             !! Surface roughness for momentum
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)         :: coastalflow_cpl !! Diffuse flow of water into the ocean, time filtered (m^3/s)
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)         :: riverflow_cpl   !! Largest rivers flowing into the ocean, time filtered (m^3/s)
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)         :: tsol_rad        !! Radiative surface temperature
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)         :: vevapp          !! Total of evaporation
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)         :: temp_sol_new    !! New soil temperature
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)         :: qsurf           !! Surface specific humidity
    REAL(r_std),DIMENSION (kjpindex,2), INTENT(out)       :: albedo          !! Albedo
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)         :: fluxsens        !! Sensible chaleur flux
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)         :: fluxlat         !! Latent chaleur flux
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)         :: emis            !! Emissivity

    !! 0.3 Modified variables
    REAL(r_std),DIMENSION (kjpindex), INTENT(inout)       :: cdrag           !! Cdrag

    !! 0.4 Optional input and output variables
    REAL(r_std),DIMENSION (kjpindex), INTENT(in), OPTIONAL:: q2m             !! Surface specific humidity
    REAL(r_std),DIMENSION (kjpindex), INTENT(in), OPTIONAL:: t2m             !! Surface air temperature
    REAL(r_std),DIMENSION (kjpindex), INTENT(out),OPTIONAL:: z0h             !! Surface roughness for heat
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(out), OPTIONAL :: veget     !! Fraction of vegetation type (unitless, 0-1) 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(out), OPTIONAL :: lai       !! Leaf area index (m^2 m^{-2}
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(out), OPTIONAL :: height    !! Vegetation Height (m)
    REAL(r_std), DIMENSION(kjpindex), OPTIONAL, INTENT(in):: coszang         !! Cosine of the solar zenith angle (unitless)
    CHARACTER(LEN=*),DIMENSION(:), OPTIONAL, INTENT(IN)   :: field_in_names  !! Names for deposit variables to be transported
                                                                             !! from chemistry model by GCM to ORCHIDEE
    REAL(r_std),DIMENSION(:,:), OPTIONAL, INTENT(IN)      :: fields_in       !! Fields for deposit variables to be transported 
                                                                             !! from chemistry model by GCM to ORCHIDEE
    CHARACTER(LEN=*),DIMENSION(:), OPTIONAL, INTENT(IN)   :: field_out_names !! Names for emission variables to be transported 
                                                                             !! to chemistry model by GCM from ORCHIDEE
    REAL(r_std),DIMENSION(:,:), OPTIONAL, INTENT(OUT)     :: fields_out      !! Fields for emission variables to be transported
                                                                             !! to chemistry model by GCM from ORCHIDEE

    !! 0.5 Local variables
    REAL(r_std),DIMENSION (kjpindex)                      :: zccanopy      !! Work array to keep ccanopy
    REAL(r_std),DIMENSION (kjpindex)                      :: zprecip_rain  !! Work array to keep precip_rain
    REAL(r_std),DIMENSION (kjpindex)                      :: zprecip_snow  !! Work array to keep precip_snow
    REAL(r_std),DIMENSION (kjpindex)                      :: zz0m, zz0h    !! Work array to keep zz0m, zz0h
    REAL(r_std),DIMENSION (kjpindex)                      :: zcdrag        !! Work array for surface drag
    REAL(r_std),DIMENSION (kjpindex)                      :: zcoastal      !! Work array to keep coastal flow
    REAL(r_std),DIMENSION (kjpindex)                      :: zriver        !! Work array to keep river out flow
    REAL(r_std),DIMENSION (kjpindex)                      :: znetco2       !! Work array to keep netco2flux
    REAL(r_std),DIMENSION (kjpindex)                      :: zcarblu       !! Work array to keep fco2_land_use
    REAL(r_std),DIMENSION (kjpindex)                      :: ztsol_rad     !! Work array to keep tsol_rad
    REAL(r_std),DIMENSION (kjpindex)                      :: zvevapp       !! Work array to keep vevapp
    REAL(r_std),DIMENSION (kjpindex)                      :: ztemp_sol_new !! Work array to keep temp_sol_new
    REAL(r_std),DIMENSION (kjpindex)                      :: zqsurf        !! Work array to keep qsurf
    REAL(r_std),DIMENSION (kjpindex,2)                    :: zalbedo       !! Work array to keep albedo
    REAL(r_std),DIMENSION (kjpindex)                      :: zfluxsens     !! Work array to keep fluxsens
    REAL(r_std),DIMENSION (kjpindex)                      :: zfluxlat      !! Work array to keep fluxlat
    REAL(r_std),DIMENSION (kjpindex)                      :: zemis         !! Work array to keep emis
    REAL(r_std),DIMENSION (kjpindex)                      :: zcoszang      !! Work array to keep coszang
    REAL(r_std),ALLOCATABLE, DIMENSION (:)                :: soilth_lev    !! Vertical soil axis for thermal scheme (m)
    REAL(r_std),DIMENSION (kjpindex)                      :: q2m_loc       !! Work array for q2m or qair
    REAL(r_std),DIMENSION (kjpindex)                      :: t2m_loc       !! Work array for t2m or temp_air
    REAL(r_std),DIMENSION (kjpindex,nvm)                  :: zveget        !! Work array to keep veget
    REAL(r_std),DIMENSION (kjpindex,nvm)                  :: zlai          !! Work array to keep lai
    REAL(r_std),DIMENSION (kjpindex,nvm)                  :: zheight       !! Work array to keep height
    INTEGER(i_std)                                        :: i, j, ik
    INTEGER(i_std)                                        :: ier
    INTEGER(i_std)                                        :: itau_sechiba
    INTEGER                                               :: old_fileout   !! old Logical Int for std IO output
    REAL,ALLOCATABLE,DIMENSION(:,:)                       :: lalo_mpi
    REAL(r_std),DIMENSION (kjpindex)                      :: landpoints    !! Local landpoints vector
    

    IF (printlev_loc >= 3) WRITE(numout,*) 'Start intersurf_main_gathered'
    CALL ipslnlf_p(new_number=numout, old_number=old_fileout)
    
    IF (lstep_init_intersurf) THEN
       ! Test if q2m and t2m are present
       IF (PRESENT(q2m) .AND. PRESENT(t2m)) THEN
          ok_q2m_t2m=.TRUE.
       ELSE
          ok_q2m_t2m=.FALSE.
       ENDIF

       ! Test if field_out_names and field_in_names are present and if so, count 
       ! the number of extra fields to exchange.
       IF (PRESENT(field_out_names)) THEN
          nb_fields_out=SIZE(field_out_names)
       ELSE
          nb_fields_out=0
       ENDIF

       IF (PRESENT(field_in_names)) THEN
          nb_fields_in=SIZE(field_in_names)
       ELSE
          nb_fields_in=0
       ENDIF
    
       CALL intersurf_initialize_gathered (kjit, iim_glo, jjm_glo, kjpindex, kindex, xrdt, &
            lrestart_read, lrestart_write, latlon, zcontfrac, zneighbours, zresolution, date0, &
            zlev,  u, v, qair, temp_air, epot_air, &
            cdrag, petAcoef, peqAcoef, petBcoef, peqBcoef, &
            precip_rain, precip_snow, lwdown, swnet, swdown, pb, &
            vevapp, fluxsens, fluxlat, coastalflow_cpl, riverflow_cpl, &
            tsol_rad, temp_sol_new, qsurf, albedo, emis, z0m,lon_scat_g, lat_scat_g, &
            q2m, t2m, zz0h, &
            field_out_names=field_out_names, fields_out=fields_out, &
            field_in_names=field_in_names, fields_in=fields_in)

       ! z0h is a optional output variable. Check first if it is present in the call from LMDZ.
       IF ( PRESENT(z0h) ) z0h(:) = zz0h(:)
       ! Return from subroutine intersurf_main_gathered
       RETURN
    END IF

    !
    !! Change to be in the orchidee context for XIOS
    !
    CALL xios_orchidee_change_context("orchidee")
    
    !
    !  Shift the time step to phase the two models
    !
    itau_sechiba = kjit + itau_offset
    
    ! Update the calendar in xios by sending the new time step
    CALL xios_orchidee_update_calendar(itau_sechiba)
    
    ! Update the calendar and all time variables in module time
    CALL time_nextstep(itau_sechiba)
    
    !
    ! 1. Just change the units of some input fields
    !
    DO ik=1, kjpindex
       
       zprecip_rain(ik) = precip_rain(ik)*xrdt
       zprecip_snow(ik) = precip_snow(ik)*xrdt
       zcdrag(ik)       = cdrag(ik)
       
    ENDDO

    !>> VOC in coupled mode 
    IF ( PRESENT(coszang) )  THEN 
       zcoszang(:) = coszang(:)
    ELSE
       zcoszang(:) = zero
    ENDIF
 
    ! Fields for deposit variables : to be transport from chemistry model by GCM to ORCHIDEE.
    DO i = 1, nb_fields_in
       WRITE(numout,*) i," Champ = ",TRIM(field_in_names(i)) 
       SELECT CASE(TRIM(field_in_names(i)))
       CASE DEFAULT 
          CALL ipslerr_p (3,'intsurf_gathered', &
               'You ask in GCM an unknown field '//TRIM(field_in_names(i))//&
               ' to give to ORCHIDEE for this specific version.',&
               'This model won''t be able to continue.', &
               '(check your tracer parameters in GCM)')
       END SELECT
    ENDDO

    !
    ! 2. modification of co2
    !
    IF ( fatmco2 ) THEN
       zccanopy(:) = atmco2
       WRITE (numout,*) 'Modification of the ccanopy value. CO2 = ',atmco2
    ELSE
       zccanopy(:) = ccanopy(:)
    ENDIF

    !
    ! 4. call sechiba for continental points only
    !
    IF ( printlev_loc >= 3 ) WRITE(numout,*) 'Before call to sechiba_main from intersurf_main_gathered'
   

    IF (ok_q2m_t2m) THEN
       t2m_loc=t2m
       q2m_loc=q2m
    ELSE
       t2m_loc=temp_air
       q2m_loc=qair
    END IF

    CALL sechiba_main (itau_sechiba, iim_g*jj_nb, kjpindex, kindex, date0_shifted, &
         lrestart_read, lrestart_write, &
         lalo, contfrac, neighbours, resolution, &
         zlev, u, v, qair, q2m_loc, t2m_loc, temp_air, epot_air, zccanopy, &
         zcdrag, petAcoef, peqAcoef, petBcoef, peqBcoef, &
         zprecip_rain ,zprecip_snow,  lwdown, swnet, swdown, zcoszang, pb, &
         zvevapp, zfluxsens, zfluxlat, zcoastal, zriver, znetco2, zcarblu, &
         ztsol_rad, ztemp_sol_new, zqsurf, zalbedo, zemis, zz0m, zz0h, &
         zveget, zlai, zheight, &
         rest_id, hist_id, hist2_id, rest_id_stom, hist_id_stom, hist_id_stom_IPCC ) 

    IF ( printlev_loc>=3 ) WRITE(numout,*) 'After call to sechiba_main'

    !
    ! 6. scatter output fields
    !
    DO ik=1, kjpindex
       z0m(ik)          = zz0m(ik)
       tsol_rad(ik)     = ztsol_rad(ik)
       temp_sol_new(ik) = ztemp_sol_new(ik)
       qsurf(ik)        = zqsurf(ik)
       albedo(ik,1)     = zalbedo(ik,1)
       albedo(ik,2)     = zalbedo(ik,2)
       fluxsens(ik)     = zfluxsens(ik)
       fluxlat(ik)      = zfluxlat(ik)
       emis(ik)         = zemis(ik)
       cdrag(ik)        = zcdrag(ik)
       ! Transform the water fluxes into Kg/m^2/s
       vevapp(ik)       = zvevapp(ik)/xrdt
    ENDDO

    ! Copy variables only if the optional variables are present in the call to the subroutine
    IF (PRESENT(veget))  veget(:,:)  = zveget(:,:) 
    IF (PRESENT(lai))    lai(:,:)    = zlai(:,:) 
    IF (PRESENT(height)) height(:,:) = zheight(:,:) 

    ! Applay time filter to distribut the river- and coastalflow over a longer time period.
    ! When coeff_rel=1(default case when tau_outflow=0), the distribution is instanteneous. 
    ! Use TAU_OUTFLOW in run.def to set the number of days of distribution.
    ! The water fluxes zriver and zcoastal coming from sechiba are at the same time transfromed 
    ! from m^3/dt_sechiba into m^3/s by dividing with the sechiba time step (xrdt).
    riverflow_cpl = coeff_rel*zriver/xrdt + (1.-coeff_rel)*riverflow_cpl0
    riverflow_cpl0 = riverflow_cpl

    coastalflow_cpl = coeff_rel*zcoastal/xrdt + (1.-coeff_rel)*coastalflow_cpl0
    coastalflow_cpl0 = coastalflow_cpl 

    
    ! z0h is a optional output variable. Check first if it is present in the call from LMDZ.
    IF ( PRESENT(z0h) ) z0h(:) = zz0h(:)
       
    CALL xios_orchidee_send_field("LandPoints" ,(/ ( REAL(ik), ik=1,kjpindex ) /))
    CALL xios_orchidee_send_field("areas", area)
    CALL xios_orchidee_send_field("contfrac",contfrac)
    CALL xios_orchidee_send_field("temp_air",temp_air)
    CALL xios_orchidee_send_field("qair",qair)
    CALL xios_orchidee_send_field("swnet",swnet)
    CALL xios_orchidee_send_field("swdown",swdown)
    CALL xios_orchidee_send_field("pb",pb)
    CALL xios_orchidee_send_field("riverflow_cpl",riverflow_cpl)
    CALL xios_orchidee_send_field("coastalflow_cpl",coastalflow_cpl)

 
    IF (ok_q2m_t2m) THEN
       CALL xios_orchidee_send_field("t2m",t2m)
       CALL xios_orchidee_send_field("q2m",q2m)
    ELSE
       CALL xios_orchidee_send_field("t2m",temp_air)
       CALL xios_orchidee_send_field("q2m",qair)
    ENDIF
    
    IF ( .NOT. almaoutput ) THEN
       !
       !  scattered during the writing
       !           
       CALL histwrite_p (hist_id, 'evap',     itau_sechiba, zvevapp, kjpindex, kindex)
       CALL histwrite_p (hist_id, 'coastalflow',itau_sechiba, zcoastal, kjpindex, kindex)
       CALL histwrite_p (hist_id, 'riverflow',itau_sechiba, zriver, kjpindex, kindex)
       CALL histwrite_p (hist_id, 'temp_sol', itau_sechiba, ztemp_sol_new, kjpindex, kindex)
       CALL histwrite_p (hist_id, 'tsol_max', itau_sechiba, ztemp_sol_new, kjpindex, kindex)
       CALL histwrite_p (hist_id, 'tsol_min', itau_sechiba, ztemp_sol_new, kjpindex, kindex)
       CALL histwrite_p (hist_id, 'fluxsens', itau_sechiba, zfluxsens, kjpindex, kindex)
       CALL histwrite_p (hist_id, 'fluxlat',  itau_sechiba, zfluxlat,  kjpindex, kindex)
       CALL histwrite_p (hist_id, 'swnet',    itau_sechiba, swnet,    kjpindex, kindex)
       CALL histwrite_p (hist_id, 'swdown',   itau_sechiba, swdown,   kjpindex, kindex)
       CALL histwrite_p (hist_id, 'alb_vis',  itau_sechiba, zalbedo(:,1), kjpindex, kindex)
       CALL histwrite_p (hist_id, 'alb_nir',  itau_sechiba, zalbedo(:,2), kjpindex, kindex)
       CALL histwrite_p (hist_id, 'tair',     itau_sechiba, temp_air, kjpindex, kindex)
       CALL histwrite_p (hist_id, 'qair',     itau_sechiba, qair, kjpindex, kindex)
       IF (ok_q2m_t2m) THEN
          CALL histwrite_p (hist_id, 't2m',      itau_sechiba, t2m, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'q2m',      itau_sechiba, q2m, kjpindex, kindex)
       ELSE
          CALL histwrite_p (hist_id, 't2m',      itau_sechiba, temp_air, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'q2m',      itau_sechiba, qair, kjpindex, kindex)
       ENDIF
       IF ( hist2_id > 0 ) THEN
          CALL histwrite_p (hist2_id, 'evap',     itau_sechiba, zvevapp, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'coastalflow',itau_sechiba, zcoastal, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'riverflow',itau_sechiba, zriver, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'temp_sol', itau_sechiba, temp_sol_new, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'tsol_max', itau_sechiba, temp_sol_new, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'tsol_min', itau_sechiba, temp_sol_new, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'fluxsens', itau_sechiba, zfluxsens, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'fluxlat',  itau_sechiba, zfluxlat,  kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'swnet',    itau_sechiba, swnet,    kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'swdown',   itau_sechiba, swdown,   kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'alb_vis',  itau_sechiba, zalbedo(:,1), kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'alb_nir',  itau_sechiba, zalbedo(:,2), kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'tair',     itau_sechiba, temp_air, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'qair',     itau_sechiba, qair, kjpindex, kindex)
          IF (ok_q2m_t2m) THEN
             CALL histwrite_p (hist2_id, 't2m',      itau_sechiba, t2m, kjpindex, kindex)
             CALL histwrite_p (hist2_id, 'q2m',      itau_sechiba, q2m, kjpindex, kindex)
          ELSE
             CALL histwrite_p (hist2_id, 't2m',      itau_sechiba, temp_air, kjpindex, kindex)
             CALL histwrite_p (hist2_id, 'q2m',      itau_sechiba, qair, kjpindex, kindex)
          ENDIF
       ENDIF
    ELSE
       CALL histwrite_p (hist_id, 'Evap', itau_sechiba, zvevapp, kjpindex, kindex)
       CALL histwrite_p (hist_id, 'SWnet',    itau_sechiba, swnet, kjpindex, kindex)
       CALL histwrite_p (hist_id, 'Qh', itau_sechiba, zfluxsens, kjpindex, kindex)
       CALL histwrite_p (hist_id, 'Qle',  itau_sechiba, zfluxlat, kjpindex, kindex)
       CALL histwrite_p (hist_id, 'AvgSurfT', itau_sechiba, ztemp_sol_new, kjpindex, kindex)
       CALL histwrite_p (hist_id, 'RadT', itau_sechiba, ztemp_sol_new, kjpindex, kindex)
       IF ( hist2_id > 0 ) THEN
          CALL histwrite_p (hist2_id, 'Evap', itau_sechiba, zvevapp, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'SWnet',    itau_sechiba, swnet, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'Qh', itau_sechiba, zfluxsens, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'Qle',  itau_sechiba, zfluxlat, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'AvgSurfT', itau_sechiba, ztemp_sol_new, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'RadT', itau_sechiba, ztemp_sol_new, kjpindex, kindex)
       ENDIF
    ENDIF
    
    ! Syncronize output but only if flag ok_histsync is set to true
    IF (ok_histsync .AND. is_omp_root) THEN
       IF ( (dw .EQ. xrdt) .AND. hist_id > 0 ) THEN
          CALL histsync(hist_id)
       ENDIF
    ENDIF
    
  
    ! Fields for emission variables : to be transport by GCM to chemistry model.
    DO i = 1, nb_fields_out
       SELECT CASE(TRIM(field_out_names(i)))
       CASE("fCO2_land") 
          fields_out(:,i)=znetco2(:)
       CASE("fCO2_land_use")
          fields_out(:,i)=zcarblu(:)
       CASE DEFAULT 
          CALL ipslerr_p (3,'intsurf_gathered', &
            &          'You ask from GCM an unknown field '//TRIM(field_out_names(i))//&
            &          ' to ORCHIDEE for this specific version.',&
            &          'This model won''t be able to continue.', &
            &          '(check your tracer parameters in GCM)')
       END SELECT
    ENDDO


    ! Write variables to restart file
    IF (lrestart_write) THEN
       CALL restput_p (rest_id, 'riverflow_cpl0', nbp_glo, 1, 1, kjit, riverflow_cpl0, 'scatter',  nbp_glo, index_g)
       CALL restput_p (rest_id, 'coastalflow_cpl0', nbp_glo, 1, 1, kjit, coastalflow_cpl0, 'scatter',  nbp_glo, index_g)
    END IF

    !
    CALL ipslnlf_p(new_number=old_fileout)
    !        

    !
    !! Finalize the XIOS orchidee context if it is the last call
    !
    IF (lrestart_write) THEN
       CALL xios_orchidee_context_finalize
    END IF
    !
    !! Change back to be in the LMDZ context for XIOS
    !
    CALL xios_orchidee_change_context("LMDZ")

    IF (printlev_loc >= 3) WRITE (numout,*) 'End intersurf_main_gathered'

  END SUBROUTINE intersurf_main_gathered

!! ================================================================================================================================
!! SUBROUTINE   : intersurf_clear
!!
!>\BRIEF         Clear intersurf module and underlaying modules
!!
!! DESCRIPTION  :  Deallocate memory and reset initialization variables to there original values. 
!!                 Call the clearing for sechiba module.
!!
!_ ================================================================================================================================
  SUBROUTINE intersurf_clear
    CALL sechiba_clear
  END SUBROUTINE intersurf_clear

END MODULE intersurf
