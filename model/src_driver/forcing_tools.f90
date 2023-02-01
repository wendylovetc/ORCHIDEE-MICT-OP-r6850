!  ==============================================================================================================================\n
!  MODULE forcing_tools	: This module concentrates on the temporal interpolation of the forcing for ORCHIDEE. 
!                         It provides basic service for the grid when this is provided in the forcing file. The main
!                         work for the grid is done in glogrid.f90. The approach of forcing_tools to handle the time
!                         aspect of the forcing is to read as many time steps as possible in memory and then
!                         interpolate that to the time interval requested by the calling program.
!                         The data is read on root_proc but then distributed over all processors according to the
!                         domain decomposition of ORCHIDEE. This allows to be more efficient in the temporal interpolation.
!                         It is important to keep in mind that forcing_tools works on time intervals. So the request for data 
!                         of ORCHIDEE as to be for an interval and the forcing file needs to have a description of the time interval
!                         over which the forcing is valid.
!                         The general description of how the attributes needed in the netCDF file for describing the cell_methods
!                         for time are provided in this document :
!                          https://forge.ipsl.jussieu.fr/orchidee/attachment/wiki/Documentation/Forcings/Description_Forcing_Files.pdf
!
!                         The most important routines of foring_tools are forcing_open and forcing_getvalues
!
!                       forcing_integration_time : Computes the interval over which the simulation should be carried out.
!                       forcing_open : Opens the forcing files and extract the main information.
!                       forcing_getvalues : Gets the forcing data for a time interval.
!                       forcing_close : Closes the forcing file
!                       forcing_printdate : A tool to print the dates in human readable form.
!                       forcing_printpoint : Print the values for a given point in time.
!                       forcing_givegridsize : Allows other programs to get the dimensions of the forcing grid.
!                       forcing_getglogrid : Allows other programs to get the spatial grid of the forcing.
!                       forcing_givegrid : Returns the description of the grid.
!                       forcing_zoomgrid : Extract a sub-region of the forcing grid.
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
!!
MODULE forcing_tools
  !
  USE defprec
  USE netcdf
  !
  USE ioipsl
  USE constantes
  USE solar
  !
  USE mod_orchidee_para
  !
  IMPLICIT NONE
  ! 
  PRIVATE
  PUBLIC :: forcing_open, forcing_close, forcing_printdate, forcing_getvalues, forcing_printpoint,&
       &    forcing_getglogrid, forcing_givegridsize, forcing_givegrid, forcing_zoomgrid, forcing_integration_time
  !
  !
  !
  INTERFACE forcing_reindex
     MODULE PROCEDURE forcing_reindex3d, forcing_reindex2dt, forcing_reindex2d, forcing_reindex1d, &
          &           forcing_reindex2to1, forcing_reindex1to2
  END INTERFACE forcing_reindex
  !
  INTERFACE forcing_printpoint
     MODULE PROCEDURE forcing_printpoint_forgrid, forcing_printpoint_forgrid2d, forcing_printpoint_gen
  END INTERFACE forcing_printpoint
  !
  INTERFACE forcing_getvalues
     MODULE PROCEDURE forcing_getvalues1d, forcing_getvalues2d
  END INTERFACE forcing_getvalues
  !
  ! This PARAMETER essentially manages the memory usage of the module as it
  ! determines how much of the forcing will be uploaded from the netCDF file into
  ! memory.
  !
  INTEGER(i_std), PARAMETER :: slab_size_max=250
  !
  ! Time variables, all in Julian days
  !
  INTEGER(i_std), PARAMETER :: nbtmethods=4
  INTEGER(i_std), SAVE :: nbtax
  INTEGER(i_std), SAVE :: nb_forcing_steps
  REAL(r_std), SAVE :: global_start_date, global_end_date, forcing_tstep_ave
  REAL(r_std), SAVE :: dt_sechiba_keep
  !
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)     :: time
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:,:)   :: time_bounds
  CHARACTER(LEN=20), SAVE, ALLOCATABLE, DIMENSION(:) :: time_axename, time_cellmethod
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)       :: preciptime
  INTEGER(i_std), SAVE, ALLOCATABLE, DIMENSION(:)    :: time_sourcefile
  INTEGER(i_std), SAVE, ALLOCATABLE, DIMENSION(:,:)  :: time_id
  LOGICAL, SAVE :: end_of_file
  !
  ! Forcing file information
  !
  INTEGER(i_std), SAVE                                :: nb_forcefile=0
  CHARACTER(LEN=100), SAVE, ALLOCATABLE, DIMENSION(:) :: forfilename
  INTEGER(i_std), SAVE, ALLOCATABLE, DIMENSION(:)     :: force_id, id_unlim
  INTEGER(i_std), SAVE, ALLOCATABLE, DIMENSION(:)     :: nb_atts, ndims, nvars
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)        :: convtosec
  INTEGER(i_std), SAVE, ALLOCATABLE, DIMENSION(:)     :: nbtime_perfile
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)      :: date0_file
  REAL(r_std), SAVE                                   :: startdate, forcingstartdate
  !
  ! Descrition of global grid
  !
  INTEGER(i_std), SAVE :: iim_glo, jjm_glo, nbpoint_glo, nbland_glo
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)   :: lon_glo, lat_glo
  INTEGER(i_std), SAVE, ALLOCATABLE, DIMENSION(:,:):: mask_glo
  INTEGER(i_std), SAVE, ALLOCATABLE, DIMENSION(:)  :: lindex_glo
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)     :: contfrac_glo
  LOGICAL, SAVE                                    :: compressed
  !
  ! Descritpion of zoomed grid
  !
  LOGICAL, SAVE :: zoom_forcing = .FALSE.
  INTEGER(i_std), SAVE :: iim_loc, jjm_loc, nbpoint_loc, nbland_loc
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)    :: lon_loc, lat_loc
  INTEGER(i_std), SAVE, ALLOCATABLE, DIMENSION(:)   :: lindex_loc
  INTEGER(i_std), SAVE, ALLOCATABLE, DIMENSION(:,:) :: mask_loc
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)    :: area_loc
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)      :: contfrac_loc
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:,:,:):: corners_loc
  ! Number of land points per proc
  INTEGER(i_std), SAVE :: nbpoint_proc
  INTEGER(i_std), SAVE, ALLOCATABLE, DIMENSION(:) :: glolindex_proc
  LOGICAL, SAVE :: landonly
  !-
  !- Heigh controls and data 
  !- 
  LOGICAL, SAVE                            :: zfixed, zsigma, zhybrid, zlevels, zheight 
  LOGICAL, SAVE                            :: zsamelev_uv 
  REAL, SAVE                               :: zlev_fixed, zlevuv_fixed 
  REAL, SAVE                               :: zhybrid_a, zhybrid_b 
  REAL, SAVE                               :: zhybriduv_a, zhybriduv_b
  LOGICAL, SAVE                            :: lwdown_cons
  !
  ! Forcing variables to be read and stored 
  !
  ! At 3000 we can fit in the slab an entire year of 3 hourly forcing.
  INTEGER(i_std), SAVE :: slab_size=-1
  INTEGER(i_std), SAVE :: current_offset=1
  INTEGER(i_std), SAVE :: position_slab(2)
  CHARACTER(LEN=20), SAVE :: calendar
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)    :: tair_slab, qair_slab
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)      :: time_tair, time_qair
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)    :: timebnd_tair, timebnd_qair
  !
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)    :: rainf_slab, snowf_slab
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)      :: time_precip
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)    :: timebnd_precip
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)      :: preciptime_slab             !! Variable needed to keep track of how much rainfall was already distributed
  !
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)    :: swdown_slab, lwdown_slab
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)      :: time_swdown, time_lwdown
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)    :: timebnd_swdown, timebnd_lwdown
  !
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)    :: u_slab, v_slab, ps_slab
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)      :: time_u, time_v, time_ps
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)    :: timebnd_u, timebnd_v, timebnd_ps
  !
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)    :: ztq_slab, zuv_slab
  INTEGER(i_std), SAVE, ALLOCATABLE, DIMENSION(:)   :: reindex_glo, reindex_loc
  INTEGER(i_std), SAVE, ALLOCATABLE, DIMENSION(:,:) :: reindex2d_loc
  INTEGER(i_std), SAVE, ALLOCATABLE, DIMENSION(:,:) :: origind
  !
  INTEGER(i_std), SAVE                              :: ncdfstart, ncdfcount
  !
CONTAINS
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_integration_time
!!
!>\BRIEF   Computes the interval over which the simulation should be carried out   
!!
!! DESCRIPTION:	 This routing will get the following parameters from the run.def : 'START_DATE', 'END_DATE' and 'DT_SECHIBA'.
!!               It allows to define the integration time of ORCHIDEE and later it will be used to verify that we have
!!               the needed data in the forcing files to perform this simulation.
!!
!! \n
!_ ==============================================================================================================================
!!
  SUBROUTINE forcing_integration_time(date_start, dt, nbdt)
    !
    !
    ! This subroutine gets the start date of the simulation, the time step and the number 
    ! of time steps we need to do until the end of the simulations.
    !
    !
    !
    REAL(r_std), INTENT(out)                     :: date_start     !! The date at which the simulation starts
    REAL(r_std), INTENT(out)                     :: dt             !! Time step length in seconds
    INTEGER(i_std), INTENT(out)                  :: nbdt           !! Number of timesteps to be executed
    !
    ! Local
    !
    CHARACTER(LEN=20) :: str_sdate(2), str_edate(2), tmpstr
    INTEGER(i_std) :: s_year, s_month, s_day, e_year, e_month, e_day
    INTEGER(i_std) :: seci, hours, minutes
    REAL(r_std) :: s_sec, e_sec, dateend, diff_sec, date_end
    INTEGER(i_std) :: i, ic
    !
    !Config Key  = START_DATE
    !Config Desc = Date at which the simulation starts
    !Config Def  = NONE
    !Config Help = The format is the same as in the CF convention : 1999-09-13 12:0:0
    str_sdate = " "
    CALL getin('START_DATE',str_sdate)
    !
    IF ( (INDEX(str_sdate(1),"-") .NE. INDEX(str_sdate(1),"-", .TRUE.)) .AND. &
         &  (INDEX(str_sdate(2),":") .NE. INDEX(str_sdate(2),":", .TRUE.)) ) THEN
       DO i=1,2
          tmpstr = str_sdate(1)
          ic = INDEX(tmpstr,"-")
          tmpstr(ic:ic) = " "
          str_sdate(1) = tmpstr
          tmpstr = str_sdate(2)
          ic = INDEX(tmpstr,":")
          tmpstr(ic:ic) = " "
          str_sdate(2) = tmpstr
       ENDDO
       READ (str_sdate(1),*) s_year, s_month, s_day
       READ (str_sdate(2),*) hours, minutes, seci
       s_sec = hours*3600. + minutes*60. + seci
    ELSE
       CALL ipslerr(3, "forcing_integration_time", "START_DATE incorrectly specified in run.def", str_sdate(1), str_sdate(2))
    ENDIF
    CALL ymds2ju (s_year, s_month, s_day, s_sec, date_start)
    CALL forcing_printdate(date_start, "This is after reading the start date")
    !
    !Config Key  = END_DATE
    !Config Desc = Date at which the simulation ends
    !Config Def  = NONE
    !Config Help =  The format is the same as in the CF convention : 1999-09-13 12:0:0
    str_edate = " "
    CALL getin('END_DATE',str_edate)
    !
    IF ( (INDEX(str_edate(1),"-") .NE. INDEX(str_edate(1),"-", .TRUE.)) .AND. &
         &  (INDEX(str_edate(2),":") .NE. INDEX(str_edate(2),":", .TRUE.)) ) THEN
       DO i=1,2
          tmpstr = str_edate(1)
          ic = INDEX(tmpstr,"-")
          tmpstr(ic:ic) = " "
          str_edate(1) = tmpstr
          tmpstr = str_edate(2)
          ic = INDEX(tmpstr,":")
          tmpstr(ic:ic) = " "
          str_edate(2) = tmpstr
       ENDDO
       READ (str_edate(1),*) e_year, e_month, e_day
       READ (str_edate(2),*) hours, minutes, seci
       e_sec = hours*3600. + minutes*60. + seci
    ELSE
       CALL ipslerr(3, "forcing_integration_time", "END_DATE incorrectly specified in run.def", str_edate(1), str_edate(2))
    ENDIF
    CALL ymds2ju (e_year, e_month, e_day, e_sec, date_end)
    !
    CALL time_diff (s_year,s_month,s_day,s_sec,e_year,e_month,e_day,e_sec,diff_sec)
    !
    !Config Key  = DT_SECHIBA
    !Config Desc = Time step length in seconds for sechiba component
    !Config Def  = 1800
    !Config Help = 
    !Config Units = [seconds]
    dt = 1800
    CALL getin('DT_SECHIBA', dt)
    dt_sechiba_keep = dt
    !
    nbdt = NINT(diff_sec/dt)
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Read the configuration options for the time interpolations.
    !
    !Config Key   = LWDOWN_CONS
    !Config Desc  = Conserve the longwave downward radiation of the forcing
    !Config Def   = n
    !Config Help  = This flag allows to conserve the downward longwave radiation
    !               provided in the forcing. It will do this by taking the closest
    !               neighbour in time from the forcing. This assumes that the forcing
    !               contains average fluxes. The default setting (LWDOWN_CONS=n) will
    !               force the model to perform a linear interpolation of the fluxes.
    !Config Units = [FLAG]
    !-
    lwdown_cons = .FALSE.
    CALL getin('LWDOWN_CONS', lwdown_cons)
    !
  END SUBROUTINE forcing_integration_time
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_open
!!
!>\BRIEF      Opens the forcing files and extract the main information.
!!
!! DESCRIPTION:	 This routine opens all the forcing files provided in the list and verifies that the grid corresponds
!!               to the coordinates provided (and which was obtained by the model from glogrid.f90.). It then zooms
!!               into the forcing as requested by the user, extracts the vertical coordinates and final reads the time axis.
!!               Some basic consistency checks are performed as for instance ensuring the that all the forcing data is available
!!               to simulate the desired period.
!!               All that information is also broadcasted to all processors.
!!               Actual forcing data is not read at this stage.
!!
!! \n
!_ ==============================================================================================================================
!
  SUBROUTINE forcing_open(filenames_in, iim, jjm, lon, lat, nbpoint_in, drvzoom_lon, drvzoom_lat, &
       &                  kindex, nbindex_perproc, wunit, landonly_arg)
    !
    ! Opens the forcing file and reads some key information and stores them in the shared variables of the
    ! module.
    !
    ! Lon, lat should come from the grid file we read before. This will give indication of the grid
    ! file is consistant with the forcing file and if we need to zoom into the forcing file.
    !
    ! Time interval of the simulation is also determined.
    !
    ! ARGUMENTS
    !
    CHARACTER(LEN=*), INTENT(in)          :: filenames_in(:)
    INTEGER(i_std), INTENT(in)            :: iim, jjm, nbpoint_in
    REAL(r_std), INTENT(in)               :: lon(iim,jjm), lat(iim,jjm)
    REAL(r_std), DIMENSION(2), INTENT(in) :: drvzoom_lon, drvzoom_lat
    INTEGER(i_std), INTENT(in)            :: kindex(nbpoint_in)
    INTEGER(i_std), INTENT(in)            :: nbindex_perproc
    INTEGER(i_std), OPTIONAL, INTENT(in)  :: wunit
    LOGICAL, OPTIONAL, INTENT(in)         :: landonly_arg
    !
    ! LOCAL
    !
    INTEGER(i_std) :: iim_tmp, jjm_tmp, nbpoint_tmp, nb_files    
    INTEGER(i_std) :: iv, it
    INTEGER(i_std) :: inl, ii, jj, ik
    INTEGER(i_std) :: land_id
    REAL(r_std)    :: dt
    INTEGER(i_std) :: nbdt
    !
    ! Check optional arguments
    !
    ! The default behaviour is to provide only land points to the calling program.
    ! But for forcing ocean model there is also the option to pass on land and ocean values.
    ! When the grid is initialized landonly_tmp=.FALSE. has to be set to obtian this behaviour.
    !
    IF ( PRESENT(landonly_arg) ) THEN
       landonly=landonly_arg
    ELSE
       landonly = .TRUE.
    ENDIF
    !
    ! How many files do we have to open ?
    !
    !
    ! All the meta information from the forcing file is ojnly needed on the root processor.
    !
    IF ( is_root_prc ) THEN
       !
       CALL forcing_filenamecheck(filenames_in, nb_files)
       IF ( PRESENT(wunit) ) THEN
          DO it=1,nb_files
             WRITE(wunit,*) "Files to be used for forcing the simulation :", it, TRIM(forfilename(it))
          ENDDO
       ENDIF
       !
       ! 0.0 Check if variables are allocated to the right size on root_proc
       !
       IF (nb_files > nb_forcefile) THEN
          IF ( ALLOCATED(force_id) ) DEALLOCATE(force_id)
          ALLOCATE(force_id(nb_files))
          IF ( ALLOCATED(id_unlim) )  DEALLOCATE(id_unlim)
          ALLOCATE(id_unlim(nb_files))
          IF ( ALLOCATED(nb_atts) ) DEALLOCATE(nb_atts)
          ALLOCATE(nb_atts(nb_files))
          IF ( ALLOCATED(ndims) ) DEALLOCATE(ndims)
          ALLOCATE(ndims(nb_files))
          IF ( ALLOCATED(nvars) ) DEALLOCATE(nvars)
          ALLOCATE( nvars(nb_files))
          IF ( ALLOCATED(nbtime_perfile) ) DEALLOCATE(nbtime_perfile)
          ALLOCATE(nbtime_perfile(nb_files))
          IF ( ALLOCATED(convtosec) ) DEALLOCATE(convtosec)
          ALLOCATE(convtosec(nb_files))
       ENDIF
       nb_forcefile = nb_files
       !
       ! Get the global grid size from the forcing file. The output is in temporary variables as in this
       ! module the values are shared.
       !
       IF ( PRESENT(wunit) ) THEN
          WRITE(wunit,*) "Getting global grid from ",  nb_forcefile, "files."
          CALL FLUSH(wunit)
       ENDIF
       CALL forcing_getglogrid(nb_forcefile, forfilename, iim_tmp, jjm_tmp, nbpoint_tmp, .FALSE., landonly)
       !
       IF ( PRESENT(wunit) ) THEN
          WRITE(wunit,*) "Getting the zoomed grid", nbpoint_tmp
          CALL FLUSH(wunit)
       ENDIF
       CALL forcing_zoomgrid(drvzoom_lon, drvzoom_lat, forfilename(1), .FALSE.)
       IF ( PRESENT(wunit) ) THEN
          WRITE(wunit,*) "Out of the zoomed grid operation"
          CALL FLUSH(wunit)
       ENDIF
       !
       ! Verification that the grid sizes coming from the calling program are consistant with what we get 
       ! from the forcing file.
       !
       IF ( (iim_loc .NE. iim) .OR. (jjm_loc .NE. jjm) ) THEN
          CALL ipslerr (3,'forcing_open',"At least one of the dimensions of the grid obtained from the",&
               &        "grid file is different from the one in the forcing file.",&
               &        "Run driver2oasis -init to generate a new grid file.")
       ENDIF
       ! Special treatment for the number of land point, as we could have a case where the forcing
       ! file does not include the land/sea mask.
       !
       IF ( nbpoint_loc .NE. nbpoint_in ) THEN
          ! We trust the number of land points obtained from the gridfile. It has the land/sea mask.
          nbpoint_loc = nbpoint_in
       ENDIF
       !
       ! Treat the time dimension now :
       !
       IF ( PRESENT(wunit) ) THEN
          WRITE(wunit,*) "Getting forcing time"
          CALL FLUSH(wunit)
       ENDIF
       CALL forcing_time(nb_forcefile, forfilename)
       !
       ! Now that we know how much time steps are in the forcing we can set some realistic slab_size
       !
       slab_size=MIN(nb_forcing_steps, slab_size_max)
       !
       !
       ! Get the vertical information from the file
       !
       CALL forcing_vertical(force_id(1))
       !
       !
       IF ( PRESENT(wunit) ) THEN
          WRITE(wunit,*) "Getting integration time"
          CALL FLUSH(wunit)
       ENDIF
       CALL forcing_integration_time(startdate, dt, nbdt)
       !
       ! Test that the time interval requested by the user correspond to the time available in the 
       ! forcing file.
       !
       IF ( startdate < time_bounds(1,1,1) .OR. startdate > time_bounds(nb_forcing_steps,1,2) ) THEN
          CALL forcing_printdate(startdate, "--> Sarte Date in forcing_open")
          CALL forcing_printdate(time_bounds(1,1,1), "--> Outer bound of forcing file.")
          CALL forcing_printdate(time_bounds(nb_forcing_steps,1,2), "--> Last date to be simulated.")
          CALL ipslerr (3,'forcing_open', 'Start time requested by the user is outside of the time interval',&
               & "covered by the forcing file.","Please verify the configuration in the run.def file.")
       ENDIF
       !
       IF ( startdate+(dt/one_day)*nbdt > time_bounds(nb_forcing_steps,1,2) .OR. &
            & startdate+(dt/one_day)*nbdt < time_bounds(1,1,1)) THEN
          CALL forcing_printdate(time_bounds(nb_forcing_steps,1,2), "Outer bound of forcing file.")
          CALL forcing_printdate(startdate+(dt/one_day)*nbdt, "Last date to be simulated.")
          WRITE(*,*) "ERROR : Final date of forcing needed is : ", startdate+(dt/one_day)*nbdt
          WRITE(*,*) "ERROR : The outer bound of the last forcing time step is :", time_bounds(nb_forcing_steps,1,2)
          CALL ipslerr (3,'forcing_open', 'End time requested by the user is outside of the time interval',&
               & "covered by the forcing file.","Please verify the configuration in the run.def file.")
       ENDIF
       !
    ENDIF
    !
    ! Broadcast the local grid (i.e. the one resulting from the zoom) to all processors
    !
    CALL bcast(iim_loc)
    CALL bcast(jjm_loc)
    CALL bcast(nbpoint_loc)
    CALL bcast(nbland_loc)
    ! Time variables needed by all procs
    CALL bcast(slab_size)
    CALL bcast(startdate)
    CALL bcast(forcingstartdate)
    CALL bcast(forcing_tstep_ave)
    !
    ! Number of points per processor
    !
    IF ( landonly ) THEN
       nbpoint_proc = nbindex_perproc
    ELSE
       nbpoint_proc = nbpoint_glo
    ENDIF
    !
    ! On the slave processes we need to allocate the memory for the data on root_prc to be bcast
    ! On the root_proc these allocations were done with CALL forcing_zoomgrid
    !
    ALLOCATE(glolindex_proc(nbpoint_proc))
    IF ( .NOT. is_root_prc ) THEN
       ALLOCATE(lon_loc(iim_loc,jjm_loc))
       ALLOCATE(lat_loc(iim_loc,jjm_loc))
       ALLOCATE(lindex_loc(nbpoint_loc)) 
       ALLOCATE(mask_loc(iim_loc,jjm_loc))
       ALLOCATE(area_loc(iim_loc,jjm_loc))
       ALLOCATE(contfrac_loc(nbpoint_loc))
       ALLOCATE(corners_loc(iim_loc,jjm_loc,4,2))
    ENDIF
    !
    ! Keep on each processor the index of each land point on the *_loc grid
    !
    IF ( landonly ) THEN
       CALL scatter(kindex, glolindex_proc)
    ELSE
       !
       ! Build a simple indexing list as the one for land cannot be used.
       !
       ik=0
       DO jj=1,jjm_loc
          DO ii=1,iim_loc
             ik=ik+1
             glolindex_proc(ik) = ik
          ENDDO
       ENDDO
    ENDIF
    !
    CALL bcast(lon_loc)
    CALL bcast(lat_loc)
    CALL bcast(lindex_loc)
    CALL bcast(mask_loc)
    CALL bcast(area_loc)
    CALL bcast(contfrac_loc)
    CALL bcast(corners_loc)
    !
  END SUBROUTINE forcing_open
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_getvalues1d
!!
!>\BRIEF   Gets the forcing data for a time interval.   
!!
!! DESCRIPTION:	The routine will get the forcing valid for the time interval provided by the caller.
!!              First it will check that the data is already in memory for that time interval. If not
!!              it will first read the data from the netCDF file.
!!              Then the forcing date will be interpolated to the requested time interval.
!!              The code calls linear interpolation for most variables except for SWdown and precipitation.
!!              These temporal interpolations can be improved later.
!!
!! \n
!_ ==============================================================================================================================
  SUBROUTINE forcing_getvalues1d(time_int, dt, zlev_tq, zlev_uv, tair, qair, rainf, snowf, &
       &                       swdown, lwdown, solarang, u, v, ps)
    !
    ! ARGUMENTS
    !
    REAL(r_std), INTENT(in)  :: time_int(2)                            !! The time interval over which the forcing is needed.
    REAL(r_std), INTENT(in)  :: dt                                     !! timestep, i.e. distance in seconds between time_int(1) and time_int(2)
    REAL(r_std), INTENT(out) :: zlev_tq(:), zlev_uv(:)
    REAL(r_std), INTENT(out) :: tair(:), qair(:), rainf(:), snowf(:)
    REAL(r_std), INTENT(out) :: swdown(:), lwdown(:), solarang(:)
    REAL(r_std), INTENT(out) :: u(:), v(:), ps(:)
    !
    ! LOCAL
    !
    INTEGER(i_std) :: i
    !
    ! Test that we have the time interval within our slab of data else we need to update it.
    ! Att : the tests are done here on time_tair as an exemple. This might need to have to be generalized.
    !
    !
    ! First case the time axis of the variable are not even yet allocated !
    IF ( .NOT. ALLOCATED(time_tair) ) THEN
       CALL forcing_readslab(time_int)
       CALL forcing_printdate(timebnd_tair(1,1), "Start of time slab just read")
       CALL forcing_printdate(timebnd_tair(slab_size,2), "End of time slab just read")
    ELSE
       ! If we have time axis (for TAIR here) we test that it is long enough in time to allow for an interpolation.
       !
       IF ( time_int(2)+forcing_tstep_ave/one_day > time_tair(slab_size) .AND. (.NOT. end_of_file) ) THEN
          CALL forcing_readslab(time_int)
          CALL forcing_printdate(timebnd_tair(1,1), "Start of time slab just read")
          CALL forcing_printdate(timebnd_tair(slab_size,2), "End of time slab just read")
       ENDIF
    ENDIF
    !
    ! Interpolate the dynamical variables to the time step at which the driver is for the moment.
    !
    CALL forcing_interpol(time_int, dt, time_u, u_slab, u)
    CALL forcing_interpol(time_int, dt, time_v, v_slab, v)
    CALL forcing_interpol(time_int, dt, time_ps, ps_slab, ps)
    !
    ! Compute the height of the first level (a routine will be needed for that !)
    ! ATT : we assume that the time axis for the height of the scalar variable is the one of TAIR
    ! and for the height of wind is the same as U.
    CALL forcing_interpol(time_int, dt, time_tair, ztq_slab, zlev_tq)
    CALL forcing_interpol(time_int, dt, time_u, zuv_slab, zlev_uv)
     !
    ! Interpolate the state variables of the lower atmospheric level
    !
    CALL forcing_interpol(time_int, dt, time_tair, tair_slab, tair)
    CALL forcing_interpol(time_int, dt, time_qair, qair_slab, qair)
    !
    ! Spread the precipitation as requested by the user
    !
    CALL forcing_spreadprec(time_int, dt, timebnd_precip, time_precip, rainf, snowf)
    !
    ! Deal with the interpolate of the radiative fluxes.
    !
    CALL forcing_solarint(time_int, dt, timebnd_swdown, time_swdown, iim_loc, jjm_loc, lon_loc, lat_loc, swdown, solarang)
    !
    ! We have the option here to conserve LWdown by taking the closest point in the forcing.
    ! So no interpolation is done.
    !
    IF ( lwdown_cons ) THEN
       CALL forcing_closest(time_int, dt, time_lwdown, lwdown_slab, lwdown)
    ELSE
       CALL forcing_interpol(time_int, dt, time_lwdown, lwdown_slab, lwdown)
    ENDIF
    !
  END SUBROUTINE forcing_getvalues1d
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_getvalues2d
!!
!>\BRIEF   Gets the forcing data in 2D field for a time interval.   
!!
!! DESCRIPTION:	The routine will get the forcing valid for the time interval provided by the caller.
!!              First it will check that the data is already in memory for that time interval. If not
!!              it will first read the data from the netCDF file.
!!              Then the forcing date will be interpolated to the requested time interval.
!!              The code calls linear interpolation for most variables except for SWdown and precipitation.
!!              These temporal interpolations can be improved later.
!!
!! \n
!_ ==============================================================================================================================
  SUBROUTINE forcing_getvalues2d(time_int, dt, zlev_tq, zlev_uv, tair, qair, rainf, snowf, &
       &                       swdown, lwdown, solarang, u, v, ps)
    !
    ! ARGUMENTS
    !
    REAL(r_std), INTENT(in)  :: time_int(2)                            !! The time interval over which the forcing is needed.
    REAL(r_std), INTENT(in)  :: dt                                     !! timestep, i.e. distance in seconds between time_int(1) and time_int(2)
    REAL(r_std), INTENT(out) :: zlev_tq(:,:), zlev_uv(:,:)
    REAL(r_std), INTENT(out) :: tair(:,:), qair(:,:), rainf(:,:), snowf(:,:)
    REAL(r_std), INTENT(out) :: swdown(:,:), lwdown(:,:), solarang(:,:)
    REAL(r_std), INTENT(out) :: u(:,:), v(:,:), ps(:,:)
    !
    REAL(r_std) :: zzlev_tq(nbpoint_loc), zzlev_uv(nbpoint_loc)
    REAL(r_std) :: ztair(nbpoint_loc), zqair(nbpoint_loc), zrainf(nbpoint_loc), zsnowf(nbpoint_loc)
    REAL(r_std) :: zswdown(nbpoint_loc), zlwdown(nbpoint_loc), zsolarang(nbpoint_loc)
    REAL(r_std) :: zu(nbpoint_loc), zv(nbpoint_loc), zps(nbpoint_loc)
    INTEGER(i_std) :: i, j, k
    !
    CALL forcing_getvalues(time_int, dt, zzlev_tq, zzlev_uv, ztair, zqair, zrainf, zsnowf, zswdown, zlwdown, zsolarang, zu, zv, zps)
    !
    k = 0
    DO j=1,jjm_loc
       DO i=1,iim_loc
          k = k + 1
          zlev_tq(i,j) = zzlev_tq(k)
          zlev_uv(i,j) = zzlev_uv(k)
          tair(i,j) = ztair(k)
          qair(i,j) = zqair(k)
          rainf(i,j) = zrainf(k)
          snowf(i,j) = zsnowf(k)
          swdown(i,j) = zswdown(k)
          lwdown(i,j) = zlwdown(k)
          solarang(i,j) = zsolarang(k)
          u(i,j) = zu(k)
          v(i,j) = zv(k)
          ps(i,j) = zps(k)
       ENDDO
    ENDDO
    !
  END SUBROUTINE forcing_getvalues2d
    
!!  =============================================================================================================================
!! SUBROUTINE: forcing_closest
!!
!>\BRIEF   This routine does not interpolate and simply uses the closes value in time. It is useful for preserving
!!         variables which are averaged in the forcing file.
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================
!
!=============================================================================================
!
  SUBROUTINE forcing_closest(time_int_in, dt, time_central_in, var_slab, var)
    !
    ! ARGUMENTS
    !
    REAL(r_std), INTENT(in)  :: time_int_in(2)
    REAL(r_std), INTENT(in)  :: dt
    REAL(r_std), INTENT(in)  :: time_central_in(:)
    REAL(r_std), INTENT(in)  :: var_slab(:,:)
    REAL(r_std), INTENT(out) :: var(:)
    !
    ! LOCAL
    !
    INTEGER(i_std) :: slabind_a, slabind_b, imin(1), i
    REAL(r_std) :: time_int(2), time_central(slab_size_max)
    REAL(r_std) :: mid_int, wa, wb, wt, wab, wae, tmp_mid_int
    LOGICAL :: mask(slab_size_max)=.FALSE.
    !
    ! Shift the input dates in order to gain in precision for the calculations 
    !
    time_int(:) = time_int_in(:)-INT(forcingstartdate)
    time_central(1:slab_size) = time_central_in(1:slab_size)-INT(forcingstartdate)
    !
    ! Create a mask so that MINLOC does not look outside of the valid interval of time_central
    !
    mask(1:slab_size) = .TRUE.
    !
    ! Select the forcing interval for which the center date is the closest to the time of 
    ! the model.
    !
    mid_int = time_int(1) + (dt/2.0)/one_day
    imin = MINLOC( ABS(time_central(1:slab_size) - mid_int), mask )
    !
    ! Verify that this is a possible date
    !
    IF ( imin(1) > 0 .AND. imin(1) <= slab_size ) THEN
       !
       slabind_a = imin(1)
       !
    ELSE
       WRITE(*,*) "imin(1) = ", imin(1), (time_int_in(1) + (dt/2.0)/one_day)
       CALL forcing_printdate(time_int_in(1), "===> Start of target time interval.")
       CALL forcing_printdate(time_int_in(2), "===> End of target time interval.")
       CALL forcing_printdate(time_central_in(imin(1)), "===> Center of forcing time interval.")
       CALL ipslerr (3,'forcing_closest', 'The target time interval has no acceptable closest',&
            & "time in the forcing slab.","")
    ENDIF
    !
    ! Transfer the data from the sloest time of the forcing data slab.
    !
    DO i=1, nbpoint_proc
       !
       var(i) = var_slab(i,slabind_a)
       !
    ENDDO
    !
    !
  END SUBROUTINE forcing_closest
  
!!  =============================================================================================================================
!! SUBROUTINE: forcing_interpol
!!
!>\BRIEF   Perform linear interpolation for the time interval requested.
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================
!
!=============================================================================================
!
  SUBROUTINE forcing_interpol(time_int_in, dt, time_central_in, var_slab, var)
    !
    ! ARGUMENTS
    !
    REAL(r_std), INTENT(in)  :: time_int_in(2)
    REAL(r_std), INTENT(in)  :: dt
    REAL(r_std), INTENT(in)  :: time_central_in(:)
    REAL(r_std), INTENT(in)  :: var_slab(:,:)
    REAL(r_std), INTENT(out) :: var(:)
    !
    ! LOCAL
    !
    INTEGER(i_std) :: slabind_a, slabind_b, imin(1), i
    REAL(r_std) :: time_int(2), time_central(slab_size_max)
    REAL(r_std) :: mid_int, wa, wb, wt, wab, wae, tmp_mid_int
    LOGICAL :: mask(slab_size_max)=.FALSE.
    !
    ! Shift the input dates in order to gain in precision for the calculations 
    !
    time_int(:) = time_int_in(:)-INT(forcingstartdate)
    time_central(1:slab_size) = time_central_in(1:slab_size)-INT(forcingstartdate)
    !
    ! Create a mask so that MINLOC does not look outside of the valid interval of time_central
    !
    mask(1:slab_size) = .TRUE.
    !
    ! Select the type of interpolation to be done.
    !
    mid_int = time_int(1) + (dt/2.0)/one_day
    imin = MINLOC( ABS(time_central(1:slab_size) - mid_int), mask )
    !
    IF ( imin(1) > 1 .AND. imin(1) < slab_size ) THEN
       !
       IF ( mid_int < time_central(imin(1)) ) THEN
          slabind_a = imin(1) - 1
          slabind_b = imin(1)
       ELSE
          slabind_a = imin(1)
          slabind_b = imin(1) + 1
       ENDIF
       !
    ELSE IF ( imin(1) == 1 ) THEN
       slabind_a = 1
       slabind_b = 2
       IF ( mid_int < time_central(slabind_a) ) THEN
          IF ( time_int(2) < time_central(slabind_a) ) THEN
             WRITE(*,*) "imin(1) = ", imin(1), (time_int_in(1) + (dt/2.0)/one_day)
             CALL forcing_printdate(time_int_in(1), "===> Start of target time interval.")
             CALL forcing_printdate(time_int_in(2), "===> End of target time interval.")
             CALL forcing_printdate(time_central_in(slabind_a), "===> Center of forcing time interval.")
             CALL ipslerr (3,'forcing_interpol', 'The target time interval lies before the first date of the slab.',&
                  & "","")
          ELSE
             mid_int = time_central(slabind_a) 
          ENDIF
       ENDIF
    ELSE IF ( imin(1) == slab_size ) THEN
       slabind_a = slab_size - 1
       slabind_b = slab_size
       IF ( mid_int > time_central(slabind_b) ) THEN
          IF ( time_int(1) > time_central(slabind_b) ) THEN
             WRITE(*,*) "imin(1) = ", imin(1), (time_int_in(1) + (dt/2.0)/one_day)
             CALL forcing_printdate(time_int_in(1), "===> Start of target time interval.")
             CALL forcing_printdate(time_int_in(2), "===> End of target time interval.")
             CALL forcing_printdate(time_central_in(slabind_b), "===> Center of forcing time interval.")
             CALL ipslerr (3,'forcing_interpol', 'The target time interval lies after the last date of the slab.',&
                  & "","")
          ELSE
             mid_int = time_central(slabind_b) 
          ENDIF
       ENDIF
    ENDIF
    !
    ! Compute the weights between the values at slabind_a and slabind_b. As with the time
    ! representation we are at the limit of precision we use 2 days to compute the distance
    ! in time between the first value (slabind_a) and the middle of the target interval.
    !
    wab = time_int(1) - time_central(slabind_a) + (dt/2.0)/one_day
    wae = time_int(2) - time_central(slabind_a) - (dt/2.0)/one_day
    wa = (wab+wae)/2.0
    wb = time_central(slabind_b) - time_central(slabind_a)
    wt = wa/wb
    !
    ! Do the weighted average of all land points with the time indices and weights computed above.
    !
    DO i=1, nbpoint_proc
       var(i) = var_slab(i,slabind_a) + wt*(var_slab(i,slabind_b) - var_slab(i,slabind_a))
       !
    ENDDO
    !
    !
  END SUBROUTINE forcing_interpol
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_spreadprec
!!
!>\BRIEF      Spreads the precipitation over the interval chosen based on the interval chosen by the user.
!!
!! DESCRIPTION:	The behaviour of this routine is controlled by the parameter SPRED_PREC_SEC in the run.def.
!!              The time in second specified by the user will be the one over which the precipitation will last
!!              where the forcing interval has rain or snow.
!!
!! \n
!_ ==============================================================================================================================
!
!=============================================================================================
!
  SUBROUTINE forcing_spreadprec(time_int, tlen, timebnd_central, time_central, rainf, snowf)
    !
    ! ARGUMENTS
    !
    REAL(r_std), INTENT(in)  :: time_int(2)         ! Time interval to which we will spread precip
    REAL(r_std), INTENT(in)  :: tlen                ! size of time interval in seconds (time step !)
    REAL(r_std), INTENT(in)  :: timebnd_central(:,:)    ! Time interval over which the read data is valid
    REAL(r_std), INTENT(in)  :: time_central(:)     ! Center of the time interval
    REAL(r_std), INTENT(out) :: rainf(:), snowf(:)
    !
    ! LOCAL
    !
    LOGICAL, SAVE :: first_call_spreadprec=.TRUE.
    REAL(r_std), SAVE :: time_to_spread=3600.0
    INTEGER(i_std) :: imin(1), i, tind(3)
    REAL(r_std) :: ft(3), dt, left, right
    INTEGER(i_std) :: offset, nb_spread
    LOGICAL :: mask(slab_size_max)=.FALSE.
    !
    IF ( first_call_spreadprec ) THEN
       !Config Key   = SPRED_PREC
       !Config Desc  = Spread the precipitation.
       !Config If    = [-]
       !Config Def   = Half of the forcing time step or uniform, depending on dt_force and dt_sechiba
       !Config Help  = Spread the precipitation over SPRED_PREC steps of the splited forcing 
       !Config         time step. This ONLY applied if the forcing time step has been splited.
       !Config         If the value indicated is greater than SPLIT_DT, SPLIT_DT is used for it.
       !Config Units = [-]
       !-
       nb_spread = -1
       CALL getin_p('SPRED_PREC', nb_spread)
       !
       ! Test if we have read the number of time steps to spread in run.def
       ! If not, then probably the time was given in seconds.
       !
       IF ( nb_spread < 0 ) THEN
          !Config Key   = SPRED_PREC_SEC
          !Config Desc  = Spread the precipitation over an interval in seconds.
          !Config Def   = 3600
          !Config Help  = Spread the precipitation over n seconds of the forcing time step
          !Config         interval. This ONLY applies when the SPRED_PREC_SEC is smaller than
          !Config         the forcing time step. Should the user set SPRED_PREC_SEC=0 we will 
          !Config         assume that the rainfall is uniformely distributed over the forcing interval.
          !Config Units = seconds
          !
          ! This is the default should 'SPRED_PREC' not be present in the run.def
          !
          time_to_spread = forcing_tstep_ave/2.0
          !
          CALL getin_p('SPRED_PREC_SEC', time_to_spread)
       ELSE
          time_to_spread = dt_sechiba_keep * nb_spread
       ENDIF
       !
       ! Do some verifications on the information read from run.def
       !
       IF ( time_to_spread > forcing_tstep_ave) THEN
          time_to_spread = forcing_tstep_ave
       ELSE IF ( time_to_spread <= 0 ) THEN
          time_to_spread = forcing_tstep_ave
       ENDIF
       !
       first_call_spreadprec = .FALSE.
       !
    ENDIF
    !
    ! First test that we have the right time interval from the forcing to spread the precipitation
    !
    IF ( time_int(1) >= timebnd_central(1,1) .AND. time_int(2) <= timebnd_central(slab_size,2)) THEN
       !
       ! Create a mask so that MINLOC does not look outside of the valid interval of time_central
       !
       mask(1:slab_size) = .TRUE.
       !
       ! To get better precision on the time difference we get a common offset to substract 
       !
       offset = INT(forcingstartdate)
       !
       ! In principle 3 time steps can contribute to the time step closest to the center of the forcing interval
       !
       imin = MINLOC( ABS(time_central(1:slab_size)-(time_int(1)+time_int(2))/2.0), mask )
       tind(1) = MAX(imin(1)-1,1)
       tind(2) = imin(1)
       tind(3) = MIN(imin(1)+1,slab_size)
       IF (imin(1)+1 > slab_size) THEN
          WRITE(*,*) "We have a problem here imin(1)+1,slab_size ", imin(1)+1,slab_size
          WRITE(*,*) "Interval : ", time_int(1),time_int(2)
       ENDIF
       !
       !
       !
       ! Do we need to take some rain from the previous time step ?
       !
       !! Time computation is not better than 1/1000 seconds
       IF ( time_int(1) < timebnd_central(tind(2),1) .AND. preciptime_slab(tind(1)) < (time_to_spread-0.001) ) THEN
          dt = ((timebnd_central(tind(2),1)-offset)-(time_int(1)-offset))*one_day
          ft(1) = MIN(time_to_spread - preciptime_slab(tind(1)), dt)/tlen
       ELSE
          ft(1) = zero
       ENDIF
       !
       ! Is there still some rain to spread from the current forcing time step ?
       !
       !! Time computation is not better than 1/1000 seconds
       IF (preciptime_slab(tind(2)) < (time_to_spread-0.001) ) THEN
          left = MAX(time_int(1), timebnd_central(tind(2),1))
          right = MIN(time_int(2),timebnd_central(tind(2),2))
          dt = ((right-offset)-(left-offset))*one_day
          ft(2) = MIN(time_to_spread - preciptime_slab(tind(2)), dt)/tlen
       ELSE
          ft(2) = zero
       ENDIF
       !
       ! Do we need to take some rain from the next time step ?
       !
       !! Time computation is not better than 1/1000 seconds
       IF ( time_int(2) > timebnd_central(tind(2),2) .AND. preciptime_slab(tind(3)) < (time_to_spread-0.001) ) THEN
          dt = ((time_int(2)-offset)-(timebnd_central(tind(2),2)-offset))*one_day
          ft(3) = MIN(time_to_spread - preciptime_slab(tind(3)), dt)/tlen
       ELSE
          ft(3) = zero
       ENDIF
       !
       ! Do the actual calculation
       !
       DO i=1, nbpoint_proc
          rainf(i) = (rainf_slab(i,tind(1)) * forcing_tstep_ave * ft(1) + &
                  &  rainf_slab(i,tind(2)) * forcing_tstep_ave * ft(2) + &
                  &  rainf_slab(i,tind(3)) * forcing_tstep_ave * ft(3))*tlen/time_to_spread
  
          snowf(i) = (snowf_slab(i,tind(1)) * forcing_tstep_ave * ft(1) + &
                  &  snowf_slab(i,tind(2)) * forcing_tstep_ave * ft(2) + &
                  &  snowf_slab(i,tind(3)) * forcing_tstep_ave * ft(3))*tlen/time_to_spread
       ENDDO
       !
       ! Update the time over which we have already spread the rainf
       !
       preciptime_slab(tind(1)) = preciptime_slab(tind(1)) + tlen * ft(1)
       preciptime_slab(tind(2)) = preciptime_slab(tind(2)) + tlen * ft(2)
       preciptime_slab(tind(3)) = preciptime_slab(tind(3)) + tlen * ft(3)
       !
    ELSE
       WRITE(numout,*) "Time interval toward which we will interpolate : ", time_int
       WRITE(numout,*) "Limits of the time slab we have : ", timebnd_central(1,1), timebnd_central(slab_size,2)
       CALL forcing_printdate(time_int(1), "Start of target time interval.")
       CALL forcing_printdate(time_int(2), "End of target time interval.")
       CALL forcing_printdate(timebnd_central(1,1), "Start of time slab we have.")
       CALL forcing_printdate(timebnd_central(slab_size,2), "End of time slab we have.")
       CALL ipslerr (3,'forcing_spreadprec', 'The sitation should not occur Why are we here ?',&
            & "","")
    ENDIF
    !
  END SUBROUTINE forcing_spreadprec
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_solarint
!!
!>\BRIEF      Interpolates incoming solar radiation to the interval requested.
!!
!! DESCRIPTION:	The interpolation here takes into account the variation of the solar zenith angle
!!              to ensure the diurnal cycle of solar radiation is as well represented as possible.
!!
!! \n
!_ ==============================================================================================================================
!
!=============================================================================================
!
  SUBROUTINE forcing_solarint(time_int_in, tlen, timebnd_in, time_cent_in, iim, jjm, lon, lat, swdown, solarangle)
    !
    ! ARGUMENTS
    !
    REAL(r_std), INTENT(in)    :: time_int_in(2)             ! Time interval for which we will compute radiation
    REAL(r_std), INTENT(in)    :: tlen                       ! size of time interval in seconds (time step !)
    REAL(r_std), INTENT(in)    :: timebnd_in(:,:)            ! Time interval over which the read data is valid
    REAL(r_std), INTENT(in)    :: time_cent_in(:)            ! Center of the time interval
    INTEGER(i_std), INTENT(in) :: iim, jjm                   ! Size of 2D domain
    REAL(r_std), INTENT(in)    :: lon(iim,jjm), lat(iim,jjm) ! Longitude and latitude
    REAL(r_std), INTENT(out)   :: swdown(:), solarangle(:)   ! interpolated downward solar radiation and corresponding
    !                                                        ! solar angle.
    !
    ! LOCAL SAVED
    !
    LOGICAL, SAVE        :: first_call_solarint=.TRUE.
    REAL(r_std), SAVE    :: solaryearstart
    INTEGER(i_std), SAVE :: split, split_max
    REAL(r_std), SAVE    :: last_time
    !
    REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)   :: mean_sinang
    REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: sinangles
    REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)     :: time_angles
    ! Dusk-dawn management
    REAL(r_std), SAVE   :: dusk_angle
    !
    ! LOCAL - temporary
    !
    REAL(r_std) :: time_int(2)
    REAL(r_std) :: timebnd(slab_size_max,2)
    REAL(r_std) :: time_cent(slab_size_max)
    INTEGER(i_std) :: year, month, day, hours, minutes
    REAL(r_std) :: sec
    !
    REAL(r_std) :: mean_sol, split_time
    REAL(r_std) :: julian, julian_tmp
    REAL(r_std) :: sinang(iim,jjm)
    INTEGER(i_std) :: is, i, ii, jj, imin(1), tmin(1)
    LOGICAL :: mask(slab_size_max)=.FALSE.
    !
    IF ( first_call_solarint ) THEN
       !
       ! Ensure the offset is on the 1st of Januray of the current years so that we do not
       ! perturbe the solar angle calculations.
       !
       CALL ju2ymds (startdate, year, month, day, sec)
       CALL ymds2ju (year, 1, 1, 0.0, solaryearstart)
       !
       last_time = -9999.0
       !
       ALLOCATE(mean_sinang(iim,jjm))
       mean_sinang(:,:) = 0.0
       !
       split = NINT(forcing_tstep_ave/tlen)
       !
       ! Allow for more space than estimated with the size of the first time step.
       !
       ALLOCATE(time_angles(split*2))
       time_angles(:) = 0.0
       !
       ALLOCATE(sinangles(iim,jjm,split*2))
       sinangles(:,:,:) = 0.0
       !
       dusk_angle=0.01
       !
       first_call_solarint = .FALSE.
       !
       split = 0
       !
    ENDIF
    !
    ! Shift the input dates in order to gain in precision for the time calculations 
    !
    time_int(:) = time_int_in(:)-INT(solaryearstart)
    time_cent(1:slab_size) = time_cent_in(1:slab_size)-INT(solaryearstart)
    timebnd(1:slab_size,1) = timebnd_in(1:slab_size,1)-INT(solaryearstart)
    timebnd(1:slab_size,2) = timebnd_in(1:slab_size,2)-INT(solaryearstart)
    !
    ! Create a mask so that MINLOC does not look outside of the valid interval of time_central
    !
    mask(1:slab_size) = .TRUE.
    !
    ! Locate the time step in the SLAB at hand
    !
    imin = MINLOC( ABS(time_cent(1:slab_size)-(time_int(1)+time_int(2))/2.0), mask )
    !
    ! Compute all the angels we will encounter for the current forcing interval
    !
    IF ( last_time .NE. timebnd(imin(1),1) ) THEN
       !
       ! Verify that we have used all the angles of the previous decomposition of the forcing
       ! time step.
       !
       IF ( split .NE. 0 ) THEN
          !
          WRITE(numout,*) "The forcing has a time step of : ", forcing_tstep_ave
          WRITE(numout,*) "The model is configured to run with a time step of : ", tlen
          WRITE(numout,*) "We are left with split = ", split, " starting from ", split_max
          !
          CALL ipslerr (3,'forcing_solarint',"The decomposition of solar downward radiation of the forcing file over the model",&
               &        "has failed. This means the average of the solar radiation over the forcing time step is not conserved.",&
               &        "This can be caused by a time step repeated twice.")
       ENDIF
       !
       ! Compute the number of time steps the model will put in the current interval of forcing data. 
       !
       split = 0
       julian_tmp = (time_int(1)+time_int(2))/2.0
       split_time = julian_tmp+split*tlen/one_day
       tmin = MINLOC( ABS(time_cent(1:slab_size) - split_time), mask)
       DO WHILE (  tmin(1) .EQ. imin(1) .AND. split_time .LE. timebnd(slab_size,2) )
          split = split + 1
          split_time = julian_tmp+split*tlen/one_day
          tmin = MINLOC( ABS(time_cent(1:slab_size) - split_time), mask)
       ENDDO
       !
       mean_sinang(:,:) = 0.0
       time_angles(:) = 0.0
       !
       DO is = 1,split
          !
          julian = julian_tmp + (is-1)*tlen/one_day
          !
          ! This call should be better at it allows to compute the difference between the
          ! current date and a reference time to higher precision. But it produces noisy
          ! SWdown fluxes !
!!          CALL solarang (julian, solaryearstart, iim, jjm, lon, lat, sinang)
          ! The old approach.
          CALL solarang (julian+INT(solaryearstart), solaryearstart, iim, jjm, lon, lat, sinang)
          !
          ! During the dusk,dawn period maintain a minimum angle to take into account the
          ! diffuse radiation which starts before the sun is over the horizon.
          !
          DO ii=1,iim
             DO jj=1,jjm
                IF ( sinang(ii,jj) > zero .AND.  sinang(ii,jj) < dusk_angle ) THEN
                   sinang(ii,jj) = dusk_angle
                ENDIF
                mean_sinang(ii,jj) = mean_sinang(ii,jj)+sinang(ii,jj)
             ENDDO
          ENDDO
          !
          ! Save the solar angle information for later use. That is when we have the target time we will
          ! look in this table the angle we have forseen.
          !
          time_angles(is) = julian
          sinangles(:,:,is) = sinang(:,:)
          !
       ENDDO
       !
       mean_sinang(:,:) = mean_sinang(:,:)/split
       last_time =  timebnd(imin(1),1)
       split_max = split
       !
    ENDIF
    !
    ! For the current timle step get the time of the closest pre-computed solar angle.
    !
    julian = (time_int(1)+time_int(2))/2.0
    tmin =  MINLOC(ABS(julian-time_angles(1:split_max)), mask)
    sinang(:,:) = sinangles(:,:,tmin(1))
    ! Remember that we have taken one value of the table for later verification
    split = split - 1
    !
    DO i=1, nbpoint_proc
       !
       jj = ((glolindex_proc(i)-1)/iim)+1
       ii = (glolindex_proc(i)-(jj-1)*iim)
       !
       IF ( mean_sinang(ii,jj) > zero ) THEN
          swdown(i) = swdown_slab(i,imin(1))*sinang(ii,jj)/mean_sinang(ii,jj)
       ELSE
          swdown(i) = zero
       ENDIF
       !
       ! Why is this ??? Should we not take the angle corresponding to this time step ?? (Jan)
       !
       solarangle(i) = mean_sinang(ii,jj)
       !
    ENDDO
    !
  END SUBROUTINE forcing_solarint
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_readslab
!!
!>\BRIEF Interface routine to read the data. This routine prepares the memory on each procesor and scatters the read data.
!!
!! DESCRIPTION:
!!
!! \n
!_ ==============================================================================================================================
!
!=============================================================================================
!
  SUBROUTINE forcing_readslab(time_int)
    !
    ! This routine serves to interface with forcing_readslab_root and ensure that
    ! the data is distributed correctly on all processors.
    !
    REAL(r_std), INTENT(in)  :: time_int(2)                            !! The time interval over which the forcing is needed.
    !
    ! Local
    !
    INTEGER(i_std)  :: is
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)    :: tair_full, qair_full
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)    :: rainf_full, snowf_full
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)    :: swdown_full, lwdown_full
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)    :: u_full, v_full
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)    :: ps_full, ztq_full, zuv_full
    !
    ! 1.0 Verify that for the slabs the memory is allocated for the variable
    ! as well as its time axis.
    !
    IF ( .NOT. ALLOCATED(tair_slab) ) ALLOCATE(tair_slab(nbpoint_proc,slab_size))
    IF ( .NOT. ALLOCATED(time_tair) ) ALLOCATE(time_tair(slab_size))
    IF ( .NOT. ALLOCATED(timebnd_tair) ) ALLOCATE(timebnd_tair(slab_size,2))
    !
    IF ( .NOT. ALLOCATED(qair_slab) ) ALLOCATE(qair_slab(nbpoint_proc,slab_size))
    IF ( .NOT. ALLOCATED(time_qair) ) ALLOCATE(time_qair(slab_size))
    IF ( .NOT. ALLOCATED(timebnd_qair) ) ALLOCATE(timebnd_qair(slab_size,2))
    !
    IF ( .NOT. ALLOCATED(rainf_slab) ) ALLOCATE(rainf_slab(nbpoint_proc,slab_size))
    IF ( .NOT. ALLOCATED(snowf_slab) ) ALLOCATE(snowf_slab(nbpoint_proc,slab_size))
    IF ( .NOT. ALLOCATED(time_precip) ) ALLOCATE(time_precip(slab_size))
    IF ( .NOT. ALLOCATED(timebnd_precip) ) ALLOCATE(timebnd_precip(slab_size,2))
    IF ( .NOT. ALLOCATED(preciptime_slab) ) ALLOCATE(preciptime_slab(slab_size))
    !
    IF ( .NOT. ALLOCATED(swdown_slab) ) ALLOCATE(swdown_slab(nbpoint_proc,slab_size))
    IF ( .NOT. ALLOCATED(time_swdown) ) ALLOCATE(time_swdown(slab_size))
    IF ( .NOT. ALLOCATED(timebnd_swdown) ) ALLOCATE(timebnd_swdown(slab_size,2))
    !
    IF ( .NOT. ALLOCATED(lwdown_slab) ) ALLOCATE(lwdown_slab(nbpoint_proc,slab_size))
    IF ( .NOT. ALLOCATED(time_lwdown) ) ALLOCATE(time_lwdown(slab_size))
    IF ( .NOT. ALLOCATED(timebnd_lwdown) ) ALLOCATE(timebnd_lwdown(slab_size,2))
    !
    IF ( .NOT. ALLOCATED(u_slab) ) ALLOCATE(u_slab(nbpoint_proc,slab_size))
    IF ( .NOT. ALLOCATED(time_u) ) ALLOCATE(time_u(slab_size))
    IF ( .NOT. ALLOCATED(timebnd_u) ) ALLOCATE(timebnd_u(slab_size,2))
    !
    IF ( .NOT. ALLOCATED(v_slab) ) ALLOCATE(v_slab(nbpoint_proc,slab_size))
    IF ( .NOT. ALLOCATED(time_v) ) ALLOCATE(time_v(slab_size))
    IF ( .NOT. ALLOCATED(timebnd_v) ) ALLOCATE(timebnd_v(slab_size,2))
    !
    IF ( .NOT. ALLOCATED(ps_slab) ) ALLOCATE(ps_slab(nbpoint_proc,slab_size))
    IF ( .NOT. ALLOCATED(time_ps) ) ALLOCATE(time_ps(slab_size))
    IF ( .NOT. ALLOCATED(timebnd_ps) ) ALLOCATE(timebnd_ps(slab_size,2))
    !
    IF ( .NOT. ALLOCATED(ztq_slab) ) ALLOCATE(ztq_slab(nbpoint_proc,slab_size))
    IF ( .NOT. ALLOCATED(zuv_slab) ) ALLOCATE(zuv_slab(nbpoint_proc,slab_size))
    !    
    !
    IF ( is_root_prc) THEN
       !
       ! Allocate ont he root processor the memory to receive the data from the file
       !
       ALLOCATE(tair_full(nbpoint_loc,slab_size))
       ALLOCATE(qair_full(nbpoint_loc,slab_size))
       ALLOCATE(rainf_full(nbpoint_loc,slab_size))
       ALLOCATE(snowf_full(nbpoint_loc,slab_size))
       ALLOCATE(swdown_full(nbpoint_loc,slab_size))
       ALLOCATE(lwdown_full(nbpoint_loc,slab_size))
       ALLOCATE(u_full(nbpoint_loc,slab_size))
       ALLOCATE(v_full(nbpoint_loc,slab_size))
       ALLOCATE(ps_full(nbpoint_loc,slab_size))
       ALLOCATE(ztq_full(nbpoint_loc,slab_size))
       ALLOCATE(zuv_full(nbpoint_loc,slab_size))
       !
       CALL forcing_readslab_root(time_int, &
            &                     tair_full, time_tair, timebnd_tair, &
            &                     qair_full, time_qair, timebnd_qair, &
            &                     rainf_full, snowf_full, time_precip, timebnd_precip, preciptime_slab, &
            &                     swdown_full, time_swdown, timebnd_swdown, &
            &                     lwdown_full, time_lwdown, timebnd_lwdown, &
            &                     u_full, time_u, timebnd_u, &
            &                     v_full, time_v, timebnd_v, &
            &                     ps_full, time_ps, timebnd_ps, &
            &                     ztq_full, zuv_full)
       !
    ELSE
       !
       ALLOCATE(tair_full(1,1))
       ALLOCATE(qair_full(1,1))
       ALLOCATE(rainf_full(1,1))
       ALLOCATE(snowf_full(1,1))
       ALLOCATE(swdown_full(1,1))
       ALLOCATE(lwdown_full(1,1))
       ALLOCATE(u_full(1,1))
       ALLOCATE(v_full(1,1))
       ALLOCATE(ps_full(1,1))
       ALLOCATE(ztq_full(1,1))
       ALLOCATE(zuv_full(1,1))
       !
    ENDIF
    !
    ! Broadcast the time information to all procs.
    !
    CALL bcast(slab_size)
    CALL bcast(time_tair)
    CALL bcast(timebnd_tair)
    CALL bcast(time_qair)
    CALL bcast(timebnd_qair)
    CALL bcast(time_precip)
    CALL bcast(timebnd_precip)
    CALL bcast(preciptime_slab)
    CALL bcast(time_swdown)
    CALL bcast(timebnd_swdown)
    CALL bcast(time_lwdown)
    CALL bcast(timebnd_lwdown)
    CALL bcast(time_u)
    CALL bcast(timebnd_u)
    CALL bcast(time_v)
    CALL bcast(timebnd_v)
    CALL bcast(time_ps)
    CALL bcast(timebnd_ps)
    !
    ! Scatter the slabs of data to all processors
    !
    IF ( landonly ) THEN
       CALL scatter(tair_full, tair_slab)
       CALL scatter(qair_full, qair_slab)
       CALL scatter(rainf_full, rainf_slab)
       CALL scatter(snowf_full, snowf_slab)
       CALL scatter(swdown_full, swdown_slab)
       CALL scatter(lwdown_full, lwdown_slab)
       CALL scatter(u_full, u_slab)
       CALL scatter(v_full, v_slab)
       CALL scatter(ps_full, ps_slab)
       CALL scatter(ztq_full, ztq_slab)
       CALL scatter(zuv_full, zuv_slab)
    ELSE
       tair_slab(:,:) = tair_full(:,:)
       qair_slab(:,:) = qair_full(:,:)
       rainf_slab(:,:) = rainf_full(:,:)
       snowf_slab(:,:) = snowf_full(:,:)
       swdown_slab(:,:) = swdown_full(:,:)
       lwdown_slab(:,:) = lwdown_full(:,:)
       u_slab(:,:) = u_full(:,:)
       v_slab(:,:) = v_full(:,:)
       ps_slab(:,:) = ps_full(:,:)
       ztq_slab(:,:) = ztq_full(:,:)
       zuv_slab(:,:) = zuv_full(:,:)
    ENDIF
    !
    ! Clean-up to free the memory on the root processor.
    !
    IF ( ALLOCATED(tair_full) ) DEALLOCATE(tair_full)
    IF ( ALLOCATED(qair_full) ) DEALLOCATE(qair_full)
    IF ( ALLOCATED(rainf_full) ) DEALLOCATE(rainf_full)
    IF ( ALLOCATED(snowf_full) ) DEALLOCATE(snowf_full)
    IF ( ALLOCATED(swdown_full) ) DEALLOCATE(swdown_full)
    IF ( ALLOCATED(lwdown_full) ) DEALLOCATE(lwdown_full)
    IF ( ALLOCATED(u_full) ) DEALLOCATE(u_full)
    IF ( ALLOCATED(v_full) ) DEALLOCATE(v_full)
    IF ( ALLOCATED(ps_full) ) DEALLOCATE(ps_full)
    IF ( ALLOCATED(ztq_full) ) DEALLOCATE(ztq_full)
    IF ( ALLOCATED(zuv_full) ) DEALLOCATE(zuv_full)
    !
  END SUBROUTINE forcing_readslab
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_readslab_root
!!
!>\BRIEF Routine which reads a slab of data from the netCDF file and writes it onto the memory.
!!
!! DESCRIPTION:	It is important to read the next slab of data while still keeping an overlap so that
!!              interpolation can continue.
!!              It also attributes a time axis to each variable.
!!
!! \n
!_ ==============================================================================================================================
!
!=============================================================================================
!
  SUBROUTINE forcing_readslab_root(time_int, &
            &                     tair, t_tair, tbnd_tair, &
            &                     qair, t_qair, tbnd_qair, &
            &                     rainf, snowf, t_prec, tbnd_prec, prectime, &
            &                     swdown, t_swdown, tbnd_swdown, &
            &                     lwdown, t_lwdown, tbnd_lwdown, &
            &                     u, t_u, tbnd_u, &
            &                     v, t_v, tbnd_v, &
            &                     ps, t_ps, tbnd_ps, &
            &                     ztq, zuv)
    !
    ! Arguments
    !
    REAL(r_std), INTENT(in)  :: time_int(2)                            !! The time interval over which the forcing is needed.
    !
    REAL(r_std), INTENT(out) :: tair(:,:)
    REAL(r_std), INTENT(out) :: t_tair(:)
    REAL(r_std), INTENT(out) :: tbnd_tair(:,:)
    !
    REAL(r_std), INTENT(out) :: qair(:,:)
    REAL(r_std), INTENT(out) :: t_qair(:)
    REAL(r_std), INTENT(out) :: tbnd_qair(:,:)
    !
    REAL(r_std), INTENT(out) :: rainf(:,:)
    REAL(r_std), INTENT(out) :: snowf(:,:)
    REAL(r_std), INTENT(out) :: t_prec(:)
    REAL(r_std), INTENT(out) :: tbnd_prec(:,:)
    REAL(r_std), INTENT(out) :: prectime(:)
    !
    REAL(r_std), INTENT(out) :: swdown(:,:)
    REAL(r_std), INTENT(out) :: t_swdown(:)
    REAL(r_std), INTENT(out) :: tbnd_swdown(:,:)
    !
    REAL(r_std), INTENT(out) :: lwdown(:,:)
    REAL(r_std), INTENT(out) :: t_lwdown(:)
    REAL(r_std), INTENT(out) :: tbnd_lwdown(:,:)
    !
    REAL(r_std), INTENT(out) :: u(:,:)
    REAL(r_std), INTENT(out) :: t_u(:)
    REAL(r_std), INTENT(out) :: tbnd_u(:,:)
    !
    REAL(r_std), INTENT(out) :: v(:,:)
    REAL(r_std), INTENT(out) :: t_v(:)
    REAL(r_std), INTENT(out) :: tbnd_v(:,:)
    !
    REAL(r_std), INTENT(out) :: ps(:,:)
    REAL(r_std), INTENT(out) :: t_ps(:)
    REAL(r_std), INTENT(out) :: tbnd_ps(:,:)
    !
    REAL(r_std), INTENT(out) :: ztq(:,:)
    REAL(r_std), INTENT(out) :: zuv(:,:)
    !
    ! Local
    !
    INTEGER(i_std) :: iret, varid
    INTEGER(i_std) :: if, it
    INTEGER(i_std) :: tstart(3), tcount(3)
    INTEGER(i_std) :: imin(1), imax(1), firstif(1)
    INTEGER(i_std) :: nctstart, nctcount, inslabpos
    INTEGER(i_std) :: start_globtime, end_globtime
    INTEGER(i_std) :: timeid_tair, timeid_qair, timeid_precip, timeid_swdown
    INTEGER(i_std) :: timeid_lwdown, timeid_u, timeid_v, timeid_ps, timeid_tmp
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: time_tmp
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: rau
    CHARACTER(LEN=80) :: cellmethod
    !
    LOGICAL, SAVE :: first_call_readslab=.TRUE.
    !
    ALLOCATE(time_tmp(slab_size,nbtax))
    ALLOCATE(rau(nbpoint_loc,slab_size))
    !
    !
    ! Catch any stupid utilisation of this routine.
    !
    IF ( .NOT. is_root_prc) THEN
       CALL ipslerr (3,'forcing_readslab_root',"Cannot run this routine o other procs than root.",&
            &        "All the information on the forcing files is only lated on the root processor."," ")
    ENDIF
    !
    !Set some variables to zero
    !
    IF ( first_call_readslab ) THEN
       !
       preciptime(:) = 0
       !
       ! If the first file is only there to provide the last time step of the previous year, we
       ! do not need to read all. We will start 2 forcing time steps before the start of the first
       ! time interval requested.
       !
       imin=MINLOC(ABS(time(:,1)-time_int(1)))
       current_offset = MAX(imin(1)-2,1)
       !
       first_call_readslab = .FALSE.
       !
    ELSE       
       !
       ! Put back the cummulated time of rainfall into the global array
       !
       preciptime(position_slab(1):position_slab(2)) = preciptime(position_slab(1):position_slab(2)) + &
            &    prectime(1:slab_size)
       !
       ! Compute new offset
       !
       current_offset = position_slab(2)-2
       !
    ENDIF
    !
    ! Check that the slab size is not too large
    !
    IF ( (current_offset-1)+slab_size > nb_forcing_steps) THEN
       slab_size = nb_forcing_steps - (current_offset-1)
    ENDIF
    !
    ! 1.1 Check that the slab we have to read exists
    !
    IF ( slab_size > 0 ) THEN
       !
       start_globtime = current_offset
       end_globtime = (current_offset-1)+slab_size
       inslabpos=1
       WRITE(*,*) ">> Reading from global position ", start_globtime, "up to ", end_globtime
       !
       DO IF=MINVAL(time_sourcefile(start_globtime:end_globtime)),MAXVAL(time_sourcefile(start_globtime:end_globtime))
          !
          ! Get position of the part of the global time axis we need to read from this file
          !
          firstif = MINLOC(ABS(time_sourcefile-if))
          ! start = distance of start_globtime or nothing + 1 to follow netCDF convention.
          nctstart =  MAX((start_globtime-firstif(1)), 0)+1
          ! count = end index - start index + 1 
          nctcount = MIN((firstif(1)-1)+nbtime_perfile(if),end_globtime)-MAX(firstif(1),start_globtime)+1
          !
          !
          ! Read time over the indexes computed above in order to position the slab correctly
          ! 
          WRITE(*,*) ">> From file ", if," reading from position ", nctstart, "up to ", (nctstart-1)+nctcount
          !
          DO it =1,nbtax
             tstart(1) = nctstart
             tcount(1) = nctcount
             iret = NF90_GET_VAR(force_id(if), time_id(if,it), time_tmp(inslabpos:inslabpos+nctcount-1,it), tstart, tcount)
             IF (iret /= NF90_NOERR) THEN
                WRITE(*,*) TRIM(NF90_STRERROR(iret))
                WRITE(*,*) "Working on file ", IF," starting at ",tstart(1)," counting ",tcount(1)
                WRITE(*,*) "The data was to be written in to section ", inslabpos,":",inslabpos+nctcount-1," of time_tmp"
                CALL ipslerr (3,'forcing_readslab_root',"Could not read the time for the current interval."," "," ")
             ENDIF
             time_tmp(inslabpos:inslabpos+nctcount-1,it) = date0_file(if,it) + &
                  time_tmp(inslabpos:inslabpos+nctcount-1,it)*convtosec(if)/one_day
          ENDDO
          !
          ! 2.0 Find and read variables.
          !
          ! 2.1 Deal with air temperature and humidity as the fist and basic case
          !
          ! 
          !
          CALL forcing_varforslab(if, "Tair", nctstart, nctcount, inslabpos, tair, cellmethod)
          CALL forcing_attributetimeaxe(cellmethod, timeid_tair)
          !
          CALL forcing_varforslab(if, "Qair", nctstart, nctcount, inslabpos, qair, cellmethod)
          CALL forcing_attributetimeaxe(cellmethod, timeid_qair)
          !
          ! 2.2 Deal with rainfall and snowfall.
          !
          CALL forcing_varforslab(if, "Rainf", nctstart, nctcount, inslabpos, rainf, cellmethod)
          CALL forcing_attributetimeaxe(cellmethod, timeid_precip)
          !
          CALL forcing_varforslab(if, "Snowf", nctstart, nctcount, inslabpos, snowf, cellmethod)
          CALL forcing_attributetimeaxe(cellmethod, timeid_tmp)
          IF ( timeid_precip .NE. timeid_tmp) THEN
             CALL ipslerr(3, 'forcing_readslab_root','Rainf and Snwof have different time axes.', &
                  &         'Please check the forcing file to ensure both variable have the same cell_method.','')
          ENDIF
          !
          !
          ! 2.3 Deal with downward shorwave and longwave radiation
          !     The SW radiation can have 2 names SWdown_aerosol or SWdown. The first one is
          !     given priority
          !
          CALL forcing_varforslab(if, "SWdown", nctstart, nctcount, inslabpos, swdown, cellmethod)
          CALL forcing_attributetimeaxe(cellmethod, timeid_swdown)
          !
          CALL forcing_varforslab(if, "LWdown", nctstart, nctcount, inslabpos, lwdown, cellmethod)
          CALL forcing_attributetimeaxe(cellmethod, timeid_lwdown)
          !
          !
          ! 2.4 Deal with wind and pressure
          !
          CALL forcing_varforslab(if, "PSurf", nctstart, nctcount, inslabpos, ps, cellmethod)
          CALL forcing_attributetimeaxe(cellmethod, timeid_ps)
          !
          CALL forcing_varforslab(if, "Wind_E", nctstart, nctcount, inslabpos, u, cellmethod)
          CALL forcing_attributetimeaxe(cellmethod, timeid_u)
          !
          CALL forcing_varforslab(if, "Wind_N", nctstart, nctcount, inslabpos, v, cellmethod)
          CALL forcing_attributetimeaxe(cellmethod, timeid_v)
          !
          ! Verify on Tair that we have a credible field.
          !
          IF (   MINVAL(tair(:,inslabpos:inslabpos+nctcount-1)) < 100.0 .OR. &
               & MAXVAL(tair(:,inslabpos:inslabpos+nctcount-1)) > 500.0 ) THEN
             WRITE(*,*) "ERROR on range of Tair : ", MINVAL(tair(:,inslabpos:inslabpos+nctcount-1)), &
                  &                                  MAXVAL(tair(:,inslabpos:inslabpos+nctcount-1))
             CALL ipslerr(3, 'forcing_readslab_root','The air temperature is not in a credible range.', &
                  &         'Please verify your forcing file.','Are variables for all points to be simulated ?')
          ENDIF
          !
          ! Do the height of the variables T&Q and U&V
          !
          ! First the T&Q level
          !
          IF ( zheight ) THEN
             ztq(:,:) = zlev_fixed
          ELSE IF ( zsigma .OR. zhybrid ) THEN
             DO it=inslabpos,inslabpos+nctcount-1
                rau(:,it) = ps(:,it)/(cte_molr*tair(:,it))
                ztq(:,it) = (ps(:,it) - (zhybrid_a + zhybrid_b*ps(:,it)))/(rau(:,it) * cte_grav)
             ENDDO
          ELSE IF ( zlevels ) THEN
             CALL forcing_varforslab(IF, "Levels", nctstart, nctcount, inslabpos, ztq, cellmethod)
          ELSE
             CALL ipslerr(3, 'forcing_readslab_root','No case for the vertical levels was specified.', &
                  &         'We cannot determine the height for T and Q.','')
          ENDIF
          !
          ! Now the U&V level
          !
          IF ( zsamelev_uv ) THEN
             zuv(:,:) = ztq(:,:)
          ELSE
             IF ( zheight ) THEN
                zuv(:,:) = zlevuv_fixed
             ELSE IF ( zsigma .OR. zhybrid ) THEN
                DO it=inslabpos,inslabpos+nctcount-1
                   rau(:,it) = ps(:,it)/(cte_molr*tair(:,it))
                   zuv(:,it) = (ps(:,it) - (zhybriduv_a + zhybriduv_b*ps(:,it)))/(rau(:,it) * cte_grav)
                ENDDO
             ELSE IF ( zlevels ) THEN
                CALL forcing_varforslab(IF, "Levels_uv", nctstart, nctcount, inslabpos, zuv, cellmethod)
             ELSE
                CALL ipslerr(3, 'forcing_readslab_root','No case for the vertical levels was specified.', &
                     &         'We cannot determine the height for U and V.','stop readdim2')
             ENDIF
          ENDIF
          
          inslabpos = inslabpos+nctcount
          
       ENDDO
       !
       ! Use the read time of the slab to place it in the global variables. We consider 
       ! that we can do that on the first axis.
       !
       imin = MINLOC(ABS(time(:,1)-time_tmp(1,1)))
       position_slab(1) = imin(1)
       imax = MINLOC(ABS(time(:,1)-time_tmp(slab_size,1)))
       position_slab(2) = imax(1)
       !
       !
       IF ( position_slab(2)-position_slab(1) .GT. slab_size ) THEN
          WRITE(*,*) "Postition_slab =",  position_slab
          WRITE(*,*) "Interval read : ", position_slab(2)-position_slab(1)
          WRITE(*,*) "Time start and end : ", time(1,1), time(slab_size,1)
          DO it =1,nbtax
             WRITE(*,*) "Checking time_tmp on idex : ", it
             WRITE(*,*) "Time_tmp start and end : ",time_tmp(1,it), time_tmp(slab_size,it)
             imin = MINLOC(ABS(time(:,1)-time_tmp(1,it)))
             imax = MINLOC(ABS(time(:,1)-time_tmp(slab_size,it)))
             WRITE(*,*) "Interval read : ", imax(1)-imin(1)+1
          ENDDO
          WRITE(*,*) "current_offset, slab_size =", current_offset, slab_size
          CALL ipslerr (3,'forcing_readslab_root',"The time slab read does not fit the number of variables read.",&
               &        "Could there be an error in the time axis ?"," ")
       ENDIF
       !
       ! Transfer the global time axis into the time variables approriate for each variable. This way
       ! the time axis for each variable will be centered on the interval of validity. This is an essential assumption 
       ! the time interpolation to be done later.
       !
       WRITE(*,*) "We have found the following axes : ", time_axename(:)
       WRITE(*,*) "For Tair we are using time axis '",TRIM(time_axename(timeid_tair)),&
            &     "' with cell method ",TRIM(time_cellmethod(timeid_tair)),"."
       t_tair(1:slab_size) = time(position_slab(1):position_slab(2), timeid_tair)
       tbnd_tair(1:slab_size,:) = time_bounds(position_slab(1):position_slab(2),timeid_tair,:)
       !
       WRITE(*,*) "For Qair we are using time axis '",TRIM(time_axename(timeid_qair)),&
            &     "' with cell method ",TRIM(time_cellmethod(timeid_qair)),"."
       t_qair(1:slab_size) = time(position_slab(1):position_slab(2), timeid_qair)
       tbnd_qair(1:slab_size,:) = time_bounds(position_slab(1):position_slab(2),timeid_qair,:)
       !
       WRITE(*,*) "For Rainf and Snowf we are using time axis '",TRIM(time_axename(timeid_precip)),&
            &     "' with cell method ",TRIM(time_cellmethod(timeid_precip)),"."
       t_prec(1:slab_size) = time(position_slab(1):position_slab(2), timeid_precip)
       tbnd_prec(1:slab_size,:) = time_bounds(position_slab(1):position_slab(2),timeid_precip,:)
       prectime(1:slab_size) = preciptime(position_slab(1):position_slab(2))
       !
       WRITE(*,*) "For SWdown we are using time axis '",TRIM(time_axename(timeid_swdown)),&
            &     "' with cell method ",TRIM(time_cellmethod(timeid_swdown)),"."
       t_swdown(1:slab_size) = time(position_slab(1):position_slab(2), timeid_swdown)
       tbnd_swdown(1:slab_size,:) = time_bounds(position_slab(1):position_slab(2),timeid_swdown,:)
       !
       WRITE(*,*) "For LWdown we are using time axis '",TRIM(time_axename(timeid_lwdown)),&
            &     "' with cell method ",TRIM(time_cellmethod(timeid_lwdown)),"."
       t_lwdown(1:slab_size) = time(position_slab(1):position_slab(2), timeid_lwdown)
       tbnd_lwdown(1:slab_size,:) = time_bounds(position_slab(1):position_slab(2),timeid_lwdown,:)
       !
       WRITE(*,*) "For Wind_E we are using time axis '",TRIM(time_axename(timeid_u)),&
            &     "' with cell method ",TRIM(time_cellmethod(timeid_u)),"."
       t_u(1:slab_size) = time(position_slab(1):position_slab(2), timeid_u)
       tbnd_u(1:slab_size,:) = time_bounds(position_slab(1):position_slab(2),timeid_u,:)
       !
       WRITE(*,*) "For Wind_N we are using time axis '",TRIM(time_axename(timeid_v)),&
            &     "' with cell method ",TRIM(time_cellmethod(timeid_v)),"."
       t_v(1:slab_size) = time(position_slab(1):position_slab(2), timeid_v)
       tbnd_v(1:slab_size,:) = time_bounds(position_slab(1):position_slab(2),timeid_v,:)
       !
       WRITE(*,*) "For PSurf we are using time axis '",TRIM(time_axename(timeid_ps)),&
            &     "' with cell method ",TRIM(time_cellmethod(timeid_ps)),"."
       t_ps(1:slab_size) = time(position_slab(1):position_slab(2), timeid_ps)
       tbnd_ps(1:slab_size,:) = time_bounds(position_slab(1):position_slab(2),timeid_ps,:)
       !
    ELSE
       CALL ipslerr (2,'forcing_readslab_root',"We have reached the end of the slabs we can read.",&
            &          "The calling program needs to manage this situation"," ")
    ENDIF
    !
    ! Have we read to the end of the files ?
    ! 
    IF ( current_offset+slab_size >= nb_forcing_steps ) THEN
       end_of_file = .TRUE.
    ELSE
       end_of_file = .FALSE.
    ENDIF
    !
    IF ( ALLOCATED(rau) ) DEALLOCATE(rau)
    IF ( ALLOCATED(time_tmp) ) DEALLOCATE(time_tmp)
    !
  END SUBROUTINE forcing_readslab_root
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_reindex3d
!!
!>\BRIEF      
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================
!
!=============================================================================================
!
  SUBROUTINE forcing_reindex3d(nbi, nbj, tin, slab_in, nbout, tout, slab_out, tstart, reindex)
    !
    ! ARGUMENTS
    !
    INTEGER(i_std), INTENT(in) :: nbi, nbj, tin, nbout, tout
    REAL(r_std), INTENT(in)    :: slab_in(nbi,nbj,tin)
    REAL(r_std), INTENT(out)   :: slab_out(nbout,tout)
    INTEGER(i_std), INTENT(in) :: tstart
    INTEGER(i_std), INTENT(in) :: reindex(nbout,2)
    !
    ! LOCAL
    !
    INTEGER(i_std) :: is, in
    !
    DO is=1,tin
       DO in=1,nbout
          slab_out(in,tstart+(is-1)) = slab_in(reindex(in,1),reindex(in,2),is)
       ENDDO
    ENDDO
    !
  END SUBROUTINE forcing_reindex3d
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_reindex2d
!!
!>\BRIEF      
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================
!
!=============================================================================================
!
  SUBROUTINE forcing_reindex2d(nbi, nbj, slab_in, nbout, slab_out, reindex)
    !
    ! ARGUMENTS
    !
    INTEGER(i_std), INTENT(in) :: nbi, nbj, nbout
    REAL(r_std), INTENT(in)    :: slab_in(nbi,nbj)
    REAL(r_std), INTENT(out)   :: slab_out(nbout)
    INTEGER(i_std), INTENT(in) :: reindex(nbout,2)
    !
    ! LOCAL
    !
    INTEGER(i_std) :: in
    !
    DO in=1,nbout
       slab_out(in) = slab_in(reindex(in,1),reindex(in,2))
    ENDDO
    !
  END SUBROUTINE forcing_reindex2d
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_reindex2dt
!!
!>\BRIEF      
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================
!
!=============================================================================================
!
  SUBROUTINE forcing_reindex2dt(nbin, tin, slab_in, nbout, tout, slab_out, tstart, reindex)
    !
    ! ARGUMENTS
    !
    INTEGER(i_std), INTENT(in) :: nbin, tin, nbout, tout
    REAL(r_std), INTENT(in)    :: slab_in(nbin,tin)
    REAL(r_std), INTENT(out)   :: slab_out(nbout,tout)
    INTEGER(i_std), INTENT(in) :: tstart
    INTEGER(i_std), INTENT(in) :: reindex(nbout)
    !
    ! LOCAL
    !
    INTEGER(i_std) :: is, in
    !
    DO is=1,tin
       DO in=1,nbout
          slab_out(in,tstart+(is-1)) = slab_in(reindex(in),is)
       ENDDO
    ENDDO
    !
  END SUBROUTINE forcing_reindex2dt
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_reindex1d
!!
!>\BRIEF      
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================
!
!=============================================================================================
!
  SUBROUTINE forcing_reindex1d(nbin, slab_in, nbout, slab_out, reindex)
    !
    ! ARGUMENTS
    !
    INTEGER(i_std), INTENT(in) :: nbin, nbout
    REAL(r_std), INTENT(in)    :: slab_in(nbin)
    REAL(r_std), INTENT(out)   :: slab_out(nbout)
    INTEGER(i_std), INTENT(in) :: reindex(nbout)
    !
    ! LOCAL
    !
    INTEGER(i_std) :: is
    !
    DO is=1,nbout
       slab_out(is) = slab_in(reindex(is))
    ENDDO
    !
  END SUBROUTINE forcing_reindex1d
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_reindex2to1
!!
!>\BRIEF      
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================
!
!=============================================================================================
!
  SUBROUTINE forcing_reindex2to1(nbi, nbj, slab_in, nbout, slab_out, reindex)
    !
    ! ARGUMENTS
    !
    INTEGER(i_std), INTENT(in) :: nbi, nbj, nbout
    REAL(r_std), INTENT(in)    :: slab_in(nbi,nbj)
    REAL(r_std), INTENT(out)   :: slab_out(nbout)
    INTEGER(i_std), INTENT(in) :: reindex(nbout)
    !
    ! LOCAL
    !
    INTEGER(i_std) :: i, j, is
    !
    DO is=1,nbout
       j = INT((reindex(is)-1)/nbi)+1
       i = (reindex(is)-(j-1)*nbi)
       slab_out(is) = slab_in(i,j)
    ENDDO
    !
  END SUBROUTINE forcing_reindex2to1
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_reindex1to2
!!
!>\BRIEF      
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================
!
!=============================================================================================
!
  SUBROUTINE forcing_reindex1to2(nbin, slab_in, nbi, nbj, slab_out, reindex)
    !
    ! ARGUMENTS
    !
    INTEGER(i_std), INTENT(in) :: nbin, nbi, nbj
    REAL(r_std), INTENT(in)    :: slab_in(nbin)
    REAL(r_std), INTENT(out)   :: slab_out(nbi, nbj)
    INTEGER(i_std), INTENT(in) :: reindex(nbin)
    !
    ! LOCAL
    !
    INTEGER(i_std) :: i, j, is
    !
    DO is=1,nbin
       j = INT((reindex(is)-1)/nbi)+1
       i = (reindex(is)-(j-1)*nbi)
       slab_out(i,j) = slab_in(is)
    ENDDO
    !
  END SUBROUTINE forcing_reindex1to2
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_close
!!
!>\BRIEF  Close all forcing files
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================
!
!=============================================================================================
!
  SUBROUTINE forcing_close()

    INTEGER(i_std) :: ierr, if

    DO if=1,nb_forcefile
       ierr = NF90_CLOSE(force_id(if))
    ENDDO

  END SUBROUTINE forcing_close
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_printdate
!!
!>\BRIEF    Print the date in a human readable format.  
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================
!
!=============================================================================================
!
  SUBROUTINE forcing_printdate(julian_day, message, wunit)
    !
    REAL(r_std), INTENT(in) :: julian_day
    CHARACTER(len=*), INTENT(in) :: message
    INTEGER, INTENT(in), OPTIONAL :: wunit
    !
    !
    !
    INTEGER(i_std) :: year, month, day, hours, minutes, seci
    REAL(r_std) :: sec
    !
    CALL ju2ymds (julian_day, year, month, day, sec)
    hours = INT(sec/3600)
    sec = sec - 3600 * hours
    minutes = INT(sec / 60)
    sec = sec - 60 * minutes
    seci = INT(sec)
    !
    IF (PRESENT(wunit)) THEN
       WRITE(wunit,'(I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2," > ", A60)') &
            &            year, month, day, hours, minutes, seci, message
    ELSE
       WRITE(*,'(I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2," > ", A60)') &
            &            year, month, day, hours, minutes, seci, message
    ENDIF
    !
  END SUBROUTINE forcing_printdate
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_printpoint_forgrid
!!
!>\BRIEF     Together with the date print some sample values. Useful for checking the forcing. 
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================
!
!=============================================================================================
!
  SUBROUTINE forcing_printpoint_forgrid(julian_day, lon_pt, lat_pt, var, message)
    !
    REAL(r_std), INTENT(in) :: julian_day
    REAL(r_std), INTENT(in) :: lon_pt, lat_pt
    REAL(r_std), INTENT(in) :: var(:)
    CHARACTER(len=*), INTENT(in) :: message
    !
    !
    !
    INTEGER(i_std) :: year, month, day, hours, minutes, seci
    REAL(r_std) :: sec
    INTEGER(i_std) :: lon_ind, lat_ind, ind
    INTEGER(i_std), DIMENSION(1) :: i,j,k
    !
    ! Check if there is anything to be done
    !
    IF ( MAX(lon_pt, lat_pt) > 360.0 ) THEN
       RETURN
    ENDIF
    !
    ! Convert time first
    !
    CALL ju2ymds (julian_day, year, month, day, sec)
    hours = INT(sec/3600)
    sec = sec - 3600 * hours
    minutes = INT(sec / 60)
    sec = sec - 60 * minutes
    seci = INT(sec)
    !
    ! Get the local to be analysed
    !
    i = MINLOC(ABS(lon_loc(:,1)-lon_pt))
    j = MINLOC(ABS(lat_loc(1,:)-lat_pt))
    ind = (j(1)-1)*iim_loc+i(1)
    k = MINLOC(ABS(lindex_loc(:)-ind))
    !
    WRITE(*,"(I2.2,':',I2.2,':',I2.2,' Loc : ', F5.1,',', F5.1,'(i=',I6,') Value = ',F12.4,A40)") &
         & hours, minutes, seci, lon_loc(i(1),1), lat_loc(1,j(1)), k(1), var(k(1)), message
    !
  END SUBROUTINE forcing_printpoint_forgrid
  !!  =============================================================================================================================
!! SUBROUTINE: forcing_printpoint_forgrid2d
!!
!>\BRIEF     Together with the date print some sample values. Useful for checking the forcing. 
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE forcing_printpoint_forgrid2d(julian_day, lon_pt, lat_pt, var, message)
    !
    REAL(r_std), INTENT(in) :: julian_day
    REAL(r_std), INTENT(in) :: lon_pt, lat_pt
    REAL(r_std), INTENT(in) :: var(:,:)
    CHARACTER(len=*), INTENT(in) :: message
    !
    !
    !
    INTEGER(i_std) :: year, month, day, hours, minutes, seci
    REAL(r_std) :: sec
    INTEGER(i_std) :: lon_ind, lat_ind
    INTEGER(i_std), DIMENSION(1) :: i,j
    !
    ! Check if there is anything to be done
    !
    IF ( MAX(lon_pt, lat_pt) > 360.0 ) THEN
       RETURN
    ENDIF
    !
    ! Convert time first
    !
    CALL ju2ymds (julian_day, year, month, day, sec)
    hours = INT(sec/3600)
    sec = sec - 3600 * hours
    minutes = INT(sec / 60)
    sec = sec - 60 * minutes
    seci = INT(sec)
    !
    ! Get the local to be analysed
    !
    i = MINLOC(ABS(lon_loc(:,1)-lon_pt))
    j = MINLOC(ABS(lat_loc(1,:)-lat_pt))
    !
    WRITE(*,"(I2.2,':',I2.2,':',I2.2,' Loc : ', F5.1,',', F5.1,'(i=',I6,') Value = ',F12.4,A40)") &
         & hours, minutes, seci, lon_loc(i(1),1), lat_loc(1,j(1)), i(1), j(1), var(i(1),j(1)), message
    
  END SUBROUTINE forcing_printpoint_forgrid2d

!!  =============================================================================================================================
!! SUBROUTINE: forcing_printpoint_gen
!!
!>\BRIEF       Together with the date print some sample values. Useful for checking the forcing. 
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================
!
!=============================================================================================
!
  SUBROUTINE forcing_printpoint_gen(julian_day, lon_pt, lat_pt, nbind, lalo_in, var, message, ktest)
    !
    REAL(r_std), INTENT(in) :: julian_day
    REAL(r_std), INTENT(in) :: lon_pt, lat_pt
    INTEGER(i_std), INTENT(in) :: nbind
    REAL(r_std), INTENT(in) :: lalo_in(:,:)
    REAL(r_std), INTENT(in) :: var(:)
    CHARACTER(len=*), INTENT(in) :: message
    INTEGER(i_std), OPTIONAL, INTENT(out) :: ktest
    !
    !
    !
    INTEGER(i_std) :: year, month, day, hours, minutes, seci
    REAL(r_std) :: sec, mindist
    REAL(r_std), ALLOCATABLE, DIMENSION(:) :: dist, refdist
    INTEGER(i_std) :: lon_ind, lat_ind, ind
    INTEGER(i_std) :: i, imin(1)
    REAL(r_std), PARAMETER :: mincos  = 0.0001
    REAL(r_std), PARAMETER :: pi = 3.141592653589793238
    REAL(r_std), PARAMETER :: R_Earth = 6378000.
    !
    ! Check if there is anything to be done
    !
    IF ( MAX(lon_pt, lat_pt) > 360.0 ) THEN
       IF ( PRESENT(ktest) ) ktest = -1
       RETURN
    ENDIF
    !
    ! Allocate memory
    !
    ALLOCATE(dist(nbind))
    ALLOCATE(refdist(nbind))
    !
    ! Convert time first
    !
    CALL ju2ymds (julian_day, year, month, day, sec)
    hours = INT(sec/3600)
    sec = sec - 3600 * hours
    minutes = INT(sec / 60)
    sec = sec - 60 * minutes
    seci = INT(sec)
    !
    ! Get the location to be analysed
    !
    DO i=1,nbind
       dist(i) = acos( sin(lat_pt*pi/180)*sin(lalo_in(i,1)*pi/180) + &
            &    cos(lat_pt*pi/180)*cos(lalo_in(i,1)*pi/180)*&
            &    cos((lalo_in(i,2)-lon_pt)*pi/180) ) * R_Earth
    ENDDO
    !
    ! Look for the next grid point closest to the one with the smalest distance.
    !
    imin = MINLOC(dist)
    DO i=1,nbind
       refdist(i) = acos( sin(lalo_in(imin(1),1)*pi/180)*sin(lalo_in(i,1)*pi/180) + &
            &       cos(lalo_in(imin(1),1)*pi/180)*cos(lalo_in(i,1)*pi/180) * &
            &       cos((lalo_in(i,2)-lalo_in(imin(1),2))*pi/180) ) * R_Earth
    ENDDO
    refdist(imin(1)) =  MAXVAL(refdist)
    mindist = MINVAL(refdist)
    !
    ! Are we closer than the closest points ?
    !
    IF ( PRESENT(ktest) ) ktest = -1
    IF ( dist(imin(1)) <= mindist ) THEN
       !
       WRITE(*,"(I2.2,':',I2.2,':',I2.2,' Loc : ', F6.1,',', F6.1,'(i=',I6,') Value = ',F12.4,A38)") &
            & hours, minutes, seci, lalo_in(imin(1),2), lalo_in(imin(1),1), imin(1), var(imin(1)), message
       !
       IF ( PRESENT(ktest) ) ktest = imin(1)
    ENDIF
    !
  END SUBROUTINE forcing_printpoint_gen
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_getglogrid
!!
!>\BRIEF       Routine to read the full grid of the forcing file.
!!
!! DESCRIPTION: The data is stored in the saved variables of the module and serve all other routines. 	  
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE forcing_getglogrid (nbfiles, filename, iim_tmp, jjm_tmp, nbpoint_tmp, closefile, landonly_arg)
    !
    ! This routine reads the global grid information from the forcing file and generates all the 
    ! information needed for this grid.
    ! 
    ! ARGUMENTS
    !
    INTEGER(i_std), INTENT(in)    :: nbfiles
    CHARACTER(LEN=*), INTENT(in)  :: filename(:)
    INTEGER(i_std), INTENT(out)   :: iim_tmp, jjm_tmp, nbpoint_tmp
    LOGICAL, INTENT(in)           :: closefile
    LOGICAL, OPTIONAL, INTENT(in) :: landonly_arg
    !
    ! LOCAL
    !
    INTEGER(i_std) :: iret, iv, if, lll
    CHARACTER(LEN=20) :: dimname, varname
    CHARACTER(LEN=60) :: lon_units, lat_units, units
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:) :: dimids, londim_ids, latdim_ids
    INTEGER(i_std) :: lon_id, lat_id, land_id, lon_nbdims, lat_nbdims, land_nbdims
    INTEGER(i_std) :: lonvar_id, latvar_id, landvar_id
    LOGICAL :: dump_mask
    INTEGER(i_std) :: ik, i, j, iff, ndimsvar
    ! Read a test variabe
    CHARACTER(len=8) :: testvarname='tair'
    INTEGER(i_std)   :: testvar_id, contfrac_id
    REAL(r_std) :: testvar_missing, contfrac_missing
    REAL(r_std), ALLOCATABLE, DIMENSION(:)   :: testvar
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: testvar2d, contfrac2d
    !
    ! 0.0 Check variables are allocated
    !
    IF ( .NOT. ALLOCATED(force_id)) ALLOCATE(force_id(nbfiles))
    IF ( .NOT. ALLOCATED(id_unlim)) ALLOCATE(id_unlim(nbfiles))
    IF ( .NOT. ALLOCATED(nb_atts)) ALLOCATE(nb_atts(nbfiles))
    IF ( .NOT. ALLOCATED(ndims)) ALLOCATE(ndims(nbfiles))
    IF ( .NOT. ALLOCATED(nvars)) ALLOCATE( nvars(nbfiles))
    !
    ! 0.1 Are we one the root proc ?
    !
    IF ( .NOT. is_root_prc ) THEN
       CALL ipslerr (3,'forcing_getglogrid'," This routine can only be called on the root processor.", " ", " ")
    ENDIF
    !
    ! The default behaviour is to provide only land points to the calling program.
    ! But for forcing ocean model there is also the option to pass on land and ocean values.
    ! When the grid is initialized landonly_tmp=.FALSE. has to be set to obtian this behaviour.
    !
    IF ( PRESENT(landonly_arg) ) THEN
       landonly=landonly_arg
    ELSE
       landonly=.TRUE.
    ENDIF
    !
    ! 1.0 Open the netCDF file and get basic dimensions
    !
    DO iff=1,nbfiles
       !
       iret = NF90_OPEN(filename(iff), NF90_NOWRITE, force_id(iff))
       IF (iret /= NF90_NOERR) THEN
          CALL ipslerr (3,'forcing_getglogrid',"Error opening the forcing file :", filename(iff), " ")
       ENDIF
       !
       iret = NF90_INQUIRE (force_id(iff), nDimensions=ndims(iff), nVariables=nvars(iff), &
            nAttributes=nb_atts(iff), unlimitedDimId=id_unlim(iff))
       !
       !
       ! 2.0 Read the dimensions found in the forcing file. Only deal with the spatial dimension as
       !     time is an unlimited dimension.
       !
       DO iv=1,ndims(iff)
          !
          iret = NF90_INQUIRE_DIMENSION (force_id(iff), iv, name=dimname, len=lll)
          IF (iret /= NF90_NOERR) THEN
             CALL ipslerr (3,'forcing_getglogrid',"Could not get size of dimensions in file : ", filename(iff), " ")
          ENDIF
          !
          SELECT CASE(lowercase(dimname))
             !
          CASE("west_east")
             CALL forcing_checkdim(iff, filename, iim_glo, lon_id, lll, iv)
          CASE("south_north")
             CALL forcing_checkdim(iff, filename, jjm_glo, lat_id, lll, iv)
          CASE("longitude")
             CALL forcing_checkdim(iff, filename, iim_glo, lon_id, lll, iv)
          CASE("latitude")
             CALL forcing_checkdim(iff, filename, jjm_glo, lat_id, lll, iv)
          CASE("lon")
             CALL forcing_checkdim(iff, filename, iim_glo, lon_id, lll, iv)
          CASE("lat")
             CALL forcing_checkdim(iff, filename, jjm_glo, lat_id, lll, iv)
          CASE("nav_lon")
             CALL forcing_checkdim(iff, filename, iim_glo, lon_id, lll, iv)
          CASE("nav_lat")
             CALL forcing_checkdim(iff, filename, jjm_glo, lat_id, lll, iv)
          CASE("x")
             CALL forcing_checkdim(iff, filename, iim_glo, lon_id, lll, iv)
          CASE("y")
             CALL forcing_checkdim(iff, filename, jjm_glo, lat_id, lll, iv)
          CASE("land")
             CALL forcing_checkdim(iff, filename, nbland_glo, land_id, lll, iv)
          END SELECT
          !
       ENDDO
    ENDDO
    !
    ! 3.0 Read the spatial coordinate variables found in the first file.
    !
    ALLOCATE(dimids(NF90_MAX_VAR_DIMS), londim_ids(NF90_MAX_VAR_DIMS), latdim_ids(NF90_MAX_VAR_DIMS))
    lonvar_id = -1
    latvar_id = -1
    landvar_id = -1
    testvar_id = -1
    contfrac_id = -1
    ! Count the number of time axis we have
    nbtax = 0
    !
    DO iv=1,nvars(1)
       !
       iret = NF90_INQUIRE_VARIABLE(force_id(1), iv, name=varname, ndims=ndimsvar, dimids=dimids)
       iret = NF90_GET_ATT(force_id(1), iv, 'units', units)
       !
       ! Check that we have the longitude
       !
       IF ( INDEX(lowercase(varname), 'lon') > 0 .AND. INDEX(lowercase(units), 'east') > 0 ) THEN
          lonvar_id = iv
          lon_units=units
          lon_nbdims = ndimsvar
          londim_ids = dimids
       ENDIF
       !
       ! Check that we have the latitude
       !
       IF ( INDEX(lowercase(varname), 'lat') > 0 .AND. INDEX(lowercase(units), 'north') > 0) THEN
          latvar_id = iv
          lat_units=units
          lat_nbdims = ndimsvar
          latdim_ids = dimids
       ENDIF
       !
       ! Check that we have the land re-indexing table
       !
       IF ( INDEX(lowercase(varname), 'land') > 0 ) THEN
          landvar_id = iv
          land_nbdims = ndimsvar
          latdim_ids = dimids
       ENDIF
       !
       ! Check if we have the contfrac variable
       !
       IF ( INDEX(lowercase(varname), 'contfrac') > 0 ) THEN
          contfrac_id = iv
          iret = NF90_GET_ATT(force_id(1), iv, "missing_value", contfrac_missing)
          IF (iret /= NF90_NOERR) THEN
             ! No missing_value found, try to read _FillValue instead
             iret = NF90_GET_ATT(force_id(1), iv, "_FillValue", contfrac_missing)
             IF (iret /= NF90_NOERR) THEN
                WRITE(*,*) TRIM(nf90_strerror(iret))
                WRITE(*,*) " >> No _FillValue or missing_value found for contfrac"
                contfrac_missing=0.0
             END IF
          ENDIF
       ENDIF
       !
       ! Find the test variable
       !
       IF ( INDEX(lowercase(varname), TRIM(testvarname)) > 0 ) THEN
          testvar_id = iv
          iret = NF90_GET_ATT(force_id(1), iv, "missing_value", testvar_missing)
          IF (iret /= NF90_NOERR) THEN
             ! No missing_value found, try to read _FillValue instead
             iret = NF90_GET_ATT(force_id(1), iv, "_FillValue", testvar_missing)
             IF (iret /= NF90_NOERR) THEN
                WRITE(*,*) TRIM(nf90_strerror(iret))
                WRITE(*,*) " >> No _FillValue or missing_value found for variable=",varname
                testvar_missing=-1
             END IF
          ENDIF
       ENDIF
       !
       ! If we come across time get the calendar and archive it.
       !
       IF ( INDEX(lowercase(units),'seconds since') > 0 .OR. &
          &  INDEX(lowercase(units),'minutes since') > 0 .OR. &
          &  INDEX(lowercase(units),'hours since') > 0) THEN 
          !
          ! Get calendar used for the time axis
          !
          iret = NF90_GET_ATT(force_id(1), iv, "calendar", calendar)
          IF (iret == NF90_NOERR) THEN
             WRITE(*,*) ">> Setting the calendar to ",calendar
          ELSE 
             WRITE(*,*) ">> Keeping proleptic Gregorian calendar" 
             calendar="proleptic_gregorian"
          ENDIF
          !
          nbtax = nbtax + 1
          !
       ENDIF
    ENDDO
    !
    ! 4.0 Verification that we have found both coordinate variables and the land point index
    !
    IF ( ( lonvar_id < 0 ) .AND. ( INDEX(lowercase(lon_units), 'east') <= 0 ) ) THEN
       CALL ipslerr (3,'forcing_getglogrid',"Have not found a valid longitude. Units = ", lon_units, " ")
    ENDIF
    IF ( ( latvar_id < 0 ) .AND. ( INDEX(lowercase(lat_units), 'north') <= 0 ) ) THEN
       CALL ipslerr (3,'forcing_getglogrid',"Have not found a valid latitude. Units = : ", lat_units, " ")
    ENDIF
    IF ( landvar_id < 0 ) THEN
       CALL ipslerr (1,'forcing_getglogrid',"No reindexing table was found. ", &
            &           "This forcing file is not compressed by gathering.", " ")
    ENDIF
    !
    ! 5.0 Allocate the space for the global variables and read them.
    !
    IF ( .NOT. ALLOCATED(lon_glo)) ALLOCATE(lon_glo(iim_glo, jjm_glo))
    IF ( .NOT. ALLOCATED(lat_glo)) ALLOCATE(lat_glo(iim_glo, jjm_glo))
    !
    IF ( lon_nbdims == 2 .AND. lat_nbdims == 2 ) THEN
       iret = NF90_GET_VAR(force_id(1), lonvar_id, lon_glo)
       iret = NF90_GET_VAR(force_id(1), latvar_id, lat_glo)
    ELSE IF ( lon_nbdims == 1 .AND. lat_nbdims == 1 ) THEN
       DO iv=1,jjm_glo
          iret = NF90_GET_VAR(force_id(1), lonvar_id, lon_glo(:,iv))
       ENDDO
       DO iv=1,iim_glo
          iret = NF90_GET_VAR(force_id(1), latvar_id, lat_glo(iv,:))
       ENDDO
    ELSE
       WRITE(*,*) "Found dimensions for lon and lat :", lon_nbdims, lat_nbdims
       CALL ipslerr (3,'forcing_getglogrid',"Longitude and Latitude variables do not have the right dimensions.", " ", " ")
    ENDIF
    !
    ! If we have a indexing variable then the data is compressed by gathering, else we have to construct it.
    !
    compressed = .FALSE.
    IF ( landvar_id > 0 ) THEN
       IF ( .NOT. ALLOCATED(lindex_glo)) ALLOCATE(lindex_glo(nbland_glo))
       iret = NF90_GET_VAR(force_id(1), landvar_id, lindex_glo)
       compressed = .TRUE.
    ENDIF
    !
    IF ( .NOT. ALLOCATED(mask_glo)) ALLOCATE(mask_glo(iim_glo, jjm_glo)) 
    !
    ! Get the land/sea mask and contfrac based on a test variable, if contfract is not available. Else
    ! we use the contfrac variable.
    !
    IF ( compressed ) THEN
       IF ( landonly ) THEN
          IF ( .NOT. ALLOCATED(contfrac_glo)) ALLOCATE(contfrac_glo(nbland_glo))
          IF ( .NOT. ALLOCATED(testvar)) ALLOCATE(testvar(nbland_glo))
          CALL forcing_contfrac1d(force_id(1), testvar_id, contfrac_id, testvar)
          nbpoint_glo = nbland_glo
       ELSE
          WRITE(*,*) "forcing_tools cannot provide data over ocean points as the"
          WRITE(*,*) "data in the file is compressed by gathering land points."
          WRITE(*,*) "Fatal error"
          CALL ipslerr (3,'forcing_getglogrid',"forcing_tools cannot provide data over ocean points as the", &
               &                               "data in the file is compressed by gathering land points.", " ")
       ENDIF
    ELSE
       IF ( .NOT. ALLOCATED(testvar2d)) ALLOCATE(testvar2d(iim_glo, jjm_glo))
       IF ( .NOT. ALLOCATED(contfrac2d)) ALLOCATE(contfrac2d(iim_glo, jjm_glo))
       CALL forcing_contfrac2d(force_id(1), testvar_id, contfrac_id, testvar2d, contfrac2d, &
            & testvar_missing, contfrac_missing, nbland_glo)
       !
       ! We have found a variable on which we can count the number of land points. We can build
       ! the indexing table and gather the information needed.
       !
       IF ( landonly ) THEN
          nbpoint_glo = nbland_glo
          IF ( .NOT. ALLOCATED(lindex_glo)) ALLOCATE(lindex_glo(nbpoint_glo))
          IF ( .NOT. ALLOCATED(contfrac_glo)) ALLOCATE(contfrac_glo(nbpoint_glo))
          IF ( .NOT. ALLOCATED(testvar)) ALLOCATE(testvar(nbpoint_glo))
          IF ( contfrac_id > 0 ) THEN
             CALL forcing_buildindex(contfrac2d, contfrac_missing, landonly, lindex_glo, contfrac_glo)
             CALL forcing_reindex(iim_glo, jjm_glo, testvar2d, nbland_glo, testvar, lindex_glo)
          ELSE
             CALL forcing_buildindex(testvar2d, testvar_missing, landonly, lindex_glo, testvar)
             contfrac_glo(:) = 1.0
          ENDIF
       ELSE
          nbpoint_glo = iim_glo*jjm_glo
          IF ( .NOT. ALLOCATED(lindex_glo)) ALLOCATE(lindex_glo(nbpoint_glo))
          IF ( .NOT. ALLOCATED(contfrac_glo)) ALLOCATE(contfrac_glo(nbpoint_glo))
          IF ( .NOT. ALLOCATED(testvar)) ALLOCATE(testvar(nbpoint_glo))
          IF ( contfrac_id > 0 ) THEN
             CALL forcing_buildindex(contfrac2d, contfrac_missing, landonly, lindex_glo, contfrac_glo)
             CALL forcing_reindex(iim_glo, jjm_glo, testvar2d, nbland_glo, testvar, lindex_glo)
          ELSE
             CALL forcing_buildindex(testvar2d, testvar_missing, landonly, lindex_glo, testvar)
             contfrac_glo(:) = 1.0
          ENDIF
       ENDIF
       !
    ENDIF
    !
    !
    ! We assume that if the forcing file is closed at the end of this subroutine this is a test
    ! or initialisation of the grids. So we will dump the mask in a netCDF file for the user to 
    ! check.
    !
    dump_mask = closefile 
    CALL forcing_checkindex(dump_mask, testvarname, testvar)
    !
    !
    ! 8.0 Close file if needed
    !
    IF ( closefile ) THEN
       CALL forcing_close()
    ENDIF
    !
    ! Prepare variables to be returnned to calling subroutine.
    !
    iim_tmp = iim_glo
    jjm_tmp = jjm_glo
    nbpoint_tmp = nbpoint_glo
    !
    ! Clean up !
    !
    IF ( ALLOCATED(dimids) ) DEALLOCATE(dimids)
    IF ( ALLOCATED(londim_ids) ) DEALLOCATE(londim_ids)
    IF ( ALLOCATED(latdim_ids) ) DEALLOCATE(latdim_ids)
    IF ( ALLOCATED(testvar) ) DEALLOCATE(testvar)
    IF ( ALLOCATED(testvar2d) ) DEALLOCATE(testvar2d)
    IF ( ALLOCATED(contfrac2d) ) DEALLOCATE(contfrac2d)
    !
  END SUBROUTINE forcing_getglogrid
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_buildindex
!!
!>\BRIEF      
!!
!! DESCRIPTION:	When the forcing file does not contain compressed variables we need
!!              to build the land index variable from the mask defined by missing variables in 
!!              a test variable.  
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE forcing_buildindex(var2d, var_missing, landonly, lindex, var)
    !
    ! When the forcing file does not contain compressed variables we need
    ! to build the land index variable from the mask defined by missing variables in 
    ! a test variable.
    !
    ! Arguments
    !
    REAL(r_std), INTENT(in) :: var2d(:,:)
    REAL(r_std), INTENT(in) :: var_missing
    LOGICAL, INTENT(in)     :: landonly
    INTEGER(i_std), INTENT(out) :: lindex(:)
    REAL(r_std), INTENT(out) :: var(:)
    !
    ! Local
    !
    INTEGER(i_std) :: i,j,k
    !
    k=0
    DO i=1,iim_glo
       DO j=1,jjm_glo
          IF ( landonly ) THEN
             IF ( var2d(i,j) /= var_missing ) THEN
                k = k + 1
                lindex(k) = (j-1)*iim_glo+i 
                var(k) = var2d(i,j)
             ENDIF
          ELSE
             ! When we take all point, no test is performed.
             k = k + 1
             lindex(k) = (j-1)*iim_glo+i 
             var(k) = var2d(i,j)
          ENDIF
       ENDDO
    ENDDO
    !
  END SUBROUTINE forcing_buildindex
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_contfrac1d
!!
!>\BRIEF      
!!
!! DESCRIPTION:	 This routine build the land/mask if needed and gets the contfrac variable from forcing file.
!!               Here we treat the case where the variables are compressed by gathering. Thus only
!!               land points are available in the file. 
!!
!! \n
!_ ==============================================================================================================================
!
!=============================================================================================
!
  SUBROUTINE forcing_contfrac1d(ifile, testvar_id, contfrac_id, testvar)
    !
    ! This routine build the land/mask if needed and gets the contfrac variable from forcing file.
    ! Here we treat the case where the variables are compressed by gathering. Thus only
    ! land points are available in the file.
    ! 
    ! ARGUMENTS
    !
    INTEGER(i_std), INTENT(in)               :: ifile
    INTEGER(i_std), INTENT(in)               :: testvar_id, contfrac_id
    REAL(r_std), DIMENSION(:), INTENT(inout) :: testvar 
    !
    ! LOCAL
    !
    INTEGER(i_std)                           :: it, iret
    INTEGER(i_std), DIMENSION(3)             :: start, count
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: contfrac2d
    !
    ! First determine the contfrac variable
    !
    IF ( contfrac_id > 0 ) THEN
       iret = NF90_INQUIRE_VARIABLE(ifile, contfrac_id, ndims=it)
       IF ( it == 1 ) THEN
          start = (/1,1,0/)
          count = (/nbpoint_glo,1,0/)
          iret = NF90_GET_VAR(ifile, contfrac_id, contfrac_glo, start, count)
          IF (iret /= NF90_NOERR) THEN
             WRITE(*,*) TRIM(nf90_strerror(iret))
             CALL ipslerr (3,'forcing_contfrac1d',"Error reading contfrac ", " ", " ")
          ENDIF
       ELSE IF ( it == 2 ) THEN
          ALLOCATE(contfrac2d(iim_glo,jjm_glo))
          start = (/1,1,0/)
          count = (/iim_glo,jjm_glo,0/)
          iret = NF90_GET_VAR(ifile, contfrac_id, contfrac2d)
          IF (iret /= NF90_NOERR) THEN
             WRITE(*,*) TRIM(nf90_strerror(iret))
             CALL ipslerr (3,'forcing_contfrac1d',"Error reading contfrac ", " ", " ")
          ENDIF
          CALL forcing_reindex(iim_glo, jjm_glo, contfrac2d, nbpoint_glo, contfrac_glo, lindex_glo)
          DEALLOCATE(contfrac2d)
       ELSE
          CALL ipslerr (3,'forcing_contfrac1d',"Contfrac has a dimension larger than 2. ", &
               "We do not know how to handle this.", " ")
       ENDIF
    ELSE
       contfrac_glo(:) = 1.0
    ENDIF
    !
    ! Read our test variable 
    !
    iret = NF90_INQUIRE_VARIABLE(ifile, testvar_id, ndims=it)
    IF ( it == 2 ) THEN
       start = (/1,1,0/)
       count = (/nbpoint_glo,1,0/)
    ELSE IF ( it == 3 ) THEN
       start = (/1,1,1/)
       count = (/nbpoint_glo,1,1/)
    ELSE
       CALL ipslerr (3,'forcing_contfrac1d',"Testvar has a dimension larger than 3.", &
            "We do not know how to handle this", " ")
    ENDIF
    iret = NF90_GET_VAR(ifile, testvar_id, testvar, start, count)
    IF (iret /= NF90_NOERR) THEN
       WRITE(*,*) TRIM(nf90_strerror(iret))
       CALL ipslerr (3,'forcing_contfrac1d',"Error reading testvar.", " ", " ")
    ENDIF
    !
  END SUBROUTINE forcing_contfrac1d
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_contfrac2d
!!
!>\BRIEF      
!!
!! DESCRIPTION: This routine build the land/mask if needed and gets the contfrac variable from forcing file.
!!              Here we treat the case where the variables is 2D. Thus we also need to identify the land points.
!!              For this we can either use the contfrac variable or the test variable.	  
!!
!! \n
!_ ==============================================================================================================================
!
!=============================================================================================
!
  SUBROUTINE forcing_contfrac2d(ifile, testvar_id, contfrac_id, testvar, contfrac, testvar_missing, contfrac_missing, nbland)
    !
    ! This routine build the land/mask if needed and gets the contfrac variable from forcing file.
    ! Here we treat the case where the variables is 2D. Thus we also need to identify the land points.
    ! For this we can either use the contfrac variable or the test variable.
    ! 
    ! ARGUMENTS
    !
    INTEGER(i_std), INTENT(in)                 :: ifile
    INTEGER(i_std), INTENT(in)                 :: testvar_id, contfrac_id
    REAL(r_std), DIMENSION(:,:), INTENT(inout) :: testvar 
    REAL(r_std), DIMENSION(:,:), INTENT(inout) :: contfrac
    REAL(r_std), INTENT(in)                    :: testvar_missing
    REAL(r_std), INTENT(in)                    :: contfrac_missing
    INTEGER(i_std), INTENT(out)                :: nbland
    !
    ! LOCAL
    !
    INTEGER(i_std)                           :: i, j, it, iret
    INTEGER(i_std), DIMENSION(4)             :: start, count
    !
    !
    nbland = 0
    !
    IF ( contfrac_id > 0 ) THEN
       !
       iret = NF90_INQUIRE_VARIABLE(ifile, contfrac_id, ndims=it)
       IF ( it == 2 ) THEN
          start = (/1,1,0,0/)
          count = (/iim_glo,jjm_glo,0,0/)
          iret = NF90_GET_VAR(ifile, contfrac_id, contfrac, start, count)
          IF (iret /= NF90_NOERR) THEN
             WRITE(*,*) TRIM(nf90_strerror(iret))
             CALL ipslerr (3,'forcing_contfrac2d',"Error reading contfrac.", " ", " ")
          ENDIF
       ELSE
          CALL ipslerr (3,'forcing_contfrac2d',"Contfrac has a dimension different of 1.", &
               "We do not know how to handle this.", " ")
       ENDIF
       !
       ! Count the number of land points.
       !
       DO i=1,iim_glo
          DO j=1,jjm_glo
             IF ( contfrac(i,j) /= contfrac_missing ) THEN
                nbland = nbland + 1
             ENDIF
          ENDDO
       ENDDO
       !
       ! If we did not find any land points on the map (i.e. iim_glo > 1 and jjm_glo > 1) then we
       ! look for fractions larger then 0.
       !
       IF ( iim_glo > 1 .AND. jjm_glo > 1 .AND. nbland < 1 ) THEN
          DO i=1,iim_glo
             DO j=1,jjm_glo
                IF ( contfrac(i,j) > 0.0 ) THEN
                   nbland = nbland + 1
                ENDIF
             ENDDO
          ENDDO
       ENDIF
       !
       ! Did we get a result ?
       !
       IF ( iim_glo > 1 .AND. jjm_glo > 1 .AND. nbland < 1 ) THEN
          CALL ipslerr (3,'forcing_contfrac2d',"Contfrac was used to count the number of land points.", &
               &          "We still have not found any land points when we looked for contfrac > 0.", " ")
       ENDIF
       !
    ELSE
       ! Just so that we have no un-initialized variable
       contfrac(:,:) = 0.0
    ENDIF
    !
    IF ( testvar_id > 0 ) THEN
       !
       iret = NF90_INQUIRE_VARIABLE(ifile, testvar_id, ndims=it)
       IF ( it == 2 ) THEN
          start = (/1,1,0,0/)
          count = (/iim_glo,jjm_glo,0,0/)
       ELSE IF ( it == 3 ) THEN
          start = (/1,1,1,0/)
          count = (/iim_glo,jjm_glo,1,0/)
       ELSE IF ( it == 4 ) THEN
          start = (/1,1,1,1/)
          count = (/iim_glo,jjm_glo,1,1/)
       ELSE
          CALL ipslerr (3,'forcing_contfrac2d',"testvar has a dimension of 1 or larger than 4.", &
               "We do not know how to handle this.", " ")
       ENDIF
       iret = NF90_GET_VAR(ifile, testvar_id, testvar, start, count)
       IF (iret /= NF90_NOERR) THEN
          WRITE(*,*) TRIM(nf90_strerror(iret))
          CALL ipslerr (3,'forcing_contfrac2d',"Error reading testvar.", " ", " ")
       ENDIF
       !
       ! IF with count frac we did not get the number of land points, let us try it here
       !
       IF ( nbland < 1 ) THEN
          DO i=1,iim_glo
             DO j=1,jjm_glo
                IF ( testvar(i,j) /= testvar_missing ) THEN
                   nbland = nbland + 1
                   ! Add infor to contfrac
                   IF ( contfrac_id < 0 ) THEN
                      contfrac(i,j) = 1.0
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDIF
       !
       !
       ! Did we get a result here ?
       !
       IF ( iim_glo > 1 .AND. jjm_glo > 1 .AND. nbland < 1 ) THEN
          CALL ipslerr (3,'forcing_contfrac2d',"Contfrac and testvar were used to count the number", &
               &          "of land points. We have not found any land points.", " ")
       ENDIF
       !
    ENDIF
    !
  END SUBROUTINE forcing_contfrac2d
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_checkindex
!!
!>\BRIEF      
!!
!! DESCRIPTION:	 For ORCHIDEE its paralelisation requires that the land points are ordered
!!               in such a way that the longitude runs fastest. That means that we go over the
!!               globle filling one line after the other.
!!               As this might not be the case in a compressed vector of land points, we need to 
!!               put all the points on the 2d grid and then scan them in the right order. 
!!               The reindexing is prepared here.  
!!
!! \n
!_ ==============================================================================================================================
!
!=============================================================================================
!
  SUBROUTINE forcing_checkindex(dump_mask, testvarname, testvar)
    !
    ! For ORCHIDEE its paralelisation requires that the land points are ordered
    ! in such a way that the longitude runs fastest. That means that we go over the
    ! globle filling one line after the other.
    ! As this might not be the case in a compressed vector of land points, we need to 
    ! put all the points on the 2d grid and then scan them in the right order. 
    ! The reindexing is prepared here.
    !
    LOGICAL          :: dump_mask
    CHARACTER(LEN=*) :: testvarname
    REAL(r_std)      :: testvar(:)
    !
    INTEGER(i_std) :: j, i, ik
    REAL(r_std), ALLOCATABLE, DIMENSION(:)   :: testvar_reind
    !
    !
    !
    ! Get the indices of the land points in the focing file
    !
    IF ( .NOT. ALLOCATED(reindex_glo)) ALLOCATE(reindex_glo(nbpoint_glo))
    IF ( .NOT. ALLOCATED(origind)) ALLOCATE(origind(iim_glo,jjm_glo))
    !
    ! Find the origine of each point in the gathered vector on the xy grid.
    !
    origind(:,:) = -1
    mask_glo(:,:) = 0
    DO ik=1,nbpoint_glo
       j = INT((lindex_glo(ik)-1)/iim_glo)+1
       i = (lindex_glo(ik)-(j-1)*iim_glo)
       origind(i,j) = ik
       mask_glo(i,j) = 1
    ENDDO
    !
    ! Prepare a reindexing array so that the vector goes in the right order : longitude runs
    ! faster than the latitude. Put then also the right information into lindex_glo.
    !
    ik=0
    DO j=1,jjm_glo
       DO i=1,iim_glo
          IF ( origind(i,j) > zero ) THEN
             ik = ik+1
             reindex_glo(ik) = origind(i,j)
             lindex_glo(ik) = (j-1)*iim_glo+i
          ENDIF
       ENDDO
    ENDDO
    !
    !
    ! Write the mask and a test variable to a file so that the user can check that all is OK
    !
    IF ( dump_mask) THEN
       !
       ! Scatter the test variable and save it in the file
       !
       WRITE(*,*) MINVAL(testvar), "<< test variable ", TRIM(testvarname), " <<", MAXVAL(testvar) 
       ALLOCATE(testvar_reind(nbpoint_glo))
       !
       CALL forcing_reindex(nbpoint_glo, testvar, nbpoint_glo, testvar_reind, reindex_glo)
       !
       CALL forcing_writetestvar("forcing_mask_glo.nc", iim_glo, jjm_glo, nbpoint_glo, &
            &                    lon_glo(:,1), lat_glo(1,:), lindex_glo, mask_glo, &
            &                    testvarname, testvar_reind)
       !
    ENDIF
    !
    ! Clean up !
    !
    IF ( ALLOCATED(testvar_reind) ) DEALLOCATE(testvar_reind)
    !
  END SUBROUTINE forcing_checkindex
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_writetestvar
!!
!>\BRIEF Write the mask and a test variable to a netCDF file.     
!!
!! DESCRIPTION:	This routine allows to test if the variables read from the forcing files is well read.
!!              Typically the forcing is compressed by gathering and thus it is a safe practice
!!              to verify that the un-compression is done correctly and that all points end-up in the
!!              right place on the global lat/lon grid.
!!
!! \n
!_ ==============================================================================================================================
!
!=============================================================================================
!
  SUBROUTINE forcing_writetestvar(ncdffile, iim, jjm, nbland, lon, lat, lindex, mask, varname, var)
    !
    ! Write the mask and a test variable to a netCDF file
    !
    ! ARGUMENTS
    !
    CHARACTER(LEN=*), INTENT(in) :: ncdffile
    INTEGER(i_std), INTENT(in)   :: iim, jjm, nbland
    REAL(r_std), INTENT(in)      :: lon(iim), lat(jjm)
    INTEGER(i_std), INTENT(in)   :: lindex(nbland)
    INTEGER(i_std), INTENT(in)   :: mask(iim,jjm)
    CHARACTER(LEN=*), INTENT(in) :: varname
    REAL(r_std), INTENT(in)      :: var(nbland)
    !
    ! Local
    !
    INTEGER(i_std) :: ik, i, j
    INTEGER(i_std) :: iret, nlonid, nlatid, varid, fid, ierr, iland
    INTEGER(i_std) :: testid
    INTEGER(i_std), DIMENSION(2) :: lolaid
    REAL(r_std) :: test_scat(iim,jjm)
    !
    !
    test_scat(:,:) = NF90_FILL_REAL
    CALL forcing_reindex(nbland, var, iim, jjm, test_scat, lindex)
    !
    iret = NF90_CREATE(ncdffile, NF90_WRITE, fid)
    IF (iret /= NF90_NOERR) THEN
       CALL ipslerr (3,'forcing_writetestvar',"Error opening the output file : ", ncdffile, " ")
    ENDIF
    !
    ! Define dimensions
    !
    iret = NF90_DEF_DIM(fid,'lon',iim,lolaid(1))
    iret = NF90_DEF_DIM(fid,'lat',jjm,lolaid(2))
    !
    !
    iret = NF90_DEF_VAR(fid,"lon",NF90_REAL4,lolaid(1),nlonid)
    iret = NF90_PUT_ATT(fid,nlonid,'axis',"X")
    iret = NF90_PUT_ATT(fid,nlonid,'standard_name',"longitude")
    iret = NF90_PUT_ATT(fid,nlonid,'units',"degree_east")
    iret = NF90_PUT_ATT(fid,nlonid,'valid_min',MINVAL(lon_glo))
    iret = NF90_PUT_ATT(fid,nlonid,'valid_max',MAXVAL(lon_glo))
    iret = NF90_PUT_ATT(fid,nlonid,'long_name',"Longitude")
    !
    iret = NF90_DEF_VAR(fid,"lat",NF90_REAL4,lolaid(2),nlatid)
    iret = NF90_PUT_ATT(fid,nlatid,'axis',"Y")
    iret = NF90_PUT_ATT(fid,nlatid,'standard_name',"latitude")
    iret = NF90_PUT_ATT(fid,nlatid,'units',"degree_north")
    iret = NF90_PUT_ATT(fid,nlatid,'valid_min',MINVAL(lat_glo))
    iret = NF90_PUT_ATT(fid,nlatid,'valid_max',MAXVAL(lat_glo))
    iret = NF90_PUT_ATT(fid,nlatid,'long_name',"Latitude")
    !
    iret = NF90_DEF_VAR(fid,"mask",NF90_REAL4,lolaid,varid)
    !
    iret = NF90_DEF_VAR(fid,TRIM(varname),NF90_REAL4,lolaid,testid)
    iret = NF90_PUT_ATT(fid,testid,'_FillValue',NF90_FILL_REAL)
    iret = NF90_PUT_ATT(fid,testid,'missing_value',NF90_FILL_REAL)
    !
    iret = NF90_ENDDEF (fid)
    IF (iret /= NF90_NOERR) THEN
       WRITE(*,*) TRIM(nf90_strerror(iret))
       CALL ipslerr (3,'forcing_writetestvar',"Error ending definitions in file : ", ncdffile, " ")
    ENDIF
    !
    ! Write variables
    !
    iret = NF90_PUT_VAR(fid,nlonid,lon)
    iret = NF90_PUT_VAR(fid,nlatid,lat)
    iret = NF90_PUT_VAR(fid,varid,REAL(mask))
    iret = NF90_PUT_VAR(fid,testid,test_scat)
    !
    ! Close file
    !
    iret = NF90_CLOSE(fid)
    IF (iret /= NF90_NOERR) THEN
       CALL ipslerr (3,'forcing_writetestvar',"Error closing the output file : ", ncdffile, " ")
    ENDIF
    !
  END SUBROUTINE forcing_writetestvar
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_zoomgrid 
!!
!>\BRIEF      We zoom into the region requested by the user.
!!
!! DESCRIPTION:	Get the area to be zoomed and the sizes of arrays we will need.
!!              This subroutine fills the *_loc variables.
!!              If requested it will dump a test vraible into a netCDF file.  
!!
!! \n
!_ ==============================================================================================================================
!
!
!=============================================================================================
!
  SUBROUTINE forcing_zoomgrid (zoom_lon, zoom_lat, filename, dumpncdf)
    !
    ! Get the area to be zoomed and the sizes of arrays we will need.
    ! This subroutine fills the *_loc variables.
    ! If requested it will dump a test vraible into a netCDF file.
    !
    ! ARGUMENTS
    !
    REAL(r_std), DIMENSION(2), INTENT(in) :: zoom_lon, zoom_lat
    CHARACTER(LEN=*), INTENT(in) :: filename
    LOGICAL, INTENT(in) :: dumpncdf
    !
    ! LOCAL
    !
    INTEGER(i_std) :: i, j, ic, jc, ik, ig
    REAL(r_std) :: dx, dy, coslat
    REAL(r_std) :: lon_tmp(iim_glo), lat_tmp(jjm_glo)
    REAL(r_std) :: longlo_tmp(iim_glo,jjm_glo)
    REAL(r_std), ALLOCATABLE, DIMENSION(:) :: lon_val, lat_val
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:,:) :: zoom_index
    !
    INTEGER(i_std) :: iret, force_id, iv
    INTEGER(i_std), DIMENSION(1) :: imin, jmin
    INTEGER(i_std), DIMENSION(2) :: start, count
    INTEGER(i_std), DIMENSION(3) :: start2d, count2d
    REAL(r_std), ALLOCATABLE, DIMENSION(:) :: readvar, zoomedvar
     REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: readvar2d
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:) :: index_glotoloc
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: lalo
    CHARACTER(LEN=8) :: testvarname="Tair"
    !
    ! 0.0 Verify we are on the root processor
    !
    IF ( .NOT. is_root_prc ) THEN
       CALL ipslerr (3,'forcing_zoomgrid'," This routine can only be called on the root processor.", " ", " ")
    ENDIF
    !
    ! 0.1 Inform the user
    !
    WRITE(*,*) "Zoom forcing : lon = ", zoom_lon
    WRITE(*,*) "Zoom forcing : lat = ", zoom_lat
    !
    ! Some forcing files have longitudes going from 0 to 360. This code works on the
    ! -180 to 180 longitude grid. So if needed we transform the longitudes of the global grid.
    !
    IF ( MAXVAL(lon_glo) <= 180.0 ) THEN
       longlo_tmp=lon_glo
    ELSE
       DO i=1,iim_glo
          DO j=1,jjm_glo
             IF ( lon_glo(i,j) <= 180.0 ) THEN
                longlo_tmp(i,j) = lon_glo(i,j)
             ELSE
                longlo_tmp(i,j) = lon_glo(i,j)-360
             ENDIF
          ENDDO
       ENDDO
    ENDIF
    !
    ! See if we need to zoom
    !
    IF (MINVAL(zoom_lon) > MINVAL(longlo_tmp) .OR. MAXVAL(zoom_lon) < MAXVAL(longlo_tmp) .OR.&
         & MINVAL(zoom_lat) > MINVAL(lat_glo) .OR. MAXVAL(zoom_lat) < MAXVAL(lat_glo) ) THEN
       zoom_forcing = .TRUE.
    ENDIF
    !
    ! Determine the size in x and y of the zoom
    !
    IF ( zoom_forcing ) THEN
       !
       ! Working with the hypothesis it is a regular global grid and bring it back to the -180 to 180 interval
       ! if needed.
       !
       lon_tmp(:) = longlo_tmp(:,1)
       lat_tmp(:) = lat_glo(1,:)
       !
       DO i=1,iim_glo
          IF ( lon_tmp(i) <= MINVAL(zoom_lon) .OR.  lon_tmp(i) >= MAXVAL(zoom_lon) ) THEN
             lon_tmp(i) = 0.0
          ELSE
             lon_tmp(i) = 1.0
          ENDIF
       ENDDO
       DO j=1,jjm_glo
          IF ( lat_tmp(j) <= MINVAL(zoom_lat) .OR. lat_tmp(j) >= MAXVAL(zoom_lat) ) THEN
             lat_tmp(j) = 0.0
          ELSE 
             lat_tmp(j) = 1.0
          ENDIF
       ENDDO
       iim_loc = NINT(SUM(lon_tmp))
       jjm_loc = NINT(SUM(lat_tmp))
    ELSE
       iim_loc = iim_glo
       jjm_loc = jjm_glo
       lon_tmp(:) = 1.0
       lat_tmp(:) = 1.0
    ENDIF
    !
    ! Determine the number of land points in the zoom
    !
    IF ( .NOT. ALLOCATED(lon_loc) ) ALLOCATE(lon_loc(iim_loc,jjm_loc))
    IF ( .NOT. ALLOCATED(lat_loc) ) ALLOCATE(lat_loc(iim_loc,jjm_loc))
    IF ( .NOT. ALLOCATED(mask_loc) ) ALLOCATE(mask_loc(iim_loc,jjm_loc))
    IF ( .NOT. ALLOCATED(zoom_index) ) ALLOCATE(zoom_index(iim_loc,jjm_loc,2))
    !
    IF ( .NOT. ALLOCATED(lon_val)) ALLOCATE(lon_val(iim_loc))
    IF ( .NOT. ALLOCATED(lat_val)) ALLOCATE(lat_val(jjm_loc))
    !
    ! Create our new lat/lon system which is in the -180/180 range and South to North and West to East
    !
    ic=0
    DO i=1,iim_glo
       IF ( lon_tmp(i) > 0 ) THEN
          ic = ic+1
          lon_val(ic) = longlo_tmp(i,1)
       ENDIF
    ENDDO
    jc=0
    DO j=1,jjm_glo
       IF ( lat_tmp(j) > 0 ) THEN
          jc = jc+1
          lat_val(jc) = lat_glo(1,j)
       ENDIF
    ENDDO
    CALL sort(lon_val, iim_loc)
    CALL sort(lat_val, jjm_loc)
    !
    ! Now find the correspondance between the zoomed & re-ordered grid and the global one.
    !
    DO i=1,iim_loc
       DO j=1,jjm_loc
          !
          imin=MINLOC(ABS(longlo_tmp(:,1)-lon_val(i)))
          jmin=MINLOC(ABS(lat_glo(1,:)-lat_val(j)))
          !
          lon_loc(i,j) = longlo_tmp(imin(1),jmin(1))
          lat_loc(i,j) = lat_glo(imin(1),jmin(1))
          mask_loc(i,j) = mask_glo(imin(1),jmin(1))
          !
          zoom_index(i,j,1) = imin(1)
          zoom_index(i,j,2) = jmin(1)
          !
       ENDDO
    ENDDO
    !
    nbpoint_loc = SUM(mask_loc)
    IF ( .NOT. zoom_forcing .AND. nbpoint_loc .NE. nbpoint_glo) THEN
       WRITE(*,*) "We have not zoomed into the forcing file still we get a number of"
       WRITE(*,*) "land points that is different from what we have in the forcing file."
       STOP "forcing_zoomgrid"
    ENDIF
    !
    IF ( .NOT. ALLOCATED(lindex_loc)) ALLOCATE(lindex_loc(nbpoint_loc))
    IF ( .NOT. ALLOCATED(reindex_loc)) ALLOCATE(reindex_loc(nbpoint_loc))
    IF ( .NOT. ALLOCATED(contfrac_loc)) ALLOCATE(contfrac_loc(nbpoint_loc))
    !
    IF ( .NOT. ALLOCATED(reindex2d_loc)) ALLOCATE(reindex2d_loc(nbpoint_loc,2))
    IF ( .NOT. ALLOCATED(index_glotoloc)) ALLOCATE(index_glotoloc(nbpoint_glo))
    IF ( .NOT. ALLOCATED(lalo)) ALLOCATE(lalo(nbpoint_loc,2))
    !
    ! Do the actual zoom on the grid
    !
    ! Set indices of all points as non existant so that we can fill in as we zoom the 
    ! indices of the points which exist.
    index_glotoloc(:) = -1
    !
    ik = 0
    !
    ! Loop only over the zoomed grid
    !
    ! Why does the inner loop need to be ic for the pralalisation ????
    !
    DO jc=1,jjm_loc
       DO ic=1,iim_loc
          !
          ! Point back from the local to the original global i*j grid
          !
          i = zoom_index(ic,jc,1)
          j = zoom_index(ic,jc,2)
          !
          IF ( origind(i,j) > 0 ) THEN
             ik = ik+1
             ! index of the points in the local grid
             lindex_loc(ik) = (jc-1)*iim_loc+ic
             !
             ! For land points, the index of global grid is saved for the this point on the local grid 
             reindex_loc(ik) = origind(i,j)
             !
             ! Keep also the i and j of the global grid for this land point on the local grid
             reindex2d_loc(ik,1) = i
             reindex2d_loc(ik,2) = j
             !
             ! Keep the reverse : on the global grid the location where we put the value of the local grid.
             index_glotoloc(origind(i,j)) = ik
             !
             contfrac_loc(ik) = contfrac_glo(origind(i,j))
             !
             lalo(ik,1) = lat_glo(i,j)
             lalo(ik,2) = longlo_tmp(i,j)
             !
          ENDIF
       ENDDO
    ENDDO
    !
    !
    nbland_loc = 0
    DO ik=1, SIZE(contfrac_loc)
       IF (contfrac_loc(ik) > 0.0) THEN
          nbland_loc = nbland_loc + 1.0
       ENDIF
    ENDDO
    !
    !
    ncdfstart = MINVAL(reindex_loc)
    reindex_loc(:) = reindex_loc(:)-ncdfstart+1
    ncdfcount =  MAXVAL(reindex_loc)
    !
    ! Compute the areas and the corners on the grid over which we will run ORCHIDEE.
    ! As this module is dedicated for regular lat/lon forcing we know that we can compute these
    ! variables without further worries.
    !
    IF ( .NOT. ALLOCATED(area_loc)) ALLOCATE(area_loc(iim_loc,jjm_loc))
    IF ( .NOT. ALLOCATED(corners_loc)) ALLOCATE(corners_loc(iim_loc,jjm_loc,4,2))
    !
    ! Treat first the longitudes
    !
    DO j=1,jjm_loc
       dx = zero
       DO i=1,iim_loc-1
          dx = dx+ABS(lon_loc(i,j)-lon_loc(i+1,j))
       ENDDO
       dx = dx/(iim_loc-1)
       DO i=1,iim_loc
          corners_loc(i,j,1,1) = lon_loc(i,j)-dx/2.0
          corners_loc(i,j,2,1) = lon_loc(i,j)+dx/2.0
          corners_loc(i,j,3,1) = lon_loc(i,j)+dx/2.0
          corners_loc(i,j,4,1) = lon_loc(i,j)-dx/2.0
       ENDDO
    ENDDO
    !
    ! Now treat the latitudes
    !
    DO i=1,iim_loc
       dy = zero
       DO j=1,jjm_loc-1
          dy = dy + ABS(lat_loc(i,j)-lat_loc(i,j+1))
       ENDDO
       dy = dy/(jjm_loc-1)
       DO j=1,jjm_loc
          corners_loc(i,j,1,2) = lat_loc(i,j)+dy/2.0
          corners_loc(i,j,2,2) = lat_loc(i,j)+dy/2.0
          corners_loc(i,j,3,2) = lat_loc(i,j)-dy/2.0
          corners_loc(i,j,4,2) = lat_loc(i,j)-dy/2.0
       ENDDO
    ENDDO
    !
    ! Compute the areas of the grid (using the simplification that the grid is regular in lon/lat).
    !
    DO i=1,iim_loc
       DO j=1,jjm_loc
          coslat = MAX( COS(lat_loc(i,j) * pi/180. ), mincos )
          dx = ABS(corners_loc(i,j,2,1) - corners_loc(i,j,1,1)) * pi/180. * R_Earth * coslat
          dy = ABS(corners_loc(i,j,1,2) - corners_loc(i,j,3,2)) * pi/180. * R_Earth
          area_loc(i,j) = dx*dy
       ENDDO
    ENDDO
    !
    ! If requested we read a variable, zoomin and dump the result into a test file.
    !
    IF ( dumpncdf ) THEN
       iret = NF90_OPEN (filename, NF90_NOWRITE, force_id)
       IF (iret /= NF90_NOERR) THEN
          CALL ipslerr (3,'forcing_zoomgrid',"Error opening the forcing file :", filename, " ")
       ENDIF
       !
       ALLOCATE(readvar(ncdfcount), readvar2d(iim_glo,jjm_glo), zoomedvar(nbpoint_loc))
       !
       iret = NF90_INQ_VARID(force_id, TRIM(testvarname), iv)
       IF (iret /= NF90_NOERR) THEN
          CALL ipslerr (3,'forcing_zoomgrid',"Could not find variable Tair in file."," "," ")
       ENDIF

       IF ( compressed ) THEN
          !
          start(1) = ncdfstart
          start(2) = 1
          count(1) = ncdfcount
          count(2) = 1
          !
          iret = NF90_GET_VAR(force_id, iv, readvar, start, count)
          IF (iret /= NF90_NOERR) THEN
             CALL ipslerr (3,'forcing_zoomgrid',"Could not read compressed variable Tair from file."," "," ")
          ENDIF
          CALL forcing_reindex(ncdfcount, readvar, nbpoint_loc, zoomedvar, reindex_loc)
          !
       ELSE
          !
          start2d(1) = 1
          start2d(2) = 1
          start2d(3) = 1
          count2d(1) = iim_glo
          count2d(2) = jjm_glo
          count2d(3) = 1
          !
          iret = NF90_GET_VAR(force_id, iv, readvar2d, start2d, count2d)
          IF (iret /= NF90_NOERR) THEN
             CALL ipslerr (3,'forcing_zoomgrid',"Could not read variable Tair from file."," "," ")
          ENDIF
          CALL forcing_reindex(iim_glo, jjm_glo, readvar2d, nbpoint_loc, zoomedvar, reindex2d_loc)
          !
       ENDIF
       !
       CALL forcing_writetestvar("forcing_mask_loc.nc", iim_loc, jjm_loc, nbpoint_loc, &
            &                    lon_loc(:,1), lat_loc(1,:), lindex_loc, mask_loc, &
            &                    TRIM(testvarname), zoomedvar)
       !
    ENDIF
    !
    ! Clean up
    !
    IF ( ALLOCATED(readvar) ) DEALLOCATE(readvar)
    IF ( ALLOCATED(readvar2d) ) DEALLOCATE(readvar2d)
    IF ( ALLOCATED(zoomedvar) ) DEALLOCATE(zoomedvar)
    IF ( ALLOCATED(index_glotoloc) ) DEALLOCATE(index_glotoloc)
    IF ( ALLOCATED(lalo) ) DEALLOCATE(lalo)
    !
  END SUBROUTINE forcing_zoomgrid
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_givegridsize
!!
!>\BRIEF      Routine which exports the size of the grid on which the model will run, i.e. the zoomed grid.
!!
!! DESCRIPTION:	This is needed to transfer the grid information from this module to the glogrid.f90 module.  
!!
!! \n
!_ ==============================================================================================================================
!
!=============================================================================================
!
  SUBROUTINE forcing_givegridsize (iim, jjm, nblindex)
    !
    ! Provides the size of the grid to be used to the calling program
    !
    ! Size of the x and y direction of the zoomed area
    INTEGER(i_std), INTENT(out) :: iim, jjm
    ! Number of land points in the zoomed area 
    INTEGER(i_std), INTENT(out) :: nblindex
    !
    IF ( .NOT. is_root_prc ) THEN
       CALL ipslerr (3,'forcing_givegridsize'," This routine can only be called on the root processor.", &
            &          "The information requested is only available on root processor.", " ")
    ENDIF
    !
    iim = iim_loc
    jjm = jjm_loc
    nblindex = nbland_loc
    !
  END SUBROUTINE forcing_givegridsize
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_
!!
!>\BRIEF      
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================
!
!=============================================================================================
!
  SUBROUTINE forcing_vertical(force_id)
    !
    !- This subroutine explores the forcing file to decide what information is available
    !- on the vertical coordinate.
    !
    INTEGER, INTENT(IN) :: force_id
    !
    INTEGER(i_std) :: iret, ireta, iretb
    !
    INTEGER(i_std) :: sigma_id = -1, sigma_uv_id = -1
    INTEGER(i_std) :: hybsiga_id = -1, hybsiga_uv_id = -1
    INTEGER(i_std) :: hybsigb_id = -1, hybsigb_uv_id = -1
    INTEGER(i_std) :: levels_id = -1, levels_uv_id = -1
    INTEGER(i_std) :: height_id = -1, height_uv_id = -1
    INTEGER(i_std) :: lev_id = -1
    !
    LOGICAL :: var_exists, vara_exists, varb_exists, varuv_exists
    LOGICAL :: foundvar
    LOGICAL :: levlegacy
    !
    !- Set all the defaults
    !
    zfixed=.FALSE.
    zsigma=.FALSE.
    zhybrid=.FALSE.
    zlevels=.FALSE.
    zheight=.FALSE.
    zsamelev_uv = .TRUE.
    levlegacy = .FALSE.
    !
    foundvar = .FALSE.
    !
    !- We have a forcing file to explore so let us see if we find any of the conventions
    !- which allow us to find the height of T,Q,U and V.
    !
    IF ( force_id > 0 ) THEN
       !
       ! Case for sigma levels
       !
       IF ( .NOT. foundvar ) THEN
          ireta = NF90_INQ_VARID(force_id, 'Sigma', sigma_id)
          IF ( (sigma_id >= 0) .AND. (ireta == NF90_NOERR) ) THEN
             foundvar = .TRUE.
             zsigma = .TRUE.
             iretb = NF90_INQ_VARID(force_id, 'Sigma_uv', sigma_uv_id)
             IF ( (sigma_uv_id >= 0) .OR. (iretb == NF90_NOERR) ) zsamelev_uv = .FALSE.
          ENDIF
       ENDIF
       !
       ! Case for Hybrid levels
       !
       IF ( .NOT. foundvar ) THEN
          var_exists = .FALSE.
          varuv_exists = .FALSE.
          ireta = NF90_INQ_VARID(force_id, 'HybSigA', hybsiga_id)
          IF ( (hybsiga_id >= 0 ) .AND. (ireta == NF90_NOERR) ) THEN
             iretb = NF90_INQ_VARID(force_id, 'HybSigB', hybsigb_id)
             IF ( (hybsigb_id >= 0 ) .AND. (iretb == NF90_NOERR) ) THEN
                var_exists=.TRUE.
             ELSE
                CALL ipslerr ( 3, 'forcing_vertical','Missing the B coefficient for', &
                     &         'Hybrid vertical levels for T and Q','stop')
             ENDIF
          ENDIF
          ireta = NF90_INQ_VARID(force_id, 'HybSigA_uv', hybsiga_uv_id)
          IF ( (hybsiga_uv_id >= 0 ) .AND. (ireta == NF90_NOERR) ) THEN
             iretb = NF90_INQ_VARID(force_id, 'HybSigB_uv', hybsigb_uv_id)
             IF ( (hybsigb_uv_id >= 0 ) .AND. (iretb == NF90_NOERR) ) THEN
                varuv_exists=.TRUE.
             ELSE
                CALL ipslerr ( 3, 'forcing_vertical','Missing the B coefficient for', &
                     &         'Hybrid vertical levels for U and V','stop')
             ENDIF
          ENDIF
          IF ( var_exists ) THEN
             foundvar = .TRUE.
             zhybrid = .TRUE.
             IF ( varuv_exists ) zsamelev_uv = .FALSE.
          ENDIF
       ENDIF
       !
       ! Case for levels (i.e. a 2d time varying field with the height in meters)
       !
       IF ( .NOT. foundvar ) THEN
          ireta = NF90_INQ_VARID(force_id, 'Levels', levels_id)
          IF ( (levels_id >= 0 ) .AND. (ireta == NF90_NOERR) ) THEN
             foundvar = .TRUE.
             zlevels = .TRUE.
             iretb = NF90_INQ_VARID(force_id, 'Levels_uv', levels_uv_id)
             IF ( (levels_uv_id >= 0 ) .AND. (iretb == NF90_NOERR) ) zsamelev_uv = .FALSE.
          ENDIF
       ENDIF
       !
       ! Case where a fixed height is provided in meters
       !
       IF ( .NOT. foundvar ) THEN
          ireta = NF90_INQ_VARID(force_id, 'Height_Lev1', height_id)
          IF ( (height_id >= 0 ) .AND. (ireta == NF90_NOERR) ) THEN
             foundvar = .TRUE.
             zheight = .TRUE.        
             iretb = NF90_INQ_VARID(force_id, 'Height_Levuv', height_uv_id)
             IF ( (height_uv_id >= 0 ) .AND. (iretb == NF90_NOERR) ) zsamelev_uv = .FALSE.
          ENDIF
       ENDIF
       !
       ! Case where a fixed height is provided in meters in the lev variable
       !
       IF ( .NOT. foundvar ) THEN
          ireta = NF90_INQ_VARID(force_id, 'lev', lev_id)
          IF ( (lev_id >= 0 ) .AND. (ireta == NF90_NOERR) ) THEN
             foundvar = .TRUE.
             zheight = .TRUE.
             levlegacy = .TRUE.
          ENDIF
       ENDIF
       !
    ENDIF
    !
    ! We found forcing variables so we need to extract the values if we are dealing with constant values (i.e. all
    ! except the case zlevels
    !
    IF ( foundvar .AND. .NOT. zlevels ) THEN
       !
       IF ( zheight ) THEN
          !
          ! Constant height
          !
          IF ( levlegacy ) THEN
             iret = NF90_GET_VAR(force_id, lev_id, zlev_fixed)
             IF ( iret /= NF90_NOERR ) THEN
                CALL ipslerr ( 3, 'forcing_vertical','Attempted to read variable lev from forcing file in legacy mode', &
                     &         'NF90_GET_VAR failed.','stop')
             ENDIF
          ELSE
             iret = NF90_GET_VAR(force_id, height_id, zlev_fixed)
             IF ( iret /= NF90_NOERR ) THEN
                CALL ipslerr ( 3, 'forcing_vertical','Attempted to read variable Height_Lev1 from forcing file', &
                     &         'NF90_GET_VAR failed.','stop')
             ENDIF
             IF ( .NOT. zsamelev_uv ) THEN
                iret = NF90_GET_VAR(force_id, height_uv_id, zlevuv_fixed)
                IF ( iret /= NF90_NOERR ) THEN
                   CALL ipslerr ( 3, 'forcing_vertical','Attempted to read variable Height_Levuv from forcing file', &
                        &         'NF90_GET_VAR failed.','stop')
                ENDIF
             ENDIF
          ENDIF
          WRITE(*,*) "forcing_vertical : case ZLEV : Read from forcing file :", zlev_fixed, zlevuv_fixed
          !
       ELSE IF ( zsigma .OR. zhybrid ) THEN
          !
          ! Sigma or hybrid levels
          !
          IF ( zsigma ) THEN
             iret = NF90_GET_VAR(force_id, sigma_id, zhybrid_b)
             zhybrid_a = zero
             IF ( .NOT. zsamelev_uv ) THEN
                iret = NF90_GET_VAR(force_id, sigma_uv_id, zhybriduv_b)
                zhybriduv_a = zero
             ENDIF
          ELSE
             ireta = NF90_GET_VAR(force_id, hybsigb_id, zhybrid_b)
             iretb = NF90_GET_VAR(force_id, hybsiga_id, zhybrid_a)
             IF ( ireta /= NF90_NOERR .OR. iretb /= NF90_NOERR) THEN
                CALL ipslerr ( 3, 'forcing_vertical','Attempted to read variable HybSigA and HybSigB from forcing file', &
                     &         'NF90_GET_VAR failed.','stop')
             ENDIF
             IF ( .NOT. zsamelev_uv ) THEN
                ireta = NF90_GET_VAR(force_id, hybsigb_uv_id, zhybriduv_b)
                iretb = NF90_GET_VAR(force_id, hybsiga_uv_id, zhybriduv_a)
                IF ( ireta /= NF90_NOERR .OR. iretb /= NF90_NOERR) THEN
                   CALL ipslerr ( 3, 'forcing_vertical',&
                        &        'Attempted to read variable HybSigA_uv and HybSigB_uv from forcing file', &
                        &        'NF90_GET_VAR failed.','stop')
                ENDIF
             ENDIF
          ENDIF
          WRITE(*,*) "forcing_vertical : case Pressure coordinates : "
          WRITE(*,*) "Read from forcing file :", zhybrid_b, zhybrid_a, zhybriduv_b, zhybriduv_a
       ELSE
          !
          ! Why are we here ???
          !
          CALL ipslerr ( 3, 'forcing_vertical','What is the option used to describe the height of', &
               &         'the atmospheric forcing ?','Please check your forcing file.')
       ENDIF
    ENDIF
    !
    !- We have no forcing file to explore or we did not find anything. So revert back to the run.def and
    !- read what has been specified by the user.
    !
    IF ( force_id < 0 .OR. .NOT. foundvar ) THEN
       !
       !-
       !Config  Key  = HEIGHT_LEV1
       !Config  Desc = Height at which T and Q are given
       !Config  Def  = 2.0
       !Config  Help = The atmospheric variables (temperature and specific
       !Config         humidity) are measured at a specific level.
       !Config         The height of this level is needed to compute
       !Config         correctly the turbulent transfer coefficients.
       !Config         Look at the description of the forcing
       !Config         DATA for the correct value.
       !-
       zlev_fixed = 2.0
       CALL getin('HEIGHT_LEV1', zlev_fixed)
       !-
       !Config  Key  = HEIGHT_LEVW
       !Config  Desc = Height at which the wind is given
       !Config  Def  = 10.0
       !Config  Help = The height at which wind is needed to compute
       !Config         correctly the turbulent transfer coefficients.
       !-
       zlevuv_fixed = 10.0
       CALL getin('HEIGHT_LEVW', zlevuv_fixed)

       zheight = .TRUE.

       IF ( ABS(zlevuv_fixed-zlev_fixed) > EPSILON(zlev_fixed)) THEN
          zsamelev_uv = .FALSE.
       ENDIF

       CALL ipslerr ( 2, 'forcing_vertical','The height of the atmospheric forcing variables', &
            &         'was not found in the netCDF file.','Thus the values in run.def were used ... or their defaults.')
    ENDIF

  END SUBROUTINE forcing_vertical
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_givegrid
!!
!>\BRIEF      Routine which exports the grid (longitude, latitude, land indices) on which the model will run, i.e. the zoomed grid.
!!
!! DESCRIPTION:	This is needed to transfer the grid information from this module to the glogrid.f90 module.  
!!
!!
!! \n
!_ ==============================================================================================================================
!
!=============================================================================================
!
  SUBROUTINE forcing_givegrid (lon, lat, mask, area, corners, lindex, contfrac, calendar_tmp)
    !
    ! This subroutine will return to the caller the grid which has been extracted from the
    ! the forcing file. It is assumed that the caller has called forcing_givegridsize before
    ! and knows the dimensions of the fields and thus has done the correct allocations.
    !
    !
    REAL(r_std), INTENT(out) :: lon(iim_loc,jjm_loc), lat(iim_loc,jjm_loc)
    REAL(r_std), INTENT(out) :: mask(iim_loc,jjm_loc)
    REAL(r_std), INTENT(out) :: area(iim_loc,jjm_loc)
    REAL(r_std), INTENT(out) :: corners(iim_loc,jjm_loc,4,2)
    INTEGER(i_std), INTENT(out) :: lindex(nbpoint_loc)
    REAL(r_std), INTENT(out) :: contfrac(nbpoint_loc)
    CHARACTER(LEN=20), INTENT(out) :: calendar_tmp
    !
    IF ( .NOT. is_root_prc ) THEN
       CALL ipslerr (3,'forcing_givegrid'," This routine can only be called on the root processor.", &
            &          "The information requested is only available on root processor.", " ")
    ENDIF
    !
    IF (nbpoint_loc .NE. nbland_loc) THEN
       WRITE(numout, *) "forcing_givegrid:: nbpoint_loc=", nbpoint_loc
       WRITE(numout, *) "forcing_givegrid:: nbland_loc=", nbland_loc
       CALL ipslerr(3,'forcing_givegrid','nbpoint_loc and nbland_loc do not match', & 
                    'The calculation of land points is not correct', &
                    'Is your forcing file valid?')
    ENDIF
    !
    lon(:,:) = lon_loc(:,:)
    lat(:,:) = lat_loc(:,:)
    !
    mask(:,:) = mask_loc(:,:)
    area(:,:) = area_loc(:,:)
    corners(:,:,:,:) = corners_loc(:,:,:,:)
    !
    !
    lindex(:) = lindex_loc(:)
    contfrac(:) = contfrac_loc(:)
    !
    calendar_tmp = calendar
    !
  END SUBROUTINE forcing_givegrid
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_checkdim
!!
!>\BRIEF      
!!
!! DESCRIPTION:	Save the dimension or check that it is equal to the previous value.
!!              Should one of the spatial dimensions be different between 2 files, then we have a big problem.
!!              They probably do not belong to the same set of forcing files.  
!!
!! \n
!_ ==============================================================================================================================
!
!=============================================================================================
!
SUBROUTINE forcing_checkdim(ifile, filenames, out_dim, out_id, in_dim, in_id)
  !
  ! Save the dimension or check that it is equal to the previous value.
  ! Should one of the spatial dimensions be different between 2 files, then we have a big problem.
  ! They probably do not belong to the same set of forcing files.
  !
  INTEGER(i_std), INTENT(in) :: ifile
  CHARACTER(LEN=*), INTENT(in) :: filenames(:)
  INTEGER(i_std), INTENT(out) :: out_dim, out_id
  INTEGER(i_std), INTENT(in) :: in_dim, in_id
  !
  IF ( ifile == 1 ) THEN
     out_dim = in_dim
     out_id = in_id
  ELSE
     IF ( out_dim /= in_dim ) THEN
        CALL ipslerr (3,'forcing_ocheckdim', 'The dimension of the file is not the same of the first file opened.', &
             &        'The offending file is : ', filenames(ifile))
     ENDIF
  ENDIF
  !
END SUBROUTINE forcing_checkdim
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_time
!!
!>\BRIEF Read the time from each file and create the time axis to be the basis for the simulation.     
!!
!! DESCRIPTION:	This is an important routine which analyses the time axis of the forcing files and
!!              stores the main information in the SAVED variables of this routine. 
!!              As this module manages a list of forcing files we also need to check that the time
!!              axis of all these files is continuous and homogeneous.
!!              The bounds are also build for all the time axes so that we know how to interpret the
!!              various variables. 
!!
!! \n
!_ ==============================================================================================================================
!
!
!=============================================================================================
!
SUBROUTINE forcing_time(nbfiles, filenames)
  !
  ! Read the time from each file and create the time axis to be the basis
  ! for the simulation.
  !
  INTEGER(i_std) :: nbfiles
  CHARACTER(LEN=*) :: filenames(nbfiles)
  !
  INTEGER(i_std) :: iv, it, iff, tcnt, itbase, itbasetmp, ittmp
  INTEGER(i_std) :: tstart, maxtime_infile
  REAL(r_std), ALLOCATABLE, DIMENSION(:) :: timeint, time_read
  REAL(r_std), ALLOCATABLE, DIMENSION(:)   :: time_infiles
  CHARACTER(LEN=20) :: axname, calendar, timevarname
  CHARACTER(LEN=60) :: timestamp, tmpatt
  INTEGER(i_std) :: tncstart(3), tnccount(3)
  !
  INTEGER(i_std) :: iret, year0, month0, day0, hours0, minutes0, seci
  INTEGER(i_std), DIMENSION(1) :: imax, imin
  REAL(r_std) :: sec0, date_int, date0_tmp
  CHARACTER :: strc
  LOGICAL :: check=.FALSE.
  !
  ! Check that we are working on the root proc.
  !
  IF ( .NOT. is_root_prc) THEN
     CALL ipslerr (3,'forcing_time',"Cannot run this routine o other procs than root.",&
          &        "All the information on the forcing files is only lated on the root processor."," ")
  ENDIF
  !
  ! Size of unlimited dimension added up through the files. If variable not allocated before by another
  ! subroutine, it needs to be done here.
  !
  IF ( .NOT. ALLOCATED(nbtime_perfile) ) ALLOCATE(nbtime_perfile(nbfiles))
  IF ( .NOT. ALLOCATED(date0_file) ) ALLOCATE(date0_file(nbfiles,nbtax))
  !
  ! Go through all files in the list in order to get the total number of time steps we have
  ! in the nbfiles files to be read
  !
  nb_forcing_steps = 0
  maxtime_infile = 0
  DO iff=1,nbfiles
     !
     iret = NF90_INQUIRE_DIMENSION(force_id(iff), id_unlim(iff), name=axname, len=nbtime_perfile(iff))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr (3,'forcing_time',"Could not get size of dimension of unlimited axis"," "," ")
     ENDIF
     nb_forcing_steps =  nb_forcing_steps + nbtime_perfile(iff)
     IF ( nbtime_perfile(iff) > maxtime_infile ) maxtime_infile = nbtime_perfile(iff)
  ENDDO
  !
  ! Allocate the variables needed with the time length just calculated.
  ! These variables are saved in the module
  !
  ALLOCATE(time_infiles(nb_forcing_steps))
  ALLOCATE(time(nb_forcing_steps, nbtax*nbtmethods), time_bounds(nb_forcing_steps,nbtax*nbtmethods,2))
  ALLOCATE(time_axename(nbtax*nbtmethods), time_cellmethod(nbtax*nbtmethods))
  ALLOCATE(preciptime(nb_forcing_steps))
  ALLOCATE(time_sourcefile(nb_forcing_steps))
  ALLOCATE(time_id(nb_forcing_steps, nbtax))
  ! Allocate local variables
  ALLOCATE(time_read(nb_forcing_steps))
  ALLOCATE(timeint(nb_forcing_steps))
  !
  ! Get through all variables to find time_id
  ! The key variables to filled up here are time (the time stamps read in the file) and
  ! time_bounds which give the validity interval for the variables.
  !
  tstart=0
  !
  IF ( check ) WRITE(*,*) "forcing_time : going through ", nbfiles, " files to get the time."
  !
  DO iff=1,nbfiles
     !
     time_id(iff,:)=-1
     !
     ! Go through the variables in the file and find the one which is a time axis.
     !
     tcnt=1
     DO iv=1,nvars(iff)
        iret = NF90_GET_ATT(force_id(iff), iv, "units", tmpatt)
        IF ( INDEX(lowercase(tmpatt),'seconds since') > 0) THEN
           time_id(iff,tcnt)=iv
           tcnt=tcnt+1
           convtosec(iff)=1.0
        ELSE IF ( INDEX(lowercase(tmpatt),'minutes since') > 0) THEN
           time_id(iff,tcnt)=iv
           tcnt=tcnt+1
           convtosec(iff)=60.0
        ELSE IF ( INDEX(lowercase(tmpatt),'hours since') > 0) THEN
           time_id(iff,tcnt)=iv
           tcnt=tcnt+1
           convtosec(iff)=3600.0
        ENDIF
     ENDDO
     IF ( ANY(time_id(iff,:) < 0) ) THEN
        CALL ipslerr (3,'forcing_time',"Incorrect numer of time axes. A time axis is missing ",&
             &        "in file :", filenames(iff))
     ENDIF
     !
     IF ( check ) WRITE(*,*) "forcing_time : Looking at time axis for file ", force_id(iff)
     !
     ! Looping through the time axes and read them.
     !
     DO tcnt=1,nbtax
        !
        iret = NF90_INQUIRE_VARIABLE(force_id(iff), time_id(iff,tcnt), name=timevarname)
        IF ( check ) WRITE(*,*) "forcing_time : in ", iff, " found variable ", timevarname
        !
        ! Get units of time axis
        !
        iret = NF90_GET_ATT(force_id(iff), time_id(iff,tcnt), "units", timestamp) 
        IF ( check ) WRITE(*,*) "forcing_time : has time stamp ", timestamp
        !
        ! Transform the start date of the netCDF file into a julian date for the model
        !
        timestamp = TRIM(timestamp(INDEX(timestamp,'since')+6:LEN_TRIM(timestamp)))
        !
        ! Temporary fix. We need a more general method to find the right format for reading
        ! the elements of the start date.
        IF (  LEN_TRIM(timestamp) == 14 ) THEN
           READ (timestamp,'(I4.4,5(a,I1))') &
                year0, strc, month0, strc, day0, &
                strc, hours0, strc, minutes0, strc, seci
        ELSE
           READ (timestamp,'(I4.4,5(a,I2.2))') &
                year0, strc, month0, strc, day0, &
                strc, hours0, strc, minutes0, strc, seci
        ENDIF
        sec0 = hours0*3600. + minutes0*60. + seci
        CALL ymds2ju (year0, month0, day0, sec0, date0_tmp)
        date0_file(iff,tcnt) = date0_tmp
        !
        ! Now get the actual dates
        !
        tncstart(1) = 1
        tnccount(1) = nbtime_perfile(iff)
        IF ( check ) WRITE(*,*) "forcing_time : number of values read : ", tnccount(1)
        iret = NF90_GET_VAR(force_id(iff), time_id(iff,tcnt), time_read, tncstart, tnccount)
        IF (iret /= NF90_NOERR) THEN
           CALL ipslerr (3,'forcing_time',"An error occured while reading time from the file."," "," ")
        ENDIF
        !
        ! Convert the variable time from seconds since to julian days
        !
        DO it=1,nbtime_perfile(iff)
           time_infiles(tstart+it) = date0_file(iff,tcnt) + time_read(it)*convtosec(iff)/one_day
        ENDDO
        if ( check ) WRITE(*,*) "File ", iff, "goes from ",  time_infiles(tstart+1), " to ", &
             time_infiles(tstart+nbtime_perfile(iff))
        !
        ! Estimate the bounds as this information is not yet in the forcing file.
        !
        date_int = (time_infiles(tstart+nbtime_perfile(iff)) - time_infiles(tstart+1))/(nbtime_perfile(iff)-1)
        forcing_tstep_ave = date_int*one_day
        !
        ! If this is the first file then simply keep the name of the time axis. Else search the same name 
        ! in what has already been read
        !
        IF ( iff == 1 ) THEN
           itbase = (tcnt-1)*nbtmethods
           time_axename(itbase+1:itbase+4) = timevarname
           time_cellmethod(itbase+1) = "reference"
           time_cellmethod(itbase+2) = "start"
           time_cellmethod(itbase+3) = "cent"
           time_cellmethod(itbase+4) = "end"
        ELSE
           !
           ! If this is not the first file then we need to find the correct axis to add to.
           ! All information have already been saved with the first file.
           !
           DO ittmp=1,nbtax
              itbasetmp=(ittmp-1)*nbtmethods
              IF ( time_axename(itbasetmp+1) == timevarname ) THEN
                 itbase = itbasetmp
              ENDIF
           ENDDO

        ENDIF
        !
        !
        ! Keep for future usage the various information on the time axis we have just read. This time axis can
        ! be understood in 3 different ways and each variable might use a different cell method for this time
        ! axis.
        !
        ! time(:,(tcnt-1)*nbtmethods+1) : corresponds to the reference time axis as it has been read from the file
        ! time(:,(tcnt-1)*nbtmethods+2) : is the time axis with a cell method which place the value at the 
        !                                beginning of the time interval
        ! time(:,(tcnt-1)*nbtmethods+3) : is the time axis corresponding to variables placed at the center of the 
        !                                time interval
        ! time(:,(tcnt-1)*nbtmethods+4) : for variables put at the end of the time interval over which they aere 
        !                                for instance averaged.
        !
        ! In variable time_cellmethod we will write the type of cell method as descirbed above so that the selection
        ! of the right axis for each variable can be made automaticaly.
        !
        ! We also keep the name of the time axis read in preparation of file where we might have to read more than one
        ! time axis.
        !
        DO it=tstart+1,tstart+nbtime_perfile(iff)
           !
           ! Reference time
           !
           time(it,itbase+1) = time_infiles(it)
           time_bounds(it,itbase+1,1) = time_infiles(it)-date_int/2.0
           time_bounds(it,itbase+1,2) = time_infiles(it)+date_int/2.0
           !
           ! Start cell method
           time(it,itbase+2) = time_infiles(it)+date_int/2.0
           time_bounds(it,itbase+2,1) = time_infiles(it)
           time_bounds(it,itbase+2,2) = time_infiles(it)+date_int
           !
           ! Centered cell method
           time(it,itbase+3) = time_infiles(it)
           time_bounds(it,itbase+3,1) = time_infiles(it)-date_int/2.0
           time_bounds(it,itbase+3,2) = time_infiles(it)+date_int/2.0
           !
           ! End cell method
           time(it,itbase+4) = time_infiles(it)-date_int/2.0
           time_bounds(it,itbase+4,1) = time_infiles(it)-date_int
           time_bounds(it,itbase+4,2) = time_infiles(it)
           !
        ENDDO
        !
        ! Keep the number of the file from which we read this time.
        !
        time_sourcefile(tstart+1:tstart+nbtime_perfile(iff))=iff
        !
        IF ( check ) WRITE(*,*) "forcing_time : finished file ", iff
        !
     ENDDO
     !
     ! Before moving to the next file advance the pointer in the time arrays.
     !
     tstart=tstart+nbtime_perfile(iff)
     !
  ENDDO
  !
  IF ( check ) WRITE(*,*) "forcing_time : All files have been processed"
  !
  ! Verify that the forcing comes in regular time intervals. If not, many of the 
  ! interpolation schemes will fail.
  ! This is only done on the first time axis ... is it enough ?
  !
  DO ittmp=1,nbtax
     itbase=(ittmp-1)*nbtmethods
     !
     date_int = (time(nb_forcing_steps,itbase+1) - time(1,itbase+1))/(nb_forcing_steps-1)
     forcing_tstep_ave = date_int*one_day
     !
     timeint(:) = 0
     DO it=1, nb_forcing_steps-1
        timeint(it) = time(it+1,itbase+1)-time(it,itbase+1)
     ENDDO
     !
     IF (  MAXVAL(timeint(1:nb_forcing_steps-1)) > date_int+0.1*date_int .OR.&
          & MINVAL(timeint(1:nb_forcing_steps-1)) < date_int-0.1*date_int) THEN
        WRITE(*,*) "The time steping of the forcing files does not seem to be regular on axis",time_axename(itbase+1),":"
        WRITE(*,*) "Average time step : ", date_int, "days = ", date_int*one_day, "sec."
        imax = MAXLOC(timeint(1:nb_forcing_steps-1))
        imin = MINLOC(timeint(1:nb_forcing_steps-1))
        WRITE(*,*) "Maximum time step : ", MAXVAL(timeint(1:nb_forcing_steps-1)), " at ", imax(1)
        WRITE(*,*) "Minimum time step : ", MINVAL(timeint(1:nb_forcing_steps-1)), " at ", imin(1)
        WRITE(*,*) "++++ Values around Maximum"
        DO it=MAX(imax(1)-5,1),MIN(imax(1)+5,nb_forcing_steps)
           WRITE(*,*) it, " from file ", time_sourcefile(it), " Value ", time(it,itbase+1)
           CALL forcing_printdate(time(it,itbase+1), "Time values.")
        ENDDO
        WRITE(*,*) "++++ Values around Minimum"
        DO it=MAX(imin(1)-5,1),MIN(imin(1)+5,nb_forcing_steps)
           WRITE(*,*) it, " from file ", time_sourcefile(it), " Value ", time(it,itbase+1)
           CALL forcing_printdate(time(it,itbase+1), "Time values.")
        ENDDO
        CALL ipslerr (3,'forcing_time', 'The time handling could be improved to allow the driver',&
             & "to cope with irregular forcing files."," ")
     ENDIF
  ENDDO
  !
  ! Print some test values
  !
  DO ittmp=1,nbtax
     itbase=(ittmp-1)*nbtmethods
     !
     WRITE(*,*) "Bounds for axis ",time_axename(itbase+1)," :"
     !
     CALL forcing_printdate(time_bounds(1,itbase+1,1), "Start time of first forcing interval.")
     CALL forcing_printdate(time_bounds(1,itbase+1,2), "End time of first forcing interval.")
     CALL forcing_printdate(time_bounds(nb_forcing_steps,itbase+1,1), "Start time of last forcing interval.")
     CALL forcing_printdate(time_bounds(nb_forcing_steps,itbase+1,2), "End time of last forcing interval.")
  ENDDO
  !
  ! Set to zero the variable for the cummulated time for rainfall
  !
  preciptime(:) = zero
  !
  ! Keep the very first date of the time axis for future use
  ! 
  forcingstartdate = time(1,1)
  !
  ! Clean-up
  !
  DEALLOCATE(timeint, time_read)
  !
END SUBROUTINE forcing_time
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_varforslab
!!
!>\BRIEF      
!!
!! DESCRIPTION:	This subroutine will read the named variable and put it in the right place in the 
!!              slab of data kept in the memory of the driver. 
!!
!! \n
!_ ==============================================================================================================================
!
!=============================================================================================
!
SUBROUTINE forcing_varforslab(fileindex, varname, timestart, timecount, inslabpos, data, cellmethod)
  !
  ! This subroutine will read the named variable and put it in the right place in the 
  ! slab of data kept in the memory of the driver.
  !
  INTEGER(i_std), INTENT(in) :: fileindex
  CHARACTER(LEN=*), INTENT(in) :: varname
  INTEGER(i_std), INTENT(in) :: timestart, timecount, inslabpos
  REAL(r_std), INTENT(inout) :: data(nbpoint_loc,slab_size)
  CHARACTER(LEN=*), INTENT(out) :: cellmethod
  !
  ! Local variables
  !
  INTEGER(i_std) :: varid, windid, windndims, it, ig, iv
  INTEGER(i_std) :: iret, ndims
  INTEGER(i_std), DIMENSION(4) :: start, count
  REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: tmp_slab
  REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:) :: tmp_slab2d
  CHARACTER(LEN=80) :: name
  LOGICAL :: windzero
  !
  ! Allocate the temporary data array if not already available
  !
  IF ( compressed ) THEN
     IF ( .NOT. ALLOCATED(tmp_slab) ) ALLOCATE(tmp_slab(ncdfcount,slab_size))
  ELSE
     IF ( .NOT. ALLOCATED(tmp_slab2d) ) ALLOCATE(tmp_slab2d(iim_glo,jjm_glo,slab_size))
  ENDIF
  !
  ! Reset the counters and flags to forget the past !
  !
  varid=-1
  windid=-1
  windzero=.FALSE.
  !
  ! Find the variable in the file
  !
  DO iv=1,nvars(fileindex)
     !
     iret = NF90_INQUIRE_VARIABLE(force_id(fileindex), iv, name=name, ndims=it)
     !
     IF ( INDEX(name, varname) > 0 ) THEN
        varid = iv
        ndims = it
     ENDIF
     IF ( (INDEX(name, "Wind") > 0) .AND. (LEN_TRIM(name) == LEN_TRIM("Wind")) ) THEN
        windid = iv
        windndims = it
     ENDIF
     !
  ENDDO
  !
  ! Treat some special cases and catch errors
  !
  IF ( varid < 0 ) THEN
     !
     ! If we requested a wind component but did not find it, it could be that only the module is available.
     ! If that is the case, then we use the module (windid) for the U component and set V top zero.
     !
     IF ( INDEX(varname, "Wind_E") > 0 ) THEN
        varid = windid
        ndims = windndims
        windzero = .FALSE.
     ELSE IF ( INDEX(varname, "Wind_N") > 0 ) THEN
        windzero = .TRUE.
     ELSE
        CALL ipslerr (3,'forcing_varforslab',"Could not find variable",varname," in file.")
     ENDIF
  ENDIF
  !
  ! If there is some variable to be read then do it
  !
  IF ( .NOT. windzero ) THEN
     !
     ! Get the attributes needed for intepretating the data
     !
     ! First get the cell method used for this variable
     iret = NF90_GET_ATT(force_id(fileindex), varid, 'cell_methods', cellmethod)
     IF (iret /= NF90_NOERR) THEN
        ! If the attribute is not found then we set a reasonable default : instantaneous and centered.
        cellmethod="time: instantaneous"
     ENDIF
     !
     !
     ! Getsize of data to be read from the netCDF file
     !
     !
     IF ( compressed ) THEN
        !
        IF ( ndims == 2 ) THEN
           start = (/ncdfstart,timestart,0,0/)
           count = (/ncdfcount,timecount,0,0/)
        ELSE IF ( ndims == 3 ) THEN
           start = (/ncdfstart,1,timestart,0/)
           count = (/ncdfcount,1,timecount,0/)
        ELSE
           CALL ipslerr (3,'forcing_varforslab',"Compressed variable : ",varname,&
                &        " does not have a compatible number of dimensions.")
        ENDIF
        !
        iret = NF90_GET_VAR(force_id(fileindex), varid, tmp_slab, start, count)
        IF (iret /= NF90_NOERR) THEN
           CALL ipslerr (3,'forcing_varforslab',"Could not read from file variable : ",varname," Compressed in the file.")
        ENDIF
        !
        ! Zoom into the data and put it in the right place in the slab of data.
        !
        CALL forcing_reindex(ncdfcount, timecount, tmp_slab, nbpoint_loc, slab_size, data, inslabpos, reindex_loc)
     ELSE
        !
        IF ( ndims == 3 ) THEN
           start = (/1,1,timestart,0/)
           count = (/iim_glo,jjm_glo,timecount,0/)
        ELSE IF (ndims == 4 ) THEN
           start = (/1,1,1,timestart/)
           count = (/iim_glo,jjm_glo,1,timecount/) 
       ELSE
           CALL ipslerr (3,'forcing_varforslab',"Full lat Lon variable : ",varname,&
                &        " does not have a compatible number of dimensions.")
        ENDIF
        !
        iret = NF90_GET_VAR(force_id(fileindex), varid, tmp_slab2d, start, count)
        IF (iret /= NF90_NOERR) THEN
           WRITE(*,*) TRIM(NF90_STRERROR(iret))
           WRITE(*,*) "File =", fileindex, "Size =", SIZE(tmp_slab2d,DIM=1), SIZE(tmp_slab2d,DIM=2), SIZE(tmp_slab2d,DIM=3)
           WRITE(*,*) "Start :", start(1:3)
           WRITE(*,*) "Count :", count(1:3)
           CALL ipslerr (3,'forcing_varforslab',"Could not read from file variable : ",varname," Not compressed.")
        ENDIF
        !
        ! Zoom into the data and put it in the right place in the slab of data.
        !
        CALL forcing_reindex(iim_glo, jjm_glo, timecount, tmp_slab2d, nbpoint_loc, slab_size, data, inslabpos, reindex2d_loc)
     ENDIF
  ELSE
     cellmethod="time: instantaneous"
     DO it=0,timecount-1
        data(:,inslabpos+it) = zero
     ENDDO
  ENDIF
  !
END SUBROUTINE forcing_varforslab
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_attributetimeaxe
!!
!>\BRIEF  Find the time axis which corresponds to the variable at hand.    
!!
!! DESCRIPTION:	 We interpret the cell_method provided in the netCDF file so that
!!               we can determine how we need to understand the values we read.
!!
!! \n
!_ ==============================================================================================================================
!
!=============================================================================================
!
  SUBROUTINE forcing_attributetimeaxe(cellmethod, timeindex)
  !
  ! We will analyse the time axis of the cell method found in the NetCDF file in order to
  ! attribute the correct time axis to this variable. 
  !
  CHARACTER(LEN=*), INTENT(in) :: cellmethod
  INTEGER(i_std), INTENT(out)  :: timeindex
  !
  INTEGER(i_std) :: itax, timepos, pos, lentime, itbase, im
  CHARACTER(LEN=20) :: TARGET, tmpstr
  CHARACTER(LEN=80) :: method
  !
  ! Clean the string to delete spaces in front of ":" and "("
  !
  method = cellmethod
  DO WHILE ( INDEX(method," :") > 0 )
     pos = INDEX(method," :")
     method = method(1:pos-1)//method(pos+1:LEN_TRIM(method))
  ENDDO
  DO WHILE ( INDEX(method,"( ") > 0 )
     pos = INDEX(method,"( ")
     method = method(1:pos)//method(pos+2:LEN_TRIM(method))
  ENDDO
  !
  ! Go through all the time axes we have to find the right one.
  !
  timeindex=0
  DO itax=1,nbtax
     !
     itbase=(itax-1)*nbtmethods
     ! find the time axis name in the cell method
     TARGET = TRIM(time_axename(itbase+1))//":"
     timepos = INDEX(method,TRIM(TARGET))
     !
     ! If we found the time axis then we look for the method with a time position description
     ! which is expected to be between parenthesis. For instance : mean(end)
     !
     IF ( timepos > 0 ) THEN
        !
        lentime=LEN_TRIM(time_axename(itbase+1))
        tmpstr = method(lentime+2:LEN_TRIM(method))
        !
        ! If there is ":" then there is information for another axis which needs to be deleted
        !
        IF ( INDEX(tmpstr,":") > 0 ) THEN
           tmpstr = tmpstr(1:INDEX(tmpstr,":")-1)
        ENDIF
        !
        ! Now that we have found a time axis see if we have between parenthesis a position
        ! on that time avis.
        !
        ! Look for a "("  
        IF ( INDEX(tmpstr, "(") > 0 ) THEN
           DO im=1,nbtmethods
              TARGET = "("//TRIM(time_cellmethod(itbase+im))
              timepos = INDEX(tmpstr,TRIM(TARGET))
              IF ( timepos > 0 ) THEN
                 timeindex = itbase+im
              ENDIF
           ENDDO
           !
           ! If there is no "(" then we have to find the centered axis.
        ELSE 
           DO im=1,nbtmethods
              IF ( INDEX(time_cellmethod(itbase+im), "cent") > 0 ) THEN
                 timeindex = itbase+im
              ENDIF
           ENDDO
        ENDIF
        !
        ! The name of the time axis was found bu no method could be identified
        !
        IF ( timeindex < 1 ) THEN
           CALL ipslerr (3,'forcing_attributetimeaxe',"Found a time axis name but could not identify method.", &
                "This should not happen !"," ")
        ENDIF
        !
     ELSE
        ! Continue in loop over nbtax
     ENDIF
  ENDDO
  !
  ! Should no corresponding time axis name be found, 
  ! then we use the first centered one.
  !
  itax=1
  DO WHILE ( timeindex < 1 ) 
     IF ( INDEX(time_cellmethod(itax), "cent") > 0 ) THEN
        timeindex = itax
     ELSE
        itax = itax + 1
     ENDIF
  ENDDO
  !
END SUBROUTINE forcing_attributetimeaxe
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_filenamecheck
!!
!>\BRIEF   Check the list of files obtained from the calling program.   
!!
!! DESCRIPTION:	A small code to check the forcing files. They have to be NetCDF (i.e. .nc termination) and 
!!              we dont keep files that appear more than once in the list.  
!!
!! \n
!_ ==============================================================================================================================
!!
SUBROUTINE forcing_filenamecheck(filenames_in, nb_files)
  !
  ! A small code to check the forcing files. They have to
  ! be NetCDF (i.e. .nc termination) and we dont keep files 
  ! that appear more than once in the list.
  !
  !
  ! INPUT
  !
  CHARACTER(LEN=*), DIMENSION(:), INTENT(in) :: filenames_in
  INTEGER(i_std), INTENT(out)                :: nb_files
  !
  ! LOCAL
  !
  INTEGER(i_std) :: ii, is, sizein
  LOGICAL        :: notfound
  !
  sizein = SIZE(filenames_in)
  IF ( sizein > 0 ) THEN
     IF ( ALLOCATED(forfilename) ) THEN
        DEALLOCATE(forfilename)
     ENDIF
     ALLOCATE(forfilename(sizein))
     nb_files=0
  ELSE
     CALL ipslerr (3,'forcing_filenamecheck',"List of forcing files is empty.","Please check your run.def file."," ")
  ENDIF
  !
  DO ii=1,sizein
     IF ( INDEX(filenames_in(ii), '.nc') > 0 ) THEN
        IF ( nb_files == 0 ) THEN
           nb_files = nb_files+1
           forfilename(nb_files)=TRIM(ADJUSTL(filenames_in(ii)))
        ELSE
           notfound=.TRUE.
           DO is=1,nb_files
              IF ( INDEX(TRIM(filenames_in(ii)), TRIM(ADJUSTL(forfilename(is)))) > 0 ) notfound=.FALSE.
           ENDDO
           IF ( notfound ) THEN
              nb_files = nb_files+1
              forfilename(nb_files)=TRIM(adjustl(filenames_in(ii)))
           ENDIF
        ENDIF
     ELSE
        !!! This is not a NetCDF file, so we ignore it
     ENDIF
  ENDDO
  !
  !
END SUBROUTINE forcing_filenamecheck
!!
!!  =============================================================================================================================
!! SUBROUTINE: lowercase, FindMinimum, Swap 
!!
!>\BRIEF      Help functions for the forcing_tools module.
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================
!
!====================================================================================================
!
!
FUNCTION lowercase(strIn) RESULT(strOut)
! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)

     IMPLICIT NONE

     CHARACTER(len=*), INTENT(in) :: strIn
     CHARACTER(len=LEN(strIn)) :: strOut
     INTEGER :: i,j

     DO i = 1, LEN(strIn)
          j = IACHAR(strIn(i:i))
          IF (j>= IACHAR("A") .AND. j<=IACHAR("Z") ) THEN
               strOut(i:i) = ACHAR(IACHAR(strIn(i:i))+32)
          ELSE
               strOut(i:i) = strIn(i:i)
          END IF
     END DO

END FUNCTION lowercase
!
! Some help function found on Internet : http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap08/sorting.f90
!
! --------------------------------------------------------------------
! INTEGER FUNCTION  FindMinimum():
!    This function returns the location of the minimum in the section
! between Start and End.
! --------------------------------------------------------------------
INTEGER FUNCTION  FindMinimum(x, Start, END)
  IMPLICIT  NONE
  REAL(r_std), DIMENSION(1:), INTENT(IN) :: x
  INTEGER(i_std), INTENT(IN)             :: Start, END
  REAL(r_std)                            :: Minimum
  INTEGER(i_std)                         :: Location
  INTEGER(i_std)                         :: i

  Minimum  = x(Start)		! assume the first is the min
  Location = Start			! record its position
  DO i = Start+1, END		! start with next elements
     IF (x(i) < Minimum) THEN	!   if x(i) less than the min?
        Minimum  = x(i)		!      Yes, a new minimum found
        Location = i            !      record its position
     ENDIF
  END DO
  FindMinimum = Location        	! return the position
END FUNCTION  FindMinimum
! --------------------------------------------------------------------
! SUBROUTINE  Swap():
!    This subroutine swaps the values of its two formal arguments.
! --------------------------------------------------------------------
SUBROUTINE  Swap(a, b)
  IMPLICIT  NONE
  REAL(r_std), INTENT(INOUT) :: a, b
  REAL(r_std)                :: Temp

  Temp = a
  a    = b
  b    = Temp
END SUBROUTINE  Swap
! --------------------------------------------------------------------
! SUBROUTINE  Sort():
!    This subroutine receives an array x() and sorts it into ascending
! order.
! --------------------------------------------------------------------
SUBROUTINE  Sort(x, Size)
  IMPLICIT  NONE
  REAL(r_std), DIMENSION(1:), INTENT(INOUT) :: x
  INTEGER(i_std), INTENT(IN)                :: Size
  INTEGER(i_std)                            :: i
  INTEGER(i_std)                            :: Location

  DO i = 1, Size-1			! except for the last
     Location = FindMinimum(x, i, Size)	! find min from this to last
     CALL  Swap(x(i), x(Location))	! swap this and the minimum
  END DO
END SUBROUTINE  Sort
!

END MODULE forcing_tools
