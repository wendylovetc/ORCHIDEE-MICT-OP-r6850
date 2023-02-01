MODULE orchoasis_tools
  !
  USE defprec
  USE netcdf
  !
  USE ioipsl_para
  USE mod_orchidee_para
  !
#ifdef OASIS
  USE mod_oasis
  !
  ! ORCHIDEE definitions
  !
  USE constantes
  USE constantes_soil
  USE constantes_mtc
  USE pft_parameters
  !
  USE control
  !
  IMPLICIT NONE
  !
  PRIVATE
  PUBLIC :: orchoasis_time, &
       &    orchoasis_printpoint, orchoasis_defvar, orchoasis_getvar, &
       &    orchoasis_putvar, orchoasis_printdate
  !
  !
  !  Global variables
  !
  !
  INTEGER(i_std), SAVE                               :: itau_offset  !! This offset is used to phase the
  ! INPUT
  INTEGER(i_std), SAVE :: il_part_id, tair_id, qair_id, zlevtq_id, zlevuv_id
  INTEGER(i_std), SAVE :: rainf_id, snowf_id, lwdown_id, swnet_id, solarang_id
  INTEGER(i_std), SAVE :: u_id, v_id, ps_id, cdrag_id
  ! OUTPUT
  INTEGER(i_std), SAVE :: vevapp_id, fluxsens_id, fluxlat_id, coastal_id, river_id
  INTEGER(i_std), SAVE :: netco2_id, carblu_id, tsolrad_id, tsolnew_id, qsurf_id
  INTEGER(i_std), SAVE :: albnir_id, albvis_id, emis_id, z0_id
  !
  LOGICAL, PARAMETER :: check_INPUTS = .FALSE.         !! (very) long print of INPUTs in intersurf 
  LOGICAL, SAVE :: check_time = .FALSE.
  !
  LOGICAL, SAVE :: landonly = .TRUE.
  !
  PUBLIC check_time
  !
CONTAINS
!
!=============================================================================================
!
  SUBROUTINE orchoasis_time(date_start, dt, nbdt)
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
       CALL ipslerr(3, "orchoasis_time", "START_DATE incorrectly specified in run.def", str_sdate(1), str_sdate(2))
    ENDIF
    CALL ymds2ju (s_year, s_month, s_day, s_sec, date_start)
    CALL orchoasis_printdate(date_start, "This is after reading the start date")
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
       CALL ipslerr(3, "orchoasis_time", "END_DATE incorrectly specified in run.def", str_edate(1), str_edate(2))
    ENDIF
    CALL ymds2ju (e_year, e_month, e_day, e_sec, date_end)
    CALL orchoasis_printdate(date_start, "This is after reading the end date")
    !
    CALL time_diff (s_year,s_month,s_day,s_sec,e_year,e_month,e_day,e_sec,diff_sec)
    !
    !Config Key  = DT_SECHIBA
    !Config Desc = Time step length in seconds for Sechiba component
    !Config Def  = 1800
    !Config Help =
    !Config Units = [seconds] 
    dt = 1800
    CALL getin('DT_SECHIBA', dt)
    !
    nbdt = NINT(diff_sec/dt)
    !
  END SUBROUTINE orchoasis_time
!
!=============================================================================================
!
  SUBROUTINE orchoasis_printpoint(julian_day, lon_pt, lat_pt, nbind, lalo, var, message, ktest)
    !
    REAL(r_std), INTENT(in) :: julian_day
    REAL(r_std), INTENT(in) :: lon_pt, lat_pt
    INTEGER(i_std), INTENT(in) :: nbind
    REAL(r_std), INTENT(in) :: lalo(:,:)
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
       dist(i) = acos( sin(lat_pt*pi/180)*sin(lalo(i,1)*pi/180) + &
            &    cos(lat_pt*pi/180)*cos(lalo(i,1)*pi/180)*cos((lalo(i,2)-lon_pt)*pi/180) ) * R_Earth
    ENDDO
    !
    ! Look for the next grid point closest to the one with the smalest distance.
    !
    imin = MINLOC(dist)
    DO i=1,nbind
       refdist(i) = acos( sin(lalo(imin(1),1)*pi/180)*sin(lalo(i,1)*pi/180) + &
            &       cos(lalo(imin(1),1)*pi/180)*cos(lalo(i,1)*pi/180) * cos((lalo(i,2)-lalo(imin(1),2))*pi/180) ) * R_Earth
    ENDDO
    refdist(imin(1)) =  MAXVAL(refdist)
    mindist = MINVAL(refdist)
    !
    ! Are we closer than the closest points ?
    !
    IF ( PRESENT(ktest) ) ktest = -1
    IF ( dist(imin(1)) <= mindist ) THEN
       !
       WRITE(numout,"(I2.2,':',I2.2,':',I2.2,' Loc : ', F6.1,',', F6.1,'(i=',I6,') Value = ',F12.4,A38)") &
            & hours, minutes, seci, lalo(imin(1),2), lalo(imin(1),1), imin(1), var(imin(1)), message
       !
       IF ( PRESENT(ktest) ) ktest = imin(1)
    ENDIF
    !
  END SUBROUTINE orchoasis_printpoint
  !
  !-------------------------------------------------------------------------------------------------------
  !-
  !- orchoasis_defvar
  !-
  !-------------------------------------------------------------------------------------------------------
  SUBROUTINE orchoasis_defvar (mpi_rank, kjpindex)
    !
    ! ARGUMENTS
    !
    INTEGER(i_std), INTENT(in) :: mpi_rank, kjpindex
    !
    ! LOCAL
    !
    INTEGER(i_std)               :: ierror, tot_error
    INTEGER(i_std), DIMENSION(3) :: ig_paral
    INTEGER(i_std), DIMENSION(2) :: var_nodims, var_shape
    CHARACTER(LEN=8) :: varname
    !
    !
    !Config Key  = COUPLING_LANDONLY
    !Config Desc = Specifies if the OASIS coupler will provide only data on land or on land and ocean.
    !Config Def  = TRUE
    !Config Help = If COUPLING_LANDONLY is set to TRUE, then ORCHIDEE expects to receive only land points
    !              from OASIS. On the contrary data on both land and ocean will be received and the arrays
    !              re-indexed so as to keep only land points.
    landonly = .TRUE.
    CALL getin('COUPLING_LANDONLY',landonly)
    !
    tot_error = 0
    !
    IF ( landonly ) THEN
       ig_paral(1) = 1
       ig_paral(2) = nbp_mpi_para_begin(mpi_rank)-1
       ig_paral(3) = kjpindex
    ELSE
       ig_paral(1) = 1
       ig_paral(2) = ij_para_begin(mpi_rank) - 1
       ig_paral(3) = ij_nb
    ENDIF

    WRITE(*,*) mpi_rank, "Start and length =",  ig_paral(2:3)

    CALL oasis_def_partition (il_part_id, ig_paral, ierror)

    IF ( landonly ) THEN
       var_nodims(1) = 1
       var_nodims(2) = 1
       var_shape(1) = 1
       var_shape(1) = kjpindex
    ELSE
       var_nodims(1) = 1
       var_nodims(2) = 1
       var_shape(1) = 1
       var_shape(1) = ij_nb
    ENDIF
    !
    ! ORCHIDEE's Input variables
    ! ===========================
    !
    !
    ! Define levels
    !
    varname="HEIGHTTQ"
    CALL oasis_def_var(zlevtq_id, varname, il_part_id, var_nodims, OASIS_In, var_shape, OASIS_Real, ierror)
    tot_error = tot_error+ierror
    !
    varname="HEIGHTUV"
    CALL oasis_def_var(zlevuv_id, varname, il_part_id, var_nodims, OASIS_In, var_shape, OASIS_Real, ierror)
    tot_error = tot_error+ierror
    !
    !
    varname="TEMPLEV1"
    CALL oasis_def_var(tair_id, varname, il_part_id, var_nodims, OASIS_In, var_shape, OASIS_Real, ierror)
    tot_error = tot_error+ierror
    varname="HUMILEV1"
    CALL oasis_def_var(qair_id, varname, il_part_id, var_nodims, OASIS_In, var_shape, OASIS_Real, ierror)
    tot_error = tot_error+ierror
    !
    ! Define precipitation variables
    !
    varname="RAINFALL"
    CALL oasis_def_var(rainf_id, varname, il_part_id, var_nodims, OASIS_In, var_shape, OASIS_Real, ierror)
    tot_error = tot_error+ierror
    varname="SNOWFALL"
    CALL oasis_def_var(snowf_id, varname, il_part_id, var_nodims, OASIS_In, var_shape, OASIS_Real, ierror)
    tot_error = tot_error+ierror
    !
    !
    ! Define precipitation variables
    !
    varname="SHORTNET"
    CALL oasis_def_var(swnet_id, varname, il_part_id, var_nodims, OASIS_In, var_shape, OASIS_Real, ierror)
    tot_error = tot_error+ierror
    varname="LONWDOWN"
    CALL oasis_def_var(lwdown_id, varname, il_part_id, var_nodims, OASIS_In, var_shape, OASIS_Real, ierror)
    tot_error = tot_error+ierror
    varname="SOLARANG"
    CALL oasis_def_var(solarang_id, varname, il_part_id, var_nodims, OASIS_In, var_shape, OASIS_Real, ierror)
    tot_error = tot_error+ierror
    !
    ! Define dynamical variables
    !
    varname="EASTWIND"
    CALL oasis_def_var(u_id, varname, il_part_id, var_nodims, OASIS_In, var_shape, OASIS_Real, ierror)
    tot_error = tot_error+ierror
    varname="NORTWIND"
    CALL oasis_def_var(v_id, varname, il_part_id, var_nodims, OASIS_In, var_shape, OASIS_Real, ierror)
    tot_error = tot_error+ierror
    varname="SURFPRES"
    CALL oasis_def_var(ps_id, varname, il_part_id, var_nodims, OASIS_In, var_shape, OASIS_Real, ierror)
    tot_error = tot_error+ierror
    varname="MODCDRAG"
    CALL oasis_def_var(cdrag_id, varname, il_part_id, var_nodims, OASIS_In, var_shape, OASIS_Real, ierror)
    tot_error = tot_error+ierror
    !
    ! ORCHIDEE's Output variables
    ! ===========================
    !
    ! Turbulent fluxes
    !
    varname="TOTEVAPS"
    CALL oasis_def_var(vevapp_id, varname, il_part_id, var_nodims, OASIS_Out, var_shape, OASIS_Real, ierror)
    tot_error = tot_error+ierror
    !
    varname="FLUXSENS"
    CALL oasis_def_var(fluxsens_id, varname, il_part_id, var_nodims, OASIS_Out, var_shape, OASIS_Real, ierror)
    tot_error = tot_error+ierror
    !
    varname="FLUXLATE"
    CALL oasis_def_var(fluxlat_id, varname, il_part_id, var_nodims, OASIS_Out, var_shape, OASIS_Real, ierror)
    tot_error = tot_error+ierror
    !
    ! Discharge to the oceans
    !
    varname="COASTFLO"
    CALL oasis_def_var(coastal_id, varname, il_part_id, var_nodims, OASIS_Out, var_shape, OASIS_Real, ierror)
    tot_error = tot_error+ierror
    !
    varname="RIVERFLO"
    CALL oasis_def_var(river_id, varname, il_part_id, var_nodims, OASIS_Out, var_shape, OASIS_Real, ierror)
    tot_error = tot_error+ierror
    !
    ! Carbon fluxes
    !
    varname="FLUNECO2"
    CALL oasis_def_var(netco2_id, varname, il_part_id, var_nodims, OASIS_Out, var_shape, OASIS_Real, ierror)
    tot_error = tot_error+ierror
    !
    varname="FLULUCO2"
    CALL oasis_def_var(carblu_id, varname, il_part_id, var_nodims, OASIS_Out, var_shape, OASIS_Real, ierror)
    tot_error = tot_error+ierror
    !
    ! Surface states
    !
    varname="TSURFRAD"
    CALL oasis_def_var(tsolrad_id, varname, il_part_id, var_nodims, OASIS_Out, var_shape, OASIS_Real, ierror)
    tot_error = tot_error+ierror
    !
    varname="TSURFNEW"
    CALL oasis_def_var(tsolnew_id, varname, il_part_id, var_nodims, OASIS_Out, var_shape, OASIS_Real, ierror)
    tot_error = tot_error+ierror
    !
    varname="QSURFNEW"
    CALL oasis_def_var(qsurf_id, varname, il_part_id, var_nodims, OASIS_Out, var_shape, OASIS_Real, ierror)
    tot_error = tot_error+ierror
    !
    varname="ALBEDVIS"
    CALL oasis_def_var(albvis_id, varname, il_part_id, var_nodims, OASIS_Out, var_shape, OASIS_Real, ierror)
    tot_error = tot_error+ierror
    !
    varname="ALBEDNIR"
    CALL oasis_def_var(albnir_id, varname, il_part_id, var_nodims, OASIS_Out, var_shape, OASIS_Real, ierror)
    tot_error = tot_error+ierror
    !
    varname="EMISLONW"
    CALL oasis_def_var(emis_id, varname, il_part_id, var_nodims, OASIS_Out, var_shape, OASIS_Real, ierror)
    tot_error = tot_error+ierror
    !
    varname="ROUGNESS"
    CALL oasis_def_var(z0_id, varname, il_part_id, var_nodims, OASIS_Out, var_shape, OASIS_Real, ierror)
    tot_error = tot_error+ierror
    !
    IF ( tot_error .NE. 0 ) THEN
       CALL ipslerr(3, "orchoasis_defvar", "Definition of one of the coupling variables failed.", &
            &          "No futher information can be given at this point.", "")
    ENDIF
    !
    CALL oasis_enddef (ierror)
    !
    IF ( ierror .NE. 0 ) THEN
       CALL ipslerr(3, "orchoasis_defvar", "End of the definition of coupling variables failed.", &
            &          "If OASIS has a way to decide the error status it should be put here.", "")
    ENDIF
    !
  END SUBROUTINE orchoasis_defvar
  !-------------------------------------------------------------------------------------------------------
  !-
  !- orchoasis_getvar
  !-
  !-------------------------------------------------------------------------------------------------------
  SUBROUTINE orchoasis_getvar (itau, dt, kjpindex, landindex, zlev_tq, zlev_uv, temp_air, qair, precip_rain, &
       &                      precip_snow, swnet, lwdown, sinang, u, v, pb, cdrag)
    !
    ! ARGUMENTS
    !
    INTEGER(i_std), INTENT(in) :: itau, kjpindex
    INTEGER(i_std), DIMENSION(kjpindex), INTENT(in) :: landindex
    REAL(r_std), INTENT(in)    :: dt
    !
    REAL(r_std), DIMENSION(kjpindex), INTENT(out) :: zlev_tq, zlev_uv, temp_air, qair, precip_rain
    REAL(r_std), DIMENSION(kjpindex), INTENT(out) :: precip_snow, swnet, lwdown, sinang, u, v, pb
    REAL(r_std), DIMENSION(kjpindex), INTENT(out) :: cdrag
    !
    ! LOCAL
    !
    INTEGER(i_std) :: oasis_info
    REAL(r_std), ALLOCATABLE, DIMENSION(:) :: buffer
    !
    ! Dimension the buffer for getting the data from OASIS according to the number of
    ! points we will get : with or without ocean points.
    !
    IF ( .NOT. ALLOCATED(buffer) ) THEN
       IF ( landonly ) THEN
          ALLOCATE(buffer(kjpindex))
       ELSE
          ALLOCATE(buffer(ij_nb))
       ENDIF
    ENDIF
    !
    ! Get first the levels
    !
    CALL oasis_get(zlevtq_id, NINT(itau*dt), buffer, oasis_info)
    IF ( landonly ) THEN
       zlev_tq = buffer
    ELSE
       zlev_tq = buffer(landindex)
    ENDIF
    CALL oasis_get(zlevuv_id, NINT(itau*dt), buffer, oasis_info)
    IF ( landonly ) THEN
       zlev_uv = buffer
    ELSE
       zlev_uv = buffer(landindex)
    ENDIF
    !
    ! Get atmospheric state variables
    !
    CALL oasis_get(tair_id, NINT(itau*dt), buffer, oasis_info)
    IF ( landonly ) THEN
       temp_air = buffer
    ELSE
       temp_air = buffer(landindex)
    ENDIF
    CALL oasis_get(qair_id, NINT(itau*dt), buffer, oasis_info)
    IF ( landonly ) THEN
       qair = buffer
    ELSE
       qair = buffer(landindex)
    ENDIF
    !
    ! Get precipitation fluxes
    !
    CALL oasis_get(rainf_id, NINT(itau*dt), buffer, oasis_info)
    IF ( landonly ) THEN
       precip_rain = buffer
    ELSE
       precip_rain = buffer(landindex)
    ENDIF
    CALL oasis_get(snowf_id, NINT(itau*dt), buffer, oasis_info)
    IF ( landonly ) THEN
       precip_snow = buffer
    ELSE
       precip_snow = buffer(landindex)
    ENDIF
    !
    ! Get Radiation fluxes
    !
    CALL oasis_get(swnet_id, NINT(itau*dt), buffer, oasis_info)
    IF ( landonly ) THEN
       swnet = buffer
    ELSE
       swnet = buffer(landindex)
    ENDIF
    CALL oasis_get(lwdown_id, NINT(itau*dt), buffer, oasis_info)
    IF ( landonly ) THEN
       lwdown = buffer
    ELSE
       lwdown = buffer(landindex)
    ENDIF
    CALL oasis_get(solarang_id, NINT(itau*dt), buffer, oasis_info)
    IF ( landonly ) THEN
       sinang = buffer
    ELSE
       sinang = buffer(landindex)
    ENDIF
    !
    ! Get dynamical variables
    !
    CALL oasis_get(u_id, NINT(itau*dt), buffer, oasis_info)
    IF ( landonly ) THEN
       u = buffer
    ELSE
       u = buffer(landindex)
    ENDIF
    CALL oasis_get(v_id, NINT(itau*dt), buffer, oasis_info)
    IF ( landonly ) THEN
       v = buffer
    ELSE
       v = buffer(landindex)
    ENDIF
    CALL oasis_get(ps_id, NINT(itau*dt), buffer, oasis_info)
    IF ( landonly ) THEN
       pb = buffer
    ELSE
       pb = buffer(landindex)
    ENDIF
    CALL oasis_get(cdrag_id, NINT(itau*dt), buffer, oasis_info)
    IF ( landonly ) THEN
       cdrag = buffer
    ELSE
       cdrag = buffer(landindex)
    ENDIF
    !
  END SUBROUTINE orchoasis_getvar
  !
  !-------------------------------------------------------------------------------------------------------
  !-
  !- orchoasis_putvar
  !-
  !-------------------------------------------------------------------------------------------------------
  SUBROUTINE orchoasis_putvar(itau, dt, kjpindex, landindex, vevapp, fluxsens, fluxlat, coastalflow, riverflow, &
       &                      netco2, carblu, tsol_rad, temp_sol_new, qsurf, albedo, emis, z0)
    !
    ! ARGUMENTS
    !
    INTEGER(i_std), INTENT(in) :: itau, kjpindex
    INTEGER(i_std), DIMENSION(kjpindex), INTENT(in) :: landindex
    REAL(r_std), INTENT(in)    :: dt
    !
    REAL(r_std), DIMENSION(kjpindex), INTENT(in) :: vevapp, fluxsens, fluxlat, coastalflow, riverflow
    REAL(r_std), DIMENSION(kjpindex), INTENT(in) :: netco2, carblu, tsol_rad, temp_sol_new, qsurf, emis, z0
    REAL(r_std), DIMENSION(kjpindex,2), INTENT(in) :: albedo
    !
    ! LOCAL
    !
    INTEGER(i_std) :: oasis_info
    REAL(r_std), ALLOCATABLE, DIMENSION(:) :: buffer
    !
    IF ( .NOT. ALLOCATED(buffer) ) THEN
       IF ( landonly ) THEN
          ALLOCATE(buffer(kjpindex))
       ELSE
          ALLOCATE(buffer(ij_nb))
       ENDIF
    ENDIF
    !
    ! Set all points to undef. This is the value which will remain over ocean points.
    !
    buffer(:) = undef_sechiba
    !
    ! Turbulent fluxes
    !
    IF ( landonly ) THEN
       buffer(:) = vevapp(:)
    ELSE
       buffer(landindex(:)) = vevapp(:)
    ENDIF
    CALL oasis_put(vevapp_id, NINT(itau*dt), buffer, oasis_info)
    IF ( landonly ) THEN
       buffer(:) = fluxsens(:)
    ELSE
       buffer(landindex(:)) = fluxsens(:)
    ENDIF
    CALL oasis_put(fluxsens_id, NINT(itau*dt), buffer, oasis_info)
    IF ( landonly ) THEN
       buffer(:) = fluxlat(:)
    ELSE
       buffer(landindex(:)) = fluxlat(:)
    ENDIF
    CALL oasis_put(fluxlat_id, NINT(itau*dt), buffer, oasis_info)
    !
    ! Water fluxes to the ocean
    !
    IF ( landonly ) THEN
       buffer(:) = coastalflow(:)
    ELSE
       buffer(landindex(:)) = coastalflow(:)
    ENDIF
    CALL oasis_put(coastal_id, NINT(itau*dt), buffer, oasis_info)
    IF ( landonly ) THEN
       buffer(:) = riverflow(:)
    ELSE
       buffer(landindex(:)) = riverflow(:)
    ENDIF
    CALL oasis_put(river_id, NINT(itau*dt), buffer, oasis_info)
    !
    ! Carbon
    !
    IF ( landonly ) THEN
       buffer(:) = netco2(:)
    ELSE
       buffer(landindex(:)) = netco2(:)
    ENDIF
    CALL oasis_put(netco2_id, NINT(itau*dt), buffer, oasis_info)
    IF ( landonly ) THEN
       buffer(:) = carblu(:)
    ELSE
       buffer(landindex(:)) = carblu(:)
    ENDIF
    CALL oasis_put(carblu_id, NINT(itau*dt), buffer, oasis_info)
    !
    ! Surface conditions
    !
    IF ( landonly ) THEN
       buffer(:) = tsol_rad(:)
    ELSE
       buffer(landindex(:)) = tsol_rad(:)
    ENDIF
    CALL oasis_put(tsolrad_id, NINT(itau*dt), buffer, oasis_info)
    IF ( landonly ) THEN
       buffer(:) = temp_sol_new(:)
    ELSE
       buffer(landindex(:)) = temp_sol_new(:)
    ENDIF
    CALL oasis_put(tsolnew_id, NINT(itau*dt), buffer, oasis_info)
    IF ( landonly ) THEN
       buffer(:) = qsurf(:)
    ELSE
       buffer(landindex(:)) = qsurf(:)
    ENDIF
    CALL oasis_put(qsurf_id, NINT(itau*dt), buffer, oasis_info)
    !
    ! Other surface conditions
    !
    IF ( landonly ) THEN
       buffer(:) = albedo(:,ivis)
    ELSE
       buffer(landindex(:)) = albedo(:,ivis)
    ENDIF
    CALL oasis_put(albvis_id, NINT(itau*dt), buffer, oasis_info)
    IF ( landonly ) THEN
       buffer(:) = albedo(:,inir)
    ELSE
       buffer(landindex(:)) = albedo(:,inir)
    ENDIF
    CALL oasis_put(albnir_id, NINT(itau*dt), buffer, oasis_info)
    IF ( landonly ) THEN
       buffer(:) = emis(:)
    ELSE
       buffer(landindex(:)) = emis(:)
    ENDIF
    CALL oasis_put(emis_id, NINT(itau*dt), buffer, oasis_info)
    IF ( landonly ) THEN
       buffer(:) = z0(:)
    ELSE
       buffer(landindex(:)) = z0(:)
    ENDIF
    CALL oasis_put(z0_id, NINT(itau*dt), buffer, oasis_info)
    !
  END SUBROUTINE orchoasis_putvar
!
!=============================================================================================
!
  SUBROUTINE orchoasis_printdate(julian_day, message)
    !
    REAL(r_std), INTENT(in) :: julian_day
    CHARACTER(len=*), INTENT(in) :: message
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
    WRITE(*,'(I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2," > ", A60)') &
         &            year, month, day, hours, minutes, seci, message
    !
  END SUBROUTINE orchoasis_printdate
  !
  !
#endif
END MODULE orchoasis_tools
