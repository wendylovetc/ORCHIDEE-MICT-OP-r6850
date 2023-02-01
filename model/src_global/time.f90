! ================================================================================================================================
!  MODULE       : time
!
!  CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
!  LICENCE      : IPSL (2006)
!  This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF   This module contians the time information for ORCHIDEE
!!
!!\n DESCRIPTION: This module contains the time information for ORCHIDEE.
!! Time variables are calculated and updated at each time step. The variables are public and can be used from other modules.
!! The time variables must not be modified outside this module. 
!!
!! Time variables in this module are given for the interval corresponding to current time-step. 
!! The variables xxx_start and xxx_end are given to describe the time interval, either on julian(julian_start, julian_end) or 
!! yymmddsec format (year_start/end, month_start/end, day_start/end, sec_start/end). The variables can be used in all other modules. 
!! Logical variables telling if it's the last or first time-step of the day, of the month or of the year are also available.  
!!
!! Subroutines in module time : 
!!   time_initialize : To initialize time variables. 
!!                     Public subroutine called from intersurf_initialize_2d, intersurf_initialize_gathered or orchideedriver.
!!                     This subroutine is called before the first call to sechiba.
!!   time_nextstep   : Update time variables in the beginning at each new time step.
!!                     Public subroutine called from intersurf_main_2d, intersurf_main_gathered or orchideedriver.
!!   time_interval   : Calculate the interval for a given time step. Local subroutine, called from time_nextstep. 
!!
!! REFERENCE(S)	: None
!!
!! SVN          :
!! $HeadURL: $
!! $Date: $
!! $Revision: $
!! \n
!_ ================================================================================================================================
MODULE time

  USE defprec
  USE constantes_var
  USE mod_orchidee_para_var, ONLY : numout
!  USE ioipsl
  USE ioipsl_para
  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: time_initialize, time_nextstep, time_get_day_of_year
!           time_is_gregorian_leap_year, time_days_in_gregorian_month, time_get_day_of_year

  REAL(r_std), PUBLIC                 :: dt_sechiba         !! Lenght of time step in sechiba (s)
!$OMP THREADPRIVATE(dt_sechiba)
  REAL(r_std), PUBLIC                 :: dt_stomate         !! Length of time step in slow processes in stomate (s)
!$OMP THREADPRIVATE(dt_stomate)
  CHARACTER(LEN=20), PUBLIC           :: calendar_str       !! The calendar type
!$OMP THREADPRIVATE(calendar_str)
  INTEGER(i_std), SAVE, PUBLIC        :: year_start         !! Year at beginng of current time step
!$OMP THREADPRIVATE(year_start)
  INTEGER(i_std), SAVE, PUBLIC        :: month_start        !! Month at beginng of current time step
!$OMP THREADPRIVATE(month_start)
  INTEGER(i_std), SAVE, PUBLIC        :: day_start          !! Day in month at beginng of current time step
!$OMP THREADPRIVATE(day_start)
  REAL(r_std), SAVE, PUBLIC           :: sec_start          !! Seconds in the day at beginng of current time step
!$OMP THREADPRIVATE(sec_start)
  INTEGER(i_std), SAVE, PUBLIC        :: year_end           !! Year at end of current time step
!$OMP THREADPRIVATE(year_end)
  INTEGER(i_std), SAVE, PUBLIC        :: month_end          !! Month at end of current time step
!$OMP THREADPRIVATE(month_end)
  INTEGER(i_std), SAVE, PUBLIC        :: day_end            !! Day in month at end of current time step
!$OMP THREADPRIVATE(day_end)
  REAL(r_std), SAVE, PUBLIC           :: sec_end            !! Seconds in the day at end of current time step
!$OMP THREADPRIVATE(sec_end)
  CHARACTER(LEN=6), SAVE, PUBLIC      :: tstepint_type      !! Position of time step in the time interval
!$OMP THREADPRIVATE(tstepint_type)
  INTEGER(i_std), PUBLIC              :: month_len          !! Lenght of current month (d)
!$OMP THREADPRIVATE(month_len)
  INTEGER(i_std), PRIVATE             :: year_length        !! Lenght of current year in time step (tstp)
!$OMP THREADPRIVATE(year_length)
  INTEGER(i_std), PUBLIC              :: year_length_in_days!! Length of current year in days
!$OMP THREADPRIVATE(year_length_in_days)
  REAL(r_std), PUBLIC                 :: one_day            !! Lenght of one day in seconds (s)
!$OMP THREADPRIVATE(one_day)
  REAL(r_std), PUBLIC                 :: one_year           !! Length of current year in days (d)
!$OMP THREADPRIVATE(one_year)
  REAL(r_std), PARAMETER, PUBLIC      :: one_hour = 3600.0  !! Lenght of hour in seconds (s)  

  LOGICAL, PUBLIC                     :: FirstTsYear        !! Flag is true for the first sechiba time step on the year.
!$OMP THREADPRIVATE(FirstTsYear)
  LOGICAL, PUBLIC                     :: LastTsYear         !! Flag is true for the last sechiba time step of the year, previously named EndOfYear
!$OMP THREADPRIVATE(LastTsYear)
  LOGICAL, PUBLIC                     :: FirstTsMonth       !! Flag is true for the first sechiba time step of the month. 
!$OMP THREADPRIVATE(FirstTsMonth)
  LOGICAL, PUBLIC                     :: LastTsMonth        !! Flag is true for the last sechiba time step of the month. 
!$OMP THREADPRIVATE(LastTsMonth)
  LOGICAL, PUBLIC                     :: FirstTsDay         !! Flag is true for the first sechiba time step of the day. 
!$OMP THREADPRIVATE(FirstTsDay)
  LOGICAL, PUBLIC                     :: LastTsDay          !! Flag is true for the last sechiba time step of the day. 
!$OMP THREADPRIVATE(LastTsDay)
  REAL(r_std), SAVE, PUBLIC           :: date0_save         !! Start date of simulation, in juilan calendar
!$OMP THREADPRIVATE(date0_save)
  REAL(r_std), PUBLIC                 :: julian_diff        !! Days since the beginning of current year
!$OMP THREADPRIVATE(julian_diff)
  REAL(r_std), PUBLIC                 :: julian_start       !! Beginning of the interval for current time step, in juilan calendar
!$OMP THREADPRIVATE(julian_start)
  REAL(r_std), PUBLIC                 :: julian_end         !! End of the interval for current time step, in juilan calendar
!$OMP THREADPRIVATE(julian_end)
  REAL(r_std), SAVE, PUBLIC           :: julian0            !! First day of this year in julian caledrier
!$OMP THREADPRIVATE(julian0)
  INTEGER(i_std), SAVE, PRIVATE       :: printlev_loc       !! Local level of text output for current module
!$OMP THREADPRIVATE(printlev_loc)

  INTEGER(i_std), SAVE, PUBLIC        :: day_of_year        !! Current day of the year
!$OMP THREADPRIVATE(day_of_year)
  INTEGER(i_std), SAVE, PUBLIC        :: year0              !! First year of the simulation 
!$OMP THREADPRIVATE(year0)
  INTEGER(i_std), SAVE, PUBLIC        :: month0             !! First month of the simulation 
!$OMP THREADPRIVATE(month0)
  INTEGER(i_std), SAVE, PUBLIC        :: day0               !! First day of the simulation 
!$OMP THREADPRIVATE(day0)
  REAL(r_std), SAVE, PUBLIC           :: sec0               !! First sec of the simulation 
!$OMP THREADPRIVATE(sec0)

CONTAINS

!!  =============================================================================================================================
!! SUBROUTINE:    time_initialize()
!!
!>\BRIEF	  Initalize time information
!!
!! DESCRIPTION:	  Initialize time information. This subroutine is called only in the intialization phase of the model and gives
!!                the basic information on the calendar and time interval.
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE time_initialize(kjit, date0_loc, dt_sechiba_loc, tstepint)

    !! 0.1 Input arguments
    INTEGER(i_std), INTENT(in)   :: kjit           !! Time step of the restart file
    REAL(r_std), INTENT(in)      :: date0_loc      !! The date at which kjit=0
    REAL(r_std), INTENT(in)      :: dt_sechiba_loc !! Time step of sechiba component
    CHARACTER(LEN=*), INTENT(in) :: tstepint       !! Position of the time stamp: beginning, centre or end of time step 
                                                   !! Possible values for tstepint : START, CENTER, END
    REAL(r_std) :: julian0_start

    !! Initialize local printlev variable
    !! It is not possible here to use the function get_printlev for problems with circular dependecies between modules
    printlev_loc=printlev
    CALL getin_p('PRINTLEV_time', printlev_loc)

    !! Save length of sechiba time step in module variable
    !! Time step for sechiba comes from the atmospheric model or from the 
    !! offline driver which reads it from parameter DT_SECHIBA
    dt_sechiba = dt_sechiba_loc

    !! Save tstepint in global variable
    tstepint_type = tstepint
    IF (tstepint_type /= "START" .AND. tstepint_type /= "END") THEN
       WRITE(numout,*) 'Unknown option for time interval, tstepint_type=', tstepint_type
       CALL ipslerr_p(3,'time_initialize', 'Unknown time iterval type.',&
            'The time stamp given can have 2 position on the interval :',' START or END')
    END IF

    !! Save the start date in the module
    date0_save = date0_loc

    !! Get the calendar from IOIPSL, it is already initialized
    CALL ioget_calendar(calendar_str)

    !! Get year lenght in days and day lenght in seconds
    CALL ioget_calendar(one_year, one_day)

    !! Get the number of sechiba time steps for 1 year (for the current year)
    CALL tlen2itau('1Y', dt_sechiba, date0_save, year_length)

    !Config Key   = DT_STOMATE
    !Config Desc  = Time step of STOMATE and other slow processes
    !Config If    = OK_STOMATE
    !Config Def   = one_day
    !Config Help  = Time step (s) of regular update of vegetation
    !Config         cover, LAI etc. This is also the time step of STOMATE.
    !Config Units = [seconds]
    dt_stomate = one_day
    CALL getin_p('DT_STOMATE', dt_stomate)

    !! Find starting year of the simulation
    julian0_start = itau2date(kjit, date0_save, dt_sechiba)
    CALL ju2ymds (julian0_start, year0, month0, day0, sec0)
 
    IF (printlev_loc >=1) THEN
       WRITE(numout,*) "time_initialize : calendar_str= ",calendar_str
       WRITE(numout,*) "time_initialize : dt_sechiba(s)= ",dt_sechiba," dt_stomate(s)= ", dt_stomate
       WRITE(numout,*) "time_initialize : date0_save= ",date0_save," kjit= ",kjit
       WRITE(numout,*) "time_initialize : year0, month0, day0, sec0 =", &
                        year0, month0, day0, sec0
    ENDIF
    
    CALL time_nextstep(kjit)
    
  END SUBROUTINE time_initialize
  
!!  =============================================================================================================================
!! SUBROUTINE:    time_nextstep()
!!
!>\BRIEF	  Update the time information for the next time step.
!!
!! DESCRIPTION:	  This subroutine will place in the public variables of the module all the time information
!!                needed by ORCHIDEE for this new time step.
!!
!! \n
!_ ==============================================================================================================================
  SUBROUTINE time_nextstep(kjit)
    
    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)   :: kjit       !! Current time step


    !! Calculate the new time interval for current time step
    CALL time_interval(kjit, year_start, month_start, day_start, sec_start, &
                             year_end,   month_end,   day_end,   sec_end)
    
    !! Calculate julian0: the julian date for the 1st of January in current year
    CALL ymds2ju(year_start,1,1,zero, julian0)

    !! Caluclate julian_diff: The diffrence of the end of current time-step and the beginning of the year in julian days
    julian_diff = julian_end - julian0

    !! Update variables for year length
    CALL ioget_calendar(one_year, one_day)
    CALL tlen2itau('1Y', dt_sechiba, date0_save, year_length)

    !! Calculate number of days in the current month, using the start of the time-interval for current time-step
    month_len = ioget_mon_len(year_start,month_start)

    !! Calculate current day of the year
    day_of_year = time_get_day_of_year(year_start, month_start, day_start)

    !! Calculate the length of the current year in days
    year_length_in_days = ioget_year_len(year_start)

    !! Calculate logical variables true if the current time-step corresponds to the end or beggining of a day, month or year
    FirstTsDay = .FALSE.
    FirstTsMonth = .FALSE.
    FirstTsYear = .FALSE.
    LastTsDay = .FALSE.
    LastTsMonth = .FALSE.
    LastTsYear = .FALSE.
 
    IF (sec_start >= -1 .AND. sec_start < dt_sechiba-1 ) THEN
       FirstTsDay = .TRUE.
       IF ( day_start == 1 ) THEN
          FirstTsMonth = .TRUE.
          IF ( month_start == 1) THEN
             FirstTsYear = .TRUE.
          END IF
       END IF
    ELSE IF (sec_start >= one_day-dt_sechiba-1 .AND. sec_start < one_day-1 ) THEN
       LastTsDay = .TRUE.
       IF ( day_start == month_len ) THEN
          LastTsMonth = .TRUE.
          IF ( month_start == 12) THEN
             LastTsYear = .TRUE.
          END IF
       END IF
    END IF

    !! Write debug information depending on printlev_loc.
    IF ( ((printlev_loc >= 4) .OR. (printlev_loc >=2 .AND. FirstTsMonth)) .OR. (printlev_loc>=1 .AND. FirstTsYear)) THEN
       WRITE(numout,*) "time_nextstep: Time interval for time step :", kjit
       WRITE(numout,"(' time_nextstep: Start of interval         ', I4.4,'-',I2.2,'-',I2.2,' ',F12.4)") &
            year_start, month_start, day_start, sec_start
       WRITE(numout,"(' time_nextstep: End of interval           ', I4.4,'-',I2.2,'-',I2.2,' ',F12.4)") &
            year_end, month_end, day_end, sec_end
       WRITE(numout,*) "time_nextstep: month_len=",month_len, " one_year=",one_year
       WRITE(numout,*) ""
    END IF
    IF (printlev_loc >= 4) THEN
       WRITE(numout,*) "time_nextstep: FirstTsDay=", FirstTsDay," FirstTsMonth=",FirstTsMonth," FirstTsYear=", FirstTsYear
       WRITE(numout,*) "time_nextstep: LastTsDay=", LastTsDay," LastTsMonth=",LastTsMonth," LastTsYear=", LastTsYear
       WRITE(numout,*) "time_nextstep: julian0=",julian0, " julian_diff=",julian_diff
       WRITE(numout,*) "time_nextstep: day_of_year=", day_of_year 
       WRITE(numout,*) ""
    END IF
    IF ((printlev_loc >= 2) .AND. FirstTsDay) THEN
       WRITE(numout,"(' Date: ', I4.4,'-',I2.2,'-',I2.2,' ',F8.4,'sec at timestep ',I12.4)") &
            year_start, month_start, day_start, sec_start, kjit
    END IF

  END SUBROUTINE time_nextstep

!!  =============================================================================================================================
!! SUBROUTINE:    time_interval()
!!
!>\BRIEF	  Computes the interval corresponging to the given time step.
!!
!! DESCRIPTION:	  This subroutine will compute the interval of time for the given time step.
!!                It will use the tstepint_type variable which was set at initilisation. 
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE time_interval(kjit, year_s, month_s, day_s, sec_s, &
                                 year_e, month_e, day_e, sec_e)

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                  :: kjit      !! Current time step

    !! 0.2 Output variables
    INTEGER(i_std), INTENT(out)                 :: year_s, month_s, day_s
    INTEGER(i_std), INTENT(out)                 :: year_e, month_e, day_e
    REAL(r_std), INTENT(out)                    :: sec_s
    REAL(r_std), INTENT(out)                    :: sec_e
    

    !! Calculate the interval for current time step
    !! tstepint_type is used to know how the interval is defined around the time-stamp of the time-step 
    SELECT CASE (TRIM(tstepint_type))
       CASE ("START ")
          !! Calculate the start date of the interval
          julian_start = itau2date(kjit, date0_save, dt_sechiba)
          !! Calculate the end date of the interval using kjti+1
          julian_end = itau2date(kjit+1, date0_save, dt_sechiba)
       CASE("END")
          !! Calculate the start date of the interval
          julian_start = itau2date(kjit-1, date0_save, dt_sechiba)
          !! Calculate the end date of the interval
          julian_end = itau2date(kjit, date0_save, dt_sechiba)
       CASE DEFAULT
          WRITE(numout,*) "time_interval: tstepint_type = ", tstepint_type
          CALL ipslerr_p(3,'time_interval', 'Unknown time iterval type.', &
               'The time stamp given can have 3 position on the interval :',' START, CENTER or END')
       END SELECT

       !! Calculate year, month, day and sec for julian_start date
       CALL ju2ymds (julian_start, year_s, month_s, day_s, sec_s)

       !! Calculate year, month, day and sec for julian_end date
       CALL ju2ymds (julian_end, year_e, month_e, day_e, sec_e)

    
  END SUBROUTINE time_interval

!!  =============================================================================================================================
!! SUBROUTINE: time_get_day_of_year
!!
!>\BRIEF	  Calculates day of the year
!!
!! DESCRIPTION:	Calculates day of the year for a given date
!!   Algorithm:
!!   http://www.herongyang.com/year/Program-Gregorian-Calendar-Algorithm.html 
!!                
!!
!! \n
!_ ==============================================================================================================================
  INTEGER FUNCTION time_get_day_of_year(year, month, day)
    ! 0.1 Input variables
    INTEGER(i_std), INTENT(in)  :: year, month, day 

    ! 0.4 Local variables
    INTEGER(i_std)              :: calc_days, item
    INTEGER(i_std)              :: days_in_month

!    WRITE(*,*) "time_get_day_of_year:: input, year, month, day==", year, month, day

!    ! Check for consistency
!    IF (day .GT. days_in_month)
!      CALL ipslerr_p(3, 'time_get_day_of_year','Incorrent date', &
!                    'Max number of days expected:'//days_in_month , &
!                    'But found:'//day )
!    ENDIF


    calc_days = 0
    DO item = 1, month-1
      calc_days=calc_days + ioget_mon_len(year, item) 
    ENDDO

    time_get_day_of_year = calc_days + day
  
  END FUNCTION time_get_day_of_year
 
END MODULE time
