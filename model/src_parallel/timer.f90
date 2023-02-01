! ==============================================================================================================================
! MODULE   : timer
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF     Timer functions to calculate MPI time consumption.
!!
!!\n DESCRIPTION  : Timer functions to calculate MPI time consumption (global time and mpi time). 
!!               - store in timer_state the state of the simulation (running, stopped, suspended)
!!               - calcutate the cpu time in cpu_timer / the real time in real_timer
!!               
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCES(S)    : None
!!
!! SVN              :
!! $HeadURL$
!! $Date$
!! $Revision$
!! \n
!_ ================================================================================================================================
MODULE timer

  USE mod_orchidee_para_var, ONLY : numout
  
  INTEGER, PARAMETER :: nb_timer=2
  INTEGER, PARAMETER :: timer_global=1
  INTEGER, PARAMETER :: timer_mpi=2
  INTEGER, PARAMETER :: stopped = 1
  INTEGER, PARAMETER :: running = 2
  INTEGER, PARAMETER :: suspended = 3
  
  DOUBLE PRECISION, DIMENSION(nb_timer),SAVE :: cpu_timer
  DOUBLE PRECISION, DIMENSION(nb_timer),SAVE :: real_timer
  INTEGER, DIMENSION(nb_timer),SAVE :: timer_state
  DOUBLE PRECISION, DIMENSION(nb_timer),SAVE :: last_cpu_time
  INTEGER, DIMENSION(nb_timer),SAVE :: last_real_time
  
  
  
  
  CONTAINS
  
  !!  =============================================================================================================================
  !! SUBROUTINE:  init_timer
  !!
  !>\BRIEF   Initialization of the timer. 
  !!
  !! DESCRIPTION:  Initialization of the timer. Need to be called at the beginning of 
  !!               the simulation by the master threads OMP on each process mpi. 
  !!
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE init_timer
  IMPLICIT NONE
    
    cpu_timer(:)=0.
    real_timer(:)=0.
    timer_state(:)=stopped
    last_cpu_time(:)=0.
    last_real_time(:)=0
    
  END SUBROUTINE init_timer
  
  
  !!  =============================================================================================================================
  !! SUBROUTINE:  start_timer
  !!
  !>\BRIEF	Start all timer variables for the type of timer which was choosen (global timer or mpi timer) 
  !!
  !! DESCRIPTION:	Start all timer variables for the type of timer which was choosen (global timer or mpi timer) 
  !!                   Need to be call by the master threads on each process mpi. 	 
  !!
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE start_timer(no_timer)
  IMPLICIT NONE
     INTEGER :: no_timer
     DOUBLE PRECISION :: x
     
     IF (timer_state(no_timer)/=stopped) THEN
       STOP 'start_timer :: timer is already running or suspended'
     ELSE
        timer_state(no_timer)=running
     ENDIF
      
     cpu_timer(no_timer)=0. 
     real_timer(no_timer)=0.
     x=Diff_real_time(no_timer)
     x=Diff_cpu_time(no_timer)
     
  END SUBROUTINE start_timer
  
  
  
  !!  =============================================================================================================================
  !! SUBROUTINE:  stop_timer
  !!
  !>\BRIEF	This subroutine will change the value of timer_state from "running" to "stopped"
  !!
  !! DESCRIPTION:	This subroutine will change the value of timer_state from "running" to "stopped"	 
  !!
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE stop_timer(no_timer)
  IMPLICIT NONE
    INTEGER :: no_timer
    
     IF (timer_state(no_timer)==running) THEN
        CALL suspend_timer(no_timer)
     ELSE IF (timer_state(no_timer)==stopped) THEN
       WRITE(numout,*) 'stop_timer :: timer is already stopped'
     ENDIF

     timer_state(no_timer)=stopped

  END SUBROUTINE stop_timer
  
  
  
  !!  =============================================================================================================================
  !! SUBROUTINE:  resume_timer
  !!
  !>\BRIEF	This subroutine will change the value of timer_state from suspended to running and update diff_cpu_time and diff_real_time
  !!
  !! DESCRIPTION:	This subroutine will change the value of timer_state from suspended to running and 
  !!                update diff_cpu_time and diff_real_time
  !!
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE resume_timer(no_timer)
  IMPLICIT NONE
    INTEGER :: no_timer
    DOUBLE PRECISION :: x
     IF (timer_state(no_timer)/=suspended) THEN
       STOP 'resume_timer :: timer is not suspended'
     ELSE
        timer_state(no_timer)=running
     ENDIF
  
     x=Diff_cpu_time(no_timer)
     x=Diff_real_time(no_timer)  
  
  END SUBROUTINE resume_timer
  
  
  
  !!  =============================================================================================================================
  !! SUBROUTINE:  suspend_timer
  !!
  !>\BRIEF	This subroutine will change the value of timer_state from running to suspended and update cpu_timer
  !!
  !! DESCRIPTION:    This subroutine will change the value of timer_state from running to suspended and update cpu_timer
  !!
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE suspend_timer(no_timer)
  
    IMPLICIT NONE
    INTEGER :: no_timer
    
     IF (timer_state(no_timer)/=running) THEN
       STOP 'suspend_timer :: timer is not running'
     ELSE
        timer_state(no_timer)=suspended
     ENDIF
  
     cpu_timer(no_timer)=cpu_timer(no_timer)+Diff_cpu_time(no_timer)
     real_timer(no_timer)=real_timer(no_timer)+Diff_real_time(no_timer)
  
  END SUBROUTINE suspend_timer
  
  
  !!  =============================================================================================================================
  !! FUNCTION:  diff_real_time
  !!
  !>\BRIEF	Calculate the increment of real time since the last call to diff_real_time
  !!
  !! DESCRIPTION:	 Calculate the increment of real time since the last call to diff_real_time
  !!
  !! \n
  !_ ==============================================================================================================================
  FUNCTION diff_real_time(no_timer)
  IMPLICIT NONE
    INTEGER :: no_timer
    DOUBLE PRECISION :: Diff_real_Time
    integer :: Last_Count,count,count_rate,count_max
    
    Last_Count=Last_real_time(no_timer)
    
    call system_clock(count,count_rate,count_max)
    if (Count>=Last_Count) then
      Diff_real_time=(1.*(Count-last_Count))/count_rate
    else
      Diff_real_time=(1.*(Count-last_Count+Count_max))/count_rate
    endif
    Last_real_time(no_timer)=Count 
    
  END FUNCTION diff_real_time
  
  !!  =============================================================================================================================
  !! FUNCTION:  Diff_Cpu_Time
  !!
  !>\BRIEF     Calculate the increment of cpu time since the last call to diff_real_time	
  !!
  !! DESCRIPTION:	 Calculate the increment of cpu time since the last call to diff_real_time	
  !!
  !! \n
  !_ ==============================================================================================================================
  function Diff_Cpu_Time(no_timer)
  implicit none
    INTEGER :: no_timer
    DOUBLE PRECISION :: Diff_Cpu_Time
    DOUBLE PRECISION :: Last_Count,Count
    
    Last_Count=Last_cpu_time(no_timer)
    
    call cpu_time(Count)
    Diff_Cpu_Time=Count-Last_Count
    Last_cpu_time(no_timer)=Count 
    
  end function Diff_Cpu_Time
  
  !!  =============================================================================================================================
  !! FUNCTION:  Get_cpu_time
  !!
  !>\BRIEF	Return the value of cpu_timer 
  !!
  !! DESCRIPTION:	This subroutine will return the value of cpu_timer 
  !!
  !! \n
  !_ ==============================================================================================================================
  FUNCTION Get_cpu_time(no_timer)
  IMPLICIT NONE
  INTEGER :: no_timer
  DOUBLE PRECISION :: Get_cpu_time
  
    IF (timer_state(no_timer)==running) THEN
      CALL suspend_timer(no_timer)
      Get_cpu_time=cpu_timer(no_timer)
      CALL resume_timer(no_timer)
    ELSE
      Get_cpu_time=cpu_timer(no_timer)
    ENDIF
    
  END FUNCTION Get_cpu_time
  
  !!  =============================================================================================================================
  !! FUNCTION:  Get_real_time
  !!
  !>\BRIEF		Return the value of real_timer
  !!
  !! DESCRIPTION: This subroutine will return the value of real_timer
  !!
  !! \n
  !_ ==============================================================================================================================
  FUNCTION Get_real_time(no_timer)
  IMPLICIT NONE
  INTEGER :: no_timer
  DOUBLE PRECISION :: Get_real_time
  
    IF (timer_state(no_timer)==running) THEN
      CALL suspend_timer(no_timer)
      Get_real_time=real_timer(no_timer)
      CALL resume_timer(no_timer)
    ELSE
      Get_real_time=real_timer(no_timer)
    ENDIF
  
  END FUNCTION Get_real_time
  
END MODULE Timer
  
