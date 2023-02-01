! ==============================================================================================================================
! MODULE   : ioipls_para
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF          Overlap of IOIPSL functions for specific parallel use in ORCHIDEE.
!!
!!\n DESCRIPTION: This module contains interfaces for some IOIPSL subroutines adapted to be used in parallel mode by ORCHIDEE. 
!!
!!                 Following interfaces are available :
!!                  - getin_p : Read a variable from run.def file. The master process will call getin in IOIPSL. 
!!                              The same result will be known by all processes after the call. 
!!                              The variable can be an integer, real, logical or character string. It can be a scalar or 
!!                              have 1 or 2 dimensions except for character string which can only be scalar or have 1 dimension. 
!!                  - restget_p :   Read a variable from restart file. The master process will call the subroutine restget in IOIPSL. 
!!                                  The variable will be distributed on the local domain for each process. 
!!                                  The variable must be a real and can have 1, 2 or 3 dimensions. It can not be a scalar.
!!                  - restput_p :   Write a variable to restart file. The master process will call the subroutine restput in IOIPSL. 
!!                                  The input variable must be given on the local domain for each process.
!!                                  The variable must be a real and can have 1, 2 or 3 dimensions. It can not be a scalar.
!!                  - histwrite_p : Write a variable to history file. The master process will call the subroutine histwrite in IOIPSL. 
!!                                  The input variable must be given on the local domain for each process. 
!!                                  The variable must be a real and can have 1, 2 or 3 dimensions. It can not be a scalar.
!!
!!                 Note that these subroutines must be called by all MPI processes and all OMP thredds because they contain 
!!                 all a MPI blocker function. 
!!                    
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

MODULE ioipsl_para
  USE ioipsl
  USE mod_orchidee_para_var
  USE mod_orchidee_transfert_para
  USE constantes_var
!-
  IMPLICIT NONE

  INTEGER, SAVE :: orch_domain_id 
!-
   INTEGER :: orch_ipslout=6, orch_ilv_cur=0, orch_ilv_max=0
!$OMP THREADPRIVATE( orch_ipslout, orch_ilv_cur, orch_ilv_max )

!-
!-
#include "src_parallel.h"
!-
  !! ==============================================================================================================================
  !! INTERFACE   : getin_p
  !!
  !>\BRIEF          interface to parallelize the call to getin in IOIPSL
  !!
  !! DESCRIPTION  :  get a variable from a text input file. Need to be call by all process
  !!
  !! \n
  !_ ================================================================================================================================
  INTERFACE getin_p
    MODULE PROCEDURE getin_p_c,getin_p_c1,   &
         getin_p_i,getin_p_i1,getin_p_i2,&
         getin_p_r,getin_p_r1,getin_p_r2,&
         getin_p_l,getin_p_l1,getin_p_l2
  END INTERFACE
!-
  !! ==============================================================================================================================
  !! INTERFACE   : restput_p
  !!
  !>\BRIEF         interface to parallelize the call to restput in IOIPSL
  !!
  !! DESCRIPTION  : allows to re-index data onto the original grid of the restart file. Need to be call by all process
  !!
  !! \n
  !_ ================================================================================================================================
  INTERFACE restput_p
     MODULE PROCEDURE &
          restput_p_r3d, restput_p_r2d, restput_p_r1d, &
          restput_p_opp_r5d, restput_p_opp_r4d, restput_p_opp_r3d, &
          restput_p_opp_r2d, restput_p_opp_r1d, restput_p_nogrid_r1d, &
          restput_p_nogrid_i_scal, restput_p_nogrid_r_scal, &
          restput_p_opp_i1d, restput_p_opp_i2d, restput_p_opp_i3d, &
          restput_p_opp_i4d, restput_p_opp_i5d
  END INTERFACE
!-
  !! ==============================================================================================================================
  !! INTERFACE   : restget_p
  !!
  !>\BRIEF    interface to parallelize the call to restget in IOIPSL     
  !!
  !! DESCRIPTION  : Transform the data from the restart file onto the model grid. 
  !!
  !! \n
  !_ ================================================================================================================================
 INTERFACE restget_p
     MODULE PROCEDURE &
          restget_p_r3d, restget_p_r2d, restget_p_r1d, &
          restget_p_opp_r5d, restget_p_opp_r4d, restget_p_opp_r3d, &
          restget_p_opp_r2d, restget_p_opp_r1d, restget_p_nogrid_r1d, &
          restget_p_nogrid_r_scal, restget_p_nogrid_i_scal, &
          restget_p_opp_i1d, restget_p_opp_i2d, restget_p_opp_i3d, &
          restget_p_opp_i4d, restget_p_opp_i5d
  END INTERFACE

  !! ==============================================================================================================================
  !! INTERFACE   : histwrite_p
  !!
  !>\BRIEF         interface to parallelize the call to histwrite in IOIPSL
  !!
  !! DESCRIPTION  : give the data to the IOIPSL system (if we don't use XIOS). Need to be call by all process
  !!
  !! \n
  !_ ================================================================================================================================

  INTERFACE histwrite_p
     MODULE PROCEDURE &
     histwrite_r1d_p,histwrite_r2d_p,histwrite_r3d_p     
  END INTERFACE

  !! ==============================================================================================================================
  !! INTERFACE   : ipslerr_p
  !!
  !>\BRIEF         information subroutine for developpers
  !!
  !! DESCRIPTION  : provide an information subroutine in 3 modes: debug, warn or error(stops)
  !!
  !! \n
  !_ ================================================================================================================================

  INTERFACE ipslerr_p
     MODULE PROCEDURE &
     ipslerr_p_str, ipslerr_p_int
  END INTERFACE
CONTAINS


  !!  =============================================================================================================================
  !! SUBROUTINE:  Init_ioipsl_para 
  !!
  !>\BRIEF	 call to IOIPSL routine : flio_dom_set 
  !!
  !! DESCRIPTION:	 will sets up the domain activity of IOIPSL. Need to be call by all process
  !!
  !! \n
  !_ ==============================================================================================================================

  SUBROUTINE Init_ioipsl_para

    IMPLICIT NONE
    
    INTEGER,DIMENSION(2) :: ddid
    INTEGER,DIMENSION(2) :: dsg
    INTEGER,DIMENSION(2) :: dsl
    INTEGER,DIMENSION(2) :: dpf
    INTEGER,DIMENSION(2) :: dpl
    INTEGER,DIMENSION(2) :: dhs
    INTEGER,DIMENSION(2) :: dhe 

    IF (is_omp_root) THEN
      ddid=(/ 1,2 /)
      dsg=(/ iim_g, jjm_g /)
      dsl=(/ iim_g, jj_nb /)
      dpf=(/ 1,jj_begin /)
      dpl=(/ iim_g, jj_end /)
      dhs=(/ ii_begin-1,0 /)
      if (mpi_rank==mpi_size-1) then
        dhe=(/0,0/)
      else
         dhe=(/ iim_g-ii_end,0 /)  
      endif
    
      call flio_dom_set(mpi_size,mpi_rank,ddid,dsg,dsl,dpf,dpl,dhs,dhe, &
                        'APPLE',orch_domain_id)
     ENDIF
     
  END SUBROUTINE Init_ioipsl_para

  !!  =============================================================================================================================
  !! SUBROUTINE:   ioconf_setatt_p
  !!
  !>\BRIEF	parallelisation of the call to IOIPSL routine ioconf_setatt 
  !!
  !! DESCRIPTION:    NONE
  !!
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE ioconf_setatt_p (attname,attvalue)
    !---------------------------------------------------------------------
    IMPLICIT NONE
    !-
    CHARACTER(LEN=*), INTENT(in) :: attname,attvalue
    !---------------------------------------------------------------------

    IF (is_root_prc) THEN 
       CALL ioconf_setatt(attname,attvalue)
    ENDIF

  END SUBROUTINE ioconf_setatt_p

  !!  =============================================================================================================================
  !! SUBROUTINE:   ipslnlf_p
  !!
  !>\BRIEF	 parallelisation of the call to IOIPSL routine ipslnlf
  !!
  !! DESCRIPTION:  The "ipslnlf" routine allows to know and modify the current logical number for the messages.
  !!
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE ipslnlf_p (new_number,old_number)
    !!--------------------------------------------------------------------
    !! The "ipslnlf" routine allows to know and modify
    !! the current logical number for the messages.
    !!
    !! SUBROUTINE ipslnlf (new_number,old_number)
    !!
    !! Optional INPUT argument
    !!
    !! (I) new_number : new logical number of the file
    !!
    !! Optional OUTPUT argument
    !!
    !! (I) old_number : current logical number of the file
    !!--------------------------------------------------------------------
    IMPLICIT NONE
    !-
    INTEGER,OPTIONAL,INTENT(IN)  :: new_number
    INTEGER,OPTIONAL,INTENT(OUT) :: old_number
    !---------------------------------------------------------------------
    IF (PRESENT(old_number)) THEN
#ifndef CPP_OMP
       CALL ipslnlf(old_number=orch_ipslout)
#endif
       old_number = orch_ipslout
    ENDIF
    IF (PRESENT(new_number)) THEN
       orch_ipslout = new_number
#ifndef CPP_OMP
       CALL ipslnlf(new_number=orch_ipslout)
#endif
    ENDIF

  END SUBROUTINE ipslnlf_p

  !!  =============================================================================================================================
  !! SUBROUTINE:   ipslerr_p_str
  !!
  !>\BRIEF         allows to handle the messages to the user.	 
  !!
  !! DESCRIPTION: NONE
  !!
  !! \n
  !_ ==============================================================================================================================
  !===
  SUBROUTINE ipslerr_p_str (plev,pcname,pstr1,pstr2,pstr3)
    !---------------------------------------------------------------------
    !! The "ipslerr_p_str" routine
    !! allows to handle the messages to the user.
    !!
    !! parallel version of IOIPSL ipslerr
    !!
    !! INPUT
    !!
    !! plev   : Category of message to be reported to the user
    !!          1 = Note to the user
    !!          2 = Warning to the user
    !!          3 = Fatal error
    !! pcname : Name of subroutine which has called ipslerr
    !! pstr1   
    !! pstr2  : Strings containing the explanations to the user
    !! pstr3
    !---------------------------------------------------------------------

#ifdef CPP_PARA
    USE mpi
#endif

    IMPLICIT NONE

    INTEGER :: plev
    CHARACTER(LEN=*) :: pcname,pstr1,pstr2,pstr3

    CHARACTER(LEN=30),DIMENSION(3) :: pemsg = &
         &  (/ "NOTE TO THE USER FROM ROUTINE ", &
         &     "WARNING FROM ROUTINE          ", &
         &     "FATAL ERROR FROM ROUTINE      " /)
    INTEGER :: ierr
    !---------------------------------------------------------------------
    IF ( (plev >= 1).AND.(plev <= 3) ) THEN
       orch_ilv_cur = plev
       orch_ilv_max = MAX(orch_ilv_max,plev)
       WRITE(orch_ipslout,'(/,A," ",A)') TRIM(pemsg(plev)),TRIM(pcname)
       WRITE(orch_ipslout,'(3(" --> ",A,/))') TRIM(pstr1),TRIM(pstr2),TRIM(pstr3)
    ENDIF
    IF (plev == 3) THEN
       WRITE(orch_ipslout,'("Fatal error from ORCHIDEE. STOP in ipslerr_p with code")')
       ! Force to pring text output using FLUSH only if cpp flag CPP_FLUSH is set in arch-XXX.fcm
#ifdef CPP_FLUSH
       CALL FLUSH(orch_ipslout)
#endif

#ifdef CPP_PARA
       CALL MPI_ABORT(MPI_COMM_ORCH, 1, ierr)
#endif     
       STOP 1
    ENDIF
    !---------------------
  END SUBROUTINE ipslerr_p_str

  !!  =============================================================================================================================
  !! SUBROUTINE:   ipslerr_p_int
  !!
  !>\BRIEF         ipslerr_p_str wrapper to allow a int argument in the last position
  !!
  !! DESCRIPTION: NONE
  !!
  !! \n
  !_ ==============================================================================================================================
  !===
  SUBROUTINE ipslerr_p_int (plev,pcname,pstr1,pstr2,pint3)
    IMPLICIT NONE

    INTEGER, INTENT(in) :: plev
    CHARACTER(LEN=*), INTENT(in) :: pcname,pstr1,pstr2
    INTEGER, INTENT(in) :: pint3

    CHARACTER(LEN=30) :: tmp_str

    WRITE(tmp_str, *) pint3
    CALL ipslerr_p_str(plev, pcname, pstr1, pstr2, TRIM(tmp_str))

  END SUBROUTINE ipslerr_p_int

  !!  =============================================================================================================================
  !! SUBROUTINE:  getin_p_c 
  !!
  !>\BRIEF      get a character variable in text input file 	 
  !!
  !! DESCRIPTION: Need to be call by all process 	 
  !!
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE getin_p_c(VarIn,VarOut)
    IMPLICIT NONE    
    CHARACTER(LEN=*),INTENT(IN) :: VarIn
    CHARACTER(LEN=*),INTENT(INOUT) :: VarOut    

    IF (is_root_prc) CALL getin(VarIn,VarOut)
    CALL bcast(VarOut)
  END SUBROUTINE getin_p_c  

  !!  =============================================================================================================================
  !! SUBROUTINE:  getin_p_c1 
  !!
  !>\BRIEF	  get a character 1D array in text input file 
  !!
  !! DESCRIPTION: Need to be call by all process 
  !!
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE getin_p_c1(VarIn,VarOut)
    IMPLICIT NONE    
    CHARACTER(LEN=*),INTENT(IN) :: VarIn
    CHARACTER(LEN=*),INTENT(INOUT) :: VarOut(:)    

    IF (is_root_prc) CALL getin(VarIn,VarOut)
    CALL bcast(VarOut)
  END SUBROUTINE getin_p_c1 

  !!  =============================================================================================================================
  !! SUBROUTINE: getin_p_i  
  !!
  !>\BRIEF	  get an integer variable in text input file 	 
  !!
  !! DESCRIPTION: Need to be call by all process 
  !!
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE getin_p_i(VarIn,VarOut)
    IMPLICIT NONE    
    CHARACTER(LEN=*),INTENT(IN) :: VarIn
    INTEGER,INTENT(INOUT) :: VarOut    

    IF (is_root_prc) CALL getin(VarIn,VarOut)
    CALL bcast(VarOut)
  END SUBROUTINE getin_p_i

  !!  =============================================================================================================================
  !! SUBROUTINE:  getin_p_i1 
  !!
  !>\BRIEF	 get an integer 1D array in text input file 
  !!
  !! DESCRIPTION:  Need to be call by all process 
  !!
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE getin_p_i1(VarIn,VarOut)
    IMPLICIT NONE    
    CHARACTER(LEN=*),INTENT(IN) :: VarIn
    INTEGER,INTENT(INOUT) :: VarOut(:)

    IF (is_root_prc) CALL getin(VarIn,VarOut)
    CALL bcast(VarOut)
  END SUBROUTINE getin_p_i1

  !!  =============================================================================================================================
  !! SUBROUTINE:  getin_p_i2 
  !!
  !>\BRIEF     get an integer 2D array in text input file 	 
  !!
  !! DESCRIPTION: Need to be call by all process 	 
  !!
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE getin_p_i2(VarIn,VarOut)
    IMPLICIT NONE    
    CHARACTER(LEN=*),INTENT(IN) :: VarIn
    INTEGER,INTENT(INOUT) :: VarOut(:,:)

    IF (is_root_prc) CALL getin(VarIn,VarOut)
    CALL bcast(VarOut)
  END SUBROUTINE getin_p_i2

  !!  =============================================================================================================================
  !! SUBROUTINE:   getin_p_r
  !!
  !>\BRIEF        get a float variable in text input file 	 	 
  !!
  !! DESCRIPTION: Need to be call by all process 
  !!
  !! \n
  !_ ==============================================================================================================================
   SUBROUTINE getin_p_r(VarIn,VarOut)
    IMPLICIT NONE    
    CHARACTER(LEN=*),INTENT(IN) :: VarIn
    REAL,INTENT(INOUT) :: VarOut

    IF (is_root_prc) CALL getin(VarIn,VarOut)
    CALL bcast(VarOut)
  END SUBROUTINE getin_p_r

  !!  =============================================================================================================================
  !! SUBROUTINE:  getin_p_r1 
  !!
  !>\BRIEF	 get a float 1D array in text input file  
  !!
  !! DESCRIPTION: Need to be call by all process 
  !!
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE getin_p_r1(VarIn,VarOut)
    IMPLICIT NONE    
    CHARACTER(LEN=*),INTENT(IN) :: VarIn
    REAL,INTENT(INOUT) :: VarOut(:)

    IF (is_root_prc) CALL getin(VarIn,VarOut)
    CALL bcast(VarOut)
  END SUBROUTINE getin_p_r1

  !!  =============================================================================================================================
  !! SUBROUTINE:  getin_p_r2 
  !!
  !>\BRIEF	 get a float 2D array in text input file  
  !!
  !! DESCRIPTION: Need to be call by all process  
  !!
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE getin_p_r2(VarIn,VarOut)
    IMPLICIT NONE    
    CHARACTER(LEN=*),INTENT(IN) :: VarIn
    REAL,INTENT(INOUT) :: VarOut(:,:)

    IF (is_root_prc) CALL getin(VarIn,VarOut)
    CALL bcast(VarOut)
  END SUBROUTINE getin_p_r2


  !!  =============================================================================================================================
  !! SUBROUTINE:  getin_p_l 
  !!
  !>\BRIEF	  get a logical variable in text input file 
  !!
  !! DESCRIPTION: Need to be call by all process 
  !!
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE getin_p_l(VarIn,VarOut)
    IMPLICIT NONE    
    CHARACTER(LEN=*),INTENT(IN) :: VarIn
    LOGICAL,INTENT(INOUT) :: VarOut

    IF (is_root_prc) CALL getin(VarIn,VarOut)
    CALL bcast(VarOut)
  END SUBROUTINE getin_p_l

  !!  =============================================================================================================================
  !! SUBROUTINE:   getin_p_l1
  !!
  !>\BRIEF      get a logical 1D array in text input file 	 
  !!
  !! DESCRIPTION: Need to be call by all process
  !!
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE getin_p_l1(VarIn,VarOut)
    IMPLICIT NONE    
    CHARACTER(LEN=*),INTENT(IN) :: VarIn
    LOGICAL,INTENT(INOUT) :: VarOut(:)

    IF (is_root_prc) CALL getin(VarIn,VarOut)
    CALL bcast(VarOut)
  END SUBROUTINE getin_p_l1

  !!  =============================================================================================================================
  !! SUBROUTINE:  getin_p_l2 
  !!
  !>\BRIEF	 get a logical 2D array in text input file 
  !!
  !! DESCRIPTION: Need to be call by all process
  !!
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE getin_p_l2(VarIn,VarOut)
    IMPLICIT NONE    
    CHARACTER(LEN=*),INTENT(IN) :: VarIn
    LOGICAL,INTENT(INOUT) :: VarOut(:,:)

    IF (is_root_prc) CALL getin(VarIn,VarOut)
    CALL bcast(VarOut)
  END SUBROUTINE getin_p_l2
!-

  !!  =============================================================================================================================
  !! SUBROUTINE:  restget_p_opp_r1d 
  !!
  !>\BRIEF	 Transform the data (real 1D) from the restart file onto the model grid with the operation MY_OPERATOR
  !!
  !! DESCRIPTION: do not use this function with non grid variable
  !!
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE restget_p_opp_r1d &
  (fid, vname_q, iim, jjm, llm, itau, def_beha, &
   var, MY_OPERATOR, nbindex, ijndex)
! DO NOT USE THIS FUNCTION WITH NON GRID VARIABLE !
    IMPLICIT NONE
!-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: iim, jjm, llm, itau
    LOGICAL def_beha
    REAL :: var(:)
    CHARACTER(LEN=*) :: MY_OPERATOR
    INTEGER :: nbindex, ijndex(nbindex)
    !-----------------------------
    REAL, ALLOCATABLE, DIMENSION(:) :: temp_g

    IF (is_root_prc) THEN 
       ALLOCATE( temp_g(iim*jjm*llm) )
    ELSE
       ALLOCATE( temp_g(1) )
    ENDIF
        
    IF (is_root_prc) THEN 
       CALL restget &
            (fid, vname_q, iim, jjm, llm, itau, def_beha, &
            temp_g, MY_OPERATOR, nbindex, ijndex)
    ENDIF
    CALL scatter(temp_g,var)
    DEALLOCATE(temp_g)
  END SUBROUTINE restget_p_opp_r1d

  !!  =============================================================================================================================
  !! SUBROUTINE:   restget_p_opp_r2d
  !!
  !>\BRIEF	Transform the data (real 2D) from the restart file onto the model grid with the operation MY_OPERATOR 
  !!
  !! DESCRIPTION: do not use this function with non grid variable.  Need to be call by all process
  !!
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE restget_p_opp_r2d &
  (fid, vname_q, iim, jjm, llm, itau, def_beha, &
   var, MY_OPERATOR, nbindex, ijndex)
    IMPLICIT NONE
    !-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: iim, jjm, llm, itau
    LOGICAL def_beha
    REAL :: var(:,:)
    CHARACTER(LEN=*) :: MY_OPERATOR
    INTEGER :: nbindex, ijndex(nbindex)
    !-----------------------------
    REAL, ALLOCATABLE, DIMENSION(:,:) :: temp_g

    IF (is_root_prc) THEN 
       ALLOCATE( temp_g(iim,jjm) )
    ELSE
      ALLOCATE( temp_g(1,1) )
    ENDIF

    IF (is_root_prc) THEN 
       CALL restget &
            (fid, vname_q, iim, jjm, llm, itau, def_beha, &
            temp_g, MY_OPERATOR, nbindex, ijndex)
    ENDIF
    CALL scatter(temp_g,var)
    DEALLOCATE(temp_g)
  END SUBROUTINE restget_p_opp_r2d

  !!  =============================================================================================================================
  !! SUBROUTINE:   restget_p_opp_r2d
  !!
  !>\BRIEF	Transform the data (real 2D) from the restart file onto the model grid with the operation MY_OPERATOR 
  !!
  !! DESCRIPTION: do not use this function with non grid variable.  Need to be call by all process
  !!
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE restget_p_opp_r3d &
  (fid, vname_q, iim, jjm, llm, itau, def_beha, &
   var, MY_OPERATOR, nbindex, ijndex)
    IMPLICIT NONE
    !-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: iim, jjm, llm, itau
    LOGICAL def_beha
    REAL :: var(:,:,:)
    CHARACTER(LEN=*) :: MY_OPERATOR
    INTEGER :: nbindex, ijndex(nbindex)
    !-----------------------------
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: temp_g

    IF (is_root_prc) THEN 
       ALLOCATE( temp_g(iim,jjm,llm) )
    ELSE
      ALLOCATE( temp_g(1,1,1) )
    ENDIF

    IF (is_root_prc) THEN 
       CALL restget &
            (fid, vname_q, iim, jjm, llm, itau, def_beha, &
            temp_g, MY_OPERATOR, nbindex, ijndex)
    ENDIF
    CALL scatter(temp_g,var)
    DEALLOCATE(temp_g)

END SUBROUTINE restget_p_opp_r3d

  !!  =============================================================================================================================
  !! SUBROUTINE:   restget_p_opp_r4d
  !!
  !>\BRIEF	Transform the data (real 4D) from the restart file onto the model grid with the operation MY_OPERATOR 
  !!
  !! DESCRIPTION: do not use this function with non grid variable.  Need to be call by all process
  !!
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE restget_p_opp_r4d &
  (fid, vname_q, iim, jjm, llm, zzm, itau, def_beha, &
   var, MY_OPERATOR, nbindex, ijndex)
    IMPLICIT NONE
    !-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: iim, jjm, llm, zzm, itau
    LOGICAL def_beha
    REAL :: var(:,:,:,:)
    CHARACTER(LEN=*) :: MY_OPERATOR
    INTEGER :: nbindex, ijndex(nbindex)
    !-----------------------------
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: temp_g

   IF (is_root_prc) THEN 
       ALLOCATE( temp_g(iim,jjm,llm,zzm) )
    ELSE
      ALLOCATE( temp_g(1,1,1,1) )
    ENDIF

    IF (is_root_prc) THEN 
       CALL restget &
            (fid, vname_q, iim, jjm, llm, zzm, itau, def_beha,  &
            temp_g, MY_OPERATOR, nbindex, ijndex)
   ENDIF
    CALL scatter(temp_g,var)
    DEALLOCATE(temp_g)

END SUBROUTINE restget_p_opp_r4d

  !!  =============================================================================================================================
  !! SUBROUTINE:   restget_p_opp_r5d
  !!
  !>\BRIEF	Transform the data (real 5D) from the restart file onto the model grid with the operation MY_OPERATOR 
  !!
  !! DESCRIPTION: do not use this function with non grid variable.  Need to be call by all process
  !!
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE restget_p_opp_r5d &
  (fid, vname_q, iim, jjm, llm, zzm, wwm, itau, def_beha, &
   var, MY_OPERATOR, nbindex, ijndex)
    IMPLICIT NONE
    !-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: iim, jjm, llm, zzm, wwm, itau
    LOGICAL def_beha
    REAL :: var(:,:,:,:,:)
    CHARACTER(LEN=*) :: MY_OPERATOR
    INTEGER :: nbindex, ijndex(nbindex)
    !-----------------------------
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: temp_g

    IF (is_root_prc) THEN 
       ALLOCATE( temp_g(iim,jjm,llm,zzm,wwm) )
    ELSE
      ALLOCATE( temp_g(1,1,1,1,1) )
    ENDIF

    IF (is_root_prc) THEN 
       CALL restget &
            (fid, vname_q, iim, jjm, llm, zzm, wwm, itau, def_beha,  &
            temp_g, MY_OPERATOR, nbindex, ijndex)
    ENDIF
    CALL scatter(temp_g,var)
    DEALLOCATE(temp_g)

END SUBROUTINE restget_p_opp_r5d

  !!  =============================================================================================================================
  !! SUBROUTINE:   restget_p_opp_i1d
  !!
  !>\BRIEF	Transform the data (integer 1D) from the restart file onto the model grid with the operation MY_OPERATOR 
  !!
  !! DESCRIPTION: do not use this function with non grid variable.  Need to be call by all process
  !!
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE restget_p_opp_i1d &
  (fid, vname_q, iim, jjm, llm, itau, def_beha, &
   var, MY_OPERATOR, nbindex, ijndex)
    IMPLICIT NONE
    !-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: iim, jjm, llm, itau
    LOGICAL def_beha
    INTEGER, INTENT(out) :: var(:)
    CHARACTER(LEN=*) :: MY_OPERATOR
    INTEGER :: nbindex, ijndex(nbindex)
    !-----------------------------
    REAL, ALLOCATABLE, DIMENSION(:) :: temp_g
    INTEGER :: ier

    ALLOCATE( temp_g(SIZE(var, DIM=1)), stat=ier ) 
    IF (ier /= 0) CALL ipslerr_p(3, 'restget_p_opp_i1d', 'Memory allocation error', vname_q, '')

    CALL restget_p &
         (fid, vname_q, iim, jjm, llm, itau, def_beha, &
         temp_g, MY_OPERATOR, nbindex, ijndex)
    var = INT(temp_g, i_std)
    
    DEALLOCATE(temp_g)
  END SUBROUTINE restget_p_opp_i1d

  !!  =============================================================================================================================
  !! SUBROUTINE:   restget_p_opp_i2d
  !!
  !>\BRIEF	Transform the data (integer 2D) from the restart file onto the model grid with the operation MY_OPERATOR 
  !!
  !! DESCRIPTION: do not use this function with non grid variable.  Need to be call by all process
  !!
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE restget_p_opp_i2d &
  (fid, vname_q, iim, jjm, llm, itau, def_beha, &
   var, MY_OPERATOR, nbindex, ijndex)
    IMPLICIT NONE
    !-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: iim, jjm, llm, itau
    LOGICAL def_beha
    INTEGER, INTENT(out) :: var(:,:)
    CHARACTER(LEN=*) :: MY_OPERATOR
    INTEGER :: nbindex, ijndex(nbindex)
    !-----------------------------
    REAL, ALLOCATABLE, DIMENSION(:,:) :: temp_g
    INTEGER :: ier

    ALLOCATE( temp_g(SIZE(var, DIM=1), SIZE(var, DIM=2)), stat=ier ) 
    IF (ier /= 0) CALL ipslerr_p(3, 'restget_p_opp_i2d', 'Memory allocation error', vname_q, '')

    CALL restget_p &
         (fid, vname_q, iim, jjm, llm, itau, def_beha, &
         temp_g, MY_OPERATOR, nbindex, ijndex)
    var = INT(temp_g, i_std)
    
    DEALLOCATE(temp_g)
  END SUBROUTINE restget_p_opp_i2d

  !!  =============================================================================================================================
  !! SUBROUTINE:   restget_p_opp_i3d
  !!
  !>\BRIEF	Transform the data (integer 3D) from the restart file onto the model grid with the operation MY_OPERATOR 
  !!
  !! DESCRIPTION: do not use this function with non grid variable.  Need to be call by all process
  !!
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE restget_p_opp_i3d &
  (fid, vname_q, iim, jjm, llm, itau, def_beha, &
   var, MY_OPERATOR, nbindex, ijndex)
    IMPLICIT NONE
    !-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: iim, jjm, llm, itau
    LOGICAL def_beha
    INTEGER, INTENT(out) :: var(:,:,:)
    CHARACTER(LEN=*) :: MY_OPERATOR
    INTEGER :: nbindex, ijndex(nbindex)
    !-----------------------------
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: temp_g
    INTEGER :: ier

    ALLOCATE( temp_g(SIZE(var, DIM=1), SIZE(var, DIM=2), SIZE(var, DIM=3)), stat=ier ) 
    IF (ier /= 0) CALL ipslerr_p(3, 'restget_p_opp_i3d', 'Memory allocation error', vname_q, '')

    CALL restget_p &
         (fid, vname_q, iim, jjm, llm, itau, def_beha, &
         temp_g, MY_OPERATOR, nbindex, ijndex)
    var = INT(temp_g, i_std)
    
    DEALLOCATE(temp_g)
  END SUBROUTINE restget_p_opp_i3d

  !!  =============================================================================================================================
  !! SUBROUTINE:   restget_p_opp_i4d
  !!
  !>\BRIEF	Transform the data (integer 4D) from the restart file onto the model grid with the operation MY_OPERATOR 
  !!
  !! DESCRIPTION: do not use this function with non grid variable.  Need to be call by all process
  !!
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE restget_p_opp_i4d &
  (fid, vname_q, iim, jjm, llm, mmm, itau, def_beha, &
   var, MY_OPERATOR, nbindex, ijndex)
    IMPLICIT NONE
    !-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: iim, jjm, llm, mmm, itau
    LOGICAL def_beha
    INTEGER, INTENT(out) :: var(:,:,:,:)
    CHARACTER(LEN=*) :: MY_OPERATOR
    INTEGER :: nbindex, ijndex(nbindex)
    !-----------------------------
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: temp_g
    INTEGER :: ier

    ALLOCATE( temp_g(SIZE(var, DIM=1), SIZE(var, DIM=2), SIZE(var, DIM=3), SIZE(var, DIM=4)), stat=ier ) 
    IF (ier /= 0) CALL ipslerr_p(3, 'restget_p_opp_i4d', 'Memory allocation error', vname_q, '')

    CALL restget_p &
         (fid, vname_q, iim, jjm, llm, mmm, itau, def_beha, &
         temp_g, MY_OPERATOR, nbindex, ijndex)
    var = INT(temp_g, i_std)
    
    DEALLOCATE(temp_g)
  END SUBROUTINE restget_p_opp_i4d

  !!  =============================================================================================================================
  !! SUBROUTINE:   restget_p_opp_i2d
  !!
  !>\BRIEF	Transform the data (integer 5D) from the restart file onto the model grid with the operation MY_OPERATOR 
  !!
  !! DESCRIPTION: do not use this function with non grid variable.  Need to be call by all process
  !!
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE restget_p_opp_i5d &
  (fid, vname_q, iim, jjm, llm, mmm, wwm, itau, def_beha, &
   var, MY_OPERATOR, nbindex, ijndex)
    IMPLICIT NONE
    !-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: iim, jjm, llm, mmm, wwm, itau
    LOGICAL def_beha
    INTEGER, INTENT(out) :: var(:,:,:,:,:)
    CHARACTER(LEN=*) :: MY_OPERATOR
    INTEGER :: nbindex, ijndex(nbindex)
    !-----------------------------
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: temp_g
    INTEGER :: ier

    ALLOCATE( temp_g(SIZE(var, DIM=1), SIZE(var, DIM=2), SIZE(var, DIM=3), SIZE(var, DIM=4), SIZE(var, DIM=5)), stat=ier ) 
    IF (ier /= 0) CALL ipslerr_p(3, 'restget_p_opp_i5d', 'Memory allocation error', vname_q, '')

    CALL restget_p &
         (fid, vname_q, iim, jjm, llm, mmm, wwm, itau, def_beha, &
         temp_g, MY_OPERATOR, nbindex, ijndex)
    var = INT(temp_g, i_std)
    
    DEALLOCATE(temp_g)
  END SUBROUTINE restget_p_opp_i5d

!!  =============================================================================================================================
!! SUBROUTINE:   restget_p_r1d
!!
!>\BRIEF	Transform the data (real 1D) from the restart file onto the model grid 	 
!!
!! DESCRIPTION:	 do not use this function with non grid variable.  Need to be call by all process
!! \n
!_ ==============================================================================================================================
  SUBROUTINE restget_p_r1d &
  (fid,vname_q,iim,jjm,llm,itau,def_beha,var)
! DO NOT USE THIS FUNCTION WITH NON GRID VARIABLE !
    IMPLICIT NONE
!-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: iim, jjm, llm, itau
    LOGICAL :: def_beha
    REAL :: var(:)
    !-------------------------
    REAL, ALLOCATABLE, DIMENSION(:) :: temp_g

    IF (is_root_prc) THEN 
       ALLOCATE( temp_g(iim*jjm*llm) )
    ELSE
       ALLOCATE( temp_g(1) )
    ENDIF

    IF (is_root_prc) THEN 
       CALL restget &
            (fid,vname_q,iim,jjm,llm,itau,def_beha,temp_g)
    ENDIF
    CALL scatter2D_mpi(temp_g,var)
    DEALLOCATE(temp_g)
  END SUBROUTINE restget_p_r1d

!!  =============================================================================================================================
!! SUBROUTINE:   restget_p_r2d
!!
!>\BRIEF	Transform the data (real 2D) from the restart file onto the model grid 	 
!!
!! DESCRIPTION:	 do not use this function with non grid variable.  Need to be call by all process
!! \n
!_ ==============================================================================================================================
  SUBROUTINE restget_p_r2d &
  (fid,vname_q,iim,jjm,llm,itau,def_beha,var)
    IMPLICIT NONE
!-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: iim, jjm, llm, itau
    LOGICAL :: def_beha
    REAL :: var(:,:)
    !-------------------------
    REAL, ALLOCATABLE, DIMENSION(:,:) :: temp_g

    IF (is_root_prc) THEN 
       ALLOCATE( temp_g(iim,jjm) )
    ELSE
       ALLOCATE( temp_g(1,1) )
    ENDIF
    IF (is_root_prc) THEN 
       CALL restget &
            (fid,vname_q,iim,jjm,llm,itau,def_beha,temp_g)
    ENDIF
    CALL scatter2D_mpi(temp_g,var)
    DEALLOCATE(temp_g)
  END SUBROUTINE restget_p_r2d

!!  =============================================================================================================================
!! SUBROUTINE:   restget_p_r3d
!!
!>\BRIEF	Transform the data (real 3D) from the restart file onto the model grid 	 
!!
!! DESCRIPTION:	 do not use this function with non grid variable.  Need to be call by all process
!! \n
!_ ==============================================================================================================================
  SUBROUTINE restget_p_r3d &
  (fid,vname_q,iim,jjm,llm,itau,def_beha,var)
    IMPLICIT NONE
!-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: iim, jjm, llm, itau
    LOGICAL def_beha
    REAL :: var(:,:,:)
    !-------------------------
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: temp_g

    IF (is_root_prc) THEN 
       ALLOCATE( temp_g(iim,jjm,llm) )
    ELSE 
       ALLOCATE( temp_g(1,1,1) )
    ENDIF
    
    IF (is_root_prc) THEN 
       CALL restget &
            (fid,vname_q,iim,jjm,llm,itau,def_beha,temp_g)
    ENDIF
    CALL scatter2D_mpi(temp_g,var)
    DEALLOCATE(temp_g)
  END SUBROUTINE restget_p_r3d

!!  =============================================================================================================================
!! SUBROUTINE:   restget_p_nogrid_r1d
!!
!>\BRIEF	Transform the data (real 1D) from the restart file onto the model grid 	 
!!
!! DESCRIPTION:	 
!! \n
!_ ==============================================================================================================================
  SUBROUTINE restget_p_nogrid_r1d &
  (fid,vname_q,itau,def_beha,def_val,var)
! 
    IMPLICIT NONE
!-
    INTEGER, INTENT(in)             :: fid
    CHARACTER(LEN=*), INTENT(in)    :: vname_q
    INTEGER, INTENT(in)             :: itau
    LOGICAL, INTENT(in)             :: def_beha
    REAL, INTENT(in)                :: def_val
    REAL, DIMENSION(:), INTENT(out) :: var
    !-------------------------
    IF (is_root_prc) THEN
       var = val_exp
       CALL restget (fid, vname_q, 1 ,1  , 1, itau, def_beha, var)
       IF(ALL(var == val_exp)) var = def_val 
    ENDIF
    CALL bcast(var)

  END SUBROUTINE restget_p_nogrid_r1d

!!  =============================================================================================================================
!! SUBROUTINE:   restget_p_nogrid_r_scal
!!
!>\BRIEF	Transform the data (real scalar) from the restart file onto the model grid 	 
!!
!! DESCRIPTION:	 
!! \n
!_ ==============================================================================================================================
  SUBROUTINE restget_p_nogrid_r_scal &
  (fid,vname_q,itau,def_beha,def_val,var)
! 
    IMPLICIT NONE
!-
    INTEGER, INTENT(in)             :: fid
    CHARACTER(LEN=*), INTENT(in)    :: vname_q
    INTEGER, INTENT(in)             :: itau
    LOGICAL, INTENT(in)             :: def_beha
    REAL, INTENT(in)                :: def_val
    REAL, INTENT(out) :: var
    !-------------------------
    REAL, DIMENSION(1) :: tmp

    tmp(1) = var
    IF (is_root_prc) THEN
       var = val_exp
       CALL restget (fid, vname_q, 1 ,1  , 1, itau, def_beha, tmp)
       var = tmp(1)
       IF(var == val_exp) var = def_val 
    ENDIF
    CALL bcast(var)

  END SUBROUTINE restget_p_nogrid_r_scal

!!  =============================================================================================================================
!! SUBROUTINE:   restget_p_nogrid_i_scal
!!
!>\BRIEF	Transform the data (integer scalar) from the restart file onto the model grid 	 
!!
!! DESCRIPTION:	 
!! \n
!_ ==============================================================================================================================
  SUBROUTINE restget_p_nogrid_i_scal &
  (fid,vname_q,itau,def_beha,def_val,varint)
! 
    IMPLICIT NONE
!-
    INTEGER, INTENT(in)             :: fid
    CHARACTER(LEN=*), INTENT(in)    :: vname_q
    INTEGER, INTENT(in)             :: itau
    LOGICAL, INTENT(in)             :: def_beha
    REAL, INTENT(in)                :: def_val
    INTEGER, INTENT(out) :: varint
    !-------------------------
    REAL :: tmp

    CALL restget_p_nogrid_r_scal(fid, vname_q, itau, def_beha, def_val, tmp)
    varint = INT(tmp)
  END SUBROUTINE restget_p_nogrid_i_scal

!!  =============================================================================================================================
!! SUBROUTINE:  restput_p_opp_r1d 
!!
!>\BRIEF       allows to re-index data (real 1D) onto the original grid of the restart file with the operation MY_OPERATOR	 
!!
!! DESCRIPTION:	  Need to be call by all process
!! \n
!_ ==============================================================================================================================
  SUBROUTINE restput_p_opp_r1d &
  (fid, vname_q, iim, jjm, llm, itau, var, MY_OPERATOR, nbindex, ijndex)
    IMPLICIT NONE
!-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: iim, jjm, llm, itau
    REAL :: var(:)
    CHARACTER(LEN=*) :: MY_OPERATOR
    INTEGER :: nbindex, ijndex(nbindex)
    !-----------------------------
    REAL, ALLOCATABLE, DIMENSION(:) :: temp_g

    IF (is_root_prc) THEN 
      ALLOCATE( temp_g(iim*jjm*llm) )
    ELSE
      ALLOCATE ( temp_g(1) )
    ENDIF
    
    CALL gather(var,temp_g)
    IF (is_root_prc) THEN
       CALL restput &
            (fid, vname_q, iim, jjm, llm, itau, temp_g, MY_OPERATOR, nbindex, ijndex)
    ENDIF

    DEALLOCATE( temp_g )
          
  END SUBROUTINE restput_p_opp_r1d

!!  =============================================================================================================================
!! SUBROUTINE:  restput_p_opp_r2d 
!!
!>\BRIEF       allows to re-index data (real 2D) onto the original grid of the restart file with the operation MY_OPERATOR	 
!!
!! DESCRIPTION:	  Need to be call by all process
!! \n
!_ ==============================================================================================================================
  SUBROUTINE restput_p_opp_r2d &
  (fid, vname_q, iim, jjm, llm, itau, var, MY_OPERATOR, nbindex, ijndex)
    IMPLICIT NONE
!-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: iim, jjm, llm, itau
    REAL :: var(:,:)
    CHARACTER(LEN=*) :: MY_OPERATOR
    INTEGER :: nbindex, ijndex(nbindex)
    !-----------------------------
    REAL, ALLOCATABLE, DIMENSION(:,:) :: temp_g

    IF (is_root_prc) THEN
      ALLOCATE( temp_g(iim,jjm) )
    ELSE
      ALLOCATE( temp_g(1,1) )
    ENDIF
          
    CALL gather(var,temp_g)
    IF (is_root_prc) THEN
       CALL restput &
            (fid, vname_q, iim, jjm, llm, itau, temp_g, MY_OPERATOR, nbindex, ijndex)
    ENDIF
    DEALLOCATE( temp_g )
          
  END SUBROUTINE restput_p_opp_r2d

!!  =============================================================================================================================
!! SUBROUTINE:  restput_p_opp_r3d 
!!
!>\BRIEF       allows to re-index data (real 3D) onto the original grid of the restart file with the operation MY_OPERATOR	 
!!
!! DESCRIPTION:	  Need to be call by all process
!! \n
!_ ==============================================================================================================================
  SUBROUTINE restput_p_opp_r3d &
  (fid, vname_q, iim, jjm, llm, itau, var, MY_OPERATOR, nbindex, ijndex)
    IMPLICIT NONE
!-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: iim, jjm, llm, itau
    REAL :: var(:,:,:)
    CHARACTER(LEN=*) :: MY_OPERATOR
    INTEGER :: nbindex, ijndex(nbindex)
    !-----------------------------
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: temp_g

    IF (is_root_prc) THEN
      ALLOCATE( temp_g(iim,jjm,llm) )
    ELSE
      ALLOCATE( temp_g(1,1,1) )
    ENDIF
          
    CALL gather(var,temp_g)
    IF (is_root_prc) THEN
       CALL restput &
            (fid, vname_q, iim, jjm, llm, itau, temp_g, MY_OPERATOR, nbindex, ijndex)
    ENDIF
    DEALLOCATE( temp_g )

          
  END SUBROUTINE restput_p_opp_r3d

!!  =============================================================================================================================
!! SUBROUTINE:  restput_p_opp_r4d 
!!
!>\BRIEF       allows to re-index data (real 4D) onto the original grid of the restart file with the operation MY_OPERATOR	 
!!
!! DESCRIPTION:	  Need to be call by all process
!! \n
!_ ==============================================================================================================================
  SUBROUTINE restput_p_opp_r4d &
  (fid, vname_q, iim, jjm, llm, zzm, itau, var, MY_OPERATOR, nbindex, ijndex)
    IMPLICIT NONE
!-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: iim, jjm, llm, zzm, itau
    REAL :: var(:,:,:,:)
    CHARACTER(LEN=*) :: MY_OPERATOR
    INTEGER :: nbindex, ijndex(nbindex)
    !-----------------------------
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: temp_g

    IF (is_root_prc) THEN
      ALLOCATE( temp_g(iim,jjm,llm,zzm) )
    ELSE
      ALLOCATE( temp_g(1,1,1,1) )
    ENDIF
          
    CALL gather(var,temp_g)
    IF (is_root_prc) THEN
       CALL restput &
            (fid, vname_q, iim, jjm, llm, zzm, itau, &
                 temp_g, MY_OPERATOR, nbindex, ijndex)
    ENDIF
    DEALLOCATE( temp_g )

          
  END SUBROUTINE restput_p_opp_r4d

!!  =============================================================================================================================
!! SUBROUTINE:  restput_p_opp_r5d 
!!
!>\BRIEF       allows to re-index data (real 5D) onto the original grid of the restart file with the operation MY_OPERATOR	 
!!
!! DESCRIPTION:	  Need to be call by all process
!! \n
!_ ==============================================================================================================================
  SUBROUTINE restput_p_opp_r5d &
  (fid, vname_q, iim, jjm, llm, zzm, wwm, itau, var, MY_OPERATOR, nbindex, ijndex)
    IMPLICIT NONE
!-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: iim, jjm, llm, zzm, wwm, itau
    REAL :: var(:,:,:,:,:)
    CHARACTER(LEN=*) :: MY_OPERATOR
    INTEGER :: nbindex, ijndex(nbindex)
    !-----------------------------
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: temp_g

    IF (is_root_prc) THEN
      ALLOCATE( temp_g(iim,jjm,llm,zzm,wwm) )
    ELSE
      ALLOCATE( temp_g(1,1,1,1,1) )
    ENDIF
          
    CALL gather(var,temp_g)
    IF (is_root_prc) THEN
       CALL restput &
            (fid, vname_q, iim, jjm, llm, zzm, wwm, itau, &
                 temp_g, MY_OPERATOR, nbindex, ijndex)
    ENDIF
    DEALLOCATE( temp_g )

          
  END SUBROUTINE restput_p_opp_r5d

!!  =============================================================================================================================
!! SUBROUTINE:  restput_p_opp_i1d 
!!
!>\BRIEF       allows to re-index data (integer 2D) onto the original grid of the restart file with the operation MY_OPERATOR	 
!!
!! DESCRIPTION:	  Need to be call by all process
!! \n
!_ ==============================================================================================================================
  SUBROUTINE restput_p_opp_i1d &
  (fid, vname_q, iim, jjm, llm, itau, var, MY_OPERATOR, nbindex, ijndex)
    IMPLICIT NONE
!-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: iim, jjm, llm, itau
    INTEGER :: var(:)
    CHARACTER(LEN=*) :: MY_OPERATOR
    INTEGER :: nbindex, ijndex(nbindex)
    !-----------------------------
    REAL, ALLOCATABLE, DIMENSION(:) :: temp_g
    INTEGER :: ier

    ALLOCATE( temp_g(SIZE(var, DIM=1)), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3, 'restput_p_opp_i1d', 'Allocation memory error ', vname_q, '')

    temp_g = REAL(var, r_std)          
    CALL restput_p &
         (fid, vname_q, iim, jjm, llm, itau, temp_g, MY_OPERATOR, nbindex, ijndex)

    DEALLOCATE( temp_g )
          
  END SUBROUTINE restput_p_opp_i1d

!!  =============================================================================================================================
!! SUBROUTINE:  restput_p_opp_i2d 
!!
!>\BRIEF       allows to re-index data (integer 2D) onto the original grid of the restart file with the operation MY_OPERATOR	 
!!
!! DESCRIPTION:	  Need to be call by all process
!! \n
!_ ==============================================================================================================================
  SUBROUTINE restput_p_opp_i2d &
  (fid, vname_q, iim, jjm, llm, itau, var, MY_OPERATOR, nbindex, ijndex)
    IMPLICIT NONE
!-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: iim, jjm, llm, itau
    INTEGER :: var(:,:)
    CHARACTER(LEN=*) :: MY_OPERATOR
    INTEGER :: nbindex, ijndex(nbindex)
    !-----------------------------
    REAL, ALLOCATABLE, DIMENSION(:,:) :: temp_g
    INTEGER :: ier

    ALLOCATE( temp_g(SIZE(var, DIM=1), SIZE(var, DIM=2)), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3, 'restput_p_opp_i2d', 'Allocation memory error', vname_q, '')

    temp_g = REAL(var, r_std)          
    CALL restput_p &
         (fid, vname_q, iim, jjm, llm, itau, temp_g, MY_OPERATOR, nbindex, ijndex)

    DEALLOCATE( temp_g )
          
  END SUBROUTINE restput_p_opp_i2d

!!  =============================================================================================================================
!! SUBROUTINE:  restput_p_opp_i3d 
!!
!>\BRIEF       allows to re-index data (integer 3D) onto the original grid of the restart file with the operation MY_OPERATOR	 
!!
!! DESCRIPTION:	  Need to be call by all process
!! \n
!_ ==============================================================================================================================
  SUBROUTINE restput_p_opp_i3d &
  (fid, vname_q, iim, jjm, llm, itau, var, MY_OPERATOR, nbindex, ijndex)
    IMPLICIT NONE
!-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: iim, jjm, llm, itau
    INTEGER :: var(:,:,:)
    CHARACTER(LEN=*) :: MY_OPERATOR
    INTEGER :: nbindex, ijndex(nbindex)
    !-----------------------------
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: temp_g
    INTEGER :: ier

    ALLOCATE( temp_g(SIZE(var, DIM=1), SIZE(var, DIM=2), SIZE(var, DIM=3)), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3, 'restput_p_opp_i2d', 'Allocation memory error', vname_q, '')

    temp_g = REAL(var, r_std)          
    CALL restput_p &
         (fid, vname_q, iim, jjm, llm, itau, temp_g, MY_OPERATOR, nbindex, ijndex)

    DEALLOCATE( temp_g )
          
  END SUBROUTINE restput_p_opp_i3d

!!  =============================================================================================================================
!! SUBROUTINE:  restput_p_opp_i4d 
!!
!>\BRIEF       allows to re-index data (integer 4D) onto the original grid of the restart file with the operation MY_OPERATOR	 
!!
!! DESCRIPTION:	  Need to be call by all process
!! \n
!_ ==============================================================================================================================
  SUBROUTINE restput_p_opp_i4d &
  (fid, vname_q, iim, jjm, llm, mmm, itau, var, MY_OPERATOR, nbindex, ijndex)
    IMPLICIT NONE
!-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: iim, jjm, llm, mmm, itau
    INTEGER :: var(:,:,:,:)
    CHARACTER(LEN=*) :: MY_OPERATOR
    INTEGER :: nbindex, ijndex(nbindex)
    !-----------------------------
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: temp_g
    INTEGER :: ier

    ALLOCATE( temp_g(SIZE(var, DIM=1), SIZE(var, DIM=2), SIZE(var, DIM=3), SIZE(var, DIM=4)), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3, 'restput_p_opp_i4d', 'Allocation memory error', vname_q, '')

    temp_g = REAL(var, r_std)          
    CALL restput_p &
         (fid, vname_q, iim, jjm, llm, mmm, itau, temp_g, MY_OPERATOR, nbindex, ijndex)

    DEALLOCATE( temp_g )
          
  END SUBROUTINE restput_p_opp_i4d

!!  =============================================================================================================================
!! SUBROUTINE:  restput_p_opp_i5d 
!!
!>\BRIEF       allows to re-index data (integer 5D) onto the original grid of the restart file with the operation MY_OPERATOR	 
!!
!! DESCRIPTION:	  Need to be call by all process
!! \n
!_ ==============================================================================================================================
  SUBROUTINE restput_p_opp_i5d &
  (fid, vname_q, iim, jjm, llm, mmm, zzm, itau, var, MY_OPERATOR, nbindex, ijndex)
    IMPLICIT NONE
!-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: iim, jjm, llm, mmm, zzm, itau
    INTEGER :: var(:,:,:,:,:)
    CHARACTER(LEN=*) :: MY_OPERATOR
    INTEGER :: nbindex, ijndex(nbindex)
    !-----------------------------
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: temp_g
    INTEGER :: ier

    ALLOCATE( temp_g(SIZE(var, DIM=1), SIZE(var, DIM=2), SIZE(var, DIM=3), &
                     SIZE(var, DIM=4), SIZE(var, DIM=5)), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3, 'restput_p_opp_i5d', 'Allocation memory error', vname_q, '')

    temp_g = REAL(var, r_std)          
    CALL restput_p &
         (fid, vname_q, iim, jjm, llm, mmm, zzm, itau, temp_g, MY_OPERATOR, nbindex, ijndex)

    DEALLOCATE( temp_g )
          
  END SUBROUTINE restput_p_opp_i5d

!!  =============================================================================================================================
!! SUBROUTINE:   restput_p_r1d
!!
!>\BRIEF	 allows to re-index data (real 1D) onto the original grid of the restart file
!!
!! DESCRIPTION:	 Need to be call by all process
!!
!! \n
!_ ==============================================================================================================================
  SUBROUTINE restput_p_r1d (fid,vname_q,iim,jjm,llm,itau,var)
    IMPLICIT NONE
!-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: iim, jjm, llm, itau
    REAL :: var(:)
    !-----------------------------
    REAL, ALLOCATABLE, DIMENSION(:) :: temp_g

    IF (is_root_prc) THEN
      ALLOCATE( temp_g(iim*jjm*llm) )
    ELSE
      ALLOCATE( temp_g(1) )
    ENDIF
    
    CALL gather2D_mpi(var,temp_g)
    IF (is_root_prc) THEN
       CALL restput (fid,vname_q,iim,jjm,llm,itau,temp_g)
    ENDIF
    DEALLOCATE( temp_g )
          
  END SUBROUTINE restput_p_r1d

!!  =============================================================================================================================
!! SUBROUTINE:   restput_p_r2d
!!
!>\BRIEF	 allows to re-index data (real 2D) onto the original grid of the restart file
!!
!! DESCRIPTION:	 Need to be call by all process
!!
!! \n
!_ ==============================================================================================================================
  SUBROUTINE restput_p_r2d (fid,vname_q,iim,jjm,llm,itau,var)
    IMPLICIT NONE
!-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: iim, jjm, llm, itau
    REAL :: var(:,:)
    !-------------------------
    REAL, ALLOCATABLE, DIMENSION(:,:) :: temp_g

    IF (is_root_prc) THEN
      ALLOCATE( temp_g(iim,jjm) )
    ELSE
      ALLOCATE( temp_g(1,1) )
    ENDIF
    
    CALL gather2D_mpi(var,temp_g)
    IF (is_root_prc) THEN
       CALL restput (fid,vname_q,iim,jjm,llm,itau,temp_g)
    ENDIF
    DEALLOCATE( temp_g )
          
  END SUBROUTINE restput_p_r2d

!!  =============================================================================================================================
!! SUBROUTINE:   restput_p_nogrid_r1d
!!
!>\BRIEF	  save reald 1D array (non-grid) data into the restart file
!!
!! DESCRIPTION:	 Need to be call by all process
!!
!! \n
!_ ==============================================================================================================================
  SUBROUTINE restput_p_nogrid_r1d (fid,vname_q,itau,var)
    IMPLICIT NONE
!-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: itau
    REAL,DIMENSION(:) :: var
    !-----------------------------

    IF (is_root_prc) THEN
       CALL restput (fid, vname_q, 1, 1, 1, itau, var)
    ENDIF
          
  END SUBROUTINE restput_p_nogrid_r1d

!!  =============================================================================================================================
!! SUBROUTINE:   restput_p_r3d
!!
!>\BRIEF	  allows to re-index data (real 3D) onto the original grid of the restart file
!!
!! DESCRIPTION:	 Need to be call by all process
!!
!! \n
!_ ==============================================================================================================================
  SUBROUTINE restput_p_r3d (fid,vname_q,iim,jjm,llm,itau,var)
    IMPLICIT NONE
!-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: iim, jjm, llm, itau
    REAL :: var(:,:,:)
    !-------------------------
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: temp_g

    IF (is_root_prc) THEN
      ALLOCATE( temp_g(iim,jjm,llm) )
    ELSE
      ALLOCATE( temp_g(1,1,1) )
    ENDIF
    
    CALL gather2D_mpi(var,temp_g)
    IF (is_root_prc) THEN
       CALL restput (fid,vname_q,iim,jjm,llm,itau,temp_g)
    ENDIF
    DEALLOCATE( temp_g )
          
  END SUBROUTINE restput_p_r3d

!!  =============================================================================================================================
!! SUBROUTINE:   restput_p_nogrid_r_scal
!!
!>\BRIEF	  save real scalar (non-grid) data into the restart file
!!
!! DESCRIPTION:	 Need to be call by all process
!!
!! \n
!_ ==============================================================================================================================
  SUBROUTINE restput_p_nogrid_r_scal (fid,vname_q,itau,var)
    IMPLICIT NONE
!-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: itau
    REAL :: var
    !-----------------------------
    REAL :: xtmp(1)

    IF (is_root_prc) THEN
       xtmp(1) = var
       CALL restput (fid, vname_q, 1, 1, 1, itau, xtmp)
    ENDIF
          
  END SUBROUTINE restput_p_nogrid_r_scal

!!  =============================================================================================================================
!! SUBROUTINE:   restput_p_nogrid_i_scal
!!
!>\BRIEF	  save integer scalar (non-grid) data into the restart file
!!
!! DESCRIPTION:	 Need to be call by all process
!!
!! \n
!_ ==============================================================================================================================
  SUBROUTINE restput_p_nogrid_i_scal (fid,vname_q,itau,var)
    IMPLICIT NONE
!-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: itau
    INTEGER :: var
    !-----------------------------
    REAL :: xtmp(1)
    REAL :: realvar

    IF (is_root_prc) THEN
       realvar = REAL(var,r_std)
       xtmp(1) = realvar
       CALL restput (fid, vname_q, 1, 1, 1, itau, xtmp)
    ENDIF
          
  END SUBROUTINE restput_p_nogrid_i_scal

!!  =============================================================================================================================
!! SUBROUTINE:   histwrite_r1d_p
!!
!>\BRIEF   give the data (real 1D) to the IOIPSL system (if we don't use XIOS). 	 
!!
!! DESCRIPTION:	 Need to be call by all process
!!
!! \n
!_ ==============================================================================================================================
  SUBROUTINE histwrite_r1d_p(pfileid,pvarname,pitau,pdata,nbindex,nindex)
    IMPLICIT NONE
!-
    INTEGER,INTENT(IN) :: pfileid, pitau, nbindex, nindex(nbindex)
    REAL,DIMENSION(:),INTENT(IN) :: pdata
    CHARACTER(LEN=*),INTENT(IN) :: pvarname
    
    REAL,DIMENSION(nbp_mpi)    :: pdata_mpi
    
    IF (pfileid > 0) THEN 
       ! Continue only if the file is initilalized
       CALL gather_omp(pdata,pdata_mpi)
       IF (is_omp_root) THEN
          CALL histwrite(pfileid,pvarname,pitau,pdata_mpi,nbp_mpi,kindex_mpi) 
       ENDIF
    END IF
      
  END SUBROUTINE histwrite_r1d_p
  
!!  =============================================================================================================================
!! SUBROUTINE:   histwrite_r2d_p
!!
!>\BRIEF	  give the data (real 2D) to the IOIPSL system (if we don't use XIOS). 	 
!!
!! DESCRIPTION:	 Need to be call by all process
!!
!! \n
!_ ==============================================================================================================================
  SUBROUTINE histwrite_r2d_p(pfileid,pvarname,pitau,pdata,nbindex,nindex)
    IMPLICIT NONE
!-
    INTEGER,INTENT(IN) :: pfileid, pitau, nbindex, nindex(nbindex)
    REAL,DIMENSION(:,:),INTENT(IN) :: pdata
    CHARACTER(LEN=*),INTENT(IN) :: pvarname

    IF (pfileid > 0) THEN 
       ! Continue only if the file is initilalized
       CALL body(size(pdata,2),nindex)
    END IF

  CONTAINS 

    SUBROUTINE body(dim,nindex)
    INTEGER :: dim
    INTEGER :: nindex(nbp_omp,dim)
    
    INTEGER :: nindex_mpi(nbp_mpi,dim)
    REAL    :: pdata_mpi(nbp_mpi,dim)
    INTEGER    :: flat_nindex_mpi(nbp_mpi * dim)
    
      CALL gather_omp(pdata,pdata_mpi)
      CALL gather_omp(nindex,nindex_mpi)
    
      IF (is_omp_root) THEN
       flat_nindex_mpi(:) = reshape(nindex_mpi,(/nbp_mpi*dim/))
       CALL histwrite(pfileid,pvarname,pitau,pdata_mpi,nbp_mpi*dim,flat_nindex_mpi) 
      ENDIF
    END SUBROUTINE body
       
  END SUBROUTINE histwrite_r2d_p

!!  =============================================================================================================================
!! SUBROUTINE:   histwrite_r3d_p
!!
!>\BRIEF      give the data (real 3D) to the IOIPSL system (if we don't use XIOS).
!!
!! DESCRIPTION:	 Need to be call by all process
!!
!! \n
!_ ==============================================================================================================================
  SUBROUTINE histwrite_r3d_p(pfileid,pvarname,pitau,pdata,nbindex,nindex)
    IMPLICIT NONE
!-
    INTEGER,INTENT(IN) :: pfileid, pitau, nbindex, nindex(nbindex)
    REAL,DIMENSION(:,:,:),INTENT(IN) :: pdata
    CHARACTER(LEN=*),INTENT(IN) :: pvarname

    CHARACTER(LEN=10) :: part_str
    CHARACTER(LEN=LEN(part_str) + LEN(pvarname) + 1) :: var_name
    REAL,DIMENSION(SIZE(pdata, 1),SIZE(pdata, 2)) :: tmparr
    INTEGER :: jv

    DO jv = 1, SIZE(pdata, 3)
       WRITE(part_str,'(I2)') jv
       IF (jv < 10) part_str(1:1) = '0'
       var_name = TRIM(pvarname)//'_'//part_str(1:LEN_TRIM(part_str))
       tmparr = pdata(:,:,jv)
       CALL histwrite_r2d_p(pfileid, var_name, pitau, tmparr, nbindex, nindex)
    ENDDO
  
    
  END SUBROUTINE histwrite_r3d_p


END MODULE ioipsl_para
