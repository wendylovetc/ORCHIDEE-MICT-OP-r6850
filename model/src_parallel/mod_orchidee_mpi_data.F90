! ==============================================================================================================================
! MODULE   : mod_orchidee_mpi_data
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF      All routines needed for the initialization of the MPI parallelization of Orchidee (coupled and offline) 
!!
!!\n DESCRIPTION  :  Definition and allocation of parallel datas for MPI.
!!                   Initialization of parallel or sequentiel IOs.
!!                   Definition of Load Balancing functions.
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


MODULE mod_orchidee_mpi_data

!-
  USE defprec
  USE ioipsl
  USE ioipsl_para
  USE xios_orchidee
  USE mod_orchidee_para_var

  IMPLICIT NONE

!-
#include "src_parallel.h"
!-
!-

CONTAINS

  !!  =============================================================================================================================
  !! SUBROUTINE:  Init_orchidee_mpi
  !!
  !>\BRIEF	 Initialization of the mpi parallelization in ORCHIDEE offline mode. 
  !!
  !! DESCRIPTION:	 Initialization of the mpi parallelization in ORCHIDEE offline mode. 
  !!
  !! \n
  !_ ============================================================================================================================== 
  SUBROUTINE Init_orchidee_mpi(communicator)

#ifdef CPP_PARA
    USE mpi
#endif

    IMPLICIT NONE
    INTEGER, OPTIONAL, INTENT(in) :: communicator
    INTEGER :: COMM
    INTEGER :: ierr

    !Config Key   = XIOS_ORCHIDEE_OK
    !Config Desc  = Use XIOS for writing diagnostics file
    !Config If    = 
    !Config Def   = y
    !Config Help  = If XIOS_ORCHIDEE_OK=y then no output with IOIPSL can be done
    !Config Units = [FLAG]
    CALL getin('XIOS_ORCHIDEE_OK',xios_orchidee_ok)

#ifdef CPP_PARA
    IF ( xios_orchidee_ok ) THEN
       CALL MPI_INIT(ierr)
       CALL xios_orchidee_comm_init(COMM)
    ELSE IF ( PRESENT(communicator) ) THEN
       COMM=communicator
    ELSE
       CALL MPI_INIT(ierr)
       COMM=MPI_COMM_WORLD
    ENDIF
    CALL MPI_COMM_SIZE(COMM,mpi_size,ierr)
    CALL MPI_COMM_RANK(COMM,mpi_rank,ierr)
    is_ok_mpi=.TRUE.
#else
    mpi_rank=0
    mpi_size=1
    is_ok_mpi=.FALSE.
    ! It is not possible to use XIOS without MPI
    WRITE(numout,*) 'XIOS cannot be run without MPI. xios_orchidee_ok is set to false.'
    xios_orchidee_ok=.FALSE.
#endif
    
    mpi_rank_root=0

    IF (mpi_rank==mpi_rank_root) THEN 
      is_mpi_root=.TRUE.
    ELSE
      is_mpi_root=.FALSE.
    ENDIF
  
    CALL Init_const_mpi(COMM)
      
  END SUBROUTINE Init_orchidee_mpi



  !!  =============================================================================================================================
  !! SUBROUTINE:  Init_orchidee_mpi_data
  !!
  !>\BRIEF	 Initialization of the mpi parallelization  
  !!
  !! DESCRIPTION: Initialization of the mpi parallelization. This subroutine is called from Init_orchidee_data_para 
  !!              which is called both in offline and coupled mode.
  !!
  !! \n
  !_ ============================================================================================================================== 
  SUBROUTINE Init_orchidee_mpi_data(arg_nbp_mpi,arg_kindex_mpi,arg_offset_mpi,COMM)
#ifdef CPP_PARA
    USE mpi
#endif

  IMPLICIT NONE
    INTEGER, INTENT(IN) :: arg_nbp_mpi
    INTEGER, INTENT(IN) :: arg_kindex_mpi(arg_nbp_mpi)
    INTEGER, INTENT(IN) :: arg_offset_mpi
    INTEGER, INTENT(IN) :: COMM

    INTEGER :: i
  
#ifdef CPP_PARA
    INTEGER :: ierr
    is_ok_mpi=.TRUE.
#else
    is_ok_mpi=.FALSE.
#endif

    ! Initialization of MPI_COMM_ORCH 
    CALL init_const_mpi(COMM)
    IF (is_ok_mpi) THEN    
#ifdef CPP_PARA
      CALL MPI_COMM_SIZE(MPI_COMM_ORCH,mpi_size,ierr)    
      CALL MPI_COMM_RANK(MPI_COMM_ORCH,mpi_rank,ierr)
#endif
    ELSE
      mpi_size=1
      mpi_rank=0
    ENDIF   
    
    IF (mpi_rank == 0) THEN
      mpi_rank_root = 0
      is_mpi_root = .true.
    ENDIF
    
    ! Test if there is enough land grid cells on each process MPI.
    ! At least 2 grid cells are needed for each process if running in parallel. Stop if this is not the case.
    IF ( arg_nbp_mpi < 1 ) THEN
       PRINT*,'Init_orchidee_mpi_data: nbp_mpi(number of grid-cells for current MPI process)=',arg_nbp_mpi
       PRINT*,'nbp_mpi=0 is not possible. It must be 1 or bigger for each MPI process. Stop now.'
       CALL ipslerr_p(3, "Init_orchidee_mpi_data","Not all MPI processes has enough land grid cells",&
            "Make the region bigger or use a lower number of MPI processes.","")
    ELSE IF (arg_nbp_mpi == 1 .AND. mpi_size > 1) THEN
       CALL ipslerr_p(1, "Init_orchidee_mpi_data","nbp=1 for current MPI process",&
            "This can be a problem in some case.", &
            "If the model crashes, then make the region bigger or use a lower number of MPI processes.")
    END IF

    ALLOCATE(nbp_mpi_para(0:mpi_size-1))
    ALLOCATE(nbp_mpi_para_begin(0:mpi_size-1))
    ALLOCATE(nbp_mpi_para_end(0:mpi_size-1))    
    ALLOCATE(jj_para_nb(0:mpi_size-1))
    ALLOCATE(jj_para_begin(0:mpi_size-1))
    ALLOCATE(jj_para_end(0:mpi_size-1))
    ALLOCATE(ii_para_begin(0:mpi_size-1))
    ALLOCATE(ii_para_end(0:mpi_size-1))    
    ALLOCATE(ij_para_nb(0:mpi_size-1))
    ALLOCATE(ij_para_begin(0:mpi_size-1))
    ALLOCATE(ij_para_end(0:mpi_size-1))
    

    nbp_mpi=arg_nbp_mpi
    ALLOCATE(kindex_mpi(nbp_mpi))
    kindex_mpi(:)=arg_kindex_mpi(:)
    
    offset_mpi=arg_offset_mpi
    
    IF (is_ok_mpi) THEN
#ifdef CPP_PARA
      CALL MPI_AllGather(nbp_mpi,1,MPI_INT_ORCH,nbp_mpi_para,1,MPI_INT_ORCH,MPI_COMM_ORCH,ierr)
#endif
    ELSE
      nbp_mpi_para(0)=nbp_mpi
    ENDIF
    
    nbp_mpi_para_begin(0)=1
    nbp_mpi_para_end(0)=nbp_mpi_para(0)
    DO i=1,mpi_size-1
      nbp_mpi_para_begin(i)=nbp_mpi_para_end(i-1)+1
      nbp_mpi_para_end(i)=nbp_mpi_para_begin(i)+nbp_mpi_para(i)-1
    ENDDO
    nbp_mpi_begin=nbp_mpi_para_begin(mpi_rank)
    nbp_mpi_end=nbp_mpi_para_end(mpi_rank)
    
    
    IF (mpi_rank==mpi_size-1) THEN
      ij_end=iim_g*jjm_g
    ELSE
      ij_end=kindex_mpi(nbp_mpi)+offset_mpi
    ENDIF

    IF (is_ok_mpi) THEN    
#ifdef CPP_PARA    
      CALL MPI_Allgather(ij_end,1,MPI_INT_ORCH,ij_para_end,1,MPI_INT_ORCH,MPI_COMM_ORCH,ierr)
#endif
    ELSE
      ij_para_end(0)=ij_end
    ENDIF
    
    ij_para_begin(0)=1
    ij_para_nb(0)=ij_para_end(0)-ij_para_begin(0)+1
    
    DO i=1,mpi_size-1
      ij_para_begin(i)=ij_para_end(i-1)+1
      ij_para_nb(i)=ij_para_end(i)-ij_para_begin(i)+1
    ENDDO
    
    DO i=0,mpi_size-1
      jj_para_begin(i)=(ij_para_begin(i)-1)/iim_g + 1
      jj_para_end(i)=(ij_para_end(i)-1)/iim_g + 1
      jj_para_nb(i)=jj_para_end(i)-jj_para_begin(i)+1
          
      ii_para_begin(i)=MOD(ij_para_begin(i)-1,iim_g)+1
      ii_para_end(i)=MOD(ij_para_end(i)-1,iim_g)+1
    ENDDO

   
    ij_nb=ij_para_nb(mpi_rank)
    ij_begin=ij_para_begin(mpi_rank)
    ij_end=ij_para_end(mpi_rank)
        
    jj_nb=jj_para_nb(mpi_rank)
    jj_begin=jj_para_begin(mpi_rank)
    jj_end=jj_para_end(mpi_rank)
    
    ii_begin=ii_para_begin(mpi_rank)
    ii_end=ii_para_end(mpi_rank)
        
      
    CALL print_mpi_data
  
    
  END SUBROUTINE Init_orchidee_mpi_data
  
  !!  =============================================================================================================================
  !! SUBROUTINE:  init_const_mpi
  !!
  !>\BRIEF	 Initialization of some constants related to the MPI parallelization 
  !!
  !! DESCRIPTION: Initialization of some constants related to the MPI parallelization
  !!
  !! \n
  !_ ============================================================================================================================== 
  SUBROUTINE init_const_mpi(COMM)

#ifdef CPP_PARA
    USE mpi
#endif

  IMPLICIT NONE
    INTEGER :: COMM

#ifdef CPP_PARA
    MPI_COMM_ORCH=COMM
    
    IF (i_std==i_4) THEN
       MPI_INT_ORCH=MPI_INTEGER4
    ELSEIF (i_std==i_8) THEN
       MPI_INT_ORCH=MPI_INTEGER8
    ELSE
       MPI_INT_ORCH=MPI_INTEGER
    ENDIF
         
    IF (r_std==r_4) THEN
       MPI_REAL_ORCH=MPI_REAL4
    ELSEIF (r_std==r_8) THEN
       MPI_REAL_ORCH=MPI_REAL8
    ELSE
       MPI_REAL_ORCH=MPI_REAL
    ENDIF
#endif

  END SUBROUTINE init_const_mpi

  !!  =============================================================================================================================
  !! SUBROUTINE:  Finalize_mpi
  !!
  !>\BRIEF	 Close the MPI parallelization 
  !!
  !! DESCRIPTION:    Close the MPI parallelization. The context XIOS will be closed before call to MPI_finalize routine
  !!
  !! \n
  !_ ============================================================================================================================== 
  SUBROUTINE Finalize_mpi

#ifdef CPP_PARA
  USE mpi
#endif

  IMPLICIT NONE
#ifdef CPP_PARA
  INTEGER :: ierr

  CALL xios_orchidee_finalize

  CALL MPI_FINALIZE(ierr)
#endif
   
  END SUBROUTINE Finalize_mpi
  
  !!  =============================================================================================================================
  !! SUBROUTINE:  print_mpi_data
  !!
  !>\BRIEF	 Print all data specific to MPI parallelization of ORCHIDEE 
  !!
  !! DESCRIPTION:  Print all data specific to MPI parallelization of ORCHIDEE 	
  !!
  !! \n
  !_ ============================================================================================================================== 
  SUBROUTINE print_mpi_data

  IMPLICIT NONE
    
    WRITE(numout,*) '==== MPI DOMAIN ===='
    WRITE(numout,*) '     ----------     '
    WRITE(numout,*) 'mpi_size',mpi_size
    WRITE(numout,*) 'mpi_rank',mpi_rank
    WRITE(numout,*) 'is_mpi_root',is_mpi_root
    WRITE(numout,*) 'mpi_rank_root',mpi_rank_root

    WRITE(numout,*) 'nbp_mpi_begin=',nbp_mpi_begin
    WRITE(numout,*) 'nbp_mpi_end  =',nbp_mpi_end
    WRITE(numout,*) 'nbp_mpi=',nbp_mpi
          
    WRITE(numout,*) 'ij_begin=',ij_begin
    WRITE(numout,*) 'ij_end=',ij_end
    WRITE(numout,*) 'ij_nb=',ij_nb
    WRITE(numout,*) 'jj_begin=',jj_begin
    WRITE(numout,*) 'jj_end=',jj_end
    WRITE(numout,*) 'jj_nb=',jj_nb	
    WRITE(numout,*) 'ii_begin=',ii_begin
    WRITE(numout,*) 'ii_end=',ii_end
    
    WRITE(numout,*) 'offset_mpi',offset_mpi
    WRITE(numout,*) 'nbp_mpi_para_begin=',nbp_mpi_para_begin
    WRITE(numout,*) 'nbp_mpi_para_end  =',nbp_mpi_para_end
    WRITE(numout,*) 'nbp_mpi_para=',nbp_mpi_para
          
    WRITE(numout,*) 'ij_para_begin=',ij_para_begin
    WRITE(numout,*) 'ij_para_end=',ij_para_end
    WRITE(numout,*) 'ij_para_nb=',ij_para_nb
    WRITE(numout,*) 'jj_para_begin=',jj_para_begin
    WRITE(numout,*) 'jj_para_end=',jj_para_end
    WRITE(numout,*) 'jj_para_nb=',jj_para_nb	
    WRITE(numout,*) 'ii_para_begin=',ii_para_begin
    WRITE(numout,*) 'ii_para_end=',ii_para_end
  
  END SUBROUTINE print_mpi_data
  
  !!  =============================================================================================================================
  !! SUBROUTINE:  Read_Load_balance
  !!
  !>\BRIEF	 Read load balance file.
  !!
  !! DESCRIPTION:	Read load balance file. This is only done in offline mode.
  !!                The load balance file contains information about the MPI partitionning on the different processes.
  !!
  !! \n
  !_ ============================================================================================================================== 
 SUBROUTINE Read_Load_balance(NbPoints,Nbpoints_loc)

    IMPLICIT NONE
    INTEGER,INTENT(IN)  :: NbPoints
    INTEGER,INTENT(OUT) :: Nbpoints_loc(0:mpi_size-1)
    INTEGER :: i,s
    INTEGER :: ierr
    
#ifdef CPP_PARA
    CHARACTER(len=255)  :: filename='Load_balance_orchidee.dat'
    INTEGER :: j
    INTEGER :: unit_number=10
#endif   

#ifdef CPP_PARA
    OPEN(UNIT=unit_number,FILE=trim(filename),STATUS='old',FORM='formatted',IOSTAT=ierr) 
#else
    ierr=1
#endif   
    Nbpoints_loc(:) = 0

    s=0
#ifdef CPP_PARA  
    IF (ierr==0) THEN
       i=0
       !- Reading for any balancing file (even with a bad structure)
       DO WHILE (i < mpi_size .AND. ierr == 0) 
          READ (unit_number,*,IOSTAT=ierr) j,Nbpoints_loc(i)
          s=s+Nbpoints_loc(i)
          i=i+1
       ENDDO
       CLOSE(unit_number)
    ENDIF
#endif   
    
    !- Correction of bad balancing file (or an empty file) => same nb of points for each procs
    IF (ierr/=0 .OR. s/=Nbpoints) THEN
       DO i=0,mpi_size-1
          Nbpoints_loc(i)=Nbpoints/mpi_size
          IF (MOD(Nbpoints,mpi_size) > i) Nbpoints_loc(i)=Nbpoints_loc(i)+1
       ENDDO
    ENDIF
    
  END SUBROUTINE Read_Load_balance
  
  !!  =============================================================================================================================
  !! SUBROUTINE:  Write_Load_balance
  !!
  !>\BRIEF	 Write the load balance file.
  !!
  !! DESCRIPTION:	Write the load balance file. This is only done in offline mode. 
  !!                The load balance file contains information about the MPI partitionning on the different processes.
  !!
  !! \n
  !_ ============================================================================================================================== 
  SUBROUTINE Write_Load_balance(times)
    IMPLICIT NONE
    REAL,INTENT(IN) :: times
  
    CHARACTER(len=255)  :: filename='Load_balance_orchidee.dat'
    INTEGER :: unit_number=10
    INTEGER :: i,ierr
    REAL :: All_Times(0:mpi_size-1)
    REAL :: average
    REAL :: efficiency
    INTEGER :: dp,S
    INTEGER :: New_nbpoints(0:mpi_size-1)
    
    IF (is_ok_mpi) THEN
#ifdef CPP_PARA
      CALL MPI_GATHER(times,1,MPI_REAL_ORCH,All_times,1,MPI_REAL_ORCH,mpi_rank_root,MPI_COMM_ORCH,ierr)
#endif
    ELSE
      All_times(:)=times
    ENDIF
    
    IF (is_mpi_root) WRITE(numout,*) 'ALL_times',All_times

    IF (is_mpi_root) THEN
     
       OPEN(UNIT=unit_number,FILE=trim(filename),STATUS='replace',FORM='formatted',IOSTAT=ierr)
       
       average=sum(All_times(:))/mpi_size
       DO i=0,mpi_size-1
          efficiency=All_times(i)/nbp_mpi_para(i)
          New_nbpoints(i)=Nbp_mpi_para(i)-(All_times(i)-average)/efficiency
       ENDDO
       
       S=sum(new_nbpoints(:))
       dp=nbp_glo-S
       
       IF ( dp > 0 ) THEN
          DO WHILE ( dp > 0 )
             New_nbpoints(MOD(dp,mpi_size))=New_nbpoints(MOD(dp,mpi_size))+1
             dp=dp-1
          ENDDO
       ELSE
          dp=-dp
          DO WHILE ( dp > 0 )
             New_nbpoints(MOD(dp,mpi_size))=New_nbpoints(MOD(dp,mpi_size))-1
             dp=dp-1
          ENDDO
       ENDIF
       

       ! If this algorithm diverge, we use previous repartition.
       IF ( ANY(New_nbpoints(:) .LE. 0) ) THEN
          New_nbpoints(:)=Nbp_mpi_para(:)
       ENDIF
       
       DO i=0,mpi_size-1
          WRITE(Unit_number,*) i,New_nbpoints(i)
       ENDDO
       CLOSE(Unit_number)
    ENDIF

  END SUBROUTINE Write_Load_Balance
  
END MODULE mod_orchidee_mpi_data

#include "mpi_dummy.h"
