! ===============================================================================================================================
! MODULE       : mod_orchidee_transfert_para
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        Contains the interfaces bcast, scatter, gather and reduce_sum
!!
!! \n DESCRIPTION : High level MPI/OpenMP parallel communication encapsulations for ORCHIDEE.
!!                  These interfaces should be called by all procesess (all MPI processes and all OMP threads)
!!          bcast      : send a variable from the root process to all processes
!!          scatter    : distribute a global field known on the root process to the local domain for each processes
!!          gather     : each process sends their local field to the root process which will recieve the field on the global domain
!!          reduce_sum : the root process will recieve the sum of the values from all processes
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

MODULE mod_orchidee_transfert_para

  USE mod_orchidee_para_var
  USE mod_orchidee_mpi_transfert
  USE mod_orchidee_omp_transfert 
  

! ===============================================================================================================================
!! INTERFACE    : bcast
!!
!>\BRIEF        Send a variable from root process to all processes
!!
!! \n DESCRIPTION : Send a variable from root process to all processes. Note that all processes must make the call.
!!
!!  Usage : CALL bcast(Var)
!!  The variable Var must be known before the call on the root process. After the call, all processes know the variable. 
!!  The variable has the same dimension on all processes. The variable can be a character string, integer, real or logical. 
!!  It can have rank 1 to 4. 
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCES(S)    : None
!!
!! \n
!_ ================================================================================================================================

  INTERFACE bcast
    MODULE PROCEDURE bcast_c, bcast_c1,                           &
                     bcast_i,bcast_i1,bcast_i2,bcast_i3,bcast_i4, &
                     bcast_r,bcast_r1,bcast_r2,bcast_r3,bcast_r4, &
             bcast_l,bcast_l1,bcast_l2,bcast_l3,bcast_l4
  END INTERFACE

! ===============================================================================================================================
!! INTERFACE    : scatter
!!
!>\BRIEF        Distribute a global field from the root process to the local domain on each processes
!!
!! \n DESCRIPTION : Distribute the global field known on the root process to the local domain for each processes.
!!                  Note that all processes must make the call.
!!
!! Usage: CALL scatter(VarIn, VarOut)
!! VarIn must be known on the root process before the call is done. 
!! After the call, VarOut contains the variable on the local domain on each process.
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCES(S)    : None
!!
!! \n
!_ ================================================================================================================================
  INTERFACE scatter
    MODULE PROCEDURE scatter_i,scatter_i1,scatter_i2,scatter_i3,scatter_i4, &
                     scatter_r,scatter_r1,scatter_r2,scatter_r3,scatter_r4, &
            scatter_l,scatter_l1,scatter_l2,scatter_l3,scatter_l4
  END INTERFACE

  
! ===============================================================================================================================
!! INTERFACE    : gather
!!
!>\BRIEF        Each process send their local field to the root process which will recieve the field on the global domain
!!
!! \n DESCRIPTION : Each process send their local field to the root process which will recieve the field on the global domain.
!!                  Note that all processes must make the call.
!!
!! Usage: CALL gather(VarIn, VarOut)
!! VarIn is the variable on the local domain on each process, known before the call. After the call, VarOut is recieved on the
!! root process containing the variable on the global domain.
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCES(S)    : None
!!
!! \n
!_ ================================================================================================================================
  INTERFACE gather
    MODULE PROCEDURE gather_i,gather_i1,gather_i2,gather_i3,gather_i4, &
                     gather_r,gather_r1,gather_r2,gather_r3,gather_r4, &
          gather_l,gather_l1,gather_l2,gather_l3,gather_l4  
  END INTERFACE
  

! ===============================================================================================================================
!! INTERFACE    : reduce_sum
!!
!>\BRIEF        The root process will recieve the sum of the values from all processes
!!
!! \n DESCRIPTION : The root process will recieve the sum of the values from all processes
!!                  Note that all processes must make the call.
!!
!! Usage: CALL reduce_sum(VarIn, VarOut)
!! VarIn is the value on each process. VarOut is the sum of these values on the root process.
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCES(S)    : None
!!x
!! \n
!_ ================================================================================================================================
  INTERFACE reduce_sum
    MODULE PROCEDURE reduce_sum_i,reduce_sum_i1,reduce_sum_i2,reduce_sum_i3,reduce_sum_i4, &
                     reduce_sum_r,reduce_sum_r1,reduce_sum_r2,reduce_sum_r3,reduce_sum_r4
  END INTERFACE 


! ===============================================================================================================================
!! INTERFACE    : all_reduce_sum
!!
!>\BRIEF        All processes will recieve the sum of the values from all processes
!!
!! \n DESCRIPTION : All processes will recieve the sum of the values from all processes
!!                  Note that all processes must make the call.
!!
!! Usage: CALL all_reduce_sum(VarIn, VarOut)
!! VarIn is the value on each process. VarOut is the sum of all these values on the all processes.
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCES(S)    : None
!!x
!! \n
!_ ================================================================================================================================
  INTERFACE allreduce_sum
    MODULE PROCEDURE allreduce_sum_r
  END INTERFACE 
   
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Definition of Broadcast 1D --> 4D !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! -- Character string -- !!

  SUBROUTINE bcast_c(var)
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(INOUT) :: Var
   
    
    IF (is_omp_root) CALL bcast_mpi(Var)
    CALL bcast_omp(Var)
    
  END SUBROUTINE bcast_c

  SUBROUTINE bcast_c1(var)
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(INOUT) :: Var(:)
   
    
    IF (is_omp_root) CALL bcast_mpi(Var)
    CALL bcast_omp(Var)
    
  END SUBROUTINE bcast_c1

!! -- Integers -- !!
  
  SUBROUTINE bcast_i(var)
  IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: Var
    
    IF (is_omp_root) CALL bcast_mpi(Var)
    CALL bcast_omp(Var)
    
  END SUBROUTINE bcast_i

  SUBROUTINE bcast_i1(var)
  IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: Var(:)
   
    
    IF (is_omp_root) CALL bcast_mpi(Var)
    CALL bcast_omp(Var)
    
  END SUBROUTINE bcast_i1


  SUBROUTINE bcast_i2(var)
  IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: Var(:,:)
   
    
    IF (is_omp_root) CALL bcast_mpi(Var)
    CALL bcast_omp(Var)
    
  END SUBROUTINE bcast_i2


  SUBROUTINE bcast_i3(var)
  IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: Var(:,:,:)
   
    
    IF (is_omp_root) CALL bcast_mpi(Var)
    CALL bcast_omp(Var)
    
  END SUBROUTINE bcast_i3


  SUBROUTINE bcast_i4(var)
  IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: Var(:,:,:,:)
   
    
    IF (is_omp_root) CALL bcast_mpi(Var)
    CALL bcast_omp(Var)
    
  END SUBROUTINE bcast_i4

 
!! -- Reals -- !!
  
  SUBROUTINE bcast_r(var)
  IMPLICIT NONE
    REAL,INTENT(INOUT) :: Var

    
    IF (is_omp_root) CALL bcast_mpi(Var)
    CALL bcast_omp(Var)
    
  END SUBROUTINE bcast_r

  SUBROUTINE bcast_r1(var)
  IMPLICIT NONE
    REAL,INTENT(INOUT) :: Var(:)
   
    
    IF (is_omp_root) CALL bcast_mpi(Var)
    CALL bcast_omp(Var)
    
  END SUBROUTINE bcast_r1


  SUBROUTINE bcast_r2(var)
  IMPLICIT NONE
    REAL,INTENT(INOUT) :: Var(:,:)
   
    
    IF (is_omp_root) CALL bcast_mpi(Var)
    CALL bcast_omp(Var)
    
  END SUBROUTINE bcast_r2


  SUBROUTINE bcast_r3(var)
  IMPLICIT NONE
    REAL,INTENT(INOUT) :: Var(:,:,:)
   
    
    IF (is_omp_root) CALL bcast_mpi(Var)
    CALL bcast_omp(Var)
    
  END SUBROUTINE bcast_r3


  SUBROUTINE bcast_r4(var)
  IMPLICIT NONE
    REAL,INTENT(INOUT) :: Var(:,:,:,:)
    
    IF (is_omp_root) CALL bcast_mpi(Var)
    CALL bcast_omp(Var)
    
  END SUBROUTINE bcast_r4 


!! -- Logicals -- !!
  
  SUBROUTINE bcast_l(var)
  IMPLICIT NONE
    LOGICAL,INTENT(INOUT) :: Var
    
    IF (is_omp_root) CALL bcast_mpi(Var)
    CALL bcast_omp(Var)
    
  END SUBROUTINE bcast_l

  SUBROUTINE bcast_l1(var)
  IMPLICIT NONE
    LOGICAL,INTENT(INOUT) :: Var(:)
   
    IF (is_omp_root) CALL bcast_mpi(Var)
    CALL bcast_omp(Var)
    
  END SUBROUTINE bcast_l1


  SUBROUTINE bcast_l2(var)
  IMPLICIT NONE
    LOGICAL,INTENT(INOUT) :: Var(:,:)
   
    
    IF (is_omp_root) CALL bcast_mpi(Var)
    CALL bcast_omp(Var)
    
  END SUBROUTINE bcast_l2


  SUBROUTINE bcast_l3(var)
  IMPLICIT NONE
    LOGICAL,INTENT(INOUT) :: Var(:,:,:)
    
    IF (is_omp_root) CALL bcast_mpi(Var)
    CALL bcast_omp(Var)
    
  END SUBROUTINE bcast_l3


  SUBROUTINE bcast_l4(var)
  IMPLICIT NONE
    LOGICAL,INTENT(INOUT) :: Var(:,:,:,:)
     
    IF (is_omp_root) CALL bcast_mpi(Var)

    CALL bcast_omp(Var)
    
  END SUBROUTINE bcast_l4


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Definition of Scatter  1D --> 4D  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE scatter_i(VarIn, VarOut)
  IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:) :: VarOut

    INTEGER,DIMENSION(nbp_mpi) :: Var_tmp
    
    
    IF (is_omp_root) CALL scatter_mpi(VarIn,Var_tmp)
    CALL scatter_omp(Var_tmp,Varout)
    
  END SUBROUTINE scatter_i


  SUBROUTINE scatter_i1(VarIn, VarOut)
  IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:) :: VarOut

    INTEGER,DIMENSION(nbp_mpi,SIZE(Varout,2)) :: Var_tmp

    
    IF (is_omp_root) CALL scatter_mpi(VarIn,Var_tmp)
    CALL scatter_omp(Var_tmp,Varout)
    
  END SUBROUTINE scatter_i1


  SUBROUTINE scatter_i2(VarIn, VarOut)
  IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    INTEGER,DIMENSION(nbp_mpi,SIZE(Varout,2),SIZE(Varout,3)) :: Var_tmp

    
    IF (is_omp_root) CALL scatter_mpi(VarIn,Var_tmp)
    CALL scatter_omp(Var_tmp,Varout)
    
  END SUBROUTINE scatter_i2


  SUBROUTINE scatter_i3(VarIn, VarOut)
  IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut

    INTEGER,DIMENSION(nbp_mpi,SIZE(Varout,2),SIZE(Varout,3),SIZE(Varout,4)) :: Var_tmp

    
    IF (is_omp_root) CALL scatter_mpi(VarIn,Var_tmp)
    CALL scatter_omp(Var_tmp,VarOut)
    
  END SUBROUTINE scatter_i3


  SUBROUTINE scatter_i4(VarIn, VarOut)
  IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:,:,:) :: VarOut

    INTEGER,DIMENSION(nbp_mpi,SIZE(Varout,2),SIZE(Varout,3),SIZE(Varout,4),SIZE(Varout,5)) :: Var_tmp

    
    IF (is_omp_root) CALL scatter_mpi(VarIn,Var_tmp)
    CALL scatter_omp(Var_tmp,VarOut)
    
  END SUBROUTINE scatter_i4

  SUBROUTINE scatter_r(VarIn, VarOut)
  IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:) :: VarOut

    REAL,DIMENSION(nbp_mpi) :: Var_tmp
    
    
    IF (is_omp_root) CALL scatter_mpi(VarIn,Var_tmp)
    CALL scatter_omp(Var_tmp,Varout)
    
  END SUBROUTINE scatter_r


  SUBROUTINE scatter_r1(VarIn, VarOut)
  IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:) :: VarOut

    REAL,DIMENSION(nbp_mpi,SIZE(Varout,2)) :: Var_tmp

    
    IF (is_omp_root) CALL scatter_mpi(VarIn,Var_tmp)
    CALL scatter_omp(Var_tmp,Varout)
    
  END SUBROUTINE scatter_r1


  SUBROUTINE scatter_r2(VarIn, VarOut)
  IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    REAL,DIMENSION(nbp_mpi,SIZE(Varout,2),SIZE(Varout,3)) :: Var_tmp

    
    IF (is_omp_root) CALL scatter_mpi(VarIn,Var_tmp)
    CALL scatter_omp(Var_tmp,Varout)
    
  END SUBROUTINE scatter_r2


  SUBROUTINE scatter_r3(VarIn, VarOut)
  IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut

    REAL,DIMENSION(nbp_mpi,SIZE(Varout,2),SIZE(Varout,3),SIZE(Varout,4)) :: Var_tmp

    
    IF (is_omp_root) CALL scatter_mpi(VarIn,Var_tmp)
    CALL scatter_omp(Var_tmp,VarOut)
    
  END SUBROUTINE scatter_r3
  

  SUBROUTINE scatter_r4(VarIn, VarOut)
  IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:,:,:) :: VarOut

    REAL,DIMENSION(nbp_mpi,SIZE(Varout,2),SIZE(Varout,3),SIZE(Varout,4),SIZE(Varout,5)) :: Var_tmp

    
    IF (is_omp_root) CALL scatter_mpi(VarIn,Var_tmp)
    CALL scatter_omp(Var_tmp,VarOut)
    
  END SUBROUTINE scatter_r4
  

  SUBROUTINE scatter_l(VarIn, VarOut)
  IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:) :: VarOut

    LOGICAL,DIMENSION(nbp_mpi) :: Var_tmp
    
    
    IF (is_omp_root) CALL scatter_mpi(VarIn,Var_tmp)
    CALL scatter_omp(Var_tmp,Varout)
    
  END SUBROUTINE scatter_l


  SUBROUTINE scatter_l1(VarIn, VarOut)
  IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:) :: VarOut

    LOGICAL,DIMENSION(nbp_mpi,SIZE(Varout,2)) :: Var_tmp

    
    IF (is_omp_root) CALL scatter_mpi(VarIn,Var_tmp)
    CALL scatter_omp(Var_tmp,Varout)
    
  END SUBROUTINE scatter_l1


  SUBROUTINE scatter_l2(VarIn, VarOut)
  IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    LOGICAL,DIMENSION(nbp_mpi,SIZE(Varout,2),SIZE(Varout,3)) :: Var_tmp

    
    IF (is_omp_root) CALL scatter_mpi(VarIn,Var_tmp)
    CALL scatter_omp(Var_tmp,Varout)
    
  END SUBROUTINE scatter_l2


  SUBROUTINE scatter_l3(VarIn, VarOut)
  IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut

    LOGICAL,DIMENSION(nbp_mpi,SIZE(Varout,2),SIZE(Varout,3),SIZE(Varout,4)) :: Var_tmp

    
    IF (is_omp_root) CALL scatter_mpi(VarIn,Var_tmp)
    CALL scatter_omp(Var_tmp,VarOut)
    
  END SUBROUTINE scatter_l3

  SUBROUTINE scatter_l4(VarIn, VarOut)
  IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:,:) :: VarOut

    LOGICAL,DIMENSION(nbp_mpi,SIZE(Varout,2),SIZE(Varout,3),SIZE(Varout,4),SIZE(Varout,5)) :: Var_tmp

    
    IF (is_omp_root) CALL scatter_mpi(VarIn,Var_tmp)
    CALL scatter_omp(Var_tmp,VarOut)
    
  END SUBROUTINE scatter_l4


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Definition of Gather  1D --> 4D   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
!!!!! --> Integers

  SUBROUTINE gather_i(VarIn, VarOut)
  IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:) :: VarOut
    
    INTEGER, DIMENSION(nbp_mpi) :: Var_tmp
    
    CALL gather_omp(VarIn,Var_tmp)
    
    IF (is_omp_root) CALL gather_mpi(Var_tmp,Varout)
  
  END SUBROUTINE gather_i


  SUBROUTINE gather_i1(VarIn, VarOut)
  IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    INTEGER, DIMENSION(nbp_mpi,SIZE(VarIn,2)) :: Var_tmp
    
    CALL gather_omp(VarIn,Var_tmp)
    
    IF (is_omp_root) CALL gather_mpi(Var_tmp,Varout)

  
  END SUBROUTINE gather_i1


  SUBROUTINE gather_i2(VarIn, VarOut)
  IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    INTEGER, DIMENSION(nbp_mpi,SIZE(VarIn,2),SIZE(VarIn,3)) :: Var_tmp
    
    CALL gather_omp(VarIn,Var_tmp)
    
    IF (is_omp_root) CALL gather_mpi(Var_tmp,VarOut)
  
  END SUBROUTINE gather_i2


  SUBROUTINE gather_i3(VarIn, VarOut)
  IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    INTEGER, DIMENSION(nbp_mpi,SIZE(VarIn,2),SIZE(VarIn,3),SIZE(VarIn,4)) :: Var_tmp
    
    CALL gather_omp(VarIn,Var_tmp)
    
    IF (is_omp_root) CALL gather_mpi(Var_tmp,VarOut)

  
  END SUBROUTINE gather_i3

  SUBROUTINE gather_i4(VarIn, VarOut)
  IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:,:,:) :: VarOut
    
    INTEGER, DIMENSION(nbp_mpi,SIZE(VarIn,2),SIZE(VarIn,3),SIZE(VarIn,4),SIZE(VarIn,5)) :: Var_tmp
    
    CALL gather_omp(VarIn,Var_tmp)
    
    IF (is_omp_root) CALL gather_mpi(Var_tmp,VarOut)

  
  END SUBROUTINE gather_i4

!!!!! --> Reals


  SUBROUTINE gather_r(VarIn, VarOut)
  IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:) :: VarOut
    
    REAL, DIMENSION(nbp_mpi) :: Var_tmp
    
    CALL gather_omp(VarIn,Var_tmp)
    
    IF (is_omp_root) CALL gather_mpi(Var_tmp,VarOut)
 
  END SUBROUTINE gather_r


  SUBROUTINE gather_r1(VarIn, VarOut)
  IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    REAL, DIMENSION(nbp_mpi,SIZE(VarIn,2)) :: Var_tmp
    
    CALL gather_omp(VarIn,Var_tmp)
    
    IF (is_omp_root) CALL gather_mpi(Var_tmp,VarOut)

  
  END SUBROUTINE gather_r1


  SUBROUTINE gather_r2(VarIn, VarOut)
  IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    REAL, DIMENSION(nbp_mpi,SIZE(VarIn,2),SIZE(VarIn,3)) :: Var_tmp
    
    CALL gather_omp(VarIn,Var_tmp)
    
    IF (is_omp_root) CALL gather_mpi(Var_tmp,VarOut)
  
  END SUBROUTINE gather_r2


  SUBROUTINE gather_r3(VarIn, VarOut)
  IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    REAL, DIMENSION(nbp_mpi,SIZE(VarIn,2),SIZE(VarIn,3),SIZE(VarIn,4)) :: Var_tmp
    
    CALL gather_omp(VarIn,Var_tmp)
    
    IF (is_omp_root) CALL gather_mpi(Var_tmp,VarOut)
  
  END SUBROUTINE gather_r3

  SUBROUTINE gather_r4(VarIn, VarOut)
  IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:,:,:) :: VarOut
    
    REAL, DIMENSION(nbp_mpi,SIZE(VarIn,2),SIZE(VarIn,3),SIZE(VarIn,4),SIZE(VarIn,5)) :: Var_tmp
    
    CALL gather_omp(VarIn,Var_tmp)
    
    IF (is_omp_root) CALL gather_mpi(Var_tmp,VarOut)
  
  END SUBROUTINE gather_r4

!!!!! --> Logiclas

  SUBROUTINE gather_l(VarIn, VarOut)
  IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:) :: VarOut
    
    LOGICAL, DIMENSION(nbp_mpi) :: Var_tmp
    
    CALL gather_omp(VarIn,Var_tmp)
    
    IF (is_omp_root) CALL gather_mpi(Var_tmp,VarOut)
  
  END SUBROUTINE gather_l


  SUBROUTINE gather_l1(VarIn, VarOut)
  IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    LOGICAL, DIMENSION(nbp_mpi,SIZE(VarIn,2)) :: Var_tmp
    
    CALL gather_omp(VarIn,Var_tmp)
    
    IF (is_omp_root) CALL gather_mpi(Var_tmp,VarOut)
  
  END SUBROUTINE gather_l1


  SUBROUTINE gather_l2(VarIn, VarOut)
  IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    LOGICAL, DIMENSION(nbp_mpi,SIZE(VarIn,2),SIZE(VarIn,3)) :: Var_tmp
    
    CALL gather_omp(VarIn,Var_tmp)
    IF (is_omp_root) CALL gather_mpi(Var_tmp,VarOut)
  
  END SUBROUTINE gather_l2


  SUBROUTINE gather_l3(VarIn, VarOut)
  IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    LOGICAL, DIMENSION(nbp_mpi,SIZE(VarIn,2),SIZE(VarIn,3),SIZE(VarIn,4)) :: Var_tmp
    
    CALL gather_omp(VarIn,Var_tmp)
    
    IF (is_omp_root) CALL gather_mpi(Var_tmp,VarOut)
  
  END SUBROUTINE gather_l3

  SUBROUTINE gather_l4(VarIn, VarOut)
  IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:,:) :: VarOut
    
    LOGICAL, DIMENSION(nbp_mpi,SIZE(VarIn,2),SIZE(VarIn,3),SIZE(VarIn,4),SIZE(VarIn,5)) :: Var_tmp
    
    CALL gather_omp(VarIn,Var_tmp)
    
    IF (is_omp_root) CALL gather_mpi(Var_tmp,VarOut)
  
  END SUBROUTINE gather_l4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Definition of reduce_sum for integers and reals   1D --> 4D   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Reduced sum subroutines for integers

  SUBROUTINE reduce_sum_i(VarIn, VarOut)
  IMPLICIT NONE
  
    INTEGER,INTENT(IN)  :: VarIn
    INTEGER,INTENT(OUT) :: VarOut
    
    INTEGER             :: Var_tmp
           
    CALL reduce_sum_omp(VarIn,Var_tmp)
    
    IF (is_omp_root) CALL reduce_sum_mpi(Var_tmp,VarOut)
  
  END SUBROUTINE reduce_sum_i  


  SUBROUTINE reduce_sum_i1(VarIn, VarOut)
  IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:)  :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:) :: VarOut
    
    INTEGER,DIMENSION(SIZE(VarIn))   :: Var_tmp
           
    CALL reduce_sum_omp(VarIn,Var_tmp)
    
    IF (is_omp_root) CALL reduce_sum_mpi(Var_tmp,VarOut)
  
  END SUBROUTINE reduce_sum_i1  


  SUBROUTINE reduce_sum_i2(VarIn, VarOut)
  IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:)  :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    INTEGER,DIMENSION(SIZE(VarIn,1),SIZE(VarIn,2)) :: Var_tmp
           
    CALL reduce_sum_omp(VarIn,Var_tmp)
    
    IF (is_omp_root) CALL reduce_sum_mpi(Var_tmp,VarOut)

  END SUBROUTINE reduce_sum_i2  
  

  SUBROUTINE reduce_sum_i3(VarIn, VarOut)
  IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:)  :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    INTEGER,DIMENSION(SIZE(VarIn,1),SIZE(VarIn,2),SIZE(VarIn,3))  :: Var_tmp
           
    CALL reduce_sum_omp(VarIn,Var_tmp)
    
    IF (is_omp_root) CALL reduce_sum_mpi(Var_tmp,VarOut)
  
  END SUBROUTINE reduce_sum_i3  


  SUBROUTINE reduce_sum_i4(VarIn, VarOut)
  IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:,:)  :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    INTEGER,DIMENSION(SIZE(VarIn,1),SIZE(VarIn,2),SIZE(VarIn,3),SIZE(VarIn,4))  :: Var_tmp
           
    CALL reduce_sum_omp(VarIn,Var_tmp)
    
    IF (is_omp_root) CALL reduce_sum_mpi(Var_tmp,VarOut)
  
  END SUBROUTINE reduce_sum_i4  


! Reduce sum subroutines for reals

  SUBROUTINE reduce_sum_r(VarIn, VarOut)
  IMPLICIT NONE
  
    REAL,INTENT(IN)  :: VarIn
    REAL,INTENT(OUT) :: VarOut
    
    REAL             :: Var_tmp
           
    CALL reduce_sum_omp(VarIn,Var_tmp)
    
    IF (is_omp_root) CALL reduce_sum_mpi(Var_tmp,VarOut)
  
  END SUBROUTINE reduce_sum_r  


  SUBROUTINE reduce_sum_r1(VarIn, VarOut)
  IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:)  :: VarIn
    REAL,INTENT(OUT),DIMENSION(:) :: VarOut
    
    REAL,DIMENSION(SIZE(VarIn))   :: Var_tmp
           
    CALL reduce_sum_omp(VarIn,Var_tmp)
    
    IF (is_omp_root) CALL reduce_sum_mpi(Var_tmp,VarOut)
  
  END SUBROUTINE reduce_sum_r1  


  SUBROUTINE reduce_sum_r2(VarIn, VarOut)
  IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:)  :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    REAL,DIMENSION(SIZE(VarIn,1),SIZE(VarIn,2)) :: Var_tmp
           
    CALL reduce_sum_omp(VarIn,Var_tmp)
    
    IF (is_omp_root) CALL reduce_sum_mpi(Var_tmp,VarOut)
  
  END SUBROUTINE reduce_sum_r2  
  

  SUBROUTINE reduce_sum_r3(VarIn, VarOut)
  IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:)  :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    REAL,DIMENSION(SIZE(VarIn,1),SIZE(VarIn,2),SIZE(VarIn,3))  :: Var_tmp
           
    CALL reduce_sum_omp(VarIn,Var_tmp)
    
    IF (is_omp_root) CALL reduce_sum_mpi(Var_tmp,VarOut)
  
  END SUBROUTINE reduce_sum_r3  


  SUBROUTINE reduce_sum_r4(VarIn, VarOut)
  IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:,:)  :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    REAL,DIMENSION(SIZE(VarIn,1),SIZE(VarIn,2),SIZE(VarIn,3),SIZE(VarIn,4))  :: Var_tmp
           
    CALL reduce_sum_omp(VarIn,Var_tmp)
    
    IF (is_omp_root) CALL reduce_sum_mpi(Var_tmp,VarOut)
  
  END SUBROUTINE reduce_sum_r4  

! AllReduce sum subroutines for reals

  SUBROUTINE allreduce_sum_r(VarIn, VarOut)
  IMPLICIT NONE
  
    REAL,INTENT(IN)  :: VarIn
    REAL,INTENT(OUT) :: VarOut
    
    REAL             :: Var_tmp
           
    CALL reduce_sum_omp(VarIn,Var_tmp)
    
    IF (is_omp_root) CALL allreduce_sum_mpi(Var_tmp,VarOut)
  
  END SUBROUTINE allreduce_sum_r  

   
END MODULE mod_orchidee_transfert_para

