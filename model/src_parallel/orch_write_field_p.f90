! ==============================================================================================================================
! MODULE   : orch_Write_Field_p
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF     Set of interfaces to create netcdf output files to test values of variables in parallel mode.
!!
!! \n DESCRIPTION  : Set of interfaces to create netcdf output files to test values of variables in parallel mode. 
!!                   IMPORTANT NOTICE: These subroutines have not been tested in hybrid run mode using MPI-OMP paralelization
!!                                     using more than 1 threds OMP. It will probably not work!
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
MODULE Orch_Write_field_p
  
  !! ==============================================================================================================================
  !! INTERFACE   :  WriteField_p
  !!
  !>\BRIEF         set of routines to write real fields (of 1d, 2d, 3d, 4d) in netcdf output file 
  !!
  !! DESCRIPTION  : set of routines to write real fields (of 1d, 2d, 3d, 4d) in netcdf output file 
  !!                CALL WriteField_p("MyVariable", variable_array) 
  !!                will create a file MyVariable.nc with all value of variable_array
  !!
  !! \n
  !_ ================================================================================================================================
  INTERFACE WriteField_p
    MODULE PROCEDURE WriteField_4d_p,WriteField_3d_p,WriteField_2d_p
  END INTERFACE WriteField_p
 
  !! ==============================================================================================================================
  !! INTERFACE   :  WriteFieldI_p
  !!
  !>\BRIEF         set of routines to write integer fields (of 1d, 2d, 3d, 4d) in netcdf output file 
  !!
  !! DESCRIPTION  : set of routines to write integer fields (of 1d, 2d, 3d, 4d) in netcdf output file 
  !!                CALL WriteFieldI_p("MyVariable", variable_array) 
  !!                will create a file MyVariable.nc with all value of variable_array
  !!
  !! \n
  !_ ================================================================================================================================
  INTERFACE WriteFieldI_p
    MODULE PROCEDURE WriteFieldI_3d_p,WriteFieldI_2d_p,WriteFieldI_1d_p
  END INTERFACE WriteFieldI_p  
  
  
CONTAINS

  SUBROUTINE init_WriteField_p
  USE mod_orchidee_para
  USE Orch_Write_Field, ONLY : Init_WriteField
  IMPLICIT NONE
    IF (is_root_prc) CALL Init_WriteField(iim_g,jjm_g,nbp_glo,index_g)
    
  END SUBROUTINE init_WriteField_p

  SUBROUTINE WriteField_4d_p(name,Field)
    USE mod_orchidee_para
    USE Orch_Write_Field, ONLY : WriteField
    IMPLICIT NONE
      CHARACTER(len=*) :: name
      REAL, DIMENSION(:,:,:,:) :: Field 
      INTEGER, DIMENSION(4) :: Dim
      
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: Field_g
      
      Dim=SHAPE(Field)
      
      ALLOCATE(Field_g(iim_g,jjm_g,DIM(3),DIM(4)))
      CALL Gather2D_mpi(Field,Field_g)

      IF (is_root_prc) CALL WriteField(name,Field_g)  
      
      DEALLOCATE(Field_g)
  END SUBROUTINE WriteField_4d_p
    
  SUBROUTINE WriteField_3d_p(name,Field)
    USE mod_orchidee_para
    USE Orch_Write_Field, ONLY : WriteField
    IMPLICIT NONE
      CHARACTER(len=*) :: name
      REAL, DIMENSION(:,:,:) :: Field 
      INTEGER, DIMENSION(3) :: Dim
      
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: Field_g
      
      
      Dim=SHAPE(Field)
      
      ALLOCATE(Field_g(iim_g,jjm_g,DIM(3)))
      CALL Gather2D_mpi(Field,Field_g)

      IF (is_root_prc) CALL WriteField(name,Field_g)  
      
      DEALLOCATE(Field_g)
  END SUBROUTINE WriteField_3d_p

  SUBROUTINE WriteField_2d_p(name,Field)
    USE mod_orchidee_para
    USE Orch_Write_Field, ONLY : WriteField
    IMPLICIT NONE
      CHARACTER(len=*) :: name
      REAL, DIMENSION(:,:) :: Field 
      INTEGER, DIMENSION(2) :: Dim
      
      REAL, ALLOCATABLE, DIMENSION(:,:) :: Field_g
      
      
      Dim=SHAPE(Field)
      
      ALLOCATE(Field_g(iim_g,jjm_g))
      CALL Gather2D_mpi(Field,Field_g)

      IF (is_root_prc) CALL WriteField_gen(name,Field_g,2,Dim)  
      
      DEALLOCATE(Field_g)
  END SUBROUTINE WriteField_2d_p

  SUBROUTINE WriteFieldI_3d_p(name,Field)
    USE mod_orchidee_para
    USE Orch_Write_Field, ONLY : WriteFieldI
    IMPLICIT NONE
      CHARACTER(len=*) :: name
      REAL, DIMENSION(:,:,:) :: Field 
      INTEGER, DIMENSION(3) :: Dim
      
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: Field_g
      
      
      Dim=SHAPE(Field)
      
      ALLOCATE(Field_g(nbp_glo,DIM(2),DIM(3)))
      CALL gather(Field,Field_g)
      
      IF (is_root_prc) CALL WriteFieldI(name,Field_g)  
      
      DEALLOCATE(Field_g)
  END SUBROUTINE WriteFieldI_3d_p

  SUBROUTINE WriteFieldI_2d_p(name,Field)
    USE mod_orchidee_para
    USE Orch_Write_Field, ONLY : WriteFieldI
    IMPLICIT NONE
      CHARACTER(len=*) :: name
      REAL, DIMENSION(:,:) :: Field 
      INTEGER, DIMENSION(2) :: Dim
      
      REAL, ALLOCATABLE, DIMENSION(:,:) :: Field_g
      
      
      Dim=SHAPE(Field)
      
      ALLOCATE(Field_g(nbp_glo,DIM(2)))
      CALL gather(Field,Field_g)
      
      IF (is_root_prc) CALL WriteFieldI(name,Field_g)  
      
      DEALLOCATE(Field_g)
  END SUBROUTINE WriteFieldI_2d_p    

  SUBROUTINE WriteFieldI_1d_p(name,Field)
    USE mod_orchidee_para
    USE Orch_Write_Field, ONLY : WriteFieldI
    IMPLICIT NONE
      CHARACTER(len=*) :: name
      REAL, DIMENSION(:) :: Field 
      INTEGER, DIMENSION(1) :: Dim
      
      REAL, ALLOCATABLE, DIMENSION(:) :: Field_g
      
      
      Dim=SHAPE(Field)
      
      ALLOCATE(Field_g(nbp_glo))
      CALL gather(Field,Field_g)
      
      IF (is_root_prc) CALL WriteFieldI(name,Field_g)  
      
      DEALLOCATE(Field_g)
  END SUBROUTINE WriteFieldI_1d_p    
    
END MODULE Orch_Write_field_p
