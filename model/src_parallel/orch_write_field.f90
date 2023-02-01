! ==============================================================================================================================
! MODULE   : orch_Write_Field
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF     Set of low level interfaces to create netcdf output files to test values of variables in sequentiel mode  
!!
!!\n DESCRIPTION  : Set of low level interfaces to create netcdf output files to test values of variables in sequentiel mode.
!!                  The interfaces in this module are only used by orch_write_field_p to create high level interfaces.
!!                  These interfaces should only be called by the master process.
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
MODULE orch_Write_Field
  
  USE mod_orchidee_para

  IMPLICIT NONE

  INTEGER, PARAMETER :: MaxWriteField = 100
  INTEGER, DIMENSION(MaxWriteField),SAVE :: FieldId
  INTEGER, DIMENSION(MaxWriteField),SAVE :: FieldVarId
  INTEGER, DIMENSION(MaxWriteField),SAVE :: FieldIndex
  CHARACTER(len=255), DIMENSION(MaxWriteField) ::  FieldName
  
  INTEGER, SAVE,DIMENSION(:), ALLOCATABLE :: Index_Write_Field
  INTEGER,SAVE :: iim
  INTEGER,SAVE :: jjm
  INTEGER,SAVE :: NbPoint
  REAL, PARAMETER :: undef_var=0.
  
  INTEGER,SAVE :: NbField = 0
  
  !! ==============================================================================================================================
  !! INTERFACE   :  WriteField
  !!
  !>\BRIEF         set of routines to write real fields (of 1d, 2d, 3d, 4d) in netcdf output file 
  !!
  !! DESCRIPTION  : set of routines to write real fields (of 1d, 2d, 3d, 4d) in netcdf output file 
  !!                CALL WriteField("MyVariable", variable_array) 
  !!                will create a file MyVariable.nc with all value of variable_array
  !!
  !! \n
  !_ ================================================================================================================================
  INTERFACE WriteField
    MODULE PROCEDURE WriteField_4d,WriteField_3d,WriteField_2d,WriteField_1d
  END INTERFACE WriteField
 
  !! ==============================================================================================================================
  !! INTERFACE   :  WriteFieldI
  !!
  !>\BRIEF         set of routines to write integer fields (of 1d, 2d, 3d, 4d) in netcdf output file 
  !!
  !! DESCRIPTION  : set of routines to write integer fields (of 1d, 2d, 3d, 4d) in netcdf output file 
  !!                CALL WriteFieldI("MyVariable", variable_array) 
  !!                will create a file MyVariable.nc with all value of variable_array
  !!
  !! \n
  !_ ================================================================================================================================
  INTERFACE WriteFieldI
    MODULE PROCEDURE WriteFieldI_3d,WriteFieldI_2d,WriteFieldI_1d
  END INTERFACE WriteFieldI

  PRIVATE :: iim,jjm,NbPoint 
  CONTAINS
  
    SUBROUTINE Init_WriteField(iim0,jjm0,NbPoint0,Index0)
    IMPLICIT NONE
      INTEGER,INTENT(in) :: iim0
      INTEGER,INTENT(in) :: jjm0
      INTEGER,INTENT(in) :: NbPoint0
      INTEGER,INTENT(in) :: Index0(NbPoint0)
    
      iim=iim0
      jjm=jjm0
      Nbpoint=Nbpoint0
      ALLOCATE(Index_Write_Field(NbPoint))
      Index_Write_Field(:)=Index0(:)
    END SUBROUTINE Init_WriteField
    
    FUNCTION GetFieldIndex(name)
    IMPLICIT NONE
      INTEGER          :: GetFieldindex
      CHARACTER(len=*) :: name
    
      CHARACTER(len=255) :: TrueName
      INTEGER            :: i
       
      
      TrueName=TRIM(ADJUSTL(name))
    
      GetFieldIndex=-1
      DO i=1,NbField
        IF (TrueName==FieldName(i)) THEN
          GetFieldIndex=i
          EXIT
        ENDIF
      ENDDO
    END FUNCTION GetFieldIndex

    SUBROUTINE WriteFieldI_3d(name,Field)
    IMPLICIT NONE
      CHARACTER(len=*) :: name
      REAL, DIMENSION(:,:,:) :: Field 
      INTEGER, DIMENSION(3) :: Dim
      INTEGER,DIMENSION(4) :: Dim_tmp
      INTEGER :: i
      
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: Field_tmp 
      
      Dim=SHAPE(Field)
      ALLOCATE(Field_tmp(iim*jjm,DIM(2),DIM(3)))
      field_tmp(:,:,:)=undef_var
      
      DO i=1,NbPoint
        field_tmp(Index_Write_Field(i),:,:)=Field(i,:,:)
      ENDDO
      
      Dim_tmp(1)=iim
      Dim_tmp(2)=jjm
      Dim_tmp(3)=DIM(2)
      Dim_tmp(4)=DIM(3)
      CALL WriteField_gen(name,Field_tmp,4,Dim_tmp)  
  
      DEALLOCATE(Field_tmp)
    END SUBROUTINE WriteFieldI_3d

    SUBROUTINE WriteFieldI_2d(name,Field)
    IMPLICIT NONE
      CHARACTER(len=*) :: name
      REAL, DIMENSION(:,:) :: Field 
      INTEGER, DIMENSION(2) :: Dim
      INTEGER,DIMENSION(3) :: Dim_tmp
      INTEGER :: i
      
      REAL, ALLOCATABLE, DIMENSION(:,:) :: Field_tmp 
      
      Dim=SHAPE(Field)
      ALLOCATE(Field_tmp(iim*jjm,DIM(2)))
      field_tmp(:,:)=undef_var
      
      DO i=1,NbPoint
        field_tmp(Index_Write_Field(i),:)=Field(i,:)
      ENDDO
      
      Dim_tmp(1)=iim
      Dim_tmp(2)=jjm
      Dim_tmp(3)=DIM(2)

      CALL WriteField_gen(name,Field_tmp,3,Dim_tmp)  
  
      DEALLOCATE(Field_tmp)
    END SUBROUTINE WriteFieldI_2d

    SUBROUTINE WriteFieldI_1d(name,Field)
    IMPLICIT NONE
      CHARACTER(len=*) :: name
      REAL, DIMENSION(:) :: Field 
      INTEGER, DIMENSION(1) :: Dim
      INTEGER,DIMENSION(2) :: Dim_tmp
      INTEGER :: i
      
      REAL, ALLOCATABLE, DIMENSION(:) :: Field_tmp 
      
      Dim=SHAPE(Field)
      ALLOCATE(Field_tmp(iim*jjm))
      field_tmp(:)=undef_var
      
      DO i=1,NbPoint
        field_tmp(Index_Write_Field(i))=Field(i)
      ENDDO
      
      Dim_tmp(1)=iim
      Dim_tmp(2)=jjm

      CALL WriteField_gen(name,Field_tmp,2,Dim_tmp)  
  
      DEALLOCATE(Field_tmp)
    END SUBROUTINE WriteFieldI_1d
        
    SUBROUTINE WriteField_4d(name,Field)
    IMPLICIT NONE
      CHARACTER(len=*) :: name
      REAL, DIMENSION(:,:,:,:) :: Field 
      INTEGER, DIMENSION(4) :: Dim
      
      Dim=SHAPE(Field)
      CALL WriteField_gen(name,Field,4,Dim)  
  
    END SUBROUTINE WriteField_4d
     
    SUBROUTINE WriteField_3d(name,Field)
    IMPLICIT NONE
      CHARACTER(len=*) :: name
      REAL, DIMENSION(:,:,:) :: Field 
      INTEGER, DIMENSION(3) :: Dim
      
      Dim=SHAPE(Field)
      CALL WriteField_gen(name,Field,3,Dim)  
  
    END SUBROUTINE WriteField_3d
    
    SUBROUTINE WriteField_2d(name,Field)
    IMPLICIT NONE
      CHARACTER(len=*) :: name
      REAL, DIMENSION(:,:) :: Field 
      INTEGER, DIMENSION(2) :: Dim
      
      Dim=SHAPE(Field)
      CALL WriteField_gen(name,Field,2,Dim)  
  
    END SUBROUTINE WriteField_2d
    
    SUBROUTINE WriteField_1d(name,Field)
    IMPLICIT NONE
      CHARACTER(len=*) :: name
      REAL, DIMENSION(:) :: Field 
      INTEGER, DIMENSION(1) :: Dim
      
      Dim=SHAPE(Field)
      CALL WriteField_gen(name,Field,1,Dim)  
  
    END SUBROUTINE WriteField_1d
        
    SUBROUTINE CreateNewField(name,NbDim,DimSize)
    USE ioipsl
    IMPLICIT NONE
    INCLUDE 'netcdf.inc'  
      CHARACTER(len=*) :: name
      INTEGER :: NbDim
      INTEGER :: DimSize(NbDim)
      INTEGER :: TabDim(NbDim+1)
      INTEGER :: status
      
      
      NbField=NbField+1
      FieldName(NbField)=TRIM(ADJUSTL(name))
      FieldIndex(NbField)=1
      
      WRITE(numout,*) 'CREATE_NEW_FIELD ',name,NbDim,DimSize
      CALL FLUSH(6)
      status = NF_CREATE(TRIM(ADJUSTL(name))//'.nc', NF_CLOBBER, FieldId(NbField))
      IF (NbDim>=1) status = NF_DEF_DIM(FieldId(NbField),'I',DimSize(1),TabDim(1))
      IF (NbDim>=2) status = NF_DEF_DIM(FieldId(NbField),'J',DimSize(2),TabDim(2))
      IF (NbDim>=3) status = NF_DEF_DIM(FieldId(NbField),'K',DimSize(3),TabDim(3))
      IF (NbDim>=4) status = NF_DEF_DIM(FieldId(NbField),'L',DimSize(4),TabDim(4))
      status = NF_DEF_DIM(FieldId(NbField),'iter',NF_UNLIMITED,TabDim(NbDim+1))
      status = NF_DEF_VAR(FieldId(NbField),FieldName(NbField),NF_DOUBLE,NbDim+1,TabDim,FieldVarId(NbField))
      status = NF_ENDDEF(FieldId(NbField))

    END SUBROUTINE CreateNewField
    
  FUNCTION int2str(int)
    IMPLICIT NONE
    INTEGER, PARAMETER :: MaxLen=10
    INTEGER,INTENT(in) :: int
    CHARACTER(len=MaxLen) :: int2str
    LOGICAL :: flag
    INTEGER :: i
    flag=.TRUE.
    
    i=int
    
    int2str=''
    DO WHILE (flag)
      int2str=CHAR(MOD(i,10)+48)//int2str
      i=i/10
      IF (i==0) flag=.FALSE.
    ENDDO
  END FUNCTION int2str


END MODULE Orch_Write_Field

    SUBROUTINE WriteField_gen(name,Field,NbDim,DimSize)
    USE orch_write_field
    IMPLICIT NONE
    INCLUDE 'netcdf.inc'
      CHARACTER(len=*) :: name
      INTEGER :: NbDim
      INTEGER,DIMENSION(NbDim) :: DimSize
      REAL,DIMENSION(*) :: Field
      
      INTEGER :: status
      INTEGER :: ind
      INTEGER :: start(NbDim+1)
      INTEGER :: COUNT(NbDim+1)
      INTEGER :: i
           
      Ind=GetFieldIndex(name)
      IF (Ind==-1) THEN
        CALL CreateNewField(name,NbDim,DimSize)
	Ind=GetFieldIndex(name)
      ELSE
        FieldIndex(Ind)=FieldIndex(Ind)+1
      ENDIF
      
      DO i=1,NbDim
        start(i)=1
        COUNT(i)=DimSize(i)
      ENDDO
      start(NbDim+1)=FieldIndex(Ind)
      COUNT(NbDim+1)=1

      status = NF_PUT_VARA_DOUBLE(FieldId(Ind),FieldVarId(Ind),start,count,Field)
      status = NF_SYNC(FieldId(Ind))
      
    END SUBROUTINE WriteField_gen
