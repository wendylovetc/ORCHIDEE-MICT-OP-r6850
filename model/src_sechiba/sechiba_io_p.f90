! ================================================================================================================================
!  MODULE       : sechiba_io_p
!
!  CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
!  LICENCE      : IPSL (2006)
!  This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF   To be used for intializing variables not read available in the restart file
!!
!!\n DESCRIPTION: This module contains the interface to setvar_p to be used for intializing variables if they were 
!!                not found in the restart file. The variable will only be modified if it was not found in the restart
!!                file (i.e. if it is eqaul val_exp).
!!
!!                Syntax : CALL setvar_p (var, val_exp, key_wd, val_put)
!!                  var : the variable to initialize; It can be an integer or a real, a scalar or have 1 or 2 dimensions
!!                  val_exp : the value set by restget_p if the variable was not found in the restart file (do not change this)
!!                  key_wd  : parameter name to be searched for in run.def
!!                  val_put : a value to be used if the kew_wd was not found in run.def. val_put must have the same or
!!                            smaller rank as var
!!
!!                Note that setvar_p must always be called by all processes because it contains call to getin_p.
!!                - The variable var, will only be modified if before the call it is equal to val_exp. Otherwise nothing is done.
!!                - If var is equal to val_exp and if key_wd is not equal "NO_KEYWORD" or "NOKEYWORD", then the value for key_wd 
!!                  is read from run.def using getin_p and used to initialize var.
!!                - If key_wd is not found in run.def or if key_wd="NO_KEYWORD" or "NOKEYWORD", then the val_put will be used to 
!!                  initialize var.
!!
!!                The interface will automatically call the right subroutine depending on the type of input variables.
!!
!! REFERENCE(S)	: None
!!
!! SVN          :
!! $HeadURL$
!! $Date$
!! $Revision$
!! \n
!_ ================================================================================================================================

MODULE sechiba_io_p

  USE defprec
!  USE constantes
!  USE ioipsl
  USE ioipsl_para
  USE mod_orchidee_para

  IMPLICIT NONE

  PRIVATE
  PUBLIC setvar_p

  INTERFACE setvar_p
    MODULE PROCEDURE i0setvar_p, i10setvar_p, i20setvar_p, i11setvar_p, i21setvar_p
    MODULE PROCEDURE r0setvar_p, r10setvar_p, r20setvar_p, r11setvar_p, r21setvar_p, r22setvar_p, r30setvar_p
  END INTERFACE

  LOGICAL, SAVE                  :: long_print_setvar_p=.FALSE.  !! change to true to have more information
!$OMP THREADPRIVATE(long_print_setvar_p)

CONTAINS 

!!  =============================================================================================================================
!! SUBROUTINE:    i0setvar_p
!!
!>\BRIEF	  Subroutine for initializing an integer scalar variable with a scalar integer.
!!
!! DESCRIPTION:	  Subroutine for initializing an integer scalar variable with a scalar integer.
!!                This subroutine must be called by all processes.
!! \n
!_ ==============================================================================================================================
SUBROUTINE i0setvar_p (var, val_exp, key_wd, val_put)

  INTEGER(i_std), INTENT(inout)                   :: var                  !! Integer scalar to modify
  INTEGER(i_std), INTENT(in)                      :: val_exp              !! Exceptional value
  CHARACTER(LEN=*), INTENT(in)                    :: key_wd               !! The Key word we will look for
  INTEGER(i_std), INTENT(in)                      :: val_put              !! Initial value to stored

  INTEGER(i_std)                                  :: val_tmp
  INTEGER(i_std)                                  :: is_key

  is_key = MAX(INDEX(key_wd, 'NO_KEYWORD'), INDEX(key_wd, 'NOKEYWORD'))
  
  IF (long_print_setvar_p) WRITE (numout,*) "i0setvar :", key_wd, val_exp, val_put

  val_tmp = val_put

  IF ( var == val_exp ) THEN 
     IF ( is_key <= 0 ) THEN
        CALL getin_p(key_wd,  val_tmp)
     ENDIF
     var = val_tmp
  END IF
  
END SUBROUTINE i0setvar_p


!!  =============================================================================================================================
!! SUBROUTINE:    i10setvar_p
!!
!>\BRIEF	  Subroutine for initializing an integer 1D array with a integer scalar variable.
!!
!! DESCRIPTION:	  Subroutine for initializing an integer 1D array with a integer scalar variable.
!!                This subroutine must be called by all processes.
!! \n
!_ ==============================================================================================================================
SUBROUTINE i10setvar_p (var, val_exp, key_wd, val_put)

  INTEGER(i_std), DIMENSION(:), INTENT(inout)     :: var                  !! 1D integer array to modify
  INTEGER(i_std), INTENT(in)                      :: val_exp              !! Exceptional value
  CHARACTER(LEN=*), INTENT(in)                :: key_wd               !! The Key word we will look for
  INTEGER(i_std), INTENT(in)                      :: val_put              !! Scalar value to stored
  
  INTEGER(i_std)                                  :: val_tmp
  INTEGER(i_std)                                  :: is_key

  is_key = MAX(INDEX(key_wd, 'NO_KEYWORD'), INDEX(key_wd, 'NOKEYWORD'))

  IF (long_print_setvar_p) WRITE (numout,*) "i10setvar :", key_wd, val_exp, val_put

  val_tmp = val_put

  IF ( ALL( var(:) == val_exp ) ) THEN
     IF ( is_key <= 0 ) THEN
       CALL getin_p(key_wd,  val_tmp)
     ENDIF
     var(:) = val_tmp
  END IF
  
END SUBROUTINE i10setvar_p


!!  =============================================================================================================================
!! SUBROUTINE:    i11setvar_p
!!
!>\BRIEF	  Subroutine for initializing an integer 1D array with another integer 1D array.
!!
!! DESCRIPTION:	  Subroutine for initializing an integer 1D array with another integer 1D array.
!!                This subroutine must be called by all processes.
!! \n
!_ ==============================================================================================================================
SUBROUTINE i11setvar_p (var, val_exp, key_wd, val_put, is_grid)
  
  INTEGER(i_std), DIMENSION(:), INTENT(inout)     :: var                 !! 1D integer array to modify
  INTEGER(i_std), INTENT(in)                      :: val_exp             !! Exceptional value
  CHARACTER(LEN=*), INTENT(in)                    :: key_wd              !! The Key word we will look for
  INTEGER(i_std), DIMENSION(:), INTENT(in)        :: val_put             !! 1D integer array to stored
  LOGICAL,        OPTIONAL                        :: is_grid             !! Parameter present indicates a setvar for a grid variable 

  INTEGER(i_std), ALLOCATABLE,DIMENSION(:)        :: val_tmp
  INTEGER(i_std), ALLOCATABLE,DIMENSION(:)        :: val_tmp_g
  INTEGER(i_std)                                  :: is_key

  is_key = MAX(INDEX(key_wd, 'NO_KEYWORD'), INDEX(key_wd, 'NOKEYWORD'))
  
  IF (long_print_setvar_p) WRITE (numout,*) "i11setvar :", key_wd, val_exp, SIZE(val_put), val_put(1)

  ALLOCATE(val_tmp(SIZE(val_put)))
  val_tmp(:) = val_put(:)

  IF ( ALL( var(:) == val_exp ) ) THEN 
     IF ( is_key <= 0 ) THEN
        IF (PRESENT(is_grid) ) THEN
           IF (is_root_prc) THEN
              ALLOCATE( val_tmp_g(nbp_glo) )
           ELSE
              ALLOCATE( val_tmp_g(1) )
           ENDIF
           CALL gather( val_tmp,val_tmp_g )
           IF (is_root_prc) &
              CALL getin(key_wd,  val_tmp_g)
           CALL scatter( val_tmp,val_tmp_g )
           DEALLOCATE( val_tmp_g )
        ELSE
           CALL getin_p(key_wd,  val_tmp)
        ENDIF
     ENDIF
     var(:) = val_tmp (:)
  END IF

  DEALLOCATE(val_tmp)
  
END SUBROUTINE i11setvar_p


!!  =============================================================================================================================
!! SUBROUTINE:    i20setvar_p
!!
!>\BRIEF	  Subroutine for initializing an integer 2D variable with a scalar integer variable.
!!
!! DESCRIPTION:	  Subroutine for initializing an integer 2D variable with a scalar integer variable.
!!                This subroutine must be called by all processes.
!! \n
!_ ==============================================================================================================================
SUBROUTINE i20setvar_p (var, val_exp, key_wd, val_put)
  
  INTEGER(i_std), DIMENSION(:,:), INTENT(inout)   :: var                  !! 2D integer array to modify
  INTEGER(i_std), INTENT(in)                      :: val_exp              !! Exceptional value
  CHARACTER(LEN=*), INTENT(in)                    :: key_wd               !! The Key word we will look for
  INTEGER(i_std), INTENT(in)                      :: val_put              !! Scalar value to be used as default

  INTEGER(i_std)                                  :: val_tmp
  INTEGER(i_std)                                  :: is_key

  is_key = MAX(INDEX(key_wd, 'NO_KEYWORD'), INDEX(key_wd, 'NOKEYWORD'))
  
  IF (long_print_setvar_p) WRITE (numout,*) "i20setvar :", key_wd, val_exp, val_put

  val_tmp = val_put

  IF ( ALL( var(:,:) == val_exp ) ) THEN 
     IF ( is_key <= 0 ) THEN
       CALL getin_p(key_wd,  val_tmp)
     ENDIF
     var(:,:) = val_tmp
  END IF
  
END SUBROUTINE i20setvar_p


!!  =============================================================================================================================
!! SUBROUTINE:    i21setvar_p
!!
!>\BRIEF	  Subroutine for initialieing an 2D integer variable with a 1D array integer.
!!
!! DESCRIPTION:	  Subroutine for initialieing an 2D integer variable with a 1D array integer.
!!                This subroutine must be called by all processes.
!!                Row or column depending size of 1D array to stored.
!!
!!                example: 1D 1,2,3     2D is 1, 2, 3,
!!                                            1, 2, 3
!!
!!                example: 1D 1,2,3     2D is 1, 1,
!!                                            2, 2,
!!                                            3, 3
!! \n
!_ ==============================================================================================================================
SUBROUTINE i21setvar_p (var, val_exp, key_wd, val_put, is_grid)
  
  INTEGER(i_std), DIMENSION(:,:), INTENT(inout)   :: var                  !! 2D integer array to modify
  INTEGER(i_std), INTENT(in)                      :: val_exp              !! Exceptional value
  CHARACTER(LEN=*), INTENT(in)                    :: key_wd               !! The Key word we will look for
  INTEGER(i_std), DIMENSION(:), INTENT(in)        :: val_put              !! 1D integer array to stored
  LOGICAL,        OPTIONAL                        :: is_grid              !! Parameter present indicates a setvar for a grid variable 
  
  INTEGER(i_std), ALLOCATABLE,DIMENSION(:)        :: val_tmp
  INTEGER(i_std)                                  :: is_key

  is_key = MAX(INDEX(key_wd, 'NO_KEYWORD'), INDEX(key_wd, 'NOKEYWORD'))

  ! test if the 1D array dimension is compatible with first or second 
  ! dimension of the 2D array

  IF (long_print_setvar_p) WRITE (numout,*) "i21setvar :", key_wd, val_exp, val_put

  ALLOCATE(val_tmp(SIZE(val_put)))
  val_tmp(:) = val_put(:)

  IF (SIZE(val_put)==SIZE(var,1)) THEN 
      !
      ! example: 1D 1.,2.,3.     2D is 1., 2., 3.,
      !                                1., 2., 3.
      !
      IF ( ALL( var(:,:) == val_exp ) ) THEN 
         IF ( is_key <= 0 ) THEN
           CALL getin_p(key_wd,  val_tmp)
         ENDIF
         var(:,:) = SPREAD(val_tmp(:),2,SIZE(var,1))
      END IF
  ELSEIF (SIZE(val_put)==SIZE(var,2)) THEN 
      !
      ! example: 1D 1.,2.,3.     2D is 1., 1.,
      !                                2., 2.,
      !                                3., 3.
      !
      IF ( ALL( var(:,:) == val_exp ) ) THEN 
         IF ( is_key <= 0 ) THEN
           CALL getin_p(key_wd,  val_tmp)
         ENDIF
         var(:,:) = SPREAD(val_tmp(:),1,SIZE(var,1))
      END IF
  ELSE 
      WRITE (numout,*) ' incompatible dimension var and val_put'
      WRITE (numout,*) ' var     ', SIZE(var,1), SIZE(var,2)
      WRITE (numout,*) ' val_put ', SIZE(val_put)
      STOP 'setvar'
  END IF

  DEALLOCATE(val_tmp)
  
END SUBROUTINE i21setvar_p


!!  =============================================================================================================================
!! SUBROUTINE:    r0setvar_p
!!
!>\BRIEF	  Subroutine for initializing a real scalar variable.
!!
!! DESCRIPTION:	  Subroutine for initializing a real scalar variable with a real scalar variable. 
!!                This subroutine must be called by all processes.
!! \n
!_ ==============================================================================================================================
SUBROUTINE r0setvar_p (var, val_exp, key_wd, val_put)
  
  REAL(r_std), INTENT(inout)                   :: var                  !! Real scalar to modify
  REAL(r_std), INTENT(in)                      :: val_exp              !! Exceptional value
  CHARACTER(LEN=*), INTENT(in)                   :: key_wd               !! The Key word we will look for
  REAL(r_std), INTENT(in)                      :: val_put              !! Initial value to stored
  
  REAL(r_std)                                  :: val_tmp
  INTEGER(i_std)                                     :: is_key

  is_key = MAX(INDEX(key_wd, 'NO_KEYWORD'), INDEX(key_wd, 'NOKEYWORD'))

  IF (long_print_setvar_p) WRITE (numout,*) "r0setvar :", key_wd, val_exp, val_put

  val_tmp = val_put

  IF ( var==val_exp ) THEN 
     IF ( is_key <= 0 ) THEN
       CALL getin_p(key_wd,  val_tmp)
     ENDIF
     var = val_tmp
  END IF
  
END SUBROUTINE r0setvar_p


!!  =============================================================================================================================
!! SUBROUTINE:    r10setvar_p
!!
!>\BRIEF	  Subroutine for initializing an real 1D array with a real scalar variable.
!!
!! DESCRIPTION:	  Subroutine for initializing an real 1D array with a real scalar variable.
!!                This subroutine must be called by all processes.
!! \n
!_ ==============================================================================================================================
SUBROUTINE r10setvar_p (var, val_exp, key_wd, val_put)
  
  REAL(r_std), DIMENSION(:), INTENT(inout)     :: var                  !! 1D real array to modify
  REAL(r_std), INTENT(in)                      :: val_exp              !! Exceptional value
  CHARACTER(LEN=*), INTENT(in)                 :: key_wd               !! The Key word we will look for
  REAL(r_std), INTENT(in)                      :: val_put              !! Scalar value to stored
   
  REAL(r_std)                                  :: val_tmp
  INTEGER(i_std)                               :: is_key

  is_key = MAX(INDEX(key_wd, 'NO_KEYWORD'), INDEX(key_wd, 'NOKEYWORD'))
 
  IF (long_print_setvar_p) WRITE (numout,*) "r10setvar :", key_wd, val_exp, val_put

  val_tmp = val_put

  IF ( ALL( var(:) == val_exp ) ) THEN 
     IF ( is_key <= 0 ) THEN
       CALL getin_p(key_wd,  val_tmp)
     ENDIF
     var(:) = val_tmp
  END IF
  
END SUBROUTINE r10setvar_p


!!  =============================================================================================================================
!! SUBROUTINE:    r11setvar_p
!!
!>\BRIEF	  Subroutine for initializing an real 1D array with another real 1D array.
!!
!! DESCRIPTION:	  Subroutine for initializing an real 1D array with another real 1D array.
!!                This subroutine must be called by all processes.
!! \n
!_ ==============================================================================================================================
SUBROUTINE r11setvar_p (var, val_exp, key_wd, val_put, is_grid)
  
  REAL(r_std), DIMENSION(:), INTENT(inout)     :: var                 !! 1D real array to modify
  REAL(r_std), INTENT(in)                      :: val_exp             !! Exceptional value
  CHARACTER(LEN=*), INTENT(in)                 :: key_wd              !! The Key word we will look for
  REAL(r_std), DIMENSION(:), INTENT(in)        :: val_put             !! 1D integer array to stored
  LOGICAL,        OPTIONAL                     :: is_grid             !! Parameter present indicates a setvar for a grid variable 

  REAL(r_std), ALLOCATABLE,DIMENSION(:)        :: val_tmp
  INTEGER(i_std)                               :: is_key

  is_key = MAX(INDEX(key_wd, 'NO_KEYWORD'), INDEX(key_wd, 'NOKEYWORD'))
   
  IF (long_print_setvar_p) WRITE (numout,*) "r11setvar :", key_wd, val_exp, SIZE(val_put), val_put(1)

  ALLOCATE(val_tmp(SIZE(val_put)))
  val_tmp(:) = val_put(:)

  IF ( ALL( var(:) == val_exp ) ) THEN 
     IF ( is_key <= 0 ) THEN
       CALL getin_p(key_wd,  val_tmp)
     ENDIF
     var(:) = val_tmp (:)
  END IF

  DEALLOCATE(val_tmp)
  
END SUBROUTINE r11setvar_p


!!  =============================================================================================================================
!! SUBROUTINE:    r20setvar_p
!!
!>\BRIEF	  Subroutine for initializing an real 2D variable with a scalar real variable.
!!
!! DESCRIPTION:	  Subroutine for initializing an real 2D variable with a scalar real variable.
!!                This subroutine must be called by all processes.
!! \n
!_ ==============================================================================================================================
SUBROUTINE r20setvar_p (var, val_exp, key_wd, val_put)
  
  REAL(r_std), DIMENSION(:,:), INTENT(inout)   :: var                  !! 2D integer array to modify
  REAL(r_std), INTENT(in)                      :: val_exp              !! Exceptional value
  CHARACTER(LEN=*), INTENT(in)                 :: key_wd               !! The Key word we will look for
  REAL(r_std), INTENT(in)                      :: val_put              !! Scalar value to stored
 
  REAL(r_std)                                  :: val_tmp  
  INTEGER(i_std)                               :: is_key

  is_key = MAX(INDEX(key_wd, 'NO_KEYWORD'), INDEX(key_wd, 'NOKEYWORD'))
 
  IF (long_print_setvar_p) WRITE (numout,*) "r20setvar :", key_wd, val_exp, val_put

  val_tmp = val_put

  IF ( ALL( var(:,:) == val_exp ) ) THEN 
     IF ( is_key <= 0 ) THEN
       CALL getin_p(key_wd,  val_tmp)
     ENDIF
     var(:,:) = val_tmp
  END IF
  
END SUBROUTINE r20setvar_p

!!  =============================================================================================================================
!! SUBROUTINE:    r21setvar_p
!!
!>\BRIEF	  Subroutine for initialieing an 2D real variable with a 1D array real.
!!
!! DESCRIPTION:	  Subroutine for initialieing an 2D real variable with a 1D array real.
!!                This subroutine must be called by all processes.
!!                Row or column depending size of 1D array to stored.
!!
!!                example: 1D 1,2,3     2D is 1, 2, 3,
!!                                            1, 2, 3
!!
!!                example: 1D 1,2,3     2D is 1, 1,
!!                                            2, 2,
!!                                            3, 3
!! \n
!_ ==============================================================================================================================
SUBROUTINE r21setvar_p (var, val_exp, key_wd, val_put, is_grid)
  
  REAL(r_std), DIMENSION(:,:), INTENT(inout)   :: var                  !! 2D real array to modify
  REAL(r_std), INTENT(in)                      :: val_exp              !! Exceptional value
  CHARACTER(LEN=*), INTENT(in)                 :: key_wd               !! The Key word we will look for
  REAL(r_std), DIMENSION(:), INTENT(in)        :: val_put              !! 1D real array to stored
  LOGICAL,        OPTIONAL                     :: is_grid              !! Parameter present indicates a setvar for a grid variable 

  REAL(r_std), ALLOCATABLE,DIMENSION(:)        :: val_tmp
  INTEGER(i_std)                               :: is_key

  is_key = MAX(INDEX(key_wd, 'NO_KEYWORD'), INDEX(key_wd, 'NOKEYWORD'))
  
  ! test if the 1D array dimension is compatible with first or second 
  ! dimension of the 2D array

  IF (long_print_setvar_p) WRITE (numout,*) "r21setvar :", key_wd, val_exp, SIZE(val_put), val_put(1)

  ALLOCATE(val_tmp(SIZE(val_put)))
  val_tmp(:) = val_put(:)

  IF (SIZE(val_put)==SIZE(var,1)) THEN 
      !
      ! example: 1D 1.,2.,3.     2D is 1., 2., 3.,
      !                                1., 2., 3.
      !
      IF ( ALL( var(:,:) == val_exp ) ) THEN 
         IF ( is_key <= 0 ) THEN
           CALL getin_p(key_wd,  val_tmp)
         ENDIF
         var(:,:) = SPREAD(val_tmp(:),2,SIZE(var,1))
      END IF
  ELSEIF (SIZE(val_put)==SIZE(var,2)) THEN 
      !
      ! example: 1D 1.,2.,3.     2D is 1., 1.,
      !                                2., 2.,
      !                                3., 3.
      !
      IF ( ALL( var(:,:) == val_exp ) ) THEN 
         IF ( is_key <= 0 ) THEN
           CALL getin_p(key_wd,  val_tmp)
         ENDIF
         var(:,:) = SPREAD(val_tmp(:),1,SIZE(var,1))
      END IF
  ELSE 
      WRITE (numout,*) ' incompatible dimension var and val_put'
      WRITE (numout,*) ' var     ', SIZE(var,1), SIZE(var,2)
      WRITE (numout,*) ' val_put ', SIZE(val_put)
      STOP 'setvar'
  END IF

  DEALLOCATE(val_tmp)
  
END SUBROUTINE r21setvar_p


!!  =============================================================================================================================
!! SUBROUTINE:    r22setvar_p
!!
!>\BRIEF	  Subroutine for initializing a 2D real variable with a real with the same size. 
!!
!! DESCRIPTION:	  Subroutine for initializing a 2D real variable with a real with the same size or by reading an scalar value
!!                from run.def if key_wd is different from "NO_KEYWORD" or "NOKEYWORD". 
!!                It is not possible to read a 2D variable from run.def.
!!                This subroutine must be called by all processes.
!! \n
!_ ==============================================================================================================================
SUBROUTINE r22setvar_p (var, val_exp, key_wd, val_put)

  REAL(r_std), DIMENSION(:,:), INTENT(inout)   :: var                 !! 2D real array to modify
  REAL(r_std), INTENT(in)                      :: val_exp             !! Exceptional value
  CHARACTER(LEN=*), INTENT(in)                 :: key_wd              !! The Key word we will look for
  REAL(r_std), DIMENSION(:,:), INTENT(in)      :: val_put             !! 2D integer array to stored
  REAL(r_std), ALLOCATABLE, DIMENSION(:,:)     :: val_tmp
  REAL(r_std)                                  :: val_scal            !! Temporary variable to read a scalar value from run.def
  INTEGER(i_std)                               :: is_key

  is_key = MAX(INDEX(key_wd, 'NO_KEYWORD'), INDEX(key_wd, 'NOKEYWORD'))

  IF (long_print_setvar_p) WRITE (numout,*) "r22setvar :", key_wd, val_exp, SIZE(val_put), val_put(1,1)

  ALLOCATE(val_tmp(SIZE(val_put,DIM=1),SIZE(val_put,DIM=2)))
  val_tmp(:,:) = val_put(:,:)

  IF ( ALL( var(:,:) == val_exp ) ) THEN 
     IF ( is_key <= 0 ) THEN
        ! This case only read a scalar value with getin
        val_scal=val_exp
        CALL getin_p(key_wd, val_scal)
        ! If a value was found in run.def, then set val_tmp to this value.
        IF (val_scal/=val_exp) val_tmp(:,:)=val_scal 
     ENDIF
     var(:,:) = val_tmp(:,:)
  END IF

  DEALLOCATE(val_tmp)
  
END SUBROUTINE r22setvar_p

!!  =============================================================================================================================
!! SUBROUTINE:    r30setvar_p
!!
!>\BRIEF	  Subroutine for initializing an real 3D variable with a scalar real variable.
!!
!! DESCRIPTION:	  Subroutine for initializing an real 3D variable with a scalar real variable.
!!                This subroutine must be called by all processes.
!! \n
!_ ==============================================================================================================================
SUBROUTINE r30setvar_p (var, val_exp, key_wd, val_put)

  REAL(r_std), DIMENSION(:,:,:), INTENT(inout) :: var                  !! 3D integer array to modify
  REAL(r_std), INTENT(in)                      :: val_exp              !! Exceptional value
  CHARACTER(LEN=*), INTENT(in)                 :: key_wd               !! The Key word we will look for
  REAL(r_std), INTENT(in)                      :: val_put              !! Scalar value to stored

  REAL(r_std)                                  :: val_tmp 
  INTEGER(i_std)                               :: is_key

  is_key = MAX(INDEX(key_wd, 'NO_KEYWORD'), INDEX(key_wd, 'NOKEYWORD'))

  IF (long_print_setvar_p) WRITE(numout,*) 'r30setvar',val_exp, val_put

  val_tmp = val_put

  IF ( ALL( var(:,:,:) == val_exp ) ) THEN
     IF ( is_key <= 0 ) THEN
       CALL getin_p(key_wd,  val_tmp)
     ENDIF
     var(:,:,:) = val_tmp
  END IF

END SUBROUTINE r30setvar_p

END MODULE sechiba_io_p
