!
! Place here all those small routines that do not fit in orchidee logic yet
!
!
!< $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_global/interpol_help.f90 $
!< $Date: 2016-06-17 13:26:43 +0200 (Fri, 17 Jun 2016) $
!< $Author: albert.jornet $
!< $Revision: 3564 $
!
!
MODULE utils 

  ! Modules used :

  USE netcdf
  USE defprec
  USE ioipsl_para

  IMPLICIT NONE

  PRIVATE
  PUBLIC nccheck, check_lai_vs_assim, show_values, compare_results
  !
  INTERFACE show_values
    MODULE PROCEDURE show_values_r_scal, show_values_r1, show_values_r2, show_values_r3, show_values_r4
  END INTERFACE
  !
  INTERFACE compare_results
    MODULE PROCEDURE compare_results_r2, compare_results_r3 
  END INTERFACE
  !
  ! Show_values parameters
  ! 
  INTEGER(i_std), PARAMETER :: PX_VAL =   1 ! Gridcell/pixel value
  INTEGER(i_std), PARAMETER :: PFT_VAL =  3 ! PFT
  INTEGER(i_std), PARAMETER :: PROC_VAL = 0 ! Processor
  ! 
CONTAINS
  !
!! ================================================================================================================================
!! SUBROUTINE 	: nccheck 
!!
!>\BRIEF        Check for netcdf exit status 
!!
!! DESCRIPTION  : Launch an orchidee error message if status variable contains a netcdf error
!!
!! RECENT CHANGE(S) : None
!!
!! REFERENCE(S)	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE nccheck(status)
    INTEGER(i_std), INTENT (IN)         :: status
    CHARACTER(LEN=200)                  :: mesg
    
    IF(status /= nf90_noerr) THEN
      
      WRITE(numout, *) trim(nf90_strerror(status))
      CALL ipslerr_p(3, 'nccheck', 'Netcdf error', 'Check out_orchide_XXXX output files', 'for more information')
    END IF  
  END SUBROUTINE nccheck
!
!! ================================================================================================================================
!! SUBROUTINE 	: check_lai_vs_assim
!!
!>\BRIEF        Check consistency between lai and assim_param
!!
!! DESCRIPTION  : Launch an orchidee error message this 2 variables are not consistent. 
!!              This will prevent any kind of error in the condveg module.
!!
!! RECENT CHANGE(S) : None
!!
!! REFERENCE(S)	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE check_lai_vs_assim(kjpindex, lalo, veget_max, lai, vcmax, caller)

    USE pft_parameters_var

    INTEGER(i_std), INTENT(in)                               :: kjpindex    !! Domain size (unitless)
    REAL(r_std),DIMENSION (kjpindex,2),   INTENT (in)        :: lalo        !! Geographical coordinates
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)        :: veget_max   !! Maximum vegetation fraction of each PFT inside 
                                                                            !! the grid box (0-1, unitless) 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)        :: vcmax       !! min+max+opt temps, vcmax, vjmax
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)        :: lai         !! Leaf area index (m^2 m^{-2})
    CHARACTER(LEN=*) :: caller

    INTEGER(i_std) :: ji, jv
    CHARACTER(LEN=150)                    :: printstr, printstr2            !! For temporary uses
!_ ================================================================================================================================

    DO jv=2, nvm
      DO ji=1,kjpindex
         !
         IF ( ( veget_max(ji,jv) .GT. min_sechiba ) ) THEN

            IF ( (lai(ji,jv) .GT. 0.01) .AND.  vcmax(ji,jv) <= 0) THEN
               WRITE(printstr, *)  'check_lai_vs_assim::coordinates lat=', lalo(ji, 1),' long=', lalo(ji, 2), 'PFT=', jv, ",npts=", ji
               WRITE(printstr2, *) 'check_lai_vs_assim:: lai=', lai(ji,jv), '  assim_param=', vcmax(ji,jv)
               CALL ipslerr_p(3, caller//'::check_lai_vs_assim', 'assim_param must be bigger than 0 to be consistent with lai', & 
                      TRIM(printstr), TRIM(printstr2))
            ENDIF
         ENDIF
      ENDDO
    ENDDO

  END SUBROUTINE check_lai_vs_assim
!
!! ================================================================================================================================
!! SUBROUTINE 	: show_values
!!
!>\BRIEF        Print values of a specific pixel and pft (if applies)
!!
!! DESCRIPTION  : Print values of a specific pixel, pft(if applies) and processor. For this purpose, define the PARAMETERS below:
!!
!!              The use of PX_VAL, PFT_VAL and PROC_VAL determines which values prints for (gridcell, pft and processor)
!!                 The rest of the values will be summed to get a unique value
!!
!!              Those subroutines are meant for developpers to help them trace values
!!
!! RECENT CHANGE(S) : None
!!
!! REFERENCE(S)	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE show_values_r_scal(subrname, varname, var)
    CHARACTER(LEN=*), INTENT(in) :: subrname
    CHARACTER(LEN=*), INTENT(in) :: varname
    REAL(r_std), INTENT(in) :: var

    IF (mpi_rank == PROC_VAL) WRITE(numout, *) "show_values_r_scal::"//TRIM(subrname)//"::"//TRIM(varname)//"=", var

  END SUBROUTINE show_values_r_scal

  SUBROUTINE show_values_r1(subrname, varname, var)
    CHARACTER(LEN=*), INTENT(in) :: subrname
    CHARACTER(LEN=*), INTENT(in) :: varname
    REAL(r_std), DIMENSION(:), INTENT(in) :: var

    IF (mpi_rank == PROC_VAL) WRITE(numout, *) "show_values_r1::"//TRIM(subrname)//"::"//TRIM(varname)//"=", var(PX_VAL)

  END SUBROUTINE show_values_r1

  SUBROUTINE show_values_r2(subrname, varname, var)
    CHARACTER(LEN=*), INTENT(in) :: subrname
    CHARACTER(LEN=*), INTENT(in) :: varname
    REAL(r_std), DIMENSION(:,:), INTENT(in) :: var

    IF (mpi_rank == PROC_VAL) WRITE(numout, *) "show_values_r2::"//TRIM(subrname)//"::"//TRIM(varname)//"=", var(PX_VAL, PFT_VAL)

  END SUBROUTINE show_values_r2

  SUBROUTINE show_values_r3(subrname, varname, var, pftind)
    CHARACTER(LEN=*), INTENT(in) :: subrname
    CHARACTER(LEN=*), INTENT(in) :: varname
    REAL(r_std), DIMENSION(:,:,:), INTENT(in) :: var
    INTEGER(i_std), INTENT(in), OPTIONAL :: pftind ! pft index found in the array

    IF (.NOT. PRESENT(pftind) .OR. (PRESENT(pftind) .AND. pftind == 2)) THEN
       IF (mpi_rank == PROC_VAL) WRITE(numout, *) "show_values_r3::"//TRIM(subrname)//"::"//TRIM(varname)//"=", var(PX_VAL, PFT_VAL,1)
    ELSE IF ( pftind == 3) THEN
       IF (mpi_rank == PROC_VAL) WRITE(numout, *) "show_values_r3::"//TRIM(subrname)//"::"//TRIM(varname)//"=", var(PX_VAL, 1, PFT_VAL)
    ELSE
        CALL ipslerr_p(3,'show_values_r3','wrong pftind value','Allowed 2 or 3 but found:',pftind)
    ENDIF

  END SUBROUTINE show_values_r3

  SUBROUTINE show_values_r4(subrname, varname, var, pftind)
    CHARACTER(LEN=*), INTENT(in) :: subrname
    CHARACTER(LEN=*), INTENT(in) :: varname
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in) :: var
    INTEGER(i_std), INTENT(in), OPTIONAL :: pftind ! pft index found in the array

    IF (.NOT. PRESENT(pftind) .OR. (PRESENT(pftind) .AND. pftind == 2)) THEN
       IF (mpi_rank == PROC_VAL) WRITE(numout, *) "show_values_r4::"//TRIM(subrname)//"::"//TRIM(varname)//"=", SUM(SUM(var(PX_VAL,PFT_VAL,:,:),DIM=2))
    ELSE IF ( pftind == 3 ) THEN
       IF (mpi_rank == PROC_VAL) WRITE(numout, *) "show_values_r4::"//TRIM(subrname)//"::"//TRIM(varname)//"=", SUM(SUM(var(PX_VAL,:,PFT_VAL,:),DIM=2))
    ELSE IF ( pftind == 4 ) THEN
       IF (mpi_rank == PROC_VAL) WRITE(numout, *) "show_values_r4::"//TRIM(subrname)//"::"//TRIM(varname)//"=", SUM(SUM(var(PX_VAL,:,:,PFT_VAL),DIM=2))
    ELSE
        CALL ipslerr_p(3,'show_values_r4','wrong pftind value','Allowed 2,3 or 4 but found:',pftind)
    ENDIF

  END SUBROUTINE show_values_r4

!
!! ================================================================================================================================
!! SUBROUTINE 	: compare_results
!!
!>\BRIEF        It checks whether all values are the same in original and test
!!
!! DESCRIPTION  : It checks for all values of original and test array are the same.
!!              Otherwise an error will be raised and the wrong values will be printed.
!!
!!              This subroutines is meant to help developpers to ensure an improved
!!              subroutine is giving the same results
!!
!!
!! RECENT CHANGE(S) : None
!!
!! REFERENCE(S)	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE compare_results_r2(calleer, subroutine_test, var_test, original, test)
    CHARACTER(LEN=*), INTENT(in) :: calleer ! From where is called
    CHARACTER(LEN=*), INTENT(in) :: subroutine_test ! Which subroutine is tested
    CHARACTER(LEN=*), INTENT(in) :: var_test ! Variable to test
    REAL(r_std), DIMENSION(:,:), INTENT(in) :: original ! Expected values
    REAL(r_std), DIMENSION(:,:), INTENT(in) :: test ! New values to test

    INTEGER(i_std) :: ji, jv ! iterators
    INTEGER(i_std) :: first, second ! 
    LOGICAL :: are_equal ! are original and test arrays equal?

    are_equal = .FALSE.
    first = SIZE(original, DIM=1) ! first dimension size
    second = SIZE(original, DIM=2) ! second dimension size

    DO ji=1, first
        DO jv=1, second
           IF (original(ji,jv) .NE. test(ji,jv)) THEN
              WRITE(numout, *) var_test//"compare_results_r2:: i, j, original vs test=", &
                                ji,jv,original(ji,jv), test(ji,jv)
              
              are_equal = .TRUE.
           ENDIF
        ENDDO ! second
    ENDDO ! first

    IF (are_equal) THEN
        CALL ipslerr_p(3, 'compare_results_r2', 'Inside '//calleer, "Compared subroutines "//subroutine_test// "for variable "//var_test, 'Produce different values')
    ENDIF

  END SUBROUTINE compare_results_r2

  SUBROUTINE compare_results_r3(calleer, subroutine_test, var_test, original, test)
    CHARACTER(LEN=*), INTENT(in) :: calleer ! From where is called
    CHARACTER(LEN=*), INTENT(in) :: subroutine_test ! Which subroutine is tested
    CHARACTER(LEN=*), INTENT(in) :: var_test ! Variable to test
    REAL(r_std), DIMENSION(:,:,:), INTENT(in) :: original ! Expected values
    REAL(r_std), DIMENSION(:,:,:), INTENT(in) :: test ! New values to test

    INTEGER(i_std) :: ji, jv, jk ! iterators
    INTEGER(i_std) :: first, second, third ! 
    LOGICAL :: are_equal ! are original and test arrays equal?

    are_equal = .FALSE.
    first = SIZE(original, DIM=1) ! first dimension size
    second = SIZE(original, DIM=2) ! second dimension size
    third = SIZE(original, DIM=3) ! third dimension size

    DO ji=1, first
        DO jv=1, second
           DO jk=1, third
              IF (original(ji,jv,jk) .NE. test(ji,jv,jk)) THEN
                 WRITE(numout, *) var_test//"compare_results_r3:: i, j, k, original vs test=", &
                                ji,jv,jk,original(ji,jv,jk), test(ji,jv,jk)
              
                 are_equal = .TRUE.
              ENDIF
           ENDDO ! third
        ENDDO ! second
    ENDDO ! first

    IF (are_equal) THEN
        CALL ipslerr_p(3, 'compare_results_r3', 'Inside '//calleer, "Compared subroutines "//subroutine_test// "for variable "//var_test, 'Produce different values')
    ENDIF

  END SUBROUTINE compare_results_r3

END MODULE utils  
