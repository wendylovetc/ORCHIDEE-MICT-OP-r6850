! ===================================================================================================\n
! MODULE        : equal_division
!
! CONTACT       : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE       : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        
!!
!!\n DESCRIPTION : Apply an integer equal division where the reminder needs to be taken
!!                 into account. 
!!                For example: 18 elements and make 4 groups -> 4, 4, 5, 5
!!
!! RECENT CHANGE(S) : None
!!
!! REFERENCE(S) :
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_sechiba/hydrol.f90 $
!! $Date: 2018-10-09 11:43:16 +0200 (Tue, 09 Oct 2018) $
!! $Revision: 5477 $
!! \n
!_ ===============================================================================================\n
MODULE equal_division 
  !---------------------------------------------------------------------
  !-
  !---------------------------------------------------------------------
  USE time
  USE ioipsl_para
  !-
  IMPLICIT NONE
  !-
  PRIVATE
  PUBLIC t_equal_division 
  !-
  ! 
  !-
  TYPE t_equal_division
    INTEGER ::  ngroups         !! Number of groups to divide
    INTEGER ::  total           !! Total Number of elements 
    INTEGER, DIMENSION(:), ALLOCATABLE :: distribution !! Number of elements per group
    INTEGER, DIMENSION(:), ALLOCATABLE :: accumulation !! Accumulation number of elements per group
  CONTAINS
    PROCEDURE :: distribute => equal_division_distribute !! Calculate equal division, private
    PROCEDURE :: get_distr => equal_division_get_distribution  !! Get the values

    PROCEDURE :: accumulate => equal_division_accumulation  !! private
    PROCEDURE :: get_accum => equal_division_get_accumulate !!

    PROCEDURE :: get => equal_division_get_group
  END TYPE t_equal_division
  !-
  INTERFACE t_equal_division
    MODULE PROCEDURE :: new_t_equal_division
  END INTERFACE t_equal_division
  !-
CONTAINS
  !-
  !===
  !-
  !! ================================================================================================================================
  !! SUBROUTINE 	: new_t_equal_division
  !!
  !>\BRIEF        Constructor for t_equal_division
  !!
  !! DESCRIPTION  : Inputs arguments define:
  !!            elements / parts
  !!                
  !!                
  !! \n
  !_ ================================================================================================================================
  TYPE(t_equal_division) FUNCTION new_t_equal_division(parts, elements) 
    ! Input Variables
    INTEGER, INTENT(in) :: parts        !! Groups number 
    INTEGER, INTENT(in) :: elements     !! Elements to distribute

    ! Output Variables

    ! Local Variables
    INTEGER :: ier
    CHARACTER(LEN=:), ALLOCATABLE :: subrname
  !_ ===============================================================================================

    subrname = "new_t_equal_division"

    ! Are inputs arguments consistent?
    IF (parts < 0) THEN
      CALL ipslerr_p(3, subrname,'parts must be positive (>0)',&
                'Value found',parts)
    ENDIF
    IF (elements < 0) THEN
      CALL ipslerr_p(3, subrname,'elements must be positive (>0)',&
                'Value found',elements)
    ENDIF
    IF (parts > elements) THEN
      CALL ipslerr_p(3, subrname,'Part must be bigger than elements',&
                'Value found:', parts)
    ENDIF

    ! Allocate vars
    ALLOCATE(new_t_equal_division%distribution(parts), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3, subrname,&
                                'Memory allocation error for distribution array',&
                                'Error code:',ier) 

    ALLOCATE(new_t_equal_division%accumulation(parts), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3, subrname,&
                                'Memory allocation error for distribution array',&
                                'Error code:',ier) 

    ! Init values
    new_t_equal_division%ngroups = parts
    new_t_equal_division%total = elements

    CALL new_t_equal_division%distribute()
    CALL new_t_equal_division%accumulate()

  END FUNCTION new_t_equal_division
  !-
  !! ================================================================================================================================
  !! SUBROUTINE 	: equal_division_distribute
  !!
  !>\BRIEF        Apply equal distribution
  !!
  !! DESCRIPTION  : Element will be distributed across each group. The reminder
  !!                 will be place from the back to the front of the array.
  !!               
  !!                
  !!                
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE equal_division_distribute(this)
    ! Input Variables
    CLASS(t_equal_division), INTENT(inout) :: this

    ! Output Variables

    ! Local Variables
    INTEGER :: remainder, ite, revidx
    REAL :: elems_per_group


    ! Equal division
    elems_per_group = this%total / this%ngroups
    this%distribution(:) = FLOOR(elems_per_group)

    ! Distribute remainder values backwards
    remainder = MODULO(this%total, this%ngroups)
    IF (remainder .GT. 0) THEN
      DO ite=1, remainder
        revidx = this%ngroups - ite + 1
        this%distribution(revidx) = this%distribution(revidx) + 1
      ENDDO
    ENDIF
    
  END SUBROUTINE equal_division_distribute

  !! ================================================================================================================================
  !! SUBROUTINE 	: equal_division_get_distribution
  !!
  !>\BRIEF        Getter
  !!
  !! DESCRIPTION  : 
  !!                
  !!                
  !! \n
  !_ ================================================================================================================================
  FUNCTION equal_division_get_distribution(this) RESULT(out_array)
    ! Input Variables
    CLASS(t_equal_division), INTENT(inout) :: this

    ! Output Variables
    INTEGER, DIMENSION(:), ALLOCATABLE :: out_array

    ! Local Variables
    out_array  = this%distribution
    
  END FUNCTION equal_division_get_distribution

  !-
  !! ================================================================================================================================
  !! SUBROUTINE 	: equal_division_accumulation
  !!
  !>\BRIEF        Calculate the accumulated at each group 
  !!
  !! DESCRIPTION  : This is an intermediate calculation 
  !!               
  !!                
  !!                
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE equal_division_accumulation(this)
    ! Input Variables
    CLASS(t_equal_division), INTENT(inout) :: this

    ! Output Variables

    ! Local Variables
    INTEGER :: remainder, ite

    this%accumulation = this%distribution
    DO ite = 2, this%ngroups
      this%accumulation(ite) = this%distribution(ite) + this%accumulation(ite - 1)
    ENDDO
    
  END SUBROUTINE equal_division_accumulation

  !! ================================================================================================================================
  !! SUBROUTINE 	: equal_division_get_accumulate
  !!
  !>\BRIEF        Getter
  !!
  !! DESCRIPTION  : 
  !!                
  !!               
  !!                
  !!                
  !! \n
  !_ ================================================================================================================================
  FUNCTION equal_division_get_accumulate(this) RESULT(out_array)
    ! Input Variables
    CLASS(t_equal_division), INTENT(inout) :: this

    ! Output Variables
    INTEGER, DIMENSION(:), ALLOCATABLE :: out_array

    ! Local Variables
    out_array  = this%accumulation
    
  END FUNCTION equal_division_get_accumulate

  !! ================================================================================================================================
  !! SUBROUTINE 	: equal_division_get_group
  !!
  !>\BRIEF        Find group number for a given value
  !!
  !! DESCRIPTION  : 
  !!                
  !!                
  !! \n
  !_ ================================================================================================================================
  FUNCTION equal_division_get_group(this, in_value) RESULT(group_number)
    ! Input Variables
    CLASS(t_equal_division), INTENT(inout)  :: this
    INTEGER, INTENT(in)                     :: in_value

    ! Output Variables
    INTEGER                                 :: group_number

    ! Local Variables
    INTEGER, PARAMETER                      :: NOT_FOUND = -1
    INTEGER                                 :: ite

    ! Init values
    group_number = NOT_FOUND

    IF (in_value .GT. this%total) THEN
      WRITE(*,*) "equal_division_get_group:: given in_value=", in_value
      WRITE(*,*) "equal_division_get_group:: max number=", this%total
      CALL ipslerr_p(3,'equal_division_get_group', &
                    'Given value is bigger than the maximum defined number', &
                    'Value found',in_value)
    ENDIF

    DO ite=1, this%ngroups
      IF ( in_value .LE. this%accumulation(ite) ) THEN
        group_number = ite
        EXIT
      ENDIF
    ENDDO
    
  END FUNCTION equal_division_get_group

END MODULE equal_division 
