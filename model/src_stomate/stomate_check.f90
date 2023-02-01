! =================================================================================================================================
! MODULE       : stomate_check
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Subroutines to check mass conservation 
!!	
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! SVN :
!! $HeadURL: $
!! $Date: $
!! $Revision: $
!! \n
!================================================================================================================================


MODULE stomate_check

  USE xios_orchidee
  USE constantes
  USE stomate_data
  USE ioipsl_para
  USE constantes_soil
  USE ieee_arithmetic 

  IMPLICIT NONE

  ! private routines
  PRIVATE

  ! public routines
  PUBLIC stomate_check_mass_values, stomate_check_cons_mass

  INTERFACE stomate_check_cons_mass
    MODULE PROCEDURE stomate_check_cons_mass_r2d, stomate_check_cons_mass_r3d
  END INTERFACE


CONTAINS 


!! ================================================================================================================================
!! SUBROUTINE 	: stomate_check_mass 
!!
!>\BRIEF       Check for biomass variable consistency  
!!
!! DESCRIPTION  : Check for biomass values: negatives, very big values and Nan
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): None 
!!
!! REFERENCES	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE stomate_check_mass_values(npts,biomass,identifier)

  INTEGER(i_std), INTENT(in)                    :: npts         !! Domain size (unitless)
  REAL(r_std), DIMENSION(:,:,:,:),INTENT(in)    :: biomass      !! Biomass components of the model tree  
  
  CHARACTER(LEN=*), INTENT(in),OPTIONAL         ::  identifier  !! string with to identify where this routine was called form
  
  ! locals
  INTEGER                                       :: nelement, npart, pft, idx
  LOGICAL                                       :: is_crash 
  CHARACTER(LEN=:), ALLOCATABLE                 :: ident_routine
! ===================================================================================================


   ident_routine='No information'
   IF (PRESENT(identifier)) THEN
     ident_routine=identifier
   ENDIF

   ! 1. check for negative biomass pools:
   IF (printlev >= 1) WRITE (numout,*) '=== BIOMASS CHECK ==='

   is_crash=.FALSE.
   DO nelement=1,nelements
     DO npart=1,nparts
        ! dont check crops and crops, it does not work yet
        DO pft=1,nvm-2
           DO idx=1,npts
             IF ((biomass(idx,pft,npart,nelement) .LT. zero).OR. & ! negative
                   (biomass(idx,pft,npart,nelement) .GT. large_value).OR. & ! infinity 
                   (isnan(biomass(idx,pft,npart,nelement))))   THEN ! NaN
                        WRITE (numout,*) 'FATAL negative biomass detected'
                        WRITE (numout,*) 'routine:'  , ident_routine
                        WRITE (numout,*) 'biomass = ', biomass(idx,pft,npart,nelement)
                        WRITE (numout,*) 'box     = ',idx
                        WRITE (numout,*) 'PFT     = ',pft
                        WRITE (numout,*) 'npart   = ',npart
                        WRITE (numout,*) '1=leaf,2=sapabov,3=sapbelo,4=heartabov,5=heatbelow,6=root,7=fruit,8=carbres,9=ilabile'
                        WRITE (numout,*) 'element = ',nelement
                        is_crash=.TRUE.
             ENDIF
           ENDDO
        ENDDO
     ENDDO
   ENDDO


   IF (is_crash) THEN
      CALL ipslerr_p ( 3, 'stomate_check_mass_values',   &
                          'Negative biomass pool was detected:', &
                          'Failed for subroutine ', &
                          identifier)
   ELSE
      IF (printlev >= 1) WRITE (numout,*) '     passed          '
   ENDIF

  END SUBROUTINE stomate_check_mass_values

!! ================================================================================================================================
!! SUBROUTINE 	: stomate_check_cons_mass_r3d 
!!
!>\BRIEF        
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): None 
!!
!! REFERENCES	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE stomate_check_cons_mass_r3d(lalo, mass_before,      &
                          mass_after,       &
                          mass_change,      &
                          identifier)

  REAL(r_std),DIMENSION(:,:),INTENT(in)              :: lalo          !! Geographical coordinates (latitude,longitude) for pixels (degrees) 
  REAL(r_std), DIMENSION(:,:,:),INTENT(in)           :: mass_before   !! Biomass old
  REAL(r_std), DIMENSION(:,:,:),INTENT(in)           :: mass_after    !! Biomass new (npts,nvm,nelements)
  REAL(r_std), DIMENSION(:,:,:),INTENT(in)           :: mass_change   !! Biomass change
  CHARACTER(LEN=*), INTENT(in),OPTIONAL              ::  identifier   !! string with to identify where this routine was called form

  ! locals
  INTEGER                                            :: nelement, pft, idx, ierr
  LOGICAL,ALLOCATABLE, DIMENSION(:,:,:)              :: is_crash
  REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)         :: mass_conserv  !! conservation 
  CHARACTER(LEN=:), ALLOCATABLE                      :: ident_routine
  ! this subroutine checks for mass conservation
  INTEGER(i_std)                                     :: npts          !! Domain size (unitless)
  INTEGER(i_std)                                     :: nvm           !! Domain size (unitless)
  INTEGER(i_std)                                     :: nelements     !! Domain size (unitless)
  
! ===================================================================================================

  npts = SIZE(mass_before, DIM=1)
  nvm = SIZE(mass_before, DIM=2)
  nelements = SIZE(mass_before, DIM=3)

  IF (PRESENT(identifier)) THEN
      ident_routine=identifier
  ELSE
       ident_routine=' no information'
  ENDIF

  ALLOCATE(is_crash(npts, nvm, nelements), stat=ierr)
  IF (ierr /= 0) CALL ipslerr_p(3, 'stomate_chekc_cons_mass_r3d', 'Memory allocation problem with ', 'is_crash variable', '')

  ALLOCATE(mass_conserv(npts, nvm, nelements), stat=ierr)
  IF (ierr /= 0) CALL ipslerr_p(3, 'stomate_check_cons_mass_r3d', 'Memory allocation problem with ', 'mass_conserv variable', '')

  is_crash(:,:,:)=.FALSE.
  ! bookkeeping:
  mass_conserv(:,:,:)          = mass_before(:,:,:) - mass_after(:,:,:) + mass_change(:,:,:)

  ! check the bookkeeping
  DO nelement=1, nelements-1
    ! take out crops and grass; till they are fixed
    DO pft=1, nvm-2
      DO idx=1, npts
         IF (ABS(mass_conserv(idx,pft,nelement)).GT. (min_stomate*100.)) THEN
                WRITE (numout,*) 'FATAL mass conservation failed (positive value = leak)'
                WRITE (numout,*) 'routine:', ident_routine
                ! the limit is given by the allocation routine, which accepts
                ! mass leaks<ABS(10E-6) 
                WRITE (numout,*) ' limit: ',min_stomate*100.
                WRITE (numout,*) ' Coordinates :', lalo(idx,:) 
                WRITE (numout,*) ' for element :', nelement, ' of nelements: ',nelements
                WRITE (numout,*) ' gridpoint: ',idx , ' of ngrids: ',npts
                WRITE (numout,*) ' PFT: ',pft , ' of npfts: ',nvm
                WRITE (numout,*) ' mismatch =', mass_conserv(idx,pft,nelement)
                WRITE (numout,*) ' mass(before) =', mass_before(idx,pft,nelement)
                WRITE (numout,*) ' mass(after) =', mass_after(idx,pft,nelement)
                is_crash(idx,pft,nelement)=.TRUE.
         ENDIF  
      ENDDO
    ENDDO
  ENDDO

  IF (ANY(is_crash(:,:,:))) THEN
     CALL ipslerr_p ( 3, 'stomate_check_cons_mass_r3d', &
                    'Mass conservation is not preserved', &
          &         'in subroutine', &
                    ident_routine)
  ENDIF

  DEALLOCATE(mass_conserv)
  DEALLOCATE(is_crash)

  END SUBROUTINE stomate_check_cons_mass_r3d

!! ================================================================================================================================
!! SUBROUTINE 	: stomate_check_cons_mass_r2d
!!
!>\BRIEF        
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): None 
!!
!! REFERENCES	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE stomate_check_cons_mass_r2d(lalo, mass_before,      &
                          mass_after,       &
                          mass_change,      &
                          identifier)


  REAL(r_std),DIMENSION(:,:),INTENT(in)              :: lalo          !! Geographical coordinates (latitude,longitude) for pixels (degrees) 
  ! biomass pools to be checked
  REAL(r_std), DIMENSION(:,:),INTENT(in)             :: mass_before   !! Biomass old
  REAL(r_std), DIMENSION(:,:),INTENT(in)             :: mass_after    !! Biomass new (npts,nvm)
  REAL(r_std), DIMENSION(:,:),INTENT(in)             :: mass_change   !! Biomass change

  CHARACTER(LEN=*), INTENT(in),OPTIONAL              ::  identifier   !! string with to identify where this routine was called form

  ! locals
  INTEGER(i_std)                                     :: pft, idx, ierr
  LOGICAL,ALLOCATABLE,DIMENSION(:,:)                 :: is_crash
  REAL(r_std),ALLOCATABLE, DIMENSION(:,:)            :: mass_conserv  !! conservation 
  CHARACTER(LEN=:), ALLOCATABLE                      :: ident_routine
  ! this subroutine checks for mass conservation
  INTEGER(i_std)                                     :: npts          !! Domain size (unitless)
  INTEGER(i_std)                                     :: nvm           !! Domain size (unitless)
! ===================================================================================================
  npts = SIZE(mass_before, DIM=1)
  nvm = SIZE(mass_before, DIM=2)

  ALLOCATE(is_crash(npts, nvm), stat=ierr)
  IF (ierr /= 0) CALL ipslerr_p(3, 'stomate_check_cons_mass_r2d', 'Memory allocation problem with ', 'is_crash variable', '')

  ALLOCATE(mass_conserv(npts, nvm), stat=ierr)
  IF (ierr /= 0) CALL ipslerr_p(3, 'stomate_check_cons_mass_r2d', 'Memory allocation problem with ', 'mass_conserv variable', '')

  IF (PRESENT(identifier)) THEN
      ident_routine=identifier
  ELSE
       ident_routine=' no information'
  ENDIF


  is_crash(:,:)=.FALSE.
  ! bookkeeping:
  mass_conserv(:,:)          = mass_before(:,:) - mass_after(:,:) + mass_change(:,:)

  ! check the bookkeeping
  ! take out crops and grass; till they are fixed
  DO pft=1,nvm-2
    DO idx=1,npts
       IF (ABS(mass_conserv(idx,pft)).GT. (min_stomate*100.)) THEN
              WRITE (numout,*) 'FATAL mass conservation failed (positive value = leak)'
              WRITE (numout,*) 'routine:', ident_routine
                ! the limit is given by the allocation routine, which accepts
                ! mass leaks<ABS(10E-6) 
              WRITE (numout,*) ' limit: '       ,min_stomate*100.
              WRITE (numout,*) ' Coordinates :' , lalo(idx,:) 
              WRITE (numout,*) ' gridpoint: '   ,idx , ' of ngrids: ',npts
              WRITE (numout,*) ' PFT: '         ,pft , ' of npfts: ',nvm
              WRITE (numout,*) ' mismatch ='    , mass_conserv(idx,pft)
              WRITE (numout,*) ' mass(before) =', mass_before(idx,pft)
              WRITE (numout,*) ' mass(after) =' , mass_after(idx,pft)
              is_crash(idx,pft)=.TRUE.
       ENDIF  
    ENDDO
  ENDDO

  IF (ANY(is_crash(:,:))) THEN
     CALL ipslerr_p ( 3, 'stomate_check_cons_mass_r2d', &
                    'Mass convervation is not preserved', &
          &         'in subroutine', &
                    ident_routine)
  ENDIF

  DEALLOCATE(mass_conserv)
  DEALLOCATE(is_crash)

  END SUBROUTINE stomate_check_cons_mass_r2d



END MODULE stomate_check
