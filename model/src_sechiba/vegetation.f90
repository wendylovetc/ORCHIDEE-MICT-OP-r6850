!  ==============================================================================================================================\n
!  MODULE 	: sechiba
! 
!  CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
!  LICENCE      : IPSL (2006)
!  This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        
!!
!!\n DESCRIPTION  : 
!! 
!! RECENT CHANGE(S): None 
!! 
!! REFERENCE(S) : None
!!   
!! SVN     :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_sechiba/sechiba.f90 $ 
!! $Date: 2018-12-13 14:18:58 +0100 (Thu, 13 Dec 2018) $
!! $Revision: 5672 $
!! \n
!_ ================================================================================================================================
 
MODULE vegetation 
 
  USE ioipsl_para

  IMPLICIT NONE

  PRIVATE
  PUBLIC vegetation_nobio_to_bg

  INTEGER(i_std), SAVE                             :: printlev_loc   !! local printlev for this module
!$OMP THREADPRIVATE(printlev_loc)

CONTAINS



!! ==============================================================================================================================\n
!! SUBROUTINE 	: vegetation_nobio_to_bg
!!
!>\BRIEF        Place all nobio to the bareground pft
!! 
!! DESCRIPTION  : 
!!
!! RECENT CHANGE(S): 
!! 
!! MAIN OUTPUT VARIABLE(S): 
!! 
!! REFERENCE(S)	: None
!! 
!! FLOWCHART    : None
!! \n
!! ================================================================================================================================ 
  SUBROUTINE vegetation_nobio_to_bg(veget_nobio, veget_bg)

    REAL(r_std), DIMENSION(:,:), INTENT(in) :: veget_nobio  !! veget_max with nobio (sum(pft)<1)
    REAL(r_std), DIMENSION(:,:), INTENT(out) :: veget_bg     !! veget_max with full veget (sum(pft)=1)

    INTEGER(i_std) :: nvm

    nvm = SIZE(veget_nobio, DIM=2)

    veget_bg(:,2:nvm) = veget_nobio(:,2:nvm)
    veget_bg(:,1) = MAX((un - SUM(veget_nobio(:,2:nvm), 2)), zero)

  END SUBROUTINE vegetation_nobio_to_bg


END MODULE vegetation 

