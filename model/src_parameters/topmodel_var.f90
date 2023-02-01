! =================================================================================================================================
! MODULE       : topmodel_var
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        
!!
!!\n DESCRIPTION: 
!!                
!! RECENT CHANGE(S):
!!
!! REFERENCE(S)	: 
!!
!! SVN          :
!! $HeadURL: $
!! $Date: 2018-06-29 14:56:38 +0200 (Fri, 29 Jun 2018) $
!! $Revision: 5316 $
!! \n
!_ ================================================================================================================================

MODULE topmodel_var

  USE defprec

  IMPLICIT NONE
!-

                         !-----------------------!
                         !  TOPMODEL CONSTANTS   !
                         !-----------------------!

  !
  ! FLAGS 
  !
!valeurs des bornes des differentes classes de WTD pour TOPMODEL
  REAL(r_std),SAVE :: WTD1_borne=0.03
  REAL(r_std),SAVE :: WTD2_borne=0.09
  REAL(r_std),SAVE :: WTD3_borne=0.15
  REAL(r_std),SAVE :: WTD4_borne=0.21
!!valeurs du shift de la distribution topo pour passer de fsat a fwet
  REAL(r_std),SAVE :: SHIFT_fsat_fwet=5.
!pss:-

END MODULE topmodel_var 
