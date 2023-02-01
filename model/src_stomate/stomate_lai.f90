! =================================================================================================================================
! MODULE       : stomate_lai 
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
!                This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Calculate lai 
!!
!! \n DESCRIPTION : Calculate lai 
!!
!! REFERENCE(S) :
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_stomate/lpj_crown.f90 $
!! $Date: 2015-11-16 14:26:03 +0100 (Mon, 16 Nov 2015) $
!! $Revision: 3026 $
!! \n
!_ ================================================================================================================================

MODULE stomate_lai 

!  USE ioipsl_para
!  USE stomate_data
!  USE constantes
  USE constantes_var
!  USE pft_parameters
  
  IMPLICIT NONE
  
  ! private & public routines

  PRIVATE
  PUBLIC setlai 
  
CONTAINS
  
!! ================================================================================================================================
!! SUBROUTINE 	: setlai
!!
!>\BRIEF        Routine to force the lai in STOMATE. The code in this routine
!! simply CALCULATES lai and is therefore not functional. The routine should be 
!! rewritten if one wants to force lai.
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): ::lai
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART : None
!! \n
!_ ================================================================================================================================
  
  SUBROUTINE setlai(biomass,sla_calc,lai)

  !! 0 Variable and parameter declaration 
  
    !! 0.1 Input variables

    REAL(r_std),DIMENSION(:,:,:,:),INTENT(in)    :: biomass !! biomass 
    REAL(r_std),DIMENSION(:,:),INTENT(in)        :: sla_calc !!  sla_calc
    
    !! 0.2 Output variables

    REAL(r_std),DIMENSION(:,:),INTENT(out)  :: lai  !! PFT leaf area index @tex $(m^{2} m^{-2})$ @endtex

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                               :: j    !! index (unitless)
!_ ================================================================================================================================
    
    !! 1. Set lai for bare soil to zero

    lai(:,ibare_sechiba) = zero

    !! 2. Multiply foliage biomass by sla to calculate lai for all PFTs and pixels

    DO j=2,SIZE(biomass,DIM=2)
       lai(:,j) = biomass(:,j,ileaf,icarbon)*sla_calc(:,j)
    ENDDO
    
  END SUBROUTINE setlai


END MODULE stomate_lai 
