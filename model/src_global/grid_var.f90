! ===============================================================================================================================
! MODULE       : grid_var
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF
!!
!! \n DESCRIPTION : This module define variables for the grid module. 
!!                The module is already USE in module grid. Therefor no need to use it seperatly if use grid is already done.

!!
!! RECENT CHANGE(S): These variables were previously in grid module. They have been moved here to avoid dependency 
!!                   problems when the variables are needed in the parallelization modules. 
!!
!!
!! SVN
!! $HeadURL$
!! $Date$
!! $Revision$
!! \n
!_ ================================================================================================================================

MODULE grid_var

  IMPLICIT NONE

  !=================================================================================
  !
  ! Horizontal grid information
  !
  !=================================================================================

  CHARACTER(LEN=20), SAVE                           :: GridType    !! Describes the grid it can be RegLonLat, RegXY or UnStruct
!$OMP THREADPRIVATE(GridType)

  CHARACTER(LEN=20), SAVE                           :: GridName    !! Name of the grid
!$OMP THREADPRIVATE(GridName)

  INTEGER, SAVE                                     :: NbSegments  !! Number of segments in 
                                                                   !! the polygone defining the grid box
!$OMP THREADPRIVATE(NbSegments)
  
  INTEGER, SAVE                                     :: NbNeighb    !! Number of neighbours
!$OMP THREADPRIVATE(NbNeighb)
  
END MODULE grid_var
