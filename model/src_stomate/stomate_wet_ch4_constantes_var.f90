! =================================================================================================================================
! MODULE       : stomate_wet_ch4_constantes_var
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       CH4_calcul main module
!!
!!\n DESCRIPTION : None
!!
!! RECENT CHANGE(S) : None
!!
!! REFERENCE(S)	: None
!!
!! SVN :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_stomate/stomate.f90 $
!! $Date: 2017-07-28 17:48:13 +0200 (Fri, 28 Jul 2017) $
!! $Revision: 4542 $
!! \n
!_ ================================================================================================================================
MODULE stomate_wet_ch4_constantes_var
  ! modules used:
  USE ioipsl_para

  IMPLICIT NONE

  ! private & public routines

  PUBLIC
!
! stomate cste WETLAND 
!
  ! Nb of vertical layers for CH4 diffussion
  INTEGER(i_std),SAVE  :: nvert = 171 
  INTEGER(i_std),SAVE  :: ns = 151
  INTEGER(i_std),SAVE  :: nday = 24
  REAL(r_std),SAVE  :: h = 0.1
  REAL(r_std),SAVE  :: rk = 1
  REAL(r_std),PARAMETER  :: rkh = 100 !rk/h**2
  REAL(r_std),SAVE  :: diffair = 7.2
  REAL(r_std),SAVE  :: pox = 0.5
  REAL(r_std),SAVE  :: dveg = 0.001
  REAL(r_std),SAVE  :: rkm = 5.0
  REAL(r_std),SAVE  :: xvmax = 20.0
  REAL(r_std),SAVE  :: oxq10 = 2.0
  REAL(r_std),PARAMETER  :: funit = 3.84 !3.84/rk
  REAL(r_std),SAVE  :: scmax = 500.
  REAL(r_std),SAVE  :: sr0pl = 600.

!valeur de WTD pour les routines de calcul de densite de flux de CH4
  REAL(r_std),SAVE  :: pwater_wet1=-3
  REAL(r_std),SAVE  :: pwater_wet2=-9
  REAL(r_std),SAVE  :: pwater_wet3=-15
  REAL(r_std),SAVE  :: pwater_wet4=-21

  REAL(r_std),SAVE  :: rpv = 0.5
  REAL(r_std),SAVE  :: iother = -1.0

!!pour l instant je les mets constantes pour toutes les latitudes
  REAL(r_std),SAVE  :: rq10 = 3.0

  REAL(r_std),SAVE , DIMENSION(3) :: alpha_CH4 = (/0.006,0.004,0.028/)

!FLAG for CH4 from wetland
  LOGICAL,SAVE                            :: CH4_WTD1, CH4_WTD2, CH4_WTD3, CH4_WTD4

!atmoshpere methane concentration/ near surface methane concentration
  REAL(r_std), SAVE            :: CH4atmo_CONC

CONTAINS


END MODULE stomate_wet_ch4_constantes_var

