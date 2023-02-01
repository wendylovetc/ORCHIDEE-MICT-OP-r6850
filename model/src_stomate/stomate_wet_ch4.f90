! =================================================================================================================================
! MODULE       : stomate_wet_ch4
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
MODULE stomate_wet_ch4

  ! modules used:
  USE stomate_wet_ch4_constantes_var
  USE stomate_wet_ch4_pt_ter_0
  USE stomate_wet_ch4_pt_ter_wet

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC stomate_wet_ch4_initialize, stomate_wet_ch4_main, stomate_wet_ch4_clear, stomate_wet_ch4_finalize, &
         stomate_wet_ch4_histdef, stomate_wet_ch4_config_parameters

  ! density flux of methane calculated for entire pixel (gCH4/dt/m**2)
  ! pour wetland avc Water Table Depth (WTD) = 0
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: ch4_flux_density_tot_0
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: ch4_flux_density_dif_0
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: ch4_flux_density_bub_0
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: ch4_flux_density_pla_0
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)  :: uo_0 !concentration dim = (kjpindex,nvert)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)  :: uold2_0 !concentration au pas de temps precedent
!pour wetland avc WTD = -x1
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: ch4_flux_density_tot_wet1
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: ch4_flux_density_dif_wet1
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: ch4_flux_density_bub_wet1
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: ch4_flux_density_pla_wet1
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)  :: uo_wet1
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)  :: uold2_wet1
!pour wetland avc WTD = -x2
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: ch4_flux_density_tot_wet2
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: ch4_flux_density_dif_wet2
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: ch4_flux_density_bub_wet2
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: ch4_flux_density_pla_wet2
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)  :: uo_wet2
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)  :: uold2_wet2
!pour wetland avc WTD = -x3
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: ch4_flux_density_tot_wet3
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: ch4_flux_density_dif_wet3
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: ch4_flux_density_bub_wet3
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: ch4_flux_density_pla_wet3
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)  :: uo_wet3
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)  :: uold2_wet3
!pour wetland avc WTD = -x4
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: ch4_flux_density_tot_wet4
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: ch4_flux_density_dif_wet4
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: ch4_flux_density_bub_wet4
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: ch4_flux_density_pla_wet4
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)  :: uo_wet4
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)  :: uold2_wet4


CONTAINS

!! ================================================================================================================================
!! SUBROUTINE   : stomate_wet_ch4_config_parameters
!!
!>\BRIEF        stomate cste WETLAND
!!
!! DESCRIPTION  : 
!!                
!! \n
!_ ================================================================================================================================
!!
  SUBROUTINE stomate_wet_ch4_config_parameters
       !Config Key   = nvert
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 171 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('NVERT',nvert)

       !Config Key   = ns
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 151 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('NS',ns)

       !Config Key   = nday
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 24
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('NDAY',nday)

       !Config Key   = h
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 0.1
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('H',h)

       !Config Key   = rk
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 1
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('RK',rk)

       !Config Key   = diffair
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 7.2 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('DIFFAIR',diffair)

       !Config Key   = pox
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 0.5 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('POX',pox)

       !Config Key   = dveg
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 0.001 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('DVEG',dveg)

       !Config Key   = rkm
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 5.0
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('RKM',rkm)

       !Config Key   = xvmax
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 20.0 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('XVMAX',xvmax)

       !Config Key   = oxq10
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 2.0 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('OXQ10',oxq10)

       !Config Key   = scmax
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 500. 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('SCMAX',scmax)

       !Config Key   = sr0pl
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 600. 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('SR0PL',sr0pl)

       !Config Key   = pwater_wet1
       !Config Desc  = depth where saturation: definition for wetland 1  
       !Config If    = CH4_CALCUL
       !Config Def   = -3 
       !Config Help  = 
       !Config Units = [cm]   
       CALL getin_p('PWATER_WET1',pwater_wet1)

       !Config Key   = pwater_wet2
       !Config Desc  = depth where saturation: definition for wetland 1  
       !Config If    = CH4_CALCUL
       !Config Def   = -9 
       !Config Help  = 
       !Config Units = [cm]   
       CALL getin_p('PWATER_WET2',pwater_wet2)

       !Config Key   = pwater_wet3
       !Config Desc  = depth where saturation: definition for wetland 1  
       !Config If    = CH4_CALCUL
       !Config Def   = -15 
       !Config Help  = 
       !Config Units = [cm]   
       CALL getin_p('PWATER_WET3',pwater_wet3)

       !Config Key   = pwater_wet4
       !Config Desc  = depth where saturation: definition for wetland 1  
       !Config If    = CH4_CALCUL
       !Config Def   = -21 
       !Config Help  = 
       !Config Units = [cm]   
       CALL getin_p('PWATER_WET4',pwater_wet4)

       !Config Key   = rpv
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 0.5 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('RPV',rpv)

       !Config Key   = iother
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = -1.0 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('IOTHER',iother)

       !Config Key   = rq10
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 3.0 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('RQ10',rq10)

       !Config Key   = alpha_CH4
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = /0.009,0.004,0.021/
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('ALPHA_CH4',alpha_CH4)

  END SUBROUTINE stomate_wet_ch4_config_parameters

!! ================================================================================================================================
!! SUBROUTINE   : stomate_wet_ch4_readrestart
!!
!>\BRIEF        Read restart variables
!!
!! DESCRIPTION  : 
!!                
!! \n
!_ ================================================================================================================================
!!
  SUBROUTINE stomate_wet_ch4_readrestart ( rest_id, itime )

    ! 0.4 Local variables
    CHARACTER(LEN=100)                   :: var_name !! Name of permafrost forcing file
    INTEGER(i_std), INTENT(in)           :: rest_id
    INTEGER(i_std), INTENT(in)           :: itime 

    INTEGER(i_std) ::  nivo
!_ ================================================================================================================================

    uo_0(:,:) = val_exp
    var_name = 'uo_0'
    CALL restget_p (rest_id, var_name, nbp_glo,   nvert, 1, itime, &
         &              .TRUE.,uo_0 , 'gather', nbp_glo, index_g)
    IF (ALL(uo_0(:,:) == val_exp)) THEN
       DO nivo=1,nvert
          IF (nivo .LE. ns) THEN
            uo_0(:,nivo) = scmax
         ELSE
            uo_0(:,nivo) = CH4atmo_CONC
         ENDIF
      ENDDO
   ENDIF

   uold2_0(:,:) = val_exp
   var_name = 'uold2_0'
   CALL restget_p (rest_id, var_name, nbp_glo,   nvert, 1, itime, &
        &              .TRUE.,uold2_0 , 'gather', nbp_glo, index_g)
   IF (ALL(uold2_0(:,:) == val_exp)) THEN
      DO nivo=1,nvert
         IF (nivo .LE. ns) THEN
            uold2_0(:,nivo) = scmax
         ELSE
            uold2_0(:,nivo) = CH4atmo_CONC
         ENDIF
      ENDDO
   ENDIF
   
   uo_wet1(:,:) = val_exp
   var_name = 'uo_wet1'
   CALL restget_p (rest_id, var_name, nbp_glo,   nvert, 1, itime, &
       &              .TRUE.,uo_wet1 , 'gather', nbp_glo, index_g)
   IF (ALL(uo_wet1(:,:) == val_exp)) THEN
      DO nivo=1,nvert
         IF (nivo .LE. ns-10) THEN
            uo_wet1(:,nivo) = scmax
         ELSE
            uo_wet1(:,nivo) = CH4atmo_CONC
         ENDIF
      ENDDO
   ENDIF
   
   uold2_wet1(:,:) = val_exp
   var_name = 'uold2_wet1'
   CALL restget_p (rest_id, var_name, nbp_glo,   nvert, 1, itime, &
        &              .TRUE.,uold2_wet1 , 'gather', nbp_glo, index_g)
   IF (ALL(uold2_wet1(:,:) == val_exp)) THEN
      DO nivo=1,nvert
         IF (nivo .LE. ns-10) THEN
            uold2_wet1(:,nivo) = scmax
         ELSE
            uold2_wet1(:,nivo) = CH4atmo_CONC
         ENDIF
      ENDDO
   ENDIF
  
   uo_wet2(:,:) = val_exp
   var_name = 'uo_wet2'
   CALL restget_p (rest_id, var_name, nbp_glo,   nvert, 1, itime, &
        &              .TRUE.,uo_wet2 , 'gather', nbp_glo, index_g)
   IF (ALL(uo_wet2(:,:) == val_exp)) THEN
      DO nivo=1,nvert
         IF (nivo .LE. ns-10) THEN
            uo_wet2(:,nivo) = scmax
         ELSE
            uo_wet2(:,nivo) = CH4atmo_CONC
         ENDIF
      ENDDO
   ENDIF

   uold2_wet2(:,:) = val_exp
   var_name = 'uold2_wet2'
   CALL restget_p (rest_id, var_name, nbp_glo,   nvert, 1, itime, &
        &              .TRUE.,uold2_wet2 , 'gather', nbp_glo, index_g)
   IF (ALL(uold2_wet2(:,:) == val_exp)) THEN
      DO nivo=1,nvert
         IF (nivo .LE. ns-10) THEN
            uold2_wet2(:,nivo) = scmax
         ELSE
            uold2_wet2(:,nivo) = CH4atmo_CONC
         ENDIF
      ENDDO
   ENDIF

   uo_wet3(:,:) = val_exp
   var_name = 'uo_wet3'
   CALL restget_p (rest_id, var_name, nbp_glo,   nvert, 1, itime, &
        &              .TRUE.,uo_wet3 , 'gather', nbp_glo, index_g)
   IF (ALL(uo_wet3(:,:) == val_exp)) THEN
      DO nivo=1,nvert
         IF (nivo .LE. ns-10) THEN
            uo_wet3(:,nivo) = scmax
         ELSE
            uo_wet3(:,nivo) = CH4atmo_CONC
         ENDIF
      ENDDO
   ENDIF
   
   uold2_wet3(:,:) = val_exp
   var_name = 'uold2_wet3'
   CALL restget_p (rest_id, var_name, nbp_glo,   nvert, 1, itime, &
        &              .TRUE.,uold2_wet3 , 'gather', nbp_glo, index_g)
   IF (ALL(uold2_wet3(:,:) == val_exp)) THEN
      DO nivo=1,nvert
         IF (nivo .LE. ns-10) THEN
            uold2_wet3(:,nivo) = scmax
         ELSE
            uold2_wet3(:,nivo) = CH4atmo_CONC
         ENDIF
      ENDDO
   ENDIF

   uo_wet4(:,:) = val_exp
   var_name = 'uo_wet4'
   CALL restget_p (rest_id, var_name, nbp_glo,   nvert, 1, itime, &
        &              .TRUE.,uo_wet4 , 'gather', nbp_glo, index_g)
   IF (ALL(uo_wet4(:,:) == val_exp)) THEN
      DO nivo=1,nvert
         IF (nivo .LE. ns-10) THEN
            uo_wet4(:,nivo) = scmax
         ELSE
            uo_wet4(:,nivo) = CH4atmo_CONC
         ENDIF
      ENDDO
   ENDIF

   uold2_wet4(:,:) = val_exp
   var_name = 'uold2_wet4'
   CALL restget_p (rest_id, var_name, nbp_glo,   nvert, 1, itime, &
        &              .TRUE.,uold2_wet4 , 'gather', nbp_glo, index_g)
   IF (ALL(uold2_wet4(:,:) == val_exp)) THEN
      DO nivo=1,nvert
         IF (nivo .LE. ns-10) THEN
            uold2_wet4(:,nivo) = scmax
         ELSE
            uold2_wet4(:,nivo) = CH4atmo_CONC
         ENDIF
      ENDDO
   ENDIF
  END SUBROUTINE stomate_wet_ch4_readrestart

!! ================================================================================================================================
!! SUBROUTINE   : stomate_wet_ch4_writerestart
!!
!>\BRIEF        Module restart write variables
!!
!! DESCRIPTION  : 
!!                
!! \n
!_ ================================================================================================================================
!!
  SUBROUTINE stomate_wet_ch4_writerestart ( rest_id, itime )

    ! 0.1 Input variables
    INTEGER(i_std), INTENT(in) :: rest_id
    INTEGER(i_std), INTENT(in) :: itime

    ! 0.4 Local variables
    CHARACTER(LEN=100)                   :: var_name !! Name of permafrost forcing file

!_ ================================================================================================================================
    var_name = 'uo_0'
    CALL restput_p (rest_id, var_name, nbp_glo, nvert, 1, itime, &
         &              uo_0, 'scatter', nbp_glo, index_g)
    
    var_name = 'uold2_0'
    CALL restput_p (rest_id, var_name, nbp_glo, nvert, 1, itime, &
         &              uold2_0, 'scatter', nbp_glo, index_g)
    
    var_name = 'uo_wet1'
    CALL restput_p (rest_id, var_name, nbp_glo, nvert, 1, itime, &
         &              uo_wet1, 'scatter', nbp_glo, index_g)
    
    var_name = 'uold2_wet1'
    CALL restput_p (rest_id, var_name, nbp_glo, nvert, 1, itime, &
         &              uold2_wet1, 'scatter', nbp_glo, index_g)
    
    var_name = 'uo_wet2'
    CALL restput_p (rest_id, var_name, nbp_glo, nvert, 1, itime, &
         &              uo_wet2, 'scatter', nbp_glo, index_g)
 
    var_name = 'uold2_wet2'
    CALL restput_p (rest_id, var_name, nbp_glo, nvert, 1, itime, &
         &              uold2_wet2, 'scatter', nbp_glo, index_g)

    var_name = 'uo_wet3'
    CALL restput_p (rest_id, var_name, nbp_glo, nvert, 1, itime, &
       &              uo_wet3, 'scatter', nbp_glo, index_g)
    
    var_name = 'uold2_wet3'
    CALL restput_p (rest_id, var_name, nbp_glo, nvert, 1, itime, &
         &              uold2_wet3, 'scatter', nbp_glo, index_g)
    
    var_name = 'uo_wet4'
    CALL restput_p (rest_id, var_name, nbp_glo, nvert, 1, itime, &
         &              uo_wet4, 'scatter', nbp_glo, index_g)
    
    var_name = 'uold2_wet4'
    CALL restput_p (rest_id, var_name, nbp_glo, nvert, 1, itime, &
         &              uold2_wet4, 'scatter', nbp_glo, index_g)

  END SUBROUTINE stomate_wet_ch4_writerestart


!! ================================================================================================================================
!! SUBROUTINE   : stomate_wet_ch4_initialize
!!
!>\BRIEF        Module initialization
!!
!! DESCRIPTION  : 
!!                
!! \n
!_ ================================================================================================================================
!!
  SUBROUTINE stomate_wet_ch4_initialize ( kjpindex, rest_id, itime )
    !
    ! 0 declarations
    !

    ! 0.1 input
    INTEGER(i_std), INTENT(in)                     :: kjpindex ! Domain size
    INTEGER(i_std), INTENT(in)                     :: rest_id
    INTEGER(i_std), INTENT(in)                     :: itime

    ! 0.2 modified fields

    ! 0.3 output

    ! 0.4 local
!_ ================================================================================================================================
   !1.4.10 wetland CH4
   !pss:+
   !appel routines pour calcul des densites de flux de CH4

        
   !Config  Key  = CH4atm_CONC
   !Config  Desc = 
   !Config If    = CH4_CALCUL
   !Config  Def  = 0.0017
   !Config  Help = 
   !               
   !Config Units = [-]
   CH4atmo_CONC=0.0017
   CALL getin_p('CH4atmo_CONC', CH4atmo_CONC)

   CH4_WTD1  = .TRUE.
   !Config  Key  = CH4_WTD1
   !Config  Desc = 
   !Config If    = CH4_CALCUL
   !Config  Def  = True
   !Config  Help = 
   !               
   !Config Units = Y/n
   CALL getin_p('CH4_WTD1', CH4_WTD1)
   !Config  Key  = CH4_WTD2
   !Config  Desc = 
   !Config If    = CH4_CALCUL
   !Config  Def  = True
   !Config  Help = 
   !               
   !Config Units = Y/n
   CH4_WTD2  = .TRUE.
   CALL getin_p('CH4_WDT2', CH4_WTD2)
   !Config  Key  = CH4_WTD3
   !Config  Desc = 
   !Config If    = CH4_CALCUL
   !Config  Def  = True
   !Config  Help = 
   !               
   !Config Units = Y/n
   CH4_WTD3  = .TRUE.
   CALL getin_p('CH4_WTD3', CH4_WTD3)
   !Config  Key  = CH4_WTD4
   !Config  Desc = 
   !Config If    = CH4_CALCUL
   !Config  Def  = True
   !Config  Help = 
   !               
   !Config Units = Y/n
   CH4_WTD4  = .TRUE.
   CALL getin_p('CH4_WTD4', CH4_WTD4)

   ! allocate
   CALL stomate_wet_ch4_init (kjpindex)

   ! read data from restart files (if any)
   CALL stomate_wet_ch4_readrestart (rest_id, itime)

  END SUBROUTINE stomate_wet_ch4_initialize

!! ================================================================================================================================
!! SUBROUTINE   : stomate_wet_ch4_init
!!
!>\BRIEF        Module variables allocation 
!!
!! DESCRIPTION  : 
!!                
!! \n
!_ ================================================================================================================================
!!
  SUBROUTINE stomate_wet_ch4_init ( kjpindex )
    !
    ! 0 declarations
    !

    ! 0.1 input

    ! Domain size
    INTEGER(i_std), INTENT(in)                                 :: kjpindex
    ! 0.2 modified fields

    ! 0.3 output

    ! 0.4 local
    INTEGER(i_std)               :: ier
!_ ================================================================================================================================

    ALLOCATE(ch4_flux_density_tot_0(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3, 'stomate_wet_ch4_initialize', 'There is an allocation variable error: ', 'ch4_flux_density_tot_0', '')
    ch4_flux_density_tot_0 = zero

    ALLOCATE(ch4_flux_density_dif_0(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3, 'stomate_wet_ch4_initialize', 'There is an allocation variable error: ', 'ch4_flux_density_dif_0', '')
    ch4_flux_density_dif_0 = zero

    ALLOCATE(ch4_flux_density_bub_0(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3, 'stomate_wet_ch4_initialize', 'There is an allocation variable error: ', 'ch4_flux_density_bub_0', '')
    ch4_flux_density_bub_0 = zero

    ALLOCATE(ch4_flux_density_pla_0(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3, 'stomate_wet_ch4_initialize', 'There is an allocation variable error: ', 'ch4_flux_density_pla_0', '')
    ch4_flux_density_pla_0 = zero

    
    ALLOCATE(ch4_flux_density_tot_wet1(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3, 'stomate_wet_ch4_initialize', 'There is an allocation variable error: ', 'ch4_flux_density_tot_wet1', '')
    ch4_flux_density_tot_wet1 = zero

    ALLOCATE(ch4_flux_density_dif_wet1(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3, 'stomate_wet_ch4_initialize', 'There is an allocation variable error: ', 'ch4_flux_density_dif_wet1', '')
    ch4_flux_density_dif_wet1 = zero

    ALLOCATE(ch4_flux_density_bub_wet1(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3, 'stomate_wet_ch4_initialize', 'There is an allocation variable error: ', 'ch4_flux_density_bub_wet1', '')
    ch4_flux_density_bub_wet1 = zero

    ALLOCATE(ch4_flux_density_pla_wet1(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3, 'stomate_wet_ch4_initialize', 'There is an allocation variable error: ', 'ch4_flux_density_pla_wet1', '')
    ch4_flux_density_pla_wet1 = zero


    ALLOCATE(ch4_flux_density_tot_wet2(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3, 'stomate_wet_ch4_initialize', 'There is an allocation variable error: ', 'ch4_flux_density_tot_wet2', '')
    ch4_flux_density_tot_wet2 = zero

    ALLOCATE(ch4_flux_density_dif_wet2(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3, 'stomate_wet_ch4_initialize', 'There is an allocation variable error: ', 'ch4_flux_density_dif_wet2', '')
    ch4_flux_density_dif_wet2 = zero

    ALLOCATE(ch4_flux_density_bub_wet2(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3, 'stomate_wet_ch4_initialize', 'There is an allocation variable error: ', 'ch4_flux_density_bub_wet2', '')
    ch4_flux_density_bub_wet2 = zero

    ALLOCATE(ch4_flux_density_pla_wet2(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3, 'stomate_wet_ch4_initialize', 'There is an allocation variable error: ', 'ch4_flux_density_pla_wet2', '')
    ch4_flux_density_pla_wet2 = zero

   
    ALLOCATE(ch4_flux_density_tot_wet3(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3, 'stomate_wet_ch4_initialize', 'There is an allocation variable error: ', 'ch4_flux_density_tot_wet3', '')
    ch4_flux_density_tot_wet3 = zero

    ALLOCATE(ch4_flux_density_dif_wet3(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3, 'stomate_wet_ch4_initialize', 'There is an allocation variable error: ', 'ch4_flux_density_dif_wet3', '')
    ch4_flux_density_dif_wet3 = zero

    ALLOCATE(ch4_flux_density_bub_wet3(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3, 'stomate_wet_ch4_initialize', 'There is an allocation variable error: ', 'ch4_flux_density_bub_wet3', '')
    ch4_flux_density_bub_wet3 = zero

    ALLOCATE(ch4_flux_density_pla_wet3(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3, 'stomate_wet_ch4_initialize', 'There is an allocation variable error: ', 'ch4_flux_density_pla_wet3', '')
    ch4_flux_density_pla_wet3 = zero

    
    ALLOCATE(ch4_flux_density_tot_wet4(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3, 'stomate_wet_ch4_initialize', 'There is an allocation variable error: ', 'ch4_flux_density_tot_wet4', '')
    ch4_flux_density_tot_wet4 = zero
    
    ALLOCATE(ch4_flux_density_dif_wet4(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3, 'stomate_wet_ch4_initialize', 'There is an allocation variable error: ', 'ch4_flux_density_dif_wet4', '')
    ch4_flux_density_dif_wet4 = zero
   
    ALLOCATE(ch4_flux_density_bub_wet4(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3, 'stomate_wet_ch4_initialize', 'There is an allocation variable error: ', 'ch4_flux_density_bub_wet4', '')
    ch4_flux_density_bub_wet4 = zero
  
    ALLOCATE(ch4_flux_density_pla_wet4(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3, 'stomate_wet_ch4_initialize', 'There is an allocation variable error: ', 'ch4_flux_density_pla_wet4', '')
    ch4_flux_density_pla_wet4 = zero
 

    ALLOCATE(uo_0(kjpindex,nvert),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3, 'stomate_wet_ch4_initialize', 'There is an allocation variable error: ', 'uo_0', '')
    uo_0 = zero

    ALLOCATE(uold2_0(kjpindex,nvert),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3, 'stomate_wet_ch4_initialize', 'There is an allocation variable error: ', 'uold2_0', '')
    uold2_0 = zero

    ALLOCATE(uo_wet1(kjpindex,nvert),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3, 'stomate_wet_ch4_initialize', 'There is an allocation variable error: ', 'uo_wet1', '')
    uo_wet1 = zero

    ALLOCATE(uold2_wet1(kjpindex,nvert),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3, 'stomate_wet_ch4_initialize', 'There is an allocation variable error: ', 'uold2_wet1', '')
    uold2_wet1 = zero

    ALLOCATE(uo_wet2(kjpindex,nvert),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3, 'stomate_wet_ch4_initialize', 'There is an allocation variable error: ', 'uo_wet2', '')
    uo_wet2 = zero

    ALLOCATE(uold2_wet2(kjpindex,nvert),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3, 'stomate_wet_ch4_initialize', 'There is an allocation variable error: ', 'uold2_wet2', '')
    uold2_wet2 = zero

    ALLOCATE(uo_wet3(kjpindex,nvert),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3, 'stomate_wet_ch4_initialize', 'There is an allocation variable error: ', 'uo_wet3', '')
    uo_wet3 = zero

    ALLOCATE(uold2_wet3(kjpindex,nvert),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3, 'stomate_wet_ch4_initialize', 'There is an allocation variable error: ', 'uold2_wet3', '')
    uold2_wet3 = zero

    ALLOCATE(uo_wet4(kjpindex,nvert),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3, 'stomate_wet_ch4_initialize', 'There is an allocation variable error: ', 'uo_wet4', '')
    uo_wet4 = zero

    ALLOCATE(uold2_wet4(kjpindex,nvert),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3, 'stomate_wet_ch4_initialize', 'There is an allocation variable error: ', 'uold2_wet4', '')
    uold2_wet4 = zero

  END SUBROUTINE stomate_wet_ch4_init

!! ================================================================================================================================
!! SUBROUTINE   : stomate_wet_ch4_main
!!
!>\BRIEF        ch4 main body
!!
!! DESCRIPTION  : 
!!                
!! \n
!_ ================================================================================================================================
!!
  SUBROUTINE stomate_wet_ch4_main ( kjpindex,       stempdiag,      tsurf_daily, tsurf_year,    &
                                    veget_cov_max,  veget,          lai, carbon,                &
                                    carbon_surf,    itime,          hist_id_stomate, hori_index) 
    !
    ! 0 declarations
    !

    ! 0.1 input

    ! Domain size
    INTEGER(i_std), INTENT(in)                                 :: kjpindex 
    REAL(r_std),DIMENSION (:,:), INTENT (in)               :: stempdiag   ! kjpindex, nslm
    ! temperature (K) at the surface
    REAL(r_std), DIMENSION(:), INTENT(in)                  :: tsurf_daily ! kjpindex

    ! temperature (K) at the surface
    REAL(r_std), DIMENSION(:), INTENT(in)                  :: tsurf_year ! kjpindex

    ! "maximal" coverage fraction of a PFT (LAI -> infinity) on nat/agri ground
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: veget_cov_max    ! kjpindex, nvm
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: veget            ! kjpindex, nvm
    REAL(r_std), DIMENSION(:,:),INTENT(in)  :: lai               !! Leaf area inex @tex $(m^2 m^{-2})$ @endtex
    !! Carbon pool integrated to over surface soils: active, slow, or passive ( kjpindex, ncarb, nvm )
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)        :: carbon_surf 
    !! Soil carbon pools per ground area: active, slow, or ( kjpindex, ncarb, nvm ) 
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)        :: carbon 
    INTEGER(i_std), INTENT(in)                                      :: itime 
    INTEGER(i_std), INTENT(in)                                      :: hist_id_stomate
    INTEGER(i_std), DIMENSION(:), INTENT(in)                        :: hori_index     !! Move to Horizontal indices

    ! 0.2 modified fields

    ! 0.3 output

    ! 0.4 local

!_ ================================================================================================================================

    IF (.NOT. CH4_calcul) CALL ipslerr_p(3, 'stomate_wet_ch4_main', & 
        'CH4 module needs to be enabled to run this subroutine', 'Set CH4_CALCUL=y', &
        'in your run.def file')


    !routine pour densite de flux d un wetland ou WTD = 0
    CALL ch4_wet_flux_density_0 (kjpindex,stempdiag,tsurf_daily,tsurf_year,veget_cov_max,veget,&
         & carbon,lai,uo_0,uold2_0, ch4_flux_density_tot_0, ch4_flux_density_dif_0,&
         & ch4_flux_density_bub_0,ch4_flux_density_pla_0, CH4atmo_CONC)  


    IF (CH4_WTD1) THEN
       !routine calcule densite de flux d un wetland ou WTD = pwater_wet1 (cf.stomate_cste_wetlands.f90) 
       CALL ch4_wet_flux_density_wet (kjpindex,stempdiag,tsurf_daily,tsurf_year,veget_cov_max,veget,&
            & carbon_surf,lai,uo_wet1,uold2_wet1,ch4_flux_density_tot_wet1, ch4_flux_density_dif_wet1, &
            & ch4_flux_density_bub_wet1,ch4_flux_density_pla_wet1, CH4atmo_CONC, pwater_wet1)  
    ENDIF
         
      IF (CH4_WTD2) THEN
         !routine calcule densite de flux d un wetland ou WTD = pwater_wet2 (cf.stomate_cste_wetlands.f90) 
         CALL ch4_wet_flux_density_wet (kjpindex,stempdiag,tsurf_daily,tsurf_year,veget_cov_max,veget,&
              & carbon_surf,lai,uo_wet2,uold2_wet2,ch4_flux_density_tot_wet2, ch4_flux_density_dif_wet2, &
              & ch4_flux_density_bub_wet2,ch4_flux_density_pla_wet2, CH4atmo_CONC, pwater_wet2)  
      ENDIF

      IF (CH4_WTD3) THEN
         !routine calcule densite de flux d un wetland ou WTD = pwater_wet3 (cf.stomate_cste_wetlands.f90) 
         CALL ch4_wet_flux_density_wet (kjpindex,stempdiag,tsurf_daily,tsurf_year,veget_cov_max,veget,&
              & carbon_surf,lai,uo_wet3,uold2_wet3,ch4_flux_density_tot_wet3, ch4_flux_density_dif_wet3, &
              & ch4_flux_density_bub_wet3,ch4_flux_density_pla_wet3, CH4atmo_CONC, pwater_wet3)  
      ENDIF

      IF (CH4_WTD4) THEN
         !routine calcule densite de flux d un wetland ou WTD = pwater_wet4 (cf.stomate_cste_wetlands.f90) 
         CALL ch4_wet_flux_density_wet (kjpindex,stempdiag,tsurf_daily,tsurf_year,veget_cov_max,veget,&
            & carbon_surf,lai,uo_wet4,uold2_wet4,ch4_flux_density_tot_wet4, ch4_flux_density_dif_wet4, &
            & ch4_flux_density_bub_wet4,ch4_flux_density_pla_wet4, CH4atmo_CONC, pwater_wet4)  
      ENDIF
      
!!!! Wetland CH4 methane
!pss:+
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_TOT_0', itime, &
                    ch4_flux_density_tot_0, kjpindex, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_DIF_0', itime, &
                    ch4_flux_density_dif_0, kjpindex, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_BUB_0', itime, &
                    ch4_flux_density_bub_0, kjpindex, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_PLA_0', itime, &
                    ch4_flux_density_pla_0, kjpindex, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_TOT_wet1', itime, &
                    ch4_flux_density_tot_wet1, kjpindex, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_DIF_wet1', itime, &
                    ch4_flux_density_dif_wet1, kjpindex, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_BUB_wet1', itime, &
                    ch4_flux_density_bub_wet1, kjpindex, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_PLA_wet1', itime, &
                    ch4_flux_density_pla_wet1, kjpindex, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_TOT_wet2', itime, &
                    ch4_flux_density_tot_wet2, kjpindex, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_DIF_wet2', itime, &
                    ch4_flux_density_dif_wet2, kjpindex, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_BUB_wet2', itime, &
                    ch4_flux_density_bub_wet2, kjpindex, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_PLA_wet2', itime, &
                    ch4_flux_density_pla_wet2, kjpindex, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_TOT_wet3', itime, &
                    ch4_flux_density_tot_wet3, kjpindex, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_DIF_wet3', itime, &
                    ch4_flux_density_dif_wet3, kjpindex, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_BUB_wet3', itime, &
                    ch4_flux_density_bub_wet3, kjpindex, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_PLA_wet3', itime, &
                    ch4_flux_density_pla_wet3, kjpindex, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_TOT_wet4', itime, &
                    ch4_flux_density_tot_wet4, kjpindex, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_DIF_wet4', itime, &
                    ch4_flux_density_dif_wet4, kjpindex, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_BUB_wet4', itime, &
                    ch4_flux_density_bub_wet4, kjpindex, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_PLA_wet4', itime, &
                    ch4_flux_density_pla_wet4, kjpindex, hori_index)

  END SUBROUTINE stomate_wet_ch4_main

!! ================================================================================================================================
!! SUBROUTINE   : stomate_wet_ch4_clear
!!
!>\BRIEF        Module variables deallocation
!!
!! DESCRIPTION  : 
!!                
!! \n
!_ ================================================================================================================================
!!
  SUBROUTINE stomate_wet_ch4_clear ( )
    !
    ! 0 declarations
    !

    ! 0.1 input

    ! 0.2 modified fields

    ! 0.3 output

    ! density flux of methane calculated for entire pixel (gCH4/dt/m**2)

    ! 0.4 local


    IF (ALLOCATED(uo_0)) DEALLOCATE(uo_0)
    IF (ALLOCATED(uold2_0)) DEALLOCATE(uold2_0)
    IF (ALLOCATED(uo_wet1)) DEALLOCATE(uo_wet1)
    IF (ALLOCATED(uold2_wet1)) DEALLOCATE(uold2_wet1)
    IF (ALLOCATED(uo_wet2)) DEALLOCATE(uo_wet2)
    IF (ALLOCATED(uold2_wet2)) DEALLOCATE(uold2_wet2)
    IF (ALLOCATED(uo_wet3)) DEALLOCATE(uo_wet3)
    IF (ALLOCATED(uold2_wet3)) DEALLOCATE(uold2_wet3)
    IF (ALLOCATED(uo_wet4)) DEALLOCATE(uo_wet4)
    IF (ALLOCATED(uold2_wet4)) DEALLOCATE(uold2_wet4)

    IF (ALLOCATED(ch4_flux_density_tot_0)) DEALLOCATE(ch4_flux_density_tot_0)
    IF (ALLOCATED(ch4_flux_density_dif_0)) DEALLOCATE(ch4_flux_density_dif_0)
    IF (ALLOCATED(ch4_flux_density_bub_0)) DEALLOCATE(ch4_flux_density_bub_0)
    IF (ALLOCATED(ch4_flux_density_pla_0)) DEALLOCATE(ch4_flux_density_pla_0)
    IF (ALLOCATED(ch4_flux_density_tot_wet1)) DEALLOCATE(ch4_flux_density_tot_wet1)
    IF (ALLOCATED(ch4_flux_density_dif_wet1)) DEALLOCATE(ch4_flux_density_dif_wet1)
    IF (ALLOCATED(ch4_flux_density_bub_wet1)) DEALLOCATE(ch4_flux_density_bub_wet1)
    IF (ALLOCATED(ch4_flux_density_pla_wet1)) DEALLOCATE(ch4_flux_density_pla_wet1)
    IF (ALLOCATED(ch4_flux_density_tot_wet2)) DEALLOCATE(ch4_flux_density_tot_wet2)
    IF (ALLOCATED(ch4_flux_density_dif_wet2)) DEALLOCATE(ch4_flux_density_dif_wet2)
    IF (ALLOCATED(ch4_flux_density_bub_wet2)) DEALLOCATE(ch4_flux_density_bub_wet2)
    IF (ALLOCATED(ch4_flux_density_pla_wet2)) DEALLOCATE(ch4_flux_density_pla_wet2)
    IF (ALLOCATED(ch4_flux_density_tot_wet3)) DEALLOCATE(ch4_flux_density_tot_wet3)
    IF (ALLOCATED(ch4_flux_density_dif_wet3)) DEALLOCATE(ch4_flux_density_dif_wet3)
    IF (ALLOCATED(ch4_flux_density_bub_wet3)) DEALLOCATE(ch4_flux_density_bub_wet3)
    IF (ALLOCATED(ch4_flux_density_pla_wet3)) DEALLOCATE(ch4_flux_density_pla_wet3)
    IF (ALLOCATED(ch4_flux_density_tot_wet4)) DEALLOCATE(ch4_flux_density_tot_wet4)
    IF (ALLOCATED(ch4_flux_density_dif_wet4)) DEALLOCATE(ch4_flux_density_dif_wet4)
    IF (ALLOCATED(ch4_flux_density_bub_wet4)) DEALLOCATE(ch4_flux_density_bub_wet4)
    IF (ALLOCATED(ch4_flux_density_pla_wet4)) DEALLOCATE(ch4_flux_density_pla_wet4)

    CALL ch4_wet_flux_density_clear_0
    CALL ch4_wet_flux_density_clear_wet

  END SUBROUTINE stomate_wet_ch4_clear


!! ================================================================================================================================
!! SUBROUTINE   : stomate_wet_ch4_finalize
!!
!>\BRIEF       Write data to restart file 
!!
!! DESCRIPTION  : 
!!                
!! \n
!_ ================================================================================================================================
!!
  SUBROUTINE stomate_wet_ch4_finalize ( rest_id, itime )

    ! 0.1 Input variables
    INTEGER(i_std), INTENT(in) :: rest_id
    INTEGER(i_std), INTENT(in) :: itime

!_ =====================================================================================================================

    !
    ! 0 declarations
    !

    ! 0.1 input

    ! 0.2 modified fields

    ! 0.3 output

    ! density flux of methane calculated for entire pixel (gCH4/dt/m**2)

    ! 0.4 local

    CALL stomate_wet_ch4_writerestart( rest_id, itime)

  END SUBROUTINE stomate_wet_ch4_finalize


!! ================================================================================================================================
!! SUBROUTINE   : stomate_wet_ch4_histdef
!!
!>\BRIEF        Define IOIPSL history output
!!
!! DESCRIPTION  :  
!!                
!! \n
!_ ================================================================================================================================
!!
  SUBROUTINE stomate_wet_ch4_histdef (iim, jjm, dt, hist_hori_id, hist_id, ave, &
                                      hist_id_stomate )
    !
    ! 0 declarations
    !

    ! 0.1 input
    INTEGER(i_std), INTENT(in)          :: iim, jjm  !! Size in x and y of the data to be handeled
    
    REAL(r_std),INTENT(in)              :: hist_id            !- Time step of history file (s)
    
    INTEGER(i_std),INTENT(in)           :: hist_hori_id !- id horizontal grid
    
    REAL(r_std),INTENT(in)              :: dt !- Time step of STOMATE (seconds)

    CHARACTER(LEN=40), DIMENSION(:), INTENT(in) ::  ave

    INTEGER(i_std),INTENT(in)           :: hist_id_stomate !- history stomate id
    ! 0.2 modified fields

    ! 0.3 output

    ! 0.4 local
!_ ================================================================================================================================

    CALL histdef (hist_id_stomate, &
         &               TRIM("CH4_FLUX_TOT_0      "), &
         &               TRIM("flux density tot of CH4 by wetlands               "), &
         &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(5), dt, hist_id)
    
    CALL histdef (hist_id_stomate, &
         &               TRIM("CH4_FLUX_DIF_0      "), &
         &               TRIM("flux density dif of CH4 by wetlands               "), &
         &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(5), dt, hist_id)
    
    CALL histdef (hist_id_stomate, &
         &               TRIM("CH4_FLUX_BUB_0      "), &
         &               TRIM("flux density bub of CH4 by wetlands               "), &
         &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(5), dt, hist_id)
    
    CALL histdef (hist_id_stomate, &
         &               TRIM("CH4_FLUX_PLA_0      "), &
         &               TRIM("flux density pla of CH4 by wetlands               "), &
         &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(5), dt, hist_id)
    
    !!pour wetland avec WTD = -x1
    CALL histdef (hist_id_stomate, &
         &               TRIM("CH4_FLUX_TOT_wet1    "), &
         &               TRIM("flux density tot of CH4 by wetlands               "), &
         &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(5), dt, hist_id)
    
   CALL histdef (hist_id_stomate, &
        &               TRIM("CH4_FLUX_DIF_wet1    "), &
        &               TRIM("flux density dif of CH4 by wetlands               "), &
        &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
        &               1,1,1, -99,32, ave(5), dt, hist_id)
   
   CALL histdef (hist_id_stomate, &
        &               TRIM("CH4_FLUX_BUB_wet1    "), &
        &               TRIM("flux density bub of CH4 by wetlands               "), &
        &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
        &               1,1,1, -99,32, ave(5), dt, hist_id)
   
   CALL histdef (hist_id_stomate, &
        &               TRIM("CH4_FLUX_PLA_wet1    "), &
        &               TRIM("flux density pla of CH4 by wetlands               "), &
        &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
        &               1,1,1, -99,32, ave(5), dt, hist_id)
   
   !pour wetland avc WTD = -x2
   CALL histdef (hist_id_stomate, &
        &               TRIM("CH4_FLUX_TOT_wet2    "), &
        &               TRIM("flux density tot of CH4 by wetlands               "), &
        &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
        &               1,1,1, -99,32, ave(5), dt, hist_id)
   
   CALL histdef (hist_id_stomate, &
        &               TRIM("CH4_FLUX_DIF_wet2    "), &
        &               TRIM("flux density dif of CH4 by wetlands               "), &
        &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
        &               1,1,1, -99,32, ave(5), dt, hist_id)
   
   CALL histdef (hist_id_stomate, &
        &               TRIM("CH4_FLUX_BUB_wet2    "), &
        &               TRIM("flux density bub of CH4 by wetlands               "), &
        &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
        &               1,1,1, -99,32, ave(5), dt, hist_id)
   
   CALL histdef (hist_id_stomate, &
        &               TRIM("CH4_FLUX_PLA_wet2    "), &
        &               TRIM("flux density pla of CH4 by wetlands               "), &
        &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
        &               1,1,1, -99,32, ave(5), dt, hist_id)

   !pour wetland avec WTD = -x3
   CALL histdef (hist_id_stomate, &
        &               TRIM("CH4_FLUX_TOT_wet3    "), &
        &               TRIM("flux density tot of CH4 by wetlands               "), &
        &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
        &               1,1,1, -99,32, ave(5), dt, hist_id)
   
   CALL histdef (hist_id_stomate, &
        &               TRIM("CH4_FLUX_DIF_wet3    "), &
        &               TRIM("flux density dif of CH4 by wetlands               "), &
        &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
        &               1,1,1, -99,32, ave(5), dt, hist_id)
   
   CALL histdef (hist_id_stomate, &
        &               TRIM("CH4_FLUX_BUB_wet3    "), &
        &               TRIM("flux density bub of CH4 by wetlands               "), &
        &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
        &               1,1,1, -99,32, ave(5), dt, hist_id)
   
   CALL histdef (hist_id_stomate, &
        &               TRIM("CH4_FLUX_PLA_wet3    "), &
        &               TRIM("flux density pla of CH4 by wetlands               "), &
        &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
        &               1,1,1, -99,32, ave(5), dt, hist_id)
   
   !wetland avc WTD = -x4
   CALL histdef (hist_id_stomate, &
        &               TRIM("CH4_FLUX_TOT_wet4    "), &
        &               TRIM("flux density tot of CH4 by wetlands               "), &
        &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
        &               1,1,1, -99,32, ave(5), dt, hist_id)
   
   CALL histdef (hist_id_stomate, &
        &               TRIM("CH4_FLUX_DIF_wet4    "), &
        &               TRIM("flux density dif of CH4 by wetlands               "), &
        &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
        &               1,1,1, -99,32, ave(5), dt, hist_id)
   
   CALL histdef (hist_id_stomate, &
        &               TRIM("CH4_FLUX_BUB_wet4    "), &
        &               TRIM("flux density bub of CH4 by wetlands               "), &
        &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
        &               1,1,1, -99,32, ave(5), dt, hist_id)
   
   CALL histdef (hist_id_stomate, &
        &               TRIM("CH4_FLUX_PLA_wet4    "), &
        &               TRIM("flux density pla of CH4 by wetlands               "), &
        &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
        &               1,1,1, -99,32, ave(5), dt, hist_id)

  END SUBROUTINE stomate_wet_ch4_histdef
  

END MODULE stomate_wet_ch4




















