!< $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/perso/albert.jornet/OMCROP/src_stomate/stomate_io.f90 $ 
!< $Date: 2016-05-09 17:27:44 +0200 (Mon, 09 May 2016) $
!< $Author: xuhui.wang $
!< $Revision: 3419 $
! IPSL (2006)
!  This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
MODULE sticslai_io
  !---------------------------------------------------------------------
  !- Not all variables saved in the start files are absolutely necessary.
  !- However, Sechiba's and Stomate's PFTs are not necessarily identical,
  !- and for that case this information needs to be saved.
  !---------------------------------------------------------------------
  USE defprec
  USE stomate_data
  USE constantes
  USE constantes_soil
  USE mod_orchidee_para
  USE ioipsl
  USE ioipsl_para 
  !-
  IMPLICIT NONE
  !-
  PRIVATE
  PUBLIC sticslai_io_readstart, sticslai_io_writestart
  !-
  ! reference temperature (K)
  !-
  REAL(r_std),ALLOCATABLE,DIMENSION(:),SAVE :: trefe
!$OMP THREADPRIVATE(trefe)
  !-
CONTAINS
  SUBROUTINE sticslai_io_writestart (npts, & 
       &  f_crop_recycle,in_cycle, f_sen_lai, st2m_max_daily, wut_cm_daily, wus_cm_daily, evapot_daily, pdbiomass, pdmasec, &
       &  masecveg, masec, dltams, gdh_daily, phoi, onarretesomcourdrp,  &
       &  nsendltams, nsendltai, nsenpfeuilverte, nsendurvie, nsenndurvie, densiteequiv, &
       &  nplt, tursla, ssla, pfeuilverte, bsenlai, &
       &  zrac, nrec, nlan, tcult, udevair, udevcult, ndrp, rfvi, nlev, nger, etatvernal, &
       &  caljvc, rfpi, upvt, utp, somcour, somcourdrp, somcourutp, tdevelop, somtemp, &
       &  somcourfauche, stpltger, R_stamflax, R_stlaxsen, R_stsenlan, stlevflo, nflo, &
       &  R_stlevdrp, R_stflodrp, R_stdrpmat, nmat, nlax, nrecbutoir, group, ndebdes, R_stdrpdes, densite, &
       &  densitelev, coeflev, densiteger, somelong, somger, humectation, nbjhumec, &
       &  somtemphumec, stpltlev, namf, stmatrec, tustress, lai, somfeuille, pdlai, &
       &  nbfeuille, reajust, ulai, pdulai, efdensite, tempeff, nstopfeuille, deltai, vmax, nsen, &
       &  laisen, pdlaisen, dltaisenat, nsencour, dltamsen, dltaisen, fgellev, &
       &  gelee, fstressgel, R_stlevamf, dernier_n, durvieI, durvie, ndebsen, somsenreste, &
       &  humrel, swfac, turfac, senfac,mafeuiljaune, msneojaune, &
       &  v_dltams, fgelflo, pdircarb, ircarb, nbgrains, pgrain, vitmoy, nbgraingel, pgraingel, &
       &  dltags, ftempremp, magrain, pdmagrain, nbj0remp, pdsfruittot, repracmax, repracmin, &
       &  kreprac, somtemprac, urac, reprac, nstoprac, c_reserve, c_leafb, gslen, drylen,  &
       &  nboxmax, box_ndays, box_lai, box_lairem, box_tdev, box_biom, box_biomrem, box_durage, box_somsenbase, &
       &  cyc_num, cyc_num_tot,rot_cmd_store, plantdate, plantdate_now )


 ! 0.1 Inputs Variables

    ! Domain size
    INTEGER(i_std),INTENT(in) :: npts
    LOGICAL, DIMENSION(npts, nvm), INTENT(IN)       :: in_cycle 
    LOGICAL, DIMENSION(npts, nvm), INTENT(IN)       :: f_sen_lai
    LOGICAL, DIMENSION(npts, nvm), INTENT(IN)       :: f_crop_recycle 
    ! daily maximum 2 meter temperatures (K)
    REAL(r_std), DIMENSION(npts), INTENT(IN)      :: st2m_max_daily
    ! daily value of soil temperature at the resolution of 1 cm, the second dimension is 3
    ! the three layers around sowing layer
    REAL(r_std), DIMENSION(npts, nvm, 3), INTENT(IN)    :: wut_cm_daily
    ! daily mean value of soil relative humidity at the resolution of 1 cm, the second dimension is 3
    ! the three layers around sowing layer
    REAL(r_std), DIMENSION(npts, nvm, 3), INTENT(IN)    :: wus_cm_daily
    ! daily potential evapotranspiration 
    REAL(r_std), DIMENSION(npts), INTENT(IN)      :: evapot_daily
    ! biomass of previous day, t/ha
    REAL(r_std), DIMENSION(npts, nvm), INTENT(IN)      :: pdbiomass
    REAL(r_std), DIMENSION(npts, nvm), INTENT(IN)      :: pdmasec
    ! vegetative biomass
    REAL(r_std), DIMENSION(npts, nvm), INTENT(IN)      :: masecveg   
    ! aboveground dry matter (t ha-1) 
    REAL(r_std), DIMENSION(npts, nvm), INTENT(IN)      :: masec
    ! growth rate of plant, it means the delta total biomass increment (t ha-1)
    REAL(r_std), DIMENSION(npts, nvm), INTENT(IN)      :: dltams
    ! daily gdh calculated according to halfhourly temperature // transmitted from stomate.f90 gdh_daily
    REAL(r_std), DIMENSION(npts, nvm), INTENT(IN)         :: gdh_daily 
    ! Photoperiod // hours
    REAL(r_std),  DIMENSION(npts), INTENT(IN)                             :: phoi 

    !  
    LOGICAL, DIMENSION(npts, nvm), INTENT(IN)           :: onarretesomcourdrp 
    !INTEGER(i_std), DIMENSION(nvm), INTENT(IN)                           :: codeulaivernal 
    !INTEGER(i_std), DIMENSION(npts, nvm), INTENT(IN)                :: nlevobs    ! the following variables ended with obs are only used for forcing simulation.  
    !INTEGER(i_std), DIMENSION(npts, nvm), INTENT(IN)                :: namfobs    ! the initial value should be always 999
    !INTEGER(i_std), DIMENSION(npts, nvm), INTENT(IN)                :: nfloobs  
    !INTEGER(i_std), DIMENSION(npts, nvm), INTENT(IN)                :: nlanobs  
    !INTEGER(i_std), DIMENSION(npts, nvm), INTENT(IN)                :: nlaxobs  
    !INTEGER(i_std), DIMENSION(npts, nvm), INTENT(IN)                :: nmatobs  
    !INTEGER(i_std), DIMENSION(npts, nvm), INTENT(IN)                :: nrecobs  
    !INTEGER(i_std), DIMENSION(npts, nvm), INTENT(IN)                :: nsenobs  
    !INTEGER(i_std), DIMENSION(npts, nvm), INTENT(IN)                :: ndrpobs  

    ! LAIdev SPECIFIC 
    REAL(r_std), DIMENSION(npts, nvm), INTENT(IN)        :: nsendltams
    REAL(r_std), DIMENSION(npts, nvm), INTENT(IN)        :: nsendltai
    REAL(r_std), DIMENSION(npts, nvm), INTENT(IN)        :: nsenpfeuilverte
    REAL(r_std), DIMENSION(npts, nvm), INTENT(IN)        :: nsendurvie
    REAL(r_std), DIMENSION(npts, nvm), INTENT(IN)        :: nsenndurvie
    REAL(r_std), DIMENSION(npts, nvm), INTENT(IN)        :: densiteequiv
    INTEGER(i_std), DIMENSION(npts, nvm), INTENT(IN)        :: nplt
    REAL(r_std), DIMENSION(npts, nvm), INTENT(IN)        :: tursla
    REAL(r_std), DIMENSION(npts, nvm), INTENT(IN)        :: ssla
    REAL(r_std), DIMENSION(npts, nvm), INTENT(IN)        :: pfeuilverte
    REAL(r_std), DIMENSION(npts, nvm), INTENT(IN)        :: bsenlai
    
    ! variables are involved in DEVELOPMENT

    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: zrac
    INTEGER(i_std), DIMENSION(npts, nvm), INTENT(IN)        :: nrec
    INTEGER(i_std), DIMENSION(npts, nvm), INTENT(IN)        :: nlan
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: tcult
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: udevair
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: udevcult
    INTEGER(i_std),    DIMENSION(npts, nvm), INTENT(IN)        :: ndrp
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: rfvi
    INTEGER(i_std),    DIMENSION(npts, nvm), INTENT(IN)        :: nlev
    INTEGER(i_std),    DIMENSION(npts, nvm), INTENT(IN)        :: nger
    logical,    DIMENSION(npts, nvm), INTENT(IN)        :: etatvernal
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: caljvc
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: rfpi
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: upvt
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: utp
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: somcour
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: somcourdrp
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: somcourutp
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: tdevelop
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: somtemp
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: somcourfauche
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: stpltger
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: R_stamflax
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: R_stlaxsen
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: R_stsenlan
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: stlevflo
    INTEGER(i_std),    DIMENSION(npts, nvm), INTENT(IN)        :: nflo
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: R_stlevdrp
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: R_stflodrp
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: R_stdrpmat
    INTEGER(i_std),    DIMENSION(npts, nvm), INTENT(IN)        :: nmat
    INTEGER(i_std),    DIMENSION(npts, nvm), INTENT(IN)        :: nlax
    INTEGER(i_std),    DIMENSION(npts, nvm), INTENT(IN)        :: nrecbutoir
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: group
    INTEGER(i_std),    DIMENSION(npts, nvm), INTENT(IN)        :: ndebdes
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: R_stdrpdes
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: densite
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: densitelev
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: coeflev
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: densiteger
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: somelong
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: somger
    logical,    DIMENSION(npts, nvm), INTENT(IN)        :: humectation
    INTEGER(i_std),    DIMENSION(npts, nvm), INTENT(IN)        :: nbjhumec
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: somtemphumec
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: stpltlev
    INTEGER(i_std),    DIMENSION(npts, nvm), INTENT(IN)        :: namf
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: stmatrec
 
    ! these variables are involved in Lai_calculation
     
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: tustress
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: lai
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: somfeuille
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: pdlai
    INTEGER(i_std),    DIMENSION(npts, nvm), INTENT(IN)        :: nbfeuille
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: reajust
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: ulai
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: pdulai
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: efdensite
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: tempeff
    INTEGER(i_std),    DIMENSION(npts, nvm), INTENT(IN)        :: nstopfeuille
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: deltai
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: vmax
    INTEGER(i_std),    DIMENSION(npts, nvm), INTENT(IN)        :: nsen
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: laisen
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: pdlaisen
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)        :: dltaisenat

    ! these variables are involved in the LAIsenescence

    INTEGER(i_std),    DIMENSION(npts, nvm), INTENT(IN)      :: nsencour
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)      :: dltamsen
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)      :: dltaisen
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)      :: fgellev
    logical,    DIMENSION(npts, nvm), INTENT(IN)      :: gelee
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)      :: fstressgel
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)      :: R_stlevamf
    INTEGER(i_std),    DIMENSION(npts, nvm), INTENT(IN)      :: dernier_n
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)      :: durvieI
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)      :: durvie
    INTEGER(i_std),    DIMENSION(npts, nvm), INTENT(IN)      :: ndebsen
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)      :: somsenreste
    INTEGER(i_std), INTENT(IN)     :: nboxmax
    INTEGER(i_std),   DIMENSION(npts, nvm, nboxmax), INTENT(IN) :: box_ndays
    REAL(r_std),   DIMENSION(npts, nvm) :: boxtemp
    REAL(r_std),   DIMENSION(npts, nvm, nboxmax), INTENT(IN) :: box_lai
    REAL(r_std),   DIMENSION(npts, nvm, nboxmax), INTENT(IN) :: box_lairem
    REAL(r_std),   DIMENSION(npts, nvm, nboxmax), INTENT(IN) :: box_tdev
    REAL(r_std),   DIMENSION(npts, nvm, nboxmax), INTENT(IN) :: box_biom
    REAL(r_std),   DIMENSION(npts, nvm, nboxmax), INTENT(IN) :: box_biomrem
    REAL(r_std),   DIMENSION(npts, nvm, nboxmax), INTENT(IN) :: box_durage
    REAL(r_std),   DIMENSION(npts, nvm, nboxmax), INTENT(IN) :: box_somsenbase

    ! these variables are involved in STRESS calculation
    
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)      :: humrel
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)      :: swfac
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)      :: turfac
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)      :: senfac

    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)      :: mafeuiljaune 
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(IN)      :: msneojaune
    ! these variables are involved in the CARBON ALLOCATION calculation

    ! grain related   
    REAL(r_std),    DIMENSION(npts, nvm, vlength)      ,INTENT(IN)       :: v_dltams
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(IN)       :: fgelflo
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(IN)       :: pdircarb
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(IN)       :: ircarb
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(IN)       :: nbgrains
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(IN)       :: pgrain
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(IN)       :: vitmoy
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(IN)       :: nbgraingel
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(IN)       :: pgraingel
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(IN)       :: dltags
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(IN)       :: ftempremp
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(IN)       :: magrain
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(IN)       :: pdmagrain
    INTEGER(i_std), DIMENSION(npts, nvm)      ,INTENT(IN)       :: nbj0remp 
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(IN)       :: pdsfruittot

    ! reprac related

    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(IN)       :: repracmax
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(IN)       :: repracmin
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(IN)       :: kreprac
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(IN)       :: somtemprac
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(IN)       :: urac
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(IN)       :: reprac
    INTEGER(i_std), DIMENSION(npts, nvm)      ,INTENT(IN)       :: nstoprac 

    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(IN)       :: c_reserve
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(IN)       :: c_leafb
    INTEGER(i_std), DIMENSION(npts, nvm)      ,INTENT(IN)       :: gslen 
    INTEGER(i_std), DIMENSION(npts, nvm)      ,INTENT(IN)       :: drylen 
    INTEGER(i_std), DIMENSION(npts, nvm)      ,INTENT(IN)       :: cyc_num
    INTEGER(i_std), DIMENSION(npts)           ,INTENT(IN)       :: cyc_num_tot
    INTEGER(i_std), DIMENSION(npts, rot_cmd_max, cyc_rot_max) ,INTENT(IN) :: rot_cmd_store
    INTEGER(i_std), DIMENSION(npts, nvm, cyc_rot_max) ,INTENT(IN) :: plantdate
    INTEGER(i_std), DIMENSION(npts, nvm) ,INTENT(IN) :: plantdate_now
!! some temporary variable in type of real for restput
    REAL(r_std), DIMENSION(npts, nvm)                           :: cyc_num_real
    REAL(r_std), DIMENSION(npts)                                :: cyc_num_tot_real
    REAL(r_std), DIMENSION(npts, rot_cmd_max, cyc_rot_max)      :: rot_cmd_store_real
    REAL(r_std), DIMENSION(npts, nvm, cyc_rot_max)              :: plantdate_real
    REAL(r_std), DIMENSION(npts, nvm)              :: plantdate_now_real

 ! 0.4 Local variables

    ! STICS--local
    REAL(r_std), DIMENSION(npts, nvm)                                 :: in_cycle_real        
    REAL(r_std), DIMENSION(npts, nvm)                                 :: f_sen_lai_real        
    REAL(r_std), DIMENSION(npts, nvm)                                 :: f_crop_recycle_real        
    REAL(r_std), DIMENSION(npts, nvm)                                 :: onarretesomcourdrp_real        
    REAL(r_std), DIMENSION(npts, nvm)                                 :: humectation_real        

    !REAL(r_std), DIMENSION(nvm)                                       :: codeulaivernal_real        
    
    !REAL(r_std), DIMENSION(npts, nvm)                                 :: nlevobs_real
    !REAL(r_std), DIMENSION(npts, nvm)                                 :: namfobs_real
    !REAL(r_std), DIMENSION(npts, nvm)                                 :: nfloobs_real  
    !REAL(r_std), DIMENSION(npts, nvm)                                 :: nlanobs_real  
    !REAL(r_std), DIMENSION(npts, nvm)                                 :: nlaxobs_real  
    !REAL(r_std), DIMENSION(npts, nvm)                                 :: nmatobs_real  
    !REAL(r_std), DIMENSION(npts, nvm)                                 :: nrecobs_real  
    !REAL(r_std), DIMENSION(npts, nvm)                                 :: nsenobs_real  
    !REAL(r_std), DIMENSION(npts, nvm)                                 :: ndrpobs_real  

    REAL(r_std), DIMENSION(npts, nvm)                                 :: etatvernal_real  
    REAL(r_std), DIMENSION(npts, nvm)                                 :: gelee_real  

    ! To store variables names for I/O
    CHARACTER(LEN=80) :: var_name
    ! string suffix indicating an index
    CHARACTER(LEN=10) :: part_str
    ! string suffix indicating biomass type crops
    CHARACTER(LEN=10) :: box_str

    !var_name = 'f_crop_init'
    !WHERE (f_crop_init)
    !   f_crop_init_real = un
    !ELSEWHERE
    !   f_crop_init_real = zero
    !ENDWHERE
    !CALL restput_p (rest_id_stomate, var_name, nbp_glo,   1, 1, itime, &
    !     &                 f_crop_init_real, 'scatter', nbp_glo, index_g)
   
    var_name = 'f_crop_recycle'
    WHERE (f_crop_recycle(:,:))
       f_crop_recycle_real = un
    ELSEWHERE
       f_crop_recycle_real = zero
    ENDWHERE
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,   nvm, 1, itime, &
         &                 f_crop_recycle_real, 'scatter', nbp_glo, index_g)

    var_name = 'in_cycle'
    WHERE (in_cycle(:,:))
       in_cycle_real = un
    ELSEWHERE
       in_cycle_real = zero
    ENDWHERE
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,   nvm, 1, itime, &
         &                 in_cycle_real, 'scatter', nbp_glo, index_g)

    var_name = 'f_sen_lai'
    WHERE (f_sen_lai(:,:))
       f_sen_lai_real = un
    ELSEWHERE
       f_sen_lai_real = zero
    ENDWHERE
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,   nvm, 1, itime, &
         &                 f_sen_lai_real, 'scatter', nbp_glo, index_g)
    
    var_name = 'st2m_max_daily'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &               st2m_max_daily, 'scatter', nbp_glo, index_g)


    CALL restput_p (rest_id_stomate, 'wut_cm_daily', nbp_glo, nvm     , 3, itime, &
            &               wut_cm_daily, 'scatter', nbp_glo, index_g)

    CALL restput_p (rest_id_stomate, 'wus_cm_daily', nbp_glo, nvm   , 3, itime, &
            &               wus_cm_daily, 'scatter', nbp_glo, index_g)

    var_name = 'evapot_daily'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &               evapot_daily, 'scatter', nbp_glo, index_g)

    var_name = 'pdbiomass'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               pdbiomass, 'scatter', nbp_glo, index_g)

    var_name = 'pdmasec'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               pdmasec, 'scatter', nbp_glo, index_g)

    var_name = 'masecveg'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               masecveg, 'scatter', nbp_glo, index_g)

    var_name = 'masec'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               masec, 'scatter', nbp_glo, index_g)

    var_name = 'dltams'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               dltams, 'scatter', nbp_glo, index_g)

    var_name = 'gdh_daily'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               gdh_daily, 'scatter', nbp_glo, index_g)

    var_name = 'phoi'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &               phoi, 'scatter', nbp_glo, index_g)

    var_name = 'onarretesomcourdrp'
    WHERE (onarretesomcourdrp(:, :))
       onarretesomcourdrp_real = un
    ELSEWHERE
       onarretesomcourdrp_real = zero
    ENDWHERE
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,   nvm, 1, itime, &
         &                 onarretesomcourdrp_real, 'scatter', nbp_glo, index_g)

    !var_name = 'codeulaivernal'
    !codeulaivernal_real = FLOAT(codeulaivernal)
    !CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
    !     &               codeulaivernal_real, 'scatter', nbp_glo, index_g)
  
    !var_name = 'nlevobs'
    !nlevobs_real = FLOAT(nlevobs)
    !CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
    !     &               nlevobs_real, 'scatter', nbp_glo, index_g)

    !var_name = 'namfobs'
    !namfobs_real = FLOAT(namfobs)
    !CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
    !     &               namfobs_real, 'scatter', nbp_glo, index_g)

    !var_name = 'nfloobs'
    !nfloobs_real = FLOAT(nfloobs)
    !CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
    !     &               nfloobs_real, 'scatter', nbp_glo, index_g)

    !var_name = 'nlanobs'
    !nlanobs_real = FLOAT(nlanobs)
    !CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
    !     &               nlanobs_real, 'scatter', nbp_glo, index_g)

    !var_name = 'nlaxobs'
    !nlaxobs_real = FLOAT(nlaxobs)
    !CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
    !     &               nlaxobs_real, 'scatter', nbp_glo, index_g)

    !var_name = 'nmatobs'
    !nmatobs_real = FLOAT(nmatobs)
    !CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
    !     &               nmatobs_real, 'scatter', nbp_glo, index_g)

    !var_name = 'nrecobs'
    !nrecobs_real = FLOAT(nrecobs)
    !CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
    !     &               nrecobs_real, 'scatter', nbp_glo, index_g)

    !var_name = 'nsenobs'
    !nsenobs_real = FLOAT(nsenobs)
    !CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
    !     &               nsenobs_real, 'scatter', nbp_glo, index_g)


    !var_name = 'ndrpobs'
    !ndrpobs_real = FLOAT(ndrpobs)
    !CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
    !     &               ndrpobs_real, 'scatter', nbp_glo, index_g)


    var_name = 'nsendltams'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               nsendltams, 'scatter', nbp_glo, index_g)

    var_name = 'nsendltai'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               nsendltai, 'scatter', nbp_glo, index_g)

    var_name = 'nsenpfeuilverte'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               nsenpfeuilverte, 'scatter', nbp_glo, index_g)

    var_name = 'nsendurvie'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               nsendurvie, 'scatter', nbp_glo, index_g)

    var_name = 'nsenndurvie'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               nsenndurvie, 'scatter', nbp_glo, index_g)

    var_name = 'densiteequiv'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               densiteequiv, 'scatter', nbp_glo, index_g)

    var_name = 'nplt'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               nplt, 'scatter', nbp_glo, index_g)

    var_name = 'tursla'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               tursla, 'scatter', nbp_glo, index_g)

    var_name = 'ssla'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               ssla, 'scatter', nbp_glo, index_g)

    var_name = 'pfeuilverte'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               pfeuilverte, 'scatter', nbp_glo, index_g)

    var_name = 'bsenlai'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               bsenlai, 'scatter', nbp_glo, index_g)

    var_name = 'zrac'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               zrac, 'scatter', nbp_glo, index_g)

    var_name = 'nrec'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               nrec, 'scatter', nbp_glo, index_g)

    var_name = 'nlan'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               nlan, 'scatter', nbp_glo, index_g)

    var_name = 'tcult'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               tcult, 'scatter', nbp_glo, index_g)

    var_name = 'udevair'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               udevair, 'scatter', nbp_glo, index_g)

    var_name = 'udevcult'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               udevcult, 'scatter', nbp_glo, index_g)

    var_name = 'ndrp'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               ndrp, 'scatter', nbp_glo, index_g)

    var_name = 'rfvi'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               rfvi, 'scatter', nbp_glo, index_g)

    var_name = 'nlev'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               nlev, 'scatter', nbp_glo, index_g)

    var_name = 'nger'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               nger, 'scatter', nbp_glo, index_g)

    var_name = 'etatvernal'
    WHERE ( etatvernal(:,:) )
       etatvernal_real = un
    ELSEWHERE
       etatvernal_real = zero
    ENDWHERE
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               etatvernal_real, 'scatter', nbp_glo, index_g)

    var_name = 'caljvc'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               caljvc, 'scatter', nbp_glo, index_g)

    var_name = 'rfpi'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               rfpi, 'scatter', nbp_glo, index_g)

    var_name = 'upvt'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               upvt, 'scatter', nbp_glo, index_g)

    var_name = 'utp'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               utp, 'scatter', nbp_glo, index_g)

    var_name = 'somcour'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               somcour, 'scatter', nbp_glo, index_g)

    var_name = 'somcourdrp'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               somcourdrp, 'scatter', nbp_glo, index_g)

    var_name = 'somcourutp'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               somcourutp, 'scatter', nbp_glo, index_g)

    var_name = 'tdevelop'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               tdevelop, 'scatter', nbp_glo, index_g)

    var_name = 'somtemp'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               somtemp, 'scatter', nbp_glo, index_g)

    var_name = 'somcourfauche'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               somcourfauche, 'scatter', nbp_glo, index_g)

    var_name = 'stpltger'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               stpltger, 'scatter', nbp_glo, index_g)

    var_name = 'R_stlaxsen'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               R_stlaxsen, 'scatter', nbp_glo, index_g)

    var_name = 'R_stamflax'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               R_stamflax, 'scatter', nbp_glo, index_g)

    var_name = 'R_stsenlan'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               R_stsenlan, 'scatter', nbp_glo, index_g)

    var_name = 'stlevflo'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               stlevflo, 'scatter', nbp_glo, index_g)

    var_name = 'nflo'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               nflo, 'scatter', nbp_glo, index_g)

    var_name = 'R_stlevdrp'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               R_stlevdrp, 'scatter', nbp_glo, index_g)

    var_name = 'R_stflodrp'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               R_stflodrp, 'scatter', nbp_glo, index_g)

    var_name = 'R_stdrpmat'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               R_stdrpmat, 'scatter', nbp_glo, index_g)

    var_name = 'nmat'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               nmat, 'scatter', nbp_glo, index_g)

    var_name = 'nlax'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               nlax, 'scatter', nbp_glo, index_g)

    var_name = 'nrecbutoir'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               nrecbutoir, 'scatter', nbp_glo, index_g)

    var_name = 'group'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               group, 'scatter', nbp_glo, index_g)

    var_name = 'ndebdes'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               ndebdes, 'scatter', nbp_glo, index_g)

    var_name = 'R_stdrpdes'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               R_stdrpdes, 'scatter', nbp_glo, index_g)

    var_name = 'densite'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               densite, 'scatter', nbp_glo, index_g)

    var_name = 'densitelev'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               densitelev, 'scatter', nbp_glo, index_g)

    var_name = 'coeflev'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               coeflev, 'scatter', nbp_glo, index_g)

    var_name = 'densiteger'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               densiteger, 'scatter', nbp_glo, index_g)

    var_name = 'somelong'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               somelong, 'scatter', nbp_glo, index_g)

    var_name = 'somger'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               somger, 'scatter', nbp_glo, index_g)

    var_name = 'humectation'
    WHERE (humectation(:, :))
       humectation_real = un
    ELSEWHERE
       humectation_real = zero
    ENDWHERE
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,   nvm, 1, itime, &
         &                 humectation_real, 'scatter', nbp_glo, index_g)

    var_name = 'nbjhumec'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               nbjhumec, 'scatter', nbp_glo, index_g)

    var_name = 'somtemphumec'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               somtemphumec, 'scatter', nbp_glo, index_g)

    var_name = 'stpltlev'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               stpltlev, 'scatter', nbp_glo, index_g)

    var_name = 'namf'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               namf, 'scatter', nbp_glo, index_g)

    var_name = 'stmatrec'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               stmatrec, 'scatter', nbp_glo, index_g)

    var_name = 'tustress'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               tustress, 'scatter', nbp_glo, index_g)

    var_name = 'lai'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               lai, 'scatter', nbp_glo, index_g)

    var_name = 'somfeuille'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               somfeuille, 'scatter', nbp_glo, index_g)

    var_name = 'pdlai'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               pdlai, 'scatter', nbp_glo, index_g)

    var_name = 'nbfeuille'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               nbfeuille, 'scatter', nbp_glo, index_g)

    var_name = 'reajust'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               reajust, 'scatter', nbp_glo, index_g)

    var_name = 'ulai'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               ulai, 'scatter', nbp_glo, index_g)

    var_name = 'pdulai'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               pdulai, 'scatter', nbp_glo, index_g)

    var_name = 'efdensite'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               efdensite, 'scatter', nbp_glo, index_g)

    var_name = 'tempeff'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               tempeff, 'scatter', nbp_glo, index_g)

    var_name = 'nstopfeuille'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               nstopfeuille, 'scatter', nbp_glo, index_g)

    var_name = 'deltai'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               deltai, 'scatter', nbp_glo, index_g)

    var_name = 'vmax'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               vmax, 'scatter', nbp_glo, index_g)

    var_name = 'nsen'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               nsen, 'scatter', nbp_glo, index_g)

    var_name = 'laisen'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               laisen, 'scatter', nbp_glo, index_g)

    var_name = 'pdlaisen'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               pdlaisen, 'scatter', nbp_glo, index_g)

    var_name = 'dltaisenat'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               dltaisenat, 'scatter', nbp_glo, index_g)

    var_name = 'nsencour'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               nsencour, 'scatter', nbp_glo, index_g)

    var_name = 'dltamsen'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               dltamsen, 'scatter', nbp_glo, index_g)

    var_name = 'dltaisen'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               dltaisen, 'scatter', nbp_glo, index_g)

    var_name = 'fgellev'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               fgellev, 'scatter', nbp_glo, index_g)

    var_name = 'gelee'
    WHERE ( gelee(:,:) )
       gelee_real = un
    ELSEWHERE
       gelee_real = zero
    ENDWHERE
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               gelee_real, 'scatter', nbp_glo, index_g)

    var_name = 'fstressgel'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               fstressgel, 'scatter', nbp_glo, index_g)

    var_name = 'laisen'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               laisen, 'scatter', nbp_glo, index_g)

    var_name = 'R_stlevamf'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               R_stlevamf, 'scatter', nbp_glo, index_g)

    var_name = 'dernier_n'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               dernier_n, 'scatter', nbp_glo, index_g)

    var_name = 'durvieI'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               durvieI, 'scatter', nbp_glo, index_g)

    var_name = 'durvie'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               durvie, 'scatter', nbp_glo, index_g)

    var_name = 'ndebsen'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               ndebsen, 'scatter', nbp_glo, index_g)

    var_name = 'somsenreste'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               somsenreste, 'scatter', nbp_glo, index_g)

    var_name = 'humrel'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               humrel, 'scatter', nbp_glo, index_g)

    var_name = 'swfac'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               swfac, 'scatter', nbp_glo, index_g)

    var_name = 'turfac'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               turfac, 'scatter', nbp_glo, index_g)

    var_name = 'senfac'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               senfac, 'scatter', nbp_glo, index_g)
  

    var_name = 'mafeuiljaune'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               mafeuiljaune, 'scatter', nbp_glo, index_g)

    var_name = 'msneojaune'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               msneojaune, 'scatter', nbp_glo, index_g)

!    var_name = 'myvar'
!    CALL restput_p(rest_id_stomate, var_name, nbp_glo, nvm, nboxmax, itime, &
!        &                myvar, 'scatter', nbp_glo, index_g)

     CALL restput_p (rest_id_stomate, 'box_ndays', nbp_glo, nvm     , nboxmax, itime, &
            &                REAL(box_ndays), 'scatter', nbp_glo, index_g)

     CALL restput_p(rest_id_stomate, 'box_lai', nbp_glo, nvm     , nboxmax, itime, &
            &                box_lai, 'scatter', nbp_glo, index_g) 

     CALL restput_p(rest_id_stomate, 'box_lairem', nbp_glo, nvm     , nboxmax, itime, &
            &                box_lairem, 'scatter', nbp_glo, index_g)

     CALL restput_p(rest_id_stomate, 'box_tdev', nbp_glo, nvm     , nboxmax, itime, &
            &                box_tdev, 'scatter', nbp_glo, index_g)

     CALL restput_p(rest_id_stomate, 'box_biom', nbp_glo, nvm     , nboxmax, itime, &
            &                box_biom, 'scatter', nbp_glo, index_g)

    CALL restput_p(rest_id_stomate, 'box_biomrem', nbp_glo, nvm     , nboxmax, itime, &
            &                box_biomrem, 'scatter', nbp_glo, index_g)

    CALL restput_p(rest_id_stomate, 'box_durage', nbp_glo, nvm     , nboxmax, itime, &
            &                box_durage, 'scatter', nbp_glo, index_g)

    CALL restput_p(rest_id_stomate, 'box_somsenbase', nbp_glo, nvm     , nboxmax, itime, &
            &                box_somsenbase, 'scatter', nbp_glo, index_g)

 
    ! STICS:: CARBON ALLOCATION
    CALL restput_p (rest_id_stomate, 'v_dltams', nbp_glo, nvm     , vlength, itime, &
         &               v_dltams, 'scatter', nbp_glo, index_g)
 
    var_name = 'fgelflo'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               fgelflo, 'scatter', nbp_glo, index_g)

    var_name = 'pdircarb'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               pdircarb, 'scatter', nbp_glo, index_g)


    var_name = 'ircarb'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               ircarb, 'scatter', nbp_glo, index_g)

    var_name = 'nbgrains'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               nbgrains, 'scatter', nbp_glo, index_g)

    var_name = 'pgrain'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               pgrain, 'scatter', nbp_glo, index_g)

    var_name = 'vitmoy'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               vitmoy, 'scatter', nbp_glo, index_g)

    var_name = 'nbgraingel'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               nbgraingel, 'scatter', nbp_glo, index_g)

    var_name = 'pgraingel'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               pgraingel, 'scatter', nbp_glo, index_g)

    var_name = 'dltags'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               dltags, 'scatter', nbp_glo, index_g)

    var_name = 'ftempremp'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               ftempremp, 'scatter', nbp_glo, index_g)

    var_name = 'magrain'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               magrain, 'scatter', nbp_glo, index_g)

    var_name = 'pdmagrain'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               pdmagrain, 'scatter', nbp_glo, index_g)


    var_name = 'nbj0remp'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               nbj0remp, 'scatter', nbp_glo, index_g)
 
    !var_name = 'nbj0remp'
    !CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
    !     &               nbj0remp, 'scatter', nbp_glo, index_g)

    var_name = 'pdsfruittot'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               pdsfruittot, 'scatter', nbp_glo, index_g)

    var_name = 'repracmax'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               repracmax, 'scatter', nbp_glo, index_g)

    var_name = 'repracmin'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               repracmin, 'scatter', nbp_glo, index_g)

    var_name = 'kreprac'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               kreprac, 'scatter', nbp_glo, index_g)

    var_name = 'somtemprac'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               somtemprac, 'scatter', nbp_glo, index_g)

    var_name = 'urac'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               urac, 'scatter', nbp_glo, index_g)

    var_name = 'reprac'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               reprac, 'scatter', nbp_glo, index_g)


    var_name = 'nstoprac'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               nstoprac, 'scatter', nbp_glo, index_g)


    var_name = 'c_leafb'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               c_leafb, 'scatter', nbp_glo, index_g)

    var_name = 'c_reserve'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               c_reserve, 'scatter', nbp_glo, index_g)
    
    var_name = 'gslen'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               gslen, 'scatter', nbp_glo, index_g)

    var_name = 'drylen'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &               drylen, 'scatter', nbp_glo, index_g)
    IF (ok_rotate) THEN
        var_name = 'cyc_num'
        CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
             &               cyc_num, 'scatter', nbp_glo, index_g)
    
        var_name = 'cyc_num_tot'
        CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
             &               cyc_num_tot, 'scatter', nbp_glo, index_g)
    
        var_name = 'rot_cmd_store'
        CALL restput_p (rest_id_stomate, var_name, nbp_glo, rot_cmd_max, cyc_rot_max, itime, &
             &               rot_cmd_store, 'scatter', nbp_glo, index_g)
    ENDIF
    var_name = 'plantdate'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, cyc_rot_max, itime, &
         &                plantdate, 'scatter', nbp_glo, index_g)

    var_name = 'plantdate_now'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1   , itime, &
         &                plantdate_now, 'scatter', nbp_glo, index_g)

  END SUBROUTINE sticslai_io_writestart
  !-
  !===
  !-
  SUBROUTINE sticslai_io_readstart (npts, & 
       & f_crop_recycle, in_cycle, f_sen_lai, st2m_max_daily, wut_cm_daily, wus_cm_daily, evapot_daily, pdbiomass, pdmasec, &
       & masecveg, masec, dltams, gdh_daily, phoi, onarretesomcourdrp,  &
       & nsendltams, nsendltai, nsenpfeuilverte, nsendurvie, nsenndurvie, densiteequiv, &
       & nplt, tursla, ssla, pfeuilverte, bsenlai, &
       & zrac, nrec, nlan, tcult, udevair, udevcult, ndrp, rfvi, nlev, nger, etatvernal, &
       & caljvc, rfpi, upvt, utp, somcour, somcourdrp, somcourutp, tdevelop, somtemp, &
       & somcourfauche, stpltger, R_stamflax, R_stlaxsen, R_stsenlan, stlevflo, nflo, &
       & R_stlevdrp, R_stflodrp, R_stdrpmat, nmat, nlax, nrecbutoir, group, ndebdes, R_stdrpdes, densite, &
       & densitelev, coeflev, densiteger, somelong, somger, humectation, nbjhumec, &
       & somtemphumec, stpltlev, namf, stmatrec, tustress, lai, somfeuille, pdlai, &
       & nbfeuille, reajust, ulai, pdulai, efdensite, tempeff, nstopfeuille, deltai, vmax, nsen, &
       & laisen, pdlaisen, dltaisenat, nsencour, dltamsen, dltaisen, fgellev, &
       & gelee, fstressgel, R_stlevamf, dernier_n, durvieI, durvie, ndebsen, somsenreste, &
       & humrel, swfac, turfac, senfac, mafeuiljaune, msneojaune,&
       & v_dltams, fgelflo, pdircarb, ircarb, nbgrains, pgrain, vitmoy, nbgraingel, pgraingel, &
       & dltags, ftempremp, magrain, pdmagrain, nbj0remp, pdsfruittot, repracmax, repracmin, &
       & kreprac, somtemprac, urac, reprac, nstoprac, c_reserve, c_leafb, gslen, drylen, &
       & nboxmax, box_ndays, box_lai, box_lairem, box_tdev, box_biom, box_biomrem,box_durage, box_somsenbase )
!!!!! end crop, xuhui
    !---------------------------------------------------------------------
    !- read start file
    !---------------------------------------------------------------------
    !-
    ! 0 declarations
    !-
    ! 0.1 input
    !-
    ! Domain size
    INTEGER(i_std),INTENT(in) :: npts
    ! Indices of the points on the map
!!!!! crops

    !LOGICAL, INTENT(OUT)       :: f_crop_init 
    LOGICAL, DIMENSION(npts, nvm), INTENT(OUT)       :: f_crop_recycle 
    LOGICAL, DIMENSION(npts, nvm), INTENT(OUT)       :: in_cycle 
    LOGICAL, DIMENSION(npts, nvm), INTENT(OUT)       :: f_sen_lai 
    ! daily maximum 2 meter temperatures (K)
    REAL(r_std), DIMENSION(npts), INTENT(OUT)      :: st2m_max_daily
    ! daily value of soil temperature at the resolution of 1 cm, the second dimension is 3
    ! the three layers around sowing layer
    REAL(r_std), DIMENSION(npts, nvm, 3), INTENT(OUT)    :: wut_cm_daily
    ! daily mean value of soil relative humidity at the resolution of 1 cm, the second dimension is 3
    ! the three layers around sowing layer
    REAL(r_std), DIMENSION(npts, nvm, 3), INTENT(OUT)    :: wus_cm_daily
    ! daily potential evapotranspiration 
    REAL(r_std), DIMENSION(npts), INTENT(OUT)      :: evapot_daily
    ! biomass of previous day, t/ha
    REAL(r_std), DIMENSION(npts, nvm), INTENT(OUT)      :: pdbiomass
    REAL(r_std), DIMENSION(npts, nvm), INTENT(OUT)      :: pdmasec
    ! vegetative biomass
    REAL(r_std), DIMENSION(npts, nvm), INTENT(OUT)      :: masecveg   
    ! aboveground dry matter (t ha-1) 
    REAL(r_std), DIMENSION(npts, nvm), INTENT(OUT)      :: masec
    ! growth rate of plant, it means the delta total biomass increment (t ha-1)
    REAL(r_std), DIMENSION(npts, nvm), INTENT(OUT)      :: dltams
    ! daily gdh calculated according to halfhourly temperature // transmitted from stomate.f90 gdh_daily
    REAL(r_std), DIMENSION(npts, nvm), INTENT(OUT)         :: gdh_daily 
    ! Photoperiod // hours
    REAL(r_std),  DIMENSION(npts), INTENT(OUT)                             :: phoi 

    !  
    LOGICAL, DIMENSION(npts, nvm), INTENT(OUT)           :: onarretesomcourdrp 
    !INTEGER(i_std), DIMENSION(nvm), INTENT(OUT)                           :: codeulaivernal 
    !INTEGER(i_std), DIMENSION(npts, nvm), INTENT(OUT)                :: nlevobs    ! the following variables ended with obs are only used for forcing simulation.  
    !INTEGER(i_std), DIMENSION(npts, nvm), INTENT(OUT)                :: namfobs    ! the initial value should be always 999
    !INTEGER(i_std), DIMENSION(npts, nvm), INTENT(OUT)                :: nfloobs  
    !INTEGER(i_std), DIMENSION(npts, nvm), INTENT(OUT)                :: nlanobs  
    !INTEGER(i_std), DIMENSION(npts, nvm), INTENT(OUT)                :: nlaxobs  
    !INTEGER(i_std), DIMENSION(npts, nvm), INTENT(OUT)                :: nmatobs  
    !INTEGER(i_std), DIMENSION(npts, nvm), INTENT(OUT)                :: nrecobs  
    !INTEGER(i_std), DIMENSION(npts, nvm), INTENT(OUT)                :: nsenobs  
    !INTEGER(i_std), DIMENSION(npts, nvm), INTENT(OUT)                :: ndrpobs  

    ! LAIdev SPECIFIC 
    REAL(r_std), DIMENSION(npts, nvm), INTENT(OUT)        :: nsendltams
    REAL(r_std), DIMENSION(npts, nvm), INTENT(OUT)        :: nsendltai
    REAL(r_std), DIMENSION(npts, nvm), INTENT(OUT)        :: nsenpfeuilverte
    REAL(r_std), DIMENSION(npts, nvm), INTENT(OUT)        :: nsendurvie
    REAL(r_std), DIMENSION(npts, nvm), INTENT(OUT)        :: nsenndurvie
    REAL(r_std), DIMENSION(npts, nvm), INTENT(OUT)        :: densiteequiv
    INTEGER(i_std), DIMENSION(npts, nvm), INTENT(OUT)        :: nplt
    REAL(r_std), DIMENSION(npts, nvm), INTENT(OUT)        :: tursla
    REAL(r_std), DIMENSION(npts, nvm), INTENT(OUT)        :: ssla
    REAL(r_std), DIMENSION(npts, nvm), INTENT(OUT)        :: pfeuilverte
    REAL(r_std), DIMENSION(npts, nvm), INTENT(OUT)        :: bsenlai
    
    ! variables are involved in DEVELOPMENT

    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: zrac
    INTEGER(i_std), DIMENSION(npts, nvm), INTENT(OUT)        :: nrec
    INTEGER(i_std), DIMENSION(npts, nvm), INTENT(OUT)        :: nlan
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: tcult
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: udevair
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: udevcult
    INTEGER(i_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: ndrp
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: rfvi
    INTEGER(i_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: nlev
    INTEGER(i_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: nger
    logical,    DIMENSION(npts, nvm), INTENT(OUT)        :: etatvernal
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: caljvc
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: rfpi
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: upvt
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: utp
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: somcour
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: somcourdrp
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: somcourutp
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: tdevelop
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: somtemp
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: somcourfauche
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: stpltger
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: R_stamflax
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: R_stlaxsen
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: R_stsenlan
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: stlevflo
    INTEGER(i_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: nflo
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: R_stlevdrp
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: R_stflodrp
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: R_stdrpmat
    INTEGER(i_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: nmat
    INTEGER(i_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: nlax
    INTEGER(i_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: nrecbutoir
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: group
    INTEGER(i_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: ndebdes
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: R_stdrpdes
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: densite
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: densitelev
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: coeflev
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: densiteger
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: somelong
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: somger
    logical,    DIMENSION(npts, nvm), INTENT(OUT)        :: humectation
    INTEGER(i_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: nbjhumec
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: somtemphumec
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: stpltlev
    INTEGER(i_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: namf
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: stmatrec
 
    ! these variables are involved in Lai_calculation
     
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: tustress
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: lai
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: somfeuille
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: pdlai
    INTEGER(i_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: nbfeuille
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: reajust
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: ulai
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: pdulai
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: efdensite
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: tempeff
    INTEGER(i_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: nstopfeuille
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: deltai
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: vmax
    INTEGER(i_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: nsen
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: laisen
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: pdlaisen
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)        :: dltaisenat

    ! these variables are involved in the LAIsenescence

    INTEGER(i_std),    DIMENSION(npts, nvm), INTENT(OUT)      :: nsencour
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)      :: dltamsen
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)      :: dltaisen
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)      :: fgellev
    logical,    DIMENSION(npts, nvm), INTENT(OUT)      :: gelee
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)      :: fstressgel
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)      :: R_stlevamf
    INTEGER(i_std),    DIMENSION(npts, nvm), INTENT(OUT)      :: dernier_n
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)      :: durvieI
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)      :: durvie
    INTEGER(i_std),    DIMENSION(npts, nvm), INTENT(OUT)      :: ndebsen
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)      :: somsenreste

    INTEGER(i_std), INTENT(IN)                             :: nboxmax   
    INTEGER(i_std),    DIMENSION(npts, nvm, nboxmax), INTENT(OUT)      :: box_ndays    
    REAL(r_std),    DIMENSION(npts, nvm, nboxmax)                   :: boxtemp
    REAL(r_std),    DIMENSION(npts, nvm, nboxmax), INTENT(OUT)      :: box_lai
    REAL(r_std),    DIMENSION(npts, nvm, nboxmax), INTENT(OUT)      :: box_lairem
    REAL(r_std),    DIMENSION(npts, nvm, nboxmax), INTENT(OUT)      :: box_tdev
    REAL(r_std),    DIMENSION(npts, nvm, nboxmax), INTENT(OUT)      :: box_biom
    REAL(r_std),    DIMENSION(npts, nvm, nboxmax), INTENT(OUT)      :: box_biomrem
    REAL(r_std),    DIMENSION(npts, nvm, nboxmax), INTENT(OUT)      :: box_durage
    REAL(r_std),    DIMENSION(npts, nvm, nboxmax), INTENT(OUT)      :: box_somsenbase
    

    ! these variables are involved in STRESS calculation
    
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)      :: humrel
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)      :: swfac
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)      :: turfac
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)      :: senfac

    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)      :: mafeuiljaune
    REAL(r_std),    DIMENSION(npts, nvm), INTENT(OUT)      :: msneojaune
    ! these variables are involved in the CARBON ALLOCATION calculation

    ! grain related   
    REAL(r_std),    DIMENSION(npts, nvm, vlength)      ,INTENT(OUT)       :: v_dltams
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(OUT)       :: fgelflo
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(OUT)       :: pdircarb
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(OUT)       :: ircarb
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(OUT)       :: nbgrains
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(OUT)       :: pgrain
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(OUT)       :: vitmoy
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(OUT)       :: nbgraingel
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(OUT)       :: pgraingel
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(OUT)       :: dltags
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(OUT)       :: ftempremp
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(OUT)       :: magrain
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(OUT)       :: pdmagrain
    INTEGER(i_std), DIMENSION(npts, nvm)      ,INTENT(OUT)       :: nbj0remp 
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(OUT)       :: pdsfruittot

    ! reprac related

    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(OUT)       :: repracmax
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(OUT)       :: repracmin
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(OUT)       :: kreprac
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(OUT)       :: somtemprac
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(OUT)       :: urac
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(OUT)       :: reprac
    INTEGER(i_std), DIMENSION(npts, nvm)      ,INTENT(OUT)       :: nstoprac

    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(OUT)       :: c_reserve
    REAL(r_std),    DIMENSION(npts, nvm)      ,INTENT(OUT)       :: c_leafb

    INTEGER(i_std), DIMENSION(npts, nvm)      ,INTENT(OUT)       :: gslen 
    INTEGER(i_std), DIMENSION(npts, nvm)      ,INTENT(OUT)       :: drylen 

!!!!! xuhui
    ! STICS--local
    CHARACTER(LEN=100)                                            :: var_name !! for restget_p
    REAL(r_std), DIMENSION(npts, nvm)                                 :: in_cycle_real        
    REAL(r_std), DIMENSION(npts, nvm)                                 :: f_sen_lai_real        
    REAL(r_std), DIMENSION(npts, nvm)                                 :: f_crop_recycle_real        
    REAL(r_std), DIMENSION(npts, nvm)                                 :: onarretesomcourdrp_real        
    REAL(r_std), DIMENSION(npts, nvm)                                 :: humectation_real        

    !REAL(r_std), DIMENSION(nvm)                                       :: codeulaivernal_real        
    
    !REAL(r_std), DIMENSION(npts, nvm)                                 :: nlevobs_real
    !REAL(r_std), DIMENSION(npts, nvm)                                 :: namfobs_real
    !REAL(r_std), DIMENSION(npts, nvm)                                 :: nfloobs_real  
    !REAL(r_std), DIMENSION(npts, nvm)                                 :: nlanobs_real  
    !REAL(r_std), DIMENSION(npts, nvm)                                 :: nlaxobs_real  
    !REAL(r_std), DIMENSION(npts, nvm)                                 :: nmatobs_real  
    !REAL(r_std), DIMENSION(npts, nvm)                                 :: nrecobs_real  
    !REAL(r_std), DIMENSION(npts, nvm)                                 :: nsenobs_real  
    !REAL(r_std), DIMENSION(npts, nvm)                                 :: ndrpobs_real  

    REAL(r_std), DIMENSION(npts, nvm)                                 :: etatvernal_real  
    REAL(r_std), DIMENSION(npts, nvm)                                 :: gelee_real  


    INTEGER(i_std)                                :: l,k,ji, jv, i, j, m      !! indices    
!!!!! xuhui

    !---------------------------------------------------------------------
    IF (printlev >= 3) WRITE(numout,*) 'Entering readstart_sticslai_io'

!ENDJCADD
!!!!! crops

    !f_crop_init_real = val_exp
    !var_name = 'f_crop_init'
    !CALL restget_p (rest_id_stomate, var_name, nbp_glo,   1, 1, itime, &
    !     &                .TRUE., f_crop_init_real, 'gather', nbp_glo, index_g)
    !IF (f_crop_init_real == val_exp) f_crop_init_real = zero
    !WHERE (f_crop_init_real == 1)
    !   f_crop_init = .TRUE.
    !ELSEWHERE
    !   f_crop_init = .FALSE.
    !ENDWHERE
   
    f_crop_recycle_real(:, :) = val_exp
    var_name = 'f_crop_recycle'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo,   nvm, 1, itime, &
         &                .TRUE., f_crop_recycle_real, 'gather', nbp_glo, index_g)
    IF (ALL(f_crop_recycle_real(:, :) == val_exp)) f_crop_recycle_real(:, :) = zero
    WHERE (f_crop_recycle_real(:, :) == un)
       f_crop_recycle = .TRUE.
    ELSEWHERE
       f_crop_recycle = .FALSE.
    ENDWHERE

    in_cycle_real(:, :) = val_exp
    var_name = 'in_cycle'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo,   nvm, 1, itime, &
         &                .TRUE., in_cycle_real, 'gather', nbp_glo, index_g)
    IF (ALL(in_cycle_real(:, :) == val_exp)) in_cycle_real(:, :) = zero
    WHERE (in_cycle_real(:, :) == un)
       in_cycle = .TRUE.
    ELSEWHERE
       in_cycle = .FALSE.
    ENDWHERE
    
    f_sen_lai_real(:, :) = val_exp
    var_name = 'f_sen_lai'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo,   nvm, 1, itime, &
         &                .TRUE., f_sen_lai_real, 'gather', nbp_glo, index_g)
    IF (ALL(f_sen_lai_real(:, :) == val_exp)) f_sen_lai_real(:, :) = un
    WHERE (f_sen_lai_real(:, :) == un)
       f_sen_lai = .TRUE.
    ELSEWHERE
       f_sen_lai = .FALSE.
    ENDWHERE
    
    st2m_max_daily(:) = val_exp
    var_name = 'st2m_max_daily'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., st2m_max_daily, 'gather', nbp_glo, index_g)
    IF (ALL(st2m_max_daily(:) == val_exp)) st2m_max_daily(:) = zero

    wut_cm_daily(:, :, :) = val_exp
    CALL restget_p (rest_id_stomate, 'wut_cm_daily', nbp_glo, nvm, 3, itime, &
         &              .TRUE., wut_cm_daily, 'gather', nbp_glo, index_g)
    IF (ALL(wut_cm_daily == val_exp)) wut_cm_daily = zero


    wus_cm_daily(:, :, :) = val_exp
    CALL restget_p (rest_id_stomate, 'wus_cm_daily', nbp_glo,  nvm , 3, itime, &
         &              .TRUE., wus_cm_daily, 'gather', nbp_glo, index_g)
    IF (ALL(wus_cm_daily == val_exp)) wus_cm_daily = zero


    evapot_daily(:) = val_exp
    var_name = 'evapot_daily'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., evapot_daily, 'gather', nbp_glo, index_g)
    IF (ALL(evapot_daily(:) == val_exp)) evapot_daily(:) = zero

    pdbiomass(:, :) = val_exp
    var_name = 'pdbiomass'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., pdbiomass, 'gather', nbp_glo, index_g)
    IF (ALL(pdbiomass(:, :) == val_exp)) pdbiomass(:, :) = zero

    pdmasec(:, :) = val_exp
    var_name = 'pdmasec'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., pdmasec, 'gather', nbp_glo, index_g)
    IF (ALL(pdmasec(:, :) == val_exp)) pdmasec(:, :) = zero

    masecveg(:, :) = val_exp
    var_name = 'masecveg'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., masecveg, 'gather', nbp_glo, index_g)
    IF (ALL(masecveg(:, :) == val_exp)) masecveg(:, :) = zero

    masec(:, :) = val_exp
    var_name = 'masec'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., masec, 'gather', nbp_glo, index_g)
    IF (ALL(masec(:, :) == val_exp)) masec(:, :) = zero

    dltams(:, :) = val_exp
    var_name = 'dltams'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., dltams, 'gather', nbp_glo, index_g)
    IF (ALL(dltams(:, :) == val_exp)) dltams(:, :) = zero

    gdh_daily(:, :) = val_exp
    var_name = 'gdh_daily'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., gdh_daily, 'gather', nbp_glo, index_g)
    IF (ALL(gdh_daily(:, :) == val_exp)) gdh_daily(:, :) = zero

    phoi(:) = val_exp
    var_name = 'phoi'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., phoi, 'gather', nbp_glo, index_g)
    IF (ALL(phoi(:) == val_exp)) phoi(:) = zero

    onarretesomcourdrp_real(:, :) = val_exp
    var_name = 'onarretesomcourdrp'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo,   nvm, 1, itime, &
         &                .TRUE., onarretesomcourdrp_real, 'gather', nbp_glo, index_g)
    IF (ALL(onarretesomcourdrp_real(:, :) == val_exp)) onarretesomcourdrp_real(:, :) = zero
    WHERE (onarretesomcourdrp_real(:, :) == un)
       onarretesomcourdrp = .TRUE.
    ELSEWHERE
       onarretesomcourdrp = .FALSE.
    ENDWHERE

    !codeulaivernal(:) = val_exp
    !var_name = 'codeulaivernal'
    !CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
    !     &              .TRUE., codeulaivernal_real, 'gather', nbp_glo, index_g)
    !IF (ALL(codeulaivernal_real(:) == val_exp)) codeulaivernal_real(:) = zero
    !codeulaivernal = INT(codeulaivernal_real)

    !nlevobs(:, :) = val_exp
    !var_name = 'nlevobs'
    !CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
    !     &              .TRUE., nlevobs_real, 'gather', nbp_glo, index_g)
    !IF (ALL(nlevobs_real(:, :) == val_exp)) nlevobs_real(:, :) = zero
    !nlevobs = INT(nlevobs_real)

    !namfobs(:, :) = val_exp
    !var_name = 'namfobs'
    !CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
    !     &              .TRUE., namfobs_real, 'gather', nbp_glo, index_g)
    !IF (ALL(namfobs_real(:, :) == val_exp)) namfobs_real(:, :) = zero
    !namfobs = INT(namfobs_real)

    !nfloobs(:, :) = val_exp
    !var_name = 'nfloobs'
    !CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
    !     &              .TRUE., nfloobs_real, 'gather', nbp_glo, index_g)
    !IF (ALL(nfloobs_real(:, :) == val_exp)) nfloobs_real(:, :) = zero
    !nfloobs = INT(nfloobs_real)

    !nlanobs(:, :) = val_exp
    !var_name = 'nlanobs'
    !CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
    !     &              .TRUE., nlanobs_real, 'gather', nbp_glo, index_g)
    !IF (ALL(nlanobs_real(:, :) == val_exp)) nlanobs_real(:, :) = zero
    !nlanobs = INT(nlanobs_real)

    !nlaxobs(:, :) = val_exp
    !var_name = 'nlaxobs'
    !CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
    !     &              .TRUE., nlaxobs_real, 'gather', nbp_glo, index_g)
    !IF (ALL(nlaxobs_real(:, :) == val_exp)) nlaxobs_real(:, :) = zero
    !nlaxobs = INT(nlaxobs_real)

    !nmatobs(:, :) = val_exp
    !var_name = 'nmatobs'
    !CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
    !     &              .TRUE., nmatobs_real, 'gather', nbp_glo, index_g)
    !IF (ALL(nmatobs_real(:, :) == val_exp)) nmatobs_real(:, :) = zero
    !nmatobs = INT(nmatobs_real)

    !nrecobs(:, :) = val_exp
    !var_name = 'nrecobs'
    !CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
    !     &              .TRUE., nrecobs_real, 'gather', nbp_glo, index_g)
    !IF (ALL(nrecobs_real(:, :) == val_exp)) nrecobs_real(:, :) = zero
    !nrecobs = INT(nrecobs_real)

    !nsenobs(:, :) = val_exp
    !var_name = 'nsenobs'
    !CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
    !     &              .TRUE., nsenobs_real, 'gather', nbp_glo, index_g)
    !IF (ALL(nsenobs_real(:, :) == val_exp)) nsenobs_real(:, :) = zero
    !nsenobs = INT(nsenobs_real)


    !ndrpobs(:, :) = val_exp
    !var_name = 'ndrpobs'
    !CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
    !     &              .TRUE., ndrpobs_real, 'gather', nbp_glo, index_g)
    !IF (ALL(ndrpobs_real(:, :) == val_exp)) ndrpobs_real(:, :) = zero
    !ndrpobs = INT(ndrpobs_real)




    nsendltams(:, :) = val_exp
    var_name = 'nsendltams'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., nsendltams, 'gather', nbp_glo, index_g)
    IF (ALL(nsendltams(:, :) == val_exp)) nsendltams(:, :) = zero

    nsendltai(:, :) = val_exp
    var_name = 'nsendltai'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., nsendltai, 'gather', nbp_glo, index_g)
    IF (ALL(nsendltai(:, :) == val_exp)) nsendltai(:, :) = zero

    nsenpfeuilverte(:, :) = val_exp
    var_name = 'nsenpfeuilverte'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., nsenpfeuilverte, 'gather', nbp_glo, index_g)
    IF (ALL(nsenpfeuilverte(:, :) == val_exp)) nsenpfeuilverte(:, :) = zero

    nsendurvie(:, :) = val_exp
    var_name = 'nsendurvie'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., nsendurvie, 'gather', nbp_glo, index_g)
    IF (ALL(nsendurvie(:, :) == val_exp)) nsendurvie(:, :) = zero

    nsenndurvie(:, :) = val_exp
    var_name = 'nsenndurvie'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., nsenndurvie, 'gather', nbp_glo, index_g)
    IF (ALL(nsenndurvie(:, :) == val_exp)) nsenndurvie(:, :) = zero

    densiteequiv(:, :) = val_exp
    var_name = 'densiteequiv'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., densiteequiv, 'gather', nbp_glo, index_g)
    IF (ALL(densiteequiv(:, :) == val_exp)) densiteequiv(:, :) = zero

    nplt(:, :) = val_exp
    var_name = 'nplt'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., nplt, 'gather', nbp_glo, index_g)
    IF (ALL(nplt(:, :) == val_exp)) nplt(:, :) = zero

    tursla(:, :) = val_exp
    var_name = 'tursla'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., tursla, 'gather', nbp_glo, index_g)
    IF (ALL(tursla(:, :) == val_exp)) tursla(:, :) = un

    ssla(:, :) = val_exp
    var_name = 'ssla'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., ssla, 'gather', nbp_glo, index_g)
    IF (ALL(ssla(:, :) == val_exp)) ssla(:, :) = zero

    pfeuilverte(:, :) = val_exp
    var_name = 'pfeuilverte'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., pfeuilverte, 'gather', nbp_glo, index_g)
    IF (ALL(pfeuilverte(:, :) == val_exp)) pfeuilverte(:, :) = zero

    bsenlai(:, :) = val_exp
    var_name = 'bsenlai'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., bsenlai, 'gather', nbp_glo, index_g)
    IF (ALL(bsenlai(:, :) == val_exp)) bsenlai(:, :) = zero

    zrac(:, :) = val_exp
    var_name = 'zrac'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., zrac, 'gather', nbp_glo, index_g)
    IF (ALL(zrac(:, :) == val_exp)) zrac(:, :) = zero

    nrec(:, :) = val_exp
    var_name = 'nrec'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., nrec, 'gather', nbp_glo, index_g)
    IF (ALL(nrec(:, :) == val_exp)) nrec(:, :) = zero

    nlan(:, :) = val_exp
    var_name = 'nlan'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., nlan, 'gather', nbp_glo, index_g)
    IF (ALL(nlan(:, :) == val_exp)) nlan(:, :) = zero

    tcult(:, :) = val_exp
    var_name = 'tcult'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., tcult, 'gather', nbp_glo, index_g)
    IF (ALL(tcult(:, :) == val_exp)) tcult(:, :) = zero

    udevair(:, :) = val_exp
    var_name = 'udevair'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., udevair, 'gather', nbp_glo, index_g)
    IF (ALL(udevair(:, :) == val_exp)) udevair(:, :) = zero

    udevcult(:, :) = val_exp
    var_name = 'udevcult'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., udevcult, 'gather', nbp_glo, index_g)
    IF (ALL(udevcult(:, :) == val_exp)) udevcult(:, :) = zero

    ndrp(:, :) = val_exp
    var_name = 'ndrp'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., ndrp, 'gather', nbp_glo, index_g)
    IF (ALL(ndrp(:, :) == val_exp)) ndrp(:, :) = zero

    rfvi(:, :) = val_exp
    var_name = 'rfvi'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., rfvi, 'gather', nbp_glo, index_g)
    IF (ALL(rfvi(:, :) == val_exp)) rfvi(:, :) = zero

    nlev(:, :) = val_exp
    var_name = 'nlev'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., nlev, 'gather', nbp_glo, index_g)
    IF (ALL(nlev(:, :) == val_exp)) nlev(:, :) = zero

    nger(:, :) = val_exp
    var_name = 'nger'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., nger, 'gather', nbp_glo, index_g)
    IF (ALL(nger(:, :) == val_exp)) nger(:, :) = zero

    etatvernal_real(:, :) = val_exp
    var_name = 'etatvernal'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., etatvernal_real, 'gather', nbp_glo, index_g)
    IF (ALL(etatvernal_real(:, :) == val_exp)) etatvernal_real(:, :) = zero
    WHERE (etatvernal_real(:,:) == un)
       etatvernal = .TRUE.
    ELSEWHERE
       etatvernal = .FALSE.
    ENDWHERE
    

    caljvc(:, :) = val_exp
    var_name = 'caljvc'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., caljvc, 'gather', nbp_glo, index_g)
    IF (ALL(caljvc(:, :) == val_exp)) caljvc(:, :) = zero

    rfpi(:, :) = val_exp
    var_name = 'rfpi'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., rfpi, 'gather', nbp_glo, index_g)
    IF (ALL(rfpi(:, :) == val_exp)) rfpi(:, :) = zero

    upvt(:, :) = val_exp
    var_name = 'upvt'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., upvt, 'gather', nbp_glo, index_g)
    IF (ALL(upvt(:, :) == val_exp)) upvt(:, :) = zero

    utp(:, :) = val_exp
    var_name = 'utp'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., utp, 'gather', nbp_glo, index_g)
    IF (ALL(utp(:, :) == val_exp)) utp(:, :) = zero

    somcour(:, :) = val_exp
    var_name = 'somcour'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., somcour, 'gather', nbp_glo, index_g)
    IF (ALL(somcour(:, :) == val_exp)) somcour(:, :) = zero

    somcourdrp(:, :) = val_exp
    var_name = 'somcourdrp'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., somcourdrp, 'gather', nbp_glo, index_g)
    IF (ALL(somcourdrp(:, :) == val_exp)) somcourdrp(:, :) = zero

    somcourutp(:, :) = val_exp
    var_name = 'somcourutp'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., somcourutp, 'gather', nbp_glo, index_g)
    IF (ALL(somcourutp(:, :) == val_exp)) somcourutp(:, :) = zero

    tdevelop(:, :) = val_exp
    var_name = 'tdevelop'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., tdevelop, 'gather', nbp_glo, index_g)
    IF (ALL(tdevelop(:, :) == val_exp)) tdevelop(:, :) = zero

    somtemp(:, :) = val_exp
    var_name = 'somtemp'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., somtemp, 'gather', nbp_glo, index_g)
    IF (ALL(somtemp(:, :) == val_exp)) somtemp(:, :) = zero

    somcourfauche(:, :) = val_exp
    var_name = 'somcourfauche'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., somcourfauche, 'gather', nbp_glo, index_g)
    IF (ALL(somcourfauche(:, :) == val_exp)) somcourfauche(:, :) = zero

    stpltger(:, :) = val_exp
    var_name = 'stpltger'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., stpltger, 'gather', nbp_glo, index_g)
    IF (ALL(stpltger(:, :) == val_exp)) THEN
       DO j= 1, nvm
          stpltger(:, j) = SP_stpltger(j)
       ENDDO
    ENDIF

    R_stlaxsen(:, :) = val_exp
    var_name = 'R_stlaxsen'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., R_stlaxsen, 'gather', nbp_glo, index_g)
    IF (ALL(R_stlaxsen(:, :) == val_exp)) THEN
       DO j= 1, nvm
          R_stlaxsen(:, j) = SP_stlaxsen(j)
       ENDDO
    ENDIF


    R_stamflax(:, :) = val_exp
    var_name = 'R_stamflax'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., R_stamflax, 'gather', nbp_glo, index_g)
    IF (ALL(R_stamflax(:, :) == val_exp)) THEN
       DO j= 1, nvm
          R_stamflax(:, j) = SP_stamflax(j)
       ENDDO
    ENDIF

    R_stsenlan(:, :) = val_exp
    var_name = 'R_stsenlan'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., R_stsenlan, 'gather', nbp_glo, index_g)
    IF (ALL(R_stsenlan(:, :) == val_exp)) THEN
       DO j= 1, nvm
          R_stsenlan(:, j) = SP_stsenlan(j)
       ENDDO
    ENDIF

    stlevflo(:, :) = val_exp
    var_name = 'stlevflo'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., stlevflo, 'gather', nbp_glo, index_g)
    IF (ALL(stlevflo(:, :) == val_exp)) THEN
       DO j= 1, nvm
          stlevflo(:, j) = SP_stlevdrp(j) - SP_stflodrp(j)
       ENDDO
    ENDIF

    nflo(:, :) = val_exp
    var_name = 'nflo'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., nflo, 'gather', nbp_glo, index_g)
    IF (ALL(nflo(:, :) == val_exp)) nflo(:, :) = zero

    R_stlevdrp(:, :) = val_exp
    var_name = 'R_stlevdrp'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., R_stlevdrp, 'gather', nbp_glo, index_g)
    IF (ALL(R_stlevdrp(:, :) == val_exp)) THEN
       DO j= 1, nvm
          R_stlevdrp(:, j) = SP_stlevdrp(j)
       ENDDO
    ENDIF

    R_stflodrp(:, :) = val_exp
    var_name = 'R_stflodrp'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., R_stflodrp, 'gather', nbp_glo, index_g)
    IF (ALL(R_stflodrp(:, :) == val_exp)) THEN
       DO j= 1, nvm
          R_stflodrp(:, j) = SP_stflodrp(j)
       ENDDO
    ENDIF

    R_stdrpmat(:, :) = val_exp
    var_name = 'R_stdrpmat'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., R_stdrpmat, 'gather', nbp_glo, index_g)
    IF (ALL(R_stdrpmat(:, :) == val_exp)) THEN
       DO j= 1, nvm
          R_stdrpmat(:, j) = SP_stdrpmat(j)
       ENDDO
    ENDIF

    nmat(:, :) = val_exp
    var_name = 'nmat'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., nmat, 'gather', nbp_glo, index_g)
    IF (ALL(nmat(:, :) == val_exp)) nmat(:, :) = zero

    nlax(:, :) = val_exp
    var_name = 'nlax'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., nlax, 'gather', nbp_glo, index_g)
    IF (ALL(nlax(:, :) == val_exp)) nlax(:, :) = zero

    nrecbutoir(:, :) = val_exp
    var_name = 'nrecbutoir'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., nrecbutoir, 'gather', nbp_glo, index_g)
    IF (ALL(nrecbutoir(:, :) == val_exp)) nrecbutoir(:, :) = 999.0

    group(:, :) = val_exp
    var_name = 'group'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., group, 'gather', nbp_glo, index_g)
    IF (ALL(group(:, :) == val_exp)) group(:, :) = zero

    ndebdes(:, :) = val_exp
    var_name = 'ndebdes'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., ndebdes, 'gather', nbp_glo, index_g)
    IF (ALL(ndebdes(:, :) == val_exp)) ndebdes(:, :) = zero

    R_stdrpdes(:, :) = val_exp
    var_name = 'R_stdrpdes'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., R_stdrpdes, 'gather', nbp_glo, index_g)
    IF (ALL(R_stdrpdes(:, :) == val_exp)) THEN
       DO j= 1, nvm
          R_stdrpdes(:, j) = SP_stdrpdes(j)
       ENDDO
    ENDIF

    densite(:, :) = val_exp
    var_name = 'densite'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., densite, 'gather', nbp_glo, index_g)
    IF (ALL(densite(:, :) == val_exp)) densite(:, :) = zero

    densitelev(:, :) = val_exp
    var_name = 'densitelev'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., densitelev, 'gather', nbp_glo, index_g)
    IF (ALL(densitelev(:, :) == val_exp)) densitelev(:, :) = zero

    coeflev(:, :) = val_exp
    var_name = 'coeflev'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., coeflev, 'gather', nbp_glo, index_g)
    IF (ALL(coeflev(:, :) == val_exp)) coeflev(:, :) = un

    densiteger(:, :) = val_exp
    var_name = 'densiteger'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., densiteger, 'gather', nbp_glo, index_g)
    IF (ALL(densiteger(:, :) == val_exp)) densiteger(:, :) = zero

    somelong(:, :) = val_exp
    var_name = 'somelong'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., somelong, 'gather', nbp_glo, index_g)
    IF (ALL(somelong(:, :) == val_exp)) somelong(:, :) = zero

    somger(:, :) = val_exp
    var_name = 'somger'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., somger, 'gather', nbp_glo, index_g)
    IF (ALL(somger(:, :) == val_exp)) somger(:, :) = zero

    humectation_real(:, :) = val_exp
    var_name = 'humectation'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo,   nvm, 1, itime, &
         &                .TRUE., humectation_real, 'gather', nbp_glo, index_g)
    IF (ALL(humectation_real(:, :) == val_exp)) humectation_real(:, :) = zero
    WHERE (humectation_real(:, :) == un)
       humectation = .TRUE.
    ELSEWHERE
       humectation = .FALSE.
    ENDWHERE

    nbjhumec(:, :) = val_exp
    var_name = 'nbjhumec'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., nbjhumec, 'gather', nbp_glo, index_g)
    IF (ALL(nbjhumec(:, :) == val_exp)) nbjhumec(:, :) = zero

    somtemphumec(:, :) = val_exp
    var_name = 'somtemphumec'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., somtemphumec, 'gather', nbp_glo, index_g)
    IF (ALL(somtemphumec(:, :) == val_exp)) somtemphumec(:, :) = zero

    stpltlev(:, :) = val_exp
    var_name = 'stpltlev'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., stpltlev, 'gather', nbp_glo, index_g)
    IF (ALL(stpltlev(:, :) == val_exp)) stpltlev(:, :) = zero

    namf(:, :) = val_exp
    var_name = 'namf'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., namf, 'gather', nbp_glo, index_g)
    IF (ALL(namf(:, :) == val_exp)) namf(:, :) = zero

    stmatrec(:, :) = val_exp
    var_name = 'stmatrec'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., stmatrec, 'gather', nbp_glo, index_g)
    IF (ALL(stmatrec(:, :) == val_exp)) stmatrec(:, :) = zero

    tustress(:, :) = val_exp
    var_name = 'tustress'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., tustress, 'gather', nbp_glo, index_g)
    IF (ALL(tustress(:, :) == val_exp)) tustress(:, :) = 1.0

    lai(:, :) = val_exp
    var_name = 'lai'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., lai, 'gather', nbp_glo, index_g)
    IF (ALL(lai(:, :) == val_exp)) lai(:, :) = zero

    somfeuille(:, :) = val_exp
    var_name = 'somfeuille'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., somfeuille, 'gather', nbp_glo, index_g)
    IF (ALL(somfeuille(:, :) == val_exp)) somfeuille(:, :) = zero

    pdlai(:, :) = val_exp
    var_name = 'pdlai'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., pdlai, 'gather', nbp_glo, index_g)
    IF (ALL(pdlai(:, :) == val_exp)) pdlai(:, :) = zero

    nbfeuille(:, :) = val_exp
    var_name = 'nbfeuille'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., nbfeuille, 'gather', nbp_glo, index_g)
    IF (ALL(nbfeuille(:, :) == val_exp)) nbfeuille(:, :) = zero

    reajust(:, :) = val_exp
    var_name = 'reajust'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., reajust, 'gather', nbp_glo, index_g)
    IF (ALL(reajust(:, :) == val_exp)) reajust(:, :) = zero

    ulai(:, :) = val_exp
    var_name = 'ulai'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., ulai, 'gather', nbp_glo, index_g)
    IF (ALL(ulai(:, :) == val_exp)) ulai(:, :) = zero

    pdulai(:, :) = val_exp
    var_name = 'pdulai'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., pdulai, 'gather', nbp_glo, index_g)
    IF (ALL(pdulai(:, :) == val_exp)) pdulai(:, :) = zero

    efdensite(:, :) = val_exp
    var_name = 'efdensite'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., efdensite, 'gather', nbp_glo, index_g)
    IF (ALL(efdensite(:, :) == val_exp)) efdensite(:, :) = zero

    tempeff(:, :) = val_exp
    var_name = 'tempeff'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., tempeff, 'gather', nbp_glo, index_g)
    IF (ALL(tempeff(:, :) == val_exp)) tempeff(:, :) = zero

    nstopfeuille(:, :) = val_exp
    var_name = 'nstopfeuille'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., nstopfeuille, 'gather', nbp_glo, index_g)
    IF (ALL(nstopfeuille(:, :) == val_exp)) nstopfeuille(:, :) = zero

    deltai(:, :) = val_exp
    var_name = 'deltai'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., deltai, 'gather', nbp_glo, index_g)
    IF (ALL(deltai(:, :) == val_exp)) deltai(:, :) = zero

    vmax(:, :) = val_exp
    var_name = 'vmax'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., vmax, 'gather', nbp_glo, index_g)
    IF (ALL(vmax(:, :) == val_exp)) vmax(:, :) = zero

    nsen(:, :) = val_exp
    var_name = 'nsen'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., nsen, 'gather', nbp_glo, index_g)
    IF (ALL(nsen(:, :) == val_exp)) nsen(:, :) = zero


    laisen(:, :) = val_exp
    var_name = 'laisen'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., laisen, 'gather', nbp_glo, index_g)
    IF (ALL(laisen(:, :) == val_exp)) laisen(:, :) = zero

    pdlaisen(:, :) = val_exp
    var_name = 'pdlaisen'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., pdlaisen, 'gather', nbp_glo, index_g)
    IF (ALL(pdlaisen(:, :) == val_exp)) pdlaisen(:, :) = zero

    dltaisenat(:, :) = val_exp
    var_name = 'dltaisenat'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., dltaisenat, 'gather', nbp_glo, index_g)
    IF (ALL(dltaisenat(:, :) == val_exp)) dltaisenat(:, :) = zero

    nsencour(:, :) = val_exp
    var_name = 'nsencour'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., nsencour, 'gather', nbp_glo, index_g)
    IF (ALL(nsencour(:, :) == val_exp)) nsencour(:, :) = zero
    nsencour = INT(nsencour)

    dltamsen(:, :) = val_exp
    var_name = 'dltamsen'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., dltamsen, 'gather', nbp_glo, index_g)
    IF (ALL(dltamsen(:, :) == val_exp)) dltamsen(:, :) = zero

    dltaisen(:, :) = val_exp
    var_name = 'dltaisen'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., dltaisen, 'gather', nbp_glo, index_g)
    IF (ALL(dltaisen(:, :) == val_exp)) dltaisen(:, :) = zero

    fgellev(:, :) = val_exp
    var_name = 'fgellev'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., fgellev, 'gather', nbp_glo, index_g)
    IF (ALL(fgellev(:, :) == val_exp)) fgellev(:, :) = un

    gelee_real(:, :) = val_exp
    var_name = 'gelee'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., gelee_real, 'gather', nbp_glo, index_g)
    IF (ALL(gelee_real(:, :) == val_exp)) gelee_real(:, :) = zero
    WHERE (gelee_real(:,:) == un)
       gelee = .TRUE.
    ELSEWHERE
       gelee = .FALSE.
    ENDWHERE

    fstressgel(:, :) = val_exp
    var_name = 'fstressgel'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., fstressgel, 'gather', nbp_glo, index_g)
    IF (ALL(fstressgel(:, :) == val_exp)) fstressgel(:, :) = zero

    R_stlevamf(:, :) = val_exp
    var_name = 'R_stlevamf'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., R_stlevamf, 'gather', nbp_glo, index_g)
    IF (ALL(R_stlevamf(:, :) == val_exp)) THEN
       DO j= 1, nvm
          R_stlevamf(:, j) = SP_stlevamf(j)
       ENDDO
    ENDIF

    dernier_n(:, :) = val_exp
    var_name = 'dernier_n'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., dernier_n, 'gather', nbp_glo, index_g)
    IF (ALL(dernier_n(:, :) == val_exp)) dernier_n(:, :) = zero

    durvieI(:, :) = val_exp
    var_name = 'durvieI'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., durvieI, 'gather', nbp_glo, index_g)
    IF (ALL(durvieI(:, :) == val_exp)) durvieI(:, :) = zero

    durvie(:, :) = val_exp
    var_name = 'durvie'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., durvie, 'gather', nbp_glo, index_g)
    IF (ALL(durvie(:, :) == val_exp)) durvie(:, :) = zero

    ndebsen(:, :) = val_exp
    var_name = 'ndebsen'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., ndebsen, 'gather', nbp_glo, index_g)
    IF (ALL(ndebsen(:, :) == val_exp)) ndebsen(:, :) = zero

    somsenreste(:, :) = val_exp
    var_name = 'somsenreste'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., somsenreste, 'gather', nbp_glo, index_g)
    IF (ALL(somsenreste(:, :) == val_exp)) somsenreste(:, :) = zero

    humrel(:, :) = val_exp
    var_name = 'humrel'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., humrel, 'gather', nbp_glo, index_g)
    IF (ALL(humrel(:, :) == val_exp)) humrel(:, :) = zero

    swfac(:, :) = val_exp
    var_name = 'swfac'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., swfac, 'gather', nbp_glo, index_g)
    IF (ALL(swfac(:, :) == val_exp)) swfac(:, :) = un

    turfac(:, :) = val_exp
    var_name = 'turfac'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., turfac, 'gather', nbp_glo, index_g)
    IF (ALL(turfac(:, :) == val_exp)) turfac(:, :) = un

    senfac(:, :) = val_exp
    var_name = 'senfac'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., senfac, 'gather', nbp_glo, index_g)
    IF (ALL(senfac(:, :) == val_exp)) senfac(:, :) = un
  

    mafeuiljaune(:, :) = val_exp
    var_name = 'mafeuiljaune'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., mafeuiljaune, 'gather', nbp_glo, index_g)
    IF (ALL(mafeuiljaune(:, :) == val_exp)) mafeuiljaune(:, :) = un
    
    msneojaune(:, :) = val_exp
    var_name = 'msneojaune'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., msneojaune, 'gather', nbp_glo, index_g)
    IF (ALL(msneojaune(:, :) == val_exp)) msneojaune(:, :) = un

    box_ndays(:,:,:) = 0
    box_lai(:,:,:) = 0.
    box_lairem(:,:,:) = 0.
    box_tdev(:,:,:) = 0.
    box_biom(:,:,:) = 0.
    box_biomrem(:,:,:) = 0.
    box_durage(:,:,:) = 0.
    box_somsenbase(:,:,:) = 0.
    CALL restget_p(rest_id_stomate, 'box_ndays', nbp_glo, nvm,  nboxmax, itime, &
             &              .TRUE., box_ndays, 'gather', nbp_glo, index_g)
    IF (ALL(box_ndays == val_exp)) box_ndays = 0

    CALL restget_p(rest_id_stomate, 'box_lai', nbp_glo, nvm,   nboxmax, itime, &
            &               .TRUE., box_lai, 'gather', nbp_glo, index_g)
    IF (ALL(box_lai == val_exp)) box_lai = 0.

    CALL restget_p(rest_id_stomate, 'box_lairem', nbp_glo, nvm,   nboxmax, itime, &
            &               .TRUE., box_lairem, 'gather', nbp_glo, index_g)
    IF (ALL(box_lairem == val_exp)) box_lairem = 0.

    CALL restget_p(rest_id_stomate, 'box_tdev', nbp_glo, nvm,   nboxmax, itime, &
            &               .TRUE., box_tdev, 'gather', nbp_glo, index_g)
    IF (ALL(box_tdev == val_exp)) box_tdev = 0.

    CALL restget_p(rest_id_stomate, 'box_biom', nbp_glo, nvm,   nboxmax, itime, &
            &               .TRUE., box_biom, 'gather', nbp_glo, index_g)
    IF (ALL(box_biom == val_exp)) box_biom = 0.

    CALL restget_p(rest_id_stomate, 'box_biomrem', nbp_glo, nvm,   nboxmax, itime, &
            &               .TRUE., box_biomrem, 'gather', nbp_glo, index_g)
    IF (ALL(box_biomrem == val_exp)) box_biomrem = 0.

    CALL restget_p(rest_id_stomate, 'box_durage', nbp_glo, nvm,   nboxmax, itime, &
            &               .TRUE., box_durage, 'gather', nbp_glo, index_g)
    IF (ALL(box_durage == val_exp)) box_durage = 0.

    CALL restget_p(rest_id_stomate, 'box_somsenbase', nbp_glo, nvm,   nboxmax, itime, &
            &               .TRUE., box_somsenbase, 'gather', nbp_glo, index_g)
    IF (ALL(box_somsenbase == val_exp)) box_somsenbase = 0.
 
    ! STICS:: CARBON ALLOCATION
   
    
    v_dltams(:,:,:) = val_exp
    CALL restget_p (rest_id_stomate, 'v_dltams', nbp_glo, nvm  , vlength, itime, &
            &                .TRUE., v_dltams, 'gather', nbp_glo, index_g)
    IF (ALL(v_dltams == val_exp))  v_dltams = zero

    fgelflo(:, :) = val_exp
    var_name = 'fgelflo'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., fgelflo, 'gather', nbp_glo, index_g)
    IF (ALL(fgelflo(:, :) == val_exp)) fgelflo(:, :) = un

    pdircarb(:, :) = val_exp
    var_name = 'pdircarb'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., pdircarb, 'gather', nbp_glo, index_g)
    IF (ALL(pdircarb(:, :) == val_exp)) pdircarb(:, :) = zero

    ircarb(:, :) = val_exp
    var_name = 'ircarb'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., ircarb, 'gather', nbp_glo, index_g)
    IF (ALL(ircarb(:, :) == val_exp)) ircarb(:, :) = zero

    nbgrains(:, :) = val_exp
    var_name = 'nbgrains'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., nbgrains, 'gather', nbp_glo, index_g)
    IF (ALL(nbgrains(:, :) == val_exp)) nbgrains(:, :) = zero

    pgrain(:, :) = val_exp
    var_name = 'pgrain'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., pgrain, 'gather', nbp_glo, index_g)
    IF (ALL(pgrain(:, :) == val_exp)) pgrain(:, :) = zero

    vitmoy(:, :) = val_exp
    var_name = 'vitmoy'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., vitmoy, 'gather', nbp_glo, index_g)
    IF (ALL(vitmoy(:, :) == val_exp)) vitmoy(:, :) = zero

    nbgraingel(:, :) = val_exp
    var_name = 'nbgraingel'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., nbgraingel, 'gather', nbp_glo, index_g)
    IF (ALL(nbgraingel(:, :) == val_exp)) nbgraingel(:, :) = zero

    pgraingel(:, :) = val_exp
    var_name = 'pgraingel'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., pgraingel, 'gather', nbp_glo, index_g)
    IF (ALL(pgraingel(:, :) == val_exp)) pgraingel(:, :) = zero

    dltags(:, :) = val_exp
    var_name = 'dltags'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., dltags, 'gather', nbp_glo, index_g)
    IF (ALL(dltags(:, :) == val_exp)) dltags(:, :) = zero

    ftempremp(:, :) = val_exp
    var_name = 'ftempremp'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., ftempremp, 'gather', nbp_glo, index_g)
    IF (ALL(ftempremp(:, :) == val_exp)) ftempremp(:, :) = zero

    magrain(:, :) = val_exp
    var_name = 'magrain'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., magrain, 'gather', nbp_glo, index_g)
    IF (ALL(magrain(:, :) == val_exp)) magrain(:, :) = zero

    pdmagrain(:, :) = val_exp
    var_name = 'pdmagrain'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., pdmagrain, 'gather', nbp_glo, index_g)
    IF (ALL(pdmagrain(:, :) == val_exp)) pdmagrain(:, :) = zero

    nbj0remp(:, :) = val_exp
    var_name = 'nbj0remp'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., nbj0remp, 'gather', nbp_glo, index_g)
    IF (ALL(nbj0remp(:, :) == val_exp)) nbj0remp(:, :) = zero

   
    !nbj0remp(:, :) = val_exp
    !var_name = 'nbj0remp'
    !CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
    !     &              .TRUE., nbj0remp, 'gather', nbp_glo, index_g)
    !IF (ALL(nbj0remp(:, :) == val_exp)) nbj0remp(:, :) = zero

    pdsfruittot(:, :) = val_exp
    var_name = 'pdsfruittot'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., pdsfruittot, 'gather', nbp_glo, index_g)
    IF (ALL(pdsfruittot(:, :) == val_exp)) pdsfruittot(:, :) = zero

    repracmax(:, :) = val_exp
    var_name = 'repracmax'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., repracmax, 'gather', nbp_glo, index_g)
    IF (ALL(repracmax(:, :) == val_exp)) repracmax(:, :) = zero

    repracmin(:, :) = val_exp
    var_name = 'repracmin'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., repracmin, 'gather', nbp_glo, index_g)
    IF (ALL(repracmin(:, :) == val_exp)) repracmin(:, :) = zero

    kreprac(:, :) = val_exp
    var_name = 'kreprac'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., kreprac, 'gather', nbp_glo, index_g)
    IF (ALL(kreprac(:, :) == val_exp)) kreprac(:, :) = zero

    somtemprac(:, :) = val_exp
    var_name = 'somtemprac'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., somtemprac, 'gather', nbp_glo, index_g)
    IF (ALL(somtemprac(:, :) == val_exp)) somtemprac(:, :) = zero

    urac(:, :) = val_exp
    var_name = 'urac'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., urac, 'gather', nbp_glo, index_g)
    IF (ALL(urac(:, :) == val_exp)) urac(:, :) = zero

    reprac(:, :) = val_exp
    var_name = 'reprac'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., reprac, 'gather', nbp_glo, index_g)
    IF (ALL(reprac(:, :) == val_exp)) reprac(:, :) = zero

    c_reserve(:, :) = val_exp
    var_name = 'c_reserve'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., c_reserve, 'gather', nbp_glo, index_g)
    IF (ALL(c_reserve(:, :) == val_exp)) c_reserve(:, :) = zero


    nstoprac(:, :) = val_exp
    var_name = 'nstoprac'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., nstoprac, 'gather', nbp_glo, index_g)
    IF (ALL(nstoprac(:, :) == val_exp)) nstoprac(:, :) = zero


    var_name = 'c_leafb'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., c_leafb, 'gather', nbp_glo, index_g)
    IF (ALL(c_leafb(:, :) == val_exp)) c_leafb(:, :) = zero

    gslen(:, :) = val_exp
    var_name = 'gslen'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., gslen, 'gather', nbp_glo, index_g)
    IF (ALL(gslen(:, :) == val_exp)) gslen(:, :) = zero

    drylen(:, :) = val_exp
    var_name = 'drylen'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &              .TRUE., drylen, 'gather', nbp_glo, index_g)
    IF (ALL(drylen(:, :) == val_exp)) drylen(:, :) = zero
!!!!! xuhui
 
    IF (printlev >= 4) WRITE(numout,*) 'Leaving readstart'
    !-----------------------
  END SUBROUTINE sticslai_io_readstart
  !-
  !===
  !-
END MODULE sticslai_io

