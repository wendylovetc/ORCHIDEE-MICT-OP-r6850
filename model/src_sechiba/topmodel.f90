! ===================================================================================================\n
! MODULE        : topmodel 
!
! CONTACT       : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE       : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        Main interface for TOP Model
!!
!!\n DESCRIPTION : 
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
MODULE topmodel 

  USE sechiba_io_p ! setvar_p
  USE ioipsl_para
  USE topmodel_var
  USE topmodel_init
  USE topmodel_subgrid


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: topmodel_initialize, topmodel_readrestart, topmodel_clear, topmodel_history, &
                topmodel_writerestart, topmodel_main, topmodel_histdef, topmodel_parameters 

  !
  ! variables used inside topmodel module : declaration and initialisation
  !

! TOPMODEL
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: fsat             !! field capacity fraction
!$OMP THREADPRIVATE(fsat)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: fwet             !! wetland fraction with WTD = 0 cm
!$OMP THREADPRIVATE(fwet)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: fwt1             !! wetland fraction with WTD entre 0 et -3cm
!$OMP THREADPRIVATE(fwt1)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: fwt2             !! wetland fraction with WTD entre -3cm et -6cm
!$OMP THREADPRIVATE(fwt2)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: fwt3             !! wetland fraction with WTD entre ... et ...
!$OMP THREADPRIVATE(fwt3)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: fwt4             !! etc.
!$OMP THREADPRIVATE(fwt4)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: drunoff          !! runoff de Dunne
!$OMP THREADPRIVATE(drunoff)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: ZMEAN            !! statistiques de la fonction de distribution des indices topo au sein de chaque maille
!$OMP THREADPRIVATE(ZMEAN)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: ZSTDT
!$OMP THREADPRIVATE(ZSTDT)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: ZSKEW
!$OMP THREADPRIVATE(ZSKEW)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: ZMIN
!$OMP THREADPRIVATE(ZMIN)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: ZMAX
!$OMP THREADPRIVATE(ZMAX)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: ZM               !! parametre TOPMODEL
!$OMP THREADPRIVATE(ZM)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: ZZPAS            !! pas des veceturs d indice topo au sein de chaque maille
!$OMP THREADPRIVATE(ZZPAS)
! vecteurs calculees par TOPMODEL pour chaque maille (contenu = f(indice seuil); fsat = f(indice seuil); etc.)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: ZTAB_FSAT        
!$OMP THREADPRIVATE(ZTAB_FSAT)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: ZTAB_WTOP
!$OMP THREADPRIVATE(ZTAB_WTOP)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: ZTAB_FWET
!$OMP THREADPRIVATE(ZTAB_FWET)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: ZTAB_WTOP_WET
!$OMP THREADPRIVATE(ZTAB_WTOP_WET)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: mcs_grid         !! Saturation dim kjpindex
!$OMP THREADPRIVATE(mcs_grid) 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: mcw_grid         !! Wilting point dim kjpindex  
!$OMP THREADPRIVATE(mcw_grid) 

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: drunoff_tot      !! Surface runoff generated Dune process
!$OMP THREADPRIVATE(drunoff_tot)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: fwet_out         !! wetland fraction
!$OMP THREADPRIVATE(fwet_out)


CONTAINS

!! ================================================================================================================================
!! SUBROUTINE 	: topmodel_initialize
!!
!>\BRIEF         Initialization subroutine for top model
!!
!! DESCRIPTION :
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE topmodel_initialize (kjit, kjpindex, nbp_glo, lalo, veget_max, zmaxh, &
                                rest_id, index_g, &
                                mcs, mcw, soiltile)
    !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                         :: kjit             !! Time step number 
    INTEGER(i_std), INTENT(in)                         :: kjpindex         !! Domain size
    INTEGER(i_std), INTENT(in)                         :: nbp_glo          !!
    REAL(r_std), INTENT(in)                            :: lalo(kjpindex,2) !! Vector of latitude and longitudes (degree)        
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: veget_max        !! Max. fraction of vegetation type (LAI -> infty)
    REAL(r_std), INTENT(in)                            :: zmaxh            !! Max hydrol depth 
    INTEGER(i_std),INTENT (in)                         :: rest_id          !! Restart file identifier
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)   :: index_g          !! Indeces of the points on the map
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)      :: mcs, mcw
    REAL(r_std),DIMENSION (kjpindex,nstm), INTENT (in) :: soiltile         !! Fraction of each soil tile within vegtot (0-1, unitless)

!_ ======================================================================================================

    ! Module allocation
    CALL topmodel_alloc(kjpindex)

    ! Read restart files
    CALL topmodel_readrestart(kjit, kjpindex, nbp_glo, rest_id, index_g)

    !
    !Config Key   = HYDROL_FWET
    !Config Desc  = Initial fwet_out if not found in restart
    !Config If    = TOPM_calcul
    !Config Def   = 0.0
    !Config Help  = The initial value of fwet_out if its value 
    !Config         is not found in the restart file. This should only be used if
    !Config         the model is started without a restart file. 
    !Config Units =
    CALL setvar_p (fwet_out, val_exp,'HYDROL_FWET', zero)
 
    ! Read input file
    CALL topmodel_io(kjpindex, lalo)
  
    !le deficit utilise pour TOPMODEL va etre calcule par rapport a la saturation
    !ZM(:)=(ZWFC(:)-ZWWILT(:))*ZD_TOP(:)/4.

    !ZM(:) = (mcs du grid_cell - mcw du grid_cell)*zmaxh/4.
    mcs_grid(:) = mcs(1)*soiltile(:,1)+mcs(2)*soiltile(:,2)+mcs(3)*soiltile(:,3)
    mcw_grid(:) = mcw(1)*soiltile(:,1)+mcw(2)*soiltile(:,2)+mcw(3)*soiltile(:,3)
    ZM(:) = ( mcs_grid(:) -  mcw_grid(:) )*zmaxh/4.


    !2 obtention des differentes fonctions necessaires a TOPMODEL en chaque grid-cell  
    CALL topmodel_init_main(kjpindex, lalo, veget_max, mcw_grid,mcs_grid,zmaxh, ZM,ZMIN, ZMAX, &
        & ZMEAN, ZSTDT, ZSKEW, ZTAB_FSAT, ZTAB_WTOP, ZTAB_FWET, ZTAB_WTOP_WET, ZZPAS)
          
  END SUBROUTINE topmodel_initialize


!! ================================================================================================================================
!! SUBROUTINE 	: topmodel_main
!!
!>\BRIEF        Main part of top model 
!!
!! DESCRIPTION :
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE topmodel_main(kjpindex, dz, humtot, profil_froz_hydro, zmaxh)

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables
  
    INTEGER(i_std), INTENT(in)  :: kjpindex         !! Domain size
    REAL(r_std),  DIMENSION(kjpindex), INTENT(in)  :: humtot !! Total Soil Moisture @tex $(kg m^{-2})$ @endtex
    REAL(r_std), DIMENSION(kjpindex,nslm), INTENT(in) :: profil_froz_hydro !! Frozen fraction for each hydrological soil layer
    REAL(r_std), INTENT(in) :: zmaxh !! Max hydrol depth 
    REAL(r_std), DIMENSION(nslm), INTENT(in)  :: dz               !! Internode thickness [dnh in vertical_soil] transformed into (mm)
 
    CALL topmodel_subgrid_main(kjpindex, ZTAB_FSAT, ZTAB_WTOP, humtot, profil_froz_hydro, fsat,&
        & ZTAB_FWET,ZTAB_WTOP_WET,fwet, zmaxh, &
        & 1000*(mcs_grid(:)-mcw_grid(:)), fwt1, fwt2, fwt3, fwt4, ZM, ZMIN, ZMAX, ZZPAS, dz)
    fwet_out(:) = fwet(:)

  END SUBROUTINE topmodel_main

!! ================================================================================================================================
!! SUBROUTINE 	: topmodel_history
!!
!>\BRIEF         Send data to history files
!!
!! DESCRIPTION :
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE topmodel_history(kjit, kjpindex, index, hist_id)

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables
  
    INTEGER(i_std), INTENT(in)                          :: kjit      !! Time step number 
    INTEGER(i_std), INTENT(in)                          :: kjpindex  !! Domain size
    INTEGER(i_std), DIMENSION (kjpindex), INTENT (in)   :: index     !! Indeces of the points on the map
    INTEGER(i_std),INTENT (in)                          :: hist_id   !!  _history_ file identifier

    CALL histwrite_p(hist_id, 'fsat', kjit, fsat, kjpindex, index)
    CALL histwrite_p(hist_id, 'fwet', kjit, fwet, kjpindex, index)
    CALL histwrite_p(hist_id, 'fwt1', kjit, fwt1, kjpindex, index)
    CALL histwrite_p(hist_id, 'fwt2', kjit, fwt2, kjpindex, index)
    CALL histwrite_p(hist_id, 'fwt3', kjit, fwt3, kjpindex, index)
    CALL histwrite_p(hist_id, 'fwt4', kjit, fwt4, kjpindex, index)
    CALL histwrite_p(hist_id, 'ZMIN', kjit, ZMIN, kjpindex, index)
    CALL histwrite_p(hist_id, 'ZMAX', kjit, ZMAX, kjpindex, index)
    CALL histwrite_p(hist_id, 'ZMEAN', kjit, ZMEAN, kjpindex, index)
    !CALL histwrite_p(hist_id, 'NB_PIXE', kjit, NB_PIXE, kjpindex, index)
    CALL histwrite_p(hist_id, 'ZSTDT', kjit, ZSTDT, kjpindex, index)
    CALL histwrite_p(hist_id, 'ZSKEW', kjit, ZSKEW, kjpindex, index)
!       CALL histwrite_p(hist_id, 'dsg', kjit, dsg, kjpindex*nvm, indexveg)
!       CALL histwrite_p(hist_id, 'dsp', kjit, dsp, kjpindex*nvm, indexveg)
!       CALL histwrite_p(hist_id, 'ZWSAT', kjit, ZWSAT, kjpindex, index)
!       CALL histwrite_p(hist_id, 'ZWWILT', kjit, ZWWILT, kjpindex, index) 
!       CALL histwrite_p(hist_id, 'ZWFC', kjit, ZWFC, kjpindex, index) 
!       CALL histwrite_p(hist_id, 'RU', kjit, ruu_ch, kjpindex, index) 
!       CALL histwrite_p(hist_id, 'mx_eau_var', kjit, mx_eau_var, kjpindex, index)
    CALL histwrite_p(hist_id, 'drunoff_tot', kjit, drunoff_tot, kjpindex, index)

  END SUBROUTINE topmodel_history

!! ================================================================================================================================
!! SUBROUTINE 	: topmodel_alloc
!!
!>\BRIEF        Module scope allocation arrays
!!
!! DESCRIPTION :
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE topmodel_alloc(kjpindex)
    !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)  :: kjpindex !! Domain size
    !! 0.4 Local variables
    INTEGER(i_std) :: ier

    
    ALLOCATE (fsat(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable fsat','','')

    ALLOCATE (fwet(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable fwet','','')

    ALLOCATE (fwt1(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable fwt1','','')
    
    ALLOCATE (fwt2(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable fwt2','','')
   
    ALLOCATE (fwt3(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable fwt3','','')

    ALLOCATE (fwt4(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable fwt4','','')

    ALLOCATE (drunoff(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable drunoff','','')
       
    ALLOCATE (ZMEAN(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable ZMEAN','','')

!    ALLOCATE (NB_PIXE(kjpindex),stat=ier)
!    IF (ier.NE.0) THEN
!        WRITE (numout,*) ' error in mx_eau_var allocation. We stop. We need kjpindex words = ',kjpindex
!        STOP 'hydrolc_init'
!    END IF
    ALLOCATE (ZSTDT(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable ZSTDT','','')
    
    ALLOCATE (ZSKEW(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable ZSKEW','','')

    ALLOCATE (ZMIN(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable ZMIN','','')
    
    ALLOCATE (ZMAX(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable ZMAX','','')
    
    ALLOCATE (ZM(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable ZM','','')

    ALLOCATE (ZZPAS(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable ZZPAS','','')

    ALLOCATE (ZTAB_FSAT(kjpindex,1000),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable ZTAB_FSAT','','')

    ALLOCATE (ZTAB_WTOP(kjpindex,1000),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable ZTAB_WTOP','','')

    ALLOCATE (ZTAB_FWET(kjpindex,1000),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable ZTAB_FWET','','')

    ALLOCATE (ZTAB_WTOP_WET(kjpindex,1000),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable ZTAB_WTOP_WET','','')

    ALLOCATE (mcw_grid(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mcw_grid','','')

    ALLOCATE (mcs_grid(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mcs_grid','','')

    ALLOCATE (fwet_out(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable fwet_out','Error code:',ier)
    fwet_out(:) = undef_sechiba

    ALLOCATE (drunoff_tot(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for drunoff_tot','','')
    drunoff_tot(:) = zero

  END SUBROUTINE topmodel_alloc

!! ================================================================================================================================
!! SUBROUTINE 	: topmodel_clear
!!
!>\BRIEF         Deallocate all module variables
!!
!! DESCRIPTION :
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE topmodel_clear()
    !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    !! 0.4 Local variables

    IF (ALLOCATED  (fsat))  DEALLOCATE (fsat)
    IF (ALLOCATED  (fwet))  DEALLOCATE (fwet)
    IF (ALLOCATED  (fwt1))  DEALLOCATE (fwt1)
    IF (ALLOCATED  (fwt2))  DEALLOCATE (fwt2)
    IF (ALLOCATED  (fwt3))  DEALLOCATE (fwt3)
    IF (ALLOCATED  (fwt4))  DEALLOCATE (fwt4)
    IF (ALLOCATED  (drunoff))  DEALLOCATE (drunoff)
    IF (ALLOCATED  (ZMEAN)) DEALLOCATE (ZMEAN)
!    IF (ALLOCATED  (NB_PIXE)) DEALLOCATE (NB_PIXE)
    IF (ALLOCATED  (ZSTDT)) DEALLOCATE (ZSTDT)
    IF (ALLOCATED  (ZSKEW)) DEALLOCATE (ZSKEW)
    IF (ALLOCATED  (ZMIN)) DEALLOCATE (ZMIN)
    IF (ALLOCATED  (ZMAX)) DEALLOCATE (ZMAX)
    IF (ALLOCATED  (ZM)) DEALLOCATE (ZM)
    IF (ALLOCATED  (ZZPAS)) DEALLOCATE (ZZPAS)
    IF (ALLOCATED  (ZTAB_FSAT)) DEALLOCATE (ZTAB_FSAT)
    IF (ALLOCATED  (ZTAB_WTOP)) DEALLOCATE (ZTAB_WTOP)
    IF (ALLOCATED  (ZTAB_FWET)) DEALLOCATE (ZTAB_FWET)
    IF (ALLOCATED  (ZTAB_WTOP_WET)) DEALLOCATE (ZTAB_WTOP_WET)

    IF ( ALLOCATED (mcs_grid)) DEALLOCATE (mcs_grid)
    IF ( ALLOCATED (mcw_grid)) DEALLOCATE (mcw_grid)

    IF ( ALLOCATED (fwet_out)) DEALLOCATE (fwet_out)
    IF ( ALLOCATED (drunoff_tot)) DEALLOCATE (drunoff_tot)

  END SUBROUTINE topmodel_clear

!! ================================================================================================================================
!! SUBROUTINE 	: topmodel_readrestart
!!
!>\BRIEF         Read data from the restart files
!!
!! DESCRIPTION :
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE topmodel_readrestart(kjit, kjpindex, nbp_glo, rest_id, index_g)

    !! 0. Variable and parameter declaration
    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                         :: kjit             !! Time step number 
    INTEGER(i_std), INTENT(in)                         :: kjpindex         !! Domain size
    INTEGER(i_std), INTENT(in)                         :: nbp_glo          !!
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)   :: index_g          !! Indeces of the points on the map
    INTEGER(i_std),INTENT (in)                         :: rest_id          !! Restart file identifier

    CALL ioconf_setatt('UNITS', '-')
    CALL ioconf_setatt('LONG_NAME','fwet pr autres routines')
    CALL restget_p (rest_id, 'fwet_out', nbp_glo, 1  , 1, kjit, .TRUE.,fwet_out , "gather", nbp_glo, index_g)

  END SUBROUTINE topmodel_readrestart

!! ================================================================================================================================
!! SUBROUTINE 	: topmodel_writerestart
!!
!>\BRIEF        Write data to restart files 
!!
!! DESCRIPTION :
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE topmodel_writerestart(kjit, kjpindex, nbp_glo, rest_id, index_g)

    !! 0. Variable and parameter declaration
    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                         :: kjit             !! Time step number 
    INTEGER(i_std), INTENT(in)                         :: kjpindex         !! Domain size
    INTEGER(i_std), INTENT(in)                         :: nbp_glo          !!
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)   :: index_g          !! Indeces of the points on the map
    INTEGER(i_std),INTENT (in)                         :: rest_id          !! Restart file identifier

    CALL restput_p(rest_id, 'fwet_out', nbp_glo,   1, 1, kjit,  fwet_out, 'scatter',  nbp_glo, index_g)
  END SUBROUTINE topmodel_writerestart


!! ================================================================================================================================
!! SUBROUTINE 	: topmodel_io
!!
!>\BRIEF         Read input file
!!
!! DESCRIPTION :
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE topmodel_io(kjpindex, lalo)
    INTEGER(i_std), INTENT(in)                         :: kjpindex         !! Domain size
    REAL(r_std), INTENT(in)                       :: lalo(kjpindex,2)          !! Vector of latitude and longitudes (degree)        


    CHARACTER(LEN=80)                                    :: filename       !! To store file names for I/O
    INTEGER(i_std)                                       :: iml, jml, lml, tml, fid
    REAL(r_std),ALLOCATABLE,DIMENSION(:,:)               :: Zminf, Zmaxf, Zmeanf, Zstdf, Zskewf
    REAL(r_std),ALLOCATABLE,DIMENSION(:)                 :: lon_temp, lat_temp
    INTEGER(i_std)                                       :: pssitau(1)
    REAL(r_std)                                          :: lev(1), pssdate, pssdt
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)             :: lat_rel, lon_rel
    INTEGER(i_std) :: ip, ix, iy, imin, jmin, ier
    REAL(r_std) :: dlon, dlonmin, dlat, dlatmin

    CHARACTER(LEN=40) :: subrname 

    subrname ="topmodel_io" 
    !  Needs to be a configurable variable
    !
    !
    !Config Key   = TOPMODEL_PARAMETERS_FILE
    !Config Desc  = Name of file from which TOPMODEL parameters file are read
    !Config Def   = TOPMODEL_param_1deg.nc
    !Config If    = TOPM_CALCUL and NOT(IMPOSE_VEG)
    !Config Help  = The name of the file to be opened to read the TOPMODEL parameters. 
    !Config         
    !Config Units = [FILE]
    !
    filename = 'TOPMODEL_param_1deg.nc'
    CALL getin_p('TOPMODEL_PARAMETERS_FILE',filename)
    !
    IF (is_root_prc) THEN
       CALL flininfo(filename,iml, jml, lml, tml, fid)
       CALL flinclo(fid)
    ENDIF
    CALL bcast(iml)
    CALL bcast(jml)
    CALL bcast(lml)
    CALL bcast(tml)
    !
    ! soils_param.nc file is 1 soit texture file.
    !
    ALLOCATE(lat_rel(iml,jml), STAT=ier)
    IF (ier /= 0) CALL ipslerr_p(3, subrname,'Memory error allocation on variable lat_rel','Error code:',ier)

    ALLOCATE(lon_rel(iml,jml), STAT=ier)
    IF (ier /= 0) CALL ipslerr_p(3, subrname,'Memory error allocation on variable lon_rel','Error code:',ier)

    ALLOCATE(Zminf(iml,jml), STAT=ier)
    IF (ier /= 0) CALL ipslerr_p(3, subrname,'Memory error allocation on variable Zminf','Error code:',ier)
    ALLOCATE(Zmaxf(iml,jml), STAT=ier)
    IF (ier /= 0) CALL ipslerr_p(3, subrname,'Memory error allocation on variable Zmaxf','Error code:',ier)
    ALLOCATE(Zmeanf(iml,jml), STAT=ier)
    IF (ier /= 0) CALL ipslerr_p(3, subrname,'Memory error allocation on variable Zmeanf','Error code:',ier)
    ALLOCATE(Zstdf(iml,jml), STAT=ier)
    IF (ier /= 0) CALL ipslerr_p(3, subrname,'Memory error allocation on variable Zstdf','Error code:',ier)
    ALLOCATE(Zskewf(iml,jml), STAT=ier)
    IF (ier /= 0) CALL ipslerr_p(3, subrname,'Memory error allocation on variable Zskewf','Error code:',ier)

    ALLOCATE (lon_temp(iml),lat_temp(jml), STAT=ier)
    IF (ier /= 0) CALL ipslerr_p(3, subrname,'Memory error allocation on variable lon_temp','Error code:',ier)
    !
    IF (is_root_prc) CALL flinopen(filename, .FALSE., iml, jml, lml, lon_rel, lat_rel, lev, tml, pssitau, pssdate, pssdt, fid)
    CALL bcast(lon_rel)
    CALL bcast(lat_rel)
    CALL bcast(pssitau)
    CALL bcast(pssdate)
    CALL bcast(pssdt)

    !
    IF (is_root_prc) CALL flinget(fid, 'Zmin', iml, jml, lml, tml, 1, 1, Zminf)
    IF (is_root_prc) CALL flinget(fid, 'Zmax', iml, jml, lml, tml, 1, 1, Zmaxf)
    IF (is_root_prc) CALL flinget(fid, 'Zmean', iml, jml, lml, tml, 1, 1, Zmeanf)
    IF (is_root_prc) CALL flinget(fid, 'Zstdev', iml, jml, lml, tml, 1, 1, Zstdf)
    IF (is_root_prc) CALL flinget(fid, 'Zskew', iml, jml, lml, tml, 1, 1, Zskewf)

    CALL bcast(Zminf)
    CALL bcast(Zmaxf)
    CALL bcast(Zmeanf)
    CALL bcast(Zstdf)
    CALL bcast(Zskewf)
    !
    IF (is_root_prc) CALL flinclo(fid)

!!!! TOPMODEL parameters 2D into 1D 
    lon_temp(:) = lon_rel(:,1)
    lat_temp(:) = lat_rel(1,:)

    DO ip = 1, kjpindex
          dlonmin = HUGE(1.)
          DO ix = 1,iml
             dlon = MIN( ABS(lalo(ip,2)-lon_temp(ix)), ABS(lalo(ip,2)+360.-lon_temp(ix)), ABS(lalo(ip,2)-360.-lon_temp(ix)) )
             IF ( dlon .LT. dlonmin ) THEN
                imin = ix
                dlonmin = dlon
             ENDIF
          ENDDO
          dlatmin = HUGE(1.)
          DO iy = 1,jml
             dlat = ABS(lalo(ip,1)-lat_temp(iy))
             IF ( dlat .LT. dlatmin ) THEN
                jmin = iy
                dlatmin = dlat
             ENDIF
          ENDDO
          ZMIN(ip) = Zminf(imin,jmin)
          ZMAX(ip) = Zmaxf(imin,jmin)
          ZMEAN(ip) = Zmeanf(imin,jmin)
          ZSTDT(ip) = Zstdf(imin,jmin)
          ZSKEW(ip) = Zskewf(imin,jmin)
    ENDDO

    DEALLOCATE (lon_temp)
    DEALLOCATE (lat_temp)
    DEALLOCATE (Zminf)
    DEALLOCATE (Zmaxf)
    DEALLOCATE (Zmeanf)
    DEALLOCATE (Zstdf)
    DEALLOCATE (Zskewf)
     
    WRITE (numout,*) 'STATS CTI OK num1!'
    WRITE (numout,*) 'psstest2'
  END SUBROUTINE topmodel_io


!! ================================================================================================================================
!! SUBROUTINE 	: topmodel_histdef
!!
!>\BRIEF         Read input file
!!
!! DESCRIPTION :
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE topmodel_histdef(iim, jjm, dt, hist_id, hori_id, dw, avescatter, fluxop )
    INTEGER(i_std), INTENT(in)     :: iim, jjm 
    REAL(r_std), INTENT(in)        :: dt        !! Time step of the counter in seconds
    INTEGER(i_std), INTENT(in)     :: hist_id !! History file identification for SECHIBA
    INTEGER(i_std), INTENT(in)     :: hori_id !! ID of the default horizontal longitude and latitude map.
    REAL(r_std), INTENT(in)        :: dw      !! Frequency of history write (sec.)
    CHARACTER(LEN=40),DIMENSION(:) :: avescatter, fluxop


    CALL histdef(hist_id, 'fsat', 'Fraction of soil saturated', '-',  &
          & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(4), dt,dw)
    CALL histdef(hist_id, 'fwet', 'Fraction of wetland', '-',  &
          & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(4), dt,dw)
    CALL histdef(hist_id, 'fwt1', 'Fraction of soil with wt [0,xcm]  & ', '-',  &
          & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(4), dt,dw)
    CALL histdef(hist_id, 'fwt2', 'Fraction of soil with wt [xcm,ycm] ', '-',  &
          & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(4), dt,dw)
    CALL histdef(hist_id, 'fwt3', 'Fraction of soil with wt [ycm,zcm] ', '-',  &
          & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(4), dt,dw)
    CALL histdef(hist_id, 'fwt4', 'Fraction of soil with wt [ucm,vcm] ', '-',  &
          & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(4), dt,dw)

    CALL histdef(hist_id, 'ZMIN', 'MIN INDICE TOPO', '-', &
          & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(3), dt,dw)      
    CALL histdef(hist_id, 'ZMAX', 'MAX INDICE TOPO', '-', &
          & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(3), dt,dw)      
    CALL histdef(hist_id, 'ZMEAN', 'MEAN INDICE TOPO', '-', &
          & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(3), dt,dw)      
! 	         CALL histdef(hist_id, 'NB_PIXE', 'NB PIXELS AVC VALEUR', '-', &
!	              & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(3), dt,dw)      
    CALL histdef(hist_id, 'ZSTDT', 'STD INDICE TOPO', '-', &
          & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(3), dt,dw)      
    CALL histdef(hist_id, 'ZSKEW', 'SKEWNESS INDICE TOPO', '-', &
          & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(3), dt,dw)      
    CALL histdef(hist_id, 'drunoff_tot', 'Surface drunoff', 'mm/d', &
          & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop(2), dt,dw) 

  END SUBROUTINE topmodel_histdef

!! ================================================================================================================================
!! SUBROUTINE 	: topmodel_parameters
!!
!>\BRIEF         Read input file
!!
!! DESCRIPTION :
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE topmodel_parameters()

!pss+
    !
    !Config Key   = SHIFT_FSAT_FWET
    !Config Desc  = shift saturation fraction to wetland fraction
    !Config If    = TOPM_calcul
    !Config Def   = 5.
    !Config Help  = 
    !Config Units = [-
    CALL getin_p('SHIFT_FSAT_FWET',SHIFT_fsat_fwet)
    !Config Key   = WTD1_borne 
    !Config Desc  = depth of subsurface saturation for wetland1 fraction
    !Config If    = TOPM_calcul 
    !Config Def   = 0.06
    !Config Help  = 
    !Config Units = [m]
    CALL getin_p('WTD1_BORNE',WTD1_borne)
    !Config Key   = WTD2_borne 
    !Config Desc  = depth of subsurface saturation for wetland1 fraction
    !Config If    = TOPM_calcul 
    !Config Def   = 0.12
    !Config Help  = 
    !Config Units = [m]
    CALL getin_p('WTD2_BORNE',WTD2_borne)
    !Config Key   = WTD3_borne 
    !Config Desc  = depth of subsurface saturation for wetland1 fraction
    !Config If    = TOPM_calcul 
    !Config Def   = 0.18
    !Config Help  = 
    !Config Units = [m]
    CALL getin_p('WTD3_BORNE',WTD3_borne)
    !Config Key   = WTD4_borne 
    !Config Desc  = depth of subsurface saturation for wetland1 fraction
    !Config If    = TOPM_calcul 
    !Config Def   = 0.24
    !Config Help  = 
    !Config Units = [m]
    CALL getin_p('WTD4_BORNE',WTD4_borne)
!pss-
  END SUBROUTINE topmodel_parameters
  
END MODULE topmodel
