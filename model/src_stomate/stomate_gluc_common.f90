! =================================================================================================================================
! MODULE       : stomate_gluc_common
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       This module contains common fuctions and subroutines used by
!              gross land use change and forestry harvest modules.
!!
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): 
!!
!! REFERENCE(S)	: None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/perso/albert.jornet/ORCHIDEE-MICT/src_stomate/stomate_lcchange.f90 $
!! $Date: 2015-07-30 15:38:45 +0200 (Thu, 30 Jul 2015) $
!! $Revision: 2847 $
!! \n
!_ ================================================================================================================================


MODULE stomate_gluc_common

  ! modules used:
  
  USE ioipsl_para
  USE stomate_data
  USE pft_parameters
  USE constantes
  USE constantes_soil_var
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC calc_cover, cross_give_receive, &
         initialize_proxy_pft, sap_take, collect_legacy_pft, &
         collect_legacy_pft_forestry, &
         add_incoming_proxy_pft, empty_pft, build_age_index, &
         prepare_balance_check, luc_balance_check
 
  
CONTAINS

!! ================================================================================================================================
!! SUBROUTINE 	: build_age_index
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
  
  SUBROUTINE build_age_index

    INTEGER, ALLOCATABLE, SAVE                  :: indall_tree(:)       !! Indices for all tree PFTs
    INTEGER, ALLOCATABLE, SAVE                  :: indold_tree(:)       !! Indices for old tree cohort only
    INTEGER, ALLOCATABLE, SAVE                  :: indagec_tree(:,:)    !! Indices for secondary tree cohorts, 
                                                                        !! note the sequence is old->young.
    INTEGER, ALLOCATABLE, SAVE                  :: indall_grass(:)      !! Indices for all grass PFTs
    INTEGER, ALLOCATABLE, SAVE                  :: indold_grass(:)      !! Indices for old grasses only
    INTEGER, ALLOCATABLE, SAVE                  :: indagec_grass(:,:)   !! Indices for secondary grass cohorts
                                                                        !! note the sequence is old->young.
    INTEGER, ALLOCATABLE, SAVE                  :: indall_pasture(:)    !! Indices for all pasture PFTs
    INTEGER, ALLOCATABLE, SAVE                  :: indold_pasture(:)    !! Indices for old pasture only
    INTEGER, ALLOCATABLE, SAVE                  :: indagec_pasture(:,:) !! Indices for secondary pasture cohorts
                                                                        !! note the sequence is old->young.
    INTEGER, ALLOCATABLE, SAVE                  :: indall_crop(:)       !! Indices for all crop PFTs
    INTEGER, ALLOCATABLE, SAVE                  :: indold_crop(:)       !! Indices for old crops only
    INTEGER, ALLOCATABLE, SAVE                  :: indagec_crop(:,:)    !! Indices for secondary crop cohorts
                                                                        !! note the sequence is old->young.
    INTEGER :: num_tree_mulagec,num_grass_mulagec,     &
               num_pasture_mulagec,num_crop_mulagec, &
               itree,itree2,igrass,igrass2,ipasture,ipasture2,icrop,icrop2
    INTEGER :: i,j,ivma,staind,endind,ivm


    !! 1. We first build all different indices that we are going to use
    !!    in handling the PFT exchanges, three types of indices are built:
    !!     - for all age classes
    !!     - include only oldest age classes
    !!     - include all age classes excpet the oldest ones
    ! We have to build these indices because we would like to extract from
    ! donating PFTs in the sequnce of old->young age classes or the revserse, 
    ! and add in the receving PFTs only in the youngest-age-class PFTs. These 
    ! indicies allow us to know where the different age classes are.

    ! calculate the total number of MTCs for each vegetation type.
    num_tree_mulagec=0          
    num_grass_mulagec=0
    num_pasture_mulagec=0
    num_crop_mulagec=0
    
    !! 1.1 Calculate the number of PFTs for different MTCs and allocate
    !! the old and all indices arrays.

    ! [Note here the sequence to identify tree,pasture,grass,crop] is
    ! critical. The similar sequence is used in the subroutine "calc_cover".
    ! Do not forget to change the sequence there if you modify here.
    DO ivma =2,nvmap
      staind=start_index(ivma)
      IF (nagec_pft(ivma)==1) THEN
        WRITE(numout,*) "Error: metaclass has only a single age group: ",ivma
        STOP
      ELSE
        IF (is_tree(staind)) THEN
          num_tree_mulagec = num_tree_mulagec+1
        ELSE IF (is_grassland_manag(staind)) THEN
          num_pasture_mulagec = num_pasture_mulagec+1
        ELSE IF (natural(staind)) THEN
          num_grass_mulagec = num_grass_mulagec+1
        ELSE
          num_crop_mulagec = num_crop_mulagec+1
        ENDIF
      ENDIF
    ENDDO
    
    !! Allocate index array
    ! allocate all index
    ALLOCATE(indall_tree(num_tree_mulagec*nagec_tree))     
    ALLOCATE(indall_grass(num_grass_mulagec*nagec_herb))     
    ALLOCATE(indall_pasture(num_pasture_mulagec*nagec_herb))     
    ALLOCATE(indall_crop(num_crop_mulagec*nagec_herb))     

    ! allocate old-ageclass index
    ALLOCATE(indold_tree(num_tree_mulagec))     
    ALLOCATE(indold_grass(num_grass_mulagec))     
    ALLOCATE(indold_pasture(num_pasture_mulagec))     
    ALLOCATE(indold_crop(num_crop_mulagec))     

    !! 1.2 Fill the oldest-age-class and all index arrays
    itree=0
    igrass=0
    ipasture=0
    icrop=0
    itree2=1
    igrass2=1
    ipasture2=1
    icrop2=1
    DO ivma =2,nvmap
      staind=start_index(ivma)
      IF (is_tree(staind)) THEN
        itree=itree+1
        indold_tree(itree) = staind+nagec_pft(ivma)-1
        DO j = 0,nagec_pft(ivma)-1
          indall_tree(itree2+j) = staind+j
        ENDDO
        itree2=itree2+nagec_pft(ivma)
      ELSE IF (natural(staind) .AND. .NOT. is_grassland_manag(staind)) THEN
        igrass=igrass+1
        indold_grass(igrass) = staind+nagec_pft(ivma)-1
        DO j = 0,nagec_pft(ivma)-1
          indall_grass(igrass2+j) = staind+j
        ENDDO
        igrass2=igrass2+nagec_pft(ivma)
      ELSE IF (is_grassland_manag(staind)) THEN
        ipasture = ipasture+1
        indold_pasture(ipasture) = staind+nagec_pft(ivma)-1
        DO j = 0,nagec_pft(ivma)-1
          indall_pasture(ipasture2+j) = staind+j
        ENDDO
        ipasture2=ipasture2+nagec_pft(ivma)
      ELSE
        icrop = icrop+1
        indold_crop(icrop) = staind+nagec_pft(ivma)-1
        DO j = 0,nagec_pft(ivma)-1
          indall_crop(icrop2+j) = staind+j
        ENDDO
        icrop2=icrop2+nagec_pft(ivma)
      ENDIF
    ENDDO
    
    !! 1.3 Allocate and fill other age class index

    ALLOCATE(indagec_tree(num_tree_mulagec,nagec_tree-1))     
    ALLOCATE(indagec_grass(num_grass_mulagec,nagec_herb-1))     
    ALLOCATE(indagec_pasture(num_pasture_mulagec,nagec_herb-1))
    ALLOCATE(indagec_crop(num_crop_mulagec,nagec_herb-1))

    ! fill the non-oldest age class index arrays when number of age classes
    ! is more than 1.
    itree=0
    igrass=0
    ipasture=0
    icrop=0
    DO ivma = 2,nvmap
      staind=start_index(ivma)
      IF (nagec_pft(ivma) > 1) THEN
        IF (is_tree(staind)) THEN
          itree=itree+1
          DO j = 1,nagec_tree-1
            indagec_tree(itree,j) = staind+nagec_tree-j-1
          ENDDO
        ELSE IF (natural(staind) .AND. .NOT. is_grassland_manag(staind)) THEN
          igrass=igrass+1
          DO j = 1,nagec_herb-1
            indagec_grass(igrass,j) = staind+nagec_herb-j-1
          ENDDO
        ELSE IF (is_grassland_manag(staind)) THEN
          ipasture=ipasture+1
          DO j = 1,nagec_herb-1
            indagec_pasture(ipasture,j) = staind+nagec_herb-j-1
          ENDDO
        ELSE
          icrop=icrop+1
          DO j = 1,nagec_herb-1
            indagec_crop(icrop,j) = staind+nagec_herb-j-1
          ENDDO
        ENDIF
      ENDIF
    ENDDO

  END SUBROUTINE build_age_index

! ================================================================================================================================
!! SUBROUTINE   : calc_cover
!!
!>\BRIEF        Calculate coverage fraction for different age classes of forest,
!!              grass, pasture and crops and also for each metaclass. Note baresoil is excluded. 
!!              
!! DESCRIPTION :
!! Note:
!! 1. "calc_cover" subroutine does not depend on how many age classes
!! there are in each MTC.
!! 2. Fraction of baresoil is excluded here. This means transformation
!! of baresoil to a vegetated PFT is excluded in gross land cover change.
!!  
!!
!! MAIN OUTPUT VARIABLE(S) :  
!!
!! \n
!_ ================================================================================================================================
  SUBROUTINE calc_cover(npts,veget_max,veget_mtc,vegagec_tree,vegagec_grass, &
                 vegagec_pasture,vegagec_crop)

   
    IMPLICIT NONE

    !! Input variables
    INTEGER, INTENT(in)                                       :: npts             !! Domain size - number of pixels (unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)         :: veget_max           !! "maximal" coverage fraction of a PFT on the ground

    !! Output variables
    REAL(r_std), DIMENSION(npts,nvmap), INTENT(inout)         :: veget_mtc        !! "maximal" coverage fraction of a PFT on the ground
    REAL(r_std), DIMENSION(npts,nagec_tree), INTENT(inout)    :: vegagec_tree     !! fraction of tree age-class groups, in sequence of old->young
    REAL(r_std), DIMENSION(npts,nagec_herb), INTENT(inout)    :: vegagec_grass    !! fraction of grass age-class groups, in sequence of old->young
    REAL(r_std), DIMENSION(npts,nagec_herb), INTENT(inout)    :: vegagec_pasture  !! fraction of pasture age-class groups, in sequence of old->young
    REAL(r_std), DIMENSION(npts,nagec_herb), INTENT(inout)    :: vegagec_crop     !! fraction of crop age-class groups, in sequence of old->young

    !! Local variables
    INTEGER(i_std)                                          :: ivma,staind,endind,j    !! indices (unitless)

    veget_mtc(:,:) = 0.
    vegagec_tree(:,:) = 0.
    vegagec_grass(:,:) = 0.
    vegagec_pasture(:,:) = 0.
    vegagec_crop(:,:) = 0.

    ! Calculate veget_max for MTCs
    DO ivma = 1,nvmap
      staind = start_index(ivma)
      IF (nagec_pft(ivma) == 1) THEN
        veget_mtc(:,ivma) = veget_max(:,staind)
      ELSE
        veget_mtc(:,ivma) = \
          SUM(veget_max(:,staind:staind+nagec_pft(ivma)-1),DIM=2)
      ENDIF
    ENDDO

    ! Calculate veget_max for each age class
    DO ivma = 2,nvmap  !here we start with 2 to exclude baresoil (always PFT1)
      staind = start_index(ivma)
      endind = staind+nagec_pft(ivma)-1

      ! Single-age-class MTC goest to oldest age class.
      IF (nagec_pft(ivma) == 1) THEN
        IF (is_tree(staind)) THEN
          vegagec_tree(:,1) = vegagec_tree(:,1)+veget_max(:,staind)
        ELSE IF (is_grassland_manag(staind)) THEN
          vegagec_pasture(:,1) = vegagec_pasture(:,1)+veget_max(:,staind)
        ELSE IF (natural(staind)) THEN
          vegagec_grass(:,1) = vegagec_grass(:,1)+veget_max(:,staind)
        ELSE
          vegagec_crop(:,1) = vegagec_crop(:,1)+veget_max(:,staind)
        ENDIF

      ELSE
        IF (is_tree(staind)) THEN
          DO j=1,nagec_tree
            vegagec_tree(:,j) = vegagec_tree(:,j)+veget_max(:,endind-j+1)
          ENDDO
        ELSE IF (is_grassland_manag(staind)) THEN
          DO j=1,nagec_herb
            vegagec_pasture(:,j) = vegagec_pasture(:,j)+veget_max(:,endind-j+1)
          ENDDO
        ELSE IF (natural(staind)) THEN
          DO j=1,nagec_herb
            vegagec_grass(:,j) = vegagec_grass(:,j)+veget_max(:,endind-j+1)
          ENDDO
        ELSE
          DO j=1,nagec_herb
            vegagec_crop(:,j) = vegagec_crop(:,j)+veget_max(:,endind-j+1)
          ENDDO
        ENDIF
      ENDIF
    ENDDO

  END SUBROUTINE calc_cover


! ================================================================================================================================
!! SUBROUTINE   : cross_give_receive
!!
!>\BRIEF        : Allocate the outgoing and receving fractions in respective
!!                PFTs.
!! \n
!! Notes:
!!  1. veget_max is subtracted when fractions are taken out, but newly added
!!     fractions in the youngest age class is not added, to avoid this newly
!!     created fractions being used again the following transitions. This is 
!!     is reasonable because the newly created youngest-age-class PFT fractions
!!     have nothing but small sapling biomass and it's unreasonable to use it
!!     for any further land use conversion activities.
!_ ================================================================================================================================
  SUBROUTINE cross_give_receive(ipts,frac_used,veget_mtc,newvegfrac,           &
                     indold_tree,indagec_crop,nagec_receive,num_crop_mulagec, &
                     veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp)


    IMPLICIT NONE

    !! 0. Input variables
    INTEGER, INTENT(in)                             :: ipts
    REAL(r_std), INTENT(in)                         :: frac_used                 !! fraction that the giving PFTs are going to collectively give
    REAL(r_std), DIMENSION(:,:), INTENT(in)         :: veget_mtc            !! "maximal" coverage fraction of a PFT on the ground
    INTEGER, DIMENSION(:), INTENT(in)               :: indold_tree          !! Indices for PFTs giving out fractions; 
                                                                            !! here use old tree cohort as an example
    INTEGER, DIMENSION(:,:), INTENT(in)             :: indagec_crop         !! Indices for secondary basic-vegetation cohorts; The youngest age classes
                                                                            !! of these vegetations are going to receive fractions. 
                                                                            !! here we use crop cohorts as an example
    INTEGER, INTENT(in)                             :: num_crop_mulagec     !! number of crop MTCs with more than one age classes
    INTEGER, INTENT(in)                             :: nagec_receive        !! number of age classes in the receiving basic types
                                                                            !! (i.e., tree, grass, pasture, crop), here we can use crop
                                                                            !! as an example, nagec_receive=nagec_herb
    REAL(r_std), DIMENSION(:,:),INTENT(in)          :: newvegfrac           !! 

    !! 1. Modified variables
    REAL(r_std), DIMENSION(:,:), INTENT(inout)      :: veget_max            !! "maximal" coverage fraction of a PFT on the ground
    REAL(r_std), DIMENSION(:,:), INTENT(inout)      :: glcc_pft             !! a temporary variable to hold the fractions each PFT is going to lose
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)    :: glcc_pftmtc          !! a temporary variable to hold the fraction of ipft->ivma, i.e., from 
                                                                            !! PFT_{ipft} to the youngest age class of MTC_{ivma}
    REAL(r_std), DIMENSION(:,:), INTENT(inout)      :: glcc_pft_tmp         !! a temporary variable to hold the fractions each PFT is going to lose

    !! Local vriables
    INTEGER  :: j,ipft, iyoung
    REAL(r_std) :: totalveg, tveg_now
    LOGICAL :: guide_newfrac
    
    guide_newfrac = gluc_newfrac_guide

    ! Out final objective is to know glcc_pftmtc, i.e., the fraction from each PFT
    ! to the youngest age group of each MTC. We separate this task into two steps:
    ! 1. we allocate the total outgoing fraction into the same age-class PFTs of 
    ! the a basic-vegetation (for example, the same age-calss PFTs of forest);
    ! 2. we further allocate the outgoing fraction of each age-class PFT to 
    ! the different receiving youngest age-class PFTs of the same basic-vegetation
    ! type, for example, the youngest age-calss PFTs of cropland.
    
    ! glcc_pft_tmp used only as a temporary variable to store the value
    glcc_pft_tmp(ipts,indold_tree) = veget_max(ipts,indold_tree)/SUM(veget_max(ipts,indold_tree))*frac_used
    glcc_pft(ipts,indold_tree) = glcc_pft(ipts,indold_tree) + glcc_pft_tmp(ipts,indold_tree)
    !we have to remove the outgoing fraction from veget_max in order to use this information for next loop
    veget_max(ipts,indold_tree) = veget_max(ipts,indold_tree) - glcc_pft_tmp(ipts,indold_tree)

    ! when receiving basic-vegetation type has a single age group, it will be considered as
    ! both old and young age group (thus recevie the fraction donation), otherwise the youngest
    ! age group is always the final element of indagec_crop.
    IF (nagec_receive == 1) THEN
      iyoung = 1
    ELSE
      iyoung = nagec_receive - 1
    ENDIF

    ! [20160728] Here we have two options: 
    ! 1. allocate the newly created young age class according to existing fractions
    ! the MTCs.
    ! 2. Use the fractions of MTCs from the current-day PFT map to guid the
    ! allocation of newly created young-age-class MTCs.

    totalveg = 0.   ! [20160130 note here totalveg is the total fraction 
                    !  of all existing MTCs that are going to recieve newly
                    ! convervted fractions.]
    tveg_now = 0.   ! total vegetation fraction in the current-day MTC
                    ! input map.
    DO j=1,num_crop_mulagec
      totalveg = totalveg + veget_mtc(ipts,agec_group(indagec_crop(j,iyoung))) 
      tveg_now = tveg_now + newvegfrac(ipts,agec_group(indagec_crop(j,iyoung)))
    ENDDO

    IF (guide_newfrac) THEN
      IF (tveg_now>min_stomate) THEN
        DO j=1,num_crop_mulagec
          ipft = indagec_crop(j,iyoung)
          glcc_pftmtc(ipts,indold_tree,agec_group(ipft)) = &
                 glcc_pftmtc(ipts,indold_tree,agec_group(ipft)) + &
                 glcc_pft_tmp(ipts,indold_tree) * newvegfrac(ipts,agec_group(ipft))/tveg_now
        ENDDO
      ELSE
        DO j=1,num_crop_mulagec
          ipft = indagec_crop(j,iyoung)
          glcc_pftmtc(ipts,indold_tree,agec_group(ipft)) = &
                glcc_pftmtc(ipts,indold_tree,agec_group(ipft)) + glcc_pft_tmp(ipts,indold_tree)/num_crop_mulagec
        ENDDO
      ENDIF

    ELSE
      IF (totalveg>min_stomate) THEN
        DO j=1,num_crop_mulagec
          ipft = indagec_crop(j,iyoung)
          glcc_pftmtc(ipts,indold_tree,agec_group(ipft)) = &
                 glcc_pftmtc(ipts,indold_tree,agec_group(ipft)) + & 
                 glcc_pft_tmp(ipts,indold_tree) * veget_mtc(ipts,agec_group(ipft))/totalveg
        ENDDO
      ELSE
        DO j=1,num_crop_mulagec
          ipft = indagec_crop(j,iyoung)
          glcc_pftmtc(ipts,indold_tree,agec_group(ipft)) = &
                glcc_pftmtc(ipts,indold_tree,agec_group(ipft)) + glcc_pft_tmp(ipts,indold_tree)/num_crop_mulagec
        ENDDO
      ENDIF

    ENDIF

  END SUBROUTINE cross_give_receive

! ================================================================================================================================
!! SUBROUTINE   : clear_forest
!!
!>\BRIEF        : Handle forest harvest before its legacy is transferred to 
!                 newly initialized youngest-age-class PFT.
!!
!>\DESCRIPTION  
!_ ================================================================================================================================
  !!++TEMP++ biomass,veget_frac are not used because the remaining biomass to be 
  !! harvested is calculated within the deforestation fire module.
  SUBROUTINE clear_forest (npts,ipts,ivm,biomass,frac,    &
                instant_loss, &
                litter, deforest_biomass_remain,&
                fuel_1hr,fuel_10hr,&
                fuel_100hr,fuel_1000hr,&
                lignin_struc,&
                bm_to_litter_pro,convflux,prod10,prod100,&
                litter_pro, fuel_1hr_pro, fuel_10hr_pro, fuel_100hr_pro, &
                fuel_1000hr_pro, lignin_content_pro)


    IMPLICIT NONE

    !! 0.1 Input variables
    INTEGER, INTENT(in)                                       :: npts
    INTEGER, INTENT(in)                                       :: ipts
    INTEGER, INTENT(in)                                       :: ivm
    REAL(r_std), INTENT(in)                                   :: instant_loss
    REAL(r_std), INTENT(in)                                   :: frac   !! the fraction of land covered by forest to be deforested
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)               :: biomass      !! biomass @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)               :: fuel_1hr
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)               :: fuel_10hr
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)               :: fuel_100hr
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)               :: fuel_1000hr
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)             :: litter   !! Vegetmax-weighted remaining litter on the ground for 
                                                                                                      !! deforestation region.
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)               :: deforest_biomass_remain  !! Vegetmax-weighted remaining biomass on the ground for 
                                                                                                      !! deforestation region.
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)         :: lignin_struc     !! ratio Lignine/Carbon in structural litter,
                                                                             !! above and below ground

    !! 0.2 Modified variables
    REAL(r_std), DIMENSION(:,:), INTENT(inout)               :: bm_to_litter_pro    !! conversion of biomass to litter 
                                                                              !! @tex ($gC m^{-2} day^{-1}$) @endtex
    REAL(r_std), DIMENSION(:), INTENT(inout)                 :: convflux         !! release during first year following land cover
                                                                                  !! change

    REAL(r_std), DIMENSION(npts,0:10), INTENT(inout)            :: prod10          !! products remaining in the 10 year-turnover
                                                                              !! pool after the annual release for each 
                                                                              !! compartment (10 + 1 : input from year of land
                                                                              !! cover change)
    REAL(r_std), DIMENSION(npts,0:100), INTENT(inout)           :: prod100         !! products remaining in the 100 year-turnover
                                                                              !! pool after the annual release for each 
                                                                              !! compartment (100 + 1 : input from year of land
                                                                              !! cover change)

    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)          :: litter_pro
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: fuel_1hr_pro
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: fuel_10hr_pro
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: fuel_100hr_pro
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: fuel_1000hr_pro
    REAL(r_std), DIMENSION(:),INTENT(inout)               :: lignin_content_pro



    !! 0.4 Local variables
    REAL(r_std)                                              :: above
      
    ! harvest of aboveground sap- and heartwood biomass after taking into
    ! account of deforestation fire
    IF (allow_deforest_fire) THEN
      above = deforest_biomass_remain(ipts,ivm,isapabove,icarbon)+ &
            deforest_biomass_remain(ipts,ivm,iheartabove,icarbon)
      convflux(ipts)  = convflux(ipts) + 0
      prod10(ipts,0)  = prod10(ipts,0) + 0.4*above
      prod100(ipts,0) = prod100(ipts,0) + 0.6*above
    ELSE
      above = (biomass(ipts,ivm,isapabove,icarbon)+ &
          biomass(ipts,ivm,ileaf,icarbon) + &
          biomass(ipts,ivm,ifruit,icarbon) + &  
          biomass(ipts,ivm,iheartabove,icarbon))*frac
      convflux(ipts)  = convflux(ipts) + instant_loss * above
      !prod10(ipts,0)  = prod10(ipts,0) + coeff_lcchange_10(ivm) * above 
      !prod100(ipts,0) = prod100(ipts,0) + coeff_lcchange_100(ivm) * above 
    ENDIF
  
    ! the transfer of dead biomass to litter
    
    bm_to_litter_pro(isapabove,:) = bm_to_litter_pro(isapabove,:) +  &
                      biomass(ipts,ivm,isapabove,:)*frac*(1-instant_loss)
    bm_to_litter_pro(iheartabove,:) = bm_to_litter_pro(iheartabove,:) +  &
                      biomass(ipts,ivm,iheartabove,:)*frac*(1-instant_loss)

    bm_to_litter_pro(isapbelow,:) = bm_to_litter_pro(isapbelow,:) +  &
                      biomass(ipts,ivm,isapbelow,:)*frac
    bm_to_litter_pro(iheartbelow,:) = bm_to_litter_pro(iheartbelow,:) + &
                      biomass(ipts,ivm,iheartbelow,:)*frac
    bm_to_litter_pro(iroot,:) = bm_to_litter_pro(iroot,:) + &
                      biomass(ipts,ivm,iroot,:)*frac
    bm_to_litter_pro(ifruit,:) = bm_to_litter_pro(ifruit,:) + &
                      biomass(ipts,ivm,ifruit,:)*frac*(1-instant_loss)
    bm_to_litter_pro(icarbres,:) = bm_to_litter_pro(icarbres,:) + &
                      biomass(ipts,ivm,icarbres,:)*frac
    bm_to_litter_pro(ileaf,:) = bm_to_litter_pro(ileaf,:) + &
                      biomass(ipts,ivm,ileaf,:)*frac*(1-instant_loss)

    !update litter_pro
    litter_pro(:,:,:) = litter_pro(:,:,:) + litter(ipts,:,ivm,:,:)*frac
    fuel_1hr_pro(:,:) = fuel_1hr_pro(:,:) + fuel_1hr(ipts,ivm,:,:)*frac
    fuel_10hr_pro(:,:) = fuel_10hr_pro(:,:) + fuel_10hr(ipts,ivm,:,:)*frac 
    fuel_100hr_pro(:,:) = fuel_100hr_pro(:,:) + fuel_100hr(ipts,ivm,:,:)*frac
    fuel_1000hr_pro(:,:) = fuel_1000hr_pro(:,:) + fuel_1000hr(ipts,ivm,:,:)*frac
    !don't forget to hanle litter lignin content
    lignin_content_pro(:)= lignin_content_pro(:) + &
      litter(ipts,istructural,ivm,:,icarbon)*frac*lignin_struc(ipts,ivm,:)

  END SUBROUTINE clear_forest

! ================================================================================================================================
!! SUBROUTINE   : harvest_industrial
!!
!>\BRIEF        : Handle forest harvest before its legacy is transferred to 
!                 newly initialized youngest-age-class PFT.
!!
!>\DESCRIPTION  
!_ ================================================================================================================================
  !!++TEMP++ biomass,veget_frac are not used because the remaining biomass to be 
  !! harvested is calculated within the deforestation fire module.
  SUBROUTINE harvest_industrial (npts,ipts,ivm,biomass,frac,    &
                litter, deforest_biomass_remain,&
                fuel_1hr,fuel_10hr,&
                fuel_100hr,fuel_1000hr,&
                lignin_struc,&
                bm_to_litter_pro,convflux,prod10,prod100,&
                litter_pro, fuel_1hr_pro, fuel_10hr_pro, fuel_100hr_pro, &
                fuel_1000hr_pro, lignin_content_pro)


    IMPLICIT NONE

    !! 0.1 Input variables
    INTEGER, INTENT(in)                                       :: npts
    INTEGER, INTENT(in)                                       :: ipts
    INTEGER, INTENT(in)                                       :: ivm
    REAL(r_std), INTENT(in)                                   :: frac   !! the fraction of land covered by forest to be deforested
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)               :: biomass      !! biomass @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)               :: fuel_1hr
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)               :: fuel_10hr
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)               :: fuel_100hr
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)               :: fuel_1000hr
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)             :: litter   !! Vegetmax-weighted remaining litter on the ground for 
                                                                                                      !! deforestation region.
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)               :: deforest_biomass_remain  !! Vegetmax-weighted remaining biomass on the ground for 
                                                                                                      !! deforestation region.
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)         :: lignin_struc     !! ratio Lignine/Carbon in structural litter,
                                                                             !! above and below ground

    !! 0.2 Modified variables
    REAL(r_std), DIMENSION(:,:), INTENT(inout)               :: bm_to_litter_pro    !! conversion of biomass to litter 
                                                                              !! @tex ($gC m^{-2} day^{-1}$) @endtex
    REAL(r_std), DIMENSION(:), INTENT(inout)                 :: convflux         !! release during first year following land cover
                                                                                  !! change

    REAL(r_std), DIMENSION(npts,0:10), INTENT(inout)            :: prod10          !! products remaining in the 10 year-turnover
                                                                              !! pool after the annual release for each 
                                                                              !! compartment (10 + 1 : input from year of land
                                                                              !! cover change)
    REAL(r_std), DIMENSION(npts,0:100), INTENT(inout)           :: prod100         !! products remaining in the 100 year-turnover
                                                                              !! pool after the annual release for each 
                                                                              !! compartment (100 + 1 : input from year of land
                                                                              !! cover change)

    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)          :: litter_pro
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: fuel_1hr_pro
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: fuel_10hr_pro
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: fuel_100hr_pro
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: fuel_1000hr_pro
    REAL(r_std), DIMENSION(:),INTENT(inout)               :: lignin_content_pro



    !! 0.4 Local variables
    REAL(r_std)                                              :: above
    REAL(r_std)                                              :: slash_frac
    REAL(r_std)                                              :: harvest_efficiency

    ! harvest_efficiency indicates the percentage of the gross area that 
    ! are effetively harvested to reach the target. 
    harvest_efficiency = 0.667
    ! harvest of aboveground sap- and heartwood biomass after taking into
    ! account of deforestation fire
    IF (allow_deforest_fire) THEN
      above = deforest_biomass_remain(ipts,ivm,isapabove,icarbon)+ &
            deforest_biomass_remain(ipts,ivm,iheartabove,icarbon)
      convflux(ipts)  = convflux(ipts) + 0
      prod10(ipts,0)  = prod10(ipts,0) + 0.4*above
      prod100(ipts,0) = prod100(ipts,0) + 0.6*above
    ELSE
      above = (biomass(ipts,ivm,isapabove,icarbon)+ &
          biomass(ipts,ivm,iheartabove,icarbon))*frac
      convflux(ipts)  = convflux(ipts) + harvest_efficiency * coeff_indwood_1(ivm) * above
      prod10(ipts,0)  = prod10(ipts,0) + harvest_efficiency * coeff_indwood_10(ivm) * above 
      prod100(ipts,0) = prod100(ipts,0) + harvest_efficiency * coeff_indwood_100(ivm) * above 
    ENDIF
  
    slash_frac = 1 - (coeff_indwood_10(ivm) + coeff_indwood_100(ivm) &
                   + coeff_indwood_1(ivm)) * harvest_efficiency
    ! the transfer of dead biomass to litter
    bm_to_litter_pro(isapabove,:) = bm_to_litter_pro(isapabove,:) +  &
                      biomass(ipts,ivm,isapabove,:)*frac*slash_frac
    bm_to_litter_pro(iheartabove,:) = bm_to_litter_pro(iheartabove,:) +  &
                      biomass(ipts,ivm,iheartabove,:)*frac*slash_frac
    bm_to_litter_pro(isapbelow,:) = bm_to_litter_pro(isapbelow,:) +  &
                      biomass(ipts,ivm,isapbelow,:)*frac
    bm_to_litter_pro(iheartbelow,:) = bm_to_litter_pro(iheartbelow,:) + &
                      biomass(ipts,ivm,iheartbelow,:)*frac
    bm_to_litter_pro(iroot,:) = bm_to_litter_pro(iroot,:) + &
                      biomass(ipts,ivm,iroot,:)*frac
    bm_to_litter_pro(ifruit,:) = bm_to_litter_pro(ifruit,:) + &
                      biomass(ipts,ivm,ifruit,:)*frac
    bm_to_litter_pro(icarbres,:) = bm_to_litter_pro(icarbres,:) + &
                      biomass(ipts,ivm,icarbres,:)*frac
    bm_to_litter_pro(ileaf,:) = bm_to_litter_pro(ileaf,:) + &
                      biomass(ipts,ivm,ileaf,:)*frac

    !update litter_pro
    litter_pro(:,:,:) = litter_pro(:,:,:) + litter(ipts,:,ivm,:,:)*frac
    fuel_1hr_pro(:,:) = fuel_1hr_pro(:,:) + fuel_1hr(ipts,ivm,:,:)*frac
    fuel_10hr_pro(:,:) = fuel_10hr_pro(:,:) + fuel_10hr(ipts,ivm,:,:)*frac 
    fuel_100hr_pro(:,:) = fuel_100hr_pro(:,:) + fuel_100hr(ipts,ivm,:,:)*frac
    fuel_1000hr_pro(:,:) = fuel_1000hr_pro(:,:) + fuel_1000hr(ipts,ivm,:,:)*frac
    !don't forget to hanle litter lignin content
    lignin_content_pro(:)= lignin_content_pro(:) + &
      litter(ipts,istructural,ivm,:,icarbon)*frac*lignin_struc(ipts,ivm,:)

  END SUBROUTINE harvest_industrial

! ================================================================================================================================
!! SUBROUTINE   : harvest_fuelwood
!!
!>\BRIEF        : Handle forest harvest before its legacy is transferred to 
!                 newly initialized youngest-age-class PFT.
!!
!>\DESCRIPTION  
!_ ================================================================================================================================
  !!++TEMP++ biomass,veget_frac are not used because the remaining biomass to be 
  !! harvested is calculated within the deforestation fire module.
  SUBROUTINE harvest_fuelwood (npts,ipts,ivm,biomass,frac,    &
                litter, deforest_biomass_remain,&
                fuel_1hr,fuel_10hr,&
                fuel_100hr,fuel_1000hr,&
                lignin_struc,&
                bm_to_litter_pro,convflux,prod10,prod100,&
                litter_pro, fuel_1hr_pro, fuel_10hr_pro, fuel_100hr_pro, &
                fuel_1000hr_pro, lignin_content_pro)


    IMPLICIT NONE

    !! 0.1 Input variables
    INTEGER, INTENT(in)                                       :: npts
    INTEGER, INTENT(in)                                       :: ipts
    INTEGER, INTENT(in)                                       :: ivm
    REAL(r_std), INTENT(in)                                   :: frac   !! the fraction of land covered by forest to be deforested
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)               :: biomass      !! biomass @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)               :: fuel_1hr
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)               :: fuel_10hr
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)               :: fuel_100hr
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)               :: fuel_1000hr
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)             :: litter   !! Vegetmax-weighted remaining litter on the ground for 
                                                                                                      !! deforestation region.
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)               :: deforest_biomass_remain  !! Vegetmax-weighted remaining biomass on the ground for 
                                                                                                      !! deforestation region.
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)         :: lignin_struc     !! ratio Lignine/Carbon in structural litter,
                                                                             !! above and below ground

    !! 0.2 Modified variables
    REAL(r_std), DIMENSION(:,:), INTENT(inout)               :: bm_to_litter_pro    !! conversion of biomass to litter 
                                                                              !! @tex ($gC m^{-2} day^{-1}$) @endtex
    REAL(r_std), DIMENSION(:), INTENT(inout)                 :: convflux         !! release during first year following land cover
                                                                                  !! change

    REAL(r_std), DIMENSION(npts,0:10), INTENT(inout)            :: prod10          !! products remaining in the 10 year-turnover
                                                                              !! pool after the annual release for each 
                                                                              !! compartment (10 + 1 : input from year of land
                                                                              !! cover change)
    REAL(r_std), DIMENSION(npts,0:100), INTENT(inout)           :: prod100         !! products remaining in the 100 year-turnover
                                                                              !! pool after the annual release for each 
                                                                              !! compartment (100 + 1 : input from year of land
                                                                              !! cover change)

    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)          :: litter_pro
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: fuel_1hr_pro
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: fuel_10hr_pro
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: fuel_100hr_pro
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: fuel_1000hr_pro
    REAL(r_std), DIMENSION(:),INTENT(inout)               :: lignin_content_pro



    !! 0.4 Local variables
    REAL(r_std)                                              :: above
      
    ! harvest of aboveground sap- and heartwood biomass after taking into
    ! account of deforestation fire
    IF (allow_deforest_fire) THEN
      above = deforest_biomass_remain(ipts,ivm,isapabove,icarbon)+ &
            deforest_biomass_remain(ipts,ivm,iheartabove,icarbon)
      convflux(ipts)  = convflux(ipts) + 0
      prod10(ipts,0)  = prod10(ipts,0) + 0.4*above
      prod100(ipts,0) = prod100(ipts,0) + 0.6*above
    ELSE
      above = (biomass(ipts,ivm,isapabove,icarbon)+ &
          biomass(ipts,ivm,iheartabove,icarbon))*frac
      ! we assume no wood goes to 10-year and 100-year pool in fuelwood collection.
      convflux(ipts)  = convflux(ipts) + above
    ENDIF
  
    ! the transfer of dead biomass to litter
    bm_to_litter_pro(isapbelow,:) = bm_to_litter_pro(isapbelow,:) +  &
                      biomass(ipts,ivm,isapbelow,:)*frac
    bm_to_litter_pro(iheartbelow,:) = bm_to_litter_pro(iheartbelow,:) + &
                      biomass(ipts,ivm,iheartbelow,:)*frac
    bm_to_litter_pro(iroot,:) = bm_to_litter_pro(iroot,:) + &
                      biomass(ipts,ivm,iroot,:)*frac
    bm_to_litter_pro(ifruit,:) = bm_to_litter_pro(ifruit,:) + &
                      biomass(ipts,ivm,ifruit,:)*frac
    bm_to_litter_pro(icarbres,:) = bm_to_litter_pro(icarbres,:) + &
                      biomass(ipts,ivm,icarbres,:)*frac
    bm_to_litter_pro(ileaf,:) = bm_to_litter_pro(ileaf,:) + &
                      biomass(ipts,ivm,ileaf,:)*frac

    !update litter_pro
    litter_pro(:,:,:) = litter_pro(:,:,:) + litter(ipts,:,ivm,:,:)*frac
    fuel_1hr_pro(:,:) = fuel_1hr_pro(:,:) + fuel_1hr(ipts,ivm,:,:)*frac
    fuel_10hr_pro(:,:) = fuel_10hr_pro(:,:) + fuel_10hr(ipts,ivm,:,:)*frac 
    fuel_100hr_pro(:,:) = fuel_100hr_pro(:,:) + fuel_100hr(ipts,ivm,:,:)*frac
    fuel_1000hr_pro(:,:) = fuel_1000hr_pro(:,:) + fuel_1000hr(ipts,ivm,:,:)*frac
    !don't forget to hanle litter lignin content
    lignin_content_pro(:)= lignin_content_pro(:) + &
      litter(ipts,istructural,ivm,:,icarbon)*frac*lignin_struc(ipts,ivm,:)

  END SUBROUTINE harvest_fuelwood
  
! ================================================================================================================================
!! SUBROUTINE   : harvest_herb
!!
!>\BRIEF        : Handle herbaceous PFT clearing before its legacy is transferred to 
!                 newly initialized youngest-age-class PFT.
!!
!>\DESCRIPTION  
!_ ================================================================================================================================
  SUBROUTINE harvest_herb (ipts,ivm,biomass,veget_frac,bm_to_litter_pro)

    IMPLICIT NONE

    !! 0.1 Input variables
    INTEGER, INTENT(in)                                       :: ipts
    INTEGER, INTENT(in)                                       :: ivm
    REAL(r_std), INTENT(in)                                   :: veget_frac   !! the fraction of land covered by herbaceous PFT to be cleared
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)               :: biomass      !! biomass @tex ($gC m^{-2}$) @endtex

    !! 0.2 Modified variables
    REAL(r_std), DIMENSION(:,:), INTENT(inout)                :: bm_to_litter_pro    


    INTEGER :: ipart, num

    ! the transfer of dead biomass to litter
    DO ipart = 1,nparts
      bm_to_litter_pro(ipart,:) = bm_to_litter_pro(ipart,:) + biomass(ipts,ivm,ipart,:)*veget_frac
    ENDDO

  END SUBROUTINE harvest_herb


! ================================================================================================================================
!! SUBROUTINE   : initialize_proxy_pft
!!
!>\BRIEF        Initialize a proxy new youngest age class PFT.
!!
!>\DESCRIPTION  Initialize a proxy new youngest age class PFT that will be 
!!              merged with existing yongest age class, or fill the empty
!!              niche of the youngest age class PFT.
!_ ================================================================================================================================
  SUBROUTINE initialize_proxy_pft(ipts,ipft_young_agec,veget_max_pro,       &
                 biomass_pro, co2_to_bm_pro, ind_pro, age_pro,              & 
                 senescence_pro, PFTpresent_pro,                            &
                 lm_lastyearmax_pro, everywhere_pro, npp_longterm_pro,      &
                 leaf_frac_pro,leaf_age_pro)

    IMPLICIT NONE

    !! 0.1 Input variables
    INTEGER, INTENT(in)                                  :: ipts              !! 
    INTEGER, INTENT(in)                                  :: ipft_young_agec   !! index of the concerned youngest-age-class PFT
    REAL(r_std), INTENT(in)                              :: veget_max_pro     !! fraction of grid cell land area that's to be occupied

    !! 0.2 Modified variables
    REAL(r_std), INTENT(inout)                           :: co2_to_bm_pro

    !! 0.3 Output variables
    REAL(r_std), DIMENSION(:,:), INTENT(out)             :: biomass_pro     !! biomass @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:), INTENT(out)               :: leaf_frac_pro   !! fraction of leaves in leaf age class 
    REAL(r_std), DIMENSION(:), INTENT(out)               :: leaf_age_pro    !! fraction of leaves in leaf age class 
    REAL(r_std), INTENT(out)     :: age_pro, ind_pro, lm_lastyearmax_pro
    REAL(r_std), INTENT(out)                             :: npp_longterm_pro 
    REAL(r_std), INTENT(out)                             :: everywhere_pro  !! is the PFT everywhere in the grid box or very 
    LOGICAL, INTENT(out)                                 :: senescence_pro  !! plant senescent (only for deciduous trees) Set
                                                                            !! to .FALSE. if PFT is introduced or killed
    LOGICAL, INTENT(out)                                 :: PFTpresent_pro  !! Is pft there (unitless)

    !! 0.4 Local variables
    !REAL(r_std), DIMENSION(npts,nvm)                     :: when_growthinit !! how many days ago was the beginning of the 
    !                                                                        !! growing season (days)

    REAL(r_std), DIMENSION(nparts,nelements)               :: bm_new          !! biomass increase @tex ($gC m^{-2}$) @endtex
    REAL(r_std) :: cn_ind,ind
    INTEGER  :: i,j,k,l

    ! -Note-
    ! This part of codes are copied from the original lcchange_main subroutine
    ! that initialize a new PFT.

    i=ipts
    j=ipft_young_agec

    !! Initialization of some variables
    leaf_frac_pro(:) = zero 
    leaf_age_pro(:) = zero 
    
    !! Initial setting of new establishment
    IF (is_tree(j)) THEN
       ! cn_sapl(j)=0.5; stomate_data.f90
       cn_ind = cn_sapl(j) 
    ELSE
       cn_ind = un
    ENDIF
    ind = veget_max_pro / cn_ind
    ind_pro = ind*veget_max_pro
    PFTpresent_pro = .TRUE.
    senescence_pro = .FALSE.
    everywhere_pro = 1.*veget_max_pro
    age_pro = zero

    ! large_value = 1.E33_r_std
    ! when_growthinit(i,j) = large_value 
    leaf_frac_pro(1) = 1.0 * veget_max_pro
    leaf_age_pro(1) = 1.0 * veget_max_pro   !This was not included in original lcchange_main subroutine
    npp_longterm_pro = npp_longterm_init * veget_max_pro
    lm_lastyearmax_pro = bm_sapl(j,ileaf,icarbon) * ind * veget_max_pro
    
    !!  Update of biomass in each each carbon stock component (leaf, sapabove, sapbelow,
    !>  heartabove, heartbelow, root, fruit, and carbres)\n
    DO k = 1, nparts ! loop over # carbon stock components, nparts = 8; stomate_constant.f90 
      DO l = 1,nelements ! loop over # elements
        biomass_pro(k,l) = ind * bm_sapl(j,k,l)
      END DO ! loop over # elements
      co2_to_bm_pro = co2_to_bm_pro + ind * bm_sapl(j,k,icarbon)
    ENDDO ! loop over # carbon stock components
    
  END SUBROUTINE initialize_proxy_pft

! ================================================================================================================================
!! SUBROUTINE   sap_take
!!
!>\BRIEF       : Take the sapling biomass of the new PFTs from the existing biomass, otherwise
!                take from co2_to_bm
!!
!>\DESCRIPTION  
!_ ================================================================================================================================
  SUBROUTINE sap_take (ipts,ivma,veget_max,biomass_pro,biomass,co2_to_bm_pro)

    INTEGER, INTENT(in)                                  :: ipts               !! 
    INTEGER, INTENT(in)                                  :: ivma
    REAL(r_std), DIMENSION(:,:), INTENT(in)              :: veget_max          !! "maximal" coverage fraction of a PFT (LAI ->
    REAL(r_std), DIMENSION(:,:), INTENT(in)              :: biomass_pro        !! biomass @tex ($gC m^{-2}$) @endtex

    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)       :: biomass            !! biomass @tex ($gC m^{-2}$) @endtex
    REAL(r_std), INTENT(inout)                           :: co2_to_bm_pro

    
    REAL(r_std), DIMENSION(nparts,nelements)             :: biomass_total      !! biomass @tex ($gC m^{-2}$) @endtex
    REAL(r_std)                             :: bm_org,bmpro_share
    INTEGER                                 :: i,ivm,ipart
    
    biomass_total(:,:) = zero
    bm_org = zero
    bmpro_share = zero

    DO i = 1,nagec_pft(ivma)
      ivm = start_index(ivma)+i-1
      IF (veget_max(ipts,ivm) .GT. min_stomate) THEN
        biomass_total = biomass_total + biomass(ipts,ivm,:,:)*veget_max(ipts,ivm)
      ENDIF
    ENDDO
  
    DO ipart = 1, nparts
      IF (biomass_total(ipart,icarbon) .GT. biomass_pro(ipart,icarbon)) THEN
        co2_to_bm_pro = co2_to_bm_pro - biomass_pro(ipart,icarbon)
        !treat each PFT of the MTC
        DO i = 1,nagec_pft(ivma)
          ivm = start_index(ivma)+i-1
          IF (veget_max(ipts,ivm) .GT. min_stomate) THEN
            bm_org = biomass(ipts,ivm,ipart,icarbon) * veget_max(ipts,ivm)
            bmpro_share = bm_org/biomass_total(ipart,icarbon) * biomass_pro(ipart,icarbon)
            biomass(ipts,ivm,ipart,icarbon) = (bm_org - bmpro_share)/veget_max(ipts,ivm)
          ENDIF
        ENDDO
      ENDIF
    ENDDO
    
  END SUBROUTINE sap_take

! ================================================================================================================================
!! SUBROUTINE   collect_legacy_pft
!!
!>\BRIEF       : Collect the legacy variables that are going to be included
!                in the newly initialized PFT.
!!
!>\DESCRIPTION  
!_ ================================================================================================================================
  SUBROUTINE collect_legacy_pft(npts, ipts, ivma, glcc_pftmtc,    &
                biomass, bm_to_litter, carbon, litter,            &
                deepC_a, deepC_s, deepC_p,                        &
                fuel_1hr, fuel_10hr, fuel_100hr, fuel_1000hr,     &
                lignin_struc, co2_to_bm, gpp_daily, npp_daily,    &
                resp_maint, resp_growth, resp_hetero, co2_fire,   &
                def_fuel_1hr_remain, def_fuel_10hr_remain,        &
                def_fuel_100hr_remain, def_fuel_1000hr_remain,    &
                deforest_litter_remain, deforest_biomass_remain,  &
                veget_max_pro, carbon_pro, lignin_struc_pro, litter_pro, &
                deepC_a_pro, deepC_s_pro, deepC_p_pro,            &
                fuel_1hr_pro, fuel_10hr_pro, fuel_100hr_pro, fuel_1000hr_pro, &
                bm_to_litter_pro, co2_to_bm_pro, gpp_daily_pro,   &
                npp_daily_pro, resp_maint_pro, resp_growth_pro,   &
                resp_hetero_pro, co2_fire_pro,                    &
                convflux,prod10,prod100)

    IMPLICIT NONE

    !! 0.1 Input variables
    INTEGER, INTENT(in)                                 :: npts               !! Domain size - number of pixels (unitless)
    INTEGER, INTENT(in)                                 :: ipts               !! Domain size - number of pixels (unitless)
    INTEGER, INTENT(in)                                 :: ivma               !! Index for metaclass
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)           :: glcc_pftmtc        !! a temporary variable to hold the fractions each PFT is going to lose
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)         :: biomass            !! biomass @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)         :: bm_to_litter       !! Transfer of biomass to litter 
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)           :: carbon             !! carbon pool: active, slow, or passive 
                                                                              !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)           :: deepC_a            !! Permafrost soil carbon (g/m**3) active
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)           :: deepC_s            !! Permafrost soil carbon (g/m**3) slow
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)           :: deepC_p            !! Permafrost soil carbon (g/m**3) passive
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)       :: litter             !! metabolic and structural litter, above and 
                                                                              !! below ground @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)         :: fuel_1hr
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)         :: fuel_10hr
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)         :: fuel_100hr
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)         :: fuel_1000hr
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)           :: lignin_struc       !! ratio Lignine/Carbon in structural litter,
                                                                              !! above and below ground
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: co2_to_bm          !! biomass uptaken 
                                                                              !! @tex ($gC m^{-2} day^{-1}$) @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: gpp_daily          !! Daily gross primary productivity  
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: npp_daily          !! Net primary productivity 
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: resp_maint         !! Maintenance respiration  
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: resp_growth        !! Growth respiration  
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: resp_hetero        !! Heterotrophic respiration  
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: co2_fire           !! Heterotrophic respiration  
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)               :: def_fuel_1hr_remain
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)               :: def_fuel_10hr_remain
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)               :: def_fuel_100hr_remain
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)               :: def_fuel_1000hr_remain
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)             :: deforest_litter_remain   !! Vegetmax-weighted remaining litter on the ground for 
                                                                                                      !! deforestation region.
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)               :: deforest_biomass_remain  !! Vegetmax-weighted remaining biomass on the ground for 
                                                                                                      !! deforestation region.

    !! 0.2 Output variables
    REAL(r_std), DIMENSION(:), INTENT(inout)              :: carbon_pro
    REAL(r_std), DIMENSION(:), INTENT(inout)              :: deepC_a_pro
    REAL(r_std), DIMENSION(:), INTENT(inout)              :: deepC_s_pro
    REAL(r_std), DIMENSION(:), INTENT(inout)              :: deepC_p_pro
    REAL(r_std), DIMENSION(:), INTENT(inout)              :: lignin_struc_pro   !! ratio Lignine/Carbon in structural litter
                                                                              !! above and below ground
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)          :: litter_pro
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: fuel_1hr_pro
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: fuel_10hr_pro
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: fuel_100hr_pro
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: fuel_1000hr_pro
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: bm_to_litter_pro
    REAL(r_std), INTENT(inout)     :: veget_max_pro, co2_to_bm_pro
    REAL(r_std), INTENT(inout)     :: gpp_daily_pro, npp_daily_pro
    REAL(r_std), INTENT(inout)     :: resp_maint_pro, resp_growth_pro
    REAL(r_std), INTENT(inout)     :: resp_hetero_pro, co2_fire_pro

    !! 0.3 Modified variables
    REAL(r_std), DIMENSION(:,:), INTENT(inout)                 :: convflux      !! release during first year following land cover
                                                                              !! change

    REAL(r_std), DIMENSION(npts,0:10,nwp), INTENT(inout)         :: prod10        !! products remaining in the 10 year-turnover
                                                                              !! pool after the annual release for each 
                                                                              !! compartment (10 + 1 : input from year of land
                                                                              !! cover change)
    REAL(r_std), DIMENSION(npts,0:100,nwp), INTENT(inout)        :: prod100       !! products remaining in the 100 year-turnover
                                                                              !! pool after the annual release for each 
                                                                              !! compartment (100 + 1 : input from year of land
                                                                              !! cover change)

    !! 0.4 Local variables
    REAL(r_std), DIMENSION(nlevs)                  :: lignin_content_pro
    REAL(r_std)                                    :: frac, instant_loss
    INTEGER                                        :: ivm


    ! All *_pro variables collect the legacy pools/fluxes of the ancestor
    ! PFTs for the receiving youngest age class. All *_pro variables 
    ! represent the quantity weighted by the fraction of ancestor contributing
    ! PFTs.
    ! Exceptions:
    ! lignin_struc_pro:: the ratio of lignin content in structural litter.

    lignin_content_pro(:)=zero

    DO ivm = 1,nvm
      frac = glcc_pftmtc(ipts,ivm,ivma)
      IF (frac>zero) THEN
        veget_max_pro = veget_max_pro+frac

        IF (is_tree(ivm)) THEN
          IF (.NOT. is_tree(start_index(ivma))) THEN

            IF (is_grassland_manag(start_index(ivma))) THEN
              instant_loss = 0.67  !forest to pasture
            ELSE IF (.NOT. natural(start_index(ivma))) THEN
              instant_loss = 0.5   !forest to crop
            ELSE
              instant_loss = 0.  !forest to natural grassland
            ENDIF

            CALL clear_forest (npts,ipts,ivm,biomass,frac,    &
                instant_loss, &
                litter, deforest_biomass_remain,&
                fuel_1hr,fuel_10hr,&
                fuel_100hr,fuel_1000hr,&
                lignin_struc,&
                bm_to_litter_pro,convflux(:,iwplcc),prod10(:,:,iwplcc),prod100(:,:,iwplcc),&
                litter_pro, fuel_1hr_pro, fuel_10hr_pro, fuel_100hr_pro, &
                fuel_1000hr_pro, lignin_content_pro)
          ENDIF
        ELSE
          !CALL harvest_herb(ipts,ivm,biomass,frac,   &
          !        bm_to_litter_pro)
          ![2016-04-19] We put the transfer of biomass to litter for herbaceous
          ! PFT directly here, because seprating them in a module harvest_herb
          ! gives some error.
          bm_to_litter_pro(:,:) = bm_to_litter_pro(:,:) + biomass(ipts,ivm,:,:)*frac

          litter_pro(:,:,:) = litter_pro(:,:,:) + litter(ipts,:,ivm,:,:)*frac
          fuel_1hr_pro(:,:) = fuel_1hr_pro(:,:) + fuel_1hr(ipts,ivm,:,:)*frac
          fuel_10hr_pro(:,:) = fuel_10hr_pro(:,:) + fuel_10hr(ipts,ivm,:,:)*frac
          fuel_100hr_pro(:,:) = fuel_100hr_pro(:,:) + fuel_100hr(ipts,ivm,:,:)*frac
          fuel_1000hr_pro(:,:) = fuel_1000hr_pro(:,:) + fuel_1000hr(ipts,ivm,:,:)*frac
          !don't forget to hanle litter lignin content
          lignin_content_pro(:)= lignin_content_pro(:) + &
            litter(ipts,istructural,ivm,:,icarbon)*lignin_struc(ipts,ivm,:)*frac
        ENDIF

        !! scalar variables to be accumulated and inherited
        !! by the destination PFT
        bm_to_litter_pro(:,:) = bm_to_litter_pro(:,:) + &
              bm_to_litter(ipts,ivm,:,:)*frac
        carbon_pro(:) = carbon_pro(:)+carbon(ipts,:,ivm)*frac
        deepC_a_pro(:) = deepC_a_pro(:)+deepC_a(ipts,:,ivm)*frac
        deepC_s_pro(:) = deepC_s_pro(:)+deepC_s(ipts,:,ivm)*frac
        deepC_p_pro(:) = deepC_p_pro(:)+deepC_p(ipts,:,ivm)*frac
        co2_to_bm_pro = co2_to_bm_pro + co2_to_bm(ipts,ivm)*frac

        gpp_daily_pro = gpp_daily_pro + gpp_daily(ipts,ivm)*frac
        npp_daily_pro = npp_daily_pro + npp_daily(ipts,ivm)*frac
        resp_maint_pro = resp_maint_pro + resp_maint(ipts,ivm)*frac
        resp_growth_pro = resp_growth_pro + resp_growth(ipts,ivm)*frac
        resp_hetero_pro = resp_hetero_pro + resp_hetero(ipts,ivm)*frac
        co2_fire_pro = co2_fire_pro + co2_fire(ipts,ivm)*frac
      ENDIF
    ENDDO

    WHERE (litter_pro(istructural,:,icarbon) .GT. min_stomate)
      lignin_struc_pro(:) = lignin_content_pro(:)/litter_pro(istructural,:,icarbon)
    ELSEWHERE
      lignin_struc_pro(:) = zero
    ENDWHERE

  END SUBROUTINE collect_legacy_pft

! ================================================================================================================================
!! SUBROUTINE   collect_legacy_pft_forestry
!!
!>\BRIEF       : Collect the legacy variables that are going to be included
!                in the newly initialized PFT.
!!
!>\DESCRIPTION  
!_ ================================================================================================================================
  SUBROUTINE collect_legacy_pft_forestry(npts, ipts, ivma, glcc_pftmtc,    &
                fuelfrac,                                         &
                biomass, bm_to_litter, carbon, litter,            &
                deepC_a, deepC_s, deepC_p,                        &
                fuel_1hr, fuel_10hr, fuel_100hr, fuel_1000hr,     &
                lignin_struc, co2_to_bm, gpp_daily, npp_daily,    &
                resp_maint, resp_growth, resp_hetero, co2_fire,   &
                def_fuel_1hr_remain, def_fuel_10hr_remain,        &
                def_fuel_100hr_remain, def_fuel_1000hr_remain,    &
                deforest_litter_remain, deforest_biomass_remain,  &
                veget_max_pro, carbon_pro, lignin_struc_pro, litter_pro, &
                deepC_a_pro, deepC_s_pro, deepC_p_pro,            &
                fuel_1hr_pro, fuel_10hr_pro, fuel_100hr_pro, fuel_1000hr_pro, &
                bm_to_litter_pro, co2_to_bm_pro, gpp_daily_pro,   &
                npp_daily_pro, resp_maint_pro, resp_growth_pro,   &
                resp_hetero_pro, co2_fire_pro,                    &
                convflux,prod10,prod100)

    IMPLICIT NONE

    !! 0.1 Input variables
    INTEGER, INTENT(in)                                 :: npts               !! Domain size - number of pixels (unitless)
    INTEGER, INTENT(in)                                 :: ipts               !! Domain size - number of pixels (unitless)
    INTEGER, INTENT(in)                                 :: ivma               !! Index for metaclass
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)           :: glcc_pftmtc        !! a temporary variable to hold the fractions each PFT is going to lose
    REAL(r_std), DIMENSION(:), INTENT(in)               :: fuelfrac           !! 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)         :: biomass            !! biomass @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)         :: bm_to_litter       !! Transfer of biomass to litter 
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)           :: carbon             !! carbon pool: active, slow, or passive 
                                                                              !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)           :: deepC_a            !! Permafrost soil carbon (g/m**3) active
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)           :: deepC_s            !! Permafrost soil carbon (g/m**3) slow
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)           :: deepC_p            !! Permafrost soil carbon (g/m**3) passive
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)       :: litter             !! metabolic and structural litter, above and 
                                                                              !! below ground @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)         :: fuel_1hr
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)         :: fuel_10hr
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)         :: fuel_100hr
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)         :: fuel_1000hr
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)           :: lignin_struc       !! ratio Lignine/Carbon in structural litter,
                                                                              !! above and below ground
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: co2_to_bm          !! biomass uptaken 
                                                                              !! @tex ($gC m^{-2} day^{-1}$) @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: gpp_daily          !! Daily gross primary productivity  
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: npp_daily          !! Net primary productivity 
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: resp_maint         !! Maintenance respiration  
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: resp_growth        !! Growth respiration  
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: resp_hetero        !! Heterotrophic respiration  
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: co2_fire           !! Heterotrophic respiration  
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)               :: def_fuel_1hr_remain
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)               :: def_fuel_10hr_remain
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)               :: def_fuel_100hr_remain
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)               :: def_fuel_1000hr_remain
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)             :: deforest_litter_remain   !! Vegetmax-weighted remaining litter on the ground for 
                                                                                                      !! deforestation region.
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)               :: deforest_biomass_remain  !! Vegetmax-weighted remaining biomass on the ground for 
                                                                                                      !! deforestation region.

    !! 0.2 Output variables
    REAL(r_std), DIMENSION(:), INTENT(inout)              :: carbon_pro
    REAL(r_std), DIMENSION(:), INTENT(inout)              :: deepC_a_pro
    REAL(r_std), DIMENSION(:), INTENT(inout)              :: deepC_s_pro
    REAL(r_std), DIMENSION(:), INTENT(inout)              :: deepC_p_pro
    REAL(r_std), DIMENSION(:), INTENT(inout)              :: lignin_struc_pro   !! ratio Lignine/Carbon in structural litter
                                                                              !! above and below ground
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)          :: litter_pro
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: fuel_1hr_pro
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: fuel_10hr_pro
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: fuel_100hr_pro
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: fuel_1000hr_pro
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: bm_to_litter_pro
    REAL(r_std), INTENT(inout)     :: veget_max_pro, co2_to_bm_pro
    REAL(r_std), INTENT(inout)     :: gpp_daily_pro, npp_daily_pro
    REAL(r_std), INTENT(inout)     :: resp_maint_pro, resp_growth_pro
    REAL(r_std), INTENT(inout)     :: resp_hetero_pro, co2_fire_pro

    !! 0.3 Modified variables
    REAL(r_std), DIMENSION(:,:), INTENT(inout)                 :: convflux      !! release during first year following land cover
                                                                              !! change

    REAL(r_std), DIMENSION(npts,0:10,nwp), INTENT(inout)         :: prod10        !! products remaining in the 10 year-turnover
                                                                              !! pool after the annual release for each 
                                                                              !! compartment (10 + 1 : input from year of land
                                                                              !! cover change)
    REAL(r_std), DIMENSION(npts,0:100,nwp), INTENT(inout)        :: prod100       !! products remaining in the 100 year-turnover
                                                                              !! pool after the annual release for each 
                                                                              !! compartment (100 + 1 : input from year of land
                                                                              !! cover change)

    !! 0.4 Local variables
    REAL(r_std), DIMENSION(nlevs)                  :: lignin_content_pro
    REAL(r_std)                                    :: frac,pfuelfrac
    INTEGER                                        :: ivm


    ! All *_pro variables collect the legacy pools/fluxes of the ancestor
    ! PFTs for the receiving youngest age class. All *_pro variables 
    ! represent the quantity weighted by the fraction of ancestor contributing
    ! PFTs.
    ! Exceptions:
    ! lignin_struc_pro:: the ratio of lignin content in structural litter.

    lignin_content_pro(:)=zero
    pfuelfrac = fuelfrac(ipts)

    DO ivm = 1,nvm
      frac = glcc_pftmtc(ipts,ivm,ivma)
      IF (frac>zero) THEN
        veget_max_pro = veget_max_pro+frac

        IF (is_tree(ivm)) THEN
          IF (is_tree(start_index(ivma))) THEN
            CALL harvest_industrial (npts,ipts,ivm,biomass,frac*(1-pfuelfrac),    &
                litter, deforest_biomass_remain,&
                fuel_1hr,fuel_10hr,&
                fuel_100hr,fuel_1000hr,&
                lignin_struc,&
                bm_to_litter_pro,convflux(:,iwphar),prod10(:,:,iwphar),prod100(:,:,iwphar),&
                litter_pro, fuel_1hr_pro, fuel_10hr_pro, fuel_100hr_pro, &
                fuel_1000hr_pro, lignin_content_pro)

            CALL harvest_fuelwood (npts,ipts,ivm,biomass,frac*pfuelfrac,    &
                litter, deforest_biomass_remain,&
                fuel_1hr,fuel_10hr,&
                fuel_100hr,fuel_1000hr,&
                lignin_struc,&
                bm_to_litter_pro,convflux(:,iwphar),prod10(:,:,iwphar),prod100(:,:,iwphar),&
                litter_pro, fuel_1hr_pro, fuel_10hr_pro, fuel_100hr_pro, &
                fuel_1000hr_pro, lignin_content_pro)
          ENDIF
        ELSE
          !CALL harvest_herb(ipts,ivm,biomass,frac,   &
          !        bm_to_litter_pro)
          ![2016-04-19] We put the transfer of biomass to litter for herbaceous
          ! PFT directly here, because seprating them in a module harvest_herb
          ! gives some error.
          bm_to_litter_pro(:,:) = bm_to_litter_pro(:,:) + biomass(ipts,ivm,:,:)*frac

          litter_pro(:,:,:) = litter_pro(:,:,:) + litter(ipts,:,ivm,:,:)*frac
          fuel_1hr_pro(:,:) = fuel_1hr_pro(:,:) + fuel_1hr(ipts,ivm,:,:)*frac
          fuel_10hr_pro(:,:) = fuel_10hr_pro(:,:) + fuel_10hr(ipts,ivm,:,:)*frac
          fuel_100hr_pro(:,:) = fuel_100hr_pro(:,:) + fuel_100hr(ipts,ivm,:,:)*frac
          fuel_1000hr_pro(:,:) = fuel_1000hr_pro(:,:) + fuel_1000hr(ipts,ivm,:,:)*frac
          !don't forget to hanle litter lignin content
          lignin_content_pro(:)= lignin_content_pro(:) + &
            litter(ipts,istructural,ivm,:,icarbon)*lignin_struc(ipts,ivm,:)*frac
        ENDIF

        !! scalar variables to be accumulated and inherited
        !! by the destination PFT
        bm_to_litter_pro(:,:) = bm_to_litter_pro(:,:) + &
              bm_to_litter(ipts,ivm,:,:)*frac
        carbon_pro(:) = carbon_pro(:)+carbon(ipts,:,ivm)*frac
        deepC_a_pro(:) = deepC_a_pro(:)+deepC_a(ipts,:,ivm)*frac
        deepC_s_pro(:) = deepC_s_pro(:)+deepC_s(ipts,:,ivm)*frac
        deepC_p_pro(:) = deepC_p_pro(:)+deepC_p(ipts,:,ivm)*frac
        co2_to_bm_pro = co2_to_bm_pro + co2_to_bm(ipts,ivm)*frac

        gpp_daily_pro = gpp_daily_pro + gpp_daily(ipts,ivm)*frac
        npp_daily_pro = npp_daily_pro + npp_daily(ipts,ivm)*frac
        resp_maint_pro = resp_maint_pro + resp_maint(ipts,ivm)*frac
        resp_growth_pro = resp_growth_pro + resp_growth(ipts,ivm)*frac
        resp_hetero_pro = resp_hetero_pro + resp_hetero(ipts,ivm)*frac
        co2_fire_pro = co2_fire_pro + co2_fire(ipts,ivm)*frac
      ENDIF
    ENDDO

    WHERE (litter_pro(istructural,:,icarbon) .GT. min_stomate)
      lignin_struc_pro(:) = lignin_content_pro(:)/litter_pro(istructural,:,icarbon)
    ELSEWHERE
      lignin_struc_pro(:) = zero
    ENDWHERE

  END SUBROUTINE collect_legacy_pft_forestry


! ================================================================================================================================
!! SUBROUTINE   : add_incoming_proxy_pft
!!
!>\BRIEF        : Merge the newly incoming proxy PFT cohort with the exisiting
!!                cohort.
!! \n
!
!_ ================================================================================================================================
  SUBROUTINE add_incoming_proxy_pft(npts, ipts, ipft, veget_max_pro,  &
       carbon_pro, litter_pro, lignin_struc_pro, bm_to_litter_pro,    &
       deepC_a_pro, deepC_s_pro, deepC_p_pro,                         &
       fuel_1hr_pro, fuel_10hr_pro, fuel_100hr_pro, fuel_1000hr_pro,  &
       biomass_pro, co2_to_bm_pro, npp_longterm_pro, ind_pro,         &
       lm_lastyearmax_pro, age_pro, everywhere_pro,                   &  
       leaf_frac_pro, leaf_age_pro, PFTpresent_pro, senescence_pro,   &
       gpp_daily_pro, npp_daily_pro, resp_maint_pro, resp_growth_pro, &
       resp_hetero_pro, co2_fire_pro,                                 &
       veget_max, carbon, litter, lignin_struc, bm_to_litter,         &
       deepC_a, deepC_s, deepC_p,                                     &
       fuel_1hr, fuel_10hr, fuel_100hr, fuel_1000hr,                  &
       biomass, co2_to_bm, npp_longterm, ind,                         &
       lm_lastyearmax, age, everywhere,                               &
       leaf_frac, leaf_age, PFTpresent, senescence,                   &
       gpp_daily, npp_daily, resp_maint, resp_growth,                 &
       resp_hetero, co2_fire)
    
    IMPLICIT NONE

    !! 0.1 Input variables
    INTEGER, INTENT(in)                                :: npts                !! Domain size - number of pixels (unitless)
    INTEGER, INTENT(in)                                :: ipts                !! Domain size - number of pixels (unitless)
    INTEGER, INTENT(in)                                :: ipft
    REAL(r_std), INTENT(in)                            :: veget_max_pro           !! The land fraction of incoming new PFTs that are 
                                                                              !! the sum of all its ancestor PFTs
    REAL(r_std), DIMENSION(:), INTENT(in)              :: carbon_pro
    REAL(r_std), DIMENSION(:), INTENT(in)              :: deepC_a_pro
    REAL(r_std), DIMENSION(:), INTENT(in)              :: deepC_s_pro
    REAL(r_std), DIMENSION(:), INTENT(in)              :: deepC_p_pro
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)          :: litter_pro
    REAL(r_std), DIMENSION(:,:), INTENT(in)            :: fuel_1hr_pro
    REAL(r_std), DIMENSION(:,:), INTENT(in)            :: fuel_10hr_pro
    REAL(r_std), DIMENSION(:,:), INTENT(in)            :: fuel_100hr_pro
    REAL(r_std), DIMENSION(:,:), INTENT(in)            :: fuel_1000hr_pro
    REAL(r_std), DIMENSION(:,:), INTENT(in)            :: bm_to_litter_pro
    REAL(r_std), DIMENSION(:), INTENT(in)              :: lignin_struc_pro    !! ratio Lignine/Carbon in structural litter
                                                                              !! above and below ground
    REAL(r_std), DIMENSION(:,:), INTENT(in)            :: biomass_pro         !! biomass @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:), INTENT(in)              :: leaf_frac_pro       !! fraction of leaves in leaf age class 
    REAL(r_std), DIMENSION(:), INTENT(in)              :: leaf_age_pro        !! fraction of leaves in leaf age class 
    REAL(r_std), INTENT(in)     :: ind_pro, age_pro, lm_lastyearmax_pro
    REAL(r_std), INTENT(in)     :: npp_longterm_pro, co2_to_bm_pro 
    REAL(r_std), INTENT(in)                            :: everywhere_pro      !! is the PFT everywhere in the grid box or very 
    LOGICAL, INTENT(in)         :: PFTpresent_pro, senescence_pro             !! Is pft there (unitless)

    REAL(r_std), INTENT(in)     :: gpp_daily_pro, npp_daily_pro
    REAL(r_std), INTENT(in)     :: resp_maint_pro, resp_growth_pro
    REAL(r_std), INTENT(in)     :: resp_hetero_pro, co2_fire_pro

    !! 0.2 Output variables

    !! 0.3 Modified variables

    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: veget_max           !! "maximal" coverage fraction of a PFT on the ground
                                                                              !! May sum to
                                                                              !! less than unity if the pixel has
                                                                              !! nobio area. (unitless, 0-1)
   
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: carbon              !! carbon pool: active, slow, or passive 
                                                                              !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: deepC_a             !! Permafrost soil carbon (g/m**3) active
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: deepC_s             !! Permafrost soil carbon (g/m**3) slow
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: deepC_p             !! Permafrost soil carbon (g/m**3) passive
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)   :: litter              !! metabolic and structural litter, above and 
                                                                              !! below ground @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: fuel_1hr
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: fuel_10hr
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: fuel_100hr
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: fuel_1000hr
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: lignin_struc        !! ratio Lignine/Carbon in structural litter,
                                                                              !! above and below ground
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: bm_to_litter        !! Transfer of biomass to litter 
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: biomass             !! Stand level biomass @tex $(gC.m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: co2_to_bm           !! CO2 taken from the atmosphere to get C to create  
                                                                              !! the seedlings @tex (gC.m^{-2}dt^{-1})$ @endtex

    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: npp_longterm        !! "Long term" mean yearly primary productivity 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: ind                 !! Number of individuals at the stand level
                                                                              !! @tex $(m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: age                 !! mean age (years)
    LOGICAL, DIMENSION(:,:), INTENT(inout)             :: PFTpresent          !! Tab indicating which PFTs are present in 
                                                                              !! each pixel
    LOGICAL, DIMENSION(:,:), INTENT(inout)             :: senescence          !! Flag for setting senescence stage (only 
                                                                              !! for deciduous trees)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: lm_lastyearmax      !! last year's maximum leaf mass for each PFT 
                                                                              !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: everywhere          !! is the PFT everywhere in the grid box or 
                                                                              !! very localized (after its introduction) (?)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: leaf_frac           !! fraction of leaves in leaf age class (unitless;0-1)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: leaf_age            !! Leaf age (days)

    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: gpp_daily           !! Daily gross primary productivity  
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: npp_daily           !! Net primary productivity 
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: resp_maint          !! Maintenance respiration  
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: resp_growth         !! Growth respiration  
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: resp_hetero         !! Heterotrophic respiration  
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: co2_fire            !! Heterotrophic respiration  
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 

    ! REAL(r_std), DIMENSION(:,:), INTENT(inout)        :: moiavail_month       !! "Monthly" moisture availability (0 to 1, 
    !                                                                           !! unitless) 
    ! REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: moiavail_week       !! "Weekly" moisture availability 
    !                                                                           !! (0 to 1, unitless)
    ! REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: gpp_week            !! Mean weekly gross primary productivity 
    !                                                                           !! @tex $(gC m^{-2} day^{-1})$ @endtex
    ! REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: ngd_minus5          !! Number of growing days (days), threshold 
    !                                                                           !! -5 deg C (for phenology)   
    ! REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: when_growthinit     !! How many days ago was the beginning of 
    !                                                                           !! the growing season (days)
    ! REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: time_hum_min        !! Time elapsed since strongest moisture 
    !                                                                           !! availability (days) 
    ! REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: gdd_midwinter       !! Growing degree days (K), since midwinter 
    !                                                                           !! (for phenology) - this is written to the
    !                                                                           !!  history files 
    ! REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: gdd_from_growthinit !! growing degree days, since growthinit 
    !                                                                           !! for crops
    ! REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: gdd_m5_dormance     !! Growing degree days (K), threshold -5 deg 
    !                                                                           !! C (for phenology)
    ! REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: ncd_dormance        !! Number of chilling days (days), since 
    !                                                                           !! leaves were lost (for phenology) 

    !! 0.4 Local variables

    INTEGER(i_std)                                     :: iele                !! Indeces(unitless)
    INTEGER(i_std)                                     :: ilit,ilev,icarb     !! Indeces(unitless)
    REAL(r_std), DIMENSION(npts,nlitt,nvm,nlevs,nelements) :: litter_old      !! metabolic and structural litter, above and 
                                                                              !! below ground @tex ($gC m^{-2}$) @endtex
    REAL(r_std) :: veget_old,veget_total
  
    
    ! Back up some variables in case they're needed later
    litter_old(:,:,:,:,:) = litter(:,:,:,:,:)

    !! General idea
    ! The established proxy vegetation has a fraction of 'veget_max_pro'; the
    ! existing iPFT has a fraction of veget_max(ipts,ipft). 
    ! Suppose we want to merge a scalar variable B, the value of B after merging
    ! is (Bi*Vi+Bj*Vj)/(Vi+Vj), where Vi is the original veget_max, Vj is the
    ! incoming veget_max. Note that in case Vi=0, this equation remains solid,
    ! i.e. the veget_max after merging is Vj and B after merging is Bj. In other
    ! words, the proxy vegetation "fills" up the empty niche of iPFT.
    ! Also note that for many scalar variables our input value is Bj*Vj, which
    ! is accumulated from multiple ancestor PFTs.
    veget_old = veget_max(ipts,ipft)
    veget_total = veget_old+veget_max_pro

    !! Different ways of handling merging depending on nature of variables:

    !! 1. Area-based scalar variables, use the equation above
    !  biomass,carbon, litter, bm_to_litter, co2_to_bm, ind,
    !  lm_lastyearmax, npp_longterm, lm_lastyearmax,
    !  lignin_struc (ratio variable depending on area-based variable)
     
    !! 2. Variables are tentatively handled like area-based variables:
    !   leaf_frac, leaf_age,

    !! 3. Variables that are overwritten by the newly initialized PFT:
    !   PFTpresent, senescence

    !! 4. Variables whose operation is uncertain and are not handled currently:
    !  when_growthinit :: how many days ago was the beginning of the growing season (days)
    !  gdd_from_growthinit :: growing degree days, since growthinit 
    !  gdd_midwinter, time_hum_min, gdd_m5_dormance, ncd_dormance, 
    !  moiavail_month, moiavail_week, ngd_minus5

    !! 5. Variables that concern with short-term fluxes that do not apply in
    !  this case:
    !  gpp_daily, npp_daily etc.

    ! Add the coming veget_max_pro into existing veget_max
    veget_max(ipts,ipft) = veget_total

    IF (veget_total .GT. min_stomate) THEN
      ! Merge scalar variables which are defined on area basis
      carbon(ipts,:,ipft) =  (veget_old * carbon(ipts,:,ipft) + &
           carbon_pro(:))/veget_total
      deepC_a(ipts,:,ipft) =  (veget_old * deepC_a(ipts,:,ipft) + &
           deepC_a_pro(:))/veget_total
      deepC_s(ipts,:,ipft) =  (veget_old * deepC_s(ipts,:,ipft) + &
           deepC_s_pro(:))/veget_total
      deepC_p(ipts,:,ipft) =  (veget_old * deepC_p(ipts,:,ipft) + &
           deepC_p_pro(:))/veget_total
      litter(ipts,:,ipft,:,:) = (veget_old * litter(ipts,:,ipft,:,:) + &
           litter_pro(:,:,:))/veget_total
      fuel_1hr(ipts,ipft,:,:) = (veget_old * fuel_1hr(ipts,ipft,:,:) + &
           fuel_1hr_pro(:,:))/veget_total
      fuel_10hr(ipts,ipft,:,:) = (veget_old * fuel_10hr(ipts,ipft,:,:) + &
           fuel_10hr_pro(:,:))/veget_total
      fuel_100hr(ipts,ipft,:,:) = (veget_old * fuel_100hr(ipts,ipft,:,:) + &
           fuel_100hr_pro(:,:))/veget_total
      fuel_1000hr(ipts,ipft,:,:) = (veget_old * fuel_1000hr(ipts,ipft,:,:) + &
           fuel_1000hr_pro(:,:))/veget_total

      WHERE (litter(ipts,istructural,ipft,:,icarbon) .GT. min_stomate) 
        lignin_struc(ipts,ipft,:) = (veget_old*litter_old(ipts,istructural,ipft,:,icarbon)* &
            lignin_struc(ipts,ipft,:) + litter_pro(istructural,:,icarbon)* &
            lignin_struc_pro(:))/(veget_total*litter(ipts,istructural,ipft,:,icarbon))
      ENDWHERE
      bm_to_litter(ipts,ipft,:,:) = (veget_old * bm_to_litter(ipts,ipft,:,:) + & 
           bm_to_litter_pro(:,:))/veget_total

      biomass(ipts,ipft,:,:) = (biomass(ipts,ipft,:,:)*veget_old + &
           biomass_pro(:,:))/veget_total
      co2_to_bm(ipts,ipft) = (veget_old*co2_to_bm(ipts,ipft) + &
           co2_to_bm_pro)/veget_total
      ind(ipts,ipft) = (ind(ipts,ipft)*veget_old + ind_pro)/veget_total
      lm_lastyearmax(ipts,ipft) = (lm_lastyearmax(ipts,ipft)*veget_old + &
           lm_lastyearmax_pro)/veget_total
      npp_longterm(ipts,ipft) = (veget_old * npp_longterm(ipts,ipft) + &
           npp_longterm_pro)/veget_total

      !CHECK: Here follows the original idea in DOFOCO, more strictly,
      ! leas mass should be considered together. The same also applies on
      ! leaf age.
      leaf_frac(ipts,ipft,:) = (leaf_frac(ipts,ipft,:)*veget_old + &
           leaf_frac_pro(:))/veget_total
      leaf_age(ipts,ipft,:) = (leaf_age(ipts,ipft,:)*veget_old + &
           leaf_age_pro(:))/veget_total
      age(ipts,ipft) = (veget_old * age(ipts,ipft) + &
           age_pro)/veget_total

      ! Everywhere deals with the migration of vegetation. Copy the
      ! status of the most migrated vegetation for the whole PFT
      everywhere(ipts,ipft) = MAX(everywhere(ipts,ipft), everywhere_pro)

      ! Overwrite the original variables with that from newly initialized
      ! proxy PFT
      PFTpresent(ipts,ipft) = PFTpresent_pro
      senescence(ipts,ipft) = senescence_pro

      ! This is to close carbon loop when writing history variables.
      gpp_daily(ipts,ipft) = (veget_old * gpp_daily(ipts,ipft) + &
           gpp_daily_pro)/veget_total
      npp_daily(ipts,ipft) = (veget_old * npp_daily(ipts,ipft) + &
           npp_daily_pro)/veget_total
      resp_maint(ipts,ipft) = (veget_old * resp_maint(ipts,ipft) + &
           resp_maint_pro)/veget_total 
      resp_growth(ipts,ipft) = (veget_old * resp_growth(ipts,ipft) + &
           resp_growth_pro)/veget_total 
      resp_hetero(ipts,ipft) = (veget_old * resp_hetero(ipts,ipft) + &
           resp_hetero_pro)/veget_total 
      co2_fire(ipts,ipft) = (veget_old * co2_fire(ipts,ipft) + &
           co2_fire_pro)/veget_total 

      ! Phenology- or time-related variables will be copied from original values if 
      ! there is already youngest-age-class PFT there, otherwise they're left
      ! untouched, because 1. to initiliaze all new PFTs here is wrong and 
      ! phenology is not explicitly considered, so we cannot assign a value
      ! to these variables. 2. We assume they will be correctly filled if 
      ! other variables are in place (e.g., non-zero leaf mass will lead to
      ! onset of growing season). In this case, merging a newly initialized PFT
      ! to an existing one is not the same as merging PFTs when they grow 
      ! old enough to exceed thresholds.
      
      ! gpp_week(ipts,ipft) = (veget_old * gpp_week(ipts,ipft) + &
      !      gpp_week_pro)/veget_total
      ! when_growthinit(ipts,ipft) = (veget_old * when_growthinit(ipts,ipft) + &
      !      when_growthinit_pro)/veget_total
      ! gdd_from_growthinit(ipts,ipft) = (veget_old * gdd_from_growthinit(ipts,ipft) + &
      !      gdd_from_growthinit_pro)/veget_total
      ! gdd_midwinter(ipts,ipft) = (veget_old * gdd_midwinter(ipts,ipft) + &
      !      gdd_midwinter_pro)/veget_total
      ! time_hum_min(ipts,ipft) = (veget_old * time_hum_min(ipts,ipft) + &
      !      time_hum_min_pro)/veget_total
      ! gdd_m5_dormance(ipts,ipft) = (veget_old * gdd_m5_dormance(ipts,ipft) + &
      !      gdd_m5_dormance_pro)/veget_total
      ! ncd_dormance(ipts,ipft) = (veget_old * ncd_dormance(ipts,ipft) + &
      !      ncd_dormance_pro)/veget_total
      ! moiavail_month(ipts,ipft) = (veget_old * moiavail_month(ipts,ipft) + &
      !      moiavail_month_pro)/veget_total
      ! moiavail_week(ipts,ipft) = (veget_old * moiavail_week(ipts,ipft) + &
      !      moiavail_week_pro)/veget_total
      ! ngd_minus5(ipts,ipft) = (veget_old * ngd_minus5(ipts,ipft) + &
      !      ngd_minus5_pro)/veget_total
    ENDIF
    
  
  END SUBROUTINE add_incoming_proxy_pft

! ================================================================================================================================
!! SUBROUTINE   : empty_pft
!!
!>\BRIEF        : Empty a PFT when,
!!                - it is exhausted because of land cover change.
!!                - it moves to the next age class
!! \n
!_ ================================================================================================================================
  SUBROUTINE empty_pft(ipts, ivm, veget_max, biomass, ind,       &
               bm_phytomer,bm_FFB,PHYbm, FFBbm,                  & !! yidi
               carbon, litter, lignin_struc, bm_to_litter,       &
               deepC_a, deepC_s, deepC_p,                        &
               fuel_1hr, fuel_10hr, fuel_100hr, fuel_1000hr,     &
               gpp_daily, npp_daily, gpp_week, npp_longterm,     &
               co2_to_bm, resp_maint, resp_growth, resp_hetero,  &
               lm_lastyearmax, leaf_frac, leaf_age, age,         &
               everywhere, PFTpresent, when_growthinit,          &
               senescence, gdd_from_growthinit, gdd_midwinter,   &
               time_hum_min, gdd_m5_dormance, ncd_dormance,      &
               moiavail_month, moiavail_week, ngd_minus5)
    
    IMPLICIT NONE

    !! 0.1 Input variables
    INTEGER, INTENT(in)                                :: ipts               !! index for grid cell
    INTEGER, INTENT(in)                                :: ivm                !! index for pft

    !! 0.2 Output variables

    !! 0.3 Modified variables

    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: veget_max           !! "maximal" coverage fraction of a PFT on the ground
                                                                              !! May sum to
                                                                              !! less than unity if the pixel has
                                                                              !! nobio area. (unitless, 0-1)
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: biomass             !! Stand level biomass @tex $(gC.m^{-2})$ @endtex
!! yidi
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: FFBbm               !! FFB mass, from sapabove            
                                                                              !! @tex $(gC.m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: PHYbm               !! PHYTOMER mass, from sapabove      
                                                                              !! @tex $(gC.m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: bm_phytomer         !! Each PHYTOMER mass, from sapabove  
                                                                              !! @tex $(gC.m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: bm_FFB              !! Fruit mass for Each PHYTOMER   
                                                                              !! @tex $(gC.m^{-2})$ @endtex
!! yidi
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: ind                 !! Number of individuals at the stand level
                                                                              !! @tex $(m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: carbon              !! carbon pool: active, slow, or passive 
                                                                              !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: deepC_a             !! Permafrost soil carbon (g/m**3) active
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: deepC_s             !! Permafrost soil carbon (g/m**3) slow
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: deepC_p             !! Permafrost soil carbon (g/m**3) passive
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)   :: litter              !! metabolic and structural litter, above and 
                                                                              !! below ground @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: fuel_1hr
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: fuel_10hr
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: fuel_100hr
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: fuel_1000hr
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: lignin_struc        !! ratio Lignine/Carbon in structural litter,
                                                                              !! above and below ground
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: bm_to_litter        !! Transfer of biomass to litter 
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: gpp_daily           !! Daily gross primary productivity  
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: npp_daily           !! Net primary productivity 
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: gpp_week            !! Mean weekly gross primary productivity 
                                                                              !! @tex $(gC m^{-2} day^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: npp_longterm        !! "Long term" mean yearly primary productivity 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: co2_to_bm           !! CO2 taken from the atmosphere to get C to create  
                                                                              !! the seedlings @tex (gC.m^{-2}dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: resp_maint          !! Maintenance respiration  
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: resp_growth         !! Growth respiration  
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: resp_hetero         !! Heterotrophic respiration  
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: lm_lastyearmax      !! last year's maximum leaf mass for each PFT 
                                                                              !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: leaf_frac           !! fraction of leaves in leaf age class (unitless;0-1)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: leaf_age            !! Leaf age (days)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: age                 !! mean age (years)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: everywhere          !! is the PFT everywhere in the grid box or 
                                                                              !! very localized (after its introduction) (?)
    LOGICAL, DIMENSION(:,:), INTENT(inout)             :: PFTpresent          !! Tab indicating which PFTs are present in 
                                                                              !! each pixel
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: when_growthinit     !! How many days ago was the beginning of 
                                                                              !! the growing season (days)
    LOGICAL, DIMENSION(:,:), INTENT(inout)             :: senescence          !! Flag for setting senescence stage (only 
                                                                              !! for deciduous trees)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: gdd_from_growthinit !! growing degree days, since growthinit 
                                                                              !! for crops
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: gdd_midwinter       !! Growing degree days (K), since midwinter 
                                                                              !! (for phenology) - this is written to the
                                                                              !!  history files 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: time_hum_min        !! Time elapsed since strongest moisture 
                                                                              !! availability (days) 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: gdd_m5_dormance     !! Growing degree days (K), threshold -5 deg 
                                                                              !! C (for phenology)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: ncd_dormance        !! Number of chilling days (days), since 
                                                                              !! leaves were lost (for phenology) 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: moiavail_month      !! "Monthly" moisture availability (0 to 1, 
                                                                              !! unitless) 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: moiavail_week       !! "Weekly" moisture availability 
                                                                              !! (0 to 1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: ngd_minus5          !! Number of growing days (days), threshold 
                                                                              !! -5 deg C (for phenology)   

    !! 0.4 Local variables
    INTEGER(i_std)                                     :: iele                !! Indeces(unitless)
    INTEGER(i_std)                                     :: ilit,ilev,icarb     !! Indeces(unitless)

    WRITE(numout,*) 'Entering empty_pft'
    WRITE(numout,*) 'yd: GLUC empty pft before,nvm=',ivm,',biomass(:,ivm,:,icarbon)',biomass(:,ivm,:,icarbon)
    veget_max(ipts,ivm) = zero
    ind(ipts,ivm) = zero
    biomass(ipts,ivm,:,:) = zero
    
    WRITE(numout,*) 'yd: GLUC empty pft after,nvm=',ivm,',biomass(:,ivm,:,icarbon)',biomass(:,ivm,:,icarbon)
!! yidi
    IF (ok_oilpalm) THEN
        WRITE(numout,*) 'yd: GLUC empty pft before, j=',ivm,'bm_phytomer(i,j,:)',bm_phytomer(ipts,ivm,:)
        WRITE(numout,*) 'yd: GLUC empty pft before, j=',ivm,'PHYbm(i,j)',PHYbm(ipts,ivm)
        IF (is_oilpalm(ivm)) THEN
           bm_phytomer(ipts,ivm,:) = zero
           bm_FFB(ipts,ivm,:) = zero
           PHYbm(ipts,ivm) = zero
           FFBbm(ipts,ivm) = zero
           WRITE(numout,*) 'yd: GLUC empty pft after, j=',ivm,'PHYbm(i,j)',PHYbm(ipts,ivm)
           WRITE(numout,*) 'yd: GLUC empty pft after, j=',ivm,'bm_phytomer(i,j,:)',bm_phytomer(ipts,ivm,:)
        ENDIF
     ENDIF
!! yidi
    litter(ipts,:,ivm,:,:) = zero
    fuel_1hr(ipts,ivm,:,:) = zero
    fuel_10hr(ipts,ivm,:,:) = zero
    fuel_100hr(ipts,ivm,:,:) = zero
    fuel_1000hr(ipts,ivm,:,:) = zero
    carbon(ipts,:,ivm) = zero 
    deepC_a(ipts,:,ivm) = zero 
    deepC_s(ipts,:,ivm) = zero 
    deepC_p(ipts,:,ivm) = zero 
    bm_to_litter(ipts,ivm,:,:) = zero
    DO ilev=1,nlevs
       lignin_struc(ipts,ivm,ilev) = zero
    ENDDO
    npp_longterm(ipts,ivm) = zero
    gpp_daily(ipts,ivm) = zero 
    gpp_week(ipts,ivm) = zero
    resp_maint(ipts,ivm) = zero
    resp_growth(ipts,ivm) = zero
    resp_hetero(ipts,ivm) = zero
    npp_daily(ipts,ivm) = zero
    co2_to_bm(ipts,ivm) = zero
    lm_lastyearmax(ipts,ivm) = zero
    age(ipts,ivm) = zero
    leaf_frac(ipts,ivm,:) = zero
    leaf_age(ipts,ivm,:) = zero
    everywhere(ipts,ivm) = zero
    when_growthinit(ipts,ivm) = zero
    gdd_from_growthinit(ipts,ivm) = zero
    gdd_midwinter(ipts,ivm) = zero
    time_hum_min(ipts,ivm) = zero
    gdd_m5_dormance(ipts,ivm) = zero
    ncd_dormance(ipts,ivm) = zero
    moiavail_month(ipts,ivm) = zero
    moiavail_week(ipts,ivm) = zero
    ngd_minus5(ipts,ivm) = zero
    PFTpresent(ipts,ivm) = .FALSE.
    senescence(ipts,ivm) = .FALSE.

  END SUBROUTINE empty_pft

  SUBROUTINE prepare_balance_check(outflux_sta,influx_sta,pool_sta,        &
                 veget_cov_max,                                            &
                 co2_to_bm,gpp_daily,npp_daily,                            &
                 biomass,litter,carbon,prod10,prod100,                     &
                 bm_to_litter,resp_maint,resp_growth,resp_hetero,          &
                 convflux,cflux_prod10,cflux_prod100,co2_fire)

    IMPLICIT NONE

    !! 0.1 Input variables
    REAL(r_std), DIMENSION(:,:), INTENT(in)            :: veget_cov_max        !! "Maximal" coverage fraction of a PFT (LAI 
                                                                                       !! -> infinity) on ground 
    REAL(r_std), DIMENSION(:,:), INTENT(in)              :: co2_to_bm            !! CO2 taken up from atmosphere when 
                                                                                       !! introducing a new PFT (introduced for 
                                                                                       !! carbon balance closure) 
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(in)            :: gpp_daily            !! Daily gross primary productivity  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(in)              :: npp_daily            !! Net primary productivity 
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in) :: biomass        !! Biomass @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in):: litter     !! Metabolic and structural litter, above 
                                                                                       !! and below ground 
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)      :: carbon               !! Carbon pool: active, slow, or passive, 
                                                                                       !! @tex $(gC m^{-2})$ @endtex  
    REAL(r_std),DIMENSION(:,:,:), INTENT(in)            :: prod10               !! Products remaining in the 10
                                                                                       !! year-turnover pool after the annual 
                                                                                       !! release for each compartment (10
                                                                                       !! + 1 : input from year of land cover 
                                                                                       !! change) @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(:,:,:), INTENT(in)           :: prod100              !! Products remaining in the 100 
                                                                                       !! year-turnover pool after the annual 
                                                                                       !! release for each compartment (100 
                                                                                       !! + 1 : input from year of land cover 
                                                                                       !! change) @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in):: bm_to_litter      !! Conversion of biomass to litter 
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(in)              :: resp_maint           !! Maintenance respiration  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    
    REAL(r_std), DIMENSION(:,:), INTENT(in)              :: resp_growth          !! Growth respiration  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(in)            :: resp_hetero          !! Heterotrophic respiration
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std),DIMENSION(:,:), INTENT(in)                 :: convflux             !! Release during first year following land 
                                                                                       !! cover change @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(:,:), INTENT(in)                 :: cflux_prod10         !! Total annual release from the 10 
                                                                                       !! year-turnover pool 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(:,:), INTENT(in)                 :: cflux_prod100        !! Total annual release from the 100 
                                                                                       !! year-turnover pool 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(in)              :: co2_fire             !! Carbon emitted into the atmosphere by 
                                                                                       !! fire (living and dead biomass)  


    !! 0.2 Output variables
    REAL(r_std), DIMENSION(:,:) :: outflux_sta
    REAL(r_std), DIMENSION(:,:) :: influx_sta
    REAL(r_std), DIMENSION(:,:) :: pool_sta


    !! 0.3 Modified variables

    !! 0.4 Local variables
    INTEGER(i_std) :: ind_biomass,ind_litter,ind_soil,ind_prod,ind_co2tobm,ind_gpp,ind_npp,&
                     ind_bm2lit,ind_resph,ind_respm,ind_respg,ind_convf,ind_cflux,ind_fire



    ind_biomass = 1
    ind_litter = 2
    ind_soil = 3
    ind_prod = 4
    pool_sta(:,ind_biomass) = SUM(SUM(biomass(:,:,:,icarbon),DIM=3) * veget_cov_max,DIM=2)
    pool_sta(:,ind_litter) = SUM(SUM(SUM(litter(:,:,:,:,icarbon),DIM=4),DIM=2) * veget_cov_max,DIM=2)
    pool_sta(:,ind_soil) = SUM(SUM(carbon(:,:,:),DIM=2) * veget_cov_max,DIM=2)
    pool_sta(:,ind_prod) = SUM(SUM(prod10,DIM=3),DIM=2) + SUM(SUM(prod100,DIM=3),DIM=2)

    ind_co2tobm = 1
    ind_gpp = 2
    ind_npp = 3
    influx_sta(:,ind_co2tobm) = SUM(co2_to_bm*veget_cov_max,DIM=2)
    influx_sta(:,ind_gpp) = SUM(gpp_daily*veget_cov_max,DIM=2)
    influx_sta(:,ind_npp) = SUM(npp_daily*veget_cov_max,DIM=2)

    ind_bm2lit = 1
    ind_respm = 2
    ind_respg = 3
    ind_resph = 4
    ind_convf = 5
    ind_cflux = 6
    ind_fire = 7
    outflux_sta(:,ind_bm2lit) = SUM(SUM(bm_to_litter(:,:,:,icarbon),DIM=3) * veget_cov_max,DIM=2)
    outflux_sta(:,ind_respm) = SUM(resp_maint*veget_cov_max,DIM=2)
    outflux_sta(:,ind_respg) = SUM(resp_growth*veget_cov_max,DIM=2)
    outflux_sta(:,ind_resph) = SUM(resp_hetero*veget_cov_max,DIM=2)
    outflux_sta(:,ind_convf) = SUM(convflux,DIM=2)
    outflux_sta(:,ind_cflux) = SUM(cflux_prod10,DIM=2) + SUM(cflux_prod100,DIM=2)
    outflux_sta(:,ind_fire) = SUM(co2_fire*veget_cov_max,DIM=2)

  END SUBROUTINE prepare_balance_check


  SUBROUTINE luc_balance_check(outflux_sta,influx_sta,pool_sta,            &
                 outflux_end,influx_end,pool_end,                          &
                 npts,lalo,identifier)

    IMPLICIT NONE

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                                 :: npts                 !! Domain size (unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(in)                    :: outflux_sta
    REAL(r_std), DIMENSION(:,:), INTENT(in)                    :: influx_sta
    REAL(r_std), DIMENSION(:,:), INTENT(in)                    :: pool_sta
    REAL(r_std), DIMENSION(:,:), INTENT(in)                    :: outflux_end
    REAL(r_std), DIMENSION(:,:), INTENT(in)                    :: influx_end
    REAL(r_std), DIMENSION(:,:), INTENT(in)                    :: pool_end
    REAL(r_std), DIMENSION(:,:), INTENT(in)                    :: lalo                 !! Geographical coordinates (latitude,longitude)
                                                                                       !! for pixels (degrees)
    CHARACTER(LEN=*), INTENT(in), OPTIONAL                     ::  identifier  !! string with to identify where this routine was called form

    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables
    REAL(r_std), DIMENSION(npts,nelements)     :: mass_before_2rd    !! Temporary variable 
    REAL(r_std), DIMENSION(npts,nelements)     :: mass_after_2rd    !! Temporary variable 
    REAL(r_std), DIMENSION(npts,nelements)     :: mass_change_2rd   !! positive
    REAL(r_std), DIMENSION(npts,nelements)     :: mass_balance_2rd    !! Temporary variable 
    INTEGER(i_std) :: ipts


    mass_before_2rd = zero
    mass_after_2rd = zero
    mass_change_2rd = zero
    mass_balance_2rd = zero

    mass_before_2rd(:,icarbon) = SUM(pool_sta(:,:),DIM=2)
    mass_after_2rd(:,icarbon) = SUM(pool_end(:,:),dim=2)

    ! mass_change_2rd is the mass increase
    mass_change_2rd(:,icarbon) = SUM(influx_end(:,:),DIM=2) - SUM(influx_sta(:,:),DIM=2) + &
                             + SUM(outflux_sta(:,:),DIM=2) - SUM(outflux_end(:,:),DIM=2)

    mass_balance_2rd(:,icarbon) = mass_before_2rd(:,icarbon) - mass_after_2rd(:,icarbon) &
                             + mass_change_2rd(:,icarbon)
     
    DO ipts = 1,npts
      IF (ABS(mass_balance_2rd(ipts,icarbon)) .GE. 1e-2) THEN
        WRITE (numout,*) ' FATAL Error'
        WRITE (numout,*) ' Mass conservation failed after ',identifier
        WRITE (numout,*) ' limit: '       , 1e-2
        WRITE (numout,*) ' Mismatch :'    , mass_balance_2rd(ipts,icarbon) 
        WRITE (numout,*) ' Coordinates :' , lalo(ipts,:) 
        WRITE (numout,*) ' gridpoint: '   , ipts , ' of ngrids: ',npts
        STOP
      ENDIF
    ENDDO
    WRITE (numout,*) 'mass balance check successful after ',identifier

  END SUBROUTINE luc_balance_check




  ! Note this subroutine does not depend on how many age classes there are
  ! in different MTCs.
  SUBROUTINE glcc_compensation_full(npts,veget_4veg,glcc,glccReal,glccDef, &
                               p2c,ipasture,g2c,igrass,f2c,itree,icrop,    &
                               IncreDeficit)

    IMPLICIT NONE

    !! 0.1 Input variables
    INTEGER, INTENT(in)                                         :: npts        !! Domain size - number of pixels (unitless)
    INTEGER, INTENT(in)    :: p2c,ipasture,g2c,igrass,f2c,itree,icrop
    REAL(r_std), DIMENSION (npts,12),INTENT(in)                 :: glcc        !! the land-cover-change (LCC) matrix in case a gross LCC is 
                                                                               !! used.

    !! 0.2 Output variables


    !! 0.3 Modified variables
    REAL(r_std), DIMENSION(npts,4), INTENT(inout)         :: veget_4veg        !! "maximal" coverage of tree/grass/pasture/crop
    REAL(r_std), DIMENSION(npts,12), INTENT(inout)        :: glccDef           !! Gross LCC deficit, negative values mean that there
                                                                               !! are not enough fractions in the source vegetations
                                                                               !! to the target ones as presribed by the LCC matrix.
    REAL(r_std), DIMENSION(npts,12), INTENT(inout)        :: glccReal          !! The "real" glcc matrix that we apply in the model
                                                                               !! after considering the consistency between presribed
                                                                               !! glcc matrix and existing vegetation fractions.
    REAL(r_std), DIMENSION(npts,4), INTENT(inout)         :: IncreDeficit      !! "Increment" deficits, negative values mean that 
                                                                               !! there are not enough fractions in the source PFTs
                                                                               !! /vegetations to target PFTs/vegetations. I.e., these
                                                                               !! fraction transfers are presribed in LCC matrix but
                                                                               !! not realized.
    
    !! 0.4 Local variables
    REAL(r_std), DIMENSION(npts)                          :: tmpdef            !! LCC deficits by summing up all the deficits to the
                                                                               !! the same target vegetation.

    !! 0. We first handle the cases where veget_4veg might be very small
    !tree
    WHERE(veget_4veg(:,itree) > min_stomate)
      glccDef(:,f2c) = veget_4veg(:,itree)-glcc(:,f2c)
      WHERE(veget_4veg(:,itree)>glcc(:,f2c))
        glccReal(:,f2c) = glcc(:,f2c)
      ELSEWHERE
        glccReal(:,f2c) = veget_4veg(:,itree)
      ENDWHERE
    ELSEWHERE
      glccReal(:,f2c) = 0.
      glccDef(:,f2c) = -1*glcc(:,f2c)
    ENDWHERE

    !pasture
    WHERE(veget_4veg(:,ipasture) > min_stomate)
      glccDef(:,p2c) = veget_4veg(:,ipasture)-glcc(:,p2c)
      WHERE(veget_4veg(:,ipasture)>glcc(:,p2c))
        glccReal(:,p2c) = glcc(:,p2c)
      ELSEWHERE
        glccReal(:,p2c) = veget_4veg(:,ipasture)
      ENDWHERE
    ELSEWHERE
      glccReal(:,p2c) = 0.
      glccDef(:,p2c) = -1*glcc(:,p2c)
    ENDWHERE

    !grass
    WHERE(veget_4veg(:,igrass) > min_stomate)
      glccDef(:,g2c) = veget_4veg(:,igrass)-glcc(:,g2c)
      WHERE(veget_4veg(:,igrass)>glcc(:,g2c))
        glccReal(:,g2c) = glcc(:,g2c)
      ELSEWHERE
        glccReal(:,g2c) = veget_4veg(:,igrass)
      ENDWHERE
    ELSEWHERE
      glccReal(:,g2c) = 0.
      glccDef(:,g2c) = -1*glcc(:,g2c)
    ENDWHERE

    !! 1. Compensation sequence: pasture,grass,forest 
    tmpdef(:) = glccDef(:,f2c)+glccDef(:,g2c)+glccDef(:,p2c)
    WHERE(glccDef(:,p2c)<0)
      WHERE(glccDef(:,g2c)<0)
        WHERE(glccDef(:,f2c)<0) ! 1 (-,-,-)
          IncreDeficit(:,icrop) = tmpdef(:)
        ELSEWHERE ! 2 (-,-,+)
          WHERE(tmpdef(:)>=min_stomate)
            glccReal(:,f2c) = glccReal(:,f2c)-glccDef(:,g2c)-glccDef(:,p2c)
          ELSEWHERE
            glccReal(:,f2c) = veget_4veg(:,itree)
            IncreDeficit(:,icrop) = tmpdef(:)
          ENDWHERE
        ENDWHERE
      ELSEWHERE
        WHERE(glccDef(:,f2c)<0) ! 3 (-,+,-)
          WHERE(tmpdef(:)>=min_stomate)
            glccReal(:,g2c) = glccReal(:,g2c)-glccDef(:,p2c)-glccDef(:,f2c)
          ELSEWHERE
            glccReal(:,g2c) = veget_4veg(:,igrass)
            IncreDeficit(:,icrop) = tmpdef(:)
          ENDWHERE
        ELSEWHERE ! 4 (-,+,+)
          WHERE(tmpdef(:)>=min_stomate)
            WHERE((glccDef(:,g2c)+glccDef(:,p2c))>=min_stomate)
              glccReal(:,g2c) = glccReal(:,g2c)-glccDef(:,p2c)
            ELSEWHERE
              glccReal(:,g2c) = veget_4veg(:,igrass)
              glccReal(:,f2c) = glccReal(:,f2c)-(glccDef(:,p2c)+glccDef(:,g2c))
            ENDWHERE
          ELSEWHERE
            glccReal(:,g2c) = veget_4veg(:,igrass)
            glccReal(:,f2c) = veget_4veg(:,itree)
            IncreDeficit(:,icrop) = tmpdef(:)
          ENDWHERE
        ENDWHERE
      ENDWHERE
    ELSEWHERE
      WHERE(glccDef(:,g2c)<0)
        WHERE(glccDef(:,f2c)<0) ! 5 (+,-,-)
          WHERE(tmpdef(:)>=min_stomate)
            glccReal(:,p2c) = glccReal(:,p2c)-glccDef(:,g2c)-glccDef(:,f2c)
          ELSEWHERE
            IncreDeficit(:,icrop) = tmpdef(:)
            glccReal(:,p2c) = veget_4veg(:,ipasture)
          ENDWHERE
        ELSEWHERE ! 6 (+,-,+)
          WHERE(tmpdef(:)>=min_stomate)
            WHERE((glccDef(:,p2c)+glccDef(:,g2c))>=min_stomate)
              glccReal(:,p2c) = glccReal(:,p2c)-glccDef(:,g2c)
            ELSEWHERE
              glccReal(:,p2c) = veget_4veg(:,ipasture)
              glccReal(:,f2c) = glccReal(:,f2c)-(glccDef(:,g2c)+glccDef(:,p2c))
            ENDWHERE
          ELSEWHERE
            IncreDeficit(:,icrop) = tmpdef(:)
            glccReal(:,p2c) = veget_4veg(:,ipasture)
            glccReal(:,f2c) = veget_4veg(:,itree)
          ENDWHERE
        ENDWHERE
      ELSEWHERE
        WHERE(glccDef(:,f2c)<0) ! 7 (+,+,-)
          WHERE(tmpdef(:)>=min_stomate)
            WHERE((glccDef(:,p2c)+glccDef(:,f2c))>=min_stomate)
              glccReal(:,p2c) = glccReal(:,p2c)-glccDef(:,f2c)
            ELSEWHERE
              glccReal(:,p2c) = veget_4veg(:,ipasture)
              glccReal(:,g2c) = glccReal(:,g2c)-(glccDef(:,f2c)+glccDef(:,p2c))
            ENDWHERE
          ELSEWHERE
            IncreDeficit(:,icrop) = tmpdef(:)
            glccReal(:,g2c) = veget_4veg(:,igrass)
            glccReal(:,p2c) = veget_4veg(:,ipasture)
          ENDWHERE
        ELSEWHERE ! 8 (+,+,+)
          !do nothing
        ENDWHERE
      ENDWHERE
    ENDWHERE
    veget_4veg(:,itree) = veget_4veg(:,itree) - glccReal(:,f2c)
    veget_4veg(:,igrass) = veget_4veg(:,igrass) - glccReal(:,g2c)
    veget_4veg(:,ipasture) = veget_4veg(:,ipasture) - glccReal(:,p2c)

  END SUBROUTINE glcc_compensation_full



  !! This subroutine implements non-full compensation, is currently
  !! abandoned.
  SUBROUTINE glcc_compensation(npts,veget_4veg,glcc,glccDef, &
                               p2c,ipasture,g2c,igrass,f2c,itree,icrop, &
                               IncreDeficit)

    IMPLICIT NONE

    !! 0.1 Input variables
    INTEGER, INTENT(in)                                         :: npts        !! Domain size - number of pixels (unitless)
    REAL(r_std), DIMENSION(npts,4), INTENT(in)                  :: veget_4veg  !! "maximal" coverage fraction of a PFT on the ground
    INTEGER, INTENT(in)    :: p2c,ipasture,g2c,igrass,f2c,itree,icrop

    !! 0.2 Output variables


    !! 0.3 Modified variables
    REAL(r_std), DIMENSION (npts,12),INTENT(inout)        :: glcc              !! the land-cover-change (LCC) matrix in case a gross LCC is 
                                                                               !! used.
    REAL(r_std), DIMENSION(npts,12), INTENT(inout)        :: glccDef           !! Gross LCC deficit, negative values mean that there
                                                                               !! are not enough fractions in the source vegetations
                                                                               !! to the target ones as presribed by the LCC matrix.
    REAL(r_std), DIMENSION(npts,4), INTENT(inout)         :: IncreDeficit      !! "Increment" deficits, negative values mean that 
                                                                               !! there are not enough fractions in the source PFTs
                                                                               !! /vegetations to target PFTs/vegetations. I.e., these
                                                                               !! fraction transfers are presribed in LCC matrix but
                                                                               !! not realized.
    
    !! 0.4 Local variables
    REAL(r_std), DIMENSION(npts)                          :: glccDef_all       !! LCC deficits by summing up all the deficits to the
                                                                               !! the same target vegetation.


    WHERE(veget_4veg(:,itree) > min_stomate)
      glccDef(:,f2c) = veget_4veg(:,itree)-glcc(:,f2c)
    ELSEWHERE
      glccDef(:,f2c) = -1*glcc(:,f2c)
      glcc(:,f2c) = 0.
    ENDWHERE

    WHERE(veget_4veg(:,ipasture) > min_stomate)
      glccDef(:,p2c) = veget_4veg(:,ipasture)-glcc(:,p2c)
    ELSEWHERE
      glccDef(:,p2c) = -1*glcc(:,p2c)
      glcc(:,p2c) = 0.
    ENDWHERE

    WHERE(veget_4veg(:,igrass) > min_stomate)
      glccDef(:,g2c) = veget_4veg(:,igrass)-glcc(:,g2c)
    ELSEWHERE
      glccDef(:,g2c) = -1*glcc(:,g2c)
      glcc(:,g2c) = 0.
    ENDWHERE

    glccDef_all(:) = glccDef(:,f2c)+glccDef(:,p2c)+glccDef(:,g2c)

    ! We allow the surpluses/deficits in p2c and g2c mutually compensating 
    ! for each other. If there are still deficits after this compensation,
    ! they will be further compensated for by the surpluses from f2c (if there are any
    ! surpluses). The ultimate deficits that cannot be compensated for 
    ! will be recorded and dropped. 

    ! Because we assume the "pasture rule" is used, i.e., the crops 
    ! are supposed to come primarily from pastures and grasses, normally
    ! we expect the deficits to occur in p2c or g2c rather than in f2c. But
    ! if it happens that f2c has deficits while p2c or g2c has surpluse,
    ! the surpluses will not be used to compensate for the f2c-deficits, 
    ! instead, we will just record and drop the f2c-deficits.

    ! In following codes for convenience we're not going to check
    ! whether surpluses in f2c are enough to compensate for deficits 
    ! in p2c or g2c or both. Instead, we just add their deficits on top
    ! of f2c. The issues of not-enough surpluses in f2c will be left for
    ! the codes after this section to handle.
    WHERE (glccDef(:,p2c) < 0.)
      glcc(:,p2c) = veget_4veg(:,ipasture)
      WHERE (glccDef(:,g2c) < 0.)
        glcc(:,g2c) = veget_4veg(:,igrass)
      ELSEWHERE
        WHERE (glccDef(:,g2c)+glccDef(:,p2c) > min_stomate)
          glcc(:,g2c) = glcc(:,g2c)-glccDef(:,p2c)
        ELSEWHERE
          glcc(:,g2c) = veget_4veg(:,igrass)
          ! whatever the case, we simply add the dificts to f2c
          glcc(:,f2c) = glcc(:,f2c)-glccDef(:,p2c)-glccDef(:,g2c)
        ENDWHERE
      ENDWHERE

    ELSEWHERE
      WHERE(glccDef(:,g2c) < 0.)
        glcc(:,g2c) = veget_4veg(:,igrass)
        WHERE(glccDef(:,p2c)+glccDef(:,g2c) > min_stomate)
          glcc(:,p2c) = glcc(:,p2c)-glccDef(:,g2c)
        ELSEWHERE
          glcc(:,p2c) = veget_4veg(:,ipasture)
          ! whatever the case, we simply add the dificts to f2c
          glcc(:,f2c) = glcc(:,f2c)-glccDef(:,p2c)-glccDef(:,g2c)
        ENDWHERE
      ELSEWHERE
        !Here p2c and g2c both show surplus, we're not going to check whether
        !glccDef(:,f2c) has negative values because we assume a "pasture rule"
        !is applied when constructing the gross LCC matrix, so deficits in 
        !f2c will just be dropped but not be compensated for by the surpluses in
        !p2c or g2c.
      ENDWHERE
    ENDWHERE

    ! 1. We calculate again the f2c-deficit because f2c-glcc is adjusted in the
    ! codes above as we allocated the deficits of p2c and g2c into f2c. 
    ! In cases where glccDef_all is less than zero, f2c-glcc will be larger
    ! than available forest veget_max and we therefore limit the f2c-glcc to
    ! available forest cover.
    ! 2. There is (probably) a second case where glccDef_all is larger then zero, 
    ! but f2c-glcc is higher than veget_tree, i.e., Originally f2c is given a 
    ! high value that there is deficit in f2c but surpluses exist for p2c and g2c.
    ! Normally we 
    ! assume this won't happen as explained above, given that a "pasture rule" was 
    ! used in constructing the gross LCC matrix. Nevertheless if this deos 
    ! happen, we will just drop the f2c deficit without being compensated 
    ! for by the surplus in p2c or g2c.
   
    ! we handle the 2nd case first
    WHERE(veget_4veg(:,itree) > min_stomate )
      WHERE(glccDef(:,f2c) < 0.)
        glcc(:,f2c) = veget_4veg(:,itree)
        WHERE (glccDef(:,p2c)+glccDef(:,g2c) > min_stomate)
          IncreDeficit(:,icrop) = glccDef(:,f2c)
        ELSEWHERE
          IncreDeficit(:,icrop) = glccDef_all(:)
        ENDWHERE
      ELSEWHERE
        WHERE(glccDef_all(:) < 0.) !handle the 1st case
          glcc(:,f2c) = veget_4veg(:,itree)
          IncreDeficit(:,icrop) = glccDef_all(:)
        ENDWHERE
      ENDWHERE
    ELSEWHERE
      WHERE(glccDef(:,p2c)+glccDef(:,g2c)>min_stomate)
        IncreDeficit(:,icrop) = glccDef(:,f2c)
      ELSEWHERE
        IncreDeficit(:,icrop) = glccDef_all(:)
      ENDWHERE
    ENDWHERE

  END SUBROUTINE glcc_compensation

END MODULE stomate_gluc_common
