! =================================================================================================================================
! MODULE       : stomate_lcchange_fh
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       This module is a copy of stomate_lcchange. It includes the forestry 
!              harvest.
!!
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): Including permafrost carbon
!!
!! REFERENCE(S)	: None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/perso/albert.jornet/ORCHIDEE-MICT/src_stomate/stomate_lcchange.f90 $
!! $Date: 2015-07-30 15:38:45 +0200 (Thu, 30 Jul 2015) $
!! $Revision: 2847 $
!! \n
!_ ================================================================================================================================


MODULE stomate_fharvest_MulAgeC

  ! modules used:
  
  USE ioipsl_para
  USE stomate_data
  USE pft_parameters
  USE constantes
  USE constantes_soil_var
  USE stomate_gluc_common
  USE stomate_glcchange_MulAgeC
  USE stomate_gluc_constants
  USE xios_orchidee
#ifdef CPP_PARA
  USE mpi
#endif

  IMPLICIT NONE

  PRIVATE
  PUBLIC fharvest_MulAgeC, Get_harvest_matrix
  
CONTAINS

! ================================================================================================================================
!! SUBROUTINE   gross_lcchange
!!
!>\BRIEF       : Apply gross land cover change.
!!
!>\DESCRIPTION  
!_ ================================================================================================================================
  SUBROUTINE fharvest_MulAgeC (npts, dt_days, harvest_matrix,newvegfrac,   &
               fuelfrac, &
               def_fuel_1hr_remain, def_fuel_10hr_remain,        &
               def_fuel_100hr_remain, def_fuel_1000hr_remain,    &
               deforest_litter_remain, deforest_biomass_remain,  &
               convflux,                   &
               glcc_pft, glcc_pftmtc,          &
               veget_max, prod10, prod100,             &
               PFTpresent, senescence, moiavail_month, moiavail_week,  &
               gpp_week, ngd_minus5, resp_maint, resp_growth,          & 
               resp_hetero, npp_daily, when_growthinit, npp_longterm,  &
               ind, lm_lastyearmax, everywhere, age,                   &
               co2_to_bm, gpp_daily, co2_fire,                         &
               time_hum_min, gdd_midwinter, gdd_from_growthinit,       &
               gdd_m5_dormance, ncd_dormance,                          &
               lignin_struc, carbon, leaf_frac,                        &
               deepC_a, deepC_s, deepC_p,                              &
               leaf_age, bm_to_litter, biomass, litter,                &
               bm_phytomer,bm_FFB,PHYbm, FFBbm,                  & !! yidi
               fuel_1hr, fuel_10hr, fuel_100hr, fuel_1000hr)
  
    IMPLICIT NONE

    !! 0.1 Input variables

    INTEGER, INTENT(in)                                  :: npts             !! Domain size - number of pixels (unitless)
    REAL(r_std), INTENT(in)                              :: dt_days          !! Time step of vegetation dynamics for stomate
    REAL(r_std), DIMENSION (npts,12),INTENT(in)          :: harvest_matrix             !! 
    REAL(r_std), DIMENSION (npts),INTENT(in)             :: fuelfrac             !! 
    REAL(r_std), DIMENSION (npts,nvmap), INTENT(in)      :: newvegfrac
                                                                             !! 

    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements), INTENT(in)                 :: def_fuel_1hr_remain
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements), INTENT(in)                 :: def_fuel_10hr_remain
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements), INTENT(in)                 :: def_fuel_100hr_remain
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements), INTENT(in)                 :: def_fuel_1000hr_remain
    REAL(r_std), DIMENSION(npts,nlitt,nvm,nlevs,nelements), INTENT(in) :: deforest_litter_remain   !! Vegetmax-weighted remaining litter on the ground for 
                                                                                                      !! deforestation region.
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(in)      :: deforest_biomass_remain  !! Vegetmax-weighted remaining biomass on the ground for 
                                                                                                      !! deforestation region.


    !! 0.2 Output variables
    REAL(r_std), DIMENSION(npts,nwp), INTENT(inout)            :: convflux         !! release during first year following land cover
                                                                             !! change
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)      :: glcc_pft         !! Loss of fraction in each PFT
    REAL(r_std), DIMENSION(npts,nvm,nvmap), INTENT(inout):: glcc_pftmtc      !! a temporary variable to hold the fractions each PFT is going to lose
                                                                             !! i.e., the contribution of each PFT to the youngest age-class of MTC

    !! 0.3 Modified variables
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)      :: veget_max        !! "maximal" coverage fraction of a PFT (LAI ->
                                                                             !! infinity) on ground (unitless) 
    REAL(r_std), DIMENSION(npts,0:10,nwp), INTENT(inout)     :: prod10           !! products remaining in the 10 year-turnover
                                                                             !! pool after the annual release for each 
                                                                             !! compartment (10 + 1 : input from year of land
                                                                             !! cover change)
    REAL(r_std), DIMENSION(npts,0:100,nwp), INTENT(inout)    :: prod100          !! products remaining in the 100 year-turnover
                                                                             !! pool after the annual release for each 
                                                                             !! compartment (100 + 1 : input from year of land
                                                                             !! cover change)
                                                                             !! pool compartments
    LOGICAL, DIMENSION(:,:), INTENT(inout)               :: PFTpresent       !! Tab indicating which PFTs are present in 
                                                                             !! each pixel
    LOGICAL, DIMENSION(:,:), INTENT(inout)               :: senescence       !! Flag for setting senescence stage (only 
                                                                             !! for deciduous trees)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)           :: moiavail_month   !! "Monthly" moisture availability (0 to 1, 
                                                                             !! unitless) 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)           :: moiavail_week    !! "Weekly" moisture availability 
                                                                             !! (0 to 1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)           :: gpp_week         !! Mean weekly gross primary productivity 
                                                                             !! @tex $(gC m^{-2} day^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)           :: ngd_minus5       !! Number of growing days (days), threshold 
                                                                             !! -5 deg C (for phenology)   
    REAL(r_std), DIMENSION(:,:), INTENT(inout)           :: resp_maint       !! Maintenance respiration  
                                                                             !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)           :: resp_growth      !! Growth respiration  
                                                                             !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)           :: resp_hetero      !! Heterotrophic respiration  
                                                                             !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)           :: npp_daily        !! Net primary productivity 
                                                                             !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)           :: when_growthinit  !! How many days ago was the beginning of 
                                                                             !! the growing season (days)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)           :: npp_longterm     !! "Long term" mean yearly primary productivity 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)           :: ind              !! Number of individuals at the stand level
                                                                             !! @tex $(m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)           :: lm_lastyearmax   !! last year's maximum leaf mass for each PFT 
                                                                             !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)           :: everywhere       !! is the PFT everywhere in the grid box or 
                                                                             !! very localized (after its introduction) (?)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)           :: age              !! mean age (years)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)           :: co2_to_bm        !! CO2 taken from the atmosphere to get C to create  
                                                                             !! the seedlings @tex (gC.m^{-2}dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)           :: gpp_daily        !! Daily gross primary productivity  
                                                                             !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)           :: co2_fire         !! Fire carbon emissions
                                                                             !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)           :: time_hum_min     !! Time elapsed since strongest moisture 
                                                                             !! availability (days) 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)           :: gdd_midwinter    !! Growing degree days (K), since midwinter 
                                                                             !! (for phenology) - this is written to the
                                                                             !!  history files 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)           :: gdd_from_growthinit !! growing degree days, since growthinit 
                                                                             !! for crops
    REAL(r_std), DIMENSION(:,:), INTENT(inout)           :: gdd_m5_dormance  !! Growing degree days (K), threshold -5 deg 
                                                                             !! C (for phenology)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)           :: ncd_dormance     !! Number of chilling days (days), since 
                                                                             !! leaves were lost (for phenology) 
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)         :: lignin_struc     !! ratio Lignine/Carbon in structural litter,
                                                                             !! above and below ground
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)         :: carbon           !! carbon pool: active, slow, or passive 
                                                                             !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)         :: deepC_a          !! Permafrost soil carbon (g/m**3) active
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)         :: deepC_s          !! Permafrost soil carbon (g/m**3) slow
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)         :: deepC_p          !! Permafrost soil carbon (g/m**3) passive
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)         :: leaf_frac        !! fraction of leaves in leaf age class (unitless;0-1)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)         :: leaf_age         !! Leaf age (days)
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)       :: bm_to_litter     !! Transfer of biomass to litter 
                                                                             !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)       :: biomass          !! Stand level biomass @tex $(gC.m^{-2})$ @endtex
!! yidi
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: FFBbm               !! FFB mass, from sapabove
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: PHYbm               !! PHYTOMER mass, from sapabove
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: bm_phytomer         !! Each PHYTOMER mass, from sapabove
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: bm_FFB              !! Fruit mass for Each PHYTOMER
!! yidi
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)     :: litter           !! metabolic and structural litter, above and 
                                                                             !! below ground @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)       :: fuel_1hr
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)       :: fuel_10hr
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)       :: fuel_100hr
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)       :: fuel_1000hr

    !! 0.4 Local variables
    REAL(r_std), DIMENSION(nvmap,nparts,nelements)             :: bm_to_litter_pro !! conversion of biomass to litter 
                                                                             !! @tex ($gC m^{-2} day^{-1}$) @endtex
    REAL(r_std), DIMENSION(nvmap,nparts,nelements)             :: biomass_pro      !! biomass @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(nvmap)                        :: veget_max_pro    !! "maximal" coverage fraction of a PFT (LAI ->
                                                                             !! infinity) on ground (unitless) 
    REAL(r_std), DIMENSION(nvmap,ncarb)                        :: carbon_pro       !! carbon pool: active, slow, or passive 
                                                                             !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(nvmap,ndeep)                        :: deepC_a_pro      !! Permafrost carbon pool: active, slow, or passive 
                                                                             !! @tex ($gC m^{-3}$) @endtex
    REAL(r_std), DIMENSION(nvmap,ndeep)                        :: deepC_s_pro      !! Permafrost carbon pool: active, slow, or passive 
                                                                             !! @tex ($gC m^{-3}$) @endtex
    REAL(r_std), DIMENSION(nvmap,ndeep)                        :: deepC_p_pro      !! Permafrost carbon pool: active, slow, or passive 
                                                                             !! @tex ($gC m^{-3}$) @endtex
    REAL(r_std), DIMENSION(nvmap,nlitt,nlevs,nelements)        :: litter_pro       !! metabolic and structural litter, above and 
                                                                             !! below ground @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(nvmap,nlitt,nelements)              :: fuel_1hr_pro
    REAL(r_std), DIMENSION(nvmap,nlitt,nelements)              :: fuel_10hr_pro
    REAL(r_std), DIMENSION(nvmap,nlitt,nelements)              :: fuel_100hr_pro
    REAL(r_std), DIMENSION(nvmap,nlitt,nelements)              :: fuel_1000hr_pro
    REAL(r_std), DIMENSION(nvmap,nlevs)                        :: lignin_struc_pro !! ratio Lignine/Carbon in structural litter
                                                                             !! above and below ground
    REAL(r_std), DIMENSION(nvmap,nleafages)                    :: leaf_frac_pro    !! fraction of leaves in leaf age class 
    REAL(r_std), DIMENSION(nvmap,nleafages)                    :: leaf_age_pro     !! fraction of leaves in leaf age class 
    LOGICAL, DIMENSION(nvmap)                :: PFTpresent_pro, senescence_pro                 !! Is pft there (unitless)
    REAL(r_std), DIMENSION(nvmap)            :: ind_pro, age_pro, lm_lastyearmax_pro, npp_longterm_pro
    REAL(r_std), DIMENSION(nvmap)            :: everywhere_pro
    REAL(r_std), DIMENSION(nvmap)            :: gpp_daily_pro, npp_daily_pro, co2_to_bm_pro
    REAL(r_std), DIMENSION(nvmap)            :: resp_maint_pro, resp_growth_pro
    REAL(r_std), DIMENSION(nvmap)            :: resp_hetero_pro, co2_fire_pro
  
    INTEGER                :: ipts,ivm,ivma,l,m,ipft_young_agec
    CHARACTER(LEN=10)      :: part_str                               !! string suffix indicating an index

    REAL(r_std), DIMENSION(npts,nvmap)       :: glcc_mtc             !! Increase in fraction of each MTC in its youngest age-class
    REAL(r_std), DIMENSION(npts,nvm)         :: glccReal_tmp         !! A temporary variable to hold glccReal
    REAL(r_std), DIMENSION(npts)             :: Deficit_pf2yf_final     !! 
    REAL(r_std), DIMENSION(npts)             :: Deficit_sf2yf_final     !! 
    REAL(r_std), DIMENSION(npts)             :: pf2yf_compen_sf2yf      !! 
    REAL(r_std), DIMENSION(npts)             :: sf2yf_compen_pf2yf      !!
    REAL(r_std), DIMENSION(npts,nvm)         :: glcc_harvest            !! Loss of fraction due to forestry harvest
    REAL(r_std), DIMENSION(npts)             :: harvest_wood            !! Loss of fraction due to forestry harvest

    WRITE(numout,*) 'Entering fharvest_MulAgeC'
    glcc_harvest(:,:) = zero
    harvest_wood(:) = zero
    glccReal_tmp(:,:) = zero

    CALL MulAgeC_fh_firstday(npts,veget_max,newvegfrac,harvest_matrix, &
                          glcc_pft,glcc_pftmtc, &
                          Deficit_pf2yf_final, Deficit_sf2yf_final, &
                          pf2yf_compen_sf2yf, sf2yf_compen_pf2yf)

    glcc_mtc(:,:) = SUM(glcc_pftmtc,DIM=2)

    ! Write out forestry harvest variables
    DO ipts = 1,npts
      DO ivm = 1,nvm
        DO ivma = 1,nvmap
          IF (is_tree(ivm) .AND. is_tree(start_index(ivma))) THEN
            glcc_harvest(ipts,ivm) = glcc_harvest(ipts,ivm) + glcc_pftmtc(ipts,ivm,ivma)
          ENDIF
        ENDDO
        harvest_wood(ipts) = harvest_wood(ipts)+ (biomass(ipts,ivm,isapabove,icarbon)+ &
               biomass(ipts,ivm,iheartabove,icarbon))*glcc_harvest(ipts,ivm)
      ENDDO
    ENDDO

    DO ipts=1,npts

      !! Initialize the _pro variables
      bm_to_litter_pro(:,:,:)=zero                                                
      biomass_pro(:,:,:)=zero
      veget_max_pro(:)=zero                                                       
      carbon_pro(:,:)=zero                                                        
      deepC_a_pro(:,:)=zero                                                       
      deepC_s_pro(:,:)=zero                                                       
      deepC_p_pro(:,:)=zero                                                       
      litter_pro(:,:,:,:)=zero                                                    
      fuel_1hr_pro(:,:,:)=zero                                                    
      fuel_10hr_pro(:,:,:)=zero                                                   
      fuel_100hr_pro(:,:,:)=zero                                                  
      fuel_1000hr_pro(:,:,:)=zero                                                 
      lignin_struc_pro(:,:)=zero                                                  

      leaf_frac_pro = zero
      leaf_age_pro = zero
      PFTpresent_pro(:) = .FALSE.
      senescence_pro(:) = .TRUE.
      ind_pro = zero
      age_pro = zero
      lm_lastyearmax_pro = zero
      npp_longterm_pro = zero
      everywhere_pro = zero
      
      gpp_daily_pro=zero                                                       
      npp_daily_pro=zero                                                       
      co2_to_bm_pro=zero                                                       
      resp_maint_pro=zero                                                      
      resp_growth_pro=zero                                                     
      resp_hetero_pro=zero                                                     
      co2_fire_pro=zero                                                        

      ! Note that we assume people don't intentionally change baresoil to 
      ! vegetated land.
      DO ivma = 2,nvmap
        ! we assume only the youngest age class receives the incoming PFT
        ! [chaoyuejoy@gmail.com 2015-08-04] This line is commented to allow
        ! the case of only single age class being handled.

        ! here we set glcc_mtc(ipts,ivma) > min_stomate as a condition,
        ! this is necessary because later on in the subroutine of 
        ! add_incoming_proxy_pft we have to merge the newly established
        ! youngest proxy with potentially exisiting youngest receiving MTC,
        ! thus have to devide a new fraction of (frac_proxy + frac_exist),
        ! but in case frac_exist = zero, we risk deviding by a very small value
        ! of frac_proxy and thus we want it to be bigger than min_stomate.
        IF ( glcc_mtc(ipts,ivma) .GT. min_stomate ) THEN

          ! 1. we accumulate the scalar variables that will be inherited
          !    note we don't handle the case of harvesting forest because
          !    we assume glcc_pftmtc(forest->forest) would be zero and this
          !    case won't occur as it's filtered by the condition of
          !    (frac>min_stomate)
          CALL collect_legacy_pft_forestry(npts, ipts, ivma, glcc_pftmtc,    &
                  fuelfrac, &
                  biomass, bm_to_litter, carbon, litter,            &
                  deepC_a, deepC_s, deepC_p,                        &
                  fuel_1hr, fuel_10hr, fuel_100hr, fuel_1000hr,     &
                  lignin_struc, co2_to_bm, gpp_daily, npp_daily,    &
                  resp_maint, resp_growth, resp_hetero, co2_fire,   &
                  def_fuel_1hr_remain, def_fuel_10hr_remain,        &
                  def_fuel_100hr_remain, def_fuel_1000hr_remain,    &
                  deforest_litter_remain, deforest_biomass_remain,  &
                  veget_max_pro(ivma), carbon_pro(ivma,:),          &
                  lignin_struc_pro(ivma,:), litter_pro(ivma,:,:,:), &
                  deepC_a_pro(ivma,:), deepC_s_pro(ivma,:), deepC_p_pro(ivma,:), &
                  fuel_1hr_pro(ivma,:,:), fuel_10hr_pro(ivma,:,:),  &
                  fuel_100hr_pro(ivma,:,:), fuel_1000hr_pro(ivma,:,:), &
                  bm_to_litter_pro(ivma,:,:), co2_to_bm_pro(ivma),  &
                  gpp_daily_pro(ivma), npp_daily_pro(ivma),         &
                  resp_maint_pro(ivma), resp_growth_pro(ivma),      &
                  resp_hetero_pro(ivma), co2_fire_pro(ivma),        &
                  convflux,prod10,prod100)

          !++TEMP++
          ! Here we substract the outgoing fraction from the source PFT.
          ! If a too small fraction remains in this source PFT, then it is
          ! exhausted, we empty it. The subroutine 'empty_pft' might be 
          ! combined with 'collect_legacy_pft', but now we just put it here.
          DO ivm = 1,nvm
            IF( glcc_pftmtc(ipts,ivm,ivma)>zero ) THEN
              veget_max(ipts,ivm) = veget_max(ipts,ivm)-glcc_pftmtc(ipts,ivm,ivma)
              IF ( veget_max(ipts,ivm)<min_stomate ) THEN
                CALL empty_pft(ipts, ivm, veget_max, biomass, ind,       &
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
              ENDIF
            ENDIF
          ENDDO

        ENDIF !IF ( glcc_mtc(ipts,ivma) .GT. min_stomate )
      ENDDO !(DO ivma = 2,nvmap)

      ! We can only establish new youngest proxy and add it to the 
      ! existing youngest-age PFT after all the harvest is done, to 
      ! avoid the dilution of harvestable biomass by the young proxy
      ! and ensure consistency. Therefore now we have to loop again
      ! over nvmap.
      DO ivma = 2,nvmap
        IF ( glcc_mtc(ipts,ivma) .GT. min_stomate ) THEN

          ipft_young_agec = start_index(ivma)

          ! 2. we establish a proxy PFT with the fraction of veget_max_pro,
          !    which is going to be either merged with existing target 
          !    `ipft_young_agec` PFT, or fill the place if no existing target PFT
          !    exits.
          CALL initialize_proxy_pft(ipts,ipft_young_agec,veget_max_pro(ivma), &
                 biomass_pro(ivma,:,:), co2_to_bm_pro(ivma), ind_pro(ivma),   &
                 age_pro(ivma),                                               & 
                 senescence_pro(ivma), PFTpresent_pro(ivma),                  &
                 lm_lastyearmax_pro(ivma), everywhere_pro(ivma),              &
                 npp_longterm_pro(ivma),                                      &
                 leaf_frac_pro(ivma,:),leaf_age_pro(ivma,:))

          CALL sap_take (ipts,ivma,veget_max,biomass_pro(ivma,:,:), &
                         biomass,co2_to_bm_pro(ivma))

          ! 3. we merge the newly initiazlized proxy PFT into existing one
          !    or use it to fill an empty PFT slot.
          CALL add_incoming_proxy_pft(npts, ipts, ipft_young_agec, veget_max_pro(ivma),&
                 carbon_pro(ivma,:), litter_pro(ivma,:,:,:), lignin_struc_pro(ivma,:), &
                 bm_to_litter_pro(ivma,:,:),    &
                 deepC_a_pro(ivma,:), deepC_s_pro(ivma,:), deepC_p_pro(ivma,:), &
                 fuel_1hr_pro(ivma,:,:), fuel_10hr_pro(ivma,:,:),               &
                 fuel_100hr_pro(ivma,:,:), fuel_1000hr_pro(ivma,:,:),           &
                 biomass_pro(ivma,:,:), co2_to_bm_pro(ivma),                    &
                 npp_longterm_pro(ivma), ind_pro(ivma),                         &
                 lm_lastyearmax_pro(ivma), age_pro(ivma), everywhere_pro(ivma), &  
                 leaf_frac_pro(ivma,:), leaf_age_pro(ivma,:),                   &
                 PFTpresent_pro(ivma), senescence_pro(ivma),                &
                 gpp_daily_pro(ivma), npp_daily_pro(ivma),                      &
                 resp_maint_pro(ivma), resp_growth_pro(ivma),                   &
                 resp_hetero_pro(ivma), co2_fire_pro(ivma),                     &
                 veget_max, carbon, litter, lignin_struc, bm_to_litter,         &
                 deepC_a, deepC_s, deepC_p,                                     &
                 fuel_1hr, fuel_10hr, fuel_100hr, fuel_1000hr,                  &
                 biomass, co2_to_bm, npp_longterm, ind,                         &
                 lm_lastyearmax, age, everywhere,                               &
                 leaf_frac, leaf_age, PFTpresent, senescence,                   &
                 gpp_daily, npp_daily, resp_maint, resp_growth,                 &
                 resp_hetero, co2_fire)
          
        ENDIF !IF ( glcc_mtc(ipts,ivma) .GT. min_stomate )
      ENDDO !(DO ivma=1,nvmap)

    ENDDO !(DO ipts=1,npts)


    CALL histwrite_p (hist_id_stomate, 'glcc_harvest', itime, &
         glcc_harvest, npts*nvm, horipft_index)
    CALL xios_orchidee_send_field ('glcc_harvest', glcc_harvest) ! kjpindex,nvm

    CALL histwrite_p (hist_id_stomate, 'harvest_wood', itime, &
         harvest_wood, npts, hori_index)
    CALL xios_orchidee_send_field ('harvest_wood', harvest_wood) ! kjpindex

    glccReal_tmp(:,:) = zero
    glccReal_tmp(:,1) = Deficit_pf2yf_final
    glccReal_tmp(:,2) = Deficit_sf2yf_final
    glccReal_tmp(:,3) = pf2yf_compen_sf2yf
    glccReal_tmp(:,4) = sf2yf_compen_pf2yf ! this is zero as currently the deficit
                                           ! in primary forest harvest is not 
                                           ! compensated for by the surplus in 
                                           ! secondary forest harvest.

    CALL histwrite_p (hist_id_stomate, 'DefiComForHarvest', itime, &
         glccReal_tmp, npts*nvm, horipft_index)
    CALL xios_orchidee_send_field ('DefiComForHarvest', glccReal_tmp) ! kjpindex,nvm

  END SUBROUTINE fharvest_MulAgeC

! ================================================================================================================================
!! SUBROUTINE   : gross_lcc_firstday
!!
!>\BRIEF        : When necessary, adjust input glcc matrix, and allocate it
!!                into different contributing age classes and receiving 
!!                youngest age classes.
!! \n
!_ ================================================================================================================================

  ! Note: it has this name because this subroutine will also be called
  ! the first day of each year to precalculate the forest loss for the
  ! deforestation fire module.
  SUBROUTINE MulAgeC_fh_firstday(npts,veget_max_org,newvegfrac,harvest_matrix, &
                          glcc_pft,glcc_pftmtc, &
                          Deficit_pf2yf_final, Deficit_sf2yf_final, &
                          pf2yf_compen_sf2yf, sf2yf_compen_pf2yf)

    IMPLICIT NONE

    !! 0.1 Input variables

    INTEGER, INTENT(in)                                     :: npts           !! Domain size - number of pixels (unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)            :: veget_max_org  !! "maximal" coverage fraction of a PFT on the ground
                                                                              !! May sum to
                                                                              !! less than unity if the pixel has
                                                                              !! nobio area. (unitless, 0-1)
    REAL(r_std), DIMENSION(npts,nvmap), INTENT(in)          :: newvegfrac     !! used to guid the allocation of new PFTs.
    REAL(r_std), DIMENSION(npts,12),INTENT(in)              :: harvest_matrix !! 
                                                                              !! 
    !! 0.2 Output variables
    REAL(r_std), DIMENSION(npts,nvm,nvmap), INTENT(inout)   :: glcc_pftmtc    !! the fractions each PFT is going to lose
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)         :: glcc_pft       !! Loss of fraction in each PFT

    REAL(r_std), DIMENSION(npts), INTENT(inout)    :: Deficit_pf2yf_final     !! 
    REAL(r_std), DIMENSION(npts), INTENT(inout)    :: Deficit_sf2yf_final     !! 
    REAL(r_std), DIMENSION(npts), INTENT(inout)    :: pf2yf_compen_sf2yf      !! 
    REAL(r_std), DIMENSION(npts), INTENT(inout)    :: sf2yf_compen_pf2yf      !!

    !! 0.3 Modified variables
    
    !! 0.4 Local variables
    REAL(r_std), DIMENSION(npts,nvmap)              :: veget_mtc           !! "maximal" coverage fraction of a PFT on the ground
    REAL(r_std), DIMENSION(npts,nvmap)              :: veget_mtc_begin     !! "maximal" coverage fraction of a PFT on the ground
    REAL(r_std), DIMENSION(npts,nagec_tree)         :: vegagec_tree        !! fraction of tree age-class groups, in sequence of old->young
    REAL(r_std), DIMENSION(npts,nagec_herb)         :: vegagec_grass       !! fraction of grass age-class groups, in sequence of old->young 
    REAL(r_std), DIMENSION(npts,nagec_herb)         :: vegagec_pasture     !! fraction of pasture age-class groups, in sequence of old->young
    REAL(r_std), DIMENSION(npts,nagec_herb)         :: vegagec_crop        !! fraction of crop age-class groups, in sequence of old->young

    
    REAL(r_std), DIMENSION(npts,4)                  :: veget_4veg      !! "maximal" coverage fraction of a PFT on the ground
    REAL(r_std), DIMENSION(npts)                    :: veget_tree      !! "maximal" coverage fraction of a PFT on the ground
    REAL(r_std), DIMENSION(npts)                    :: veget_grass     !! "maximal" coverage fraction of a PFT on the ground
    REAL(r_std), DIMENSION(npts)                    :: veget_pasture   !! "maximal" coverage fraction of a PFT on the ground
    REAL(r_std), DIMENSION(npts)                    :: veget_crop      !! "maximal" coverage fraction of a PFT on the ground

    REAL(r_std), DIMENSION(npts,nvm)         :: veget_max              !! "maximal" coverage fraction of a PFT on the ground
    REAL(r_std), DIMENSION(npts,nvm)         :: veget_max_tmp          !! "maximal" coverage fraction of a PFT on the ground
    REAL(r_std), DIMENSION(npts,nvm)         :: veget_max_old          !! "maximal" coverage fraction of a PFT on the ground
    REAL(r_std), DIMENSION(npts,nvm)         :: glcc_pft_tmp           !! Loss of fraction in each PFT

    REAL(r_std), DIMENSION(npts,nvm,nvmap)  :: glcc_pftmtc_harvest    !! a temporary variable to hold the fractions each PFT is going to lose
    REAL(r_std), DIMENSION(npts,12)         :: glccRealHarvest    !! real matrix applied for forestry harvest.

    REAL(r_std), DIMENSION(npts,nvm)         :: glccReal_tmp         !! A temporary variable to hold glccReal

    LOGICAL, SAVE  :: MulAgeC_fh_firstday_done = .FALSE.

    INTEGER :: ivma, pf2yf, sf2yf, ivm

    REAL(r_std), DIMENSION(npts,12)         :: FHmatrix_remainA        !! 
    REAL(r_std), DIMENSION(npts,12)         :: FHmatrix_remainB        !! 
    REAL(r_std), DIMENSION(npts,12)         :: glccRemain              !! 

    INTEGER :: ipts,IndStart_f,IndEnd_f
    CHARACTER(LEN=10)      :: part_str                               !! string suffix indicating an index

    !Some more local configurations
    LOGICAL                                 :: allow_youngest_forest_WoodHarvest = .FALSE.
    
 
    ! we make copies of original input veget_max (which is veget_max_org
    ! in the subroutine parameter list).
    ! veget_max will be modified through different operations in order to 
    ! check for various purposes, e.g., whether input harvest matrix 
    ! is compatible with existing veget_max and how to allocate it etc.
    ! veget_max_old will not be modified
    veget_max(:,:) = veget_max_org(:,:)
    veget_max_old(:,:) = veget_max_org(:,:)

    !********************** block to handle forestry harvest ****************
    !! 1. Handle the forestry harvest process

    !! 1.0 Some preparation 

    pf2yf=1   !primary to young forest conversion because of harvest
    sf2yf=2   !old secondary to young forest conversion because of harvest
    
    Deficit_pf2yf_final(:) = zero 
    Deficit_sf2yf_final(:) = zero

    ! Note in the naming of pf2yf_compen_sf2yf and sf2yf_compen_pf2yf, active
    ! tense is used. I.e., pf2yf_compen_sf2yf means the fraction which pf2yf
    ! compenstates for sf2yf
    pf2yf_compen_sf2yf(:) = zero  !primary->young conversion that compensates 
                               !the secondary->young conversion because of deficit
                               !in the latter
    sf2yf_compen_pf2yf(:) = zero  !seondary->young conversion that compensates
                               !the primary->young conversion because of the deficit
                               !in the latter

    ! we now have to fill the transtion of forest->forest in the process of harvest
    ! into our target matrix glcc_pftmtc. Thus we will initiliaze them first.
    glcc_pft(:,:) = 0.
    glcc_pft_tmp(:,:) = 0.
    glcc_pftmtc(:,:,:) = 0.
    glccRemain(:,:) = harvest_matrix(:,:)

    !! 2.1 Handle secondary forest harvest

    CALL calc_cover(npts,veget_max,veget_mtc,vegagec_tree,vegagec_grass, &
           vegagec_pasture,vegagec_crop)
    
    ! In following call of calc_cover, veget_mtc will be updated each time,
    ! but we don't want this, so we put its initial value into veget_mtc_begin
    ! in order to retrieve this initial value later. veget_mtc is used in the 
    ! subroutine type_conversion to guild the allocation of newly created
    ! lands of a certain vegetation type (e.g., pasture) into its compoenent
    ! meta-classes (e.g. C3 and C4 pasture).
    veget_mtc_begin = veget_mtc

    ! Allocate harvest-caused out-going primary and secondary forest fraction
    ! into different primary and secondary (all other younger age classes) forest PFTs. 

    ! [Note] 
    ! Below we used the tempelate of type_conversion but in fact we need
    ! only glcc_pft, which means the fraction loss in each PFT. We then need to
    ! use glcc_pft to fill glcc_pftmtc (our final target matrix), assuming that
    ! the loss of forest PFT will go to the youngest age class of its forest MTC.
    ! Therefore we assume no changes of forset species (meta-class) in forestry 
    ! harvest and following tree-replanting.
    ! Although glcc_pftmtc and glcc_pft_tmp will be automatically filled when 
    ! we use the tempelate `type_conversion` by calling it as below, but it makes
    ! no sense because they will be reset later.

    !! 1.1 Secondary forest harvest.

    ! We first handle within the secondary forest age classes, in the sequence 
    ! of old->young

    !IndStart_f = 2             ! note the indecies for tree age class are in the
    IndStart_f = 1             ! note the indecies for tree age class are in the
                               ! sequence of from old to young, thus index=2 means the
                               ! 2nd oldest age class.
    !IndEnd_f = nagec_tree-1    ! nagec_tree-2: The 3rd youngest age class
    IndEnd_f = nagec_tree-5    ! nagec_tree-2: The 3rd youngest age class
                               ! nagec_tree-1: The 2nd youngest age class
                               ! nagec_tree: The youngest age class

    IF (IndStart_f .LE. 0 .OR. IndEnd_f .LE. 0) THEN
      write(numout,*) 'Forest harvest: Age class index cannot be negative or zero!'
      STOP
    ENDIF

    DO ipts=1,npts
      !sf2yf
      CALL type_conversion(ipts,sf2yf,harvest_matrix,veget_mtc,newvegfrac, &
                       indold_tree,indagec_tree,indagec_pasture,num_pasture_mulagec,&
                       IndEnd_f,nagec_herb,                    &
                       vegagec_tree,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain, &
                       .TRUE., iagec_start=IndStart_f)
    ENDDO
    FHmatrix_remainA(:,:) = glccRemain

    !! 1.2 Use primary forest harvest to compensate the deficit in secondary 
    !!       forest harvest. Note such compensation happens automatically here.

    CALL calc_cover(npts,veget_max,veget_mtc,vegagec_tree,vegagec_grass, &
           vegagec_pasture,vegagec_crop)

    ! retrieve the initial veget_mtc value, as explained above.
    veget_mtc = veget_mtc_begin

    ! we check whether the required harvest of secondary forest
    ! is met by the existing secondary forest fractions. Otherwise
    ! we use the oldest-age-class forest to compenstate it. 
    DO ipts=1,npts
      IF (FHmatrix_remainA(ipts,sf2yf) .GT. zero) THEN
        ! in this case, the existing secondary forest fraction
        ! is not enough for secondary forest harvest, we have to
        ! use primary (oldest age class) foret to compensate it.

        IndStart_f = 1             ! Oldest age class
        IndEnd_f = 1               ! Oldest age class

        !sf2yf
        CALL type_conversion(ipts,sf2yf,FHmatrix_remainA,veget_mtc,newvegfrac, &
                         indold_tree,indagec_tree,indagec_pasture,num_pasture_mulagec,&
                         IndEnd_f,nagec_herb,                    &
                         vegagec_tree,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                         glccRemain, &
                         .TRUE., iagec_start=IndStart_f)
        pf2yf_compen_sf2yf(ipts) = FHmatrix_remainA(ipts,sf2yf) - glccRemain(ipts,sf2yf)
        
        !!++Temp++
        ! Normally we don't allow youngest foret cohort to be harvested,
        ! but when external input driving data are not consistent, we harvest
        ! even the youngest forest cohort to make sure it's consistent with
        ! the input data and comparable with the case of SinAgeC. However, except
        ! for the thse two reasons, we should not harvest the youngest age class
        ! as the wood there is not feasible for usage in most cases.
        IF (allow_youngest_forest_WoodHarvest) THEN
          FHmatrix_remainA(ipts,sf2yf) = glccRemain(ipts,sf2yf)
          IF (FHmatrix_remainA(ipts,sf2yf) .GT. zero) THEN
    
            IndStart_f = nagec_tree    ! youngest age class
            IndEnd_f = nagec_tree      ! youngest age class

            IF (IndStart_f .LE. 0 .OR. IndEnd_f .LE. 0) THEN
              write(numout,*) 'Age class index cannot be negative or zero!'
              STOP
            ENDIF

            !sf2yf
            CALL type_conversion(ipts,sf2yf,FHmatrix_remainA,veget_mtc,newvegfrac, &
                             indold_tree,indagec_tree,indagec_pasture,num_pasture_mulagec,&
                             IndEnd_f,nagec_herb,                    &
                             vegagec_tree,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                             glccRemain, &
                             .TRUE., iagec_start=IndStart_f)
          ENDIF
        ENDIF !(IF allow_youngest_forest_WoodHarvest)

      ENDIF
    ENDDO
    FHmatrix_remainB(:,:) = glccRemain

    !! 1.3 Handle primary forest harvest
    
    CALL calc_cover(npts,veget_max,veget_mtc,vegagec_tree,vegagec_grass, &
           vegagec_pasture,vegagec_crop)
    veget_mtc = veget_mtc_begin

    ! We check first whether there is still deficit in the required secondary
    ! harvest even after being compensated for by the primary forests. 
    ! If yes, that means all existing harvestable forests
    ! are depleted, thus required primary harvest will be suppressed.
    ! Otherwise we will handle primary forest harvest starting from the
    ! oldest-age-class forest

    DO ipts=1,npts
      IF (FHmatrix_remainB(ipts,sf2yf) .GT. min_stomate) THEN
        ! in this case, all forest fraction is depleted in handling
        ! required secondary forest harvest. We thus suppress the 
        ! the required primary forest harvest.
        Deficit_sf2yf_final(ipts) = -1 * FHmatrix_remainB(ipts,sf2yf)
        Deficit_pf2yf_final(ipts) = -1 * FHmatrix_remainB(ipts,pf2yf)
        

      ELSE
        ! there are still forest can be used for required primary forest harvest.
        ! we assume primary harvest occurs in the oldest age class. Here,
        ! we will start from the oldest forest and then move to the younger ones,
        ! until the second youngest when the youngest forest cohort is not allowed
        ! for harvesting, and to youngest one when it is allowed.

        IndStart_f = 1             ! Oldest age class
        
        !!++Temp++
        IF (allow_youngest_forest_WoodHarvest) THEN
          IndEnd_f = nagec_tree      ! youngest age class
        ELSE
          IndEnd_f = nagec_tree-1    ! 2nd youngest age class
        ENDIF

        IF (IndStart_f .LE. 0 .OR. IndEnd_f .LE. 0) THEN
          write(numout,*) 'Age class index cannot be negative or zero!'
          STOP
        ENDIF

        !pf2yf
        CALL type_conversion(ipts,pf2yf,FHmatrix_remainB,veget_mtc,newvegfrac, &
                         indold_tree,indagec_tree,indagec_pasture,num_pasture_mulagec,&
                         IndEnd_f,nagec_herb,                    &
                         vegagec_tree,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                         glccRemain, &
                         .TRUE., iagec_start=IndStart_f)
      ENDIF
      
      IF (glccRemain(ipts,pf2yf) .GT. min_stomate) THEN
        Deficit_pf2yf_final(ipts) = -1 * glccRemain(ipts,pf2yf)
      ENDIF
    ENDDO
    glccRealHarvest = harvest_matrix - glccRemain

    ! Because we use the template of `type_conversion, now the glcc_pft_tmp
    ! and glcc_pftmtc have wrong information (because harvest loss is assigned
    ! on the newly created youngest-age-class pasture/crop MTCs). So they have 
    ! to be re-initialized to zero. Only the information in glcc_pft is what
    ! we need, as explained above.
    glcc_pft_tmp(:,:) = 0.
    glcc_pftmtc(:,:,:) = 0.
    ! Here we need to put glcc_pft into glcc_pftmtc for forestry harvest.
    ! The same MTC will be maintained when forest is harvested.
    DO ivm =1,nvm
      IF (is_tree(ivm)) THEN
        glcc_pftmtc(:,ivm,pft_to_mtc(ivm)) = glcc_pft(:,ivm)
      ENDIF 
    ENDDO
    glcc_pftmtc_harvest = glcc_pftmtc
    !****************** end block to handle forestry harvest ****************

    IF (.NOT. MulAgeC_fh_firstday_done) THEN

      glccReal_tmp = zero
      glccReal_tmp(:,1:12) = glccRealHarvest
      CALL histwrite_p (hist_id_stomate, 'glccRealHarvest', itime, &
           glccReal_tmp, npts*nvm, horipft_index)
      CALL xios_orchidee_send_field ('glccRealHarvest', glccReal_tmp) ! kjpindex,nvm

      DO ivma = 1, nvmap
        WRITE(part_str,'(I2)') ivma
        IF (ivma < 10) part_str(1:1) = '0'
        CALL histwrite_p (hist_id_stomate, 'glcc_pftmtc_H_'//part_str(1:LEN_TRIM(part_str)), &
             itime, glcc_pftmtc_harvest(:,:,ivma), npts*nvm, horipft_index)
      ENDDO
      CALL xios_orchidee_send_field ('glcc_pftmtc_H', glcc_pftmtc_harvest) ! kjpindex,nvm,nvmap
     
      MulAgeC_fh_firstday_done = .TRUE.
    ENDIF

  END SUBROUTINE MulAgeC_fh_firstday


  SUBROUTINE Get_harvest_matrix(npts,veget_max,newvegfrac,harvest_matrix, &
                          biomass, harvest_biomass, area_land_m2, &
                          glcc_pft,glcc_pftmtc, &
                          Deficit_pf2yf_final, Deficit_sf2yf_final,   &
                          pf2yf_compen_sf2yf, sf2yf_compen_pf2yf)

    IMPLICIT NONE

    !! 0.1 Input variables

    INTEGER, INTENT(in)                                     :: npts            !! Domain size - number of pixels (unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)            :: veget_max       !! 
    REAL(r_std), DIMENSION(npts), INTENT(in)                :: area_land_m2

    REAL(r_std), DIMENSION(npts,nvmap), INTENT(in)          :: newvegfrac
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(in)  :: biomass   !! Stand level biomass @tex $(gC.m^{-2})$ @endtex

    REAL(r_std), DIMENSION(npts,12),INTENT(in)              :: harvest_biomass  !! the 1st dimension is for industrial wood, 2nd dimension for fuelwood.

    !! 0.2 Output variables
    REAL(r_std), DIMENSION(npts,12),INTENT(inout)           :: harvest_matrix !! 


    !! 0.3 Modified variables
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)         :: glcc_pft       !! Loss of fraction in each PFT
    REAL(r_std), DIMENSION(npts,nvm,nvmap), INTENT(inout)   :: glcc_pftmtc    !! a temporary variable to hold the fractions each PFT is going to lose

    !! 0.4 Local variables
    REAL(r_std), DIMENSION(npts)             :: mean_abwood_dens        !! Stand level biomass @tex $(gC.m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts)             :: total_abwood            !! Stand level biomass @tex $(gC.m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts)             :: total_abwood_harvest    !! Stand level biomass @tex $(gC.m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts)             :: treefrac                !! "maximal" coverage fraction of a PFT on the ground

    INTEGER                :: ipts,ivm,ivma,num,ierr

    ! Note: these parameters are useless here but to conform with the form they
    ! are needed to call the subroutine MulAgeC_fh_firstday.
    REAL(r_std), DIMENSION(npts,nvmap)       :: glcc_mtc                !! Increase in fraction of each MTC in its youngest age-class
    REAL(r_std), DIMENSION(npts,nvm)         :: glccReal_tmp            !! A temporary variable to hold glccReal
    REAL(r_std), DIMENSION(npts)             :: Deficit_pf2yf_final     !! 
    REAL(r_std), DIMENSION(npts)             :: Deficit_sf2yf_final     !! 
    REAL(r_std), DIMENSION(npts)             :: pf2yf_compen_sf2yf      !! 
    REAL(r_std), DIMENSION(npts)             :: sf2yf_compen_pf2yf      !!
    REAL(r_std)                              :: ratio,input_total,model_total
    REAL(r_std)                              :: input_total_allproc,model_total_allproc


    input_total = SUM(harvest_biomass(:,1:2))

    treefrac(:) = zero
    total_abwood(:) = zero
    DO ivm = 1,nvm
      IF (is_tree(ivm) ) THEN
        total_abwood(:) = total_abwood(:) + veget_max(:,ivm) * &
          (biomass(:,ivm,isapabove,icarbon)+biomass(:,ivm,iheartabove,icarbon))
        treefrac(:) = treefrac(:) + veget_max(:,ivm)
      ENDIF
    ENDDO

    mean_abwood_dens(:) = zero
    WHERE( treefrac(:) .GT. min_stomate )
      mean_abwood_dens(:) = total_abwood(:)/treefrac(:)
    ENDWHERE

    ! harvest_matrix(:,:) = zero
    ! WHERE( mean_abwood_dens(:) .GT. min_stomate .AND. area_land_m2(:) .GT. min_stomate )
    !   harvest_matrix(:,1) = SUM(harvest_biomass(:,1:2),DIM=2)/(mean_abwood_dens(:)*area_land_m2(:))
    ! ENDWHERE

    CALL MulAgeC_fh_firstday(npts,veget_max,newvegfrac,harvest_matrix,   &
                          glcc_pft,glcc_pftmtc,  &
                          Deficit_pf2yf_final, Deficit_sf2yf_final,   &
                          pf2yf_compen_sf2yf, sf2yf_compen_pf2yf)

    CALL ratio_io_wood_harvest(npts, glcc_pft, area_land_m2, biomass, harvest_biomass, & ! In
                               ratio ) ! Out

    ! We adjust the harvest_matrix if the output total harvet wood
    ! geos beyond 10% difference with the input harvest biomass.
    num=1
    DO WHILE ( (ABS(ratio-1) .GT. 0.01) .AND. (ratio.GT.min_stomate) .AND. num .LT. 20)

      ! Note here we treat all harvest as primary harvest
      harvest_matrix(:,1) = harvest_matrix(:,1)/ratio
      CALL MulAgeC_fh_firstday(npts,veget_max,newvegfrac,harvest_matrix,   &
                            glcc_pft,glcc_pftmtc, &
                            Deficit_pf2yf_final, Deficit_sf2yf_final,   &
                            pf2yf_compen_sf2yf, sf2yf_compen_pf2yf)

      CALL ratio_io_wood_harvest(npts, glcc_pft, area_land_m2, biomass, harvest_biomass, & ! In
                               ratio ) ! Out

      num=num+1
      WRITE(numout,*) 'Get_harvest_matrix:: Wood harvest matrix iteration: ',num

    ENDDO

  END SUBROUTINE Get_harvest_matrix


! ================================================================================================================================
!! SUBROUTINE   ratio_io_wood_harvest
!!
!>\BRIEF       : Calculate ratio for input/output wood harvest 
!!
!>\DESCRIPTION  
!_ ================================================================================================================================
  SUBROUTINE ratio_io_wood_harvest(npts,    glcc_pft,   area_land_m2,   biomass, harvest_biomass, & ! In
                                   ratio ) ! Out

    !! 0.1 Input variables
    INTEGER, INTENT(in)                                             :: npts             !! Domain size - number of pixels (unitless)
    REAL(r_std), DIMENSION(npts), INTENT(in)                        :: area_land_m2
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(in)   :: biomass          !! Stand level biomass @tex $(gC.m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)                    :: glcc_pft         !! Loss of fraction in each PFT
    REAL(r_std), DIMENSION(npts,12), INTENT(in)                     :: harvest_biomass  !! the 1st dimension is for industrial wood, 2nd dimension for fuelwood.

    !! 0.2 Output variables
    REAL(r_std), INTENT(out)                                        :: ratio            !! Input/Output wood harvest ratio 

    !! 0.4 Local variables
    REAL(r_std), DIMENSION(npts)             :: total_abwood_harvest    !! Stand level biomass @tex $(gC.m^{-2})$ @endtex
    REAL(r_std)                              :: input_total,model_total
    REAL(r_std)                              :: input_total_allproc,model_total_allproc
    INTEGER(i_std)                           :: ivm


    input_total = SUM(harvest_biomass(:,1:2))

    ! calculate model_total
    total_abwood_harvest(:) = zero
    DO ivm = 1,nvm
      IF (is_tree(ivm) ) THEN
        total_abwood_harvest(:) = total_abwood_harvest(:) + glcc_pft(:,ivm) *area_land_m2(:) * &
          (biomass(:,ivm,isapabove,icarbon)+biomass(:,ivm,iheartabove,icarbon))
      ENDIF
    ENDDO

    model_total = SUM( total_abwood_harvest ) ! Change to
    WRITE(numout,*) 'ratio_io_wood_harvest:: model_total: ',model_total

    CALL allreduce_sum(input_total, input_total_allproc)
    WRITE(numout,*) 'ratio_io_wood_harvest:: Regional total input wood harvest,', input_total_allproc
    IF ( input_total_allproc .GT. min_stomate) THEN
      CALL allreduce_sum(model_total, model_total_allproc)
      ratio = model_total_allproc/input_total_allproc
      WRITE(numout,*) 'ratio_io_wood_harvest:: Initial ratio is: ',ratio
      WRITE(numout,*) 'ratio_io_wood_harvest:: Initial regional total output wood harvest,', model_total_allproc
    ENDIF

  END SUBROUTINE ratio_io_wood_harvest

END MODULE stomate_fharvest_MulAgeC
