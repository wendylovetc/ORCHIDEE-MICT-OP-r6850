! Remove definition of type_conversion

! =================================================================================================================================
! MODULE       : stomate_glcc_bioe1
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       This module implements gross land use change with age classes.
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


MODULE stomate_glcc_bioe1

  ! modules used:
  
  USE ioipsl_para
  USE stomate_data
  USE pft_parameters
  USE constantes
  USE constantes_soil_var
  USE stomate_gluc_common
  USE stomate_gluc_constants
  USE stomate_glcchange_MulAgeC
  USE xios_orchidee

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC glcc_bioe1_firstday, glcc_bioe1 
  
CONTAINS

! ================================================================================================================================
!! SUBROUTINE   gross_lcchange
!!
!>\BRIEF       : Apply gross land cover change.
!!
!>\DESCRIPTION  
!_ ================================================================================================================================
  SUBROUTINE glcc_bioe1 (npts, dt_days, newvegfrac,  &
               glccSecondShift,&
               def_fuel_1hr_remain, def_fuel_10hr_remain,        &
               def_fuel_100hr_remain, def_fuel_1000hr_remain,    &
               deforest_litter_remain, deforest_biomass_remain,  &
               convflux, cflux_prod10, cflux_prod100,                  &
               glcc_pft, glcc_pftmtc,          &
               veget_max, prod10, prod100, flux10, flux100,            &
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
    REAL(r_std), DIMENSION (npts,12),INTENT(in)        :: glccSecondShift     !! the land-cover-change (LCC) matrix in case a gross LCC is 
                                                                              !! used.
    REAL(r_std), DIMENSION (npts,nvmap),INTENT(in)       :: newvegfrac             !! 
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
    REAL(r_std), DIMENSION(npts,nwp), INTENT(out)            :: convflux         !! release during first year following land cover
                                                                             !! change
    REAL(r_std), DIMENSION(npts,nwp), INTENT(out)            :: cflux_prod10     !! total annual release from the 10 year-turnover
                                                                             !! pool @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts,nwp), INTENT(out)            :: cflux_prod100    !! total annual release from the 100 year-
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
    REAL(r_std), DIMENSION(npts,10,nwp), INTENT(inout)       :: flux10           !! annual release from the 10/100 year-turnover 
                                                                             !! pool compartments
    REAL(r_std), DIMENSION(npts,100,nwp), INTENT(inout)      :: flux100          !! annual release from the 10/100 year-turnover
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
    REAL(r_std), DIMENSION(npts,nvm)         :: glccReal_tmp         !! A temporary variable

    WRITE(numout,*) 'Entering glcc_MulAgeC'
    glccReal_tmp(:,:) = zero

    CALL glcc_bioe1_firstday(npts,veget_max,newvegfrac,   &
                          glccSecondShift,glcc_pft,glcc_pftmtc)

    glcc_mtc(:,:) = SUM(glcc_pftmtc,DIM=2)

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

        ! here we set glcc_mtc(ipts,ivma) > min_stomate as a condition,
        ! this is necessary because later on in the subroutine of 
        ! add_incoming_proxy_pft we have to merge the newly established
        ! youngest proxy with potentially exisiting youngest receiving MTC,
        ! thus have to devide a new fraction of (frac_proxy + frac_exist),
        ! but in case frac_exist = zero, we risk deviding by a very small value
        ! of frac_proxy and thus we want it to be bigger than min_stomate.
        IF ( glcc_mtc(ipts,ivma) .GT. min_stomate ) THEN

          ! 1. we accumulate the scalar variables that will be inherited.
          !    note that in the subroutine of `collect_legacy_pft`, all
          !    zero transitions will be automatically skipped. 
          CALL collect_legacy_pft(npts, ipts, ivma, glcc_pftmtc,    &
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
      ENDDO ! (DO ivma = 2,nvmap)

      ! We can only establish new youngest proxy and add it to the 
      ! existing youngest-age PFT after all wood product extraction is done, to 
      ! avoid the dilution of extractable biomass by the young proxy
      ! and ensure consistency. Therefore now we have to loop again
      ! over nvmap.
      DO ivma = 2,nvmap
        !IF ( glcc_mtc(ipts,ivma) .GT. min_stomate ) THEN

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
          
        !ENDIF !IF ( glcc_mtc(ipts,ivma) .GT. min_stomate )
      ENDDO !(DO ivma=1,nvmap)

    ENDDO !(DO ipts=1,npts)

    !! Update 10 year-turnover pool content following flux emission
    !!     (linear decay (10%) of the initial carbon input)
    DO  l = 0, 8
      m = 10 - l
      cflux_prod10(:,:) =  cflux_prod10(:,:) + flux10(:,m,:)
      prod10(:,m,:)     =  prod10(:,m-1,:)   - flux10(:,m-1,:)
      flux10(:,m,:)     =  flux10(:,m-1,:)
    ENDDO
    
    cflux_prod10(:,:) = cflux_prod10(:,:) + flux10(:,1,:) 
    flux10(:,1,:)     = 0.1 * prod10(:,0,:)
    prod10(:,1,:)     = prod10(:,0,:)
    
    !! 2.4.3 update 100 year-turnover pool content following flux emission\n
    DO l = 0, 98
       m = 100 - l
       cflux_prod100(:,:)  =  cflux_prod100(:,:) + flux100(:,m,:)
       prod100(:,m,:)      =  prod100(:,m-1,:)   - flux100(:,m-1,:)
       flux100(:,m,:)      =  flux100(:,m-1,:)
    ENDDO
    
    cflux_prod100(:,:)  = cflux_prod100(:,:) + flux100(:,1,:) 
    flux100(:,1,:)      = 0.01 * prod100(:,0,:)
    prod100(:,1,:)      = prod100(:,0,:)
    prod10(:,0,:)        = zero
    prod100(:,0,:)       = zero 

    ! Write out history files
    CALL histwrite_p (hist_id_stomate, 'glcc_pft', itime, &
         glcc_pft, npts*nvm, horipft_index)
    CALL xios_orchidee_send_field ('glcc_pft', glcc_pft)

    DO ivma = 1, nvmap
      WRITE(part_str,'(I2)') ivma
      IF (ivma < 10) part_str(1:1) = '0'
      CALL histwrite_p (hist_id_stomate, 'glcc_pftmtc_'//part_str(1:LEN_TRIM(part_str)), &
           itime, glcc_pftmtc(:,:,ivma), npts*nvm, horipft_index)
    ENDDO
    CALL xios_orchidee_send_field ('glcc_pftmtc', glcc_pftmtc) ! kjpindex,nvm,nvmap

  END SUBROUTINE glcc_bioe1


! ================================================================================================================================
!! SUBROUTINE   : glcc_bioe1_firstday
!!
!>\BRIEF        : When necessary, adjust input glcc matrix, and allocate it
!!                into different contributing age classes and receiving 
!!                youngest age classes.
!! \n
!_ ================================================================================================================================

  ! Note: it has this name because this subroutine will also be called
  ! the first day of each year to precalculate the forest loss for the
  ! deforestation fire module.
  SUBROUTINE glcc_bioe1_firstday(npts,veget_max_org,newvegfrac, &
                          glccSecondShift, glcc_pft, glcc_pftmtc)

    IMPLICIT NONE

    !! 0.1 Input variables

    INTEGER, INTENT(in)                                     :: npts           !! Domain size - number of pixels (unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)            :: veget_max_org  !! "maximal" coverage fraction of a PFT on the ground
                                                                              !! May sum to
                                                                              !! less than unity if the pixel has
                                                                              !! nobio area. (unitless, 0-1)
    REAL(r_std), DIMENSION(npts,nvmap), INTENT(in)          :: newvegfrac     !! used to guid the allocation of new PFTs.
                                                                              !! 
    REAL(r_std), DIMENSION (npts,12),INTENT(in)        :: glccSecondShift     !! the land-cover-change (LCC) matrix in case a gross LCC is 
                                                                              !! used.

    !! 0.2 Output variables
    REAL(r_std), DIMENSION(npts,nvm,nvmap), INTENT(inout)   :: glcc_pftmtc    !! a temporary variable to hold the fractions each PFT is going to lose
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)         :: glcc_pft       !! Loss of fraction in each PFT

    !! 0.3 Modified variables
    
    !! 0.4 Local variables
    REAL(r_std), DIMENSION(npts,nvmap)              :: veget_mtc           !! "maximal" coverage fraction of a PFT on the ground
    REAL(r_std), DIMENSION(npts,nvmap)              :: veget_mtc_begin     !! "maximal" coverage fraction of a PFT on the ground
    REAL(r_std), DIMENSION(npts,nagec_tree)         :: vegagec_tree        !! fraction of tree age-class groups, in sequence of old->young
    REAL(r_std), DIMENSION(npts,nagec_herb)         :: vegagec_grass       !! fraction of grass age-class groups, in sequence of old->young 
    REAL(r_std), DIMENSION(npts,nagec_herb)         :: vegagec_pasture     !! fraction of pasture age-class groups, in sequence of old->young
    REAL(r_std), DIMENSION(npts,nagec_herb)         :: vegagec_crop        !! fraction of crop age-class groups, in sequence of old->young
    REAL(r_std), DIMENSION(npts,nagec_bioe1)        :: vegagec_bioe1       !! fraction of crop age-class groups, in sequence of old->young

    
    REAL(r_std), DIMENSION(npts,4)                  :: veget_4veg      !! "maximal" coverage fraction of a PFT on the ground
    REAL(r_std), DIMENSION(npts)                    :: veget_tree      !! "maximal" coverage fraction of a PFT on the ground
    REAL(r_std), DIMENSION(npts)                    :: veget_grass     !! "maximal" coverage fraction of a PFT on the ground
    REAL(r_std), DIMENSION(npts)                    :: veget_pasture   !! "maximal" coverage fraction of a PFT on the ground
    REAL(r_std), DIMENSION(npts)                    :: veget_crop      !! "maximal" coverage fraction of a PFT on the ground
    REAL(r_std), DIMENSION(npts)                    :: veget_bioe1     !! "maximal" coverage fraction of a PFT on the ground

    REAL(r_std), DIMENSION(npts,nvm)         :: veget_max              !! "maximal" coverage fraction of a PFT on the ground
    REAL(r_std), DIMENSION(npts,nvm)         :: veget_max_old          !! "maximal" coverage fraction of a PFT on the ground
    REAL(r_std), DIMENSION(npts,nvm)         :: glcc_pft_tmp           !! Loss of fraction in each PFT

    REAL(r_std), DIMENSION(npts,nvm,nvmap)  :: glcc_pftmtc_SecShift   !! a temporary variable to hold the fractions each PFT is going to lose
    REAL(r_std), DIMENSION(npts,12)      :: glccRealSecShift   !! real matrix applied for secondary shifting cultivation.
    REAL(r_std), DIMENSION(npts,12)      :: glccDefSecShift    !! deficit for the glccSecondShift
    REAL(r_std), DIMENSION(npts,12)         :: glccRemain      !! 
    REAL(r_std), DIMENSION(npts,12)         :: glccSecondShift_remain      !! 

    REAL(r_std), DIMENSION(npts,nvm)         :: glccReal_tmp         !! A temporary variable 

    LOGICAL, SAVE  :: glcc_bioe1_firstday_done = .FALSE.

    ! Different indexes for convenient local uses
    INTEGER :: f_to_bioe1=1, g_to_bioe1=2, p_to_bioe1=3, c_to_bioe1=4, & 
               bioe1_to_f=5, bioe1_to_g=6, bioe1_to_p=7, bioe1_to_c=8

    INTEGER :: ivma

    INTEGER :: ipts,IndStart_f,IndEnd_f
    CHARACTER(LEN=10)      :: part_str                               !! string suffix indicating an index

    !Some more local configurations
    LOGICAL                                 :: allow_youngest_forest_convert = .TRUE.
    
    ! Initialization
    glcc_pftmtc = zero
    glcc_pft = zero
    glcc_pft_tmp = zero

    !!! ** Land cover change processes start here ** !!!
    ! we make copies of original input veget_max (which is veget_max_org
    ! in the subroutine parameter list).
    ! veget_max will be modified through different operations in order to 
    ! check various purposes, e.g., whether input glcc matrix 
    ! is compatible with existing veget_max and how to allocate it etc.
    ! veget_max_old will not be modified
    veget_max(:,:) = veget_max_org(:,:)
    veget_max_old(:,:) = veget_max_org(:,:)

    !! 3. Treat secondary-agriculture shifting cultivation transition matrix
    !! [The primary-agriculture shifting cultivation will be treated together
    !!  with the netLCC transitions, with the conversion sequence of oldest->
    !!  youngest is applied.]
    ! When we prepare the driving data, secondary-agriculture shifting cultivation
    ! is intended to include the "constant transitions" over time. Ideally, we
    ! should start applying this secondary-agriculture shifting cultivation with
    ! the "secondary forest" in the model. Here we tentatively start with the 3rd
    ! youngest age class and move to the 2ne youngest age class. But if the prescribed
    ! transition fraction is not met, we then move further to 4th youngest age class
    ! and then move to the oldest age class sequentially.

    CALL calc_cover_bioe1(npts,veget_max,veget_mtc_begin,vegagec_tree,vegagec_grass, &
           vegagec_pasture,vegagec_crop,vegagec_bioe1)
    veget_mtc = veget_mtc_begin
 
    !! 3.1 We start treating secondary-agriculture cultivation from the 3rd youngest 
    !! age class and then move to the younger age class. 
    ! Because it's rather complicated to calculate which transtion fraction between
    ! which vegetation types should occur in case there is deficit occuring
    ! for the overall donation vegetation type, we will just start from some 
    ! priority and leave the unrealized parts into the latter section.

    ! For this purpose, we should first make a copy of glccSecondShift into 
    ! glccRemain. glccRemain will tell us the transition fractions that have to
    ! be treated starting from `IndStart_f+1` oldest age class and moving torward older
    ! age class.
    glccRemain(:,:) = glccSecondShift(:,:)

    ! Now we will call type_conversion for each of the 12 transitions, starting
    ! from `IndStart_f` age class moving to the 2nd youngest age class. We use glccRemain
    ! to track the transtion fractions we should leave for the second case. 
    ! To make the code more flexible, we will store the start and end indecies
    ! in variables.

    !*[Note: we do above process only for forest now, as we assume the conversion
    !  of crop/pasture/grass to other types will start always from the oldest
    !  age class]

    IndStart_f = nagec_tree-1  ! note the indecies and vegetfrac for tree age class
                               ! is from old to young
    IndEnd_f = nagec_tree    ! nagec_tree-2: The 3rd youngest age class
                               ! nagec_tree-1: The 2nd youngest age class
                               ! nagec_tree: The youngest age class

    IF (IndStart_f .LE. 0 .OR. IndEnd_f .LE. 0) THEN
      write(numout,*) 'glcc_bioe1_firstday: Age class index cannot be negative or zero!'
      STOP
    ENDIF

    DO ipts=1,npts
      !f_to_bioe1
      CALL type_conversion(ipts,f_to_bioe1,glccSecondShift,veget_mtc,newvegfrac,       &
                       indold_tree,indagec_tree,indagec_bioe1,num_bioe1_mulagec,     &
                       IndEnd_f,nagec_bioe1,                    &
                       vegagec_tree,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp, &
                       glccRemain, &
                       .TRUE., iagec_start=IndStart_f)
      !g_to_bioe1
      CALL type_conversion(ipts,g_to_bioe1,glccSecondShift,veget_mtc,newvegfrac,       &
                       indold_grass,indagec_grass,indagec_bioe1,num_bioe1_mulagec,     &
                       nagec_herb,nagec_bioe1,                    &
                       vegagec_grass,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain, &
                       .TRUE.)
      !p_to_bioe1
      CALL type_conversion(ipts,p_to_bioe1,glccSecondShift,veget_mtc,newvegfrac,       &
                       indold_pasture,indagec_pasture,indagec_bioe1,num_bioe1_mulagec,     &
                       nagec_herb,nagec_bioe1,                    &
                       vegagec_pasture,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain, &
                       .TRUE.)
      !c_to_bioe1
      CALL type_conversion(ipts,c_to_bioe1,glccSecondShift,veget_mtc,newvegfrac,       &
                       indold_crop,indagec_crop,indagec_bioe1,num_bioe1_mulagec,     &
                       nagec_herb,nagec_bioe1,                    &
                       vegagec_crop,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain, &
                       .TRUE.)
      !bioe1_to_f
      CALL type_conversion(ipts,bioe1_to_f,glccSecondShift,veget_mtc,newvegfrac,       &
                       indold_bioe1,indagec_bioe1,indagec_tree,num_tree_mulagec,     &
                       nagec_bioe1,nagec_tree,                    &
                       vegagec_bioe1,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain, &
                       .TRUE.)
      !bioe1_to_g
      CALL type_conversion(ipts,bioe1_to_g,glccSecondShift,veget_mtc,newvegfrac,       &
                       indold_bioe1,indagec_bioe1,indagec_grass,num_grass_mulagec,     &
                       nagec_bioe1,nagec_herb,                    &
                       vegagec_bioe1,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain, &
                       .TRUE.)
      !bioe1_to_p
      CALL type_conversion(ipts,bioe1_to_p,glccSecondShift,veget_mtc,newvegfrac,       &
                       indold_bioe1,indagec_bioe1,indagec_pasture,num_pasture_mulagec,     &
                       nagec_bioe1,nagec_herb,                    &
                       vegagec_bioe1,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain, &
                       .TRUE.)
      !bioe1_to_c
      CALL type_conversion(ipts,bioe1_to_c,glccSecondShift,veget_mtc,newvegfrac,       &
                       indold_bioe1,indagec_bioe1,indagec_crop,num_crop_mulagec,     &
                       nagec_bioe1,nagec_herb,                    &
                       vegagec_bioe1,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain, &
                       .TRUE.)
    ENDDO
    glccSecondShift_remain(:,:) = glccRemain(:,:)

    !! 3.2 We treat the remaing unrealized transtions from forest. Now we will
    !! start with the 3rd oldest age class and then move to the oldest age class. 

    CALL calc_cover_bioe1(npts,veget_max,veget_mtc,vegagec_tree,vegagec_grass, &
           vegagec_pasture,vegagec_crop,vegagec_bioe1)
    veget_mtc = veget_mtc_begin
 
    IndStart_f = nagec_tree  ! note the indecies and vegetfrac for tree age class
                               ! is from old to young. 
                               ! nagec_tree -3: The 4th youngest age class.

    IndEnd_f = 1               ! oldest-age class forest.

    IF (IndStart_f .LE. 0 .OR. IndEnd_f .LE. 0) THEN
      write(numout,*) 'glcc_bioe1_firstday: Age class index cannot be negative or zero!'
      STOP
    ENDIF

    ! we start with the 3rd youngest age class and move up to the oldest age
    ! class in the sequence of young->old, as indicated by the .FALSE. parameter
    ! when calling the subroutine type_conversion.
    DO ipts=1,npts
      !f_to_bioe1
      CALL type_conversion(ipts,f_to_bioe1,glccSecondShift_remain,veget_mtc,newvegfrac,       &
                       indold_tree,indagec_tree,indagec_bioe1,num_bioe1_mulagec,     &
                       IndEnd_f,nagec_bioe1,                    &
                       vegagec_tree,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp, &
                       glccRemain, &
                       .FALSE., iagec_start=IndStart_f)
    ENDDO

    IF (allow_youngest_forest_convert) THEN
      !!++Temp++!!
      !! this block of 3.3 could be commented to remove this process as desribed
      !! below.

      ! [2016-04-20] This is temporarily added: Normally we assume the youngest
      ! forest age cohort will not be cut because in a shifting cultivation, they
      ! are grown to let the land recover from agricultural process. (Or at least)
      ! we can set the threshold of youngest age cohort to be very small. But there
      ! are two reasons we allow the youngest forest cohort to be cut for shifting
      ! cultivation purpose: a). Farmers may decide to harvest the wood of a forest
      ! and then convert to crop. We don't simulate explicitly this process because 
      ! this will depend on input land change matrix and land use data assumptions.
      ! However,we can implicitly account for this by assuming "farmers plant young
      ! trees after harvesting the wood, and immediately convert this young trees
      ! to crops. b). For the sake of conserving the total sum of veget_max before
      ! and after the transition as one, we need to allow the youngest forest cohort
      ! eligible for cutting.

      !! 3.3 We treat the remaing unrealized transtions from forest, allowing 
      !! the youngest forest cohort to be cut. For this purpose, we will
      !! start with the 2nd youngest age class and then move to the youngest one. 

      glccSecondShift_remain(:,:) = glccRemain(:,:)

      CALL calc_cover_bioe1(npts,veget_max,veget_mtc,vegagec_tree,vegagec_grass, &
           vegagec_pasture,vegagec_crop,vegagec_bioe1)
      veget_mtc = veget_mtc_begin
 
      ! Note: the setting of index here must be consistent with those of 3.1 and 3.2
      IndStart_f = nagec_tree-1  ! note the indecies and vegetfrac for tree age class
                                 ! is from old to young. 
                                 ! nagec_tree -1: The 2nd youngest age class.

      IndEnd_f = nagec_tree      ! youngest class forest.

      IF (IndStart_f .LE. 0 .OR. IndEnd_f .LE. 0) THEN
        write(numout,*) 'glcc_bioe1_firstday: Age class index cannot be negative or zero!'
        STOP
      ENDIF

      ! we start with the 3rd youngest age class and move up to the oldest age
      ! class in the sequence of young->old, as indicated by the .FALSE. parameter
      ! when calling the subroutine type_conversion.
      DO ipts=1,npts
        !f_to_bioe1
        CALL type_conversion(ipts,f_to_bioe1,glccSecondShift_remain,veget_mtc,newvegfrac,       &
                         indold_tree,indagec_tree,indagec_bioe1,num_bioe1_mulagec,     &
                         IndEnd_f,nagec_bioe1,                    &
                         vegagec_tree,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp, &
                         glccRemain, &
                         .FALSE., iagec_start=IndStart_f)
      ENDDO
      !! End of ++Temp++ Section 3.3
    ENDIF

    ! Final handling of some output variables.
    ! we separate the glcc_pftmtc_SecShift
    glcc_pftmtc_SecShift = glcc_pftmtc

    ! we put the remaining glccRemain into the deficit
    glccDefSecShift = -1 * glccRemain
    glccRealSecShift = glccSecondShift - glccRemain

    !*****end block to handle glcc involving bioenergy vegtation type *******

    IF (.NOT. glcc_bioe1_firstday_done) THEN

      glccReal_tmp = zero

      glccReal_tmp(:,1:12) = glccRealSecShift
      CALL histwrite_p (hist_id_stomate, 'glccRealSecShift', itime, &
           glccReal_tmp, npts*nvm, horipft_index)
      CALL xios_orchidee_send_field ('glccRealSecShift', glccReal_tmp) ! kjpindex,nvm

      glccReal_tmp(:,1:12) = glccDefSecShift
      CALL histwrite_p (hist_id_stomate, 'glccDefSecShift', itime, &
           glccReal_tmp, npts*nvm, horipft_index)
      CALL xios_orchidee_send_field ('glccDefSecShift', glccReal_tmp) ! kjpindex,nvm

      DO ivma = 1, nvmap
        WRITE(part_str,'(I2)') ivma
        IF (ivma < 10) part_str(1:1) = '0'
        CALL histwrite_p (hist_id_stomate, 'glcc_pftmtc_SF_'//part_str(1:LEN_TRIM(part_str)), &
             itime, glcc_pftmtc_SecShift(:,:,ivma), npts*nvm, horipft_index)
      ENDDO
      CALL xios_orchidee_send_field ('glcc_pftmtc_SecShift', glcc_pftmtc_SecShift) ! kjpindex,nvm,nvmap

      glcc_bioe1_firstday_done = .TRUE.
    ENDIF

  END SUBROUTINE glcc_bioe1_firstday

  SUBROUTINE calc_cover_bioe1(npts,veget_max,veget_mtc,vegagec_tree,vegagec_grass, &
                 vegagec_pasture,vegagec_crop,vegagec_bioe1)

   
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
    REAL(r_std), DIMENSION(npts,nagec_bioe1), INTENT(inout)    :: vegagec_bioe1     !! fraction of bioenergy tree age-class groups, in sequence of old->young

    !! Local variables
    INTEGER(i_std)                                          :: ivma,staind,endind,j    !! indices (unitless)

    veget_mtc(:,:) = 0.
    vegagec_tree(:,:) = 0.
    vegagec_grass(:,:) = 0.
    vegagec_pasture(:,:) = 0.
    vegagec_crop(:,:) = 0.
    vegagec_bioe1(:,:) = 0.

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
        WRITE(numout,*) "Error: metaclass has only a single age group: ",ivma
        STOP
      ELSE
        IF (is_tree(staind)) THEN
          DO j=1,nagec_tree
            vegagec_tree(:,j) = vegagec_tree(:,j)+veget_max(:,endind-j+1)
          ENDDO
        ELSE IF (is_bioe1(staind)) THEN
          DO j=1,nagec_bioe1
            vegagec_bioe1(:,j) = vegagec_bioe1(:,j)+veget_max(:,endind-j+1)
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

  END SUBROUTINE calc_cover_bioe1

END MODULE stomate_glcc_bioe1
