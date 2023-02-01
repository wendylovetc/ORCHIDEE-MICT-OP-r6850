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


MODULE stomate_glcchange_SinAgeC

  ! modules used:
  
  USE ioipsl_para
  USE stomate_data
  USE pft_parameters
  USE constantes
  USE constantes_soil_var
  USE stomate_gluc_common
  USE xios_orchidee
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC glcc_SinAgeC_firstday, glcc_SinAgeC, type_conversion
  
CONTAINS

! ================================================================================================================================
!! SUBROUTINE   gross_lcchange
!!
!>\BRIEF       : Apply gross land cover change.
!!
!>\DESCRIPTION  
!_ ================================================================================================================================
  SUBROUTINE glcc_SinAgeC (npts, dt_days, harvest_matrix,newvegfrac,   &
               glccSecondShift,glccPrimaryShift,glccNetLCC,&
               def_fuel_1hr_remain, def_fuel_10hr_remain,        &
               def_fuel_100hr_remain, def_fuel_1000hr_remain,    &
               deforest_litter_remain, deforest_biomass_remain,  &
               convflux, cflux_prod10, cflux_prod100,                  &
               glccReal, IncreDeficit, glcc_pft, glcc_pftmtc,          &
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
    REAL(r_std), DIMENSION (npts,12),INTENT(in)        :: glccPrimaryShift    !! the land-cover-change (LCC) matrix in case a gross LCC is 
                                                                              !! used.
    REAL(r_std), DIMENSION (npts,12),INTENT(in)        :: glccNetLCC          !! the land-cover-change (LCC) matrix in case a gross LCC is 
                                                                              !! used.
    REAL(r_std), DIMENSION (npts,12),INTENT(in)          :: harvest_matrix             !! 
                                                                             !! 

    REAL(r_std), DIMENSION (npts,nvmap),INTENT(in)        :: newvegfrac             !! 
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
    REAL(r_std), DIMENSION(npts,12), INTENT(inout)       :: glccReal         !! The "real" glcc matrix that we apply in the model
                                                                             !! after considering the consistency between presribed
                                                                             !! glcc matrix and existing vegetation fractions.
    REAL(r_std), DIMENSION(npts,12), INTENT(inout)        :: IncreDeficit     !! "Increment" deficits, negative values mean that 
                                                                             !! there are not enough fractions in the source PFTs
                                                                             !! /vegetations to target PFTs/vegetations. I.e., these
                                                                             !! fraction transfers are presribed in LCC matrix but
                                                                             !! not realized.
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
    REAL(r_std), DIMENSION(npts,nvm)         :: glccReal_tmp         !! A temporary variable to hold glccReal
    REAL(r_std), DIMENSION(npts)             :: Deficit_pf2yf_final     !! 
    REAL(r_std), DIMENSION(npts)             :: Deficit_sf2yf_final     !! 
    REAL(r_std), DIMENSION(npts)             :: pf2yf_compen_sf2yf      !! 
    REAL(r_std), DIMENSION(npts)             :: sf2yf_compen_pf2yf      !!
    REAL(r_std), DIMENSION(npts,nvm)         :: glcc_harvest            !! Loss of fraction due to forestry harvest

    WRITE(numout,*) 'Entering glcc_SinAgeC'
    glcc_harvest(:,:) = zero
    glccReal_tmp(:,:) = zero

    CALL glcc_SinAgeC_firstday(npts,veget_max,newvegfrac,harvest_matrix,  &
                          glccSecondShift,glccPrimaryShift,glccNetLCC,&
                          glccReal,glcc_pft,glcc_pftmtc,IncreDeficit,  &
                          Deficit_pf2yf_final, Deficit_sf2yf_final,   &
                          pf2yf_compen_sf2yf, sf2yf_compen_pf2yf)

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

    !convflux        = convflux/one_year*dt_days
    !cflux_prod10    = cflux_prod10/one_year*dt_days
    !cflux_prod100   = cflux_prod100/one_year*dt_days

    ! Write out history files
    CALL histwrite_p (hist_id_stomate, 'glcc_pft', itime, &
         glcc_pft, npts*nvm, horipft_index)
    CALL xios_orchidee_send_field ('glcc_pft', glcc_pft) ! kjpindex,nvm

    glccReal_tmp(:,1:12) = glccReal
    CALL histwrite_p (hist_id_stomate, 'glccReal', itime, &
         glccReal_tmp, npts*nvm, horipft_index)
    CALL xios_orchidee_send_field ('glccReal', glccReal_tmp) ! kjpindex,nvm

    ! ! Write out forestry harvest variables
    ! DO ipts = 1,npts
    !   DO ivm = 1,nvm
    !     DO ivma = 1,nvmap
    !       IF (is_tree(ivm) .AND. is_tree(start_index(ivma))) THEN
    !         glcc_harvest(ipts,ivm) = glcc_harvest(ipts,ivm) + glcc_pftmtc(ipts,ivm,ivma)
    !       ENDIF
    !     ENDDO
    !   ENDDO
    ! ENDDO
    ! CALL histwrite_p (hist_id_stomate, 'glcc_harvest', itime, &
    !      glcc_harvest, npts*nvm, horipft_index)

    glccReal_tmp(:,:) = zero
    glccReal_tmp(:,1:12) = IncreDeficit
    CALL histwrite_p (hist_id_stomate, 'IncreDeficit', itime, &
         glccReal_tmp, npts*nvm, horipft_index)
    CALL xios_orchidee_send_field ('IncreDeficit', glccReal_tmp) ! kjpindex,nvm

    ! glccReal_tmp(:,:) = zero
    ! glccReal_tmp(:,1) = Deficit_pf2yf_final
    ! glccReal_tmp(:,2) = Deficit_sf2yf_final    ! is always zero in case of
    !                                            ! single age class
    ! glccReal_tmp(:,3) = pf2yf_compen_sf2yf     ! alawys zero for SinAgeC
    ! glccReal_tmp(:,4) = sf2yf_compen_pf2yf     ! always zero for SinAgeC

    ! CALL histwrite_p (hist_id_stomate, 'DefiComForHarvest', itime, &
    !      glccReal_tmp, npts*nvm, horipft_index)

    DO ivma = 1, nvmap
      WRITE(part_str,'(I2)') ivma
      IF (ivma < 10) part_str(1:1) = '0'
      CALL histwrite_p (hist_id_stomate, 'glcc_pftmtc_'//part_str(1:LEN_TRIM(part_str)), &
           itime, glcc_pftmtc(:,:,ivma), npts*nvm, horipft_index)
    ENDDO
    CALL xios_orchidee_send_field ('glcc_pftmtc', glcc_pftmtc) ! kjpindex,nvm,nvmap

  END SUBROUTINE glcc_SinAgeC


! ================================================================================================================================
!! SUBROUTINE   : glcc_SinAgeC_firstday
!!
!>\BRIEF        : When necessary, adjust input glcc matrix, and allocate it
!!                into different contributing age classes and receiving 
!!                youngest age classes.
!! \n
!_ ================================================================================================================================

  ! Note: it has this name because this subroutine will also be called
  ! the first day of each year to precalculate the forest loss for the
  ! deforestation fire module.
  SUBROUTINE glcc_SinAgeC_firstday(npts,veget_max_org,newvegfrac,harvest_matrix,&
                          glccSecondShift,glccPrimaryShift,glccNetLCC,&
                          glccReal,glcc_pft,glcc_pftmtc,IncreDeficit, &
                          Deficit_pf2yf_final, Deficit_sf2yf_final,   &
                          pf2yf_compen_sf2yf, sf2yf_compen_pf2yf)

    IMPLICIT NONE

    !! 0.1 Input variables

    INTEGER, INTENT(in)                                     :: npts           !! Domain size - number of pixels (unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)            :: veget_max_org  !! "maximal" coverage fraction of a PFT on the ground
                                                                              !! May sum to
                                                                              !! less than unity if the pixel has
                                                                              !! nobio area. (unitless, 0-1)
    REAL(r_std), DIMENSION(npts,12),INTENT(in)              :: harvest_matrix !! 
                                                                              !! 
    REAL(r_std), DIMENSION (npts,nvmap),INTENT(in)        :: newvegfrac             !! 
    REAL(r_std), DIMENSION (npts,12),INTENT(in)        :: glccSecondShift     !! the land-cover-change (LCC) matrix in case a gross LCC is 
                                                                              !! used.
    REAL(r_std), DIMENSION (npts,12),INTENT(in)        :: glccPrimaryShift    !! the land-cover-change (LCC) matrix in case a gross LCC is 
                                                                              !! used.
    REAL(r_std), DIMENSION (npts,12),INTENT(in)        :: glccNetLCC          !! the land-cover-change (LCC) matrix in case a gross LCC is 
                                                                              !! used.

    !! 0.2 Output variables
    REAL(r_std), DIMENSION(npts,nvm,nvmap), INTENT(inout)   :: glcc_pftmtc    !! a temporary variable to hold the fractions each PFT is going to lose
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)         :: glcc_pft       !! Loss of fraction in each PFT
    REAL(r_std), DIMENSION(npts,12), INTENT(inout)          :: glccReal       !! The "real" glcc matrix that we apply in the model
                                                                              !! after considering the consistency between presribed
                                                                              !! glcc matrix and existing vegetation fractions.
    REAL(r_std), DIMENSION(npts,12), INTENT(inout)           :: IncreDeficit   !! "Increment" deficits, negative values mean that 
                                                                              !! there are not enough fractions in the source PFTs
                                                                              !! /vegetations to target PFTs/vegetations. I.e., these
                                                                              !! fraction transfers are presribed in LCC matrix but
                                                                              !! not realized.
    REAL(r_std), DIMENSION(npts), INTENT(inout)    :: Deficit_pf2yf_final     !! 
    REAL(r_std), DIMENSION(npts), INTENT(inout)    :: Deficit_sf2yf_final     !! 
    REAL(r_std), DIMENSION(npts), INTENT(inout)    :: pf2yf_compen_sf2yf      !! 
    REAL(r_std), DIMENSION(npts), INTENT(inout)    :: sf2yf_compen_pf2yf      !!
     

    !! 0.3 Modified variables
    
    !! 0.4 Local variables
    REAL(r_std), DIMENSION (npts,12)                :: glcc                !! the land-cover-change (LCC) matrix in case a gross LCC is 
                                                                           !! used.
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

    REAL(r_std), DIMENSION(npts,nvm)         :: veget_max         !! "maximal" coverage fraction of a PFT on the ground
    REAL(r_std), DIMENSION(npts,nvm)         :: veget_max_tmp     !! "maximal" coverage fraction of a PFT on the ground
    REAL(r_std), DIMENSION(npts,nvm)         :: veget_max_old     !! "maximal" coverage fraction of a PFT on the ground
    REAL(r_std), DIMENSION(npts,nvm)         :: glcc_pft_tmp      !! Loss of fraction in each PFT

    ! Different indexes for convenient local uses
    ! We define the rules for gross land cover change matrix:
    ! 1 forest->grass
    ! 2 forest->pasture 
    ! 3 forest->crop
    ! 4 grass->forest
    ! 5 grass->pasture
    ! 6 grass->crop
    ! 7 pasture->forest
    ! 8 pasture->grass
    ! 9 pasture->crop
    ! 10 crop->forest
    ! 11 crop->grass
    ! 12 crop->pasture
    INTEGER :: f2g=1, f2p=2, f2c=3
    INTEGER :: g2f=4, g2p=5, g2c=6, p2f=7, p2g=8, p2c=9, c2f=10, c2g=11, c2p=12

    INTEGER, ALLOCATABLE                  :: indall_tree(:)       !! Indices for all tree PFTs
    INTEGER, ALLOCATABLE                  :: indold_tree(:)       !! Indices for old tree cohort only
    INTEGER, ALLOCATABLE                  :: indagec_tree(:,:)    !! Indices for secondary tree cohorts, 
                                                                  !! note the sequence is old->young.
    INTEGER, ALLOCATABLE                  :: indall_grass(:)      !! Indices for all grass PFTs
    INTEGER, ALLOCATABLE                  :: indold_grass(:)      !! Indices for old grasses only
    INTEGER, ALLOCATABLE                  :: indagec_grass(:,:)   !! Indices for secondary grass cohorts
                                                                  !! note the sequence is old->young.
    INTEGER, ALLOCATABLE                  :: indall_pasture(:)    !! Indices for all pasture PFTs
    INTEGER, ALLOCATABLE                  :: indold_pasture(:)    !! Indices for old pasture only
    INTEGER, ALLOCATABLE                  :: indagec_pasture(:,:) !! Indices for secondary pasture cohorts
                                                                  !! note the sequence is old->young.
    INTEGER, ALLOCATABLE                  :: indall_crop(:)       !! Indices for all crop PFTs
    INTEGER, ALLOCATABLE                  :: indold_crop(:)       !! Indices for old crops only
    INTEGER, ALLOCATABLE                  :: indagec_crop(:,:)    !! Indices for secondary crop cohorts
                                                                  !! note the sequence is old->young.
    INTEGER :: num_tree_sinagec,num_tree_mulagec,num_grass_sinagec,num_grass_mulagec,     &
               num_pasture_sinagec,num_pasture_mulagec,num_crop_sinagec,num_crop_mulagec, &
               itree,itree2,igrass,igrass2,ipasture,ipasture2,icrop,icrop2,pf2yf,sf2yf
    INTEGER :: i,j,ivma,staind,endind,ivm


    REAL(r_std), DIMENSION(npts,12)         :: glccDef            !! Gross LCC deficit, negative values mean that there
                                                                  !! are not enough fractions in the source vegetations
                                                                  !! to the target ones as presribed by the LCC matrix.
    REAL(r_std), DIMENSION(npts)            :: Deficit_pf2yf      !! 
    REAL(r_std), DIMENSION(npts)            :: Deficit_sf2yf      !! 
    REAL(r_std), DIMENSION(npts)            :: Surplus_pf2yf      !! 
    REAL(r_std), DIMENSION(npts)            :: Surplus_sf2yf      !! 
    REAL(r_std), DIMENSION(npts,12)         :: glccRemain      !! 
    REAL(r_std), DIMENSION(npts,12)         :: HmatrixReal        !! 
    INTEGER :: ipts
    

    !! 1. We first build all different indices that we are going to use
    !!    in handling the PFT exchanges, three types of indices are built:
    !!     - for all age classes
    !!     - include only oldest age classes
    !!     - include all age classes excpet the oldest ones
    ! We have to build these indices because we would like to extract from
    ! donating PFTs in the sequnce of old->young age classes, and add in the
    ! receving PFTs only in the youngest-age-class PFTs. These indicies allow
    ! us to know where the different age classes are.

    num_tree_sinagec=0          ! number of tree PFTs with only one single age class 
                                ! considered as the oldest age class
    num_tree_mulagec=0          ! number of tree PFTs having multiple age classes
    num_grass_sinagec=0
    num_grass_mulagec=0
    num_pasture_sinagec=0
    num_pasture_mulagec=0
    num_crop_sinagec=0
    num_crop_mulagec=0
    
    !! 1.1 Calculate the number of PFTs for different MTCs and allocate
    !! the old and all indices arrays.

    ! [Note here the sequence to identify tree,pasture,grass,crop] is
    ! critical. The similar sequence is used in the subroutine "calc_cover".
    ! Do not forget to change the sequence there if you modify here.
    DO ivma =2,nvmap
      staind=start_index(ivma)
      IF (nagec_pft(ivma)==1) THEN
        IF (is_tree(staind)) THEN
          num_tree_sinagec = num_tree_sinagec+1
        ELSE IF (is_grassland_manag(staind)) THEN
          num_pasture_sinagec = num_pasture_sinagec+1
        ELSE IF (natural(staind)) THEN
          num_grass_sinagec = num_grass_sinagec+1
        ELSE
          num_crop_sinagec = num_crop_sinagec+1
        ENDIF

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
    ALLOCATE(indall_tree(num_tree_sinagec+num_tree_mulagec*nagec_tree))     
    ALLOCATE(indall_grass(num_grass_sinagec+num_grass_mulagec*nagec_herb))     
    ALLOCATE(indall_pasture(num_pasture_sinagec+num_pasture_mulagec*nagec_herb))     
    ALLOCATE(indall_crop(num_crop_sinagec+num_crop_mulagec*nagec_herb))     

    ! allocate old-ageclass index
    ALLOCATE(indold_tree(num_tree_sinagec+num_tree_mulagec))     
    ALLOCATE(indold_grass(num_grass_sinagec+num_grass_mulagec))     
    ALLOCATE(indold_pasture(num_pasture_sinagec+num_pasture_mulagec))     
    ALLOCATE(indold_crop(num_crop_sinagec+num_crop_mulagec))     

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

    ! [chaoyuejoy@gmail.com 2015-08-05]
    ! note that we treat the case of (num_tree_mulagec==0) differently. In this
    ! case there is no distinction of age groups among tree PFTs. But we still
    ! we want to use the "gross_lcchange" subroutine. In this case we consider 
    ! them as having a single age group. In the subroutines
    ! of "type_conversion" and "cross_give_receive", only the youngest-age-group
    ! PFTs of a given MTC or vegetation type could receive the incoming fractions.
    ! To be able to handle this case with least amount of code change, we assign the index 
    ! of PFT between youngest and second-oldes (i.e., indagec_tree etc) the same as
    ! those of oldest tree PFTs (or all tree PFTs because in this cases these two indices
    ! are identical) . So that this case could be correctly handled in the subrountines
    ! of "type_conversion" and "cross_give_receive". This treatment allows use
    ! of gross land cover change subroutine with only one single age class. This single
    ! age class is "simultanously the oldest and youngest age class". At the same
    ! time, we also change the num_tree_mulagec as the same of num_crop_sinagec.
    ! The similar case also applies in grass,pasture and crop.

    IF (num_tree_mulagec .EQ. 0) THEN
      ALLOCATE(indagec_tree(num_tree_sinagec,1))
      indagec_tree(:,1) = indall_tree(:)
      num_tree_mulagec = num_tree_sinagec
    ELSE
      ALLOCATE(indagec_tree(num_tree_mulagec,nagec_tree-1))     
    END IF

    IF (num_grass_mulagec .EQ. 0) THEN
      ALLOCATE(indagec_grass(num_grass_sinagec,1))
      indagec_grass(:,1) = indall_grass(:)
      num_grass_mulagec = num_grass_sinagec
    ELSE
      ALLOCATE(indagec_grass(num_grass_mulagec,nagec_herb-1))     
    END IF

    IF (num_pasture_mulagec .EQ. 0) THEN
      ALLOCATE(indagec_pasture(num_pasture_sinagec,1))
      indagec_pasture(:,1) = indall_pasture(:)
      num_pasture_mulagec = num_pasture_sinagec
    ELSE
      ALLOCATE(indagec_pasture(num_pasture_mulagec,nagec_herb-1))
    END IF

    IF (num_crop_mulagec .EQ. 0) THEN
      ALLOCATE(indagec_crop(num_crop_sinagec,1))
      indagec_crop(:,1) = indall_crop(:)
      num_crop_mulagec = num_crop_sinagec
    ELSE
      ALLOCATE(indagec_crop(num_crop_mulagec,nagec_herb-1))
    END IF

    ! fill the non-oldest age class index arrays when number of age classes
    ! is more than 1.
    ! [chaoyuejoy@gmail.com, 2015-08-05]
    ! Note the corresponding part of code  will be automatically skipped
    ! when nagec_tree ==1 and/or nagec_herb ==1, i.e., the assginment
    ! in above codes when original num_*_mulagec variables are zero will be retained.
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


    ! we make copies of original input veget_max
    ! veget_max will be modified through different operations in order to 
    ! check various purposes, e.g., whether input glcc is compatible with
    ! existing veget_max and how to allocate it etc.
    ! veget_max_old will not be modified
    veget_max(:,:) = veget_max_org(:,:)
    veget_max_old(:,:) = veget_max_org(:,:)

    !! 2. Calcuate the fractions covered by tree, grass, pasture and crops
    !!    for each age class

    !************************************************************************!
    !****block to calculate fractions for basic veg types and age classes ***!
    ! Note:
    ! 1. "calc_cover" subroutine does not depend on how many age classes
    ! there are in each MTC.
    ! 2. Fraction of baresoil is excluded here. This means transformation
    ! of baresoil to a vegetated PFT is excluded in gross land cover change.
    veget_mtc(:,:) = 0.
    vegagec_tree(:,:) = 0.
    vegagec_grass(:,:) = 0.
    vegagec_pasture(:,:) = 0.
    vegagec_crop(:,:) = 0.

    CALL calc_cover(npts,veget_max,veget_mtc,vegagec_tree,vegagec_grass, &
           vegagec_pasture,vegagec_crop)
    ! In following call of calc_cover, veget_mtc will be updated each time,
    ! but we don't want this, so we put its initial value into veget_mtc_begin
    ! in order to retrieve this initial value later.
    veget_mtc_begin = veget_mtc
 
    veget_tree(:) = SUM(vegagec_tree(:,:),DIM=2)
    veget_grass(:) = SUM(vegagec_grass(:,:),DIM=2)
    veget_pasture(:) = SUM(vegagec_pasture(:,:),DIM=2)
    veget_crop(:) = SUM(vegagec_crop(:,:),DIM=2)

    !****end block to calculate fractions for basic veg types and age classes ***!
    !****************************************************************************!

    !********************** block to handle forestry harvest ****************
    !! 2B. Here we handle the forestry wood harvest
    ! Rules:
    ! 1. We take first from second oldest forest, then oldest forest
    
    pf2yf=1   !primary to young forest conversion because of harvest
    sf2yf=2   !old secondary to young forest conversion because of harvest
    
    !! Note that Deficit_pf2yf and Deficit_sf2yf are temporary, intermediate
    !! variables. The final deficits after mutual compensation are stored in 
    !! Deficit_pf2yf_final and Deficit_sf2yf_final.
    Deficit_pf2yf(:) = zero 
    Deficit_sf2yf(:) = zero
    Deficit_pf2yf_final(:) = zero 
    Deficit_sf2yf_final(:) = zero

    !! Note that both Surplus_pf2yf and Surplus_sf2yf and temporary intermediate
    !! variables, the final surplus after mutual compensation are not outputed.
    Surplus_pf2yf(:) = zero
    Surplus_sf2yf(:) = zero

    !! Note in the naming of pf2yf_compen_sf2yf and sf2yf_compen_pf2yf, active
    !! tense is used.
    pf2yf_compen_sf2yf(:) = zero  !primary->young conversion that compensates 
                               !the secondary->young conversion because of deficit
                               !in the latter
    sf2yf_compen_pf2yf(:) = zero  !seondary->young conversion that compensates
                               !the primary->young conversion because of the deficit
                               !in the latter
    

    !! Define the "real" harvest matrix after considering the mutual compenstation
    !! between primary->young and secondary->young transitions.
    HmatrixReal(:,:) = zero  !Harvest matrix real, used to hold the 
                                       !harvest matrix after considering the mutual
                                       !compensation between primary and old secondary
                                       !forest

    ! we sum together harvest from primary and secondary forest and consider 
    ! as all happening on parimary forest.
    HmatrixReal(:,1) = harvest_matrix(:,pf2yf) + harvest_matrix(:,sf2yf)

    ! Check the availability of forest fractions for harvest
    WHERE (veget_tree(:) .LE. HmatrixReal(:,1)) 
      Deficit_pf2yf_final(:) = veget_tree(:)-HmatrixReal(:,1)
      HmatrixReal(:,1) = veget_tree(:)
    ENDWHERE

    glccRemain(:,:) = HmatrixReal(:,:)
    glcc_pft(:,:) = 0.
    glcc_pft_tmp(:,:) = 0.
    glcc_pftmtc(:,:,:) = 0.

    !! Allocate harvest-caused out-going primary and secondary forest fraction
    !! into different primary and secondary forest PFTs. 
    ! [Note: here we need only glcc_pft, but not glcc_pft_tmp and glcc_pftmtc.
    ! The latter two variables will be set to zero again when handling LCC in 
    ! later sections.]
    DO ipts=1,npts
      !pf2yf
      CALL type_conversion(ipts,pf2yf,HmatrixReal,veget_mtc,newvegfrac,  &
                       indold_tree,indagec_tree,indagec_crop,num_crop_mulagec, &
                       1,nagec_herb,               &
                       vegagec_tree,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain)
    ENDDO

    ! Because we use the container of type_conversion, now the glcc_pft_tmp
    ! and glcc_pftmtc have wrong information (because harvest loss is assigned
    ! on the newly created youngest-age-class pasture/crop MTCs). So they have 
    ! to be re-initialized to zero. Only the information in glcc_pft is what
    ! we need.
    glcc_pft_tmp(:,:) = 0.
    glcc_pftmtc(:,:,:) = 0.
    !Here we need to put glcc_pft into glcc_pftmtc for forestry harvest.
    !The same MTC will be maintained when forest is harvested.
    DO ivm =1,nvm
      IF (is_tree(ivm)) THEN
        glcc_pftmtc(:,ivm,pft_to_mtc(ivm)) = glcc_pft(:,ivm)
      ENDIF 
    ENDDO
    !****************** end block to handle forestry harvest ****************
    veget_max_tmp(:,:) = veget_max(:,:)


    !************************************************************************!
    !****block to calculate fractions for basic veg types and age classes ***!
    ! Note:
    ! 1. "calc_cover" subroutine does not depend on how many age classes
    ! there are in each MTC.
    ! 2. Fraction of baresoil is excluded here. This means transformation
    ! of baresoil to a vegetated PFT is excluded in gross land cover change.
    veget_mtc(:,:) = 0.
    vegagec_tree(:,:) = 0.
    vegagec_grass(:,:) = 0.
    vegagec_pasture(:,:) = 0.
    vegagec_crop(:,:) = 0.


    CALL calc_cover(npts,veget_max,veget_mtc,vegagec_tree,vegagec_grass, &
           vegagec_pasture,vegagec_crop)
    veget_mtc = veget_mtc_begin
 
    veget_tree(:) = SUM(vegagec_tree(:,:),DIM=2)
    veget_grass(:) = SUM(vegagec_grass(:,:),DIM=2)
    veget_pasture(:) = SUM(vegagec_pasture(:,:),DIM=2)
    veget_crop(:) = SUM(vegagec_crop(:,:),DIM=2)
    itree=1
    igrass=2
    ipasture=3
    icrop=4
    veget_4veg(:,itree) = veget_tree(:)
    veget_4veg(:,igrass) = veget_grass(:)
    veget_4veg(:,ipasture) = veget_pasture(:)
    veget_4veg(:,icrop) = veget_crop(:)
    !****end block to calculate fractions for basic veg types and age classes ***!
    !****************************************************************************!

    !! 3. Decompose the LCC matrix to different PFTs
    !! We do this through several steps:
    !  3.1 Check whether input LCC matrix is feasible with current PFT fractions
    !      (i.e., the fractions of forest,grass,pasture and crops)
    !      and if not, adjust the transfer matrix by compensating the deficits
    !      using the surpluses.
    !  3.2 Allocate the decreasing fractions of tree/grass/pasture/crop to their
    !      respective age classes, in the sequences of old->young.
    !  3.3 Allocate the incoming fractions of tree/grass/pasture/crop to their
    !      respective youngest age classes. The incoming fractions are distributed
    !      according to the existing fractions of youngest-age-class PFTs of the
    !      same receiving vegetation type. If none of them exists, the incoming
    !      fraction is distributed equally.

    !!  3.1 Adjust LCC matrix if it's not feasible with current PFT fractions

    !++code freezing++
    !codes below handle the mutual compenstation of transition matrices
    !among different land cover types. This is desgined for consistency 
    !with activated DGVM.    

    ! glcc(:,:) = glccSecondShift+glccPrimaryShift+glccNetLCC
    ! glccReal(:,:) = 0.
    ! glccDef(:,:) = 0.

    ! !to crop - sequence: p2c,g2c,f2c
    ! CALL glcc_compensation_full(npts,veget_4veg,glcc,glccReal,glccDef, &
    !                        p2c,ipasture,g2c,igrass,f2c,itree,icrop, &
    !                        IncreDeficit)

    ! !to pasture - sequence: g2p,c2p,f2p
    ! CALL glcc_compensation_full(npts,veget_4veg,glcc,glccReal,glccDef, &
    !                        g2p,igrass,c2p,icrop,f2p,itree,ipasture, &
    !                        IncreDeficit)

    ! !to grass - sequence: p2g,c2g,f2g
    ! CALL glcc_compensation_full(npts,veget_4veg,glcc,glccReal,glccDef, &
    !                        p2g,ipasture,c2g,icrop,f2g,itree,igrass, &
    !                        IncreDeficit)

    ! !to forest - sequence: c2f,p2f,g2f
    ! CALL glcc_compensation_full(npts,veget_4veg,glcc,glccReal,glccDef, &
    !                        c2f,icrop,p2f,ipasture,g2f,igrass,itree, &
    !                        IncreDeficit)

    ! !!  3.2 & 3.3 Allocate LCC matrix to different PFTs/age-classes

    ! ! because we use veget_max as a proxy variable and it has been changed
    ! ! when we derive the glccReal, so here we have to recover its original
    ! ! values, which is veget_max_tmp after the forestry harvest.
    ! veget_max(:,:) = veget_max_tmp(:,:)

    ! ! Calculate again fractions for different age-classes.
    ! veget_mtc(:,:) = 0.
    ! vegagec_tree(:,:) = 0.
    ! vegagec_grass(:,:) = 0.
    ! vegagec_pasture(:,:) = 0.
    ! vegagec_crop(:,:) = 0.

    ! CALL calc_cover(npts,veget_max,veget_mtc,vegagec_tree,vegagec_grass, &
    !        vegagec_pasture,vegagec_crop)


    !++end codes freezing ++

    IncreDeficit(:,:) = 0.
    glcc(:,:) = glccSecondShift+glccPrimaryShift+glccNetLCC
    glccReal(:,:) = glcc(:,:)
    glccRemain(:,:) = glcc(:,:)

    !  We allocate in the sequences of old->young. Within the same age-class
    !  group, we allocate in proportion with existing PFT fractions.
    DO ipts=1,npts
      !f2c
      CALL type_conversion(ipts,f2c,glccReal,veget_mtc,newvegfrac,       &
                       indold_tree,indagec_tree,indagec_crop,num_crop_mulagec,     &
                       nagec_tree,nagec_herb,                    &
                       vegagec_tree,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp, &
                       glccRemain)
      !f2p
      CALL type_conversion(ipts,f2p,glccReal,veget_mtc,newvegfrac,       &
                       indold_tree,indagec_tree,indagec_pasture,num_pasture_mulagec,     &
                       nagec_tree,nagec_herb,                    &
                       vegagec_tree,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain)
      !f2g
      CALL type_conversion(ipts,f2g,glccReal,veget_mtc,newvegfrac,       &
                       indold_tree,indagec_tree,indagec_grass,num_grass_mulagec,     &
                       nagec_tree,nagec_herb,                    &
                       vegagec_tree,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain)
      !g2c
      CALL type_conversion(ipts,g2c,glccReal,veget_mtc,newvegfrac,       &
                       indold_grass,indagec_grass,indagec_crop,num_crop_mulagec,     &
                       nagec_herb,nagec_herb,                    &
                       vegagec_grass,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain)
      !g2p
      CALL type_conversion(ipts,g2p,glccReal,veget_mtc,newvegfrac,       &
                       indold_grass,indagec_grass,indagec_pasture,num_pasture_mulagec,     &
                       nagec_herb,nagec_herb,                    &
                       vegagec_grass,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain)
      !g2f
      CALL type_conversion(ipts,g2f,glccReal,veget_mtc,newvegfrac,       &
                       indold_grass,indagec_grass,indagec_tree,num_tree_mulagec,     &
                       nagec_herb,nagec_tree,                    &
                       vegagec_grass,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain)
      !p2c
      CALL type_conversion(ipts,p2c,glccReal,veget_mtc,newvegfrac,       &
                       indold_pasture,indagec_pasture,indagec_crop,num_crop_mulagec,     &
                       nagec_herb,nagec_herb,                    &
                       vegagec_pasture,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain)
      !p2g
      CALL type_conversion(ipts,p2g,glccReal,veget_mtc,newvegfrac,       &
                       indold_pasture,indagec_pasture,indagec_grass,num_grass_mulagec,     &
                       nagec_herb,nagec_herb,                    &
                       vegagec_pasture,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain)
      !p2f
      CALL type_conversion(ipts,p2f,glccReal,veget_mtc,newvegfrac,       &
                       indold_pasture,indagec_pasture,indagec_tree,num_tree_mulagec,     &
                       nagec_herb,nagec_tree,                    &
                       vegagec_pasture,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain)
      !c2p
      CALL type_conversion(ipts,c2p,glccReal,veget_mtc,newvegfrac,       &
                       indold_crop,indagec_crop,indagec_pasture,num_pasture_mulagec,     &
                       nagec_herb,nagec_herb,                    &
                       vegagec_crop,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain)
      !c2g
      CALL type_conversion(ipts,c2g,glccReal,veget_mtc,newvegfrac,       &
                       indold_crop,indagec_crop,indagec_grass,num_grass_mulagec,     &
                       nagec_herb,nagec_herb,                    &
                       vegagec_crop,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain)
      !c2f
      CALL type_conversion(ipts,c2f,glccReal,veget_mtc,newvegfrac,       &
                       indold_crop,indagec_crop,indagec_tree,num_tree_mulagec,     &
                       nagec_herb,nagec_tree,                    &
                       vegagec_crop,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain)
    ENDDO
   
    WHERE (glccRemain .GT. zero) 
      glccReal = glcc - glccRemain
      IncreDeficit = -1 * glccRemain
    ENDWHERE 

  END SUBROUTINE glcc_SinAgeC_firstday



! ================================================================================================================================
!! SUBROUTINE   : type_conversion
!>\BRIEF        : Allocate outgoing into different age classes and incoming into
!!                yongest age-class of receiving MTCs.
!!
!! REMARK       : The current dummy variables give an example of converting forests
!!                to crops.
!! \n
!_ ================================================================================================================================
  SUBROUTINE type_conversion(ipts,f2c,glccReal,veget_mtc,newvegfrac,       &
                     indold_tree,indagec_tree,indagec_crop,num_crop_mulagec,     &
                     nagec_giving,nagec_receive,                    &
                     vegagec_tree,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp, &
                     glccRemain, &
                     iagec_start)

    IMPLICIT NONE

    !! Input variables
    INTEGER, INTENT(in)                             :: ipts,f2c
    REAL(r_std), DIMENSION(:,:), INTENT(in)         :: glccReal             !! The "real" glcc matrix that we apply in the model
                                                                            !! after considering the consistency between presribed
                                                                            !! glcc matrix and existing vegetation fractions.
    REAL(r_std), DIMENSION(:,:), INTENT(in)         :: veget_mtc            !! "maximal" coverage fraction of a PFT on the ground
    REAL(r_std), DIMENSION(:,:),INTENT(in)          :: newvegfrac           !! 
    INTEGER, DIMENSION(:), INTENT(in)               :: indold_tree          !! Indices for PFTs giving out fractions; 
                                                                            !! here use old tree cohort as an example
    INTEGER, DIMENSION(:,:), INTENT(in)             :: indagec_tree         !! Indices for PFTs giving out fractions; 
                                                                            !! here use old tree cohort as an example
    INTEGER, DIMENSION(:,:), INTENT(in)             :: indagec_crop         !! Indices for secondary basic-vegetation cohorts; The youngest age classes
                                                                            !! of these vegetations are going to receive fractions. 
                                                                            !! here we use crop cohorts as an example
    INTEGER, INTENT(in)                             :: num_crop_mulagec     !! number of crop MTCs with more than one age classes
    INTEGER, INTENT(in)                             :: nagec_giving         !! number of age classes in the giving basic types
                                                                            !! (i.e., tree, grass, pasture, crop), here we can use tree
                                                                            !! as an example, nagec=nagec_tree
    INTEGER, INTENT(in)                             :: nagec_receive        !! number of age classes in the receiving basic types
                                                                            !! (i.e., tree, grass, pasture, crop), here we can use crop
                                                                            !! as an example, nagec=nagec_herb
    INTEGER, OPTIONAL, INTENT(in)                   :: iagec_start          !! starting index for iagec, this is added in order to handle
                                                                            !! the case of secondary forest harvest.

    !! 1. Modified variables
    REAL(r_std), DIMENSION(:,:), INTENT(inout)      :: vegagec_tree         !! fraction of tree age-class groups, in sequence of old->young
    REAL(r_std), DIMENSION(:,:), INTENT(inout)      :: veget_max            !! "maximal" coverage fraction of a PFT on the ground
    REAL(r_std), DIMENSION(:,:), INTENT(inout)      :: glcc_pft             !! a temporary variable to hold the fractions each PFT is going to lose
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)    :: glcc_pftmtc          !! a temporary variable to hold the fraction of ipft->ivma, i.e., from 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)      :: glcc_pft_tmp         !! Loss of fraction in each PFT
    REAL(r_std), DIMENSION(:,:), INTENT(inout)      :: glccRemain           !! The remaining glcc matrix after applying the conversion. I.e., it will
                                                                            !! record the remaining unrealized transition fraction in case the donation
                                                                            !! vegetation is not enough compared with prescribed transition fraction.
                                                                            !! This variable should be initialized the same as glccReal before it's fed
                                                                            !! into this function.

    !! Local vriables
    INTEGER  :: j,iagec,iagec_start_proxy
    REAL(r_std) :: frac_begin,frac_used
                                                                            !! PFT_{ipft} to the youngest age class of MTC_{ivma}
    IF (.NOT. PRESENT(iagec_start)) THEN
      iagec_start_proxy=1
    ELSE
      iagec_start_proxy=iagec_start
    ENDIF
    
    ! This subroutine handles the conversion from one basic-vegetation type
    ! to another, by calling the subroutine cross_give_receive, which handles
    ! allocation of giving-receiving fraction among the giving age classes
    ! and receiving basic-vegetation young age classes.
    ! We allocate in the sequences of old->young. Within the same age-class
    ! group, we allocate in proportion with existing PFT fractions. The same
    ! also applies in the receiving youngest-age-class PFTs, i.e., the receiving
    ! total fraction is allocated according to existing fractions of 
    ! MTCs of the same basic vegetation type, otherwise it will be equally
    ! distributed.

    frac_begin = glccReal(ipts,f2c)
    !DO WHILE (frac_begin>min_stomate)
      DO iagec=iagec_start_proxy,nagec_giving
        IF (vegagec_tree(ipts,iagec)>frac_begin) THEN
          frac_used = frac_begin
        ELSE IF (vegagec_tree(ipts,iagec)>min_stomate) THEN
          frac_used = vegagec_tree(ipts,iagec)
        ELSE
          frac_used = 0.
        ENDIF
        
        IF (frac_used>min_stomate) THEN
          IF (iagec==1) THEN
            ! Note that vegagec_tree is fractions of tree age-class groups in the 
            ! the sequence of old->young, so iagec==1 means that we're handling 
            ! first the oldest-age-group tree PFTs.
            CALL cross_give_receive(ipts,frac_used,veget_mtc,newvegfrac,      &
                     indold_tree,indagec_crop,nagec_receive,num_crop_mulagec, &
                      veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp)
          ELSE
            ! Note also the sequence of indagec_tree is from old->young, so by
            ! increasing iagec, we're handling progressively the old to young
            ! tree age-class PFTs.
            CALL cross_give_receive(ipts,frac_used,veget_mtc,newvegfrac,      &
                     indagec_tree(:,iagec-1),indagec_crop,nagec_receive,num_crop_mulagec, &
                      veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp)
          ENDIF
          frac_begin = frac_begin-frac_used
          vegagec_tree(ipts,iagec)=vegagec_tree(ipts,iagec)-frac_used
          glccRemain(ipts,f2c) = glccRemain(ipts,f2c) - frac_used
        ENDIF
      ENDDO
    !ENDDO

  END SUBROUTINE type_conversion

END MODULE stomate_glcchange_SinAgeC
