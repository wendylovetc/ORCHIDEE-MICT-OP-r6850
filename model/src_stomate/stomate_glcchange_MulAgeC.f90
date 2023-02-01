! =================================================================================================================================
! MODULE       : stomate_glcchange_MulAgeC
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


MODULE stomate_glcchange_MulAgeC

  ! modules used:
  
  USE ioipsl_para
  USE stomate_data
  USE pft_parameters
  USE constantes
  USE constantes_soil_var
  USE stomate_gluc_common
  USE stomate_gluc_constants
  USE xios_orchidee

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC glcc_MulAgeC_firstday, glcc_MulAgeC, age_class_distr, type_conversion
  
CONTAINS

! ================================================================================================================================
!! SUBROUTINE   : age_class_distr
!!
!>\BRIEF        Redistribute biomass, litter, soilcarbon and water across
!!              the age classes
!!
!! DESCRIPTION  : Following growth, the trees from an age class may have become
!! too big to belong to this age class. The biomass, litter, soilcarbon and
!! soil water then need to be moved from one age class to the next age class.
!!
!! RECENT CHANGE(S) : 
!!
!! MAIN OUTPUT VARIABLE(S) :  
!!
!! REFERENCES   : None
!!
!! FLOWCHART    : 
!! \n
!_ ================================================================================================================================

  SUBROUTINE age_class_distr(npts, lalo, resolution, bound_spa, &  
       biomass, veget_max, ind, &
       bm_phytomer,bm_FFB,PHYbm, FFBbm,                  & !! yidi
       lm_lastyearmax, leaf_frac, co2_to_bm, &
       fuel_1hr, fuel_10hr, fuel_100hr, fuel_1000hr, &
       everywhere, litter, carbon, lignin_struc, &
       deepC_a, deepC_s, deepC_p, &
       bm_to_litter, PFTpresent, when_growthinit,&
       senescence, npp_longterm, gpp_daily, leaf_age, age, &
       gdd_from_growthinit, gdd_midwinter, time_hum_min, hum_min_dormance, &
       gdd_m5_dormance, &
       ncd_dormance, moiavail_month, moiavail_week, ngd_minus5, &
       gpp_week, resp_maint, resp_growth, npp_daily, resp_hetero)

    IMPLICIT NONE
    
  !! 0. Variable and parameter declaration 
    
    !! 0.1 Input variables

    INTEGER, INTENT(in)                                :: npts                !! Domain size - number of pixels (unitless)
    REAL(r_std),DIMENSION(npts,2),INTENT(in)                   :: lalo                 !! Geographical coordinates (latitude,longitude)
                                                                                       !! for pixels (degrees)
    REAL(r_std), DIMENSION(npts,2), INTENT(in)                 :: resolution           !! Resolution at each grid point (m)  
                                                                                       !! [1=E-W, 2=N-S] 

    !! 0.2 Output variables


    !! 0.3 Modified variables

    LOGICAL, DIMENSION(:,:), INTENT(inout)             :: PFTpresent          !! Tab indicating which PFTs are present in 
                                                                              !! each pixel
    LOGICAL, DIMENSION(:,:), INTENT(inout)             :: senescence          !! Flag for setting senescence stage (only 
                                                                              !! for deciduous trees)
     REAL(r_std), DIMENSION(:,:), INTENT(inout)        :: moiavail_month      !! "Monthly" moisture availability (0 to 1, 
                                                                              !! unitless) 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: moiavail_week       !! "Weekly" moisture availability 
                                                                              !! (0 to 1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: gpp_week            !! Mean weekly gross primary productivity 
                                                                              !! @tex $(gC m^{-2} day^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: ngd_minus5          !! Number of growing days (days), threshold 
                                                                              !! -5 deg C (for phenology)   
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: resp_maint          !! Maintenance respiration  
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: resp_growth         !! Growth respiration  
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: npp_daily           !! Net primary productivity 
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: resp_hetero           !! Net primary productivity 
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: when_growthinit     !! How many days ago was the beginning of 
                                                                              !! the growing season (days)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: npp_longterm        !! "Long term" mean yearly primary productivity 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: ind                 !! Number of individuals at the stand level
                                                                              !! @tex $(m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: veget_max           !! "maximal" coverage fraction of a PFT on the ground
                                                                              !! May sum to
                                                                              !! less than unity if the pixel has
                                                                              !! nobio area. (unitless, 0-1)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: lm_lastyearmax      !! last year's maximum leaf mass for each PFT 
                                                                              !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: everywhere          !! is the PFT everywhere in the grid box or 
                                                                              !! very localized (after its introduction) (?)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: co2_to_bm           !! CO2 taken from the atmosphere to get C to create  
                                                                              !! the seedlings @tex (gC.m^{-2}dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: gpp_daily           !! Daily gross primary productivity  
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: time_hum_min        !! Time elapsed since strongest moisture 
                                                                              !! availability (days) 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: hum_min_dormance    !! minimum moisture during dormance 
                                                                              !! (0-1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: gdd_midwinter       !! Growing degree days (K), since midwinter 
                                                                              !! (for phenology) - this is written to the
                                                                              !!  history files 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: gdd_from_growthinit !! growing degree days, since growthinit 
                                                                              !! for crops
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: gdd_m5_dormance     !! Growing degree days (K), threshold -5 deg 
                                                                              !! C (for phenology)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: ncd_dormance        !! Number of chilling days (days), since 
                                                                              !! leaves were lost (for phenology) 
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: lignin_struc        !! ratio Lignine/Carbon in structural litter,
                                                                              !! above and below ground
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: carbon              !! carbon pool: active, slow, or passive 
                                                                              !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: deepC_a             !! Permafrost soil carbon (g/m**3) active
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: deepC_s             !! Permafrost soil carbon (g/m**3) slow
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: deepC_p             !! Permafrost soil carbon (g/m**3) passive
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: age                 !! Age (years)    
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: leaf_frac           !! fraction of leaves in leaf age class (unitless;0-1)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: leaf_age            !! Leaf age (days)
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: bm_to_litter        !! Transfer of biomass to litter 
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: biomass             !! Stand level biomass @tex $(gC.m^{-2})$ @endtex
!! yidi
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: FFBbm               !! FFB mass, from sapabove
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: PHYbm               !! PHYTOMER mass, from sapabove
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: bm_phytomer         !! Each PHYTOMER mass, from sapabove
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: bm_FFB              !! Fruit mass for Each PHYTOMER
!! yidi
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)   :: litter              !! metabolic and structural litter, above and 
                                                                              !! below ground @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: fuel_1hr
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: fuel_10hr
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: fuel_100hr
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: fuel_1000hr
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: bound_spa           !! Spatial age class boundaries.

    !! 0.4 Local variables

    INTEGER(i_std)                                     :: ipts,ivm,igroup     !! Indeces(unitless)
    INTEGER(i_std)                                     :: iele,ipar,ipft      !! Indeces(unitless)
    INTEGER(i_std)                                     :: iagec,imbc,icirc    !! Indeces(unitless)
    INTEGER(i_std)                                     :: ilit,ilev,icarb     !! Indeces(unitless)
    INTEGER(i_std)                                     :: ivma                !! Indeces(unitless)
    REAL(r_std)                                        :: share_expanded      !! Share of the veget_max of the existing vegetation
                                                                              !! within a PFT over the total veget_max following 
                                                                              !! expansion of that PFT (unitless, 0-1)
                                                                              !! @tex $(ind m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nmbcomp,nelements) :: check_intern        !! Contains the components of the internal
                                                                              !! mass balance chech for this routine
                                                                              !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)         :: closure_intern      !! Check closure of internal mass balance
                                                                              !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)         :: pool_start          !! Start and end pool of this routine 
                                                                              !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)         :: pool_end            !! Start and end pool of this routine 
                                                                              !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(nelements)                  :: temp_start          !! Start and end pool of this routine 
                                                                              !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(nelements)                  :: temp_end            !! Start and end pool of this routine 
                                                                              !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(nlitt,nlevs)                :: litter_weight_expanded !! The fraction of litter on the expanded
                                                                              !! PFT.
                                                                              !! @tex $-$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)                   :: woodmass            !! Woodmass of individuals (gC) 
    REAL(r_std), DIMENSION(npts,nvm)                   :: reverse_soilcarbon          !! 
    REAL(r_std), DIMENSION(npts,nvm)                   :: agec_indicator      !! 
    CHARACTER(LEN=80)                                  :: data_filename

!_ ================================================================================================================================
  
    IF (printlev.GE.3) WRITE(numout,*) 'Entering age class distribution'

    !CALL getin_p('AgeC_Threshold_File',data_filename)
    !CALL slowproc_read_data(npts, lalo, resolution, bound_spa, data_filename, 'matrix')

    IF (.NOT. use_bound_spa) THEN
      DO ipts = 1,npts
        bound_spa(ipts,:) = age_class_bound(:)
      ENDDO
    ENDIF

 !! 1. Initialize

    woodmass(:,:) = biomass(:,:,isapabove,icarbon)+biomass(:,:,isapbelow,icarbon) &
                    +biomass(:,:,iheartabove,icarbon)+biomass(:,:,iheartbelow,icarbon) 
    reverse_soilcarbon(:,:) = -1 *(SUM(carbon(:,:,:),DIM=2) + &
                      SUM(SUM(litter(:,:,:,:,icarbon),DIM=2),DIM=3))

    !! 1.2 Initialize check for mass balance closure
    !  The mass balance is calculated at the end of this routine
    !  in section 3. Initial biomass and wood product pool all other
    !  relevant pools were just set to zero.
    pool_start(:,:,:) = zero
    DO iele = 1,nelements
       
       ! co2_to_bm
       pool_start(:,:,iele) = pool_start(:,:,iele) + co2_to_bm(:,:)

       ! Biomass pool + bm_to_litter
       DO ipar = 1,nparts
          pool_start(:,:,iele) = pool_start(:,:,iele) + &
               (biomass(:,:,ipar,iele) + bm_to_litter(:,:,ipar,iele)) * &
               veget_max(:,:)
       ENDDO

       ! Litter pool (gC m-2) *  (m2 m-2) 
       DO ilit = 1,nlitt
          DO ilev = 1,nlevs
             pool_start(:,:,iele) = pool_start(:,:,iele) + &
                  litter(:,ilit,:,ilev,iele) * veget_max(:,:)
          ENDDO
       ENDDO

       ! Soil carbon (gC m-2) *  (m2 m-2)
       DO icarb = 1,ncarb
          pool_start(:,:,iele) = pool_start(:,:,iele) + &
               carbon(:,icarb,:) * veget_max(:,:)
       ENDDO

    ENDDO


 !! 2. Handle the merge of PFTs when one age class moves to the next one.

    !  Following growth, the value of age-class indicator variable 
    !  from an age class may have become too big to stay
    !  in this age class. The biomass, litter, reverse_soilcarbon and soil
    !  water then need to be moved from one age class to the next age class.
    DO ipts = 1,npts
      ! This loops over all the MTCs that we have ignoring age classes
      DO ivma=1,nvmap
        ivm=start_index(ivma)

        ! If we only have a single age class for this
        ! PFT, we can skip it.
        IF(nagec_pft(ivma) .EQ. 1)CYCLE

        IF(is_tree(ivm)) THEN
          agec_indicator(:,:) = woodmass(:,:)
        ELSE
          agec_indicator(:,:) = reverse_soilcarbon(:,:)
        ENDIF ! is_tree(ivm)

        CALL check_merge_same_MTC(ipts, ivma, agec_indicator, bound_spa, &
                biomass, veget_max, ind, &
                bm_phytomer,bm_FFB,PHYbm, FFBbm,npts,                  & !! yidi
                lm_lastyearmax, leaf_frac, co2_to_bm, &
                fuel_1hr, fuel_10hr, fuel_100hr, fuel_1000hr, &
                everywhere, litter, carbon, lignin_struc, &
                deepC_a, deepC_s, deepC_p, &
                bm_to_litter, PFTpresent, when_growthinit,&
                senescence, npp_longterm, gpp_daily, leaf_age, age, &
                gdd_from_growthinit, gdd_midwinter, time_hum_min, hum_min_dormance, &
                gdd_m5_dormance, &
                ncd_dormance, moiavail_month, moiavail_week, ngd_minus5, &
                gpp_week, resp_maint, resp_growth, npp_daily, resp_hetero)

      ENDDO ! Looping over MTCs
    ENDDO ! loop over #pixels - domain size


!! 3. Mass balance closure
    
    !! 3.1 Calculate components of the mass balance
    pool_end(:,:,:) = zero
    DO iele = 1,nelements

       ! co2_to_bm
       pool_end(:,:,iele) = pool_end(:,:,iele) + co2_to_bm(:,:)

       ! Biomass pool + bm_to_litter
       DO ipar = 1,nparts
          pool_end(:,:,iele) = pool_end(:,:,iele) + &
               (biomass(:,:,ipar,iele) + bm_to_litter(:,:,ipar,iele)) * &
               veget_max(:,:)
       ENDDO

       ! Litter pool (gC m-2) *  (m2 m-2) 
       DO ilit = 1,nlitt
          DO ilev = 1,nlevs
             pool_end(:,:,iele) = pool_end(:,:,iele) + &
                  litter(:,ilit,:,ilev,iele) * veget_max(:,:)
          ENDDO
       ENDDO

       ! Soil carbon (gC m-2) *  (m2 m-2)
       DO icarb = 1,ncarb
          pool_end(:,:,iele) = pool_end(:,:,iele) + &
               carbon(:,icarb,:) * veget_max(:,:)
       ENDDO
    ENDDO

    !! 3.2 Calculate mass balance
    check_intern(:,:,iatm2land,icarbon) = zero 
    check_intern(:,:,iland2atm,icarbon) = -un * zero
    check_intern(:,:,ilat2out,icarbon) = zero
    check_intern(:,:,ilat2in,icarbon) = -un * zero
    check_intern(:,:,ipoolchange,icarbon) = -un * (pool_end(:,:,icarbon) - pool_start(:,:,icarbon))
    closure_intern = zero
    DO imbc = 1,nmbcomp
       closure_intern(:,:,icarbon) = closure_intern(:,:,icarbon) + check_intern(:,:,imbc,icarbon)
    ENDDO

    !! 3.3 Write outcome of the check
    !  Sum over ivm because of age class redistribution
    DO ipts = 1,npts
       IF (SUM(closure_intern(ipts,:,icarbon)) .LT. min_stomate .AND. &
            SUM(closure_intern(ipts,:,icarbon)) .GT. -min_stomate) THEN
          IF (ld_massbal) WRITE(numout,*) 'Mass balance closure: age_class_distr', ipts
       ELSE
          WRITE(numout,*) 'Error: mass balance is not closed in age_class_distr'
          WRITE(numout,*) '   Difference, ipts, ', ipts, SUM(closure_intern(ipts,:,icarbon)) 
       ENDIF
    ENDDO

    IF (printlev.GE.4) WRITE(numout,*) 'Leaving age class distribution'
    
  END SUBROUTINE age_class_distr


  SUBROUTINE check_merge_same_MTC(ipts, ivma, woodmass, bound_spa, &
       biomass, veget_max, ind, &
       bm_phytomer,bm_FFB,PHYbm, FFBbm,npts,                  & !! yidi
       lm_lastyearmax, leaf_frac, co2_to_bm, &
       fuel_1hr, fuel_10hr, fuel_100hr, fuel_1000hr, &
       everywhere, litter, carbon, lignin_struc, &
       deepC_a, deepC_s, deepC_p, &
       bm_to_litter, PFTpresent, when_growthinit,&
       senescence, npp_longterm, gpp_daily, leaf_age, age, &
       gdd_from_growthinit, gdd_midwinter, time_hum_min,hum_min_dormance, &
       gdd_m5_dormance, &
       ncd_dormance, moiavail_month, moiavail_week, ngd_minus5, &
       gpp_week, resp_maint, resp_growth, npp_daily, resp_hetero)

    IMPLICIT NONE
    
  !! 0. Variable and parameter declaration 
    
    !! 0.1 Input variables

    INTEGER, INTENT(in)                                :: npts                !! Domain size - number of pixels (unitless)
    INTEGER, INTENT(in)                                :: ipts                !! Domain size - number of pixels (unitless)
    INTEGER, INTENT(in)                                :: ivma                !! 
    REAL(r_std), DIMENSION(:,:), INTENT(in)            :: woodmass            !! Woodmass of individuals (gC) 
    REAL(r_std), DIMENSION(:,:), INTENT(in)            :: bound_spa           !!

    !! 0.2 Output variables


    !! 0.3 Modified variables

    LOGICAL, DIMENSION(:,:), INTENT(inout)             :: PFTpresent          !! Tab indicating which PFTs are present in 
                                                                              !! each pixel
    LOGICAL, DIMENSION(:,:), INTENT(inout)             :: senescence          !! Flag for setting senescence stage (only 
                                                                              !! for deciduous trees)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)        :: moiavail_month       !! "Monthly" moisture availability (0 to 1, 
                                                                              !! unitless) 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: moiavail_week       !! "Weekly" moisture availability 
                                                                              !! (0 to 1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: gpp_week            !! Mean weekly gross primary productivity 
                                                                              !! @tex $(gC m^{-2} day^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: ngd_minus5          !! Number of growing days (days), threshold 
                                                                              !! -5 deg C (for phenology)   
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: resp_maint          !! Maintenance respiration  
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: resp_growth         !! Growth respiration  
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: npp_daily           !! Net primary productivity 
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: resp_hetero         !! Net primary productivity 
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: when_growthinit     !! How many days ago was the beginning of 
                                                                              !! the growing season (days)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: npp_longterm        !! "Long term" mean yearly primary productivity 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: ind                 !! Number of individuals at the stand level
                                                                              !! @tex $(m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: veget_max           !! "maximal" coverage fraction of a PFT on the ground
                                                                              !! May sum to
                                                                              !! less than unity if the pixel has
                                                                              !! nobio area. (unitless, 0-1)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: lm_lastyearmax      !! last year's maximum leaf mass for each PFT 
                                                                              !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: everywhere          !! is the PFT everywhere in the grid box or 
                                                                              !! very localized (after its introduction) (?)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: co2_to_bm           !! CO2 taken from the atmosphere to get C to create  
                                                                              !! the seedlings @tex (gC.m^{-2}dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: gpp_daily           !! Daily gross primary productivity  
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: time_hum_min        !! Time elapsed since strongest moisture 
                                                                              !! availability (days) 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: hum_min_dormance    !! minimum moisture during dormance 
                                                                              !! (0-1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: gdd_midwinter       !! Growing degree days (K), since midwinter 
                                                                              !! (for phenology) - this is written to the
                                                                              !!  history files 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: gdd_from_growthinit !! growing degree days, since growthinit 
                                                                              !! for crops
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: gdd_m5_dormance     !! Growing degree days (K), threshold -5 deg 
                                                                              !! C (for phenology)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: ncd_dormance        !! Number of chilling days (days), since 
                                                                              !! leaves were lost (for phenology) 
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: lignin_struc        !! ratio Lignine/Carbon in structural litter,
                                                                              !! above and below ground
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: carbon              !! carbon pool: active, slow, or passive 
                                                                              !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: deepC_a             !! Permafrost soil carbon (g/m**3) active
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: deepC_s             !! Permafrost soil carbon (g/m**3) slow
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: deepC_p             !! Permafrost soil carbon (g/m**3) passive
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: age                 !! Age (years)    
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: leaf_frac           !! fraction of leaves in leaf age class (unitless;0-1)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: leaf_age            !! Leaf age (days)
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: bm_to_litter        !! Transfer of biomass to litter 
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: biomass             !! Stand level biomass @tex $(gC.m^{-2})$ @endtex
!! yidi
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: FFBbm               !! FFB mass, from sapabove
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: PHYbm               !! PHYTOMER mass, from sapabove
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: bm_phytomer         !! Each PHYTOMER mass, from sapabove
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: bm_FFB              !! Fruit mass for Each PHYTOMER
!! yidi
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)   :: litter              !! metabolic and structural litter, above and 
                                                                              !! below ground @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: fuel_1hr
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: fuel_10hr
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: fuel_100hr
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: fuel_1000hr

    !! 0.4 Local variables

    INTEGER(i_std)                                     :: iele,ipar,ipft      !! Indeces(unitless)
    INTEGER(i_std)                                     :: iagec,imbc,icirc    !! Indeces(unitless)
    INTEGER(i_std)                                     :: ilit,ilev,icarb     !! Indeces(unitless)
    REAL(r_std)                                        :: share_expanded      !! Share of the veget_max of the existing vegetation
                                                                              !! within a PFT over the total veget_max following 
                                                                              !! expansion of that PFT (unitless, 0-1)
                                                                              !! @tex $(ind m^{-2})$ @endtex
    REAL(r_std), DIMENSION(nlitt,nlevs)                :: litter_weight_expanded !! The fraction of litter on the expanded
                                                                              !! PFT.
!! yidi
   REAL(r_std), DIMENSION(npts,nvm,nelements)          :: tot_biomass    !! Total living biomass
!! yidi


!_ ================================================================================================================================

    !! 1 Check if the trees still belong to this age class
    !  Note that the term age class is used but that the classes used in the
    !  code are not defined on an age criterion. Instead the biomass or
    !  or soil carbon pool is used.
    WRITE(numout,*) 'Entering checkmerge'
    WRITE(numout,*) 'ivma',ivma
    WRITE(numout,*) 'start_index(ivma)',start_index(ivma)
  IF (is_tree(start_index(ivma))) THEN
    
    WRITE(numout,*) 'bound_spa',bound_spa(ipts,ivma)
    WRITE(numout,*) 'woodmass(ivma)',woodmass(ipts,ivma)
    DO iagec = nagec_pft(ivma),1,-1

       WRITE(numout,*) 'iagec',iagec
       !start from oldest age class and then move to younger age classes.
       ipft = start_index(ivma)+iagec-1

       WRITE(numout,*) 'ipft',ipft
       !  Check whether woodmass exceeds boundaries of
       !  the age class. For forest PFTs, bound_spa stores the upper boundary
       !  value .
       IF(ld_agec)THEN
          WRITE(numout,*) 'Checking to merge for: '
          WRITE(numout,*) 'ipft,iagec,ipts: ',ipft,iagec,ipts
          WRITE(numout,*) 'nagec_pft,woodmass,age_class_bound: ',nagec_pft(ivma),&
               woodmass(ipts,ipft),bound_spa(ipts,ipft)
       ENDIF
!! default
       !IF ( (iagec .EQ. nagec_pft(ivma)) .AND. &
       !     woodmass(ipts,ipft) .GT. bound_spa(ipts,ipft) ) THEN
!! yidi
       tot_biomass(ipts,ipft,:) = biomass(ipts,ipft,ileaf,:) + biomass(ipts,ipft,isapabove,:) + &
             &                    biomass(ipts,ipft,isapbelow,:) +&
             &                    biomass(ipts,ipft,iheartabove,:) +biomass(ipts,ipft,iheartbelow,:) + &
             &                    biomass(ipts,ipft,iroot,:) + biomass(ipts,ipft,ifruit,:) + &
             &                    biomass(ipts,ipft,icarbres,:)
       WRITE(numout,*) 'tot_biomass(ipts,ipft,icarbon)',tot_biomass(ipts,ipft,:)
       IF ( (iagec .EQ. nagec_pft(ivma)) .AND. &
            tot_biomass(ipts,ipft,icarbon) .GT. bound_spa(ipts,ipft) ) THEN
  
          ! If these conditions are satisfied our woodmass is
          ! very unrealist
          WRITE(numout,*) 'WARNING: age class indicator exceeds: ', &
               bound_spa(ipts,ipft) 
  
       !ELSEIF ( iagec .NE. nagec_pft(ivma) .AND. iagec .NE. nagec_pft(ivma)-1 .AND. &
       !     woodmass(ipts,ipft) .GT. bound_spa(ipts,ipft) ) THEN
       !ELSEIF ( iagec .NE. nagec_pft(ivma) .AND. iagec .NE. nagec_pft(ivma)-1 .AND. &
       !     tot_biomass(ipts,ipft,icarbon) .GT. bound_spa(ipts,ipft) ) THEN
       ELSEIF ( iagec .NE. nagec_pft(ivma) .AND. &
            tot_biomass(ipts,ipft,icarbon) .GT. bound_spa(ipts,ipft) ) THEN
!! yidi
          IF(ld_agec)THEN
             WRITE(numout,*) 'Merging biomass'
             WRITE(numout,*) 'ipts,ipft,iagec: ',ipts,ipft,iagec
             WRITE(numout,*) 'age_class_bound: ',bound_spa(ipts,ipft)
             WRITE(numout,*) 'woodmass: ',woodmass(ipts,ipft)

          ENDIF

          !! 2 Merge biomass
          !  Biomass of two age classes needs to be merged. The established
          !  vegetation is stored in ipft+1, the new vegetation is stored in
          !  ipft
          share_expanded = veget_max(ipts,ipft+1) / &
               ( veget_max(ipts,ipft+1) + veget_max(ipts,ipft) )
          ! We also need a scaling factor which includes the litter
          DO ilev=1,nlevs
             DO ilit=1,nlitt
                IF(litter(ipts,ilit,ipft,ilev,icarbon) .GE. min_stomate)THEN
                   litter_weight_expanded(ilit,ilev)=litter(ipts,ilit,ipft+1,ilev,icarbon) * veget_max(ipts,ipft+1)/ &
                        (litter(ipts,ilit,ipft+1,ilev,icarbon) * veget_max(ipts,ipft+1) + &
                        litter(ipts,ilit,ipft,ilev,icarbon) * veget_max(ipts,ipft))
                ELSE
                   litter_weight_expanded(ilit,ilev)=zero
                ENDIF
             END DO
          ENDDO

          


          IF (((ipft .GE. 2) .AND. (ipft .LE. 7)) .OR. (ipft .GE. 62)) THEN
            WRITE(numout,*) 'yd: check merge before,nvm=',ipft,',biomass(:,ipft,:,icarbon)',biomass(:,ipft,:,icarbon)
            WRITE(numout,*) 'yd: check merge before,nvm=',ipft+1,',biomass(:,ipft+1,:,icarbon)',biomass(:,ipft+1,:,icarbon)
            WRITE(numout,*) 'yd: check merge before,nvm=',ipft,',biomass(:,ipft,all,icarbon)',SUM(biomass(:,ipft,1:8,icarbon))
          ENDIF
          ! Merge the biomass and ind of the two age classes
          biomass(ipts,ipft+1,:,:) = share_expanded * biomass(ipts,ipft+1,:,:) + &
               (un - share_expanded) * biomass(ipts,ipft,:,:)
          ind(ipts,ipft+1) = share_expanded * ind(ipts,ipft+1) + &
               (un - share_expanded) * ind(ipts,ipft)
          
          !! 3 Empty the age class that was merged and update veget_max
          ind(ipts,ipft) = zero
          biomass(ipts,ipft,:,:) = zero
          veget_max(ipts,ipft+1) = veget_max(ipts,ipft+1) + veget_max(ipts,ipft)
          veget_max(ipts,ipft) = zero
!! yidi
          IF (ok_oilpalm) THEN
              IF (is_oilpalm(ipft)) THEN 
                 bm_FFB(ipts,ipft+1,:) = share_expanded * bm_FFB(ipts,ipft+1,:) + &
                     (un - share_expanded) * bm_FFB(ipts,ipft,:)
                 bm_phytomer(ipts,ipft+1,:) = share_expanded * bm_phytomer(ipts,ipft+1,:) + &
                     (un - share_expanded) * bm_phytomer(ipts,ipft,:)
                 FFBbm(ipts,ipft+1) = share_expanded * FFBbm(ipts,ipft+1) + &
                     (un - share_expanded) * FFBbm(ipts,ipft)
                 PHYbm(ipts,ipft+1) = share_expanded * PHYbm(ipts,ipft+1) + &
                     (un - share_expanded) * PHYbm(ipts,ipft)
                 bm_phytomer(ipts,ipft,:) = zero
                 bm_FFB(ipts,ipft,:) = zero
                 FFBbm(ipts,ipft) = zero
                 PHYbm(ipts,ipft) = zero
                 WRITE(numout,*) 'yd: check merge after,nvm=',ipft,',bmphytomer(:,j,1):', bm_phytomer(1,ipft,1) !! yidi output
                 WRITE(numout,*) 'yd: check merge after,nvm=',ipft,',bmFFB(:,j,1):', bm_FFB(1,ipft,1) !! yidi output
                 WRITE(numout,*) 'yd: check merge after,nvm=',ipft,', FFBbm(:,j):',FFBbm(1,ipft) !! yidi output
                 WRITE(numout,*) 'yd: check merge after,nvm=',ipft,', PHYbm(:,j):',PHYbm(1,ipft) !! yidi output
              ENDIF
          ENDIF
 
          IF (((ipft .GE. 2) .AND. (ipft .LE. 7)) .OR. (ipft .GE. 62)) THEN
            WRITE(numout,*) 'yd: check merge after,nvm=',ipft,',biomass(:,ipft,:,icarbon)',biomass(:,ipft,:,icarbon)
            WRITE(numout,*) 'yd: check merge after,nvm=',ipft+1,',biomass(:,ipft+1,:,icarbon)',biomass(:,ipft+1,:,icarbon)
            WRITE(numout,*) 'yd: check merge after,nvm=',ipft,',biomass(:,ipft,all,icarbon)',SUM(biomass(:,ipft,1:8,icarbon))
          ENDIF
!! yidi
          !! 4 Calculate the PFT characteristics of the merged PFT
          !  Take the weighted mean of the existing vegetation and the new 
          !  vegetation joining this PFT. 
          !  Note that co2_to_bm is in gC. m-2 dt-1 ,
          !  so we should also take the weighted mean (rather than sum if
          !  this where absolute values).

          lm_lastyearmax(ipts,ipft+1) = share_expanded * lm_lastyearmax(ipts,ipft+1) + &
               (un - share_expanded) * lm_lastyearmax(ipts,ipft)
          lm_lastyearmax(ipts,ipft) = zero
          age(ipts,ipft+1) = share_expanded * age(ipts,ipft+1) + &
               (un - share_expanded) * age(ipts,ipft)
          age(ipts,ipft) = zero

          !CHECK: more strictly this should be considered together with leaf mass
          leaf_frac(ipts,ipft+1,:) = share_expanded * leaf_frac(ipts,ipft+1,:) + &
               (un - share_expanded) * leaf_frac(ipts,ipft,:)
          leaf_frac(ipts,ipft,:) = zero
          leaf_age(ipts,ipft+1,:) = share_expanded * leaf_age(ipts,ipft+1,:) + &
               (un - share_expanded) * leaf_age(ipts,ipft,:)
          leaf_age(ipts,ipft,:) = zero
          co2_to_bm(ipts,ipft+1) = share_expanded * co2_to_bm(ipts,ipft+1) + &
               (un - share_expanded) * co2_to_bm(ipts,ipft)
          co2_to_bm(ipts,ipft) = zero

          ! Everywhere deals with the migration of vegetation. Copy the
          ! status of the most migrated vegetation for the whole PFT
          everywhere(ipts,ipft+1) = MAX(everywhere(ipts,ipft), everywhere(ipts,ipft+1))
          everywhere(ipts,ipft) = zero

          ! The new soil&litter pools are the weighted mean of the newly 
          ! established vegetation for that PFT and the soil&litter pools
          ! of the original vegetation that already exists in that PFT.
          ! Since it is not only the amount of vegetation present (veget_max) but also
          ! the amount of structural litter (litter) that is important, we have to 
          ! weight by both items here.
          DO ilev=1,nlevs
             lignin_struc(ipts,ipft+1,ilev) = litter_weight_expanded(istructural,ilev) * lignin_struc(ipts,ipft+1,ilev) + &
                  (un - litter_weight_expanded(istructural,ilev)) * lignin_struc(ipts,ipft,ilev) 
             lignin_struc(ipts,ipft,ilev) = zero
          ENDDO
          litter(ipts,:,ipft+1,:,:) = share_expanded * litter(ipts,:,ipft+1,:,:) + &
               (un - share_expanded) * litter(ipts,:,ipft,:,:)
          litter(ipts,:,ipft,:,:) = zero

          fuel_1hr(ipts,ipft+1,:,:) = share_expanded * fuel_1hr(ipts,ipft+1,:,:) + &
               (un - share_expanded) * fuel_1hr(ipts,ipft,:,:)
          fuel_1hr(ipts,ipft,:,:) = zero

          fuel_10hr(ipts,ipft+1,:,:) = share_expanded * fuel_10hr(ipts,ipft+1,:,:) + &
               (un - share_expanded) * fuel_10hr(ipts,ipft,:,:)
          fuel_10hr(ipts,ipft,:,:) = zero

          fuel_100hr(ipts,ipft+1,:,:) = share_expanded * fuel_100hr(ipts,ipft+1,:,:) + &
               (un - share_expanded) * fuel_100hr(ipts,ipft,:,:)
          fuel_100hr(ipts,ipft,:,:) = zero

          fuel_1000hr(ipts,ipft+1,:,:) = share_expanded * fuel_1000hr(ipts,ipft+1,:,:) + &
               (un - share_expanded) * fuel_1000hr(ipts,ipft,:,:)
          fuel_1000hr(ipts,ipft,:,:) = zero

          carbon(ipts,:,ipft+1) =  share_expanded * carbon(ipts,:,ipft+1) + &
               (un - share_expanded) * carbon(ipts,:,ipft)
          carbon(ipts,:,ipft) = zero 

          deepC_a(ipts,:,ipft+1) =  share_expanded * deepC_a(ipts,:,ipft+1) + &
               (un - share_expanded) * deepC_a(ipts,:,ipft)
          deepC_a(ipts,:,ipft) = zero 

          deepC_s(ipts,:,ipft+1) =  share_expanded * deepC_s(ipts,:,ipft+1) + &
               (un - share_expanded) * deepC_s(ipts,:,ipft)
          deepC_s(ipts,:,ipft) = zero 

          deepC_p(ipts,:,ipft+1) =  share_expanded * deepC_p(ipts,:,ipft+1) + &
               (un - share_expanded) * deepC_p(ipts,:,ipft)
          deepC_p(ipts,:,ipft) = zero 

          bm_to_litter(ipts,ipft+1,:,:) = share_expanded * bm_to_litter(ipts,ipft+1,:,:) + & 
               (un - share_expanded) * bm_to_litter(ipts,ipft,:,:)
          bm_to_litter(ipts,ipft,:,:) = zero

          ! Copy variables that depend on veget_max 
          when_growthinit(ipts,ipft+1) = share_expanded * when_growthinit(ipts,ipft+1) + &
               (un - share_expanded) * when_growthinit(ipts,ipft)
          when_growthinit(ipts,ipft) = zero
          gdd_from_growthinit(ipts,ipft+1) = share_expanded * &
               gdd_from_growthinit(ipts,ipft+1) + &
               (un - share_expanded) * gdd_from_growthinit(ipts,ipft)
          gdd_from_growthinit(ipts,ipft) = zero
          gdd_midwinter(ipts,ipft+1) = share_expanded * gdd_midwinter(ipts,ipft+1) + &
               (un - share_expanded) * gdd_midwinter(ipts,ipft)
          gdd_midwinter(ipts,ipft) = zero
          time_hum_min(ipts,ipft+1) = share_expanded * time_hum_min(ipts,ipft+1) + &
               (un - share_expanded) * time_hum_min(ipts,ipft)
          time_hum_min(ipts,ipft) = zero
          hum_min_dormance(ipts,ipft+1) = share_expanded * hum_min_dormance(ipts,ipft+1) + &
               (un - share_expanded) * hum_min_dormance(ipts,ipft)
          hum_min_dormance(ipts,ipft) = zero
          gdd_m5_dormance(ipts,ipft+1) = share_expanded * gdd_m5_dormance(ipts,ipft+1) + &
               (un - share_expanded) * gdd_m5_dormance(ipts,ipft)
          gdd_m5_dormance(ipts,ipft) = zero
          ncd_dormance(ipts,ipft+1) = share_expanded * ncd_dormance(ipts,ipft+1) + &
               (un - share_expanded) * ncd_dormance(ipts,ipft)
          ncd_dormance(ipts,ipft) = zero
          moiavail_month(ipts,ipft+1) = share_expanded * moiavail_month(ipts,ipft+1) + &
               (un - share_expanded) * moiavail_month(ipts,ipft)
          moiavail_month(ipts,ipft) = zero
          moiavail_week(ipts,ipft+1) = share_expanded * moiavail_week(ipts,ipft+1) + &
               (un - share_expanded) * moiavail_week(ipts,ipft)
          moiavail_week(ipts,ipft) = zero
          ngd_minus5(ipts,ipft+1) = share_expanded * ngd_minus5(ipts,ipft+1) + &
               (un - share_expanded) * ngd_minus5(ipts,ipft)
          ngd_minus5(ipts,ipft) = zero
    
          ! Copy remaining properties 
          PFTpresent(ipts,ipft+1) = PFTpresent(ipts,ipft)
          PFTpresent(ipts,ipft) = .FALSE.
          senescence(ipts,ipft+1) = senescence(ipts,ipft)
          senescence(ipts,ipft) = .FALSE.
          npp_longterm(ipts,ipft+1) = share_expanded * npp_longterm(ipts,ipft+1) + &
               (un - share_expanded) * npp_longterm(ipts,ipft)
          npp_longterm(ipts,ipft) = zero
          gpp_daily(ipts,ipft+1) = share_expanded * gpp_daily(ipts,ipft+1) + &
               (un - share_expanded) * gpp_daily(ipts,ipft)
          gpp_daily(ipts,ipft) = zero 
          gpp_week(ipts,ipft+1) = share_expanded * gpp_week(ipts,ipft+1) + &
               (un - share_expanded) * gpp_week(ipts,ipft)
          gpp_week(ipts,ipft) = zero
          resp_maint(ipts,ipft+1) = share_expanded * resp_maint(ipts,ipft+1) + &
               (un - share_expanded) * resp_maint(ipts,ipft) 
          resp_maint(ipts,ipft) = zero
          resp_growth(ipts,ipft+1) = share_expanded * resp_growth(ipts,ipft+1) + &
               (un - share_expanded) * resp_growth(ipts,ipft) 
          resp_growth(ipts,ipft) = zero
          npp_daily(ipts,ipft+1) = share_expanded * npp_daily(ipts,ipft+1) + &
               (un - share_expanded) * npp_daily(ipts,ipft) 
          npp_daily(ipts,ipft) = zero
          resp_hetero(ipts,ipft+1) = share_expanded * resp_hetero(ipts,ipft+1) + &
               (un - share_expanded) * resp_hetero(ipts,ipft) 
          resp_hetero(ipts,ipft) = zero

       ENDIF
    ENDDO
  ENDIF
  ! ! concerned MTC is grass/pasture/crop
  ! ELSE
  !   DO iagec = 1,nagec_pft(ivma),1

  !      ! As the soil C gets smaller when forest-generating crop gets older,
  !      ! we start from young age class and then move to older age classes.
  !      ! If the soil C of ipft is smaller than the threshold, then it should 
  !      ! go to the next age class.
  !      ipft = start_index(ivma)+iagec-1

  !      !  Check whether woodmass exceeds boundaries of
  !      !  the age class.
  !      IF(ld_agec)THEN
  !         WRITE(numout,*) 'Checking to merge for: '
  !         WRITE(numout,*) 'ipft,iagec,ipts: ',ipft,iagec,ipts
  !         WRITE(numout,*) 'nagec_pft,woodmass,age_class_bound: ',nagec_pft(ivma),&
  !              woodmass(ipts,ipft),bound_spa(ipts,ipft)
  !      ENDIF

  !      !IF ( (iagec .EQ. 1) .AND. &
  !      !     woodmass(ipts,ipft) .GT. bound_spa(ipts,ipft) ) THEN
  !      ! 
  !      !   ! If this is satisfied than we're having a quite large
  !      !   ! soil C in the newly initiated crop
  !      !   WRITE(numout,*) 'WARNING: age class indicator exceeds: ', &
  !      !        bound_spa(ipts,ipft) 
  ! 
  !      !ELSEIF ( (iagec .NE. nagec_pft(ivma)) .AND. &
  !      !     woodmass(ipts,ipft) .LT. bound_spa(ipts,ipft)) THEN

  !      ! If the soil C is smaller than the threshold and the concerned
  !      ! ipft is not the oldest age class, then it should move to the
  !      ! next (older) age class. So we have to set the soil C threshold
  !      ! for crop as:

  !      ! youngest:   0.9 of maximum end-spinup forest soil C
  !      ! 2nd young:  0.75 of maximum end-spniup forest soil C
  !      ! old:        0.55 of maximum end-spniup forest soil C
  !      ! oldest:     the oldest one should not be less than zero. 
  !      IF ( (iagec .NE. nagec_pft(ivma)) .AND. &
  !           woodmass(ipts,ipft) .LT. bound_spa(ipts,ipft) .AND. veget_max(ipts,ipft) .GT. min_stomate) THEN
  !         IF(ld_agec)THEN
  !            WRITE(numout,*) 'Merging biomass'
  !            WRITE(numout,*) 'ipts,ipft,iagec: ',ipts,ipft,iagec
  !            WRITE(numout,*) 'age_class_bound: ',bound_spa(ipts,ipft)
  !            WRITE(numout,*) 'woodmass: ',woodmass(ipts,ipft)

  !         ENDIF

  !         !! 2 Merge biomass
  !         !  Biomass of two age classes needs to be merged. The established
  !         !  vegetation is stored in ipft+1, the new vegetation is stored in
  !         !  ipft
  !         share_expanded = veget_max(ipts,ipft+1) / &
  !              ( veget_max(ipts,ipft+1) + veget_max(ipts,ipft) )
  !         ! We also need a scaling factor which includes the litter
  !         DO ilev=1,nlevs
  !            DO ilit=1,nlitt
  !               IF(litter(ipts,ilit,ipft,ilev,icarbon) .GE. min_stomate)THEN
  !                  litter_weight_expanded(ilit,ilev)=litter(ipts,ilit,ipft+1,ilev,icarbon) * veget_max(ipts,ipft+1)/ &
  !                       (litter(ipts,ilit,ipft+1,ilev,icarbon) * veget_max(ipts,ipft+1) + &
  !                       litter(ipts,ilit,ipft,ilev,icarbon) * veget_max(ipts,ipft))
  !               ELSE
  !                  litter_weight_expanded(ilit,ilev)=zero
  !               ENDIF
  !            END DO
  !         ENDDO

  !         ! Merge the biomass and ind of the two age classes
  !         biomass(ipts,ipft+1,:,:) = share_expanded * biomass(ipts,ipft+1,:,:) + &
  !              (un - share_expanded) * biomass(ipts,ipft,:,:)
  !         ind(ipts,ipft+1) = share_expanded * ind(ipts,ipft+1) + &
  !              (un - share_expanded) * ind(ipts,ipft)
  !         
  !         !! 3 Empty the age class that was merged and update veget_max
  !         ind(ipts,ipft) = zero
  !         biomass(ipts,ipft,:,:) = zero
  !         veget_max(ipts,ipft+1) = veget_max(ipts,ipft+1) + veget_max(ipts,ipft)
  !         veget_max(ipts,ipft) = zero
 
  !         !! 4 Calculate the PFT characteristics of the merged PFT
  !         !  Take the weighted mean of the existing vegetation and the new 
  !         !  vegetation joining this PFT. 
  !         !  Note that co2_to_bm is in gC. m-2 dt-1 ,
  !         !  so we should also take the weighted mean (rather than sum if
  !         !  this where absolute values).
  !         lm_lastyearmax(ipts,ipft+1) = share_expanded * lm_lastyearmax(ipts,ipft+1) + &
  !              (un - share_expanded) * lm_lastyearmax(ipts,ipft)
  !         lm_lastyearmax(ipts,ipft) = zero
  !         !age(ipts,ipft+1) = share_expanded * age(ipts,ipft+1) + &
  !         !     (un - share_expanded) * age(ipts,ipft)
  !         !age(ipts,ipft) = zero

  !         !CHECK: more strictly this should be considered together with leaf mass
  !         leaf_frac(ipts,ipft+1,:) = share_expanded * leaf_frac(ipts,ipft+1,:) + &
  !              (un - share_expanded) * leaf_frac(ipts,ipft,:)
  !         leaf_frac(ipts,ipft,:) = zero
  !         leaf_age(ipts,ipft+1,:) = share_expanded * leaf_age(ipts,ipft+1,:) + &
  !              (un - share_expanded) * leaf_age(ipts,ipft,:)
  !         leaf_age(ipts,ipft,:) = zero
  !         co2_to_bm(ipts,ipft+1) = share_expanded * co2_to_bm(ipts,ipft+1) + &
  !              (un - share_expanded) * co2_to_bm(ipts,ipft)
  !         co2_to_bm(ipts,ipft) = zero

  !         ! Everywhere deals with the migration of vegetation. Copy the
  !         ! status of the most migrated vegetation for the whole PFT
  !         everywhere(ipts,ipft+1) = MAX(everywhere(ipts,ipft), everywhere(ipts,ipft+1))
  !         everywhere(ipts,ipft) = zero

  !         ! The new soil&litter pools are the weighted mean of the newly 
  !         ! established vegetation for that PFT and the soil&litter pools
  !         ! of the original vegetation that already exists in that PFT.
  !         ! Since it is not only the amount of vegetation present (veget_max) but also
  !         ! the amount of structural litter (litter) that is important, we have to 
  !         ! weight by both items here.
  !         DO ilev=1,nlevs
  !            lignin_struc(ipts,ipft+1,ilev) = litter_weight_expanded(istructural,ilev) * lignin_struc(ipts,ipft+1,ilev) + &
  !                 (un - litter_weight_expanded(istructural,ilev)) * lignin_struc(ipts,ipft,ilev) 
  !            lignin_struc(ipts,ipft,ilev) = zero
  !         ENDDO
  !         litter(ipts,:,ipft+1,:,:) = share_expanded * litter(ipts,:,ipft+1,:,:) + &
  !              (un - share_expanded) * litter(ipts,:,ipft,:,:)
  !         litter(ipts,:,ipft,:,:) = zero

  !         fuel_1hr(ipts,ipft+1,:,:) = share_expanded * fuel_1hr(ipts,ipft+1,:,:) + &
  !              (un - share_expanded) * fuel_1hr(ipts,ipft,:,:)
  !         fuel_1hr(ipts,ipft,:,:) = zero

  !         fuel_10hr(ipts,ipft+1,:,:) = share_expanded * fuel_10hr(ipts,ipft+1,:,:) + &
  !              (un - share_expanded) * fuel_10hr(ipts,ipft,:,:)
  !         fuel_10hr(ipts,ipft,:,:) = zero

  !         fuel_100hr(ipts,ipft+1,:,:) = share_expanded * fuel_100hr(ipts,ipft+1,:,:) + &
  !              (un - share_expanded) * fuel_100hr(ipts,ipft,:,:)
  !         fuel_100hr(ipts,ipft,:,:) = zero

  !         fuel_1000hr(ipts,ipft+1,:,:) = share_expanded * fuel_1000hr(ipts,ipft+1,:,:) + &
  !              (un - share_expanded) * fuel_1000hr(ipts,ipft,:,:)
  !         fuel_1000hr(ipts,ipft,:,:) = zero

  !         carbon(ipts,:,ipft+1) =  share_expanded * carbon(ipts,:,ipft+1) + &
  !              (un - share_expanded) * carbon(ipts,:,ipft)
  !         carbon(ipts,:,ipft) = zero 

  !         deepC_a(ipts,:,ipft+1) =  share_expanded * deepC_a(ipts,:,ipft+1) + &
  !              (un - share_expanded) * deepC_a(ipts,:,ipft)
  !         deepC_a(ipts,:,ipft) = zero 

  !         deepC_s(ipts,:,ipft+1) =  share_expanded * deepC_s(ipts,:,ipft+1) + &
  !              (un - share_expanded) * deepC_s(ipts,:,ipft)
  !         deepC_s(ipts,:,ipft) = zero 

  !         deepC_p(ipts,:,ipft+1) =  share_expanded * deepC_p(ipts,:,ipft+1) + &
  !              (un - share_expanded) * deepC_p(ipts,:,ipft)
  !         deepC_p(ipts,:,ipft) = zero 

  !         bm_to_litter(ipts,ipft+1,:,:) = share_expanded * bm_to_litter(ipts,ipft+1,:,:) + & 
  !              (un - share_expanded) * bm_to_litter(ipts,ipft,:,:)
  !         bm_to_litter(ipts,ipft,:,:) = zero

  !         ! Copy variables that depend on veget_max 
  !         when_growthinit(ipts,ipft+1) = share_expanded * when_growthinit(ipts,ipft+1) + &
  !              (un - share_expanded) * when_growthinit(ipts,ipft)
  !         when_growthinit(ipts,ipft) = zero
  !         gdd_from_growthinit(ipts,ipft+1) = share_expanded * &
  !              gdd_from_growthinit(ipts,ipft+1) + &
  !              (un - share_expanded) * gdd_from_growthinit(ipts,ipft)
  !         gdd_from_growthinit(ipts,ipft) = zero
  !         gdd_midwinter(ipts,ipft+1) = share_expanded * gdd_midwinter(ipts,ipft+1) + &
  !              (un - share_expanded) * gdd_midwinter(ipts,ipft)
  !         gdd_midwinter(ipts,ipft) = zero
  !         time_hum_min(ipts,ipft+1) = share_expanded * time_hum_min(ipts,ipft+1) + &
  !              (un - share_expanded) * time_hum_min(ipts,ipft)
  !         time_hum_min(ipts,ipft) = zero
  !         hum_min_dormance(ipts,ipft+1) = share_expanded * hum_min_dormance(ipts,ipft+1) + &
  !              (un - share_expanded) * hum_min_dormance(ipts,ipft)
  !         hum_min_dormance(ipts,ipft) = zero
  !         gdd_m5_dormance(ipts,ipft+1) = share_expanded * gdd_m5_dormance(ipts,ipft+1) + &
  !              (un - share_expanded) * gdd_m5_dormance(ipts,ipft)
  !         gdd_m5_dormance(ipts,ipft) = zero
  !         ncd_dormance(ipts,ipft+1) = share_expanded * ncd_dormance(ipts,ipft+1) + &
  !              (un - share_expanded) * ncd_dormance(ipts,ipft)
  !         ncd_dormance(ipts,ipft) = zero
  !         moiavail_month(ipts,ipft+1) = share_expanded * moiavail_month(ipts,ipft+1) + &
  !              (un - share_expanded) * moiavail_month(ipts,ipft)
  !         moiavail_month(ipts,ipft) = zero
  !         moiavail_week(ipts,ipft+1) = share_expanded * moiavail_week(ipts,ipft+1) + &
  !              (un - share_expanded) * moiavail_week(ipts,ipft)
  !         moiavail_week(ipts,ipft) = zero
  !         ngd_minus5(ipts,ipft+1) = share_expanded * ngd_minus5(ipts,ipft+1) + &
  !              (un - share_expanded) * ngd_minus5(ipts,ipft)
  !         ngd_minus5(ipts,ipft) = zero
  !   
  !         ! Copy remaining properties 
  !         PFTpresent(ipts,ipft+1) = PFTpresent(ipts,ipft)
  !         PFTpresent(ipts,ipft) = .FALSE.
  !         senescence(ipts,ipft+1) = senescence(ipts,ipft)
  !         senescence(ipts,ipft) = .FALSE.
  !         npp_longterm(ipts,ipft+1) = share_expanded * npp_longterm(ipts,ipft+1) + &
  !              (un - share_expanded) * npp_longterm(ipts,ipft)
  !         npp_longterm(ipts,ipft) = zero
  !         gpp_daily(ipts,ipft+1) = share_expanded * gpp_daily(ipts,ipft+1) + &
  !              (un - share_expanded) * gpp_daily(ipts,ipft)
  !         gpp_daily(ipts,ipft) = zero 
  !         gpp_week(ipts,ipft+1) = share_expanded * gpp_week(ipts,ipft+1) + &
  !              (un - share_expanded) * gpp_week(ipts,ipft)
  !         gpp_week(ipts,ipft) = zero
  !         resp_maint(ipts,ipft+1) = share_expanded * resp_maint(ipts,ipft+1) + &
  !              (un - share_expanded) * resp_maint(ipts,ipft) 
  !         resp_maint(ipts,ipft) = zero
  !         resp_growth(ipts,ipft+1) = share_expanded * resp_growth(ipts,ipft+1) + &
  !              (un - share_expanded) * resp_growth(ipts,ipft) 
  !         resp_growth(ipts,ipft) = zero
  !         npp_daily(ipts,ipft+1) = share_expanded * npp_daily(ipts,ipft+1) + &
  !              (un - share_expanded) * npp_daily(ipts,ipft) 
  !         npp_daily(ipts,ipft) = zero
  !         resp_hetero(ipts,ipft+1) = share_expanded * resp_hetero(ipts,ipft+1) + &
  !              (un - share_expanded) * resp_hetero(ipts,ipft) 
  !         resp_hetero(ipts,ipft) = zero

  !      ENDIF
  !   ENDDO

  ! ENDIF

  END SUBROUTINE check_merge_same_MTC



! ================================================================================================================================
!! SUBROUTINE   gross_lcchange
!!
!>\BRIEF       : Apply gross land cover change.
!!
!>\DESCRIPTION  
!_ ================================================================================================================================
  SUBROUTINE glcc_MulAgeC (npts, dt_days, newvegfrac,  &
               glccSecondShift,glccPrimaryShift,glccNetLCC,&
               def_fuel_1hr_remain, def_fuel_10hr_remain,        &
               def_fuel_100hr_remain, def_fuel_1000hr_remain,    &
               deforest_litter_remain, deforest_biomass_remain,  &
               convflux,                    &
               glccReal, IncreDeficit, glcc_pft, glcc_pftmtc,          &
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
               bm_phytomer,bm_FFB,PHYbm, FFBbm,                        & !! yidi
               fuel_1hr, fuel_10hr, fuel_100hr, fuel_1000hr)
  
    IMPLICIT NONE

    !! 0.1 Input variables

    INTEGER, INTENT(in)                                  :: npts             !! Domain size - number of pixels (unitless)
    REAL(r_std), INTENT(in)                              :: dt_days          !! Time step of vegetation dynamics for stomate
    REAL(r_std), DIMENSION (npts,12),INTENT(in)          :: glccSecondShift  !! the land-cover-change (LCC) matrix in case a gross LCC is 
                                                                             !! used.
    REAL(r_std), DIMENSION (npts,12),INTENT(in)          :: glccPrimaryShift !! the land-cover-change (LCC) matrix in case a gross LCC is 
                                                                             !! used.
    REAL(r_std), DIMENSION (npts,12),INTENT(in)          :: glccNetLCC       !! the land-cover-change (LCC) matrix in case a gross LCC is 
                                                                             !! used.
    REAL(r_std), DIMENSION (npts,nvmap),INTENT(in)       :: newvegfrac       !! veget max fraction matrix to guide the allocation of newly
                                                                             !! created lands of a given vegetation type.

    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements), INTENT(in)       :: def_fuel_1hr_remain
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements), INTENT(in)       :: def_fuel_10hr_remain
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements), INTENT(in)       :: def_fuel_100hr_remain
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nelements), INTENT(in)       :: def_fuel_1000hr_remain
    REAL(r_std), DIMENSION(npts,nlitt,nvm,nlevs,nelements), INTENT(in) :: deforest_litter_remain   !! Vegetmax-weighted remaining litter on the ground for 
                                                                                                      !! deforestation region.
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(in)      :: deforest_biomass_remain  !! Vegetmax-weighted remaining biomass on the ground for 
                                                                                                      !! deforestation region.


    !! 0.2 Output variables
    REAL(r_std), DIMENSION(npts,nwp), INTENT(inout)         :: convflux         !! release during first year following land cover
                                                                              !! change
    REAL(r_std), DIMENSION(npts,12), INTENT(inout)        :: glccReal         !! The "real" glcc matrix that we apply in the model
                                                                              !! after considering the consistency between presribed
                                                                              !! glcc matrix and existing vegetation fractions.
    REAL(r_std), DIMENSION(npts,12), INTENT(inout)        :: IncreDeficit     !! Originally "Increment" deficits, negative values mean that 
                                                                              !! there are not enough fractions in the source PFTs
                                                                              !! /vegetations to target PFTs/vegetations. I.e., these
                                                                              !! fraction transfers are presribed in LCC matrix but
                                                                              !! not realized. Now the glccDeficit for all land cover changes
                                                                              !! except forestry harvest.
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)       :: glcc_pft         !! Loss of fraction in each PFT
    REAL(r_std), DIMENSION(npts,nvm,nvmap), INTENT(inout) :: glcc_pftmtc      !! a temporary variable to hold the fractions each PFT is going to lose
                                                                              !! i.e., the contribution of each PFT to the youngest age-class of MTC

    !! 0.3 Modified variables
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)       :: veget_max        !! "maximal" coverage fraction of a PFT (LAI ->
                                                                              !! infinity) on ground (unitless) 
    REAL(r_std), DIMENSION(npts,0:10,nwp), INTENT(inout)  :: prod10           !! products remaining in the 10 year-turnover
                                                                              !! pool after the annual release for each 
                                                                              !! compartment (10 + 1 : input from year of land
                                                                              !! cover change)
    REAL(r_std), DIMENSION(npts,0:100,nwp), INTENT(inout) :: prod100              !! products remaining in the 100 year-turnover
                                                                              !! pool after the annual release for each 
                                                                              !! compartment (100 + 1 : input from year of land
                                                                              !! cover change)
    LOGICAL, DIMENSION(:,:), INTENT(inout)                :: PFTpresent       !! Tab indicating which PFTs are present in 
                                                                              !! each pixel
    LOGICAL, DIMENSION(:,:), INTENT(inout)                :: senescence       !! Flag for setting senescence stage (only 
                                                                              !! for deciduous trees)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: moiavail_month   !! "Monthly" moisture availability (0 to 1, 
                                                                              !! unitless) 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: moiavail_week    !! "Weekly" moisture availability 
                                                                              !! (0 to 1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: gpp_week         !! Mean weekly gross primary productivity 
                                                                              !! @tex $(gC m^{-2} day^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: ngd_minus5       !! Number of growing days (days), threshold 
                                                                              !! -5 deg C (for phenology)   
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: resp_maint       !! Maintenance respiration  
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: resp_growth      !! Growth respiration  
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: resp_hetero      !! Heterotrophic respiration  
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: npp_daily        !! Net primary productivity 
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: when_growthinit  !! How many days ago was the beginning of 
                                                                              !! the growing season (days)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: npp_longterm     !! "Long term" mean yearly primary productivity 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: ind              !! Number of individuals at the stand level
                                                                              !! @tex $(m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: lm_lastyearmax   !! last year's maximum leaf mass for each PFT 
                                                                              !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: everywhere       !! is the PFT everywhere in the grid box or 
                                                                              !! very localized (after its introduction) (?)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: age              !! mean age (years)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: co2_to_bm        !! CO2 taken from the atmosphere to get C to create  
                                                                              !! the seedlings @tex (gC.m^{-2}dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: gpp_daily        !! Daily gross primary productivity  
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: co2_fire         !! Fire carbon emissions
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: time_hum_min     !! Time elapsed since strongest moisture 
                                                                              !! availability (days) 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: gdd_midwinter    !! Growing degree days (K), since midwinter 
                                                                              !! (for phenology) - this is written to the
                                                                              !!  history files 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: gdd_from_growthinit !! growing degree days, since growthinit 
                                                                              !! for crops
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: gdd_m5_dormance  !! Growing degree days (K), threshold -5 deg 
                                                                              !! C (for phenology)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: ncd_dormance     !! Number of chilling days (days), since 
                                                                              !! leaves were lost (for phenology) 
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)          :: lignin_struc     !! ratio Lignine/Carbon in structural litter,
                                                                              !! above and below ground
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)          :: carbon           !! carbon pool: active, slow, or passive 
                                                                              !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)          :: deepC_a          !! Permafrost soil carbon (g/m**3) active
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)          :: deepC_s          !! Permafrost soil carbon (g/m**3) slow
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)          :: deepC_p          !! Permafrost soil carbon (g/m**3) passive
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)          :: leaf_frac        !! fraction of leaves in leaf age class (unitless;0-1)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)          :: leaf_age         !! Leaf age (days)
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)        :: bm_to_litter     !! Transfer of biomass to litter 
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)        :: biomass          !! Stand level biomass @tex $(gC.m^{-2})$ @endtex
!! yidi
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: FFBbm               !! FFB mass, from sapabove
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: PHYbm               !! PHYTOMER mass, from sapabove
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: bm_phytomer         !! Each PHYTOMER mass, from sapabove
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: bm_FFB              !! Fruit mass for Each PHYTOMER
!! yidi
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)      :: litter           !! metabolic and structural litter, above and 
                                                                              !! below ground @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)        :: fuel_1hr
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)        :: fuel_10hr
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)        :: fuel_100hr
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)        :: fuel_1000hr

    !! 0.4 Local variables
    REAL(r_std), DIMENSION(nvmap,nparts,nelements)        :: bm_to_litter_pro !! conversion of biomass to litter 
                                                                              !! @tex ($gC m^{-2} day^{-1}$) @endtex
    REAL(r_std), DIMENSION(nvmap,nparts,nelements)        :: biomass_pro      !! biomass @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(nvmap)                         :: veget_max_pro    !! "maximal" coverage fraction of a PFT (LAI ->
                                                                              !! infinity) on ground (unitless) 
    REAL(r_std), DIMENSION(nvmap,ncarb)                   :: carbon_pro       !! carbon pool: active, slow, or passive 
                                                                              !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(nvmap,ndeep)                   :: deepC_a_pro      !! Permafrost carbon pool: active, slow, or passive 
                                                                              !! @tex ($gC m^{-3}$) @endtex
    REAL(r_std), DIMENSION(nvmap,ndeep)                   :: deepC_s_pro      !! Permafrost carbon pool: active, slow, or passive 
                                                                              !! @tex ($gC m^{-3}$) @endtex
    REAL(r_std), DIMENSION(nvmap,ndeep)                   :: deepC_p_pro      !! Permafrost carbon pool: active, slow, or passive 
                                                                              !! @tex ($gC m^{-3}$) @endtex
    REAL(r_std), DIMENSION(nvmap,nlitt,nlevs,nelements)   :: litter_pro       !! metabolic and structural litter, above and 
                                                                              !! below ground @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(nvmap,nlitt,nelements)         :: fuel_1hr_pro
    REAL(r_std), DIMENSION(nvmap,nlitt,nelements)         :: fuel_10hr_pro
    REAL(r_std), DIMENSION(nvmap,nlitt,nelements)         :: fuel_100hr_pro
    REAL(r_std), DIMENSION(nvmap,nlitt,nelements)         :: fuel_1000hr_pro
    REAL(r_std), DIMENSION(nvmap,nlevs)                   :: lignin_struc_pro !! ratio Lignine/Carbon in structural litter
                                                                              !! above and below ground
    REAL(r_std), DIMENSION(nvmap,nleafages)               :: leaf_frac_pro    !! fraction of leaves in leaf age class 
    REAL(r_std), DIMENSION(nvmap,nleafages)               :: leaf_age_pro     !! fraction of leaves in leaf age class 
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

    WRITE(numout,*) 'Entering glcc_MulAgeC'
    glccReal_tmp(:,:) = zero

    CALL glcc_MulAgeC_firstday(npts,veget_max,newvegfrac,   &
                          glccSecondShift,glccPrimaryShift,glccNetLCC,&
                          glccReal,glcc_pft,glcc_pftmtc,IncreDeficit)

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

        ! here we set (glcc_mtc(ipts,ivma) GT. min_stomate) as a condition,
        ! this is necessary because later on in the subroutine of 
        ! `add_incoming_proxy_pft` we have to merge the newly established
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

          ! Here we substract the outgoing fraction from the source PFT.
          ! If a too small fraction remains in this source PFT, then it is
          ! considered exhausted and we empty it. The subroutine `empty_pft`
          ! might be combined with `collect_legacy_pft` later.
          DO ivm = 1,nvm
            ! In the above we limit the collection of legacy pools to
            ! (glcc_mtc(ipts,ivma) GT. min_stomate). Here we tentatively use
            ! a lower threshold of `min_stomate*0.1`.
            IF( glcc_pftmtc(ipts,ivm,ivma)>min_stomate*0.1 ) THEN
              ! this is the key line to implement reduction of fraction of source
              ! PFT.
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

      ! We establish new youngest proxy and add it to the 
      ! existing youngest-age PFT. 
      DO ivma = 2,nvmap
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

        ! we take as a priority from exsiting PFTs of the same meta-class
        ! the sapling biomass needed to initialize the youngest-age-class
        ! PFT, to avoid a too much high amount of CO2 dragged down from
        ! the air.
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
        
      ENDDO !(DO ivma=1,nvmap)

    ENDDO !(DO ipts=1,npts)

    ! Write out history files
    CALL histwrite_p (hist_id_stomate, 'glcc_pft', itime, &
         glcc_pft, npts*nvm, horipft_index)
    CALL xios_orchidee_send_field ('glcc_pft', glcc_pft) ! kjpindex,nvm

    glccReal_tmp(:,1:12) = glccReal
    CALL histwrite_p (hist_id_stomate, 'glccReal', itime, &
         glccReal_tmp, npts*nvm, horipft_index)
    CALL xios_orchidee_send_field ('glccReal', glccReal_tmp) ! kjpindex,nvm

    glccReal_tmp(:,:) = zero
    glccReal_tmp(:,1:12) = IncreDeficit
    CALL histwrite_p (hist_id_stomate, 'IncreDeficit', itime, &
         glccReal_tmp, npts*nvm, horipft_index)
    CALL xios_orchidee_send_field ('IncreDeficit', glccReal_tmp) ! kjpindex,nvm

    DO ivma = 1, nvmap
      WRITE(part_str,'(I2)') ivma
      IF (ivma < 10) part_str(1:1) = '0'
      CALL histwrite_p (hist_id_stomate, 'glcc_pftmtc_'//part_str(1:LEN_TRIM(part_str)), &
           itime, glcc_pftmtc(:,:,ivma), npts*nvm, horipft_index)
    ENDDO
    CALL xios_orchidee_send_field ('glcc_pftmtc', glcc_pftmtc) ! kjpindex,nvm,nvmap
  
    WRITE(numout,*) 'End glcc_MulAgeC'
  END SUBROUTINE glcc_MulAgeC


! ================================================================================================================================
!! SUBROUTINE   : glcc_MulAgeC_firstday
!!
!>\BRIEF        : When necessary, adjust input glcc matrix, and allocate it
!!                into different contributing age classes and receiving 
!!                youngest age classes.
!! \n
!_ ================================================================================================================================

  ! Note: it has this name because this subroutine will also be called
  ! the first day of each year to precalculate the forest loss for the
  ! deforestation fire module.
  SUBROUTINE glcc_MulAgeC_firstday(npts,veget_max_org,newvegfrac, &
                          glccSecondShift,glccPrimaryShift,glccNetLCC,&
                          glccReal,glcc_pft,glcc_pftmtc,IncreDeficit)

    IMPLICIT NONE

    !! 0.1 Input variables

    INTEGER, INTENT(in)                                    :: npts             !! Domain size - number of pixels (unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)           :: veget_max_org    !! "maximal" coverage fraction of a PFT on the ground
                                                                               !! May sum to
                                                                               !! less than unity if the pixel has
                                                                               !! nobio area. (unitless, 0-1)
    REAL(r_std), DIMENSION(npts,nvmap), INTENT(in)         :: newvegfrac       !! used to guid the allocation of new PFTs.
                                                                               !! 
    REAL(r_std), DIMENSION (npts,12),INTENT(in)            :: glccSecondShift  !! the land-cover-change (LCC) matrix in case a gross LCC is 
                                                                               !! used.
    REAL(r_std), DIMENSION (npts,12),INTENT(in)            :: glccPrimaryShift !! the land-cover-change (LCC) matrix in case a gross LCC is 
                                                                               !! used.
    REAL(r_std), DIMENSION (npts,12),INTENT(in)            :: glccNetLCC       !! the land-cover-change (LCC) matrix in case a gross LCC is 
                                                                               !! used.

    !! 0.2 Output variables
    REAL(r_std), DIMENSION(npts,nvm,nvmap), INTENT(inout)  :: glcc_pftmtc      !! a temporary variable to hold the fractions each PFT is going to lose
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)        :: glcc_pft         !! Loss of fraction in each PFT

    REAL(r_std), DIMENSION(npts,12), INTENT(inout)         :: glccReal         !! Originally the "real" glcc matrix that we apply in the model
                                                                               !! after considering the consistency between presribed
                                                                               !! glcc matrix and existing vegetation fractions. Now the glcc
                                                                               !! by summing SecShift,NetLCC and PriShift
    REAL(r_std), DIMENSION(npts,12), INTENT(inout)         :: IncreDeficit     !! Originally "Increment" deficits, negative values mean that 
                                                                               !! there are not enough fractions in the source PFTs
                                                                               !! /vegetations to target PFTs/vegetations. I.e., these
                                                                               !! fraction transfers are presribed in LCC matrix but
                                                                               !! not realized. Now the glccDeficit for all land cover changes
                                                                               !! except forestry harvest.

    !! 0.3 Modified variables
    
    !! 0.4 Local variables
    REAL(r_std), DIMENSION(npts,nvmap)              :: veget_mtc           !! "maximal" coverage fraction of a PFT on the ground
    REAL(r_std), DIMENSION(npts,nvmap)              :: veget_mtc_begin     !! "maximal" coverage fraction of a PFT on the ground
    REAL(r_std), DIMENSION(npts,nagec_tree)         :: vegagec_tree        !! fraction of tree age-class groups, in sequence of old->young
    REAL(r_std), DIMENSION(npts,nagec_herb)         :: vegagec_grass       !! fraction of grass age-class groups, in sequence of old->young 
    REAL(r_std), DIMENSION(npts,nagec_herb)         :: vegagec_pasture     !! fraction of pasture age-class groups, in sequence of old->young
    REAL(r_std), DIMENSION(npts,nagec_herb)         :: vegagec_crop        !! fraction of crop age-class groups, in sequence of old->young

    
    REAL(r_std), DIMENSION(npts,4)                  :: veget_4veg          !! "maximal" coverage fraction of a PFT on the ground
    REAL(r_std), DIMENSION(npts)                    :: veget_tree          !! "maximal" coverage fraction of a PFT on the ground
    REAL(r_std), DIMENSION(npts)                    :: veget_grass         !! "maximal" coverage fraction of a PFT on the ground
    REAL(r_std), DIMENSION(npts)                    :: veget_pasture       !! "maximal" coverage fraction of a PFT on the ground
    REAL(r_std), DIMENSION(npts)                    :: veget_crop          !! "maximal" coverage fraction of a PFT on the ground

    REAL(r_std), DIMENSION(npts,nvm)                :: veget_max           !! "maximal" coverage fraction of a PFT on the ground
    REAL(r_std), DIMENSION(npts,nvm)                :: veget_max_tmp       !! "maximal" coverage fraction of a PFT on the ground
    REAL(r_std), DIMENSION(npts,nvm)                :: veget_max_old       !! "maximal" coverage fraction of a PFT on the ground
    REAL(r_std), DIMENSION(npts,nvm)                :: glcc_pft_tmp        !! Loss of fraction in each PFT

    REAL(r_std), DIMENSION(npts,nvm,nvmap)          :: glcc_pftmtc_SecShift   !! a temporary variable to hold the fractions each PFT is going to lose
    REAL(r_std), DIMENSION(npts,nvm,nvmap)          :: glcc_pftmtc_NetPri     !! a temporary variable to hold the fractions each PFT is going to lose

    REAL(r_std), DIMENSION(npts,12)                 :: glccRealSecShift       !! real matrix applied for secondary shifting cultivation.
    REAL(r_std), DIMENSION(npts,12)                 :: glccRealNetPriShift    !! real matrix applied for NetLCC+primary shifting cultivation.    

    REAL(r_std), DIMENSION(npts,12)                 :: glccDefSecShift        !! deficit for the glccSecondShift
    REAL(r_std), DIMENSION(npts,12)                 :: glccDefNetPriShift     !! deficit for the glccNetLCC + glccPriShift

    REAL(r_std), DIMENSION(npts,nvm)                :: glccReal_tmp           !! A temporary variable to hold glccReal

    LOGICAL, SAVE  :: glcc_MulAgeC_firstday_done = .FALSE.

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

    INTEGER :: ivma

    REAL(r_std), DIMENSION(npts,12)         :: glccRemain      
    REAL(r_std), DIMENSION(npts,12)         :: glccSecondShift_remain      

    INTEGER :: ipts,IndStart_f,IndEnd_f
    CHARACTER(LEN=10)      :: part_str                               !! string suffix indicating an index

    !Some more local configurations
    LOGICAL                                 :: allow_youngest_forest_SecShift = .TRUE.
    

    WRITE(numout,*) 'Entering glcc_MulAgeC_firstday'
    ! check for equal bi-directional transition in glccSecondShift
    DO ipts = 1,npts
      IF (ABS(glccSecondShift(ipts,f2g)-glccSecondShift(ipts,g2f)) .GE. min_stomate*100) THEN
        WRITE(numout,*) 'Error in input transition matrix for shifting cultivation'
        WRITE(numout,*) 'transition between f2g and g2f is not equal !'
        WRITE(numout,*) 'Grid cell number: ', ipts
        STOP
      ENDIF

      IF (ABS(glccSecondShift(ipts,f2c)-glccSecondShift(ipts,c2f)) .GE. min_stomate*100) THEN
        WRITE(numout,*) 'Error in input transition matrix for shifting cultivation'
        WRITE(numout,*) 'transition between f2c and c2f is not equal !'
        WRITE(numout,*) 'Grid cell number: ', ipts
        STOP
      ENDIF

      IF (ABS(glccSecondShift(ipts,f2p)-glccSecondShift(ipts,p2f)) .GE. min_stomate*100) THEN
        WRITE(numout,*) 'Error in input transition matrix for shifting cultivation'
        WRITE(numout,*) 'transition between f2p and p2f is not equal !'
        WRITE(numout,*) 'Grid cell number: ', ipts
        STOP
      ENDIF

      IF (ABS(glccSecondShift(ipts,g2p)-glccSecondShift(ipts,p2g)) .GE. min_stomate*100) THEN
        WRITE(numout,*) 'Error in input transition matrix for shifting cultivation'
        WRITE(numout,*) 'transition between g2p and p2g is not equal !'
        WRITE(numout,*) 'Grid cell number: ', ipts
        STOP
      ENDIF

      IF (ABS(glccSecondShift(ipts,g2c)-glccSecondShift(ipts,c2g)) .GE. min_stomate*100) THEN
        WRITE(numout,*) 'Error in input transition matrix for shifting cultivation'
        WRITE(numout,*) 'transition between g2c and c2g is not equal !'
        WRITE(numout,*) 'Grid cell number: ', ipts
        STOP
      ENDIF

      IF (ABS(glccSecondShift(ipts,p2c)-glccSecondShift(ipts,c2p)) .GE. min_stomate*100) THEN
        WRITE(numout,*) 'Error in input transition matrix for shifting cultivation'
        WRITE(numout,*) 'transition between p2c and c2p is not equal !'
        WRITE(numout,*) 'Grid cell number: ', ipts
        STOP
      ENDIF
    ENDDO

    ! check for equal bi-directional transition in glccPrimaryShift
    DO ipts = 1,npts
      IF (ABS(glccPrimaryShift(ipts,f2g)-glccPrimaryShift(ipts,g2f)) .GE. min_stomate*100) THEN
        WRITE(numout,*) 'Error in input transition matrix for shifting cultivation'
        WRITE(numout,*) 'transition between f2g and g2f is not equal !'
        WRITE(numout,*) 'Grid cell number: ', ipts
        STOP
      ENDIF

      IF (ABS(glccPrimaryShift(ipts,f2c)-glccPrimaryShift(ipts,c2f)) .GE. min_stomate*100) THEN
        WRITE(numout,*) 'Error in input transition matrix for shifting cultivation'
        WRITE(numout,*) 'transition between f2c and c2f is not equal !'
        WRITE(numout,*) 'Grid cell number: ', ipts
        STOP
      ENDIF

      IF (ABS(glccPrimaryShift(ipts,f2p)-glccPrimaryShift(ipts,p2f)) .GE. min_stomate*100) THEN
        WRITE(numout,*) 'Error in input transition matrix for shifting cultivation'
        WRITE(numout,*) 'transition between f2p and p2f is not equal !'
        WRITE(numout,*) 'Grid cell number: ', ipts
        STOP
      ENDIF

      IF (ABS(glccPrimaryShift(ipts,g2p)-glccPrimaryShift(ipts,p2g)) .GE. min_stomate*100) THEN
        WRITE(numout,*) 'Error in input transition matrix for shifting cultivation'
        WRITE(numout,*) 'transition between g2p and p2g is not equal !'
        WRITE(numout,*) 'Grid cell number: ', ipts
        STOP
      ENDIF

      IF (ABS(glccPrimaryShift(ipts,g2c)-glccPrimaryShift(ipts,c2g)) .GE. min_stomate*100) THEN
        WRITE(numout,*) 'Error in input transition matrix for shifting cultivation'
        WRITE(numout,*) 'transition between g2c and c2g is not equal !'
        WRITE(numout,*) 'Grid cell number: ', ipts
        STOP
      ENDIF

      IF (ABS(glccPrimaryShift(ipts,p2c)-glccPrimaryShift(ipts,c2p)) .GE. min_stomate*100) THEN
        WRITE(numout,*) 'Error in input transition matrix for shifting cultivation'
        WRITE(numout,*) 'transition between p2c and c2p is not equal !'
        WRITE(numout,*) 'Grid cell number: ', ipts
        STOP
      ENDIF
    ENDDO

    ! Initialization
    glccReal = zero
    glcc_pftmtc = zero
    glcc_pft = zero
    glcc_pft_tmp = zero

    !!! ** Land cover change processes start here ** !!!
    ! we make copies of original input veget_max (which is veget_max_org
    ! in the subroutine parameter list).
    ! veget_max will be modified through different operations in order to 
    ! check for various purposes, e.g., whether input glcc matrix 
    ! is compatible with existing veget_max and how to allocate it etc.
    ! veget_max_old will not be modified
    veget_max(:,:) = veget_max_org(:,:)
    veget_max_old(:,:) = veget_max_org(:,:)

    !! 3. Treat secondary-agriculture shifting cultivation transition matrix.
    !! The primary-agriculture shifting cultivation will be treated together
    !! with the netLCC transitions, with the conversion sequence of oldest->
    !! youngest is applied.

    ! When we prepare the driving data, secondary-agriculture shifting cultivation
    ! is intended to include the "constant transitions" over time. Ideally, we
    ! should start applying this secondary-agriculture shifting cultivation with
    ! the "secondary forest" in the model. Here we tentatively start with the 3rd
    ! youngest age class and move to the 2ne youngest age class. But if the prescribed
    ! transition fraction is not met, we then move further to 4th youngest age class
    ! and then move to the oldest age class sequentially.

    ! Note for the first call, we have to pass veget_mtc_begin instead of veget_mtc
    ! in order to keep the original veget_mtc before any convserion is made. The
    ! veget_mtc is used to in type_conversion to guide the allocation of newly
    ! created fraction of a certain mtc to its componenet youngest PFTs.
    CALL calc_cover(npts,veget_max,veget_mtc_begin,vegagec_tree,vegagec_grass, &
           vegagec_pasture,vegagec_crop)
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
      write(numout,*) 'glcc_MulAgeC: Age class index cannot be negative or zero!'
      STOP
    ENDIF

    DO ipts=1,npts
      !f2c
      CALL type_conversion(ipts,f2c,glccSecondShift,veget_mtc,newvegfrac,       &
                       indold_tree,indagec_tree,indagec_crop,num_crop_mulagec,     &
                       IndEnd_f,nagec_herb,                    &
                       vegagec_tree,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp, &
                       glccRemain, &
                       .TRUE., iagec_start=IndStart_f)
      !f2p
      CALL type_conversion(ipts,f2p,glccSecondShift,veget_mtc,newvegfrac,       &
                       indold_tree,indagec_tree,indagec_pasture,num_pasture_mulagec,     &
                       IndEnd_f,nagec_herb,                    &
                       vegagec_tree,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp, &
                       glccRemain, &
                       .TRUE., iagec_start=IndStart_f)
      !f2g
      CALL type_conversion(ipts,f2g,glccSecondShift,veget_mtc,newvegfrac,       &
                       indold_tree,indagec_tree,indagec_grass,num_grass_mulagec,     &
                       IndEnd_f,nagec_herb,                    &
                       vegagec_tree,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain, &
                       .TRUE., iagec_start=IndStart_f)
      !g2c
      CALL type_conversion(ipts,g2c,glccSecondShift,veget_mtc,newvegfrac,       &
                       indold_grass,indagec_grass,indagec_crop,num_crop_mulagec,     &
                       nagec_herb,nagec_herb,                    &
                       vegagec_grass,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain, &
                       .TRUE.)
      !g2p
      CALL type_conversion(ipts,g2p,glccSecondShift,veget_mtc,newvegfrac,       &
                       indold_grass,indagec_grass,indagec_pasture,num_pasture_mulagec,     &
                       nagec_herb,nagec_herb,                    &
                       vegagec_grass,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain, &
                       .TRUE.)
      !g2f
      CALL type_conversion(ipts,g2f,glccSecondShift,veget_mtc,newvegfrac,       &
                       indold_grass,indagec_grass,indagec_tree,num_tree_mulagec,     &
                       nagec_herb,nagec_tree,                    &
                       vegagec_grass,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain, &
                       .TRUE.)
      !p2c
      CALL type_conversion(ipts,p2c,glccSecondShift,veget_mtc,newvegfrac,       &
                       indold_pasture,indagec_pasture,indagec_crop,num_crop_mulagec,     &
                       nagec_herb,nagec_herb,                    &
                       vegagec_pasture,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain, &
                       .TRUE.)
      !p2g
      CALL type_conversion(ipts,p2g,glccSecondShift,veget_mtc,newvegfrac,       &
                       indold_pasture,indagec_pasture,indagec_grass,num_grass_mulagec,     &
                       nagec_herb,nagec_herb,                    &
                       vegagec_pasture,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain, &
                       .TRUE.)
      !p2f
      CALL type_conversion(ipts,p2f,glccSecondShift,veget_mtc,newvegfrac,       &
                       indold_pasture,indagec_pasture,indagec_tree,num_tree_mulagec,     &
                       nagec_herb,nagec_tree,                    &
                       vegagec_pasture,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain, &
                       .TRUE.)
      !c2p
      CALL type_conversion(ipts,c2p,glccSecondShift,veget_mtc,newvegfrac,       &
                       indold_crop,indagec_crop,indagec_pasture,num_pasture_mulagec,     &
                       nagec_herb,nagec_herb,                    &
                       vegagec_crop,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain, &
                       .TRUE.)
      !c2g
      CALL type_conversion(ipts,c2g,glccSecondShift,veget_mtc,newvegfrac,       &
                       indold_crop,indagec_crop,indagec_grass,num_grass_mulagec,     &
                       nagec_herb,nagec_herb,                    &
                       vegagec_crop,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain, &
                       .TRUE.)
      !c2f
      CALL type_conversion(ipts,c2f,glccSecondShift,veget_mtc,newvegfrac,       &
                       indold_crop,indagec_crop,indagec_tree,num_tree_mulagec,     &
                       nagec_herb,nagec_tree,                    &
                       vegagec_crop,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain, &
                       .TRUE.)
    ENDDO
    glccSecondShift_remain(:,:) = glccRemain(:,:)

    !! 3.2 We treat the remaing unrealized transtions from forest. Now we will
    !! start with the 3rd oldest age class and then move to the oldest age class. 

    CALL calc_cover(npts,veget_max,veget_mtc,vegagec_tree,vegagec_grass, &
           vegagec_pasture,vegagec_crop)
    veget_mtc = veget_mtc_begin
 
    IndStart_f = nagec_tree  ! note the indecies and vegetfrac for tree age class
                               ! is from old to young. 
                               ! nagec_tree -3: The 4th youngest age class.

    IndEnd_f = 1               ! oldest-age class forest.

    IF (IndStart_f .LE. 0 .OR. IndEnd_f .LE. 0) THEN
      write(numout,*) 'glcc_MulAgeC: Age class index cannot be negative or zero!'
      STOP
    ENDIF

    ! we start with the 3rd youngest age class and move up to the oldest age
    ! class in the sequence of young->old, as indicated by the .FALSE. parameter
    ! when calling the subroutine type_conversion.
    DO ipts=1,npts
      !f2c
      CALL type_conversion(ipts,f2c,glccSecondShift_remain,veget_mtc,newvegfrac,       &
                       indold_tree,indagec_tree,indagec_crop,num_crop_mulagec,     &
                       IndEnd_f,nagec_herb,                    &
                       vegagec_tree,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp, &
                       glccRemain, &
                       .FALSE., iagec_start=IndStart_f)
      !f2p
      CALL type_conversion(ipts,f2p,glccSecondShift_remain,veget_mtc,newvegfrac,       &
                       indold_tree,indagec_tree,indagec_pasture,num_pasture_mulagec,     &
                       IndEnd_f,nagec_herb,                    &
                       vegagec_tree,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp, &
                       glccRemain, &
                       .FALSE., iagec_start=IndStart_f)
      !f2g
      CALL type_conversion(ipts,f2g,glccSecondShift_remain,veget_mtc,newvegfrac,       &
                       indold_tree,indagec_tree,indagec_grass,num_grass_mulagec,     &
                       IndEnd_f,nagec_herb,                    &
                       vegagec_tree,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain, &
                       .FALSE., iagec_start=IndStart_f)
    ENDDO

    IF (allow_youngest_forest_SecShift) THEN
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

      CALL calc_cover(npts,veget_max,veget_mtc,vegagec_tree,vegagec_grass, &
             vegagec_pasture,vegagec_crop)
      veget_mtc = veget_mtc_begin
 
      ! Note: the setting of index here must be consistent with those of 3.1 and 3.2
      IndStart_f = nagec_tree-1  ! note the indecies and vegetfrac for tree age class
                                 ! is from old to young. 
                                 ! nagec_tree -1: The 2nd youngest age class.

      IndEnd_f = nagec_tree      ! youngest class forest.

      IF (IndStart_f .LE. 0 .OR. IndEnd_f .LE. 0) THEN
        write(numout,*) 'glcc_MulAgeC: Age class index cannot be negative or zero!'
        STOP
      ENDIF

      ! we start with the 3rd youngest age class and move up to the oldest age
      ! class in the sequence of young->old, as indicated by the .FALSE. parameter
      ! when calling the subroutine type_conversion.
      DO ipts=1,npts
        !f2c
        CALL type_conversion(ipts,f2c,glccSecondShift_remain,veget_mtc,newvegfrac,       &
                         indold_tree,indagec_tree,indagec_crop,num_crop_mulagec,     &
                         IndEnd_f,nagec_herb,                    &
                         vegagec_tree,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp, &
                         glccRemain, &
                         .TRUE., iagec_start=IndStart_f)
        !f2p
        CALL type_conversion(ipts,f2p,glccSecondShift_remain,veget_mtc,newvegfrac,       &
                         indold_tree,indagec_tree,indagec_pasture,num_pasture_mulagec,     &
                         IndEnd_f,nagec_herb,                    &
                         vegagec_tree,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp, &
                         glccRemain, &
                         .TRUE., iagec_start=IndStart_f)
        !f2g
        CALL type_conversion(ipts,f2g,glccSecondShift_remain,veget_mtc,newvegfrac,       &
                         indold_tree,indagec_tree,indagec_grass,num_grass_mulagec,     &
                         IndEnd_f,nagec_herb,                    &
                         vegagec_tree,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                         glccRemain, &
                         .TRUE., iagec_start=IndStart_f)
      ENDDO
      !! End of ++Temp++ Section 3.3
    ENDIF

    ! Final handling of some output variables.
    ! we separate the glcc_pftmtc_SecShift
    glcc_pftmtc_SecShift = glcc_pftmtc

    ! we put the remaining glccRemain into the deficit
    glccDefSecShift = -1 * glccRemain
    glccRealSecShift = glccSecondShift - glccRemain

    !*****end block to handle secondary-agriculture shifting cultivation *******


    !! 4. Treat the transtions involving the oldest age classes, which include
    !!    the first-time primary-agriculture cultivation and the net land cover
    !!    transtions

    CALL calc_cover(npts,veget_max,veget_mtc,vegagec_tree,vegagec_grass, &
           vegagec_pasture,vegagec_crop)
    veget_mtc = veget_mtc_begin
 

    ! the variable "glccReal" is originally for storing the realized maxtrix
    ! after considering the constraining and compensation of existing vegetation
    ! fractions. But as this case is not allowed at the moment, we will just
    ! simply put it as the sum of glccPrimaryShift and glccNetLCC
    glccReal(:,:) = glccPrimaryShift+glccNetLCC

    ! We copy the glccReal to glccRemain in order to track the remaining
    ! prescribed transtion fraction after applying each transition by calling
    ! the subroutine "type_conversion". For the moment this is mainly to fufill
    ! the parameter requirement of the type_conversion subroutine.
    glccRemain(:,:) = glccReal(:,:)

    ! We allocate in the sequences of old->young. Within the same age-class
    ! group, we allocate in proportion with existing PFT fractions.
    DO ipts=1,npts
      !f2c
      CALL type_conversion(ipts,f2c,glccReal,veget_mtc,newvegfrac,       &
                       indold_tree,indagec_tree,indagec_crop,num_crop_mulagec,     &
                       nagec_tree,nagec_herb,                    &
                       vegagec_tree,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp, &
                       glccRemain, &
                       .TRUE.)
      !f2p
      CALL type_conversion(ipts,f2p,glccReal,veget_mtc,newvegfrac,       &
                       indold_tree,indagec_tree,indagec_pasture,num_pasture_mulagec,     &
                       nagec_tree,nagec_herb,                    &
                       vegagec_tree,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp, &
                       glccRemain, &
                       .TRUE.)
      !f2g
      CALL type_conversion(ipts,f2g,glccReal,veget_mtc,newvegfrac,       &
                       indold_tree,indagec_tree,indagec_grass,num_grass_mulagec,     &
                       nagec_tree,nagec_herb,                    &
                       vegagec_tree,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain, &
                       .TRUE.)
      !g2c
      CALL type_conversion(ipts,g2c,glccReal,veget_mtc,newvegfrac,       &
                       indold_grass,indagec_grass,indagec_crop,num_crop_mulagec,     &
                       1,nagec_herb,                    &
                       vegagec_grass,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain, &
                       .FALSE., iagec_start=nagec_herb)
      !g2p
      CALL type_conversion(ipts,g2p,glccReal,veget_mtc,newvegfrac,       &
                       indold_grass,indagec_grass,indagec_pasture,num_pasture_mulagec,     &
                       1,nagec_herb,                    &
                       vegagec_grass,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain, &
                       .FALSE., iagec_start=nagec_herb)
      !g2f
      CALL type_conversion(ipts,g2f,glccReal,veget_mtc,newvegfrac,       &
                       indold_grass,indagec_grass,indagec_tree,num_tree_mulagec,     &
                       1,nagec_tree,                    &
                       vegagec_grass,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain, &
                       .FALSE., iagec_start=nagec_herb)
      !p2c
      CALL type_conversion(ipts,p2c,glccReal,veget_mtc,newvegfrac,       &
                       indold_pasture,indagec_pasture,indagec_crop,num_crop_mulagec,     &
                       1,nagec_herb,                    &
                       vegagec_pasture,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain, &
                       .FALSE., iagec_start=nagec_herb)
      !p2g
      CALL type_conversion(ipts,p2g,glccReal,veget_mtc,newvegfrac,       &
                       indold_pasture,indagec_pasture,indagec_grass,num_grass_mulagec,     &
                       1,nagec_herb,                    &
                       vegagec_pasture,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain, &
                       .FALSE., iagec_start=nagec_herb)
      !p2f
      CALL type_conversion(ipts,p2f,glccReal,veget_mtc,newvegfrac,       &
                       indold_pasture,indagec_pasture,indagec_tree,num_tree_mulagec,     &
                       1,nagec_tree,                    &
                       vegagec_pasture,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain, &
                       .FALSE., iagec_start=nagec_herb)
      !c2p
      CALL type_conversion(ipts,c2p,glccReal,veget_mtc,newvegfrac,       &
                       indold_crop,indagec_crop,indagec_pasture,num_pasture_mulagec,     &
                       1,nagec_herb,                    &
                       vegagec_crop,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain, &
                       .FALSE., iagec_start=nagec_herb)
      !c2g
      CALL type_conversion(ipts,c2g,glccReal,veget_mtc,newvegfrac,       &
                       indold_crop,indagec_crop,indagec_grass,num_grass_mulagec,     &
                       1,nagec_herb,                    &
                       vegagec_crop,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain, &
                       .FALSE., iagec_start=nagec_herb)
      !c2f
      CALL type_conversion(ipts,c2f,glccReal,veget_mtc,newvegfrac,       &
                       indold_crop,indagec_crop,indagec_tree,num_tree_mulagec,     &
                       1,nagec_tree,                    &
                       vegagec_crop,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp,&
                       glccRemain, &
                       .FALSE., iagec_start=nagec_herb)
    ENDDO

    glccDefNetPriShift = -1 * glccRemain
    glccRealNetPriShift = glccPrimaryShift + glccNetLCC - glccRemain
    glcc_pftmtc_NetPri = glcc_pftmtc - glcc_pftmtc_SecShift
    glccReal = glccRealSecShift + glccRealNetPriShift
    ! Note here IncreDeficit  includes the deficit from secondary<->agriculgure shifting
    ! cultivation and the primary<->agriculture+NetLCC transitions.
    IncreDeficit = glccDefSecShift + glccDefNetPriShift

    IF (.NOT. glcc_MulAgeC_firstday_done) THEN

      glccReal_tmp = zero

      glccReal_tmp(:,1:12) = glccRealSecShift
      CALL histwrite_p (hist_id_stomate, 'glccRealSecShift', itime, &
           glccReal_tmp, npts*nvm, horipft_index)
      CALL xios_orchidee_send_field ('glccRealSecShift', glccReal_tmp) ! kjpindex,nvm

      glccReal_tmp(:,1:12) = glccRealNetPriShift
      CALL histwrite_p (hist_id_stomate, 'glccRealNetPriShift', itime, &
           glccReal_tmp, npts*nvm, horipft_index)
      CALL xios_orchidee_send_field ('glccRealNetPriShift', glccReal_tmp) ! kjpindex,nvm

      glccReal_tmp(:,1:12) = glccDefSecShift
      CALL histwrite_p (hist_id_stomate, 'glccDefSecShift', itime, &
           glccReal_tmp, npts*nvm, horipft_index)
      CALL xios_orchidee_send_field ('glccDefSecShift', glccReal_tmp) ! kjpindex,nvm

      glccReal_tmp(:,1:12) = glccDefNetPriShift
      CALL histwrite_p (hist_id_stomate, 'glccDefNetPriShift', itime, &
           glccReal_tmp, npts*nvm, horipft_index)
      CALL xios_orchidee_send_field ('glccDefNetPriShift', glccReal_tmp) ! kjpindex,nvm

      DO ivma = 1, nvmap
        WRITE(part_str,'(I2)') ivma
        IF (ivma < 10) part_str(1:1) = '0'
        CALL histwrite_p (hist_id_stomate, 'glcc_pftmtc_SF_'//part_str(1:LEN_TRIM(part_str)), &
             itime, glcc_pftmtc_SecShift(:,:,ivma), npts*nvm, horipft_index)
      ENDDO
      CALL xios_orchidee_send_field ('glcc_pftmtc_SF', glcc_pftmtc_SecShift) ! kjpindex,nvm,nvmap

      DO ivma = 1, nvmap
        WRITE(part_str,'(I2)') ivma
        IF (ivma < 10) part_str(1:1) = '0'
        CALL histwrite_p (hist_id_stomate, 'glcc_pftmtc_NPF_'//part_str(1:LEN_TRIM(part_str)), &
             itime, glcc_pftmtc_NetPri(:,:,ivma), npts*nvm, horipft_index)
      ENDDO
      CALL xios_orchidee_send_field ('glcc_pftmtc_NPF', glcc_pftmtc_NetPri) ! kjpindex,nvm,nvmap
     
      glcc_MulAgeC_firstday_done = .TRUE.
    ENDIF
    WRITE(numout,*) 'End glcc_MulAgeC_firstday'
  END SUBROUTINE glcc_MulAgeC_firstday



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
                     iagec_tree_end,nagec_receive,                    &
                     vegagec_tree,veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp, &
                     glccRemain, &
                     old_to_young, iagec_start)

    IMPLICIT NONE

    !! Input variables
    INTEGER, INTENT(in)                             :: ipts,f2c
    REAL(r_std), DIMENSION(:,:), INTENT(in)         :: glccReal             !! The "real" glcc matrix that we apply in the model
                                                                            !! after considering the consistency between presribed
                                                                            !! glcc matrix and existing vegetation fractions.
    REAL(r_std), DIMENSION(:,:), INTENT(in)          :: newvegfrac
    REAL(r_std), DIMENSION(:,:), INTENT(in)         :: veget_mtc            !! "maximal" coverage fraction of a PFT on the ground
    INTEGER, DIMENSION(:), INTENT(in)               :: indold_tree          !! Indices for PFTs giving out fractions; 
                                                                            !! here use old tree cohort as an example. When iagec_start and 
                                                                            !! a 'old_to_young' or 'young_to_old' sequence is prescribed,
                                                                            !! this index can be possibly skipped.
    INTEGER, DIMENSION(:,:), INTENT(in)             :: indagec_tree         !! Indices for other cohorts except the oldest one giving out fractions; 
                                                                            !! here we use an example of other forest chorts except the oldest one.
    INTEGER, DIMENSION(:,:), INTENT(in)             :: indagec_crop         !! Indices for secondary basic-vegetation cohorts; The youngest age classes
                                                                            !! of these vegetations are going to receive fractions. 
                                                                            !! here we use crop cohorts as an example
    INTEGER, INTENT(in)                             :: num_crop_mulagec     !! number of crop MTCs with more than one age classes
    INTEGER, INTENT(in)                             :: iagec_tree_end       !! End index of age classes in the giving basic types
                                                                            !! (i.e., tree, grass, pasture, crop)
    INTEGER, INTENT(in)                             :: nagec_receive        !! number of age classes in the receiving basic types
                                                                            !! (i.e., tree, grass, pasture, crop), here we can use crop
                                                                            !! as an example, nagec=nagec_herb
    LOGICAL, INTENT(in)                             :: old_to_young         !! an logical variable indicating whether we should handle donation 
                                                                            !! vegetation in a sequence of old->young or young->old. TRUE is for
                                                                            !! old->young. If TRUE, the index will be in increasing sequence of 
                                                                            !! (iagec_start,iagec_tree_end).
    INTEGER, OPTIONAL, INTENT(in)                   :: iagec_start          !! starting index for iagec, this is added in order to handle
                                                                            !! the case of secondary forest clearing.

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

    WRITE(numout,*) 'Entering type_conversion'
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
      IF (old_to_young) THEN
        ! note that both indagec_tree and vegagec_tree are in sequence of old->young
        ! thus iagec_start_proxy must be smaller than iagec_tree_end
        DO iagec=iagec_start_proxy,iagec_tree_end,1
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
              CALL cross_give_receive(ipts,frac_used,veget_mtc,newvegfrac,              &
                       indold_tree,indagec_crop,nagec_receive,num_crop_mulagec, &
                        veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp)
            ELSE
              ! Note also the sequence of indagec_tree is from old->young, so by
              ! increasing iagec, we're handling progressively the old to young
              ! tree age-class PFTs.
              CALL cross_give_receive(ipts,frac_used,veget_mtc,newvegfrac,              &
                       indagec_tree(:,iagec-1),indagec_crop,nagec_receive,num_crop_mulagec, &
                        veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp)
            ENDIF
            frac_begin = frac_begin-frac_used
            vegagec_tree(ipts,iagec)=vegagec_tree(ipts,iagec)-frac_used
            glccRemain(ipts,f2c) = glccRemain(ipts,f2c) - frac_used
          ENDIF
        ENDDO
      ELSE ! in the sequence of young->old
        DO iagec=iagec_start_proxy,iagec_tree_end,-1
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
              CALL cross_give_receive(ipts,frac_used,veget_mtc,newvegfrac,              &
                       indold_tree,indagec_crop,nagec_receive,num_crop_mulagec, &
                        veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp)
            ELSE
              ! Note also the sequence of indagec_tree is from old->young, so by
              ! increasing iagec, we're handling progressively the old to young
              ! tree age-class PFTs.
              CALL cross_give_receive(ipts,frac_used,veget_mtc,newvegfrac,              &
                       indagec_tree(:,iagec-1),indagec_crop,nagec_receive,num_crop_mulagec, &
                        veget_max,glcc_pft,glcc_pftmtc,glcc_pft_tmp)
            ENDIF
            frac_begin = frac_begin-frac_used
            vegagec_tree(ipts,iagec)=vegagec_tree(ipts,iagec)-frac_used
            glccRemain(ipts,f2c) = glccRemain(ipts,f2c) - frac_used
          ENDIF
        ENDDO
      ENDIF
    !ENDDO

  END SUBROUTINE type_conversion

END MODULE stomate_glcchange_MulAgeC
