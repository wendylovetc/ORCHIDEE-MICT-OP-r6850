! =================================================================================================================================
! MODULE          : stomate_npp
!
! CONTACT         : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE         : IPSL (2006)
!                 This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF          This modules calculates NPP: Maintenance and growth respiration
!!
!!\n DESCRIPTION: We calculate first the maintenance respiration. This is substracted from the
!!                allocatable biomass (and from the present biomass if the GPP is too low).\n
!!                Of the rest, a part is lost as growth respiration, while the other part is
!!                effectively allocated.
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S)	:
!!
!! SVN          :
!! $HeadURL$
!! $Date$
!! $Revision$
!! \n
!_ ================================================================================================================================

MODULE stomate_npp

  ! modules used:
  USE xios_orchidee
  USE ioipsl_para
  USE stomate_data
  USE constantes
  USE constantes_soil
  USE pft_parameters
  USE crop_alloc

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC npp_calc,npp_calc_clear

  LOGICAL, SAVE                                              :: firstcall_npp = .TRUE.         !! first call
!$OMP THREADPRIVATE(firstcall_npp)

CONTAINS

!! ================================================================================================================================
!! SUBROUTINE 	: npp_calc_clear
!!
!>\BRIEF        : Set the flag ::firstcall_npp to .TRUE. and as such activate section 
!! 1.1 of the subroutine npp_calc (see below).\n
!_ ================================================================================================================================

  SUBROUTINE npp_calc_clear
    firstcall_npp=.TRUE.
  END SUBROUTINE npp_calc_clear





!! ================================================================================================================================
!! SUBROUTINE	: npp_calc
!!
!>\BRIEF        Calculate NPP as the difference between GPP and respiration (= growth + maintenance respiration).
!!              Update biomass of all compartments after calculating respiration and allocation.
!!
!!
!! DESCRIPTION  : NPP is calculated from three components: Gross Primary Productivity (GPP), maintenance respiration 
!! and growth respiration (all in @tex $ gC.m^{-2}dt^{-1} $ @endtex), following the convention that positive fluxes denote 
!! fluxes plants to the atmosphere. GPP is the input variable from which, in the end, NPP or total allocatable biomass 
!! @tex $(gC.m^{-2}dt^{-1}))$ @endtex is calculated. Net primary production is then calculated as:\n	
!! NPP = GPP - growth_resp - maint-resp   [eq. 1]\n   
!!	
!! The calculation of maintenance respiration is done in routine stomate_resp.f90. Maintenance respiration is calculated for 
!! the whole plant and is therefore removed from the total allocatable biomass. In order to prevent all allocatable biomass 
!! from being used for maintenance respiration, a limit fraction of total allocatable biomass, tax_max, is defined (in 
!! variables declaration). If maintenance respiration exceeds tax_max (::bm_tax_max), the maximum allowed allocatable biomass
!! will be respired and the remaining respiration, required in excess of tax_max, is taken out from tissues already present in
!! the plant (biomass).\n  
!! 
!! After total allocatable biomass has been updated by removing maintenance respiration, total allocatable biomass is distributed 
!! to all plant compartments according to the f_alloc fractions calculated in stomate_alloc.f90.\n
!!
!! Growth respiration is calculated as a fraction of allocatable biomass for each part of the plant. The fraction coefficient 
!! ::frac_growth_resp is defined in stomate_constants.f90 and is currently set to be the same for all plant compartments. 
!! Allocatable biomass of all plant compartments are updated by removing what is lost through growth respiration. Net allocatable
!! biomass (total allocatable biomass after maintenance and growth respiration) is added to the current biomass for  each plant 
!! compartment.
!!
!! Finally, leaf age and plant age are updated. Leaf age is described with the concept of "leaf age classes". A number of leaf 
!! age classes (nleafages) is defined in stomate_constants.f90. Each leaf age class contains a fraction (::leaf_frac) of the 
!! total leaf biomass. When new biomass is added to leaves, the age of the biomass in the youngest leaf age class is decreased. 
!! The fractions of leaves in the other leaf ages classes are also updated as the total biomass has increased. Plant age is 
!! updated first by increasing the age of the previous biomass by one time step, and then by adjusting this age as the average 
!! of the ages of the previous and the new biomass.
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::npp
!!
!! REFERENCE(S)	: 
!! - F.W.T.Penning De Vries, A.H.M. Brunsting, H.H. Van Laar. 1974. Products, requirements and efficiency of biosynthesis a 
!! quantitative approach. Journal of Theoretical Biology, Volume 45, Issue 2, June 1974, Pages 339-377.
!! 
!! FLOWCHART : 
!! \latexonly
!! \includegraphics[scale=0.14]{stomate_npp_flow.jpg}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE npp_calc (npts, dt, &
       PFTpresent, &
       t2m, tsoil, lai, rprof,  &
!!!@xchen lai for leaf age class1-4 for caculating PC in Wu's science paper
       lai_cl1, lai_cl2, lai_cl3,lai_cl4,delta_turnover, var_date, &
       gpp, f_alloc, bm_alloc, resp_maint_part,&
       biomass, leaf_age, leaf_frac, age, &
       alloc_phy,FtoPHY, & !! yidi
       bmFFB_alloc, bmphy_alloc, bm_phytomer, bm_FFB, PHYbm, FFBbm, & !! yidi      
       resp_maint, resp_growth, npp,  &
!!!! crop variables
       ! for crop bm_alloc
       in_cycle, deltai, dltaisen, ssla, pgrain, deltgrain, reprac, &
       nger, nlev, ndrp, nlax, nmat, nrec, &
       c_reserve, c_leafb, slai, tday_counter, veget_max, &
!!!! end crop variables, xuhui
!gmjc
       sla_calc, sla_age1,N_limfert)
!end gmjc   
   
!! 0 Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                                :: npts             !! Domain size - number of pixels (unitless)
    REAL(r_std), INTENT(in)                                   :: dt               !! Time step (days)
    LOGICAL, DIMENSION(npts,nvm), INTENT(in)                  :: PFTpresent       !! PFT exists (true/false)
    REAL(r_std), DIMENSION(npts), INTENT(in)                  :: t2m              !! Temperature at 2 meter (K)
    REAL(r_std), DIMENSION(npts,nslm), INTENT(in)             :: tsoil            !! Soil temperature of each soil layer (K)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: lai              !! PFT leaf area index (unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: rprof            !! PFT root depth as calculated in stomate.f90
                                                                                  !! from root profile parameter humcste (m) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: gpp              !! PFT gross primary productivity 
                                                                                  !! @tex $(gC.m^{-2}dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nparts), INTENT(in)       :: f_alloc          !! Fraction of total allocatable biomass that 
                                                                                  !! goes into each plant part (unitless)
!! yidi
    REAL(r_std), DIMENSION(npts,nvm,nphs), INTENT(in)         :: alloc_phy        !! Fraction of sapabove allocatable biomass that 
                                                                                  !! goes into each phytomer (unitless)
    REAL(r_std), DIMENSION(npts,nvm,nphs), INTENT(in)         :: FtoPHY           !! Fraction of allocatable biomass for each photomer 
                                                                                  !! that goes into each FFB (unitless)
!! yidi end

    REAL(r_std), DIMENSION(npts,nvm,nparts), INTENT(in)       :: resp_maint_part  !! Maintenance respiration of different plant 
                                                                                  !! parts @tex $(gC.m^{-2}dt^{-1})$ @endtex
!!!!! crop var
    LOGICAL, DIMENSION(npts,nvm), INTENT(in)                 :: in_cycle
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: deltai
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: dltaisen
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)              :: ssla
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: pgrain
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: deltgrain
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: reprac
    INTEGER(i_std), DIMENSION(npts,nvm), INTENT(in)              :: nger
    INTEGER(i_std), DIMENSION(npts,nvm), INTENT(in)              :: nlev
    INTEGER(i_std), DIMENSION(npts,nvm), INTENT(in)              :: ndrp
    INTEGER(i_std), DIMENSION(npts,nvm), INTENT(in)              :: nlax
    INTEGER(i_std), DIMENSION(npts,nvm), INTENT(in)              :: nmat
    INTEGER(i_std), DIMENSION(npts,nvm), INTENT(in)              :: nrec
    INTEGER(i_std), INTENT(in)                                   :: tday_counter
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)                 :: veget_max
!!!!! xuhui

    !! 0.2 Output variables

    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)             :: resp_maint       !! PFT maintenance respiration 
                                                                                  !! @tex $(gC.m^{-2}dt^{-1})$ @endtex		    
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)             :: resp_growth      !! PFT growth respiration 
                                                                                  !! @tex $(gC.m^{-2}dt^{-1})$ @endtex				
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)             :: npp              !! PFT net primary productivity 
                                                                                  !! @tex $(gC.m^{-2}dt^{-1})$ @endtex		
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(out) :: bm_alloc    !! PFT biomass increase, i.e. NPP per plant part 
                                                                                  !! @tex $(gC.m^{-2}dt^{-1})$ @endtex
!! yidi
    REAL(r_std), DIMENSION(npts,nvm,nphs), INTENT(out)         :: bmFFB_alloc     !! OP-FFB biomass increase, i.e. NPP per FFB 
                                                                                  !! @tex $(gC.m^{-2}dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nphs), INTENT(out)         :: bmphy_alloc     !! OP-phy biomass increase, i.e. NPP per phytomer (FFB excluded) 
                                                                                  !! @tex $(gC.m^{-2}dt^{-1})$ @endtex
   !REAL(r_std), DIMENSION(npts,nvm,nphs), INTENT(in)      :: phytomer_age     !! age of the phytomers (days)
   !REAL(r_std), DIMENSION(npts,nvm),INTENT(out)              :: bm_alloc_sapabold !! sapabove mass after npp allocation, before FFB substraction  
                                                                                  !! @tex $(gC.m^{-2})$ @endtex			
    REAL(r_std), DIMENSION(npts,nvm),INTENT(inout)            :: FFBbm            !! FFB mass, from sapabove   
                                                                                  !! @tex $(gC.m^{-2})$ @endtex			
    REAL(r_std), DIMENSION(npts,nvm,nphs),INTENT(inout)       :: bm_FFB           !! FRUIT mass for each PHYTOMER, from sapabove   
                                                                                  !! @tex $(gC.m^{-2})$ @endtex			
    REAL(r_std), DIMENSION(npts,nvm),INTENT(inout)            :: PHYbm            !! PHYTOMER mass, from sapabove   
                                                                                  !! @tex $(gC.m^{-2})$ @endtex			
    REAL(r_std), DIMENSION(npts,nvm,nphs),INTENT(inout)       :: bm_phytomer      !! Each PHYTOMER mass, from sapabove   
                                                                                  !! @tex $(gC.m^{-2})$ @endtex			
!! yidi

    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)             :: lai_cl1          !! @xchen to add lai of each leaf class age and for                                                                                  !!calculating PC in wu's science paper
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)             :: lai_cl2
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)             :: lai_cl3        
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)             :: lai_cl4        
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: delta_turnover            !! Intermediate variable for turnover ??
    INTEGER(i_std),INTENT(in)                            :: var_date               !! xchen Date (days) 


    !! 0.3 Modified variables

    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: biomass   !! PFT total biomass of each plant part 
                                                                                  !! @tex $(gC.m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nleafages), INTENT(inout) :: leaf_age         !! PFT age of different leaf age classes (days)
    REAL(r_std), DIMENSION(npts,nvm,nleafages), INTENT(inout) :: leaf_frac        !! PFT fraction of total leaves in leaf age 
                                                                                  !! class (unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: age              !! PFT age (years)
!gmjc
    ! N fertilization limitation factor for grassland Vcmax and SLA
    REAL(r_std), DIMENSION(:,:), INTENT(in)                    :: N_limfert
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: sla_calc
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: sla_age1
!end gmjc
!!!!! crop var
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)             :: c_reserve
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)             :: c_leafb
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)             :: slai
!!!!! xuhui

    !! 0.4 Local variables
!gmjc
    ! weighted leaf age (leaf_age * fraction of the age class)
    REAL(r_std), DIMENSION(npts,nvm)                      :: leaf_age_w
    REAL(r_std), DIMENSION(npts,nvm)                      :: sla_age2
    REAL(r_std), DIMENSION(npts,nvm)                      :: sla_age3
    REAL(r_std), DIMENSION(npts,nvm)                      :: sla_age4
    ! SLA max and SLA min affected by N fertilization
    REAL(r_std), DIMENSION(npts,nvm)                      :: sla_max_Nfert
    REAL(r_std), DIMENSION(npts,nvm)                      :: sla_min_Nfert
!end gmjc
!BEGINNVADD
    ! npp above (gC/m**2) for the whole plant
    REAL(r_std), DIMENSION(npts,nvm)                           :: npp_above
    ! npp below (gC/m**2) for the whole plant
    REAL(r_std), DIMENSION(npts,nvm)                           :: npp_below
   ! NPP per plant part
    REAL(r_std), DIMENSION(npts,nvm,nparts)                    :: npp_part
!ENDNVADD
    REAL(r_std), DIMENSION(npts,nvm)                          :: t_root           !! Root temperature (convolution of root and 
                                                                                  !! soil temperature profiles)(K)
    REAL(r_std), DIMENSION(npts,nvm,nparts)                   :: coeff_maint      !! PFT maintenance respiration coefficients of 
                                                                                  !! different plant compartments at 0 deg C 
                                                                                  !! @tex $(g.g^{-1}dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nparts)                       :: t_maint          !! Temperature which is pertinent for maintenance
                                                                                  !! respiration, which is air/root temperature for
                                                                                  !! above/below-ground compartments (K)
    REAL(r_std), DIMENSION(npts)                              :: rpc              !! Scaling factor for integrating vertical soil 
                                                                                  !! profiles (unitless)
    REAL(r_std), DIMENSION(npts)                              :: tl               !! Long term annual mean temperature (C)
    REAL(r_std), DIMENSION(npts)                              :: slope            !! Slope of maintenance respiration coefficient
                                                                                  !! (1/K)
    REAL(r_std), DIMENSION(npts,nparts)                       :: resp_growth_part !! Growth respiration of different plant parts
                                                                                  !! @tex $(gC.m^{-2}dt^{-1})$ @endtex		
    REAL(r_std), DIMENSION(npts,nvm)                          :: bm_alloc_tot     !! Allocatable biomass for the whole plant
                                                                                  !! @tex $(gC.m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts)                              :: bm_add           !! Biomass increase @tex $(gC.m^{-2})$ @endtex		
    REAL(r_std), DIMENSION(npts)                              :: bm_new           !! New biomass @tex $(gC.m^{-2})$ @endtex	
    REAL(r_std), DIMENSION(npts,nvm)                          :: leaf_mass_young  !! Leaf mass in youngest age class 
                                                                                  !! @tex $(gC.m^{-2})$ @endtex		
    REAL(r_std), DIMENSION(npts,nvm)                          :: lm_old           !! Leaf mass after maintenance respiration 
                                                                                  !! @tex $(gC.m^{-2})$ @endtex			
    REAL(r_std), DIMENSION(npts,nvm)                          :: bm_create        !! Biomass created when biomass<0 because of dark
                                                                                  !! respiration @tex $(gC.m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts)                              :: bm_tax_max       !! Maximum part of allocatable biomass used for 
                                                                                  !! respiration @tex $(gC.m^{-2})$ @endtex	
    REAL(r_std), DIMENSION(npts)                              :: bm_pump          !! Biomass that remains to be taken away 
                                                                                  !! @tex $(gC.m^{-2})$ @endtex
    INTEGER(i_std)                                            :: i,j,k,l,m,p,jv   !! Indeces(unitless)
    INTEGER(i_std)                                            :: ier              !! Error handling

!! yidi
    CHARACTER(LEN=10)                                         :: part_str         !! string suffix indicating an index
!! yidi
!_ ================================================================================================================================

    IF (printlev>=3) WRITE(numout,*) 'Entering npp##############################'
    
 !! 1. Initializations
    
    !! 1.1 First call
    IF ( firstcall_npp ) THEN

       !! 1.1.2 Output message
       !  Write message including value used for tax_max 	
       IF (printlev >= 2) THEN
          WRITE(numout,*) 'npp:'

          WRITE(numout,*) '   > max. fraction of allocatable biomass used for'// &
               ' maint. resp.:', tax_max
       END IF
       firstcall_npp = .FALSE.

    ENDIF ! End if first call

    !! 1.2 Set output variables to zero
    bm_alloc(:,:,:,:) = zero
    resp_maint(:,:) = zero
    resp_growth(:,:) = zero
    npp(:,:) = zero
!BEGINNVADD
    npp_above(:,:) = zero
    npp_below(:,:) = zero
    npp_part(:,:,:) = zero
    !!@xchen
    lai_cl1(:,:) = zero
    lai_cl2(:,:) = zero                                                       
    lai_cl3(:,:) = zero                                                       
    lai_cl4(:,:) = zero                                                       
    !! yidi
    !FFBbm(:,:) = zero
    !PHYbm(:,:) = zero
   ! FFBharvest(:,:) =zero
    bmFFB_alloc(:,:,:) = zero
    bmphy_alloc(:,:,:) = zero
    !bm_FFB(:,:,:) = zero
    !bm_phytomer(:,:,:) = zero
   ! bm_alloc_sapabold(:,:)=zero
!ENDNVADD

!!  ydtemp
    !WRITE(numout,*) 'yd: npp1 bmFFB_alloc(1,62,1) initial',bmFFB_alloc(1,62,1) !! yd output

    !! 1.3 Total allocatable biomass
    ! total allocatable biomass during this time step determined from GPP.
    ! GPP was calculated as CO2 assimilation in enerbil.f90
    !WRITE(numout,*) 'zd bm_alloc 1','bm_alloc_tot(1,10)',bm_alloc_tot(1,10),'gpp(1,10)',gpp(1,10),'dt',dt
    !!xchen201702
    IF ((var_date .GT. 195) .AND. (var_date .LT. 330)) THEN 
    bm_alloc_tot(:,:) = gpp(:,:) * dt
    ELSE
    bm_alloc_tot(:,:) = gpp(:,:) * dt   
    ENDIF
    
    !WRITE(numout,*) 'yd: npp2  bm_alloc_tot(1,2:7)',bm_alloc_tot(1,2:7)  !!!!!yd output
    WRITE(numout,*) 'yd: npp2  bm_alloc_tot(1,62:67)',bm_alloc_tot(1,62:67)  !!!!!yd output

    !DO j = 2, nvm
       !IF (ok_oilpalm .AND. is_oilpalm(j)) THEN
         !WRITE(numout,*) 'yd: npp2  nvm=',j,'bm_alloc_tot(1,j)',bm_alloc_tot(1,j)  !!!!!yd output
         !WRITE(numout,*) 'yd: npp2 NVM=',j,'before npp alloc: f_alloc(i,j,isapabove) ',f_alloc(i,j,isapabove) 
         !WRITE(numout,*) 'yd: npp2 NVM=',j,'before npp alloc: alloc_phy(i,j,:) ',alloc_phy(i,j,:) 
         !WRITE(numout,*) 'yd: npp2 NVM=',j,'before npp alloc: SUM(alloc_phy(i,j,:)) ',SUM(alloc_phy(i,j,:)) 
         !WRITE(numout,*) 'yd: npp2 NVM=',j,'before npp alloc: FtoPHY(i,j,:) ',FtoPHY(i,j,:) 
      ! ENDIF
    !ENDDO
!    WRITE(numout,*) 'biomass(1,12:14,:,icarbon) before PFT loop:',biomass(1,12:14,:,icarbon)
    !! 3. Calculate maintenance and growth respiration
    ! First, total maintenance respiration for the whole plant is calculated by summing maintenance 
    ! respiration of the different plant compartments. Then, maintenance respiration is subtracted 
    ! from whole-plant allocatable biomass (up to a maximum fraction of the total allocatable biomass). 
    ! Growth respiration is then calculated for each plant compartment as a fraction of remaining 
    ! allocatable biomass for this compartment. NPP is calculated by substracting total autotrophic 
    ! respiration from GPP i.e. NPP = GPP - maintenance resp - growth resp.
    !WRITE(numout,*) 'zd nppcalc1 biomass(1,10,ileaf,icarbon)',biomass(1,10,ileaf,icarbon),'bm_pump(1)',bm_pump(1),'resp_maint_part(1,10,ileaf)',resp_maint_part(1,10,ileaf),'resp_maint(1,10)',resp_maint(1,10)
    DO j = 2,nvm	! Loop over # of PFTs

       !! 3.1 Maintenance respiration of the different plant parts
       !      Maintenance respiration of the different plant parts is calculated in 
       !      stomate_resp.f90 as a function of the plant's temperature, 
       !      the long term temperature and plant coefficients 
       !      VPP killer:
       resp_maint(:,j) = zero

       !  Following the calculation of hourly maintenance respiration, verify that 
       !  the PFT has not been killed after calcul of resp_maint_part in stomate.
       DO k= 1, nparts
          WHERE (PFTpresent(:,j))
             resp_maint(:,j) = resp_maint(:,j) + resp_maint_part(:,j,k)
          ENDWHERE
       ENDDO
       
       !! 3.2 Substract maintenance respiration from allocatable biomass
       !      The total maintenance respiration calculated in 3.2 is substracted  from the newly 
       !      produced allocatable biomass (bm_alloc_tot). However, ensure that not all allocatable 
       !      biomass is removed by setting a maximum to the fraction of allocatable biomass used 
       !      for maintenance respiration: tax_max. If the maintenance respiration is larger than 
       !      tax_max,the amount tax_max is taken from allocatable biomass, and the remaining of 
       !      maintenance respiration is taken from the tissues themselves (biomass). We suppose 
       !      that respiration is not dependent on leaf age -> therefore the leaf age structure is 
       !      not changed. 
       !      The maximum fraction of allocatable biomass used for respiration is defined as tax_max. 
       !      The value of tax_max is set in the declarations section (0.4 Local variables) of this
       !      routine
       bm_tax_max(:) = tax_max * bm_alloc_tot(:,j)

!       IF (j.EQ.10) WRITE(numout,*) 'zd bm_alloc 3','bm_alloc_tot(1,10)',bm_alloc_tot(1,10)
       DO i = 1, npts	! Loop over # of pixels

          ! If there is enough allocatable biomass to cover maintenance respiration, 
	  ! then biomass associated with maintenance respiration is removed from allocatable biomass
          IF ( bm_alloc_tot(i,j) .GT. zero ) THEN
               IF (is_oilpalm(j)) THEN
                  WRITE(numout,*) 'yd: npp2.1  nvm=',j,'bm_alloc_tot(1,j)',bm_alloc_tot(1,j)  !!!!!yd output
               ENDIF
               IF ( ( resp_maint(i,j) * dt ) .LT. bm_tax_max(i) )  THEN	
	
                  bm_alloc_tot(i,j) = bm_alloc_tot(i,j) - resp_maint(i,j) * dt

                  WRITE(numout,*) 'yd: npp2.2  nvm=',j,'bm_alloc_tot(1,j)',bm_alloc_tot(1,j)  !!!!!yd output
                  ! If there is not enough allocatable biomass to cover maintenance respiration, the  
                  ! - maximum allowed allocatable biomass (bm_tax_max) is removed from allocatable biomass.
               ELSE
	     
                  bm_alloc_tot(i,j) = bm_alloc_tot(i,j) - bm_tax_max(i)

                  IF (is_oilpalm(j)) THEN
                     WRITE(numout,*) 'yd: npp2.3  nvm=',j,'bm_alloc_tot(1,j)',bm_alloc_tot(1,j)  !!!!!yd output
                  ENDIF
                  ! ::bm_pump is the amount of maintenance respiration that exceeds the maximum allocatable biomass
                  ! This amount of biomass still needs to be respired and will be removed from tissues biomass of each
                  ! plant compartment
                  bm_pump(i) = resp_maint(i,j) * dt - bm_tax_max(i)

                  ! The biomass is removed from each plant compartment tissues as the ratio of the maintenance		
                  ! respiration of the plant compartment to the total maintenance respiration (resp_maint_part/resp_maint)
                  biomass(i,j,ileaf,icarbon) = biomass(i,j,ileaf,icarbon) - &
                       bm_pump(i) * resp_maint_part(i,j,ileaf) / resp_maint(i,j)
                  biomass(i,j,isapabove,icarbon) = biomass(i,j,isapabove,icarbon) - &
                       bm_pump(i) * resp_maint_part(i,j,isapabove) / resp_maint(i,j)
                  biomass(i,j,isapbelow,icarbon) = biomass(i,j,isapbelow,icarbon) - &
                       bm_pump(i) * resp_maint_part(i,j,isapbelow) / resp_maint(i,j)
                  biomass(i,j,iroot,icarbon) = biomass(i,j,iroot,icarbon) - &
                       bm_pump(i) * resp_maint_part(i,j,iroot) / resp_maint(i,j)
                  biomass(i,j,ifruit,icarbon) = biomass(i,j,ifruit,icarbon) - &
                       bm_pump(i) * resp_maint_part(i,j,ifruit) / resp_maint(i,j)
                  biomass(i,j,icarbres,icarbon) = biomass(i,j,icarbres,icarbon) - &
                       bm_pump(i) * resp_maint_part(i,j,icarbres) / resp_maint(i,j)
!! yidi
 !                 WRITE(numout,*) 'yd: npp3, nvm=2:7,(biomass(isapabove)',biomass(i,2:7,isapabove,icarbon)
 !                 WRITE(numout,*) 'yd: npp3, nvm=2,  (biomass(all)',biomass(i,2,:,icarbon)
                  IF (ok_oilpalm) THEN
                     IF (is_oilpalm(j)) THEN
                         WRITE(numout,*) 'yd: npp3, nvm=',j,'(biomass(isapabove)',biomass(i,j,isapabove,icarbon)
                         WRITE(numout,*) 'yd: npp3, nvm=',j,'(biomass(all)',biomass(i,j,:,icarbon)
                         IF (biomass(i,j,isapabove,icarbon) .GT. min_stomate) THEN
                            DO p = 1, nphs
                               bm_phytomer(i,j,p) = bm_phytomer(i,j,p) - &
                                                   bm_pump(i) * resp_maint_part(i,j,isapabove) / resp_maint(i,j) * &
                                                   bm_phytomer(i,j,p) / biomass(i,j,isapabove,icarbon)
                               !WRITE(numout,*) 'yd: npp3, nvm=',j,'nphs=',p,'bm_phytomer(i,j,p)',bm_phytomer(i,j,p)
                            ENDDO
                            PHYbm(i,j) = PHYbm(i,j) - &
                                         bm_pump(i) * resp_maint_part(i,j,isapabove) / resp_maint(i,j) * &
                                         PHYbm(i,j) / biomass(i,j,isapabove,icarbon)
                            WRITE(numout,*) 'yd: npp3, nvm=',j,'PHYbm(i,j)',PHYbm(i,j)
                            WRITE(numout,*) 'yd: npp3, nvm=',j,'bm_phytomer(i,j,:)',bm_phytomer(i,j,:)
                        ENDIF
                     ENDIF
!! may not separate is_oilpalm and is_oilpalm_ffbharvestsince bm_FFB=0 
                     IF (is_oilpalm_ffbharvest(j)) THEN
                         !WRITE(numout,*) 'yd: npp4, nvm=',j,'(biomass(isapabove)',biomass(i,j,isapabove,icarbon)
                         IF (biomass(i,j,isapabove,icarbon) .GT. min_stomate) THEN
                            DO p = 1, nphs
                               bm_FFB(i,j,p) = bm_FFB(i,j,p) - &
                                                   bm_pump(i) * resp_maint_part(i,j,isapabove) / resp_maint(i,j) * &
                                                   bm_FFB(i,j,p) / biomass(i,j,isapabove,icarbon)
                               !WRITE(numout,*) 'yd: npp4, nvm=',j,'nphs=',p,'bm_FFB(i,j,p)',bm_FFB(i,j,p)
                            ENDDO
                            FFBbm(i,j) = FFBbm(i,j) - &
                                         bm_pump(i) * resp_maint_part(i,j,isapabove) / resp_maint(i,j) * &
                                         FFBbm(i,j) / biomass(i,j,isapabove,icarbon)
                            WRITE(numout,*) 'yd: npp4, nvm=',j,'FFBbm(i,j)',FFBbm(i,j)
                            WRITE(numout,*) 'yd: npp4, nvm=',j,'bm_FFB(i,j,:)',bm_FFB(i,j,:)
                        ENDIF
                     ENDIF
                  ENDIF

!! yidi
!BEGINNVADD
          npp_part (i,j,:)= npp_part (i,j,:) - &
                                bm_pump(i) * resp_maint_part(i,j,:) / resp_maint(i,j)
!ENDNVADD
             ENDIF
          ELSE
             biomass(i,j,icarbres,icarbon) = biomass(i,j,icarbres,icarbon) + & 
                  bm_alloc_tot(i,j) - resp_maint(i,j) * dt 
             bm_alloc_tot(i,j) = 0. 
          ENDIF	! End if there is enough allocatable biomass to cover maintenance respiration

       ENDDO   ! Fortran95: WHERE - ELSEWHERE construct
!       IF (j.EQ.10) WRITE(numout,*) 'zd bm_alloc 4','bm_alloc_tot(1,10)',bm_alloc_tot(1,10)
       !IF ((j.GE.2) .AND.(j.LE.7)) THEN
       !   WRITE(numout,*) 'yd: npp3,after maintenance respiration nvm=2:7,(biomass(isapabove)',biomass(i,2:7,isapabove,icarbon)
       !   WRITE(numout,*) 'yd: npp3,after maintenance respiration nvm=2,  (biomass(all)',biomass(i,2,:,icarbon)
       !ENDIF
       !! 3.3 Allocate allocatable biomass to different plant compartments.
       !      The amount of allocatable biomass of each compartment is a fraction according f_alloc of total 
       !      allocatable biomass (the f_alloc of the different plant parts are calculated in stomate_alloc.f90)
!       IF (j.EQ.10) WRITE(numout,*) 'zd bm_alloc 5','bm_alloc(1,10,:,icarbon)',bm_alloc(1,10,:,icarbon),'f_alloc(1,10,:)',f_alloc(1,10,:)
       IF (ok_LAIdev(j)) THEN 
!          WRITE(numout,*) 'slai before npp_alloc: ',j, 'pft, ', slai(1,j)
!          WRITE(numout,*) 'ssla before npp_alloc: ',j, 'pft, ', ssla(1,j)
!          WRITE(numout,*) 'biomass(1,j,:,icarbon): ',j, 'pft, ',biomass(1,j,:,icarbon)
!          WRITE(numout,*) 'bm_alloc(1,j,:,icarbon) before npp_alloc: ',j, 'pft, ', bm_alloc(1,j,:,icarbon)
          DO i = 1, npts
             ! we call the crop bm allocation subroutine        
             if ( in_cycle(i,j) ) then
                 write(numout,*) '(i,j)', i, j
                 write(numout,*) 'in_cycle(i,j)',in_cycle(i,j)
                 if (veget_max(i,j) .gt. 0.) then 
                     ! there will be error in crop_alloc if vegetmax for a certain crop is 0
                     call crop_bmalloc(in_cycle(i, j),         &
                               deltai(i, j),           &
                               dltaisen(i, j),         &
                               ssla(i, j),             &
                               pgrain(i, j),           &
                               deltgrain(i, j),          &
                               reprac(i, j),           &
                               nger(i, j),             &
                               nlev(i,j),             &
                               ndrp(i,j),             &
                               nlax(i,j),             &    ! input 
                               nmat(i,j),             &
                               nrec(i, j),            &
!                               f_crop_recycle(i,j),   & 
                               bm_alloc_tot(i,j),     &    ! input
                               biomass(i, j, :, icarbon),&    ! input
                               c_reserve(i,j),        &    ! inout
                               c_leafb(i,j),          &    ! inout
                               bm_alloc(i, j, :, icarbon),     &    ! inout
                               SP_densitesem(j),      &
                               SP_pgrainmaxi(j),      &
                               SP_tigefeuil(j),       &
                               SP_slamax(j),          &
                               slai(i,j),             &
                               tday_counter)               ! parameter
                 endif ! j==11
              endif ! in_cycle
   
          ENDDO ! npts
       ELSE  ! natural vegetation (is not crop)
        !   IF ((j .GE. 2) .AND. (j .LE. 7)) THEN
        !      WRITE(numout,*) 'yd: npp5  nvm=',j,', gpp(:,j)*dt', gpp(:,j)*dt
        !      WRITE(numout,*) 'yd: npp5  resp_maint(:,j)', resp_maint(:,j)
        !      WRITE(numout,*) 'yd: npp5  bm_alloc_tot(:,j)', bm_alloc_tot(:,j)
        !      WRITE(numout,*) 'yd: npp5  bm_alloc(:,j,:,icarbon) before npp_alloc', bm_alloc(:,j,:,icarbon)
        !   ENDIF
           IF ( (printlev>=5) .AND. is_oilpalm(j) ) THEN
               WRITE(numout,*) 'yd: npp5  nvm=',j,'gpp(:,j)*dt', gpp(:,j)*dt
               WRITE(numout,*) 'yd: npp5  bm_alloc_tot(:,j)', bm_alloc_tot(:,j)
               WRITE(numout,*) 'yd: npp5  resp_maint(:,j)', resp_maint(:,j)
               WRITE(numout,*) 'ye: npp5  bm_alloc(:,j,:,icarbon) before npp_alloc', bm_alloc(:,j,:,icarbon)
           ENDIF
           DO k = 1, nparts
              bm_alloc(:,j,k,icarbon) = f_alloc(:,j,k) * bm_alloc_tot(:,j)
           ENDDO
           !WRITE(numout,*) 'yd: npp5.25 before npp alloc: nvm=',j,' f_alloc(i,j,isapabove) ',f_alloc(i,j,isapabove) 
!! yidi    

         !  IF ((j .GE. 2) .AND. (j .LE. 7)) THEN
         !     WRITE(numout,*) 'yd: npp 6 bm_alloc(:,j,:,icarbon),NVM=',j,' after npp_alloc', bm_alloc(:,j,:,icarbon)
         !  ENDIF
           IF (is_oilpalm(j)) THEN
              WRITE(numout,*) 'yd: npp 6 bm_alloc(:,j,:,icarbon),NVM=',j,' after npp_alloc', bm_alloc(:,j,:,icarbon)
           ENDIF

           IF (ok_oilpalm) THEN
              IF (is_oilpalm(j)) THEN
              !WRITE(numout,*) 'bm_alloc(1,j,:,icarbon) before npp_alloc: ',j,'nvm, ',bm_alloc(1,j,:,icarbon)

                 DO p = 1, nphs
                    bmFFB_alloc(:,j,p) = bm_alloc_tot(:,j)* f_alloc(:,j,isapabove) &
                            * alloc_phy(:,j,p) * FtoPHY(:,j,p)
                    bmphy_alloc(:,j,p) = bm_alloc_tot(:,j)* f_alloc(:,j,isapabove) &
                            * alloc_phy(:,j,p) * (1-FtoPHY(:,j,p))
                    !WRITE(numout,*) 'yd: npp5.5:NVM=',j,' before npp alloc: alloc_phy(i,j,p) ',alloc_phy(i,j,p),'nphs',p 
                    !WRITE(numout,*) 'yd: npp5.5:NVM=',j,' test1: ',alloc_phy(i,j,p)*bm_alloc(i,j,isapabove,icarbon),'nphs',p 
                    !WRITE(numout,*) 'yd: npp5.5:NVM=',j,' test2: ',alloc_phy(i,j,p)*FtoPHY(i,j,p),'nphs',p 
                 !   WRITE(numout,*) 'yd: npp6 after alloc: bmFFB_alloc(1,j,p) ',j,'nvm',p,'nphs', bmFFB_alloc(1,j,p) 
                 !   WRITE(numout,*) 'yd: npp6 after alloc: bmphy_alloc(1,j,p) ',j,'nvm',p,'nphs', bmphy_alloc(1,j,p) 
                 ENDDO
                 WRITE(numout,*) 'yd: npp6 after alloc: bmFFB_alloc(1,j) ',j,'nvm', bmFFB_alloc(1,j,:) 
                 WRITE(numout,*) 'yd: npp6 after alloc: bmphy_alloc(1,j) ',j,'nvm', bmphy_alloc(1,j,:) 
                 !WRITE(numout,*) 'yd: npp6 after alloc:npts=',i,', bmFFB_alloc(:,j) ',j,'nvm', bmFFB_alloc(:,j,:) 
                 !WRITE(numout,*) 'yd: npp6 after alloc:npts=',i,', bmphy_alloc(:,j) ',j,'nvm', bmphy_alloc(:,j,:) 
              ENDIF
           ENDIF       

    !          bm_alloc_sapabold(:,j)=bm_alloc(:,j,isapabove,icarbon)
    !          bm_alloc(:,j,isapabove,icarbon)=biomass(:,j,isapabove,icarbon)-FFBbm(:,j)
    !          WRITE(numout,*) 'bm_alloc_sapabold(1,j) after npp_alloc: ',j,'pft, ',bm_alloc_sapabold(1,j)
    !          WRITE(numout,*) 'bm_alloc(1,j,:,icarbon) after npp_alloc minus FFB: ',j,'pft, ',bm_alloc(1,j,:,icarbon)
    !          WRITE(numout,*) 'FFBbm(1,j) after npp_alloc: ',j,'pft, ',FFBbm(1,j)
!! yidi        

!!xchen1111
!!           IF ((var_date .GT. 210) .AND. (var_date .LT. 330)) THEN                      
!!             bm_alloc(:,j,ileaf,icarbon) = 7.5*f_alloc(:,j,ileaf) * bm_alloc_tot(:,j) !!+ delta_turnover(:,j)
!!             bm_alloc(:,j,isapabove,icarbon) = f_alloc(:,j,isapabove) * bm_alloc_tot(:,j) !!-0.9*f_alloc(:,j,ileaf) * bm_alloc_tot(:,j)!!- delta_turnover(:,j)/3.0
!!             bm_alloc(:,j,iroot,icarbon) = f_alloc(:,j,iroot) * bm_alloc_tot(:,j) !!-0.4*f_alloc(:,j,ileaf) * bm_alloc_tot(:,j)!!- delta_turnover(:,j)/3.0
!!             bm_alloc(:,j,ifruit,icarbon) = f_alloc(:,j,ifruit) * bm_alloc_tot(:,j) !!- 0.2*f_alloc(:,j,ileaf) * bm_alloc_tot(:,j) !!- delta_turnover(:,j)/3.0
!!           ELSE
!!             bm_alloc(:,j,ileaf,icarbon) = 0.5*f_alloc(:,j,ileaf) * bm_alloc_tot(:,j) !!+ delta_turnover(:,j)
!!             bm_alloc(:,j,isapabove,icarbon) = f_alloc(:,j,isapabove) * bm_alloc_tot(:,j) !!- f_alloc(:,j,ileaf) * bm_alloc_tot(:,j)/6!!- delta_turnover(:,j)/3.0
!!             bm_alloc(:,j,iroot,icarbon) = f_alloc(:,j,iroot) * bm_alloc_tot(:,j) !!- f_alloc(:,j,ileaf) * bm_alloc_tot(:,j)/6!!- delta_turnover(:,j)/3.0
!!             bm_alloc(:,j,ifruit,icarbon) = f_alloc(:,j,ifruit) * bm_alloc_tot(:,j) !!- f_alloc(:,j,ileaf) * bm_alloc_tot(:,j)/6!!- delta_turnover(:,j)/3.0
!!           ENDIF
!!           DO k = 1, nparts                                                     
!!              f_alloc(:,j,k) =  bm_alloc(:,j,k,icarbon)/bm_alloc_tot(:,j)      
!!           ENDDO
       ENDIF
!       IF (j.EQ.10) WRITE(numout,*) 'zd bm_alloc 6','bm_alloc(1,10,:,icarbon)',bm_alloc(1,10,:,icarbon)

       
       !! 3.4 Calculate growth respiration of each plant compartment.
       !      Growth respiration of a plant compartment is a fraction of the allocatable biomass remaining after
       !      maintenance respiration losses have been taken into account. The fraction of allocatable biomass 
       !      removed for growth respiration is the same for all plant compartments and is defined by the parameter
       !      frac_growth_resp in stomate_constants.f90. Allocatable biomass ::bm_alloc is updated as a result of 
       !      the removal of growth resp.

       !!! xuhui: note that we should exclude those negative bm_alloc induced by crop remobilizations
!       WRITE(numout,*) 'bm_alloc(1,j,:,icarbon) after growth_resp: ',j ,'pft, ',bm_alloc(1,j,:,icarbon)
       DO i=1, npts
           resp_growth_part(i,:) = zero
           WHERE (bm_alloc(i,j,:,icarbon) .GT. zero)
               resp_growth_part(i,:) = frac_growthresp(j) * bm_alloc(i,j,:,icarbon) / dt
               bm_alloc(i,j,:,icarbon) = ( un - frac_growthresp(j) ) * bm_alloc(i,j,:,icarbon)
           ENDWHERE

!! yidi
        !IF ((j .GE. 2) .AND. (j .LE. 7)) THEN
        !   WRITE(numout,*) 'yd: npp7 bm_alloc(:,j,:,icarbon) after growthresp', bm_alloc(:,j,:,icarbon)
        !   WRITE(numout,*) 'yd: npp7 resp_growth_part(:)', resp_growth_part(:,j)
        !ENDIF
        IF (is_oilpalm(j)) THEN
           WRITE(numout,*) 'yd: npp7 bm_alloc(:,j,:,icarbon) after growthresp', bm_alloc(:,j,:,icarbon)
           WRITE(numout,*) 'yd: npp7 resp_growth_part(:)', resp_growth_part(:,j)
        ENDIF 
           IF (ok_oilpalm) THEN
              IF (is_oilpalm(j)) THEN
                  DO p = 1, nphs
                     IF (bm_alloc(i,j,isapabove,icarbon) .GT. zero) THEN
                        bmFFB_alloc(i,j,p) =  ( un - frac_growthresp(j) ) * bmFFB_alloc(i,j,p)
                        bmphy_alloc(i,j,p)  = ( un - frac_growthresp(j) ) * bmphy_alloc(i,j,p)
                     !WRITE(numout,*) 'yd:npp7 bmphy_alloc(1,j,p) after growthresp',j,'nvm',p,'nphs', bmphy_alloc(1,j,p) 
                     !WRITE(numout,*) 'yd:npp7 bmFFB_alloc(1,j,p) after growthresp',j,'nvm',p,'nphs', bmFFB_alloc(1,j,p) 
                     ENDIF
                  ENDDO
                  WRITE(numout,*) 'yd:npp7 bmphy_alloc(i,j,:) after growthresp',i,'npts,',j,'nvm', bmphy_alloc(i,j,:) 
                  WRITE(numout,*) 'yd:npp7 bmFFB_alloc(i,j,:) after growthresp',i,'npts,',j,'nvm', bmFFB_alloc(i,j,:) 
              ENDIF
           ENDIF
!! yidi
       ENDDO
    
!       WRITE(numout,*) 'bm_alloc(1,j,:,icarbon) after growth_resp: ',j ,'pft, ',bm_alloc(1,j,:,icarbon)
!       IF (j.EQ.10) WRITE(numout,*) 'zd bm_alloc 7','bm_alloc(1,10,:,icarbon)',bm_alloc(1,10,:,icarbon),'frac_growthresp(10)',frac_growthresp(10)
 
       !! 3.5 Total growth respiration 
       !      Calculate total growth respiration of the plant as the sum of growth respiration of all plant parts	
       resp_growth(:,j) = zero

       DO k = 1, nparts
          resp_growth(:,j) = resp_growth(:,j) + resp_growth_part(:,k)
       ENDDO

    ENDDO ! # End Loop over # of PFTs
    !WRITE(numout,*) 'zd nppcalc2 biomass(1,10,ileaf,icarbon)',biomass(1,10,ileaf,icarbon),'bm_pump(1)',bm_pump(1),'resp_maint_part(1,10,ileaf)',resp_maint_part(1,10,ileaf),'resp_maint(1,10)',resp_maint(1,10)
!    WRITE(numout,*) 'frac_growthresp: ', frac_growthresp
!    WRITE(numout,*) 'bm_alloc(1,12:14,:,icarbon) after update: ', bm_alloc(1,12:14,:,icarbon)
    
 !! 4. Update the biomass with newly allocated biomass after respiration
 
    !  Save the old leaf biomass for later. "old" leaf mass is leaf mass after maintenance respiration in the case 
    !  where maintenance respiration has required taking biomass from tissues in section 3.3 
    lm_old(:,:) = biomass(:,:,ileaf,icarbon)
    biomass(:,:,:,:) = biomass(:,:,:,:) + bm_alloc(:,:,:,:)

!! yidi
    DO  j = 1, nvm
      IF (ok_oilpalm) THEN
         IF (is_oilpalm(j)) THEN
            DO p = 1, nphs
  !!  bmFFB_alloc=0, FtoPHY=0 when is_oilpalm_ffbharvest(j) = FALSE
               bm_FFB(:,j,p) = bm_FFB(:,j,p) + bmFFB_alloc(:,j,p)  
               FFBbm(:,j) = FFBbm(:,j) + bmFFB_alloc(:,j,p)
               bm_phytomer(:,j,p) = bm_phytomer(:,j,p) + bmphy_alloc(:,j,p)
               PHYbm(:,j) = PHYbm(:,j) + bmphy_alloc(:,j,p)
              ! WRITE(numout,*) 'yd npp8: bm_phytomer(1,j,p) after npp update: ',j,'nvm',p,'nphs',bm_phytomer(1,j,:)
              ! WRITE(numout,*) 'yd npp8: bm_FFB(1,j,p) after npp update: ',j,'nvm',p,'nphs',bm_FFB(1,j,:)
            ENDDO
               WRITE(numout,*) 'yd npp8: bm_phytomer(1,j,:) after npp update: ',j,'nvm',bm_phytomer(1,j,:)
               WRITE(numout,*) 'yd npp8: bm_FFB(1,j,:) after npp update: ',j,'nvm',bm_FFB(1,j,:)
               WRITE(numout,*) 'yd npp8: PHYbm(1,j) after npp update: ',j,'nvm',PHYbm(1,j)
               WRITE(numout,*) 'yd npp8: FFBbm(1,j) after npp update: ',j,'nvm',FFBbm(1,j)
         ENDIF
      ENDIF
    ENDDO
!! yidi

    !WRITE(numout,*) 'yd: npp8: biomass(1,2:7,ileaf,icarbon) after bmalloc',biomass(1,2:7,ileaf,icarbon)
    WRITE(numout,*) 'yd: npp8: biomass(1,62:67,ileaf,icarbon) after bmalloc',biomass(1,62:67,ileaf,icarbon)
!BEGINNVADD
    npp_part (:,:,:)=  npp_part (:,:,:) + bm_alloc(:,:,:,icarbon)
!ENDNVADD
    !WRITE(numout,*) 'yd: npp8: biomass(1,2,:,icarbon) after update: ', biomass(1,2,:,icarbon)
    WRITE(numout,*) 'yd: npp8: biomass(1,62,:,icarbon) after update: ', biomass(1,62,:,icarbon)
 !! 5. Deal with negative biomasses
    
    !  Biomass can become negative in some rare cases, as the GPP can be negative. This corresponds to  
    !  situations that can be seen as the 'creation' of a seed ('virtual photosynthesis'). In this case, we set
    !  biomass to some small value min_stomate. For carbon budget to remain balanced, this creation of matter (carbon)
    !  is taken into account by decreasing the autotrophic respiration by the same amount that has been added to biomass
    !  for it to become positive. In this case, maintenance respiration can become negative in extreme cases (deserts)!!
    DO k = 1, nparts	! Loop over # of plant parts

       DO j = 2,nvm	! Loop over # of PFTs

          WHERE ( biomass(:,j,k,icarbon) .LT. zero )
            
             bm_create(:,j) = min_stomate - biomass(:,j,k,icarbon)

             ! Set biomass to min_stomate 
             biomass(:,j,k,icarbon) = biomass(:,j,k,icarbon) + bm_create(:,j)

             resp_maint(:,j) = resp_maint(:,j) - bm_create(:,j) / dt

          ENDWHERE

       ENDDO	! Loop over # of PFTs

    ENDDO	! Loop over # plant parts

!! yidi
    !WRITE(numout,*) 'yd:npp9 biomass(1,2,:,icarbon) after deal with negative biomass: ',biomass(1,2,:,icarbon)
    WRITE(numout,*) 'yd:npp9 biomass(1,62,:,icarbon) after deal with negative biomass: ',biomass(1,62,:,icarbon)
    DO j = 2,nvm ! Loop over # of PFTs 
       IF (ok_oilpalm) THEN
          IF (is_oilpalm(j)) THEN
             WRITE(numout,*) 'yd:npp9 negative before bm_FFB(1,j,1)',j,'nvm', bm_FFB(1,j,1) 
             WRITE(numout,*) 'yd:npp9 negative before bm_phytomer(1,j,1)',j,'nvm', bm_phytomer(1,j,1) 
             WRITE(numout,*) 'yd:npp9 negative before PHYbm(1,j)',j,'nvm', PHYbm(1,j) 
             WRITE(numout,*) 'yd:npp9 negative before FFBbm(1,j)',j,'nvm', FFBbm(1,j) 
             WHERE (bm_FFB(:,j,p) .LT. zero)
                bm_create(:,j) = min_stomate - bm_FFB(:,j,p)
                bm_FFB(:,j,p) = bm_FFB(:,j,p) + bm_create(:,j)
             ENDWHERE
                
             WHERE (bm_phytomer(:,j,p) .LT. zero)
                bm_create(:,j) = min_stomate - bm_phytomer(:,j,p)
                bm_phytomer(:,j,p) = bm_phytomer(:,j,p) + bm_create(:,j)
             ENDWHERE
             WHERE (PHYbm(:,j) .LT. zero)
                bm_create(:,j) = min_stomate - PHYbm(:,j)
                PHYbm(:,j) = PHYbm(:,j) + bm_create(:,j)
             ENDWHERE
             WHERE (FFBbm(:,j) .LT. zero)
                bm_create(:,j) = min_stomate - FFBbm(:,j)
                FFBbm(:,j) = FFBbm(:,j) + bm_create(:,j)
             ENDWHERE
              
          !   WRITE(numout,*) 'yd:npp10 negative after bm_FFB(1,j,1)',j,'nvm', bm_FFB(1,j,1) 
          !   WRITE(numout,*) 'yd:npp10 negative after bm_phytomer(1,j,1)',j,'nvm', bm_phytomer(1,j,1) 
          ENDIF
       ENDIF
    ENDDO ! Loop over # of PFTs
    WRITE(numout,*) 'yd:npp10 bm_phytomer(1,62:67,1) after deal with negative biomass: ',bm_phytomer(1,62:67,1)
    WRITE(numout,*) 'yd:npp10 bm_FFB(1,62:67,1) after deal with negative biomass: ',bm_FFB(1,62:67,1)
    DO j =62,nvm
       IF (biomass(1,j,isapabove,icarbon) .GT. min_stomate) THEN
          WRITE(numout,*) 'yd:npp10 bm_phytomer(1,j,:) after deal with negative biomass: ',bm_phytomer(1,j,:)
          WRITE(numout,*) 'yd:npp10 bm_FFB(1,j,:) after deal with negative biomass: ',bm_FFB(1,j,:)
       ENDIF
    ENDDO
    WRITE(numout,*) 'yd:npp10 PHYbm(1,62:67) after deal with negative biomass: ',PHYbm(1,62:67)
    WRITE(numout,*) 'yd:npp10 FFBbm(1,62:67) after deal with negative biomass: ',FFBbm(1,62:67)
    !WRITE(numout,*) 'yd:npp10 biomass(1,62:67,:,1) after deal with negative biomass: ',biomass(1,62:67,:,1)
    WRITE(numout,*) 'yd:npp10 biomass(1,62:67,isapabove,1) after : ',biomass(1,62:67,isapabove,1)
    !WRITE(numout,*) 'yd:npp10 biomass(1,2:7,isapabove,1) after : ',biomass(1,2:7,isapabove,1)
    !WRITE(numout,*) 'yd:npp10 biomass(1,2,:,1) after : ',biomass(1,2,:,1)
    WRITE(numout,*) 'yd:npp10 biomass(1,62,:,1) after : ',biomass(1,62,:,1)
    !WRITE(numout,*) 'yd:npp10 biomass(1,3,:,1) after : ',biomass(1,3,:,1)
    WRITE(numout,*) 'yd:npp10 biomass(1,63,:,1) after : ',biomass(1,63,:,1)
             
!! yidi

    !WRITE(numout,*) 'zd nppcalc4 biomass(1,10,ileaf,icarbon)',biomass(1,10,ileaf,icarbon)

!    WRITE(numout,*) 'biomass(1,12:14,:,icarbon) after negative removal: ', biomass(1,12:14,:,icarbon)
 !! 6. Calculate NPP (See Eq 1 in header)
    
    !  Calculate the NPP @tex $(gC.m^{-2}dt^{-1})$ @endtex as the difference between GPP
    !  and autotrophic respiration (maintenance and growth respirations)
    DO j = 2,nvm	! Loop over # PFTs
       npp(:,j) = gpp(:,j) - resp_growth(:,j) - resp_maint(:,j)
    ENDDO	! Loop over # PFTs

    WRITE(numout,*) 'npp(:,2)*dt', npp(:,2)*dt
    
    DO j =62,nvm
       WRITE(numout,*) 'yd: npp11 caluclate nvm=',j,',npp(:,j)*dt', npp(:,j)*dt
       WRITE(numout,*) 'yd: npp11 nvm=',j,',resp_maint(:,j)', resp_maint(:,j)
    ENDDO
    !DO j =2,7
    !   WRITE(numout,*) 'yd: npp11 caluclate nvm=',j,',npp(:,j)*dt', npp(:,j)*dt
    !   WRITE(numout,*) 'yd: npp11 nvm=',j,',resp_maint(:,j)', resp_maint(:,j)
    !   !WRITE(numout,*) 'yd: npp11 nvm=',j,',resp_maint_part(:,j,nparts)', resp_maint_part(:,j,:)
    !ENDDO
 !! 7. Update leaf age

    !  Leaf age is needed for calculation of turnover and vmax in stomate_turnover.f90 and stomate_vmax.f90 routines. 
    !  Leaf biomass is distributed according to its age into several "age classes" with age class=1 representing the
    !  youngest class, and consisting of the most newly allocated leaf biomass 
    
    !! 7.1 Update quantity and age of the leaf biomass in the youngest class
    !      The new amount of leaf biomass in the youngest age class (leaf_mass_young) is the sum of :
    !      - the leaf biomass that was already in the youngest age class (leaf_frac(:,j,1) * lm_old(:,j)) with the 
    !        leaf age given in leaf_age(:,j,1) 
    !      - and the new biomass allocated to leaves (bm_alloc(:,j,ileaf)) with a leaf age of zero.
    DO j = 2,nvm
       leaf_mass_young(:,j) = leaf_frac(:,j,1) * lm_old(:,j) + bm_alloc(:,j,ileaf,icarbon)
    ENDDO

    ! The age of the updated youngest age class is the average of the ages of its 2 components: bm_alloc(leaf) of age
    ! '0', and leaf_frac*lm_old(=leaf_mass_young-bm_alloc) of age 'leaf_age(:,j,1)' 
    DO j = 2,nvm
       WHERE ( ( bm_alloc(:,j,ileaf,icarbon) .GT. zero ) .AND. &
         ( leaf_mass_young(:,j) .GT. min_sechiba ) )

          leaf_age(:,j,1) = MAX ( zero, &
               & leaf_age(:,j,1) * &
               & ( leaf_mass_young(:,j) - bm_alloc(:,j,ileaf,icarbon) ) / &
               & leaf_mass_young(:,j) )
          
       ENDWHERE
    ENDDO
    !WRITE(numout,*) 'yd: npp12 leaf_age(:,2:7,1) before', leaf_age(:,2:7,1) !! yidi
    WRITE(numout,*) 'yd: npp12 leaf_age(:,62:67,1) before', leaf_age(:,62:67,1) !! yidi
    !! 7.2 Update leaf age
    !      Update fractions of leaf biomass in each age class (fraction in youngest class increases)

    !! 7.2.1 Update age of youngest leaves
    !        For age class 1 (youngest class), because we have added biomass to the youngest class, we need to update
    !        the fraction of total leaf biomass that belongs to the youngest age class : updated mass in class divided
    !        by new total leaf mass
    DO j = 2,nvm
       WHERE ( biomass(:,j,ileaf,icarbon) .GT. min_stomate )

          leaf_frac(:,j,1) = leaf_mass_young(:,j) / biomass(:,j,ileaf,icarbon)

       ENDWHERE
    ENDDO

    !! 7.2.2 Update age of other age classes
    !        Because the total leaf biomass has changed, we need to update the fraction of leaves in each age class:
    !        mass in leaf age class (from previous fraction of leaves in this class and previous total leaf biomass) 
    !        divided by new total mass
    DO m = 2, nleafages	! Loop over # leaf age classes

       DO j = 2,nvm	! Loop over # PFTs
          WHERE ( biomass(:,j,ileaf,icarbon) .GT. min_stomate )

             leaf_frac(:,j,m) = leaf_frac(:,j,m) * lm_old(:,j) / biomass(:,j,ileaf,icarbon)

          ENDWHERE
       ENDDO

    ENDDO	! Loop over # leaf age classes
!gmjc varied sla for managed grassland
    leaf_age_w = 0.0
    DO j = 2,nvm
      IF (is_grassland_manag(j)) THEN
         sla_max_Nfert(:,j)=sla_max(j)
         sla_min_Nfert(:,j)=sla_min(j)

      ELSE
         sla_max_Nfert(:,j)=sla_max(j)
         sla_min_Nfert(:,j)=sla_min(j)
      ENDIF
      IF (is_oilpalm(j)) THEN
        WRITE(numout,*) 'yd: npp12.5 sla_age1',sla_age1(:,j),'nvm=',j
        WRITE(numout,*) 'yd: npp12.5 sla_age2',sla_age2(:,j),'nvm=',j
        WRITE(numout,*) 'yd: npp12.5 sla_age3',sla_age3(:,j),'nvm=',j
        WRITE(numout,*) 'yd: npp12.5 sla_age4',sla_age4(:,j),'nvm=',j
        WRITE(numout,*) 'yd: npp12.5 sla_max(j)',sla_max(j)
        WRITE(numout,*) 'yd: npp12.5 sla_calc(:,j)',sla_calc(:,j)
      ENDIF

      WHERE ( ( bm_alloc(:,j,ileaf,icarbon) .GT. 0.0 ) .AND. &
             ( leaf_mass_young(:,j) .GT. 0.0 ) )

       sla_age1(:,j) = (sla_age1(:,j) * &
                       (leaf_mass_young(:,j)-bm_alloc(:,j,ileaf,icarbon)) + &
                       sla_max(j) * bm_alloc(:,j,ileaf,icarbon)) / leaf_mass_young(:,j)
       sla_age2(:,j) = sla_max(j)*0.9
       sla_age3(:,j) = sla_max(j)*0.85
       sla_age4(:,j) = sla_max(j)*0.8

       sla_calc(:,j) = sla_age1(:,j) * leaf_frac(:,j,1) + &
                       sla_age2(:,j) * leaf_frac(:,j,2) +  &
                       sla_age3(:,j) * leaf_frac(:,j,3) + & 
                       sla_age4(:,j) * leaf_frac(:,j,4)
      ENDWHERE

      IF (is_oilpalm(j)) THEN
        WRITE(numout,*) 'yd: npp12.75 sla_age1',sla_age1(:,j),'nvm=',j
        WRITE(numout,*) 'yd: npp12.75 sla_age2',sla_age2(:,j),'nvm=',j
        WRITE(numout,*) 'yd: npp12.75 sla_age3',sla_age3(:,j),'nvm=',j
        WRITE(numout,*) 'yd: npp12.75 sla_age4',sla_age4(:,j),'nvm=',j
        WRITE(numout,*) 'yd: npp12.75 sla_calc(:,j)',sla_calc(:,j),'nvm=',j
      ENDIF
      leaf_age_w(:,j) = 0.0
      DO m = 1, nleafages
        leaf_age_w(:,j) = leaf_age_w(:,j)+ leaf_age(:,j,m)*leaf_frac(:,j,m)
      END DO

      ! sla_calc can not be greater than sla_max or less than sla_min, and sla
      ! will be at maximum when age< 10
      WHERE (sla_calc(:,j) .GT. sla_max(j))
        sla_calc(:,j) = sla_max(j)
      ELSE WHERE (leaf_age_w(:,j) .LT. 5.0)
        sla_calc(:,j) = sla_max(j)
      ELSE WHERE (sla_calc(:,j) .LT. sla_min(j))
        sla_calc(:,j) = sla_min(j)
      ENDWHERE
    END DO
      IF (is_oilpalm(j)) THEN
        WRITE(numout,*) 'yd: npp12.85 sla_calc',sla_calc(:,j),'nvm=',j
      ENDIF
!end gmjc

 !! 8. Update whole-plant age 
    
    WRITE(numout,*) 'age(:,62:67) after update', age(:,62:67) !! yidi
    !! 8.1 PFT age
    !      At every time step, increase age of the biomass that was already present at previous time step. 
    !      Age is expressed in years, and the time step 'dt' in days so age increase is: dt divided by number 
    !      of days in a year.
    WHERE ( PFTpresent(:,:) )

       age(:,:) = age(:,:) + dt/one_year

    ELSEWHERE

       age(:,:) = zero

    ENDWHERE

    !WRITE(numout,*) 'yd: npp13 leaf_age(:,2:7,1) after', leaf_age(:,2:7,1) !! yidi
    WRITE(numout,*) 'yd: npp13 leaf_age(:,62:67,1) after', leaf_age(:,62:67,1) !! yidi
    !WRITE(numout,*) 'age(:,2:7) after update', age(:,2:7) !! yidi
    WRITE(numout,*) 'age(:,62:67) after update', age(:,62:67) !! yidi
    !! 8.2 Age of grasses and crops
    !  For grasses and crops, biomass with age 0 has been added to the whole plant with age 'age'. New biomass is the sum of 
    !  the current total biomass in all plant parts (bm_new), bm_new(:) = SUM( biomass(:,j,:), DIM=2 ). The biomass that has 
    !  just been added is the sum of the allocatable biomass of all plant parts (bm_add), its age is zero. bm_add(:) = 
    !  SUM( bm_alloc(:,j,:), DIM=2 ). Before allocation, the plant biomass is bm_new-bm_add, its age is "age(:,j)". The age of
    !  the new biomass is the average of the ages of previous and added biomass.
    !  For trees, age is treated in "establish" if vegetation is dynamic, and in turnover routines if it is static (in this 
    !  case, only the age of the heartwood is accounted for).
    DO j = 2,nvm

       IF ( .NOT. is_tree(j) ) THEN

          bm_new(:) = biomass(:,j,ileaf,icarbon) + biomass(:,j,isapabove,icarbon) + &
               biomass(:,j,iroot,icarbon) + biomass(:,j,ifruit,icarbon)
          bm_add(:) = bm_alloc(:,j,ileaf,icarbon) + bm_alloc(:,j,isapabove,icarbon) + &
               bm_alloc(:,j,iroot,icarbon) + bm_alloc(:,j,ifruit,icarbon)

          WHERE ( ( bm_new(:) .GT. min_sechiba ) .AND. ( bm_add(:) .GT. zero ) )
             age(:,j) = age(:,j) * ( bm_new(:) - bm_add(:) ) / bm_new(:)
          ENDWHERE

       ENDIF

    ENDDO

 !! 9. Write history files

    CALL xios_orchidee_send_field("BM_ALLOC_LEAF",bm_alloc(:,:,ileaf,icarbon))
    CALL xios_orchidee_send_field("BM_ALLOC_SAP_AB",bm_alloc(:,:,isapabove,icarbon))
    CALL xios_orchidee_send_field("BM_ALLOC_SAP_BE",bm_alloc(:,:,isapbelow,icarbon))
    CALL xios_orchidee_send_field("BM_ALLOC_ROOT",bm_alloc(:,:,iroot,icarbon))
    CALL xios_orchidee_send_field("BM_ALLOC_FRUIT",bm_alloc(:,:,ifruit,icarbon))
    CALL xios_orchidee_send_field("BM_ALLOC_RES",bm_alloc(:,:,icarbres,icarbon))
!!xchen add to check the balance of biomass for each parts
    !!CALL xios_orchidee_send_field("BM_LEAF_BAL",biomass(:,:,ileaf,icarbon))  
    !!CALL xios_orchidee_send_field("BM_SAP_BAL",biomass(:,:,isapabove,icarbon))
    !!CALL xios_orchidee_send_field("BM_SAP_BAL",biomass(:,:,isapbelow,icarbon))
    !!CALL xios_orchidee_send_field("BM_ROOT_BAL",biomass(:,:,iroot,icarbon))  
    !!CALL xios_orchidee_send_field("BM_FRUIT_BAL",biomass(:,:,ifruit,icarbon))
    !!CALL xios_orchidee_send_field("BM_RES_BAL",biomass(:,:,icarbres,icarbon))
    !!CALL xios_orchidee_send_field("BM_TOT_BAL",biomass(:,:,:,icarbon))   
!! yidi output wroung for test 
   ! CALL xios_orchidee_send_field("BM_ALLOC_SAPABOVE_PRE",bm_alloc_sapabold(:,:))
   ! CALL xios_orchidee_send_field("BM_ALLOC_SAPABOVE_PRIOR",bm_alloc(:,:,isapabove,icarbon))
    IF (ok_oilpalm) THEN
       CALL xios_orchidee_send_field("BM_ALLOC_PHYTOMER",bmphy_alloc(:,:,:))
       CALL xios_orchidee_send_field("BM_ALLOC_FFB",bmFFB_alloc(:,:,:))
       CALL xios_orchidee_send_field("BM_ALLOC_FFB_AGE1",bmFFB_alloc(:,62,:))
       CALL xios_orchidee_send_field("BM_ALLOC_PHYTOMER_AGE1",bmphy_alloc(:,62,:))
       !CALL xios_orchidee_send_field("BM_PHYEACH",bm_phytomer(:,:,:))
       !CALL xios_orchidee_send_field("BM_FFBEACH",bm_FFB(:,:,:))
       !CALL xios_orchidee_send_field("BM_FFB",FFBbm(:,:))
       !CALL xios_orchidee_send_field("BM_PHYTOMER",PHYbm(:,:))
    ENDIF
!! yidi
!gmjc
    npp_above(:,:) = npp_part(:,:,ileaf)+npp_part(:,:,isapabove)+&
                     npp_part(:,:,ifruit)+npp_part(:,:,icarbres)/2.
    npp_below(:,:) = npp_part(:,:,iroot)+npp_part(:,:,isapbelow)+&
                     npp_part(:,:,icarbres)/2.
    CALL xios_orchidee_send_field("SLA_CALC",sla_calc)
    CALL xios_orchidee_send_field("NPP_ABOVE",npp_above)
    CALL xios_orchidee_send_field("NPP_BELOW",npp_below)
    !!@xchen lai for leaf age class1-4 for caculating PC in Wu's science paper
    lai_cl1(:,:) = biomass(:,:,ileaf,icarbon)*leaf_frac(:,:,1)*sla_age1(:,:)
    lai_cl2(:,:) = biomass(:,:,ileaf,icarbon)*leaf_frac(:,:,2)*sla_age2(:,:)    
    lai_cl3(:,:) = biomass(:,:,ileaf,icarbon)*leaf_frac(:,:,3)*sla_age3(:,:)    
    lai_cl4(:,:) = biomass(:,:,ileaf,icarbon)*leaf_frac(:,:,4)*sla_age4(:,:)    
    !! yidi
    DO j=62,67
        WRITE(numout,*) 'yd: npp14 nvm=',j, 'leaf_frac',leaf_frac(:,j,1:4)
        WRITE(numout,*) 'yd: npp14 nvm=',j, 'biomass(:,:,ileaf,icarbon)',biomass(:,:,ileaf,icarbon)
        WRITE(numout,*) 'yd: npp14 sla_age1',sla_age1(:,j),'nvm=',j
        WRITE(numout,*) 'yd: npp14 sla_age2',sla_age2(:,j),'nvm=',j
        WRITE(numout,*) 'yd: npp14 sla_age3',sla_age3(:,j),'nvm=',j
        WRITE(numout,*) 'yd: npp14 sla_age4',sla_age4(:,j),'nvm=',j
        WRITE(numout,*) 'yd: npp14 lai_cl4', lai_cl4(:,j),'nvm=',j
        WRITE(numout,*) 'yd: npp14 lai_cl3', lai_cl3(:,j),'nvm=',j
        WRITE(numout,*) 'yd: npp14 lai_cl2', lai_cl2(:,j),'nvm=',j
        WRITE(numout,*) 'yd: npp14 lai_cl1', lai_cl1(:,j),'nvm=',j
    ENDDO
    !! yidi
    !!CALL xios_orchidee_send_field("LAI_CL1",biomass(:,:,ileaf,icarbon)*leaf_frac(:,:,1)*sla_age1(:,:))
    !!CALL xios_orchidee_send_field("LAI_CL2",biomass(:,:,ileaf,icarbon)*leaf_frac(:,:,2)*sla_age2(:,:))           
    !!CALL xios_orchidee_send_field("LAI_CL3",biomass(:,:,ileaf,icarbon)*leaf_frac(:,:,3)*sla_age3(:,:))                  
    !!CALL xios_orchidee_send_field("LAI_CL4",biomass(:,:,ileaf,icarbon)*leaf_frac(:,:,4)*sla_age4(:,:))
    !!CALL xios_orchidee_send_field("LAI_CL1",lai_cl1(:,:))
    !!CALL xios_orchidee_send_field("LAI_CL2",lai_cl2(:,:))
    !!CALL xios_orchidee_send_field("LAI_CL3",lai_cl3(:,:))
    !!CALL xios_orchidee_send_field("LAI_CL4",lai_cl4(:,:))

    CALL xios_orchidee_send_field("LEAF_BM_CL1",biomass(:,:,ileaf,icarbon)*leaf_frac(:,:,1))
    CALL xios_orchidee_send_field("LEAF_BM_CL2",biomass(:,:,ileaf,icarbon)*leaf_frac(:,:,2))
    CALL xios_orchidee_send_field("LEAF_BM_CL3",biomass(:,:,ileaf,icarbon)*leaf_frac(:,:,3))
    CALL xios_orchidee_send_field("LEAF_BM_CL4",biomass(:,:,ileaf,icarbon)*leaf_frac(:,:,4))
    CALL histwrite_p (hist_id_stomate, 'LEAF_BM_CL1', itime, &
         biomass(:,:,ileaf,icarbon)*leaf_frac(:,:,1), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LEAF_BM_CL2', itime, &
         biomass(:,:,ileaf,icarbon)*leaf_frac(:,:,2), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LEAF_BM_CL3', itime, &
         biomass(:,:,ileaf,icarbon)*leaf_frac(:,:,3), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LEAF_BM_CL4', itime, &
         biomass(:,:,ileaf,icarbon)*leaf_frac(:,:,4), npts*nvm, horipft_index)
    
   
    


    CALL xios_orchidee_send_field("BM_LEAF_BAL",biomass(:,:,ileaf,icarbon)) 
    CALL xios_orchidee_send_field("BM_SAP_AB_BAL",biomass(:,:,isapabove,icarbon))  
    CALL xios_orchidee_send_field("BM_SAP_BE_BAL",biomass(:,:,isapbelow,icarbon))  
    CALL xios_orchidee_send_field("BM_ROOT_BAL",biomass(:,:,iroot,icarbon))     
    CALL xios_orchidee_send_field("BM_FRUIT_BAL",biomass(:,:,ifruit,icarbon))   
    CALL xios_orchidee_send_field("BM_RES_BAL",biomass(:,:,icarbres,icarbon))
    !!CALL xios_orchidee_send_field("BM_TOT_BAL",biomass(:,:,:,icarbon))        

!end gmjc

    ! Save in history file the variables describing the biomass allocated to the plant parts
    CALL histwrite_p (hist_id_stomate, 'SLA_CROP',itime, &
         ssla, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'BM_ALLOC_LEAF', itime, &
         bm_alloc(:,:,ileaf,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'BM_ALLOC_SAP_AB', itime, &
         bm_alloc(:,:,isapabove,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'BM_ALLOC_SAP_BE', itime, &
         bm_alloc(:,:,isapbelow,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'BM_ALLOC_ROOT', itime, &
         bm_alloc(:,:,iroot,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'BM_ALLOC_FRUIT', itime, &
         bm_alloc(:,:,ifruit,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'BM_ALLOC_RES', itime, &
         bm_alloc(:,:,icarbres,icarbon), npts*nvm, horipft_index)
!gmjc

!! yidi !! wrong bm_alloc_phytomer in lpj !! yidi output
    IF (ok_oilpalm) THEN
       DO jv = 1, nvm
          IF (is_oilpalm(jv)) THEN  !don't bother to write if there are pfts that don't exist
             WRITE(part_str,'(I2)') jv
             IF (jv < 10) part_str(1:1) = '0'
             CALL histwrite_p(hist_id_stomate,'BM_ALLOC_PHY_'//part_str(1:LEN_TRIM(part_str)), &
                     itime, bm_phytomer(:,jv,:), npts*nphs, horiphytomer_index)
             CALL histwrite_p(hist_id_stomate,'BM_ALLOC_FFB_'//part_str(1:LEN_TRIM(part_str)), &
                     itime, bm_FFB(:,jv,:), npts*nphs, horiphytomer_index)
          ENDIF
       ENDDO
    ENDIF
!! yidi
    CALL histwrite_p(hist_id_stomate, 'SLA_CALC', itime, &
                    sla_calc(:,:), npts*nvm, horipft_index)

!end gmjc
!BEGINNVADD
    CALL histwrite_p(hist_id_stomate, 'NPP_ABOVE', itime, &
                    npp_above, npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate, 'NPP_BELOW', itime, &
                   npp_below, npts*nvm, horipft_index)
!ENDNVADD
    IF (printlev>=3) WRITE(numout,*) 'Leaving npp########################'

  END SUBROUTINE npp_calc

END MODULE stomate_npp
