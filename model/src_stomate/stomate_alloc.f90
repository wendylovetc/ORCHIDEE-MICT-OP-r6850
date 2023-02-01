! =================================================================================================================================
! MODULE       : stomate_alloc
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
!                This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Allocate net primary production to: carbon reserves, aboveground sapwood,
!! belowground sapwood, root, fruits and leaves.      
!!
!!\n DESCRIPTION: None
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

MODULE stomate_alloc

  ! Modules used:

  USE ioipsl_para
  USE pft_parameters
  USE stomate_data
  USE constantes
  USE constantes_soil

  IMPLICIT NONE

  ! Private & public routines

  PRIVATE
  PUBLIC alloc,alloc_clear

 ! Variables shared by all subroutines in this module

  LOGICAL, SAVE                                             :: firstcall_alloc = .TRUE.  !! Is this the first call? (true/false) 
!$OMP THREADPRIVATE(firstcall_alloc)
CONTAINS


!! ================================================================================================================================
!! SUBROUTINE   : alloc_clear
!!
!>\BRIEF          Set the flag ::firstcall_alloc to .TRUE. and as such activate section 
!! 1.1 of the subroutine alloc (see below).\n
!!
!_ ================================================================================================================================

  SUBROUTINE alloc_clear
    firstcall_alloc = .TRUE.
  END SUBROUTINE alloc_clear



!! ================================================================================================================================
!! SUBROUTINE 	: alloc
!!
!>\BRIEF         Allocate net primary production (= photosynthesis
!! minus autothrophic respiration) to: carbon reserves, aboveground sapwood,
!! belowground sapwood, root, fruits and leaves following Friedlingstein et al. (1999).
!!
!! DESCRIPTION (definitions, functional, design, flags):\n
!! The philosophy underlying the scheme is that allocation patterns result from
!! evolved responses that adjust carbon investments to facilitate capture of most
!! limiting resources i.e. light, water and mineral nitrogen. The implemented scheme 
!! calculates the limitation of light, water and nitrogen. However, nitrogen is not a 
!! prognostic variable of the model and therefore soil temperature and soil moisture 
!! are used as a proxy for soil nitrogen availability.\n
!! Sharpe & Rykiel (1991) proposed a generic relationship between the allocation of 
!! carbon to a given plant compartment and the availability of a particular resource:\n
!! \latexonly 
!!   \input{alloc1.tex}
!! \endlatexonly
!! \n
!! where A is the allocation of biomass production (NPP) to a given compartment (either 
!! leaves, stem, or roots). Xi and Yj are resource availabilities (e.g. light, water, 
!! nutrient). For a given plant compartment, a resource can be of type X or Y. An increase 
!! in a X-type resource will increase the allocation to compartment A. An increase in a 
!! Y-type resource will, however, lead to a decrease in carbon allocation to that compartment. 
!! In other words, Y-type resources are those for which uptake increases with increased 
!! investment in the compartment in question. X-type resources, as a consequence of 
!! trade-offs, are the opposite. For example, water is a Y-type resource for root allocation. 
!! Water-limited conditions should promote carbon allocation to roots, which enhance water 
!! uptake and hence minimize plant water stress. Negative relationships between investment 
!! and uptake arise when increased investment in one compartment leads, as required for 
!! conservation of mass, to decreased investment in a component involved in uptake of 
!! that resource.\n
!!
!! The implemented scheme allocates carbon to the following components:\n
!! - Carbon reserves;\n
!! - Aboveground sapwood;\n
!! - Belowground sapwood;\n
!! - Roots;\n
!! - Fruits/seeds and\n
!! - Leaves.
!! \n
!!
!! The allocation to fruits and seeds is simply a 10% "tax" of the total biomass 
!! production.\n
!! Following carbohydrate use to support budburst and initial growth, the 
!! carbohydrate reserve is refilled. The daily amount of carbon allocated to the 
!! reserve pool is proportional to leaf+root allocation (::LtoLSR and ::RtoLSR).\n
!! Sapwood and root allocation (respectively ::StoLSR and ::RtoLSR) are proportional 
!! to the estimated light and soil (water and nitrogen) stress (::Limit_L and 
!! ::Limit_NtoW). Further, Sapwood allocation is separated in belowground sapwood 
!! and aboveground sapwood making use of the parameter (:: alloc_sap_above_tree
!! or ::alloc_sap_above_grass). For trees partitioning between above and 
!! belowground compartments is a function of PFT age.\n 
!! Leaf allocation (::LtoLSR) is calculated as the residual of root and sapwood 
!! allocation (LtoLSR(:) = 1. - RtoLSR(:) - StoLSR(:).\n
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): :: f_alloc; fraction of NPP that is allocated to the 
!! six different biomass compartments (leaves, roots, above and belowground wood, 
!! carbohydrate reserves and fruits). DIMENSION(npts,nvm,nparts).
!!
!! REFERENCE(S)	:
!! - Friedlingstein, P., G. Joel, C.B. Field, and Y. Fung (1999), Towards an allocation
!! scheme for global terrestrial carbon models, Global Change Biology, 5, 755-770.\n
!! - Sharpe, P.J.H., and Rykiel, E.J. (1991), Modelling integrated response of plants 
!! to multiple stresses. In: Response of Plants to Multiple Stresses (eds Mooney, H.A., 
!! Winner, W.E., Pell, E.J.), pp. 205-224, Academic Press, San Diego, CA.\n
!! - Krinner G, Viovy N, de Noblet-Ducoudr N, Ogee J, Polcher J, Friedlingstein P,
!! Ciais P, Sitch S, Prentice I C (2005) A dynamic global vegetation model for studies
!! of the coupled atmosphere-biosphere system. Global Biogeochemical Cycles, 19, GB1015,
!! doi: 10.1029/2003GB002199.\n
!! - Malhi, Y., Doughty, C., and Galbraith, D. (2011). The allocation of ecosystem net primary productivity in tropical forests, 
!! Philosophical Transactions of the Royal Society B-Biological Sciences, 366, 3225-3245, DOI 10.1098/rstb.2011.0062.\n
!!
!! FLOWCHART    : 
!! \latexonly 
!!   \includegraphics[scale=0.5]{allocflow.jpg}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE alloc (npts, dt, &
       lai,last_lai,sup_flag, veget_cov_max, senescence, when_growthinit, &
       moiavail_week, tsoil_month, soilhum_month, &
       biomass, age, leaf_age, leaf_frac,rprof, f_alloc, &
       alloc_phy, FtoPHY, phytomer_age, & !! yidi
       deltai, ssla, & !added for crops, xuhui
!gmjc
       sla_calc, when_growthinit_cut, &
!!xchen1111
       delta_turnover, var_date, VPD, swdown,swdown_week, LtoLSR, StoLSR, RtoLSR)
!end gmjc
 !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                                 :: npts			!! Domain size - number of grid cells 
                                                                                        !! (unitless)
    REAL(r_std), INTENT(in)                                    :: dt			!! Time step of the simulations for stomate
                                                                                        !! (days)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: lai			!! PFT leaf area index 
                                                                                        !! @tex $(m^2 m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: veget_cov_max		!! PFT "Maximal" coverage fraction of a PFT 
                                                                                        !! (= ind*cn_ind) 
                                                                                        !! @tex $(m^2 m^{-2})$ @endtex
    LOGICAL, DIMENSION(npts,nvm), INTENT(in)                   :: senescence		!! Is the PFT senescent?  - only for 
                                                                                        !! deciduous trees (true/false)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: when_growthinit	!! Days since beginning of growing season 
                                                                                        !! (days)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: moiavail_week		!! PFT moisture availability - integrated
                                                                                        !! over a week (0-1, unitless) 
    REAL(r_std), DIMENSION(npts,nslm), INTENT(in)              :: tsoil_month		!! PFT soil temperature - integrated over 
                                                                                        !! a month (K) 
    REAL(r_std), DIMENSION(npts,nslm), INTENT(in)              :: soilhum_month		!! PFT soil humidity - integrated over a 
                                                                                        !! month (0-1, unitless) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: age			!! PFT age (days)
    !  for STICS
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: deltai
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: ssla
    ! end for STICS
!! xchen add vpd for gpp simulation in Amazon                                   
    REAL(r_std), DIMENSION(npts),     INTENT (in)       :: VPD              !! Vapor Pressure Deficit (kPa)
    REAL(r_std), DIMENSION(npts), INTENT(IN)        :: swdown               
    REAL(r_std), DIMENSION(npts), INTENT(IN)        :: swdown_week
    REAL(r_std), DIMENSION(npts), INTENT(out)                               :: LtoLSR                !! Ratio between leaf-allocation and 
                                                                                        !! (leaf+sapwood+root)-allocation 
                                                                                        !! (0-1, unitless)
    REAL(r_std), DIMENSION(npts), INTENT(out)                               :: StoLSR                !! Ratio between sapwood-allocation and 
                                                                                        !! (leaf+sapwood+root)-allocation 
                                                                                        !! (0-1, unitless)
    REAL(r_std), DIMENSION(npts), INTENT(out)                               :: RtoLSR                !! Ratio between root-allocation and 
                                                                                        !! (leaf+sapwood+root)-allocation 
                                                                                        !! (0-1, unitless)
    !! 0.2 Output variables

    !! 0.3 Modified variables
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)               :: last_lai		!! PFT leaf area index, last time step 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)               :: sup_flag		!! flag of whether use sup_LtoLSR 

    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: biomass		!! PFT total biomass 
                                                                                        !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nleafages), INTENT(inout)  :: leaf_age		!! PFT age of different leaf classes 
                                                                                        !! (days)
    REAL(r_std), DIMENSION(npts,nvm,nleafages), INTENT(inout)  :: leaf_frac		!! PFT fraction of leaves in leaf age class
                                                                                        !! (0-1, unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: rprof
!! [DISPENSABLE] PFT rooting depth - not
                                                                                        !!
                                                                                        !calculated
                                                                                        !in
                                                                                        !the
                                                                                        !current
                                                                                        !version
                                                                                        !of
                                                                                        !!
                                                                                        !the
                                                                                        !model
                                                                                        !(m)

    REAL(r_std), DIMENSION(npts,nvm,nparts), INTENT(out)       :: f_alloc		!! PFT fraction of NPP that is allocated to
                                                                                        !! the different components (0-1, unitless)
!!yidi for OP
    REAL(r_std), DIMENSION(npts,nvm,nphs), INTENT(out)         :: alloc_phy		!! Fraction of sapabove allocatable 
                                                                                        !! biomass that goes into 
                                                                                        !! each phytomer (unitless) 

    REAL(r_std), DIMENSION(npts,nvm,nphs), INTENT(in)          :: phytomer_age          !! age of phytomers (days)

                                                                                             
    REAL(r_std), DIMENSION(npts,nvm,nphs), INTENT(out)         :: FtoPHY                !! Fraction of allocatable 
                                                                                        !! biomass for each photomer
                                                                                        !! that goes into each FFB (unitless)
   ! REAL(r_std), DIMENSION(npts,nvm), INTENT(out)              :: FtoPHY_sum            !! sum of phytomerFFB-allocation 
!! yidi

    !! local
!gmjc
    ! how many days ago was the beginning of the cut
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)        :: when_growthinit_cut
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)     :: sla_calc
!!xchen1111                                                                     
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: delta_turnover            !! Intermediate variable for turnover ??
    INTEGER(i_std),INTENT(in)                            :: var_date               !! xchen Date (days) 

!end gmjc
    !! 0.4 Local variables
    REAL(r_std), DIMENSION(nvm)                                :: lai_happy		!! Lai threshold below which carbohydrate 
                                                                                        !! reserve may be used 
                                                                                        !! @tex $(m^2 m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts)                               :: limit_L		!! Lights stress (0-1, unitless)
    REAL(r_std), DIMENSION(npts)                               :: limit_N		!! Total nitrogen stress (0-1, unitless)
    REAL(r_std), DIMENSION(npts)                               :: limit_N_temp		!! Stress from soil temperature on nitrogen
                                                                                        !! mineralisation (0-1, unitless)
    REAL(r_std), DIMENSION(npts)                               :: limit_N_hum		!! Stress from soil humidity on nitrogen 
                                                                                        !! mineralisation (0-1, unitless)
    REAL(r_std), DIMENSION(npts)                               :: limit_W		!! Soil water stress (0-1, unitless)
    REAL(r_std), DIMENSION(npts)                               :: limit_WorN		!! Most limiting factor in the soil: 
                                                                                        !! nitrogen or water (0-1, unitless)
    REAL(r_std), DIMENSION(npts)                               :: limit			!! Most limiting factor: amongst limit_N, 
                                                                                        !! limit_W and limit_L (0-1, unitless)
    REAL(r_std), DIMENSION(npts)                               :: t_nitrogen		!! Preliminairy soil temperature stress 
                                                                                        !! used as a proxy for nitrogen stress (K)
    REAL(r_std), DIMENSION(npts)                               :: h_nitrogen		!! Preliminairy soil humidity stress used 
                                                                                        !! as a proxy for nitrogen stress 
                                                                                        !! (unitless)  
    REAL(r_std), DIMENSION(npts)                               :: rpc			!! Scaling factor for integrating vertical
                                                                                        !!  soil profiles (unitless) 	
    !!REAL(r_std), DIMENSION(npts)                               :: LtoLSR		!! Ratio between leaf-allocation and 
                                                                                        !! (leaf+sapwood+root)-allocation 
                                                                                        !! (0-1, unitless)
    !!REAL(r_std), DIMENSION(npts)                               :: StoLSR		!! Ratio between sapwood-allocation and 
                                                                                        !! (leaf+sapwood+root)-allocation 
                                                                                        !! (0-1, unitless)
    !!REAL(r_std), DIMENSION(npts)                               :: RtoLSR		!! Ratio between root-allocation and 
                                                                                        !! (leaf+sapwood+root)-allocation 
                                                                                        !! (0-1, unitless)
    REAL(r_std), DIMENSION(npts)                               :: carb_rescale		!! Rescaling factor for allocation factors
                                                                                        !! if carbon is allocated to carbohydrate 
                                                                                        !! reserve (0-1, unitless)
    REAL(r_std), DIMENSION(npts)                               :: use_reserve		!! Mass of carbohydrate reserve used to 
                                                                                        !! support growth
!!xchen1111
    REAL(r_std), DIMENSION(npts)                               :: var_use_reserve           !! Mass of carbohydrate reserve used to  
                                                                                        !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts)                               :: transloc_leaf		!! Fraction of carbohydrate reserve used 
                                                                                        !! (::use_reserve) to support leaf growth 
                                                                                        !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts)                               :: leaf_mass_young	!! Leaf biomass in youngest leaf age class 
                                                                                        !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)                           :: lm_old		!! Variable to store leaf biomass from 
                                                                                        !! previous time step 
                                                                                        !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std)                                                :: reserve_time		!! Maximum number of days during which 
                                                                                        !! carbohydrate reserve may be used (days)
    REAL(r_std), DIMENSION(npts,nvm)                           :: lai_around		!! lai on natural part of the grid cell, or
                                                                                        !! of agricultural PFTs 
                                                                                        !! @tex $(m^2 m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)                           :: veget_cov_max_nat 	!! Vegetation cover of natural PFTs on the 
                                                                                        !! grid cell (agriculture masked) 
                                                                                        !! (0-1, unitless)
    REAL(r_std), DIMENSION(npts)                               :: natveg_tot		!! Total natural vegetation cover on 
                                                                                        !! natural part of the grid cell 
                                                                                        !! (0-1, unitless)
    REAL(r_std), DIMENSION(npts)                               :: lai_nat		!! Average LAI on natural part of the grid
                                                                                        !! cell @tex $(m^2 m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts)                               :: zdiff_min		!! [DISPENSABLE] intermediate array for 
                                                                                        !! looking for minimum
    REAL(r_std), DIMENSION(npts)                               :: alloc_sap_above	!! Prescribed fraction of sapwood 
                                                                                        !! allocation to above ground sapwood 
                                                                                        !! (0-1, unitless)
    REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: z_soil                !! Variable to store depth of the different
                                                                                        !! soil layers (m)
!$OMP THREADPRIVATE(z_soil)
    INTEGER(i_std)                                             :: i,j,l,m,p            	!! Indices (unitless)
!gmjc
    ! rescaling factor for LtoLSR of grass,PFT=10
    REAL(r_std), DIMENSION(npts)                               :: alloc_leaf
    REAL(r_std), DIMENSION(npts)                               :: old_LtoLSR     !!xchen2017 Ratio between leaf-allocation and 
    REAL(r_std), DIMENSION(npts)                               :: xzLtoLSR       !!xchen2017 Ratio between leaf-allocation and 
    REAL(r_std), DIMENSION(npts)                               :: old_RtoLSR     !! 
    REAL(r_std), DIMENSION(npts)                               :: old_StoLSR     !! 
    REAL(r_std), DIMENSION(npts)                               :: new_StoLSR     !! 
    REAL(r_std), DIMENSION(npts)                               :: new_LtoLSR     !! 
    REAL(r_std), DIMENSION(npts)                               :: new_LtoLSR2     !! 
    REAL(r_std), DIMENSION(npts)                               :: new_RtoLSR     !! 
    REAL(r_std), DIMENSION(npts)                               :: new_RtoLSR2     !! 
    REAL(r_std), DIMENSION(npts)                               :: new_StoLSR2     !! 
!end gmjc
    INTEGER(i_std)                                             :: ier                   !! Error handling
    REAL(r_std), DIMENSION(npts,nvm)             :: lai_cl1          !! @xchen to add lai of each leaf class age and for                                                                                  !!calculating PC in wu's science paper
    REAL(r_std), DIMENSION(npts,nvm)             :: lai_cl2                     
    REAL(r_std), DIMENSION(npts,nvm)             :: lai_cl3                     
    REAL(r_std), DIMENSION(npts,nvm)             :: lai_cl4                     
    REAL(r_std), DIMENSION(npts,nvm)                      :: sla_age2           
    REAL(r_std), DIMENSION(npts,nvm)                      :: sla_age3           
    REAL(r_std), DIMENSION(npts,nvm)                      :: sla_age4   
!! yidi
    REAL(r_std), DIMENSION(npts,nvm)              :: ffbharvest_age_crit        !! critical ffbharvest age (days)
    REAL(r_std), DIMENSION(npts,nvm)              :: phytomer_age_crit          !! critical phytomer age (days)
!! yidi
!_ ================================================================================================================================

    IF (printlev>=3) WRITE(numout,*) 'Entering alloc###########'

!! 1. Initialize

    !! 1.1 First call only
    IF ( firstcall_alloc ) THEN

       !
       ! 1.1.0 Initialization
       !
       L0(2:nvm) = un - R0(2:nvm) - S0(2:nvm)  
       IF ((MINVAL(L0(2:nvm)) .LT. zero) .OR. (MAXVAL(S0(2:nvm)) .EQ. un)) THEN
          CALL ipslerr_p (3,'in module stomate_alloc', &
               &           'Something wrong happened', &
               &           'L0 negative or division by zero if S0 = 1', &
               &           '(Check your parameters.)')
       ENDIF
       

       !! 1.1.1 Copy the depth of the different soil layers (number of layers=nslm) 
       !        previously calculated as variable diaglev in routines sechiba.f90 and slowproc.f90  
       ALLOCATE(z_soil(0:nslm), stat=ier)
       IF ( ier /= 0 ) CALL ipslerr_p(3,'stomate_alloc','Pb in allocate of z_soil','','')
       z_soil(0) = zero
       z_soil(1:nslm) = diaglev(1:nslm)

       !! 1.1.2 Print flags and parameter settings
       IF (printlev >= 2) THEN
          WRITE(numout,*) 'alloc:'
          WRITE(numout,'(a,$)') '    > We'
          IF ( .NOT. ok_minres ) WRITE(numout,'(a,$)') ' do NOT'
          WRITE(numout,*) 'try to reach a minumum reservoir when severely stressed.'
          WRITE(numout,*) '   > Time delay (days) to build leaf mass (::tau_leafinit): ', &
               tau_leafinit(:)
          WRITE(numout,*) '   > Curvature of root mass with increasing soil depth (::z_nitrogen): ', &
               z_nitrogen
          WRITE(numout,*) '   > Sap allocation above the ground / total sap allocation (0-1, unitless): '
          WRITE(numout,*) '       grasses (::alloc_sap_above_grass) :', alloc_sap_above_grass
          WRITE(numout,*) '   > Default root alloc fraction (1; ::R0): ', R0(:)
          WRITE(numout,*) '   > Default sapwood alloc fraction (1; ::S0): ', S0(:)
          WRITE(numout,*) '   > Default fruit allocation (1, ::f_fruit): ', f_fruit
          WRITE(numout,*) '   > Minimum (min_LtoLSR)/maximum (::max_LtoLSR)leaf alloc fraction (0-1, unitless): ',&
               min_LtoLSR,max_LtoLSR
          WRITE(numout,*) '   > Maximum time (days) the carbon reserve can be used:'
          WRITE(numout,*) '       trees (reserve_time_tree):',reserve_time_tree
          WRITE(numout,*) '       grasses (reserve_time_grass):',reserve_time_grass
       END IF
       firstcall_alloc = .FALSE.

    ENDIF


    !! 1.2 Every call
    !! 1.2.1 Reset output variable (::f_alloc)
    f_alloc(:,:,:) = zero
    f_alloc(:,:,icarbres) = un
    alloc_phy(:,:,:) = zero
    FtoPHY(:,:,:) = zero

 
    !! 1.2.2 Proxy for soil nitrogen stress
    !        Nitrogen availability and thus N-stress can not be calculated by the model. Water and
    !        temperature stress are used as proxy under the assumption that microbial activity is 
    !        determined by soil temperature and water availability. In turn, microbial activity is 
    !         assumed to be an indicator for nitrogen mineralisation and thus its availability. 

    !! 1.2.2.1 Convolution of nitrogen stress with root profile 
    !          Here we calculate preliminary soil temperature and soil humidity stresses that will be used 
    !          as proxies for nitrogen stress. Their calculation follows the nitrogen-uptake capacity of roots.
    !          The capacity of roots to take up nitrogen is assumed to decrease exponentially with 
    !          increasing soil depth. The curvature of the exponential function describing the 
    !          nitrogen-uptake capacity of roots (= root mass * uptake capacity) is given by 
    !          ::z_nitrogen. Strictly speaking its unit is meters (m). Despite its units this parameter 
    !          has no physical meaning.
    !          Because the roots are described by an exponential function but the soil depth is limited to
    !          ::z_soil(nslm), the root profile is truncated at ::z_soil(nslm). For numerical reasons, 
    !          the total capacity of the soil profile for nitrogen uptake should be 1. To this aim a scaling
    !          factor (::rpc) is calculated as follows:\n
    !          \latexonly
    !            \input{alloc2.tex} 
    !          \endlatexonly
    !          Then temperature (::t_nitrogen) and humidity (::h_nitrogen) proxies for nitrogen stress are 
    !          calculated using mean weighted (weighted by nitrogen uptake capacity) soil temperature (::tsoil_month)
    !          or soil moisture (::soil_hum_month) (calculated in stomate_season.f90). 
    !          \latexonly
    !            \input{alloc3.tex} 
    !          \endlatexonly
    !          \latexonly
    !            \input{alloc4.tex} 
    !          \endlatexonly	
    !          \n 
		  
    ! Scaling factor for integration
    rpc(:) = un / ( un - EXP( -z_soil(nslm) / z_nitrogen ) )

    ! Integrate over # soil layers
    t_nitrogen(:) = zero

    DO l = 1, nslm ! Loop over # soil layers

       t_nitrogen(:) = &
            t_nitrogen(:) + tsoil_month(:,l) * rpc(:) * &
            ( EXP( -z_soil(l-1)/z_nitrogen ) - EXP( -z_soil(l)/z_nitrogen ) )

    ENDDO ! Loop over # soil layers

 
!!$    !! 1.2.2.2 Convolution for soil moisture
!!$    !          Scaling factor for integration
!!$    rpc(:) = 1. / ( 1. - EXP( -z_soil(nslm) / z_nitrogen ) )

    ! Integrate over # soil layers
    h_nitrogen(:) = zero

    DO l = 1, nslm ! Loop over # soil layers

       h_nitrogen(:) = &
            h_nitrogen(:) + soilhum_month(:,l) * rpc(:) * &
            ( EXP( -z_soil(l-1)/z_nitrogen ) - EXP( -z_soil(l)/z_nitrogen ) )

    ENDDO ! Loop over # soil layers


    !! 1.2.3 Separate between natural and agrigultural LAI
    !        The model distinguishes different natural PFTs but does not contain information
    !        on whether these PFTs are spatially separated or mixed. In line with the DGVM the 
    !        models treats the natural PFT's as mixed. Therefore, the average LAI over the
    !        natural PFTs is calculated to estimate light stress. Agricultural PFTs are spatially
    !        separated. 
    natveg_tot(:) = zero
    lai_nat(:) = zero

    DO j = 2, nvm ! Loop over # PFTs

       !WRITE(numout,*) 'yd:alloc0_vegetcovmax nvm=',j
       IF ( natural(j) .AND. .NOT. pasture(j)) THEN
          ! Mask agricultural vegetation
          veget_cov_max_nat(:,j) = veget_cov_max(:,j)
       ELSE
          ! Mask natural vegetation
          veget_cov_max_nat(:,j) = zero
       ENDIF

       ! Sum up fraction of natural space covered by vegetation
       natveg_tot(:) = natveg_tot(:) + veget_cov_max_nat(:,j)

       ! Sum up lai
       lai_nat(:) = lai_nat(:) + veget_cov_max_nat(:,j) * lai(:,j)

    ENDDO ! Loop over # PFTs

    DO j = 2, nvm ! Loop over # PFTs

       !WRITE(numout,*) 'yd:alloc0_lai nvm=',j
       IF ( natural(j) .AND. .NOT. pasture(j)) THEN

          ! Use the mean LAI over all natural PFTs when estimating light stress 
          ! on a specific natural PFT
          lai_around(:,j) = lai_nat(:)
       ELSE

          ! Use the actual LAI (specific for that PFT) when estimating light
          ! stress on a specific agricultural PFT
          lai_around(:,j) = lai(:,j)
       ENDIF

    ENDDO ! Loop over # PFTs

    !WRITE(numout,*) 'yd:alloc0_lai ,ffballocmax(:,62)',ffballocmax(62)
    !WRITE(numout,*) 'yd:alloc0_lai ,ffballocmax',ffballocmax(62:67)
    !WRITE(numout,*) 'yd:alloc0_lai nvm=62,lai_around(:,62)',lai_around(:,62)

    !! 1.2.4 Calculate LAI threshold below which carbohydrate reserve is used.
    !        Lai_max is a PFT-dependent parameter specified in stomate_constants.f90
!JCMODIF
    lai_happy(:) = lai_max(:) * lai_max_to_happy(:)
!    DO j = 2, nvm
!      lai_happy(:,j) = lai_max(j) * lai_max_to_happy(j)
!    ENDDO
!ENDJCMODIF
 !! 2. Use carbohydrate reserve to support growth and update leaf age

    ! Save old leaf mass, biomass got last updated in stomate_phenology.f90
    lm_old(:,:) = biomass(:,:,ileaf,icarbon)

    WRITE(numout,*) 'yd:alloc0_nvm=62,biomass(:,62,:,icarbon)',biomass(:,62,:,icarbon)
    DO j = 2, nvm ! Loop over # PFTs 
       IF (is_oilpalm(j)) THEN
          WRITE(numout,*) 'yd:alloc1 nvm=',j,'f_alloc(i,j,isapabove)', f_alloc(:,j,isapabove)
       ENDIF
       !! 2.1 Calculate demand for carbohydrate reserve to support leaf and root growth.
       !      Maximum time (days) since start of the growing season during which carbohydrate 
       !      may be used
       IF ( is_tree(j) ) THEN
          reserve_time = reserve_time_tree
       ELSE
          reserve_time = reserve_time_grass
       ENDIF

       IF ( ok_LAIdev(j) ) THEN
           ! For crops, at the moment, we do not relocate C from leaf to root
           ! & carbres, but future processes should be added during maturity
           ! stage, transloc_leaf represents the overall allocated biomass for
           ! leaf
           ! This is likely a bug that lead to decrease of leaf age for
           ! crops. Xuhui
           DO i = 1,npts

              IF ( ( biomass(i,j,ileaf,icarbon) .GT. zero ) .AND. &
                    ( .NOT. senescence(i,j) ) .AND. &
                    ( lai(i,j) .LT. lai_happy(j) ) .AND. &
                    ( when_growthinit(i,j) .LT. reserve_time ) ) THEN

                   !transloc_leaf(i) = deltai(i,j)/ssla(i,j)*10000.
                   transloc_leaf(i) = 0.0 !deltai(i,j)/ssla(i,j)*10000.
    !              biomass(i,j,ileaf) = biomass(i,j,ileaf)
    !              biomass(i,j,iroot) = biomass(i,j,iroot)
    !              biomass(i,j,icarbres) = biomass(i,j,icarbres)
              ENDIF
           ENDDO
       ELSE

           ! Growth is only supported by the use of carbohydrate reserves if the following 
           ! conditions are  statisfied:\n
           ! - PFT is not senescent;\n
           ! - LAI must be low (i.e. below ::lai_happy) and\n
           ! - Day of year of the simulation is in the beginning of the growing season.
           IF ((j .EQ. 2) .OR. (j .EQ. 62)) THEN
              WRITE(numout,*) 'xuhui: stomate_alloc debug for PFT '
              WRITE(numout,*) 'when_growthinit(:,j),', when_growthinit(:,j)
              WRITE(numout,*) 'lai(:,j),', lai(:,j)
              WRITE(numout,*) 'lai_happy(j),', lai_happy(j)
           ENDIF
           WHERE ( ( biomass(:,j,ileaf,icarbon) .GT. zero ) .AND. & 
                ( .NOT. senescence(:,j) ) .AND. &
                ( lai(:,j) .LT. lai_happy(j) ) .AND. &
                ( when_growthinit(:,j) .LT. reserve_time ) ) 
    
              ! Determine the mass from the carbohydrate reserve that can be used @tex $(gC m^{-2})$ @endtex. 
              ! Satisfy the demand or use everything that is available 
              ! (i.e. ::biomass(:,j,icarbres)). Distribute the demand evenly over the time 
              ! required (::tau_leafinit) to develop a minimal canopy from reserves (::lai_happy). 
              use_reserve(:) = &
                   MIN( biomass(:,j,icarbres,icarbon), &
    !JCMODIF
    !               deux * dt/tau_leafinit(j) * lai_happy(j)/ sla(j) )
                   deux * dt/tau_leafinit(j) * lai_happy(j)/ sla_calc(:,j) ) !!+ delta_turnover(:,j) @chen
    !ENDJCMODIF
              ! Distribute the reserve over leaves and fine roots.
    	  ! The part of the reserve going to the leaves is the ratio of default leaf allocation to default root and leaf allocation.
    	  ! The remaining of the reserve is alocated to the roots.
              transloc_leaf(:) = L0(j)/(L0(j)+R0(j)) * use_reserve(:)
              biomass(:,j,ileaf,icarbon) = biomass(:,j,ileaf,icarbon) + transloc_leaf(:)
              biomass(:,j,iroot,icarbon) = biomass(:,j,iroot,icarbon) + ( use_reserve(:) - transloc_leaf(:) )
    
              ! Adjust the carbohydrate reserve mass by accounting for the reserves allocated to leaves and roots during
              ! this time step
              biomass(:,j,icarbres,icarbon) = biomass(:,j,icarbres,icarbon) - use_reserve(:)
    
           ELSEWHERE
    
              transloc_leaf(:) = zero
    
           ENDWHERE
    !gmjc modify leaf_age and leaf_frac when_growthinit
           IF (is_grassland_manag(j)) THEN
              WHERE ( when_growthinit(:,j) .EQ. 0)
                 leaf_age(:,j,2)=10
                 leaf_age(:,j,3)=20
                 leaf_age(:,j,4)=30
              ENDWHERE
              WHERE ( ( biomass(:,j,ileaf,icarbon) .GT. 0.0 ) .AND. &
                    ( when_growthinit_cut(:,j) .LT. reserve_time_cut )   .AND. &
    !                 ( lai(:,j) .LT. lai_max_calc(:,j)*0.5 ) )!1.25
                    ( lai(:,j) .LT. lai_happy(j) ) )
    !             ( turnover_daily(:,j,ileaf) .LT. turnover_max_cut) )
    ! determine mass to put on
                    use_reserve(:) = &
                        MIN( biomass(:,j,icarbres,icarbon), &
                        !                   2._r_std * dt/tau_leafinit_cut *
                        !                   lai_happy_cut  / sla(j) )
                        2._r_std * dt/tau_leafinit_cut * lai_happy_cut /sla_calc(:,j) )
    ! grow leaves and fine roots
    
                    transloc_leaf(:) = use_reserve(:)

                    biomass(:,j,ileaf,icarbon) = biomass(:,j,ileaf,icarbon) + transloc_leaf(:)
                    biomass(:,j,iroot,icarbon) = biomass(:,j,iroot,icarbon) + ( use_reserve(:) - transloc_leaf(:) )
    
    ! decrease reserve mass
    
                    biomass(:,j,icarbres,icarbon) = biomass(:,j,icarbres,icarbon) - use_reserve(:)
    
              ENDWHERE
           ENDIF
    
    !end gmjc    
       ENDIF
       !! 2.2 Update leaf age
       !! 2.2.1 Decrease leaf age in youngest class
       !        Adjust the mass of the youngest leaves by the newly grown leaves
       leaf_mass_young(:) = leaf_frac(:,j,1) * lm_old(:,j) + transloc_leaf(:)

       WHERE ( ( transloc_leaf(:) .GT. min_stomate ) .AND. ( leaf_mass_young(:) .GT. min_stomate ) )
          
          ! Adjust leaf age by the ratio of leaf_mass_young (t-1)/leaf_mass_young (t) 
          leaf_age(:,j,1) = MAX( zero, leaf_age(:,j,1) * ( leaf_mass_young(:) - transloc_leaf(:) ) / &
               leaf_mass_young(:) )

       ENDWHERE

       !! 2.2.2 Update leaf mass fraction for the different age classes
       !        Mass fraction in the youngest age class is calculated as the ratio between
       !        the new mass in the youngest class and the total leaf biomass 
       !        (inc. the new leaves)
       WHERE ( biomass(:,j,ileaf,icarbon) .GT. min_stomate )
          
          leaf_frac(:,j,1) = leaf_mass_young(:) / biomass(:,j,ileaf,icarbon)

       ENDWHERE

       ! Mass fraction in the other classes is calculated as the ratio bewteen
       ! the current mass in that age and the total leaf biomass 
       ! (inc. the new leaves)\n
       DO m = 2, nleafages ! Loop over # leaf age classes

          WHERE ( biomass(:,j,ileaf,icarbon) .GT. min_stomate )

             leaf_frac(:,j,m) = leaf_frac(:,j,m) * lm_old(:,j) / biomass(:,j,ileaf,icarbon)

          ENDWHERE

       ENDDO ! Loop over # leaf age classes
    ENDDO ! loop over # PFTs
                 
    DO j = 2,nvm                                                              
       sla_age2(:,j) = sla_max(j)*0.9 !!0.9                                           
       sla_age3(:,j) = sla_max(j)*0.85 !!0.85                                          
       sla_age4(:,j) = sla_max(j)*0.8 !!0.8                                          
    ENDDO                                                                       
    lai_cl2(:,:) = biomass(:,:,ileaf,icarbon)*leaf_frac(:,:,2)*sla_age2(:,:)    
    lai_cl3(:,:) = biomass(:,:,ileaf,icarbon)*leaf_frac(:,:,3)*sla_age3(:,:)    
    lai_cl4(:,:) = biomass(:,:,ileaf,icarbon)*leaf_frac(:,:,4)*sla_age4(:,:)    
    lai_cl1(:,:) = lai(:,:) - lai_cl2(:,:) -lai_cl3(:,:) -lai_cl4(:,:) 
 !! 3. Calculate allocatable fractions of biomass production (NPP) 
     
    ! Calculate fractions of biomass production (NPP) to be allocated to the different 
    ! biomass components.\n
    ! The fractions of NPP allocated (0-1, unitless) to the different compartments depend on the
    ! availability of light, water, and nitrogen.
    WRITE(numout,*) 'swdown_week in alloc:', swdown_week

    DO j = 2, nvm !! yidi modify  !xchen@ for tropical evergreen forest

      !WRITE(numout,*) 'yd: alloc2 testnvm1 , nvm=',j 
      !IF (pheno_model(j) == "none" .AND. is_tropical(j)) THEN  ! is tropical evergreen forest? nvm=2-7,62-67
      IF (is_oilpalm(j)) THEN  ! is tropical evergreen forest? nvm=2-7,62-67
        WRITE(numout,*) 'yd: alloc2 testnvm2 , nvm=',j 
        IF (.NOT. ok_LAIdev(j) ) THEN ! bypass this part for crops
           !WRITE(numout,*) 'yd: alloc2 testnvm3 , nvm=',j 
           ! Reset values
           RtoLSR(:) = zero
           LtoLSR(:) = zero
           StoLSR(:) = zero
    !gmjc
           alloc_leaf(:)=1.0
    !end gmjc
           ! For trees, partitioning between above and belowground sapwood biomass is a function 
           ! of age. An older tree gets more allocation to the aboveground sapwoood than a younger tree.
           ! For the other PFTs it is prescribed. 
           ! ::alloc_min, ::alloc_max and ::demi_alloc are specified in stomate_constants.f90
           IF ( is_tree(j) ) THEN
    
              alloc_sap_above(:) = alloc_min(j)+(alloc_max(j)-alloc_min(j))*(un-EXP(-age(:,j)/demi_alloc(j)))
           
           ELSE
              
              alloc_sap_above(:) = alloc_sap_above_grass
           ENDIF
    
    
           !! 3.1 Calculate light stress, water stress and proxy for nitrogen stress.\n
           !      For the limiting factors a low value indicates a strong limitation
           WRITE(numout,*) 'yd: alloc2.1 alloc_sap_above',alloc_sap_above,'nvm=',j
           WRITE(numout,*) 'yd: alloc2.1 biomass(:,j,ileaf,icarbon)', biomass(:,j,ileaf,icarbon),'nvm=',j
           WRITE(numout,*) 'yd: alloc2.1 lai_around', lai_around(:,j),'nvm=',j
           WRITE(numout,*) 'yd: alloc2.1 max_possible_lai', max_possible_lai,'nvm=',j
           WRITE(numout,*) 'yd: alloc2.1 leaf_frac', leaf_frac(:,j,1:4),'nvm=',j
           WRITE(numout,*) 'yd: alloc2.1 lai_cl4', lai_cl4(:,j),'nvm=',j
           WRITE(numout,*) 'yd: alloc2.1 lai_cl3', lai_cl3(:,j),'nvm=',j
           WRITE(numout,*) 'yd: alloc2.1 lai_cl2', lai_cl2(:,j),'nvm=',j
           WRITE(numout,*) 'yd: alloc2.1 lai_cl1', lai_cl1(:,j),'nvm=',j
           WRITE(numout,*) 'yd: alloc2.1 lai', lai(:,j),'nvm=',j
           WRITE(numout,*) 'yd: alloc2.1 lai_max', lai_max(j),'nvm=',j
           WRITE(numout,*) 'yd: alloc2.1 lsr_coeff', lsr_coeff
           WRITE(numout,*) 'yd: alloc2.1 last_lai', last_lai(:,j)
           WRITE(numout,*) 'yd: alloc2.1 sup_flag', sup_flag(:,j),'nvm=',j
           WHERE ( biomass(:,j,ileaf,icarbon) .GT. min_stomate )
    
              !! 3.1.1 Light stress
              !        Light stress is a function of the mean lai on the natural part of the grid box 
              !        and of the PFT-specific LAI for agricultural crops. In line with the DGVM, natural
              !        PFTs in the same gridbox are treated as if they were spatially mixed whereas 
              !        agricultural PFTs are considered to be spatially separated.
              !        The calculation of the lights stress depends on the extinction coefficient (set to 0.5)
              !        and of a mean LAI.

              WHERE( lai_around(:,j) < max_possible_lai )
    
                 limit_L(:) = MAX( 0.1_r_std, EXP( -ext_coeff(j) * lai_around(:,j) ) )
              
              ELSEWHERE
                 
                 limit_L(:) = 0.1_r_std
              
              ENDWHERE
    

              !! 3.1.2 Water stress
              !        Water stress is calculated as the weekly moisture availability.
              !        Weekly moisture availability is calculated in stomate_season.f90.
              limit_W(:) = MAX( 0.1_r_std, MIN( un, moiavail_week(:,j) ) )
    
    
              !! 3.1.3 Proxy for nitrogen stress
              !         The proxy for nitrogen stress depends on monthly soil water availability
              !         (::soilhum_month) and monthly soil temperature (::tsoil_month). See section
              !         1.2.2 for details on how ::t_nitrogen and ::h_nitrogen were calculated.\n
              !         Currently nitrogen-stress is calculated for both natural and agricultural PFTs.
              !         Due to intense fertilization of agricultural PFTs this is a strong 
              !         assumption for several agricultural regions in the world (US, Europe, India, ...)
              !         Water stress on nitrogen mineralisation 
              limit_N_hum(:) = MAX( undemi, MIN( un, h_nitrogen(:) ) )
    
    	  ! Temperature stress on nitrogen mineralisation using a Q10 decomposition model
              ! where Q10 was set to 2
              limit_N_temp(:) = 2.**((t_nitrogen(:) - ZeroCelsius - Nlim_tref )/Nlim_Q10)
              limit_N_temp(:) = MAX( 0.1_r_std, MIN( un, limit_N_temp(:) ) )
    
              ! Combine water and temperature factors to get total nitrogen stress
              limit_N(:) = MAX( 0.1_r_std, MIN( un, limit_N_hum(:) * limit_N_temp(:) ) )
    
              ! Take the most limiting factor among soil water and nitrogen 
              limit_WorN(:) = MIN( limit_W(:), limit_N(:) )
    
              ! Take the most limiting factor among aboveground (i.e. light) and belowground 
              ! (i.e. water & nitrogen) limitations
              limit(:) = MIN( limit_WorN(:), limit_L(:) )
    
              !! 3.2 Calculate ratio between allocation to leaves, sapwood and roots
              !      Partitioning between belowground and aboveground biomass components is assumed 
              !      to be proportional to the ratio of belowground and aboveground stresses.\n
              !      \latexonly 
              !        \input{alloc1.tex} 
              !      \endlatexonly
              !      Root allocation is the default root allocation corrected by a normalized ratio of aboveground stress to total stress.
    	  !      The minimum root allocation is 0.15.
              !RtoLSR(:) = &
              !     MAX( .15_r_std, &
              !     R0(j) * trois * limit_L(:) / ( limit_L(:) + deux * limit_WorN(:) ) )
    	  !      The minimum root allocation is modified to 0.1 for oil palm.
              RtoLSR(:) = &
                   MAX( min_R, &
                   R0(j) * trois * limit_L(:) / ( limit_L(:) + deux * limit_WorN(:) ) )
               
              ! Sapwood allocation is the default sapwood allocation corrected by a normalized ratio of belowground stress to total stress.
              StoLSR(:) = S0(j) * 3. * limit_WorN(:) / ( 2._r_std * limit_L(:) + limit_WorN(:) )
    
              ! Leaf allocation is calculated as the remaining allocation fraction
    	  ! The range of variation of leaf allocation is constrained by ::min_LtoLSR and ::max_LtoLSR. 
    !JCMODIF
              LtoLSR(:) = 1.0 - RtoLSR(:) - StoLSR(:)
!!xchen201701
              old_LtoLSR(:) = LtoLSR(:)
              old_RtoLSR(:) = RtoLSR(:)
              old_StoLSR(:) = StoLSR(:)

              !LtoLSR(:) = MIN((ABS(swdown_week(:))*EXP(-0.5*(lai_cl4(:,j)))/30.0)**3.0*LtoLSR(:), 0.90)              
              !old_LtoLSR2(:) = LtoLSR(:) !! yidi
              !LtoLSR(:) = MAX(min_LtoLSR, MIN(max_LtoLSR, LtoLSR(:)))
!!VPD and water stress
!!xchen          

!! yidi      
              !LtoLSR(:) = 1.0 - RtoLSR(:) - StoLSR(:)
              xzLtoLSR(:) = MIN((ABS(swdown_week(:))*EXP(-0.5*(lai_cl4(:,j)))/lsr_coeff(1))**lsr_coeff(2)*LtoLSR(:), max_LtoLSR)

              !LtoLSR(:) = min_LtoLSR+(max_LtoLSR-min_LtoLSR)*(ABS(swdown_week(:))*(1-EXP(-lsr_coeff(3)*(lai_cl4(:,j))))/lsr_coeff(1))**lsr_coeff(2)*LtoLSR(:)
              LtoLSR(:) = min_LtoLSR+(max_LtoLSR-min_LtoLSR)*(ABS(swdown_week(:))*EXP(-lsr_coeff(3)*(lai_cl4(:,j)))/lsr_coeff(1))**lsr_coeff(2)*LtoLSR(:)
              new_LtoLSR(:) = LtoLSR(:)
              LtoLSR(:) = MAX(min_LtoLSR, MIN(max_LtoLSR, LtoLSR(:)))
              !lai_diff(:,j)=lai_max(j)-lai_round(:,j)
              WHERE ((sup_flag(:,j) .GT. val_exp) )
                 sup_flag(:,j)=zero
              ENDWHERE

              WHERE ((sup_flag(:,j) .GT. min_stomate) .AND. (lai_around(:,j) .GE. lai_max(j)))
                 sup_flag(:,j)=zero
              ENDWHERE

              WHERE (last_lai(:,j) .LT. lai_max(j) .AND. lai_around(:,j) .LT. lai_max(j) .AND. last_lai(:,j) .GT. lai_around(:,j)) 
                 !LtoLSR(:) = sup_LtoLSR
                 sup_flag(:,j)=1.0 
              ENDWHERE
              WHERE (sup_flag(:,j) .GT. min_stomate)
                 LtoLSR(:) = sup_LtoLSR
              ENDWHERE
              last_lai(:,j) = lai(:,j)
!! yidi
!! xzchen
              !StoLSR(:) = MAX(0.05, StoLSR(:) - 0.5*(LtoLSR(:)-old_LtoLSR(:)))
              !StoLSR(:) = 0.41*StoLSR(:)
              !!!StoLSR(:) = 1.0 - LtoLSR(:) - RtoLSR(:)
              !RtoLSR = 1.0 -  LtoLSR(:) - StoLSR(:)
    	 
!! xzchen
!! yidi/default
              
              RtoLSR(:) = MAX(min_RtoLSR, RtoLSR(:) - rs_coeff(1)*ABS(LtoLSR(:)-old_LtoLSR(:)))
              RtoLSR(:) = MIN(max_RtoLSR, RtoLSR(:)) 
              StoLSR(:) = 1.0 - LtoLSR(:) - RtoLSR(:)
              !RtoLSR = 1.0 -  LtoLSR(:) - StoLSR(:)
              new_RtoLSR(:) = RtoLSR(:)
              new_StoLSR(:) = StoLSR(:)

           ! Check whether allocation needs to be adjusted. If LAI exceeds maximum LAI 
           ! (::lai_max), no addition carbon should be allocated to leaf biomass. Allocation is 
           ! then partioned between root and sapwood biomass.
             WHERE ( (biomass(:,j,ileaf,icarbon) .GT. min_stomate) .AND. (lai(:,j) .GT. lai_max(j)) )
    
!! xzchen
              !StoLSR(:) = StoLSR(:) + 0.41*LtoLSR(:)
              !LtoLSR(:) = zero
              !RtoLSR = 1.0 -  LtoLSR(:) - StoLSR(:) 
!! xzchen

              StoLSR(:) = StoLSR(:) + rs_coeff(2)*LtoLSR(:)
              LtoLSR(:) = zero
              RtoLSR(:) = MAX(min_RtoLSR, RtoLSR(:) + (1-rs_coeff(2))*LtoLSR(:) )
              RtoLSR(:) = MIN(max_RtoLSR, RtoLSR(:)) 
              StoLSR(:) = 1.0 - LtoLSR(:) - RtoLSR(:)
              !RtoLSR = 1.0 -  LtoLSR(:) - StoLSR(:) 
             
              new_RtoLSR2(:) = RtoLSR(:)
              new_LtoLSR2(:) = LtoLSR(:)
              new_StoLSR2(:) = StoLSR(:)
!! yidi
             ENDWHERE
           ENDWHERE  

              !WRITE(numout,*) 'yd: alloc2.2 limit_WorN', limit_WorN(:)
              !WRITE(numout,*) 'yd: alloc2.2 limit_L', limit_L(:)
              WRITE(numout,*) 'yd: alloc2.2 xzLtoLSR', xzLtoLSR(:)
              WRITE(numout,*) 'yd: alloc2.2 nvm=',j,' old_LtoLSR', old_LtoLSR(:)
              WRITE(numout,*) 'yd: alloc2.2 old_RtoLSR', old_RtoLSR(:),' nvm=',j                                            
              WRITE(numout,*) 'yd: alloc2.2 old_StoLSR', old_StoLSR(:)
              WRITE(numout,*) 'yd: alloc2.2 new_LtoLSR', new_LtoLSR(:)
              WRITE(numout,*) 'yd: alloc2.3 new_RtoLSR', new_RtoLSR(:)
              WRITE(numout,*) 'yd: alloc2.3 new_StoLSR', new_StoLSR(:)
              WRITE(numout,*) 'yd: alloc2.2 new_LtoLSR2', new_LtoLSR2(:)
              WRITE(numout,*) 'yd: alloc2.3 new_RtoLSR2', new_RtoLSR2(:)
              WRITE(numout,*) 'yd: alloc2.3 new_StoLSR2', new_StoLSR2(:)
              WRITE(numout,*) 'yd: alloc2.3 swdown_week', swdown_week(:)

           !IF (is_oilpalm(j) .OR. ((j .GE. 2) .AND. (j .LE. 7))) THEN
           IF (is_oilpalm(j)) THEN
              WRITE(numout,*) 'yd: alloc3 LtoLSR', LtoLSR(:)
              WRITE(numout,*) 'yd: alloc3 RtoLSR', RtoLSR(:)                                            
              WRITE(numout,*) 'yd: alloc3 StoLSR', StoLSR(:)
              WRITE(numout,*) 'yd: alloc2.1 last_lai', last_lai(:,j)
              WRITE(numout,*) 'yd: alloc2.3 sup_flag', sup_flag(:,j),'nvm=',j
              WRITE(numout,*) 'yd: alloc3 nvm=',j, ' biomass(i,j,:,icarbon)=', biomass(1,j,:,icarbon) 
           ENDIF                                            
           !CALL xios_orchidee_send_field("LtoLSR",LtoLSR(:))                         
           !CALL xios_orchidee_send_field("RtoLSR",RtoLSR(:))                           
           !CALL xios_orchidee_send_field("StoLSR",StoLSR(:)) 

           !! 3.3 Calculate the allocation fractions.
           !      The allocation fractions (::f_alloc) are an output variable (0-1, unitless). f_alloc
           !      has three dimensions (npts,nvm,nparts). Where ::npts is the number of grid cells, ::nvm is the 
           !      number of PFTs and ::nparts the number of biomass components. Currently six biomass compartments
           !      are distinguished: (1) Carbon reserves, (2) Aboveground sapwood, (3) Belowground 
           !      sapwood, (4) Roots, (5) fruits/seeds and (6) Leaves.@tex $(gC m^{-2})$ @endtex \n
 
           WRITE(numout,*) 'yd: alloc4 testnvm4 , nvm=',j 
           DO i = 1, npts ! Loop over grid cells 
              WRITE(numout,*) 'yd: alloc4 1) testnpts , npts=',i 
              WRITE(numout,*) 'yd: alloc4 1), nvm=',j, ' biomass(i,j,ileaf,icarbon)=', biomass(1,j,ileaf,icarbon) 
              IF ( biomass(i,j,ileaf,icarbon) .GT. min_stomate ) THEN
          
                 WRITE(numout,*) 'yd: alloc4 2) biomass/senescence, nvm=', j 
                 WRITE(numout,*) 'yd: alloc4 1) senescence(i,j)=', senescence(i,j) 
                 WRITE(numout,*) 'yd:alloc nvm=',j,'f_alloc(i,j,isapabove)', f_alloc(i,j,isapabove)
                 IF ( senescence(i,j) ) THEN
                    
                    WRITE(numout,*) 'yd: alloc4 3.3,icarbres 2) senescence(i,j)=', senescence(i,j) 
                    !! 3.3.1 Allocate all C to carbohydrate reserve
                    !        If the PFT is senescent allocate all C to carbohydrate reserve, 
                    !        then the allocation fraction to reserves is 1. 
                    f_alloc(i,j,icarbres) = un
    
                 ELSE
    
                    WRITE(numout,*) 'yd: alloc4 3.3,icarbres nvm=',j,'3) senescence(i,j)=', senescence(i,j) 
                    !! 3.3.2 Allocation during the growing season  
                    f_alloc(i,j,ifruit) = f_fruit
                 
                    ! Allocation to the carbohydrate reserve is proportional to leaf and root 
                    ! allocation. If carbon is allocated to the carbohydrate reserve, rescaling 
                    ! of allocation factors is required to ensure carbon mass preservation.\n
                    ! Carbon is allocated to the carbohydrate reserve when the pool size of the 
                    ! reserve is less than the carbon needed to grow a canopy twice the size of 
                    ! the maximum LAI (::lai_max). Twice the size was used as a threshold because
                    ! the reserves needs to be sufficiently to grow a canopy and roots. In case 
                    ! the carbohydrate pool is full, there is no need to rescale the other
                    ! allocation factors.
    	        ! If there is no rescaling of the allocation factors (carbres=1, no carbon put
    		! to reserve), then fraction remaining after fruit allocation (1-fruit_alloc) 
    		! is distributed between leaf, root and sap (sap carbon also distributed between   
    		! sap_above and sap_below with factor alloc_sap_above).
    		! If carbon is allocated to the carbohydrate reserve, all these factors are 
    		! rescaled through carb_rescale, and an allocation fraction for carbohydrate pool 
    		! appears. carb_rescale depends on the parameter (::ecureuil). 
                    ! (::ecureuil) is the fraction of primary leaf and root allocation put into 
                    ! reserve, it is specified in stomate_constants.f90 and is either 0 or 1.
    !JCMODIF
                    IF ( ( biomass(i,j,icarbres,icarbon)*sla_calc(i,j) ) .LT. 2*lai_max(j) ) THEN
    !ENDJCMODIF       
                       carb_rescale(i) = un / ( un + ecureuil(j) * ( LtoLSR(i) + RtoLSR(i) ) )
                    ELSE
                       carb_rescale(i) = un
                    ENDIF

                    !WRITE(numout,*) 'yd:alloc5 nvm=',j,'f_alloc(i,j,isapabove)', f_alloc(i,j,isapabove)
                    f_alloc(i,j,ileaf) = LtoLSR(i) * ( un - f_alloc(i,j,ifruit) ) * carb_rescale(i)
                    f_alloc(i,j,isapabove) = StoLSR(i) * alloc_sap_above(i) * &
                         ( un - f_alloc(i,j,ifruit) ) * carb_rescale(i)
                    f_alloc(i,j,isapbelow) = StoLSR(i) * ( un - alloc_sap_above(i) ) * &
                         ( un - f_alloc(i,j,ifruit) ) * carb_rescale(i)
                    f_alloc(i,j,iroot) = RtoLSR(i) * (un - f_alloc(i,j,ifruit) ) * carb_rescale(i)
                    f_alloc(i,j,icarbres) = ( un - carb_rescale(i) ) * ( un - f_alloc(i,j,ifruit) )
                    WRITE(numout,*) 'yd:alloc5 nvm=',j,'f_alloc(i,j,:)', f_alloc(i,j,:)
                    WRITE(numout,*) 'yd:alloc5 testnpts , npts=',i 
              !! yidi for OP
                    IF (ok_oilpalm) THEN
                        FtoPHY(i,j,:) = zero
                        alloc_phy(i,j,:)= zero
                        phytomer_age_crit(:,j) = phytomeragecrit(j)
                        ffbharvest_age_crit(:,j) = ffbharvestagecrit(j)
                        IF (is_oilpalm_ffbharvest(j)) THEN
                           WRITE(numout,*) 'yd: alloc6; nvm=',j !! yd: output
                           WRITE(numout,*) 'yd: alloc6; ffb_coeff',ffb_coeff(:) !! yd: output
                          ! WRITE(numout,*) 'yd: alloc6; ffbharvestagecrit',ffbharvestagecrit(j) !! yd: output
                           WRITE(numout,*) 'yd: alloc6; ftophymin',ftophy_min(j) !! yd: output
                           WRITE(numout,*) 'yd: alloc6; ftophymax',ftophy_max(j) !! yd: output
                           DO p=1,nphs
                              IF (phytomer_age (i,j,p) .LE. ffbharvest_age_crit(i,j))  THEN
                                 !FtoPHY(i,j,p)=MIN(un - exp( - ffb_coeff(1) * phytomer_age(i,j,p)) , ffb_coeff(2))
                                 !FtoPHY(i,j,p)=MAX(FtoPHY(i,j,p),0.0)
                                 FtoPHY(i,j,p)=ftophy_min(j)+(ftophy_max(j)-ftophy_min(j))*(un - exp( - ffb_coeff(1) * phytomer_age(i,j,p)))
                                 !FtoPHY(i,j,p)=FtoPHY(i,j,p)*(un - exp( - age(i,j) / ffb_coeff(2)))
                              ELSE
                                 FtoPHY(i,j,p) = 0.0
                              ENDIF
                           ENDDO
                              WRITE(numout,*) 'yd: alloc6; nvm=',j !! yd: output
                              WRITE(numout,*) 'yd: alloc6; FtoPHY(i,j,:)',FtoPHY(i,j,:) !! yd: output
                              WRITE(numout,*) 'yd: alloc6; phytomer_age(:,j,:)',phytomer_age(:,j,:) !! yd: output
                        ENDIF

                        IF (is_oilpalm(j)) THEN
                           WRITE(numout,*) 'yd: alloc6; nvm=',j !! yd: output
                           WRITE(numout,*) 'yd: alloc6; phytomeragecrit',phytomeragecrit(j) !! yd: output
                           WRITE(numout,*) 'yd: alloc6; phy_coeff',phy_coeff !! yd: output
                           DO p=1,nphs
                              IF (phytomer_age (i,j,p) .LE. phytomer_age_crit(i,j) .AND. phytomer_age (i,j,p) .GT. min_stomate) THEN 
                                 !alloc_phy(i,j,p)= (phy_coeff(1) * phytomer_age(i,j,p) / phytomeragecrit(j))**phy_coeff(2)
                                 alloc_phy(i,j,p)= phyalloc_min(j)+(phyalloc_max(j)-phyalloc_min(j))*((phy_coeff(1) * phytomer_age(i,j,p) / phytomeragecrit(j))**phy_coeff(2))
                                 alloc_phy(i,j,p)= alloc_phy(i,j,p)*(un - exp( - age(i,j) / phy_coeff(3)))
                                 alloc_phy(i,j,p)=MAX(MIN(alloc_phy(i,j,p), phy_allocmax),0.0)
                              ELSE
                                 alloc_phy(i,j,p) = 0.0
                              ENDIF
                           ENDDO
                              WRITE(numout,*) 'yd: alloc6; nvm=',j !! yd: output
                              WRITE(numout,*) 'yd: alloc6; alloc_phy(i,j,:)',alloc_phy(i,j,:) !! yd: output
                              WRITE(numout,*) 'yd: alloc6 phytomer_age(:,j,:)',phytomer_age(:,j,:) !! yd: output
                        ENDIF
                    ENDIF

              !! yidi for OP 
                 ENDIF  ! Is senescent?
    
              ENDIF  ! There are leaves
    
           ENDDO  ! Loop over # pixels - domain size
        ENDIF ! Is crop?

      ELSE !! for the other pft, remain the default (yidi)              

        IF (.NOT. ok_LAIdev(j) ) THEN ! bypass this part for crops

           ! Reset values
           RtoLSR(:) = zero
           LtoLSR(:) = zero
           StoLSR(:) = zero
    !gmjc
           alloc_leaf(:)=1.0
    !end gmjc
           ! For trees, partitioning between above and belowground sapwood biomass is a function 
           ! of age. An older tree gets more allocation to the aboveground sapwoood than a younger tree.
           ! For the other PFTs it is prescribed. 
           ! ::alloc_min, ::alloc_max and ::demi_alloc are specified in stomate_constants.f90
           IF ( is_tree(j) ) THEN
    
              alloc_sap_above(:) = alloc_min(j)+(alloc_max(j)-alloc_min(j))*(un-EXP(-age(:,j)/demi_alloc(j)))
           
           ELSE
              
              alloc_sap_above(:) = alloc_sap_above_grass
    !!end gmjc       
           ENDIF
    
    
           !! 3.1 Calculate light stress, water stress and proxy for nitrogen stress.\n
           !      For the limiting factors a low value indicates a strong limitation
           WHERE ( biomass(:,j,ileaf,icarbon) .GT. min_stomate )
    
              !! 3.1.1 Light stress
              !        Light stress is a function of the mean lai on the natural part of the grid box 
              !        and of the PFT-specific LAI for agricultural crops. In line with the DGVM, natural
              !        PFTs in the same gridbox are treated as if they were spatially mixed whereas 
              !        agricultural PFTs are considered to be spatially separated.
              !        The calculation of the lights stress depends on the extinction coefficient (set to 0.5)
              !        and of a mean LAI.
              WHERE( lai_around(:,j) < max_possible_lai )
    
                 limit_L(:) = MAX( 0.1_r_std, EXP( -ext_coeff(j) * lai_around(:,j) ) )
              
              ELSEWHERE
                 
                 limit_L(:) = 0.1_r_std
              
              ENDWHERE
    
              !! 3.1.2 Water stress
              !        Water stress is calculated as the weekly moisture availability.
              !        Weekly moisture availability is calculated in stomate_season.f90.
              limit_W(:) = MAX( 0.1_r_std, MIN( un, moiavail_week(:,j) ) )
    
    
              !! 3.1.3 Proxy for nitrogen stress
              !         The proxy for nitrogen stress depends on monthly soil water availability
              !         (::soilhum_month) and monthly soil temperature (::tsoil_month). See section
              !         1.2.2 for details on how ::t_nitrogen and ::h_nitrogen were calculated.\n
              !         Currently nitrogen-stress is calculated for both natural and agricultural PFTs.
              !         Due to intense fertilization of agricultural PFTs this is a strong 
              !         assumption for several agricultural regions in the world (US, Europe, India, ...)
              !         Water stress on nitrogen mineralisation 
              limit_N_hum(:) = MAX( undemi, MIN( un, h_nitrogen(:) ) )
    
    	  ! Temperature stress on nitrogen mineralisation using a Q10 decomposition model
              ! where Q10 was set to 2
              limit_N_temp(:) = 2.**((t_nitrogen(:) - ZeroCelsius - Nlim_tref )/Nlim_Q10)
              limit_N_temp(:) = MAX( 0.1_r_std, MIN( un, limit_N_temp(:) ) )
    
              ! Combine water and temperature factors to get total nitrogen stress
              limit_N(:) = MAX( 0.1_r_std, MIN( un, limit_N_hum(:) * limit_N_temp(:) ) )
    
              ! Take the most limiting factor among soil water and nitrogen 
              limit_WorN(:) = MIN( limit_W(:), limit_N(:) )
    
              ! Take the most limiting factor among aboveground (i.e. light) and belowground 
              ! (i.e. water & nitrogen) limitations
              limit(:) = MIN( limit_WorN(:), limit_L(:) )
    
              !! 3.2 Calculate ratio between allocation to leaves, sapwood and roots
              !      Partitioning between belowground and aboveground biomass components is assumed 
              !      to be proportional to the ratio of belowground and aboveground stresses.\n
              !      \latexonly 
              !        \input{alloc1.tex} 
              !      \endlatexonly
              !      Root allocation is the default root allocation corrected by a normalized ratio of aboveground stress to total stress.
    	  !      The minimum root allocation is 0.15.
              RtoLSR(:) = &
                   MAX( .15_r_std, &
                   R0(j) * trois * limit_L(:) / ( limit_L(:) + deux * limit_WorN(:) ) )
    
              ! Sapwood allocation is the default sapwood allocation corrected by a normalized ratio of belowground stress to total stress.
              StoLSR(:) = S0(j) * 3. * limit_WorN(:) / ( 2._r_std * limit_L(:) + limit_WorN(:) )
    
              ! Leaf allocation is calculated as the remaining allocation fraction
    	  ! The range of variation of leaf allocation is constrained by ::min_LtoLSR and ::max_LtoLSR. 
    !JCMODIF
              LtoLSR(:) = un - RtoLSR(:) - StoLSR(:)
              LtoLSR(:) = MAX( min_LtoLSR, MIN( max_LtoLSR, LtoLSR(:) ) )
    
              ! Roots allocation is recalculated as the residual carbon after leaf allocation has been calculated. 
              RtoLSR(:) = un - LtoLSR(:) - StoLSR(:)
    
           ENDWHERE
    	   
           ! Check whether allocation needs to be adjusted. If LAI exceeds maximum LAI 
           ! (::lai_max), no addition carbon should be allocated to leaf biomass. Allocation is 
           ! then partioned between root and sapwood biomass.
           WHERE ( (biomass(:,j,ileaf,icarbon) .GT. min_stomate) .AND. (lai(:,j) .GT. lai_max(j)) )
    
              StoLSR(:) = StoLSR(:) + LtoLSR(:)
              LtoLSR(:) = zero
    
           ENDWHERE
    
           !IF ((j .GE. 2) .AND. (j .LE. 7)) THEN 
           !   WRITE(numout,*) 'yd:alloc end nvm=',j,'f_alloc(1,j,:)', f_alloc(1,j,:)
           !   WRITE(numout,*) 'yd: alloc3 LtoLSR', LtoLSR(:)
           !   WRITE(numout,*) 'yd: alloc3 RtoLSR', RtoLSR(:)                                            
           !   WRITE(numout,*) 'yd: alloc3 StoLSR', StoLSR(:)
           !   WRITE(numout,*) 'yd: alloc3 nvm=',j, ' biomass(i,j,:,icarbon)=', biomass(1,j,:,icarbon) 
           !ENDIF
           !! 3.3 Calculate the allocation fractions.
           !      The allocation fractions (::f_alloc) are an output variable (0-1, unitless). f_alloc
           !      has three dimensions (npts,nvm,nparts). Where ::npts is the number of grid cells, ::nvm is the 
           !      number of PFTs and ::nparts the number of biomass components. Currently six biomass compartments
           !      are distinguished: (1) Carbon reserves, (2) Aboveground sapwood, (3) Belowground 
           !      sapwood, (4) Roots, (5) fruits/seeds and (6) Leaves.@tex $(gC m^{-2})$ @endtex \n
           DO i = 1, npts ! Loop over grid cells 
    
              IF ( biomass(i,j,ileaf,icarbon) .GT. min_stomate ) THEN
          
                 IF ( senescence(i,j) ) THEN
                    
                    !! 3.3.1 Allocate all C to carbohydrate reserve
                    !        If the PFT is senescent allocate all C to carbohydrate reserve, 
                    !        then the allocation fraction to reserves is 1. 
                    f_alloc(i,j,icarbres) = un
    
                 ELSE
    
                    !! 3.3.2 Allocation during the growing season  
                    f_alloc(i,j,ifruit) = f_fruit
    
    
                    ! Allocation to the carbohydrate reserve is proportional to leaf and root 
                    ! allocation. If carbon is allocated to the carbohydrate reserve, rescaling 
                    ! of allocation factors is required to ensure carbon mass preservation.\n
                    ! Carbon is allocated to the carbohydrate reserve when the pool size of the 
                    ! reserve is less than the carbon needed to grow a canopy twice the size of 
                    ! the maximum LAI (::lai_max). Twice the size was used as a threshold because
                    ! the reserves needs to be sufficiently to grow a canopy and roots. In case 
                    ! the carbohydrate pool is full, there is no need to rescale the other
                    ! allocation factors.
    	        ! If there is no rescaling of the allocation factors (carbres=1, no carbon put
    		! to reserve), then fraction remaining after fruit allocation (1-fruit_alloc) 
    		! is distributed between leaf, root and sap (sap carbon also distributed between   
    		! sap_above and sap_below with factor alloc_sap_above).
    		! If carbon is allocated to the carbohydrate reserve, all these factors are 
    		! rescaled through carb_rescale, and an allocation fraction for carbohydrate pool 
    		! appears. carb_rescale depends on the parameter (::ecureuil). 
                    ! (::ecureuil) is the fraction of primary leaf and root allocation put into 
                    ! reserve, it is specified in stomate_constants.f90 and is either 0 or 1.
    !JCMODIF
                    IF ( ( biomass(i,j,icarbres,icarbon)*sla_calc(i,j) ) .LT. 2*lai_max(j) ) THEN
    !ENDJCMODIF
                       carb_rescale(i) = un / ( un + ecureuil(j) * ( LtoLSR(i) + RtoLSR(i) ) )
                    ELSE
                       carb_rescale(i) = un
                    ENDIF
    
                    IF ((j .GE. 2) .AND. (j .LE. 7)) THEN 
                       WRITE(numout,*) 'yd:alloc4 before nvm=',j,'f_alloc(1,j,:)', f_alloc(1,j,:)
                    ENDIF
                    f_alloc(i,j,ileaf) = LtoLSR(i) * ( un - f_alloc(i,j,ifruit) ) * carb_rescale(i)
                    f_alloc(i,j,isapabove) = StoLSR(i) * alloc_sap_above(i) * &
                         ( un - f_alloc(i,j,ifruit) ) * carb_rescale(i)
                    f_alloc(i,j,isapbelow) = StoLSR(i) * ( un - alloc_sap_above(i) ) * &
                         ( un - f_alloc(i,j,ifruit) ) * carb_rescale(i)
                    f_alloc(i,j,iroot) = RtoLSR(i) * (un - f_alloc(i,j,ifruit) ) * carb_rescale(i)
                    f_alloc(i,j,icarbres) = ( un - carb_rescale(i) ) * ( un - f_alloc(i,j,ifruit) )
                    !IF ( (printlev>=5) .AND. (j .GE. 10) ) THEN
                    !    WRITE(numout,*) 'nvm=',j,'f_alloc(i,10,ileaf)', f_alloc(i,10,ileaf)
                    !ENDIF
    
                    IF ((j .GE. 2) .AND. (j .LE. 7)) THEN 
                       WRITE(numout,*) 'yd:alloc5 after nvm=',j,'f_alloc(1,j,:)', f_alloc(1,j,:)
                    ENDIF
                 ENDIF  ! Is senescent?
    
              ENDIF  ! There are leaves
    
           ENDDO  ! Loop over # pixels - domain size
        ENDIF ! Is crop?

      ENDIF ! is tropical evergreen forest?  nvm=2-7 and 62-67 yidi
    ENDDO  ! loop over # PFTs
    
    DO j= 1,nvm
      IF (is_oilpalm(j)) THEN
         WRITE(numout,*) 'yd:alloc end nvm=',j,'f_alloc(1,j,:)', f_alloc(1,j,:)
         WRITE(numout,*) 'yd:alloc end nvm=',j,'alloc_phy(1,j,:)', alloc_phy(1,j,:)
         WRITE(numout,*) 'yd:alloc end nvm=',j,'SUM(alloc_phy(1,j,:))', SUM(alloc_phy(1,j,:))
         WRITE(numout,*) 'yd:alloc end nvm=',j,'SUM(FtoPHY(1,j,:))', SUM(FtoPHY(1,j,:))
      ENDIF
    ENDDO
    DO j=2,7
      WRITE(numout,*) 'yd:alloc end nvm=',j,'f_alloc(1,j,:)', f_alloc(1,j,:)
    ENDDO
    IF (printlev>=3) WRITE(numout,*) 'Leaving alloc#############'

  END SUBROUTINE alloc

END MODULE stomate_alloc
