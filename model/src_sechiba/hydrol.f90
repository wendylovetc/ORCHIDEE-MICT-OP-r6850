! ===================================================================================================\n
! MODULE        : hydrol
!
! CONTACT       : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE       : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        This module computes the soil moisture processes on continental points. 
!!
!!\n DESCRIPTION : contains hydrol_main, hydrol_initialize, hydrol_finalise, hydrol_init, 
!!                 hydrol_var_init, hydrol_alma,
!!                 hydrol_snow, hydrol_vegupd, hydrol_canop, hydrol_flood, hydrol_soil.
!!                 The assumption in this module is that very high vertical resolution is
!!                 needed in order to properly resolve the vertical diffusion of water in
!!                 the soils. Furthermore we have taken into account the sub-grid variability
!!                 of soil properties and vegetation cover by allowing the co-existence of
!!                 different soil moisture columns in the same grid box.
!!                 This routine was originaly developed by Patricia deRosnay.
!!
!! RECENT CHANGE(S) : None
!!
!! REFERENCE(S) :
!! - de Rosnay, P., J. Polcher, M. Bruen, and K. Laval, Impact of a physically based soil
!! water flow and soil-plant interaction representation for modeling large-scale land surface
!! processes, J. Geophys. Res, 107 (10.1029), 2002. \n
!! - de Rosnay, P. and Polcher J. (1998) Modeling root water uptake in a complex land surface scheme coupled 
!! to a GCM. Hydrology and Earth System Sciences, 2(2-3):239-256. \n
!! - de Rosnay, P., M. Bruen, and J. Polcher, Sensitivity of surface fluxes to the number of layers in the soil 
!! model used in GCMs, Geophysical research letters, 27 (20), 3329 - 3332, 2000. \n
!! - d’Orgeval, T., J. Polcher, and P. De Rosnay, Sensitivity of the West African hydrological
!! cycle in ORCHIDEE to infiltration processes, Hydrol. Earth Syst. Sci. Discuss, 5, 2251 - 2292, 2008. \n
!! - Carsel, R., and R. Parrish, Developing joint probability distributions of soil water retention
!! characteristics, Water Resources Research, 24 (5), 755 - 769, 1988. \n
!! - Mualem, Y., A new model for predicting the hydraulic conductivity of unsaturated porous
!! media, Water Resources Research, 12 (3), 513 - 522, 1976. \n
!! - Van Genuchten, M., A closed-form equation for predicting the hydraulic conductivity of
!! unsaturated soils, Soil Science Society of America Journal, 44 (5), 892 - 898, 1980. \n
!! - Campoy, A., Ducharne, A., Cheruy, F., Hourdin, F., Polcher, J., and Dupont, J.-C., Response 
!! of land surface fluxes and precipitation to different soil bottom hydrological conditions in a
!! general circulation model,  J. Geophys. Res, in press, 2013. \n
!! - Gouttevin, I., Krinner, G., Ciais, P., Polcher, J., and Legout, C. , 2012. Multi-scale validation
!! of a new soil freezing scheme for a land-surface model with physically-based hydrology.
!! The Cryosphere, 6, 407-430, doi: 10.5194/tc-6-407-2012. \n
!!
!! SVN          :
!! $HeadURL$
!! $Date$
!! $Revision$
!! \n
!_ ===============================================================================================\n
MODULE hydrol

  USE ioipsl
  USE xios_orchidee
  USE constantes
  USE time, ONLY : one_day, dt_sechiba, julian_diff
  USE constantes_soil
  USE pft_parameters
  USE sechiba_io_p
  USE grid
  USE explicitsnow

!pss:+USE TOPMODEL routines
  USE ioipsl_para
  USE topmodel
!pss:-
  USE interpol_help


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: hydrol_main, hydrol_initialize, hydrol_finalize, hydrol_clear, hydrol_rotation_update

  !
  ! variables used inside hydrol module : declaration and initialisation
  !
  LOGICAL, SAVE                                   :: doponds=.FALSE.           !! Reinfiltration flag (true/false)
!$OMP THREADPRIVATE(doponds)

  REAL(r_std), SAVE                               :: froz_frac_corr            !! Coefficient for water frozen fraction correction
!$OMP THREADPRIVATE(froz_frac_corr)
  REAL(r_std), SAVE                               :: max_froz_hydro            !! Coefficient for water frozen fraction correction
!$OMP THREADPRIVATE(max_froz_hydro)
  REAL(r_std), SAVE                               :: smtot_corr                !! Coefficient for water frozen fraction correction
!$OMP THREADPRIVATE(smtot_corr)
  LOGICAL, SAVE                                   :: do_rsoil=.FALSE.          !! Flag to calculate rsoil for bare soile evap 
                                                                               !! (true/false)
!$OMP THREADPRIVATE(do_rsoil)
  LOGICAL, SAVE                                   :: ok_dynroot                !! Flag to activate dynamic root profile to optimize soil  
                                                                               !! moisture usage, similar to Beer et al.2007
!$OMP THREADPRIVATE(ok_dynroot)
  CHARACTER(LEN=80) , SAVE                        :: var_name                  !! To store variables names for I/O
!$OMP THREADPRIVATE(var_name)
  !
  REAL(r_std), PARAMETER                          :: allowed_err =  2.0E-8_r_std
  REAL(r_std), PARAMETER                          :: EPS1 = EPSILON(un)      !! A small number
  ! one dimension array allocated, computed, saved and got in hydrol module
  ! Values per soil type
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: nvan                !! Van Genuchten coeficients n (unitless)
                                                                          ! RK: 1/n=1-m
!$OMP THREADPRIVATE(nvan)                                                 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: avan                !! Van Genuchten coeficients a
                                                                         !!  @tex $(mm^{-1})$ @endtex
!$OMP THREADPRIVATE(avan)                                                
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: mcr                 !! Residual volumetric water content
                                                                         !!  @tex $(m^{3} m^{-3})$ @endtex
!$OMP THREADPRIVATE(mcr)                                                 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: mcs_mineral         !! Saturated volumetric water content
                                                                         !!  @tex $(m^{3} m^{-3})$ @endtex
!$OMP THREADPRIVATE(mcs_mineral)                                                  
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: ks                  !! Hydraulic conductivity at saturation
                                                                         !!  @tex $(mm d^{-1})$ @endtex
!$OMP THREADPRIVATE(ks)                                                  
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: pcent               !! Fraction of saturated volumetric soil moisture above 
                                                                         !! which transpir is max (0-1, unitless)
!$OMP THREADPRIVATE(pcent)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: mcfc_mineral        !! Volumetric water content at field capacity
                                                                         !!  @tex $(m^{3} m^{-3})$ @endtex 
!$OMP THREADPRIVATE(mcfc_mineral)                                                 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: mcw_mineral         !! Volumetric water content at wilting point
                                                                         !!  @tex $(m^{3} m^{-3})$ @endtex 
!$OMP THREADPRIVATE(mcw_mineral)                                                 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: mcs                
!$OMP THREADPRIVATE(mcs)                                                  
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: mcfc                !! Volumetric water content at field capacity
!$OMP THREADPRIVATE(mcfc)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: mcw                 !! Volumetric water content at wilting point
!$OMP THREADPRIVATE(mcw)                                                  
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: VG_m
!$OMP THREADPRIVATE(VG_m)                                                  
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: VG_n
!$OMP THREADPRIVATE(VG_n)                                                  
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: VG_alpha
!$OMP THREADPRIVATE(VG_alpha)                                                  
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: VG_psi_fc
!$OMP THREADPRIVATE(VG_psi_fc)                                                  
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: VG_psi_wp
!$OMP THREADPRIVATE(VG_psi_wp)                                                  

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: mc_awet             !! Vol. wat. cont. above which albedo is cst
                                                                         !!  @tex $(m^{3} m^{-3})$ @endtex 
!$OMP THREADPRIVATE(mc_awet)                                             
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: mc_adry             !! Vol. wat. cont. below which albedo is cst
                                                                         !!  @tex $(m^{3} m^{-3})$ @endtex 
!$OMP THREADPRIVATE(mc_adry)                                             
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: tot_watveg_beg   !! Total amount of water on vegetation at start of time 
                                                                         !! step @tex $(kg m^{-2})$ @endtex 
!$OMP THREADPRIVATE(tot_watveg_beg)                                      
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: tot_watveg_end   !! Total amount of water on vegetation at end of time step
                                                                         !!  @tex $(kg m^{-2})$ @endtex 
!$OMP THREADPRIVATE(tot_watveg_end)                                      
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: tot_watsoil_beg  !! Total amount of water in the soil at start of time step
                                                                         !!  @tex $(kg m^{-2})$ @endtex 
!$OMP THREADPRIVATE(tot_watsoil_beg)                                     
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: tot_watsoil_end  !! Total amount of water in the soil at end of time step
                                                                         !!  @tex $(kg m^{-2})$ @endtex 
!$OMP THREADPRIVATE(tot_watsoil_end)                                     
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: snow_beg         !! Total amount of snow at start of time step
                                                                         !!  @tex $(kg m^{-2})$ @endtex 
!$OMP THREADPRIVATE(snow_beg)                                            
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: snow_end         !! Total amount of snow at end of time step
                                                                         !!  @tex $(kg m^{-2})$ @endtex 
!$OMP THREADPRIVATE(snow_end)                                           
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: delsoilmoist     !! Change in soil moisture @tex $(kg m^{-2})$ @endtex 
!$OMP THREADPRIVATE(delsoilmoist)                                         
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: delintercept     !! Change in interception storage
                                                                         !!  @tex $(kg m^{-2})$ @endtex 
!$OMP THREADPRIVATE(delintercept)                                        
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: delswe           !! Change in SWE @tex $(kg m^{-2})$ @endtex 
!$OMP THREADPRIVATE(delswe)
  REAL(r_std),ALLOCATABLE, SAVE, DIMENSION (:)       :: undermcr         !! Nb of tiles under mcr for a given time step
!$OMP THREADPRIVATE(undermcr) 

  ! array allocated, computed, saved and got in hydrol module
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: mask_veget       !! zero/one when veget fraction is zero/higher (1)
!$OMP THREADPRIVATE(mask_veget)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: mask_soiltile    !! zero/one where soil tile is zero/higher (1)
!$OMP THREADPRIVATE(mask_soiltile)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: humrelv          !! Water stress index for transpiration
                                                                         !! for each soiltile x PFT couple (0-1, unitless)
!$OMP THREADPRIVATE(humrelv)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: vegstressv       !! Water stress index for vegetation growth
                                                                         !! for each soiltile x PFT couple (0-1, unitless)
!$OMP THREADPRIVATE(vegstressv)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:,:):: us               !! Water stress index for transpiration 
                                                                         !! (by soil layer and PFT) (0-1, unitless)
!$OMP THREADPRIVATE(us)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: precisol         !! Throughfall+Totmelt per PFT 
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(precisol)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: throughfall      !! Throughfall per PFT 
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(throughfall)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: precisol_ns      !! Throughfall per soiltile 
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(precisol_ns)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: ae_ns            !! Bare soil evaporation per soiltile
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(ae_ns)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: evap_bare_lim_ns !! Limitation factor (beta) for bare soil evaporation 
                                                                         !! per soiltile (used to deconvoluate vevapnu)  
                                                                         !!  (0-1, unitless) 
!$OMP THREADPRIVATE(evap_bare_lim_ns)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: free_drain_coef  !! Coefficient for free drainage at bottom
                                                                         !!  (0-1, unitless) 
!$OMP THREADPRIVATE(free_drain_coef) 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: zwt_force        !! Prescribed water table depth (m)
!$OMP THREADPRIVATE(zwt_force) 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: frac_bare_ns     !! Evaporating bare soil fraction per soiltile
                                                                         !!  (0-1, unitless)
!$OMP THREADPRIVATE(frac_bare_ns)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: rootsink         !! Transpiration sink by soil layer and soiltile
                                                                         !! @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(rootsink)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: subsnowveg       !! Sublimation of snow on vegetation 
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(subsnowveg)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: subsnownobio     !! Sublimation of snow on other surface types  
                                                                         !! (ice, lakes,...) @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(subsnownobio)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: snowmelt          !! Quantite de glace fondue
!$OMP THREADPRIVATE(snowmelt)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: icemelt          !! Quantite de glace fondue
!$OMP THREADPRIVATE(icemelt)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: subsinksoil      !! Excess of sublimation as a sink for the soil
                                                                         !! @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(subsinksoil)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: vegtot           !! Total Total fraction of grid-cell covered by PFTs
                                                                         !! (bare soil + vegetation) (1; 1)
!$OMP THREADPRIVATE(vegtot)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: resdist          !! Soiltile values from previous time-step (1; 1)
!$OMP THREADPRIVATE(resdist)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: vegtot_old       !! Total Total fraction of grid-cell covered by PFTs
                                                                         !! from previous time-step (1; 1)
!$OMP THREADPRIVATE(vegtot_old)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: mx_eau_var       !! Maximum water content of the soil @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(mx_eau_var)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: irrig_fin        !! final application of irrigation water
!$OMP THREADPRIVATE(irrig_fin)

  ! arrays used by cwrr scheme
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: nroot            !! Normalized root length fraction in each soil layer 
                                                                         !! (0-1, unitless)
                                                                         !! DIM = kjpindex * nvm * nslm
!$OMP THREADPRIVATE(nroot)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: kfact_root       !! Factor to increase Ks towards the surface
                                                                         !! (unitless)
                                                                         !! DIM = kjpindex * nslm * nstm
!$OMP THREADPRIVATE(kfact_root)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: kfact            !! Factor to reduce Ks with depth (unitless)
                                                                         !! DIM = nslm * nscm
!$OMP THREADPRIVATE(kfact)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: zz               !! Depth of nodes [znh in vertical_soil] transformed into (mm) 
!$OMP THREADPRIVATE(zz)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: dz               !! Internode thickness [dnh in vertical_soil] transformed into (mm)
!$OMP THREADPRIVATE(dz)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: dh               !! Layer thickness [dlh in vertical_soil] transformed into (mm)
!$OMP THREADPRIVATE(dh)
  INTEGER(i_std), SAVE                               :: itopmax          !! Number of layers where the node is above 0.1m depth
!$OMP THREADPRIVATE(itopmax)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: mc_lin   !! 50 Vol. Wat. Contents to linearize K and D, for each texture 
                                                                 !!  @tex $(m^{3} m^{-3})$ @endtex
                                                                 !! DIM = imin:imax * nscm
!$OMP THREADPRIVATE(mc_lin)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: k_lin    !! 50 values of unsaturated K, for each soil layer and texture
                                                                 !!  @tex $(mm d^{-1})$ @endtex
                                                                 !! DIM = imin:imax * nslm * nscm
!$OMP THREADPRIVATE(k_lin)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: d_lin    !! 50 values of diffusivity D, for each soil layer and texture
                                                                 !!  @tex $(mm^2 d^{-1})$ @endtex
                                                                 !! DIM = imin:imax * nslm * nscm
!$OMP THREADPRIVATE(d_lin)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: a_lin    !! 50 values of the slope in K=a*mc+b, for each soil layer and texture
                                                                 !!  @tex $(mm d^{-1})$ @endtex
                                                                 !! DIM = imin:imax * nslm * nscm
!$OMP THREADPRIVATE(a_lin)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: b_lin    !! 50 values of y-intercept in K=a*mc+b, for each soil layer and texture
                                                                 !!  @tex $(m^{3} m^{-3})$ @endtex
                                                                 !! DIM = imin:imax * nslm * nscm
!$OMP THREADPRIVATE(b_lin)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: humtot   !! Total Soil Moisture @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(humtot)
  LOGICAL, ALLOCATABLE, SAVE, DIMENSION (:)          :: resolv   !! Mask of land points where to solve the diffusion equation
                                                                 !! (true/false)
!$OMP THREADPRIVATE(resolv)

!! for output
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: kk_moy   !! Mean hydraulic conductivity over soiltiles (mm/d)
!$OMP THREADPRIVATE(kk_moy)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: kk       !! Hydraulic conductivity for each soiltiles (mm/d)
!$OMP THREADPRIVATE(kk)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: avan_mod_tab  !! VG parameter a modified from  exponantial profile
                                                                      !! @tex $(mm^{-1})$ @endtex !! DIMENSION (nslm,nscm)
!$OMP THREADPRIVATE(avan_mod_tab)  
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: nvan_mod_tab  !! VG parameter n  modified from  exponantial profile
                                                                      !! (unitless) !! DIMENSION (nslm,nscm)  
!$OMP THREADPRIVATE(nvan_mod_tab)
  
!! linarization coefficients of hydraulic conductivity K (hydrol_soil_coef)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: k        !! Hydraulic conductivity K for each soil layer
                                                                 !!  @tex $(mm d^{-1})$ @endtex
                                                                 !! DIM = (:,nslm)
!$OMP THREADPRIVATE(k)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: a        !! Slope in K=a*mc+b(:,nslm)
                                                                 !!  @tex $(mm d^{-1})$ @endtex
                                                                 !! DIM = (:,nslm)
!$OMP THREADPRIVATE(a)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: b        !! y-intercept in K=a*mc+b
                                                                 !!  @tex $(m^{3} m^{-3})$ @endtex
                                                                 !! DIM = (:,nslm)
!$OMP THREADPRIVATE(b)
!! linarization coefficients of hydraulic diffusivity D (hydrol_soil_coef)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: d        !! Diffusivity D for each soil layer
                                                                 !!  @tex $(mm^2 d^{-1})$ @endtex
                                                                 !! DIM = (:,nslm)
!$OMP THREADPRIVATE(d)
!! matrix coefficients (hydrol_soil_tridiag and hydrol_soil_setup), see De Rosnay (1999), p155-157
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: e        !! Left-hand tridiagonal matrix coefficients
!$OMP THREADPRIVATE(e)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: f        !! Left-hand tridiagonal matrix coefficients
!$OMP THREADPRIVATE(f)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: g1       !! Left-hand tridiagonal matrix coefficients
!$OMP THREADPRIVATE(g1)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: ep       !! Right-hand matrix coefficients
!$OMP THREADPRIVATE(ep)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: fp       !! Right-hand atrix coefficients
!$OMP THREADPRIVATE(fp)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: gp       !! Right-hand atrix coefficients
!$OMP THREADPRIVATE(gp)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: rhs      !! Right-hand system
!$OMP THREADPRIVATE(rhs)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: srhs     !! Temporarily stored rhs
!$OMP THREADPRIVATE(srhs)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: tmat             !! Left-hand tridiagonal matrix
!$OMP THREADPRIVATE(tmat)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: stmat            !! Temporarily stored tmat
  !$OMP THREADPRIVATE(stmat)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: water2infilt     !! Water to be infiltrated
                                                                         !! @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(water2infilt)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmc              !! Total moisture content per soiltile 
                                                                         !!  @tex $(kg m^{-2})$ @endtex 
!$OMP THREADPRIVATE(tmc)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmcr             !! Total moisture content at residual per soiltile 
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(tmcr)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmcs             !! Total moisture content at saturation per soiltile 
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(tmcs)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmcfc            !! Total moisture content at field capacity per soiltile 
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(tmcfc)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmcw             !! Total moisture content at wilting point per soiltile 
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(tmcw)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmc_litter       !! Total moisture in the litter per soiltile
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(tmc_litter)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: tmc_litt_mea     !! Total moisture in the litter over the grid
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(tmc_litt_mea)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmc_litter_wilt  !! Total moisture of litter at wilt point per soiltile 
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(tmc_litter_wilt)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmc_litter_field !! Total moisture of litter at field cap. per soiltile 
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(tmc_litter_field)
!!! A CHANGER DANS TOUT HYDROL: tmc_litter_res et sat ne devraient pas dependre de ji - tdo
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmc_litter_res   !! Total moisture of litter at residual moisture per soiltile 
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(tmc_litter_res)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmc_litter_sat   !! Total moisture of litter at saturation per soiltile 
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(tmc_litter_sat)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmc_litter_awet  !! Total moisture of litter at mc_awet per soiltile 
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(tmc_litter_awet)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmc_litter_adry  !! Total moisture of litter at mc_adry per soiltile 
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(tmc_litter_adry)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: tmc_litt_wet_mea !! Total moisture in the litter over the grid below which 
                                                                         !! albedo is fixed constant
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(tmc_litt_wet_mea)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: tmc_litt_dry_mea !! Total moisture in the litter over the grid above which 
                                                                         !! albedo is constant
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(tmc_litt_dry_mea)
  LOGICAL, SAVE                                      :: tmc_init_updated = .FALSE. !! Flag allowing to determine if tmc is initialized.
!$OMP THREADPRIVATE(tmc_init_updated)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: v1               !! Temporary variable (:)
!$OMP THREADPRIVATE(v1)

  !! par type de sol :
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: ru_ns            !! Surface runoff per soiltile
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(ru_ns)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: dr_ns            !! Drainage per soiltile
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(dr_ns)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tr_ns            !! Transpiration per soiltile
!$OMP THREADPRIVATE(tr_ns)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: vegetmax_soil    !! (:,nvm,nstm) percentage of each veg. type on each soil 
                                                                         !! of each grid point 
!$OMP THREADPRIVATE(vegetmax_soil)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: mc               !! Total volumetric water content at the calculation nodes
                                                                         !! (eg : liquid + frozen)
                                                                         !!  @tex $(m^{3} m^{-3})$ @endtex
!$OMP THREADPRIVATE(mc)

   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:) :: mc_read_prev       !! Soil moisture from file at previous timestep in the file
!$OMP THREADPRIVATE(mc_read_prev)
   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:) :: mc_read_next       !! Soil moisture from file at next time step in the file
!$OMP THREADPRIVATE(mc_read_next)
   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:) :: mask_mc_interp     !! Mask of valid data in soil moisture nudging file
!$OMP THREADPRIVATE(mask_mc_interp)
   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)   :: snowdz_read_prev   !! snowdz read from file at previous timestep in the file
!$OMP THREADPRIVATE(snowdz_read_prev)
   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)   :: snowdz_read_next   !! snowdz read from file at next time step in the file
!$OMP THREADPRIVATE(snowdz_read_next)
   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)   :: snowrho_read_prev  !! snowrho read from file at previous timestep in the file
!$OMP THREADPRIVATE(snowrho_read_prev)
   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)   :: snowrho_read_next  !! snowrho read from file at next time step in the file
!$OMP THREADPRIVATE(snowrho_read_next)
   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)   :: snowtemp_read_prev !! snowtemp read from file at previous timestep in the file
!$OMP THREADPRIVATE(snowtemp_read_prev)
   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)   :: snowtemp_read_next !! snowtemp read from file at next time step in the file
!$OMP THREADPRIVATE(snowtemp_read_next)
   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)   :: mask_snow_interp   !! Mask of valid data in snow nudging file
!$OMP THREADPRIVATE(mask_snow_interp)

   REAL(r_std),ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: mcl              !! Liquid water content
                                                                         !!  @tex $(m^{3} m^{-3})$ @endtex
!$OMP THREADPRIVATE(mcl)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: soilmoist        !! (:,nslm) Mean of each soil layer's moisture
                                                                         !! across soiltiles
                                                                         !!  @tex $(kg m^{-2})$ @endtex 
!$OMP THREADPRIVATE(soilmoist)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: soilmoist_liquid !! (:,nslm) Mean of each soil layer's liquid moisture
                                                                         !! across soiltiles
                                                                         !!  @tex $(kg m^{-2})$ @endtex 
!$OMP THREADPRIVATE(soilmoist_liquid)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: soil_wet_ns      !! Soil wetness above mcw (0-1, unitless)
!$OMP THREADPRIVATE(soil_wet_ns)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: soil_wet_litter  !! Soil wetness aove mvw in the litter (0-1, unitless)
!$OMP THREADPRIVATE(soil_wet_litter)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: qflux            !! Diffusive water fluxes between soil layers
!$OMP THREADPRIVATE(qflux)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: profil_froz_hydro     !! Frozen fraction for each hydrological soil layer
!$OMP THREADPRIVATE(profil_froz_hydro)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: profil_froz_hydro_ns  !! As  profil_froz_hydro per soiltile
!$OMP THREADPRIVATE(profil_froz_hydro_ns)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: temp_hydro            !! Temp profile on hydrological levels
!$OMP THREADPRIVATE(temp_hydro)

!gmjc top 5 layer soil moisture for grazing
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmc_trampling
!$OMP THREADPRIVATE(tmc_trampling)
!end gmjc
  REAL(r_std), ALLOCATABLE, SAVE,DIMENSION(:)    :: refSOC_1d            !!initialize soil organic carbon only used to calculate thermal insulating effect
!$OMP THREADPRIVATE(refSOC_1d)
  LOGICAL, SAVE                                  :: use_refSOC_hydrol = .TRUE. !! consider impacts of soil organic carbon on hydrologic parameters
!$OMP THREADPRIVATE(use_refSOC_hydrol)


CONTAINS

!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_initialize
!!
!>\BRIEF         Allocate module variables, read from restart file or initialize with default values
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

  SUBROUTINE hydrol_initialize ( kjit,           kjpindex,  index,         rest_id,   &
                                 njsc,           soiltile,  veget,         veget_max,        &
                                 humrel,         vegstress, drysoil_frac,                    &
                                 shumdiag_perma,    qsintveg,                        &
                                 evap_bare_lim,  snow,      snow_age,      snow_nobio,       &
!                                 evap_bare_lim,  evap_bare_lim_pft, snow,      snow_age,      snow_nobio,       &
                                 snow_nobio_age, snowrho,   snowtemp,                        &
                                 snowgrain,      snowdz,    snowheat,      &
                                 totfrac_nobio,  precip_rain, precip_snow, returnflow,       &
                                 reinfiltration, irrigation,vevapwet,                        &
                                 transpir,       vevapnu,   vevapsno,      vevapflo,         &
                                 floodout,       runoff,    drainage,                        &
                                 mc_layh,        mcl_layh,  soilmoist_out,                   &
                                 mc_layh_s,      mcl_layh_s,                                 &
!gmjc
                                 tmc_topgrass, humcste_use, altmax, swdown) !!xchen add swdown
!end gmjc

    !! 0. Variable and parameter declaration
    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                         :: kjit             !! Time step number 
    INTEGER(i_std), INTENT(in)                         :: kjpindex         !! Domain size
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)   :: index            !! Indeces of the points on the map
    INTEGER(i_std),INTENT (in)                         :: rest_id          !! Restart file identifier
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)   :: njsc             !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
    REAL(r_std),DIMENSION (kjpindex,nstm), INTENT (in) :: soiltile         !! Fraction of each soil tile within vegtot (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: veget            !! Fraction of vegetation type           
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: veget_max        !! Max. fraction of vegetation type (LAI -> infty)

    !! 0.2 Output variables
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)   :: humrel         !! Relative humidity
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)   :: vegstress      !! Veg. moisture stress (only for vegetation growth)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)       :: drysoil_frac   !! function of litter wetness
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (out)  :: shumdiag_perma !! Percent of porosity filled with water (mc/mcs) used for the thermal computations
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)   :: qsintveg       !! Water on vegetation due to interception
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)        :: evap_bare_lim  !! Limitation factor for bare soil evaporation
!    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(out)    :: evap_bare_lim_pft  !! Limitation factor for bare soil evaporation
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)       :: snow           !! Snow mass [Kg/m^2]
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)       :: snow_age       !! Snow age
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (out):: snow_nobio     !! Water balance on ice, lakes, .. [Kg/m^2]
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (out):: snow_nobio_age !! Snow age on ice, lakes, ...
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(out) :: snowrho        !! Snow density
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(out) :: snowtemp       !! Snow temperature
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(out) :: snowgrain      !! Snow grainsize
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(out) :: snowdz         !! Snow layer thickness
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(out) :: snowheat       !! Snow heat content
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (out)  :: mc_layh        !! Volumetric moisture content for each layer in hydrol (liquid+ice) m3/m3
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (out)  :: mcl_layh       !! Volumetric moisture content for each layer in hydrol (liquid) m3/m3
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (out)  :: soilmoist_out  !! Total soil moisture content for each layer in hydrol (liquid+ice), mm
    REAL(r_std),DIMENSION (kjpindex,nslm,nstm), INTENT (out)  :: mc_layh_s        !! Volumetric moisture content for each layer in hydrol (liquid+ice) m3/m3
    REAL(r_std),DIMENSION (kjpindex,nslm,nstm), INTENT (out)  :: mcl_layh_s       !! Volumetric moisture content for each layer in hydrol (liquid) m3/m3
    REAL(r_std),DIMENSION (kjpindex)                     :: soilwetdummy   !! Temporary variable never used

!gmjc
    REAL(r_std),DIMENSION (kjpindex),INTENT(out)         :: tmc_topgrass
!end gmjc
    REAL(r_std),DIMENSION (kjpindex,nvm),INTENT(out)         :: humcste_use
    REAL(r_std),DIMENSION (kjpindex,nvm),INTENT(in)         :: altmax
    INTEGER(i_std)                                      :: ier
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: swdown           !!xchen Down-welling surface short-wave flux  

    !! 0.4 Local variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: totfrac_nobio    !! Total fraction of ice+lakes+...
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: precip_rain      !! Rain precipitation
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: precip_snow      !! Snow precipitation
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: returnflow       !! Routed water which comes back into the soil (from the 
                                                                             !! bottom) 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: reinfiltration   !! Routed water which comes back into the soil (at the 
                                                                             !! top) 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: irrigation       !! Water from irrigation returning to soil moisture  
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: vevapwet         !! Interception loss
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: transpir         !! Transpiration
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)   :: vevapnu          !! Bare soil evaporation
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)   :: vevapsno         !! Snow evaporation
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)   :: vevapflo         !! Floodplain evaporation
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)       :: floodout     !! Flux out of floodplains
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)       :: runoff       !! Complete runoff
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)       :: drainage     !! Drainage

    INTEGER(i_std)                                       :: jsl
    REAL(r_std), DIMENSION(kjpindex,1) :: refsoc1d_yat !! debug array
    CHARACTER(LEN=150) :: filename
!_ ================================================================================================================================

    !Config Key   = use_refSOC_hydrol
    !Config Desc  = 
    !Config Def   = 
    !Config If    = 
    !Config Help  = 
    !Config         
    !Config Units = 
    !
    CALL getin_p('use_refSOC_hydrol',use_refSOC_hydrol)

    IF (use_refSOC_hydrol) THEN
      ALLOCATE (refSOC_1d(kjpindex),stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3,'hydrol_initialize', 'Error in allocation of refSOC_1d','','')
  
      CALL restget_p (rest_id, 'refSOC_1d', nbp_glo, 1, 1, kjit, .TRUE., refSOC_1d, "gather", nbp_glo, index_g)

      IF (ALL(refSOC_1d(:) == val_exp)) THEN
          ! Open file with SOC map

           !Config Key   = SOIL_REFSOC_1d_FILE
           !Config Desc  = File with climatological soil temperature
           !Config If    = READ_REFTEMP
           !Config Def   = reftemp.nc
           !Config Help  = 
           !Config Units = [FILE]
           filename = 'refSOC_1d.nc'
           CALL getin_p('SOIL_REFSOC_1d_FILE',filename)

           CALL read_refSOCfile_1dfile(kjpindex, filename, 'soil_organic_carbon_1d', lalo, refsoc1d_yat)
           refSOC_1d(:) = refsoc1d_yat(:,1)
      ENDIF
    ENDIF

!pss+       
    !Config Key   = TOPM_CALCUL
    !Config Desc  = Enable or disable TOPMODEL module
    !Config Def   = False
    !Config If    = OK_SECHIBA
    !Config Help  = Enable or disable TOPMODEL module.
    !Config         
    !Config Units = [-] 
    TOPM_calcul  = .FALSE.
    CALL getin_p('TOPM_CALCUL', TOPM_calcul)
!pss+

    CALL hydrol_init (kjit, kjpindex, index, rest_id, veget_max, soiltile, &
         humrel, vegstress, snow,       snow_age,   snow_nobio, snow_nobio_age, qsintveg, &
         snowdz, snowgrain, snowrho,    snowtemp,   snowheat, &
         drysoil_frac, evap_bare_lim)
!         snowflx,snowcap,   cgrnd_snow, dgrnd_snow, drysoil_frac, evap_bare_lim, evap_bare_lim_pft, &
    
    CALL hydrol_var_init (kjpindex, veget, veget_max, &
         soiltile, njsc, mx_eau_var, shumdiag_perma, &
         drysoil_frac, qsintveg, mc_layh, mcl_layh, &
         mc_layh_s, mcl_layh_s, &
!gmjc
         tmc_topgrass, humcste_use,altmax, swdown) !!xchen add swdown
!end gmjc 

!pss+       
    IF (TOPM_calcul) CALL topmodel_initialize (kjit, kjpindex, nbp_glo, lalo, veget_max, zmaxh, &
                                rest_id, index_g, &
                                mcs, mcw, soiltile)
!pss:-

      !! Initialize alma output variables if they were not found in the restart file. This is done in the end of 
      !! hydrol_initialize so that all variables(humtot,..) that will be used are initialized.
      IF (ALL(tot_watveg_beg(:)==val_exp) .OR.  ALL(tot_watsoil_beg(:)==val_exp) .OR. ALL(snow_beg(:)==val_exp)) THEN
            ! The output variable soilwetdummy is not calculated at first call to hydrol_alma.
            CALL hydrol_alma(kjpindex, index, .TRUE., qsintveg, snow, snow_nobio, soilwetdummy)
       END IF

    !! Calculate itopmax indicating the number of layers where the node is above 0.1m depth
    itopmax=1
    DO jsl = 1, nslm
       ! znh : depth of nodes
       IF (znh(jsl) <= 0.1) THEN
          itopmax=jsl
       END IF
    END DO
    IF (printlev>=3) WRITE(numout,*) "Number of layers where the node is above 0.1m depth: itopmax=",itopmax

    ! Copy soilmoist into a local variable to be sent to thermosoil
    soilmoist_out(:,:) = soilmoist(:,:)

  END SUBROUTINE hydrol_initialize


!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_main
!!
!>\BRIEF         
!!
!! DESCRIPTION :
!! - called every time step
!! - initialization and finalization part are not done in here
!!
!! - 1 computes snow  ==> hydrol_snow
!! - 2 computes vegetations reservoirs  ==> hydrol_vegupd
!! - 3 computes canopy  ==> hydrol_canop
!! - 4 computes surface reservoir  ==> hydrol_flood
!! - 5 computes soil hydrology ==> hydrol_soil
!!
!! IMPORTANT NOTICE : The water fluxes are used in their integrated form, over the time step 
!! dt_sechiba, with a unit of kg m^{-2}.
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE hydrol_main (kjit, kjpindex, &
       & index, indexveg, indexsoil, indexlayer, indexnslm, &
       & temp_sol_new, floodout, runoff, drainage, frac_nobio, totfrac_nobio, vevapwet, veget, veget_max, njsc, &
       & qsintmax, qsintveg, vevapnu, vevapnu_pft, vevapsno, vevapflo, snow, snow_age, snow_nobio, snow_nobio_age,  &
       & tot_melt, transpir, precip_rain, precip_snow, returnflow, reinfiltration, irrigation, vegstress_old, transpot, &
       & humrel, vegstress, drysoil_frac, evapot, evapot_penm, evap_bare_lim, flood_frac, flood_res, &
!       & humrel, vegstress, drysoil_frac, evapot, evapot_penm, evap_bare_lim, evap_bare_lim_pft, flood_frac, flood_res, &
       & shumdiag,shumdiag_perma, k_litt, litterhumdiag, soilcap, soiltile, fraclut, reinf_slope, rest_id, hist_id, hist2_id, soil_deficit, is_crop_soil, &
       & contfrac, stempdiag, &
       & temp_air, pb, u, v, tq_cdrag, pgflux, &
       & snowrho,snowtemp, snowgrain,snowdz,snowheat,snowliq,&
       & grndflux,gtemp, tot_bare_soil, &
       & soilflxresid, mc_layh, mcl_layh, soilmoist_out, &
       & mc_layh_s, mcl_layh_s, & 
       & lambda_snow, cgrnd_snow, dgrnd_snow, frac_snow_veg, temp_sol_add, &
!gmjc
       & tmc_topgrass, humcste_use, swdown) !!xchen add swdown
!end gmjc

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables
  
    INTEGER(i_std), INTENT(in)                         :: kjit             !! Time step number 
    INTEGER(i_std), INTENT(in)                         :: kjpindex         !! Domain size
    INTEGER(i_std),INTENT (in)                         :: rest_id,hist_id  !! _Restart_ file and _history_ file identifier
    INTEGER(i_std),INTENT (in)                         :: hist2_id         !! _history_ file 2 identifier
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)   :: index            !! Indeces of the points on the map
    INTEGER(i_std),DIMENSION (kjpindex*nvm), INTENT (in):: indexveg        !! Indeces of the points on the 3D map for veg
    INTEGER(i_std),DIMENSION (kjpindex*nstm), INTENT (in):: indexsoil      !! Indeces of the points on the 3D map for soil
    INTEGER(i_std),DIMENSION (kjpindex*nslm), INTENT (in):: indexlayer     !! Indeces of the points on the 3D map for soil layers
    INTEGER(i_std),DIMENSION (kjpindex*nslm), INTENT (in):: indexnslm      !! Indeces of the points on the 3D map for of diagnostic soil layers

    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: precip_rain      !! Rain precipitation
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: precip_snow      !! Snow precipitation
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: returnflow       !! Routed water which comes back into the soil (from the 
                                                                           !! bottom) 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: reinfiltration   !! Routed water which comes back into the soil (at the 
                                                                           !! top) 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: irrigation       !! Water from irrigation returning to soil moisture  
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: vegstress_old    !! vegstress of previous step
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: transpot         !! potential transpiratio

    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: temp_sol_new     !! New soil temperature

    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)   :: njsc             !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (in) :: frac_nobio     !! Fraction of ice, lakes, ...
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: totfrac_nobio    !! Total fraction of ice+lakes+...
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: soilcap          !! Soil capacity
    REAL(r_std),DIMENSION (kjpindex,nstm), INTENT (in) :: soiltile         !! Fraction of each soil tile within vegtot (0-1, unitless)
    LOGICAL, DIMENSION (nstm), INTENT (in)             :: is_crop_soil     !! whether soil tile is under cropland
    REAL(r_std),DIMENSION (kjpindex,nlut), INTENT (in) :: fraclut          !! Fraction of each landuse tile (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: vevapwet         !! Interception loss
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: veget            !! Fraction of vegetation type           
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: veget_max        !! Max. fraction of vegetation type (LAI -> infty)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: qsintmax         !! Maximum water on vegetation for interception
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: transpir         !! Transpiration
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: reinf_slope      !! Slope coef
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: evapot           !! Soil Potential Evaporation
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: evapot_penm      !! Soil Potential Evaporation Correction
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: flood_frac       !! flood fraction
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: contfrac         !! Fraction of continent in the grid
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in) :: stempdiag        !! Diagnostic temp profile from thermosoil
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)       :: temp_air         !! Air temperature
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)       :: u,v              !! Horizontal wind speed
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)       :: tq_cdrag         !! Surface drag coefficient (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)       :: pb               !! Surface pressure
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)       :: pgflux           !! Net energy into snowpack
    REAL(r_std),DIMENSION (kjpindex),INTENT(inout)     :: soilflxresid     !! Energy flux to the snowpack
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)       :: gtemp            !! First soil layer temperature
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)       :: tot_bare_soil    !! Total evaporating bare soil fraction 
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)       :: lambda_snow      !! Coefficient of the linear extrapolation of surface temperature 
    REAL(r_std),DIMENSION (kjpindex,nsnow), INTENT (in):: cgrnd_snow       !! Integration coefficient for snow numerical scheme
    REAL(r_std),DIMENSION (kjpindex,nsnow), INTENT (in):: dgrnd_snow       !! Integration coefficient for snow numerical scheme
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)       :: frac_snow_veg    !! Snow cover fraction on vegetation
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: swdown           !!xchen Down-welling surface short-wave flux    

    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: vegstress        !! Veg. moisture stress (only for vegetation growth)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: drysoil_frac     !! function of litter wetness
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (out):: shumdiag         !! Relative soil moisture in each soil layer 
                                                                           !! with respect to (mcfc-mcw)
                                                                           !! (unitless; can be out of 0-1)
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (out):: shumdiag_perma   !! Percent of porosity filled with water (mc/mcs) used for the thermal computations 
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: k_litt           !! litter approximate conductivity
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: litterhumdiag    !! litter humidity
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: tot_melt         !! Total melt    
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: floodout         !! Flux out of floodplains
    
!pss:+
!    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)     :: drunoff_tot    !! Dunne runoff 
!    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: fwet_out   !! fwet: to change the balance of energy according to wetland fraction
!pss:-

    !! 0.3 Modified variables

    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(inout):: qsintveg         !! Water on vegetation due to interception
    REAL(r_std),DIMENSION (kjpindex), INTENT(inout)    :: evap_bare_lim    !! Limitation factor (beta) for bare soil evaporation    
!    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(inout) :: evap_bare_lim_pft    !! Limitation factor (beta) for bare soil evaporation    
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(out) :: soil_deficit    !! water deficit to reach IRRIG_FULFILL of holding capacity
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(inout):: humrel           !! Relative humidity
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)   :: vevapnu          !! Bare soil evaporation
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout) :: vevapnu_pft      !! Bare soil evaporation
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)   :: vevapsno         !! Snow evaporation
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)   :: vevapflo         !! Floodplain evaporation
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)   :: flood_res        !! flood reservoir estimate
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)   :: snow             !! Snow mass [kg/m^2]
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)   :: snow_age         !! Snow age
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (inout) :: snow_nobio  !! Water balance on ice, lakes, .. [Kg/m^2]
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (inout) :: snow_nobio_age !! Snow age on ice, lakes, ...
    !! We consider that any water on the ice is snow and we only peforme a water balance to have consistency.
    !! The water balance is limite to + or - 10^6 so that accumulation is not endless

    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: runoff       !! Complete surface runoff
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: drainage     !! Drainage
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(inout) :: snowrho      !! Snow density
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(inout) :: snowtemp     !! Snow temperature
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(inout) :: snowgrain    !! Snow grainsize
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(inout) :: snowdz       !! Snow layer thickness
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(inout) :: snowheat     !! Snow heat content
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(inout) :: snowliq      !! Snow liquid content (m)
    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)       :: grndflux     !! Net flux into soil W/m2
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT(out)     :: mc_layh      !! Volumetric moisture content for each layer in hydrol(liquid + ice) [m3/m3)]
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT(out)     :: mcl_layh     !! Volumetric moisture content for each layer in hydrol(liquid) [m3/m3]
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT(out)     :: soilmoist_out!! Total soil moisture content for each layer in hydrol(liquid + ice) [mm]
    REAL(r_std),DIMENSION (kjpindex,nslm,nstm), INTENT(out)  :: mc_layh_s      !! Volumetric moisture content for each layer in hydrol(liquid + ice) [m3/m3)]
    REAL(r_std),DIMENSION (kjpindex,nslm,nstm), INTENT(out)  :: mcl_layh_s     !! Volumetric moisture content for each layer in hydrol(liquid) [m3/m3]
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)       :: temp_sol_add !! additional surface temperature due to the melt of first layer
                                                                           !! at the present time-step @tex ($K$) @endtex

!gmjc
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)       :: tmc_topgrass
!end gmjc
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: humcste_use
    !! 0.4 Local variables

    INTEGER(i_std)                                     :: jst              !! Index of soil tiles (unitless, 1-3)
    INTEGER(i_std)                                     :: jsl              !! Index of soil layers (unitless)
    INTEGER(i_std)                                     :: ji, jv
    REAL(r_std)                                        :: tempfrac         !! temporary fraction for irrigation
    REAL(r_std),DIMENSION (kjpindex)                   :: soilwet          !! A temporary diagnostic of soil wetness
    REAL(r_std),DIMENSION (kjpindex)                   :: snowdepth        !! Depth of snow layer, only for diagnostics with ok_explicitsnow=n
    REAL(r_std),DIMENSION (kjpindex)                   :: snowdepth_diag   !! Depth of snow layer containing default values, only for diagnostics
    REAL(r_std),DIMENSION (kjpindex, nsnow)            :: snowdz_diag      !! Depth of snow layer on all layers containing default values, only for diagnostics
    REAL(r_std),DIMENSION (kjpindex)                   :: njsc_tmp         !! Temporary REAL value for njsc to write it
    REAL(r_std),DIMENSION (kjpindex,nvm)               :: irrig_demand_ratio !! irrigation demand for different PFTs
    REAL(r_std), DIMENSION (kjpindex)                  :: snowmelt         !! Snow melt [mm/dt_sechiba]
    REAL(r_std), DIMENSION (kjpindex,nstm)             :: tmc_top          !! Moisture content in the itopmax upper layers, per tile
    REAL(r_std), DIMENSION (kjpindex)                  :: humtot_top       !! Moisture content in the itopmax upper layers, for diagnistics
    REAL(r_std), DIMENSION(kjpindex)                   :: histvar          !! Temporary variable when computations are needed
    REAL(r_std), DIMENSION (kjpindex,nvm)              :: frac_bare        !! Fraction(of veget_max) of bare soil in each vegetation type

    REAL(r_std), DIMENSION(kjpindex)                   :: twbr             !! Grid-cell mean of TWBR Total Water Budget Residu[kg/m2/dt]
    REAL(r_std),DIMENSION (kjpindex,nslm)              :: land_nroot       !! To ouput the grid-cell mean of nroot
    REAL(r_std),DIMENSION (kjpindex,nslm)              :: land_dlh         !! To ouput the soil layer thickness on all grid points [m]
    REAL(r_std),DIMENSION (kjpindex,nslm)              :: land_mcs         !! To ouput the mean of mcs
    REAL(r_std),DIMENSION (kjpindex,nslm)              :: land_mcfc        !! To ouput the mean of mcfc 
    REAL(r_std),DIMENSION (kjpindex,nslm)              :: land_mcw         !! To ouput the mean of mcw 
    REAL(r_std),DIMENSION (kjpindex,nslm)              :: land_mcr         !! To ouput the mean of mcr 
    REAL(r_std),DIMENSION (kjpindex)                   :: land_tmcs        !! To ouput the grid-cell mean of tmcs 
    REAL(r_std),DIMENSION (kjpindex)                   :: land_tmcfc       !! To ouput the grid-cell mean of tmcfc
    REAL(r_std),DIMENSION (kjpindex)                   :: drain_upd        !! Change in drainage due to decrease in vegtot
                                                                           !! on mc [kg/m2/dt]
    REAL(r_std),DIMENSION (kjpindex)                   :: runoff_upd       !! Change in runoff due to decrease in vegtot
                                                                           !! on water2infilt[kg/m2/dt]
    REAL(r_std), DIMENSION (kjpindex,nlut)             :: humtot_lut       !! Moisture content on landuse tiles, for diagnostics
    REAL(r_std), DIMENSION (kjpindex,nlut)             :: humtot_top_lut   !! Moisture content in upper layers on landuse tiles, for diagnostics
    REAL(r_std), DIMENSION (kjpindex,nlut)             :: mrro_lut         !! Total runoff from landuse tiles, for diagnostics

!_ ================================================================================================================================
    !! 1. Update vegtot_old and recalculate vegtot
    vegtot_old(:) = vegtot(:)
    soil_deficit = zero

    DO ji = 1, kjpindex
       vegtot(ji) = SUM(veget_max(ji,:))
    ENDDO

    !! 2. Applay nudging for soil moisture and/or snow variables
    IF (ok_nudge_mc .OR. ok_nudge_snow) THEN
       CALL hydrol_nudge(kjit, kjpindex, mc, snowdz, snowrho, snowtemp, soiltile)
    END IF

    !! 3. Shared time step
    IF (printlev>=3) WRITE (numout,*) 'hydrol pas de temps = ',dt_sechiba

!pss:+ 
    !! 3.0 Calculate wetland fractions
    IF (TOPM_calcul) THEN
      CALL topmodel_main(kjpindex, dz, humtot, profil_froz_hydro, zmaxh)
    ENDIF
!pss:-

    ! 
    !! 3.1 Calculate snow processes with explicit method or bucket snow model
    IF (ok_explicitsnow) THEN
       ! Explicit snow model
       IF (printlev>=3) WRITE (numout,*) ' ok_explicitsnow : use multi-snow layer '
       
       CALL explicitsnow_main(kjpindex,    precip_rain,   precip_snow,    temp_air,   pb,       &
                              u,           v,             temp_sol_new,   soilcap,    pgflux,   &
                              frac_nobio,  totfrac_nobio, gtemp,                                &
                              lambda_snow, cgrnd_snow,    dgrnd_snow,     contfrac,             & 
                              vevapsno,    snow_age,      snow_nobio_age, snow_nobio, snowrho,  &
                              snowgrain,   snowdz,        snowtemp,       snowheat,   snow,     &
                              temp_sol_add,                                                     &
                              snowliq,     subsnownobio,  grndflux,       snowmelt,   tot_melt, &
                              soilflxresid,subsinksoil)
 
    ELSE
       ! Bucket snow model
       CALL hydrol_snow(kjpindex, precip_rain, precip_snow, temp_sol_new, soilcap, &
            frac_nobio, totfrac_nobio, vevapsno, snow, snow_age, snow_nobio, snow_nobio_age, &
            tot_melt, snowdepth,snowmelt)
    END IF
        
    !
    !! 3.2 computes vegetations reservoirs  ==>hydrol_vegupd
! Modif temp vuichard
    CALL hydrol_vegupd(kjpindex, veget, veget_max, soiltile, qsintveg, frac_bare, drain_upd, runoff_upd)

    !! Calculate kfact_root
    !! An exponential factor is used to increase ks near the surface depending on the amount of roots in the soil 
    !! through a geometric average over the vegets
    !! This comes from the PhD thesis of d'Orgeval, 2006, p82; d'Orgeval et al. 2008, Eqs. 3-4
    !! (Calibrated against Hapex-Sahel measurements)
    !! Since rev 2916: veget_max/2 is used instead of veget
    kfact_root(:,:,:) = un
    DO jsl = 1, nslm
       DO jv = 2, nvm
          jst = pref_soil_veg(jv)
          DO ji = 1, kjpindex
             IF(soiltile(ji,jst) .GT. min_sechiba) THEN
                kfact_root(ji,jsl,jst) = kfact_root(ji,jsl,jst) * &
                     & MAX((MAXVAL(ks_usda)/ks(njsc(ji)))**(- vegetmax_soil(ji,jv,jst)/2 * (humcste(jv)*zz(jsl)/mille - un)/deux), &
                     un) 
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    !
    !! 3.3 computes canopy  ==>hydrol_canop
    CALL hydrol_canop(kjpindex, precip_rain, vevapwet, veget_max, veget, qsintmax, qsintveg,precisol,tot_melt)

    !
    !! 3.4 computes surface reservoir  ==>hydrol_flood
    CALL hydrol_flood(kjpindex,  vevapflo, flood_frac, flood_res, floodout)

    !
    !! 3.5 computes soil hydrology ==>hydrol_soil
    !!! calculating ratio of irrigation for each pft at each point
    irrig_demand_ratio(:,:) = zero
!    irrig_totfrac(:) = zero
    DO ji = 1,kjpindex
        DO jv = 2,nvm
            IF (veget_max(ji,jv) .GT. zero) THEN
                IF (irrig_drip) THEN
                    tempfrac = veget(ji,jv)/veget_max(ji,jv)
                    IF ( ok_LAIdev(jv) .AND. (vegstress_old(ji,jv) .LT. irrig_threshold(jv)) .AND.  &
                        & (transpot(ji,jv)*tempfrac + evapot(ji)*(1-tempfrac) .GT. precip_rain(ji)) ) THEN
        !                irrig_totfrac(ji) = irrig_totfrac(ji) + veget_max(ji,jv)
                        irrig_demand_ratio(ji,jv) = MIN( irrig_dosmax, irrig_fulfill(jv) * &
                                                    & ( transpot(ji,jv)*tempfrac &
                                                    & + evapot(ji)*(1-tempfrac) &
                                                    & - precip_rain(ji) ) ) * veget_max(ji,jv)
                        !!!! reconsider if evapot or evapot_corr to be used as irrigation demand
                        !!!! if re-infiltration is considered in sechiba, it
                        !should also be considered here, xuhui
                    ENDIF ! since irrigated ratio is the same across pfts on the same grid, no need to consider
                ELSE ! flooding
                    IF ( ok_LAIdev(jv) .AND. (vegstress_old(ji,jv) .LT. irrig_threshold(jv)) ) THEN 
                        irrig_demand_ratio(ji,jv) = MIN( irrig_dosmax, MAX( zero, soil_deficit(ji,jv) ) ) * veget_max(ji,jv)
                    ENDIF
                    !!!! if re-infiltration is considered in sechiba, it
                        !should also be considered here, xuhui
                ENDIF
            ENDIF
        ENDDO
        IF ( SUM(irrig_demand_ratio(ji,:)) .GT. zero ) THEN
            irrig_demand_ratio(ji,:) = irrig_demand_ratio(ji,:) / SUM(irrig_demand_ratio(ji,:))
        ENDIF
    ENDDO

    !!! end ratio_irrig, Xuhui

    CALL hydrol_soil(kjpindex, veget_max, soiltile, njsc, reinf_slope,  &
         transpir, vevapnu, vevapnu_pft, evapot, evapot_penm, runoff, &
         drainage, returnflow, reinfiltration, irrigation, irrig_demand_ratio, &
         tot_melt,evap_bare_lim, shumdiag, shumdiag_perma, &
         k_litt, litterhumdiag, humrel, vegstress, drysoil_frac, irrig_fin, &
         is_crop_soil, stempdiag,snow,snowdz, tot_bare_soil, &
         u, v, tq_cdrag, &
         mc_layh, mcl_layh, &
         mc_layh_s, mcl_layh_s, &
!gmjc
         tmc_topgrass, swdown) !!xchen add swdown
!end gmjc

    DO ji = 1,kjpindex
        DO jv = 2,nvm
           IF (.NOT. natural(jv) .AND. ok_LAIdev(jv)) THEN
               IF (.NOT. is_crop_soil(pref_soil_veg(jv))) THEN
                   STOP 'hydrol irrig'
               ENDIF
               !! soil_deficit(ji,jv) = MAX( zero, irrig_fulfill(jv)*tmcs(ji,4) - tmc(ji,4) ) !mm
               ! note that since crop soil may not necessarily be the fourth
               ! soil colum, this needs to be changed !xuhui 151214
               soil_deficit(ji,jv) = MAX( zero, irrig_fulfill(jv)*tmcs(ji,pref_soil_veg(jv)) - tmc(ji,pref_soil_veg(jv)) ) !mm
           ENDIF
        ENDDO
    ENDDO


    ! The update fluxes come from hydrol_vegupd
    drainage(:) =  drainage(:) +  drain_upd(:)
    runoff(:) =  runoff(:) +  runoff_upd(:)


    !! 4 write out file  ==> hydrol_alma/histwrite(*)
    !
    ! If we use the ALMA standards
    CALL hydrol_alma(kjpindex, index, .FALSE., qsintveg, snow, snow_nobio, soilwet)
    

    ! Calculate the moisture in the upper itopmax layers corresponding to 0.1m (humtot_top): 
    ! For ORCHIDEE with nslm=11 and zmaxh=2, itopmax=6. 
    ! We compute tmc_top as tmc but only for the first itopmax layers. Then we compute a humtot with this variable.
    DO jst=1,nstm
       DO ji=1,kjpindex
          tmc_top(ji,jst) = dz(2) * ( trois*mc(ji,1,jst) + mc(ji,2,jst) )/huit
          DO jsl = 2, itopmax
             tmc_top(ji,jst) = tmc_top(ji,jst) + dz(jsl) * (trois*mc(ji,jsl,jst)+mc(ji,jsl-1,jst))/huit &
                  + dz(jsl+1) * (trois*mc(ji,jsl,jst)+mc(ji,jsl+1,jst))/huit
          ENDDO
       ENDDO
    ENDDO
 
    ! We average the values of each soiltile and multiply by vegtot to transform to a grid-cell mean
    humtot_top(:) = zero
    DO jst=1,nstm
       DO ji=1,kjpindex
          humtot_top(ji) = humtot_top(ji) + soiltile(ji,jst) * tmc_top(ji,jst) * vegtot(ji)
       ENDDO
    ENDDO

    ! Calculate the Total Water Budget Residu (in kg/m2 over dt_sechiba)
    ! All the delstocks and fluxes below are averaged over the mesh
    ! snow_nobio included in delswe
    ! Does not include the routing reservoirs, although the flux to/from routing are integrated
    DO ji=1,kjpindex
       twbr(ji) = (delsoilmoist(ji) + delintercept(ji) + delswe(ji)) &
            - ( precip_rain(ji) + precip_snow(ji) + irrigation(ji) + floodout(ji) &
            + returnflow(ji) + reinfiltration(ji) ) &
            + ( runoff(ji) + drainage(ji) + SUM(vevapwet(ji,:)) &
            + SUM(transpir(ji,:)) + vevapnu(ji) + vevapsno(ji) + vevapflo(ji) ) 
    ENDDO
    ! Transform unit from kg/m2/dt to kg/m2/s (or mm/s)
    CALL xios_orchidee_send_field("twbr",twbr/dt_sechiba)
    CALL xios_orchidee_send_field("undermcr",undermcr) ! nb of tiles undermcr at end of timestep

    ! Calculate land_nroot : grid-cell mean of nroot 
    ! Do not treat PFT1 because it has no roots
    land_nroot(:,:) = zero
    DO jsl=1,nslm
       DO jv=2,nvm
          DO ji=1,kjpindex
               IF ( vegtot(ji) > min_sechiba ) THEN 
               land_nroot(ji,jsl) = land_nroot(ji,jsl) + veget_max(ji,jv) * nroot(ji,jv,jsl) / vegtot(ji) 
            END IF
          END DO
       ENDDO
    ENDDO
    CALL xios_orchidee_send_field("nroot",land_nroot)   

    DO jsl=1,nslm
       land_dlh(:,jsl)=dlh(jsl)
    ENDDO
    CALL xios_orchidee_send_field("dlh",land_dlh)

    ! Particular soil moisture values, spatially averaged over the grid-cell
    ! (a) total SM in kg/m2
    !     we average the total values of each soiltile and multiply by vegtot to transform to a grid-cell mean (over total land)
    land_tmcs(:) = zero
    land_tmcfc(:) = zero
    DO jst=1,nstm
       DO ji=1,kjpindex
          land_tmcs(ji) = land_tmcs(ji) + soiltile(ji,jst) * tmcs(ji,jst) * vegtot(ji)
          land_tmcfc(ji) = land_tmcfc(ji) + soiltile(ji,jst) * tmcfc(ji,jst) * vegtot(ji)
       ENDDO
    ENDDO
    CALL xios_orchidee_send_field("tmcs",land_tmcs) ! in kg/m2
    CALL xios_orchidee_send_field("tmcfc",land_tmcfc) ! in kg/m2

    ! (b) volumetric moisture content by layers in m3/m3
    !     mcs etc are identical in all layers (no normalization by vegtot to be comparable to mc)
    DO jsl=1,nslm
       land_mcs(:,jsl) = mcs(:)
       land_mcfc(:,jsl) = mcfc(:)
       land_mcw(:,jsl) = mcw(:)
       land_mcr(:,jsl) = mcr(njsc(:))
    ENDDO
    CALL xios_orchidee_send_field("mcs",land_mcs) ! in m3/m3
    CALL xios_orchidee_send_field("mcfc",land_mcfc) ! in m3/m3 
    CALL xios_orchidee_send_field("mcw",land_mcw) ! in m3/m3 
    CALL xios_orchidee_send_field("mcr",land_mcr) ! in m3/m3 
          
    CALL xios_orchidee_send_field("water2infilt",water2infilt)   

    CALL xios_orchidee_send_field("mc",mc)
    CALL xios_orchidee_send_field("kfact_root",kfact_root)
    CALL xios_orchidee_send_field("vegetmax_soil",vegetmax_soil)

    CALL xios_orchidee_send_field("evapnu_soil",ae_ns/dt_sechiba)
    CALL xios_orchidee_send_field("drainage_soil",dr_ns/dt_sechiba)
    CALL xios_orchidee_send_field("transpir_soil",tr_ns/dt_sechiba)
    CALL xios_orchidee_send_field("runoff_soil",ru_ns/dt_sechiba)
    CALL xios_orchidee_send_field("humrel",humrel)     
    CALL xios_orchidee_send_field("irrig_fin",irrig_fin*one_day/dt_sechiba)
    CALL xios_orchidee_send_field("drainage",drainage/dt_sechiba) ! [kg m-2 s-1]
    CALL xios_orchidee_send_field("runoff",runoff/dt_sechiba) ! [kg m-2 s-1]
    CALL xios_orchidee_send_field("precisol",precisol/dt_sechiba)
    CALL xios_orchidee_send_field("throughfall",throughfall/dt_sechiba)
    CALL xios_orchidee_send_field("precip_rain",precip_rain/dt_sechiba)
    CALL xios_orchidee_send_field("precip_snow",precip_snow/dt_sechiba)
    CALL xios_orchidee_send_field("qsintmax",qsintmax)
    CALL xios_orchidee_send_field("qsintveg",qsintveg)
    CALL xios_orchidee_send_field("qsintveg_tot",SUM(qsintveg(:,:),dim=2))
    histvar(:)=(precip_rain(:)-SUM(throughfall(:,:),dim=2))
    CALL xios_orchidee_send_field("prveg",histvar/dt_sechiba)

    IF ( do_floodplains ) THEN
       CALL xios_orchidee_send_field("floodout",floodout/dt_sechiba)
    END IF

    CALL xios_orchidee_send_field("snowmelt",snowmelt/dt_sechiba)
    CALL xios_orchidee_send_field("tot_melt",tot_melt/dt_sechiba)

    CALL xios_orchidee_send_field("soilmoist",soilmoist)
    CALL xios_orchidee_send_field("soilmoist_liquid",soilmoist_liquid)
    CALL xios_orchidee_send_field("humtot_frozen",SUM(soilmoist(:,:),2)-SUM(soilmoist_liquid(:,:),2))
    CALL xios_orchidee_send_field("tmc",tmc)
    CALL xios_orchidee_send_field("humtot",humtot)
    CALL xios_orchidee_send_field("humtot_top",humtot_top)


    
    ! Prepare diagnostic snow variables depending on snow scheme
    IF (ok_explicitsnow) THEN
       !  Add XIOS default value where no snow
       DO ji=1,kjpindex
          IF (snow(ji) > 0) THEN
             snowdz_diag(ji,:) = snowdz(ji,:)
             snowdepth_diag(ji) = SUM(snowdz(ji,:))*(1-totfrac_nobio(ji))*frac_snow_veg(ji)
          ELSE
             snowdz_diag(ji,:) = xios_default_val
             snowdepth_diag(ji) = xios_default_val             
          END IF
       END DO
       CALL xios_orchidee_send_field("snowdz",snowdz_diag)
       CALL xios_orchidee_send_field("snowdepth",snowdepth_diag)
    ELSE
       ! Add XIOS default value where no snow
       DO ji=1,kjpindex
          IF (snow(ji) > 0) THEN
             snowdz_diag(ji,:) = snowdepth(ji)
             snowdepth_diag(ji) = snowdepth(ji)*(1-totfrac_nobio(ji))
          ELSE
             snowdz_diag(ji,:) = xios_default_val
             snowdepth_diag(ji) = xios_default_val             
          END IF
       END DO
       CALL xios_orchidee_send_field("snowdz",snowdz_diag(:,1))
       CALL xios_orchidee_send_field("snowdepth",snowdepth_diag)
    END IF

    CALL xios_orchidee_send_field("frac_bare",frac_bare)

    CALL xios_orchidee_send_field("soilwet",soilwet)
    CALL xios_orchidee_send_field("delsoilmoist",delsoilmoist)
    CALL xios_orchidee_send_field("delswe",delswe)
    CALL xios_orchidee_send_field("delintercept",delintercept)  

    IF (ok_freeze_cwrr) THEN
       CALL xios_orchidee_send_field("profil_froz_hydro",profil_froz_hydro)
       CALL xios_orchidee_send_field("temp_hydro",temp_hydro)
       CALL xios_orchidee_send_field("profil_froz_hydro_ns", profil_froz_hydro_ns)
    END IF
    
    CALL xios_orchidee_send_field("kk_moy",kk_moy) ! in mm/d


    !! Calculate diagnostic variables on Landuse tiles for LUMIP/CMIP6
    humtot_lut(:,:)=0
    humtot_top_lut(:,:)=0
    mrro_lut(:,:)=0
    DO jv=1,nvm
       jst=pref_soil_veg(jv) ! soil tile index
       IF (natural(jv)) THEN
          humtot_lut(:,id_psl) = humtot_lut(:,id_psl) + tmc(:,jst)*veget_max(:,jv)
          humtot_top_lut(:,id_psl) = humtot_top_lut(:,id_psl) + tmc_top(:,jst)*veget_max(:,jv)
          mrro_lut(:,id_psl) = mrro_lut(:,id_psl) + (dr_ns(:,jst)+ru_ns(:,jst))*veget_max(:,jv)
       ELSE
          humtot_lut(:,id_crp) = humtot_lut(:,id_crp) + tmc(:,jst)*veget_max(:,jv)
          humtot_top_lut(:,id_crp) = humtot_top_lut(:,id_crp) + tmc_top(:,jst)*veget_max(:,jv)
          mrro_lut(:,id_crp) = mrro_lut(:,id_crp) + (dr_ns(:,jst)+ru_ns(:,jst))*veget_max(:,jv)
       ENDIF
    END DO

    WHERE (fraclut(:,id_psl)>min_sechiba)
       humtot_lut(:,id_psl) = humtot_lut(:,id_psl)/fraclut(:,id_psl)
       humtot_top_lut(:,id_psl) = humtot_top_lut(:,id_psl)/fraclut(:,id_psl)
       mrro_lut(:,id_psl) = mrro_lut(:,id_psl)/fraclut(:,id_psl)/dt_sechiba
    ELSEWHERE
       humtot_lut(:,id_psl) = val_exp
       humtot_top_lut(:,id_psl) = val_exp
       mrro_lut(:,id_psl) = val_exp
    END WHERE
    WHERE (fraclut(:,id_crp)>min_sechiba)
       humtot_lut(:,id_crp) = humtot_lut(:,id_crp)/fraclut(:,id_crp)
       humtot_top_lut(:,id_crp) = humtot_top_lut(:,id_crp)/fraclut(:,id_crp)
       mrro_lut(:,id_crp) = mrro_lut(:,id_crp)/fraclut(:,id_crp)/dt_sechiba
    ELSEWHERE
       humtot_lut(:,id_crp) = val_exp
       humtot_top_lut(:,id_crp) = val_exp
       mrro_lut(:,id_crp) = val_exp
    END WHERE

    humtot_lut(:,id_pst) = val_exp
    humtot_lut(:,id_urb) = val_exp
    humtot_top_lut(:,id_pst) = val_exp
    humtot_top_lut(:,id_urb) = val_exp
    mrro_lut(:,id_pst) = val_exp
    mrro_lut(:,id_urb) = val_exp

    CALL xios_orchidee_send_field("humtot_lut",humtot_lut)
    CALL xios_orchidee_send_field("humtot_top_lut",humtot_top_lut)
    CALL xios_orchidee_send_field("mrro_lut",mrro_lut)


    IF ( .NOT. almaoutput ) THEN
       CALL histwrite_p(hist_id, 'frac_bare', kjit, frac_bare, kjpindex*nvm, indexveg)

       DO jst=1,nstm
          ! var_name= "mc_1" ... "mc_3"
          WRITE (var_name,"('moistc_',i1)") jst
          CALL histwrite_p(hist_id, TRIM(var_name), kjit,mc(:,:,jst), kjpindex*nslm, indexlayer)

          ! var_name= "kfactroot_1" ... "kfactroot_3"
          WRITE (var_name,"('kfactroot_',i1)") jst
          CALL histwrite_p(hist_id, TRIM(var_name), kjit, kfact_root(:,:,jst), kjpindex*nslm, indexlayer)

          ! var_name= "vegetsoil_1" ... "vegetsoil_3"
          WRITE (var_name,"('vegetsoil_',i1)") jst
          CALL histwrite_p(hist_id, TRIM(var_name), kjit,vegetmax_soil(:,:,jst), kjpindex*nvm, indexveg)
       ENDDO
       CALL histwrite_p(hist_id, 'precip_soil', kjit, precisol_ns, kjpindex*nstm, indexsoil)
       CALL histwrite_p(hist_id, 'evapnu_soil', kjit, ae_ns, kjpindex*nstm, indexsoil)
       CALL histwrite_p(hist_id, 'drainage_soil', kjit, dr_ns, kjpindex*nstm, indexsoil)
       CALL histwrite_p(hist_id, 'transpir_soil', kjit, tr_ns, kjpindex*nstm, indexsoil)
       CALL histwrite_p(hist_id, 'runoff_soil', kjit, ru_ns, kjpindex*nstm, indexsoil)
       CALL histwrite_p(hist_id, 'humtot_soil', kjit, tmc, kjpindex*nstm, indexsoil)
       ! mrso is a perfect duplicate of humtot
       CALL histwrite_p(hist_id, 'humtot', kjit, humtot, kjpindex, index)
       CALL histwrite_p(hist_id, 'mrso', kjit, humtot, kjpindex, index)
       CALL histwrite_p(hist_id, 'mrsos', kjit, humtot_top, kjpindex, index)
       njsc_tmp(:)=njsc(:)
       CALL histwrite_p(hist_id, 'soilindex', kjit, njsc_tmp, kjpindex, index)
       CALL histwrite_p(hist_id, 'humrel',   kjit, humrel,   kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'vegstress',   kjit, vegstress, kjpindex*nvm, indexveg)
       ! CROP variable
       IF (ANY(ok_LAIdev)) CALL histwrite_p(hist_id, 'soil_deficit', kjit, soil_deficit, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'drainage', kjit, drainage, kjpindex, index)
       ! NB! According to histdef in intersurf, the variables 'runoff' and 'mrros' have different units
       CALL histwrite_p(hist_id, 'runoff', kjit, runoff, kjpindex, index)
       CALL histwrite_p(hist_id, 'mrros', kjit, runoff, kjpindex, index)
       histvar(:)=(runoff(:)+drainage(:))
       CALL histwrite_p(hist_id, 'mrro', kjit, histvar, kjpindex, index)
       CALL histwrite_p(hist_id, 'precisol', kjit, precisol, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'rain', kjit, precip_rain, kjpindex, index)

       histvar(:)=(precip_rain(:)-SUM(throughfall(:,:),dim=2))
       CALL histwrite_p(hist_id, 'prveg', kjit, histvar, kjpindex, index)

       CALL histwrite_p(hist_id, 'snowf', kjit, precip_snow, kjpindex, index)
       CALL histwrite_p(hist_id, 'qsintmax', kjit, qsintmax, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'qsintveg', kjit, qsintveg, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'irrig_fin', kjit, irrig_fin, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'snowmelt',kjit,snowmelt,kjpindex,index)
       CALL histwrite_p(hist_id, 'shumdiag_perma',kjit,shumdiag_perma,kjpindex*nslm,indexnslm)

!pss:+ ! write out wetland fraction and CTI parameters
       IF (TOPM_calcul) CALL topmodel_history(kjit, kjpindex, index, hist_id)
!pss:-

       IF ( river_routing .AND. do_floodplains ) THEN
          CALL histwrite_p(hist_id, 'floodout', kjit, floodout, kjpindex, index)
       ENDIF
       !
       IF ( hist2_id > 0 ) THEN
          DO jst=1,nstm
             ! var_name= "mc_1" ... "mc_3"
             WRITE (var_name,"('moistc_',i1)") jst
             CALL histwrite_p(hist2_id, TRIM(var_name), kjit,mc(:,:,jst), kjpindex*nslm, indexlayer)

             ! var_name= "kfactroot_1" ... "kfactroot_3"
             WRITE (var_name,"('kfactroot_',i1)") jst
             CALL histwrite_p(hist2_id, TRIM(var_name), kjit, kfact_root(:,:,jst), kjpindex*nslm, indexlayer)

             ! var_name= "vegetsoil_1" ... "vegetsoil_3"
             WRITE (var_name,"('vegetsoil_',i1)") jst
             CALL histwrite_p(hist2_id, TRIM(var_name), kjit,vegetmax_soil(:,:,jst), kjpindex*nvm, indexveg)
          ENDDO
          CALL histwrite_p(hist2_id, 'evapnu_soil', kjit, ae_ns, kjpindex*nstm, indexsoil)
          CALL histwrite_p(hist2_id, 'drainage_soil', kjit, dr_ns, kjpindex*nstm, indexsoil)
          CALL histwrite_p(hist2_id, 'transpir_soil', kjit, tr_ns, kjpindex*nstm, indexsoil)
          CALL histwrite_p(hist2_id, 'runoff_soil', kjit, ru_ns, kjpindex*nstm, indexsoil)
          CALL histwrite_p(hist2_id, 'humtot_soil', kjit, tmc, kjpindex*nstm, indexsoil)
          ! mrso is a perfect duplicate of humtot
          CALL histwrite_p(hist2_id, 'humtot', kjit, humtot, kjpindex, index)
          CALL histwrite_p(hist2_id, 'mrso', kjit, humtot, kjpindex, index)
          CALL histwrite_p(hist2_id, 'mrsos', kjit, humtot_top, kjpindex, index)
          njsc_tmp(:)=njsc(:)
          CALL histwrite_p(hist2_id, 'soilindex', kjit, njsc_tmp, kjpindex, index)
          CALL histwrite_p(hist2_id, 'humrel',   kjit, humrel,   kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'drainage', kjit, drainage, kjpindex, index)
          ! NB! According to histdef in intersurf, the variables 'runoff' and 'mrros' have different units
          CALL histwrite_p(hist2_id, 'runoff', kjit, runoff, kjpindex, index)
          CALL histwrite_p(hist2_id, 'mrros', kjit, runoff, kjpindex, index)
          histvar(:)=(runoff(:)+drainage(:))
          CALL histwrite_p(hist2_id, 'mrro', kjit, histvar, kjpindex, index)

          IF ( river_routing .AND. do_floodplains ) THEN
             CALL histwrite_p(hist2_id, 'floodout', kjit, floodout, kjpindex, index)
          ENDIF
          CALL histwrite_p(hist2_id, 'precisol', kjit, precisol, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'rain', kjit, precip_rain, kjpindex, index)
          CALL histwrite_p(hist2_id, 'snowf', kjit, precip_snow, kjpindex, index)
          CALL histwrite_p(hist2_id, 'snowmelt',kjit,snowmelt,kjpindex,index)
          CALL histwrite_p(hist2_id, 'qsintmax', kjit, qsintmax, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'qsintveg', kjit, qsintveg, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'irrig_fin', kjit, irrig_fin, kjpindex*nvm, indexveg)

       ENDIF
    ELSE
       CALL histwrite_p(hist_id, 'Snowf', kjit, precip_snow, kjpindex, index)
       CALL histwrite_p(hist_id, 'Rainf', kjit, precip_rain, kjpindex, index)
       CALL histwrite_p(hist_id, 'Qs', kjit, runoff, kjpindex, index)
       CALL histwrite_p(hist_id, 'Qsb', kjit, drainage, kjpindex, index)
       CALL histwrite_p(hist_id, 'Qsm', kjit, snowmelt, kjpindex, index)
       CALL histwrite_p(hist_id, 'DelSoilMoist', kjit, delsoilmoist, kjpindex, index)
       CALL histwrite_p(hist_id, 'DelSWE', kjit, delswe, kjpindex, index)
       CALL histwrite_p(hist_id, 'DelIntercept', kjit, delintercept, kjpindex, index)
       !
       CALL histwrite_p(hist_id, 'SoilMoist', kjit, soilmoist, kjpindex*nslm, indexlayer)
       CALL histwrite_p(hist_id, 'SoilWet', kjit, soilwet, kjpindex, index)
       !
       CALL histwrite_p(hist_id, 'RootMoist', kjit, tot_watsoil_end, kjpindex, index)
       CALL histwrite_p(hist_id, 'SubSnow', kjit, vevapsno, kjpindex, index)
       !
       IF (.NOT. ok_explicitsnow) CALL histwrite_p(hist_id, 'SnowDepth', kjit, snowdepth, kjpindex, index)
       !
       IF ( hist2_id > 0 ) THEN
          CALL histwrite_p(hist2_id, 'Snowf', kjit, precip_snow, kjpindex, index)
          CALL histwrite_p(hist2_id, 'Rainf', kjit, precip_rain, kjpindex, index)
          CALL histwrite_p(hist2_id, 'Qs', kjit, runoff, kjpindex, index)
          CALL histwrite_p(hist2_id, 'Qsb', kjit, drainage, kjpindex, index)
          CALL histwrite_p(hist2_id, 'Qsm', kjit, snowmelt, kjpindex, index)
          CALL histwrite_p(hist2_id, 'DelSoilMoist', kjit, delsoilmoist, kjpindex, index)
          CALL histwrite_p(hist2_id, 'DelSWE', kjit, delswe, kjpindex, index)
          CALL histwrite_p(hist2_id, 'DelIntercept', kjit, delintercept, kjpindex, index)
          !
          CALL histwrite_p(hist2_id, 'SoilMoist', kjit, soilmoist, kjpindex*nslm, indexlayer)
          CALL histwrite_p(hist2_id, 'SoilWet', kjit, soilwet, kjpindex, index)
          !
          CALL histwrite_p(hist2_id, 'RootMoist', kjit, tot_watsoil_end, kjpindex, index)
          CALL histwrite_p(hist2_id, 'SubSnow', kjit, vevapsno, kjpindex, index)
          !
          IF (.NOT. ok_explicitsnow) CALL histwrite_p(hist2_id, 'SnowDepth', kjit, snowdepth, kjpindex, index)
       ENDIF
    ENDIF

    IF (ok_freeze_cwrr) THEN
       CALL histwrite_p(hist_id, 'profil_froz_hydro', kjit,profil_froz_hydro , kjpindex*nslm, indexlayer)
       
       CALL histwrite_p(hist_id, 'temp_hydro', kjit,temp_hydro , kjpindex*nslm, indexlayer)
    ENDIF
    CALL histwrite_p(hist_id, 'kk_moy', kjit, kk_moy,kjpindex*nslm, indexlayer) ! averaged over soiltiles
    DO jst=1,nstm
       WRITE (var_name,"('profil_froz_hydro_',i1)") jst
       CALL histwrite_p(hist_id, TRIM(var_name), kjit, profil_froz_hydro_ns(:,:,jst), kjpindex*nslm, indexlayer)
    ENDDO

    ! Copy soilmoist into a local variable to be sent to thermosoil
    soilmoist_out(:,:) = soilmoist(:,:)

    IF (printlev>=3) WRITE (numout,*) ' hydrol_main Done '

  END SUBROUTINE hydrol_main


!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_finalize
!!
!>\BRIEF         
!!
!! DESCRIPTION : This subroutine writes the module variables and variables calculated in hydrol to restart file
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE hydrol_finalize( kjit,           kjpindex, rest_id,  vegstress,  &
                              qsintveg,       humrel,                         &
                              snow,           snow_age, snow_nobio,           &
                              snow_nobio_age, snowrho,  snowtemp,             &
                              snowdz,         snowheat,                       &
                              snowgrain, drysoil_frac, evap_bare_lim)
!                              snowcap,        snowgrain, drysoil_frac, evap_bare_lim, evap_bare_lim_pft)

    !! 0. Variable and parameter declaration
    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                           :: kjit           !! Time step number 
    INTEGER(i_std), INTENT(in)                           :: kjpindex       !! Domain size
    INTEGER(i_std),INTENT (in)                           :: rest_id        !! Restart file identifier
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: vegstress      !! Veg. moisture stress (only for vegetation growth)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: qsintveg       !! Water on vegetation due to interception
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: humrel
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: snow           !! Snow mass [Kg/m^2]
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: snow_age       !! Snow age
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (in) :: snow_nobio     !! Water balance on ice, lakes, .. [Kg/m^2]
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (in) :: snow_nobio_age !! Snow age on ice, lakes, ...
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in)  :: snowrho        !! Snow density
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in)  :: snowtemp       !! Snow temperature
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in)  :: snowdz         !! Snow layer thickness
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in)  :: snowheat       !! Snow heat content
    REAL(r_std),DIMENSION (kjpindex),INTENT(in)          :: drysoil_frac   !! function of litter wetness
    REAL(r_std),DIMENSION (kjpindex),INTENT(in)          :: evap_bare_lim
!    REAL(r_std),DIMENSION (kjpindex,nvm),INTENT(in)       :: evap_bare_lim_pft

    REAL(r_std),DIMENSION (kjpindex,nsnow), INTENT(in)   :: snowgrain      !! Snow grain size

    !! 0.4 Local variables
    INTEGER(i_std)                                       :: jst, jsl
   
!_ ================================================================================================================================


    IF (printlev>=3) WRITE (numout,*) ' we have to complete restart file with HYDROLOGIC variables '

    CALL restput_p(rest_id, 'moistc', nbp_glo,  nslm, nstm, kjit, mc, 'scatter',  nbp_glo, index_g)
    CALL restput_p(rest_id, 'moistcl', nbp_glo,  nslm, nstm, kjit, mcl, 'scatter',  nbp_glo, index_g)
     
    CALL restput_p(rest_id, 'us', nbp_glo,nvm, nstm, nslm, kjit,us,'scatter',nbp_glo,index_g)
    
    CALL restput_p(rest_id, 'free_drain_coef', nbp_glo,   nstm, 1, kjit,  free_drain_coef, 'scatter',  nbp_glo, index_g)
    CALL restput_p(rest_id, 'zwt_force', nbp_glo,   nstm, 1, kjit,  zwt_force, 'scatter',  nbp_glo, index_g)
    CALL restput_p(rest_id, 'water2infilt', nbp_glo,   nstm, 1, kjit,  water2infilt, 'scatter',  nbp_glo, index_g)
    CALL restput_p(rest_id, 'ae_ns', nbp_glo,   nstm, 1, kjit,  ae_ns, 'scatter',  nbp_glo, index_g)
    CALL restput_p(rest_id, 'vegstress', nbp_glo,   nvm, 1, kjit,  vegstress, 'scatter',  nbp_glo, index_g)
    CALL restput_p(rest_id, 'snow', nbp_glo,   1, 1, kjit,  snow, 'scatter',  nbp_glo, index_g)
    CALL restput_p(rest_id, 'snow_age', nbp_glo,   1, 1, kjit,  snow_age, 'scatter',  nbp_glo, index_g)
    CALL restput_p(rest_id, 'snow_nobio', nbp_glo,   nnobio, 1, kjit,  snow_nobio, 'scatter', nbp_glo, index_g)
    CALL restput_p(rest_id, 'snow_nobio_age', nbp_glo,   nnobio, 1, kjit,  snow_nobio_age, 'scatter', nbp_glo, index_g)
    CALL restput_p(rest_id, 'qsintveg', nbp_glo, nvm, 1, kjit,  qsintveg, 'scatter',  nbp_glo, index_g)
    CALL restput_p(rest_id, 'evap_bare_lim_ns', nbp_glo, nstm, 1, kjit,  evap_bare_lim_ns, 'scatter',  nbp_glo, index_g)
    CALL restput_p(rest_id, 'evap_bare_lim', nbp_glo, 1, 1, kjit,  evap_bare_lim, 'scatter',  nbp_glo, index_g)
    CALL restput_p(rest_id, 'resdist', nbp_glo, nstm, 1, kjit,  resdist, 'scatter',  nbp_glo, index_g)  
    CALL restput_p(rest_id, 'vegtot_old', nbp_glo, 1, 1, kjit,  vegtot_old, 'scatter',  nbp_glo, index_g)            
    CALL restput_p(rest_id, 'drysoil_frac', nbp_glo,   1, 1, kjit, drysoil_frac, 'scatter', nbp_glo, index_g)
    CALL restput_p(rest_id, 'humrel', nbp_glo,   nvm, 1, kjit,  humrel, 'scatter',  nbp_glo, index_g)
    IF (use_refSOC_hydrol)  CALL restput_p (rest_id, 'refSOC_1d', nbp_glo, 1, 1, kjit, refSOC_1d, 'scatter', nbp_glo, index_g)

    CALL restput_p(rest_id, 'tot_watveg_beg', nbp_glo,  1, 1, kjit,  tot_watveg_beg, 'scatter',  nbp_glo, index_g)
    CALL restput_p(rest_id, 'tot_watsoil_beg', nbp_glo, 1, 1, kjit,  tot_watsoil_beg, 'scatter',  nbp_glo, index_g)
    CALL restput_p(rest_id, 'snow_beg', nbp_glo,        1, 1, kjit,  snow_beg, 'scatter',  nbp_glo, index_g)
    
    ! Write variables for explictsnow module to restart file
    IF (ok_explicitsnow) THEN

      CALL explicitsnow_finalize ( kjit,     kjpindex, rest_id,    snowrho,   &
                                    snowtemp, snowdz,   snowheat,   snowgrain)
    END IF

    !pss:+
    IF (TOPM_calcul) CALL topmodel_writerestart(kjit, kjpindex, nbp_glo, rest_id, index_g)
    !pss:-

  END SUBROUTINE hydrol_finalize


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_init
!!
!>\BRIEF        Initializations and memory allocation   
!!
!! DESCRIPTION  :
!! - 1 Some initializations
!! - 2 make dynamic allocation with good dimension
!! - 2.1 array allocation for soil textur
!! - 2.2 Soil texture choice
!! - 3 Other array allocation
!! - 4 Open restart input file and read data for HYDROLOGIC process
!! - 5 get restart values if none were found in the restart file
!! - 6 Vegetation array      
!! - 7 set humrelv from us
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!!_ hydrol_init

  SUBROUTINE hydrol_init(kjit, kjpindex, index, rest_id, veget_max, soiltile, &
         humrel, vegstress, snow,       snow_age,   snow_nobio, snow_nobio_age, qsintveg, &
         snowdz, snowgrain, snowrho,    snowtemp,   snowheat, &
         drysoil_frac, evap_bare_lim)
!         snowflx,snowcap,   cgrnd_snow, dgrnd_snow, drysoil_frac, evap_bare_lim, evap_bare_lim_pft, &

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT (in)                         :: kjit               !! Time step number 
    INTEGER(i_std), INTENT (in)                         :: kjpindex           !! Domain size
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)    :: index              !! Indeces of the points on the map
    INTEGER(i_std), INTENT (in)                         :: rest_id            !! _Restart_ file identifier 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)   :: veget_max          !! Carte de vegetation max
    REAL(r_std),DIMENSION (kjpindex,nstm), INTENT (in)  :: soiltile           !! Fraction of each soil tile within vegtot (0-1, unitless)

    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)  :: humrel             !! Stress hydrique, relative humidity
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)  :: vegstress          !! Veg. moisture stress (only for vegetation growth)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: snow               !! Snow mass [Kg/m^2]
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: snow_age           !! Snow age
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (out) :: snow_nobio       !! Snow on ice, lakes, ...
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (out) :: snow_nobio_age   !! Snow age on ice, lakes, ...
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)    :: qsintveg         !! Water on vegetation due to interception
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(out)    :: snowdz           !! Snow depth
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(out)    :: snowgrain        !! Snow grain size
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(out)    :: snowheat         !! Snow heat content
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(out)    :: snowtemp         !! Snow temperature
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(out)    :: snowrho          !! Snow density
    REAL(r_std),DIMENSION (kjpindex),INTENT(out)          :: drysoil_frac     !! function of litter wetness
    REAL(r_std),DIMENSION (kjpindex),INTENT(out)          :: evap_bare_lim
!    REAL(r_std),DIMENSION (kjpindex,nvm),INTENT(out)          ::evap_bare_lim_pft

    !! 0.3 Modified variable

    !! 0.4 Local variables

    INTEGER(i_std)                                     :: ier                   !! Error code
    INTEGER(i_std)                                     :: ji                    !! Index of land grid cells (1)
    INTEGER(i_std)                                     :: jv                    !! Index of PFTs (1)
    INTEGER(i_std)                                     :: jst                   !! Index of soil tiles (1)
    INTEGER(i_std)                                     :: jsl                   !! Index of soil layers (1)
    INTEGER(i_std)                                     :: jsc                   !! Index of soil texture (1)
    INTEGER(i_std), PARAMETER                          :: error_level = 3       !! Error level for consistency check
                                                                                !! Switch to 2 tu turn fatal errors into warnings  
    REAL(r_std), ALLOCATABLE, DIMENSION (:)            :: free_drain_max        !! Temporary var for initialization of free_drain_coef 
    REAL(r_std), ALLOCATABLE, DIMENSION (:)            :: zwt_default           !! Temporary variable for initialization of zwt_force
    LOGICAL                                            :: zforce                !! To test if we force the WT in any of the soiltiles

!_ ================================================================================================================================

    !! 1 Some initializations
    !
    !Config Key   = DO_PONDS
    !Config Desc  = Should we include ponds 
    !Config Def   = n
    !Config If    = HYDROL_CWRR
    !Config Help  = This parameters allows the user to ask the model
    !Config         to take into account the ponds and return 
    !Config         the water into the soil moisture. If this is 
    !Config         activated, then there is no reinfiltration 
    !Config         computed inside the hydrol module.
    !Config Units = [FLAG]
    !
    doponds = .FALSE.
    CALL getin_p('DO_PONDS', doponds)

    !Config Key   = FROZ_FRAC_CORR 
    !Config Desc  = Coefficient for the frozen fraction correction
    !Config Def   = 1.0
    !Config If    = HYDROL_CWRR and OK_FREEZE
    !Config Help  =
    !Config Units = [-]
    froz_frac_corr = 1.0
    CALL getin_p("FROZ_FRAC_CORR", froz_frac_corr)

    !Config Key   = MAX_FROZ_HYDRO
    !Config Desc  = Coefficient for the frozen fraction correction
    !Config Def   = 1.0
    !Config If    = HYDROL_CWRR and OK_FREEZE
    !Config Help  =
    !Config Units = [-]
    max_froz_hydro = 1.0
    CALL getin_p("MAX_FROZ_HYDRO", max_froz_hydro)

    !Config Key   = SMTOT_CORR
    !Config Desc  = Coefficient for the frozen fraction correction
    !Config Def   = 2.0
    !Config If    = HYDROL_CWRR and OK_FREEZE
    !Config Help  =
    !Config Units = [-]
    smtot_corr = 2.0
    CALL getin_p("SMTOT_CORR", smtot_corr)

    !Config Key   = DO_RSOIL
    !Config Desc  = Should we reduce soil evaporation with a soil resistance
    !Config Def   = n
    !Config If    = HYDROL_CWRR
    !Config Help  = This parameters allows the user to ask the model
    !Config         to calculate a soil resistance to reduce the soil evaporation
    !Config Units = [FLAG]
    !
    do_rsoil = .FALSE.
    CALL getin_p('DO_RSOIL', do_rsoil) 

    !Config Key   = OK_DYNROOT
    !Config Desc  = Calculate dynamic root profile to optimize soil moisture usage  
    !Config Def   = y
    !Config If    = HYDROL_CWRR
    !Config Help  = 
    !Config Units = [FLAG]
    ok_dynroot = .TRUE.
    CALL getin_p('OK_DYNROOT',ok_dynroot)

    !! 2 make dynamic allocation with good dimension

    !! 2.1 array allocation for soil texture

    ALLOCATE (nvan(nscm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable nvan','','')

    ALLOCATE (avan(nscm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable avan','','')

    ALLOCATE (mcr(nscm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mcr','','')

    ALLOCATE (mcs_mineral(nscm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mcs_mineral','','')

    ALLOCATE (ks(nscm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable ks','','')

    ALLOCATE (pcent(nscm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable pcent','','')

    ALLOCATE (mcfc_mineral(nscm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mcfc_mineral','','')

    ALLOCATE (mcw_mineral(nscm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mcw_mineral','','')
    
    ALLOCATE (mcs(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mcs','','')
    
    ALLOCATE (mcw(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mcw','','')
    
    ALLOCATE (mcfc(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mcfc','','')
    
    ALLOCATE (VG_m(nscm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable VG_m','','')

    ALLOCATE (VG_n(nscm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable VG_n','','')

    ALLOCATE (VG_alpha(nscm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable VG_alpha','','')

    ALLOCATE (VG_psi_fc(nscm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable VG_psi_fc','','')

    ALLOCATE (VG_psi_wp(nscm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable VG_psi_wp','','')

    ALLOCATE (mc_awet(nscm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mc_awet','','')

    ALLOCATE (mc_adry(nscm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mc_adry','','')
       
    !!__2.2 Soil texture choose

    SELECTCASE (nscm)
    CASE (3)
              
       nvan(:) = nvan_fao(:)       
       avan(:) = avan_fao(:)
       mcr(:) = mcr_fao(:)
       mcs_mineral(:) = mcs_fao(:)
       ks(:) = ks_fao(:)
       pcent(:) = pcent_fao(:)
       mcfc_mineral(:) = mcf_fao(:)
       mcw_mineral(:) = mcw_fao(:)
       mc_awet(:) = mc_awet_fao(:)
       mc_adry(:) = mc_adry_fao(:)
    CASE (12)
       
       nvan(:) = nvan_usda(:)
       avan(:) = avan_usda(:)
       mcr(:) = mcr_usda(:)
       mcs_mineral(:) = mcs_usda(:)
       ks(:) = ks_usda(:)
       pcent(:) = pcent_usda(:)
       mcfc_mineral(:) = mcf_usda(:)
       mcw_mineral(:) = mcw_usda(:)
       mc_awet(:) = mc_awet_usda(:)
       mc_adry(:) = mc_adry_usda(:)
       VG_m(:)=VG_m_usda(:)
       VG_n(:)=VG_n_usda(:)
       VG_alpha(:)=VG_alpha_usda(:)
       VG_psi_fc(:)=VG_psi_fc_usda(:)
       VG_psi_wp(:)=VG_psi_wp_usda(:)
       
    CASE DEFAULT
       WRITE (numout,*) 'Unsupported soil type classification. Choose between zobler and usda according to the map'
       CALL ipslerr_p(3,'hydrol_init','Unsupported soil type classification. ',&
            'Choose between zobler and usda according to the map','')
    ENDSELECT

    !! 2.3 Read in the run.def the parameters values defined by the user

    !Config Key   = CWRR_N_VANGENUCHTEN
    !Config Desc  = Van genuchten coefficient n
    !Config If    = HYDROL_CWRR
    !Config Def   = 1.89, 1.56, 1.31
    !Config Help  = This parameter will be constant over the entire 
    !Config         simulated domain, thus independent from soil
    !Config         texture.   
    !Config Units = [-]
    CALL getin_p("CWRR_N_VANGENUCHTEN",nvan)

    !! Check parameter value (correct range)
    IF ( ANY(nvan(:) <= zero) ) THEN
       CALL ipslerr_p(error_level, "hydrol_init.", &
            &     "Wrong parameter value for CWRR_N_VANGENUCHTEN.", &
            &     "This parameter should be positive. ", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !Config Key   = CWRR_A_VANGENUCHTEN
    !Config Desc  = Van genuchten coefficient a
    !Config If    = HYDROL_CWRR
    !Config Def   = 0.0075, 0.0036, 0.0019
    !Config Help  = This parameter will be constant over the entire 
    !Config         simulated domain, thus independent from soil
    !Config         texture.   
    !Config Units = [1/mm]  
    CALL getin_p("CWRR_A_VANGENUCHTEN",avan)

    !! Check parameter value (correct range)
    IF ( ANY(avan(:) <= zero) ) THEN
       CALL ipslerr_p(error_level, "hydrol_init.", &
            &     "Wrong parameter value for CWRR_A_VANGENUCHTEN.", &
            &     "This parameter should be positive. ", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !Config Key   = VWC_RESIDUAL
    !Config Desc  = Residual soil water content
    !Config If    = HYDROL_CWRR
    !Config Def   = 0.065, 0.078, 0.095
    !Config Help  = This parameter will be constant over the entire 
    !Config         simulated domain, thus independent from soil
    !Config         texture.   
    !Config Units = [m3/m3]  
    CALL getin_p("VWC_RESIDUAL",mcr)

    !! Check parameter value (correct range)
    IF ( ANY(mcr(:) < zero) .OR. ANY(mcr(:) > 1.)  ) THEN
       CALL ipslerr_p(error_level, "hydrol_init.", &
            &     "Wrong parameter value for VWC_RESIDUAL.", &
            &     "This parameter is ranged between 0 and 1. ", &
            &     "Please, check parameter value in run.def. ")
    END IF

    
    !Config Key   = VWC_SAT
    !Config Desc  = Saturated soil water content
    !Config If    = HYDROL_CWRR
    !Config Def   = 0.41, 0.43, 0.41
    !Config Help  = This parameter will be constant over the entire 
    !Config         simulated domain, thus independent from soil
    !Config         texture.   
    !Config Units = [m3/m3]  
    CALL getin_p("VWC_SAT",mcs_mineral)

    !! Check parameter value (correct range)
    IF ( ANY(mcs_mineral(:) < zero) .OR. ANY(mcs_mineral(:) > 1.) .OR. ANY(mcs_mineral(:) <= mcr(:)) ) THEN
       CALL ipslerr_p(error_level, "hydrol_init.", &
            &     "Wrong parameter value for VWC_SAT.", &
            &     "This parameter should be greater than VWC_RESIDUAL and less than 1. ", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !Config Key   = CWRR_KS 
    !Config Desc  = Hydraulic conductivity Saturation
    !Config If    = HYDROL_CWRR 
    !Config Def   = 1060.8, 249.6, 62.4
    !Config Help  = This parameter will be constant over the entire 
    !Config         simulated domain, thus independent from soil
    !Config         texture.   
    !Config Units = [mm/d]   
    CALL getin_p("CWRR_KS",ks)

    !! Check parameter value (correct range)
    IF ( ANY(ks(:) <= zero) ) THEN
       CALL ipslerr_p(error_level, "hydrol_init.", &
            &     "Wrong parameter value for CWRR_KS.", &
            &     "This parameter should be positive. ", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !Config Key   = WETNESS_TRANSPIR_MAX
    !Config Desc  = Soil moisture above which transpir is max
    !Config If    = HYDROL_CWRR
    !Config Def   = 0.5, 0.5, 0.5
    !Config Help  = This parameter is independent from soil texture for
    !Config         the time being.
    !Config Units = [-]    
    CALL getin_p("WETNESS_TRANSPIR_MAX",pcent)

    !! Check parameter value (correct range)
    IF ( ANY(pcent(:) <= zero) .OR. ANY(pcent(:) > 1.) ) THEN
       CALL ipslerr_p(error_level, "hydrol_init.", &
            &     "Wrong parameter value for WETNESS_TRANSPIR_MAX.", &
            &     "This parameter should be positive and less or equals than 1. ", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !Config Key   = VWC_FC 
    !Config Desc  = Volumetric water content field capacity
    !Config If    = HYDROL_CWRR
    !Config Def   = 0.32, 0.32, 0.32
    !Config Help  = This parameter is independent from soil texture for
    !Config         the time being.
    !Config Units = [m3/m3]   
    CALL getin_p("VWC_FC",mcfc_mineral)

    !! Check parameter value (correct range)
    IF ( ANY(mcfc_mineral(:) > mcs_mineral(:)) ) THEN
       CALL ipslerr_p(error_level, "hydrol_init.", &
            &     "Wrong parameter value for VWC_FC.", &
            &     "This parameter should be less than VWC_SAT. ", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !Config Key   = VWC_WP
    !Config Desc  = Volumetric water content Wilting pt
    !Config If    = HYDROL_CWRR
    !Config Def   = 0.10, 0.10, 0.10 
    !Config Help  = This parameter is independent from soil texture for
    !Config         the time being.
    !Config Units = [m3/m3]   
    CALL getin_p("VWC_WP",mcw_mineral)

    !! Check parameter value (correct range)
    IF ( ANY(mcw_mineral(:) > mcfc_mineral(:)) .OR. ANY(mcw_mineral(:) < mcr(:)) ) THEN
       CALL ipslerr_p(error_level, "hydrol_init.", &
            &     "Wrong parameter value for VWC_WP.", &
            &     "This parameter should be greater or equal than VWC_RESIDUAL and less or equal than VWC_SAT.", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !Config Key   = VWC_MIN_FOR_WET_ALB
    !Config Desc  = Vol. wat. cont. above which albedo is cst
    !Config If    = HYDROL_CWRR
    !Config Def   = 0.25, 0.25, 0.25
    !Config Help  = This parameter is independent from soil texture for
    !Config         the time being.
    !Config Units = [m3/m3]  
    CALL getin_p("VWC_MIN_FOR_WET_ALB",mc_awet)

    !! Check parameter value (correct range)
    IF ( ANY(mc_awet(:) < 0) ) THEN
       CALL ipslerr_p(error_level, "hydrol_init.", &
            &     "Wrong parameter value for VWC_MIN_FOR_WET_ALB.", &
            &     "This parameter should be positive. ", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !Config Key   = VWC_MAX_FOR_DRY_ALB
    !Config Desc  = Vol. wat. cont. below which albedo is cst
    !Config If    = HYDROL_CWRR
    !Config Def   = 0.1, 0.1, 0.1
    !Config Help  = This parameter is independent from soil texture for
    !Config         the time being.
    !Config Units = [m3/m3]   
    CALL getin_p("VWC_MAX_FOR_DRY_ALB",mc_adry)

    !! Check parameter value (correct range)
    IF ( ANY(mc_adry(:) < 0) .OR. ANY(mc_adry(:) > mc_awet(:)) ) THEN
       CALL ipslerr_p(error_level, "hydrol_init.", &
            &     "Wrong parameter value for VWC_MAX_FOR_DRY_ALB.", &
            &     "This parameter should be positive and not greater than VWC_MIN_FOR_WET_ALB.", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !! 3 Other array allocation


    ALLOCATE (mask_veget(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mask_veget','','')

    ALLOCATE (irrig_fin(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable irrig_fin','','')

    ALLOCATE (mask_soiltile(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mask_soiltile','','')

    ALLOCATE (humrelv(kjpindex,nvm,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable humrelv','','')

    ALLOCATE (vegstressv(kjpindex,nvm,nstm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable vegstressv','','')

    ALLOCATE (us(kjpindex,nvm,nstm,nslm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable us','','')

    ALLOCATE (precisol(kjpindex,nvm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable precisol','','')

    ALLOCATE (throughfall(kjpindex,nvm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable throughfall','','')

    ALLOCATE (precisol_ns(kjpindex,nstm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable precisol_nc','','')

    ALLOCATE (free_drain_coef(kjpindex,nstm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable free_drain_coef','','')

    ALLOCATE (zwt_force(kjpindex,nstm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable zwt_force','','')

    ALLOCATE (frac_bare_ns(kjpindex,nstm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable frac_bare_ns','','')

    ALLOCATE (water2infilt(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable water2infilt','','')

    ALLOCATE (ae_ns(kjpindex,nstm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable ae_ns','','')

    ALLOCATE (evap_bare_lim_ns(kjpindex,nstm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable evap_bare_lim_ns','','')

    ALLOCATE (rootsink(kjpindex,nslm,nstm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable rootsink','','')

    ALLOCATE (subsnowveg(kjpindex),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable subsnowveg','','')

    ALLOCATE (subsnownobio(kjpindex,nnobio),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable subsnownobio','','')

    ALLOCATE (snowmelt(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in snowmelt allocation. We stop. We need kjpindex words = ',kjpindex
        STOP 'hydrol_init'
    END IF

    ALLOCATE (icemelt(kjpindex),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable icemelt','','')

    ALLOCATE (subsinksoil(kjpindex),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable subsinksoil','','')

    ALLOCATE (mx_eau_var(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mx_eau_var','','')

    ALLOCATE (vegtot(kjpindex),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable vegtot','','')

    ALLOCATE (vegtot_old(kjpindex),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable vegtot_old','','')

    ALLOCATE (resdist(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable resdist','','')

    ALLOCATE (humtot(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable humtot','','')

    ALLOCATE (resolv(kjpindex),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable resolv','','')

    ALLOCATE (k(kjpindex,nslm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable k','','')

    ALLOCATE (kk_moy(kjpindex,nslm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable kk_moy','','')
    kk_moy(:,:) = 276.48
    
    ALLOCATE (kk(kjpindex,nslm,nstm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable kk','','')
    kk(:,:,:) = 276.48
    
    ALLOCATE (avan_mod_tab(nslm,nscm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable avan_mod_tab','','')
    
    ALLOCATE (nvan_mod_tab(nslm,nscm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable nvan_mod_tab','','')

    ALLOCATE (a(kjpindex,nslm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable a','','')

    ALLOCATE (b(kjpindex,nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable b','','')

    ALLOCATE (d(kjpindex,nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable d','','')

    ALLOCATE (e(kjpindex,nslm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable e','','')

    ALLOCATE (f(kjpindex,nslm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable f','','')

    ALLOCATE (g1(kjpindex,nslm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable g1','','')

    ALLOCATE (ep(kjpindex,nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable ep','','')

    ALLOCATE (fp(kjpindex,nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable fp','','')

    ALLOCATE (gp(kjpindex,nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable gp','','')

    ALLOCATE (rhs(kjpindex,nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable rhs','','')

    ALLOCATE (srhs(kjpindex,nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable srhs','','')

    ALLOCATE (tmc(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tmc','','')

    ALLOCATE (tmcs(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tmcs','','')

    ALLOCATE (tmcr(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tmcr','','')

    ALLOCATE (tmcfc(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tmcfc','','')

    ALLOCATE (tmcw(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tmcw','','')

    ALLOCATE (tmc_litter(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tmc_litter','','')
!gmjc top 5 layer mc for grazing
    ALLOCATE (tmc_trampling(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tmc_trampling','','')
!end gmjc
    ALLOCATE (tmc_litt_mea(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tmc_litt_mea','','')

    ALLOCATE (tmc_litter_res(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tmc_litter_res','','')

    ALLOCATE (tmc_litter_wilt(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tmc_litter_wilt','','')

    ALLOCATE (tmc_litter_field(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tmc_litter_field','','')

    ALLOCATE (tmc_litter_sat(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tmc_litter_sat','','')

    ALLOCATE (tmc_litter_awet(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tmc_litter_awet','','')

    ALLOCATE (tmc_litter_adry(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tmc_litter_adry','','')

    ALLOCATE (tmc_litt_wet_mea(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tmc_litt_wet_mea','','')

    ALLOCATE (tmc_litt_dry_mea(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tmc_litt_dry_mea','','')

    ALLOCATE (v1(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable v1','','')

    ALLOCATE (ru_ns(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable ru_ns','','')
    ru_ns(:,:) = zero

    ALLOCATE (dr_ns(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable dr_ns','','')
    dr_ns(:,:) = zero

    ALLOCATE (tr_ns(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tr_ns','','')

    ALLOCATE (vegetmax_soil(kjpindex,nvm,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable vegetmax_soil','','')

    ALLOCATE (mc(kjpindex,nslm,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mc','','')


    ! Variables for nudging of soil moisture
    IF (ok_nudge_mc) THEN
       ALLOCATE (mc_read_prev(kjpindex,nslm,nstm),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mc_read_prev','','')
       ALLOCATE (mc_read_next(kjpindex,nslm,nstm),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mc_read_next','','')
       ALLOCATE (mask_mc_interp(kjpindex,nslm,nstm),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mask_mc_interp','','')
    END IF

    ! Variables for nudging of snow variables
    IF (ok_nudge_snow) THEN
       ALLOCATE (snowdz_read_prev(kjpindex,nsnow),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable snowdz_read_prev','','')
       ALLOCATE (snowdz_read_next(kjpindex,nsnow),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable snowdz_read_next','','')
       
       ALLOCATE (snowrho_read_prev(kjpindex,nsnow),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable snowrho_read_prev','','')
       ALLOCATE (snowrho_read_next(kjpindex,nsnow),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable snowrho_read_next','','')
       
       ALLOCATE (snowtemp_read_prev(kjpindex,nsnow),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable snowtemp_read_prev','','')
       ALLOCATE (snowtemp_read_next(kjpindex,nsnow),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable snowtemp_read_next','','')
       
       ALLOCATE (mask_snow_interp(kjpindex,nsnow),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mask_snow_interp','','')
    END IF

    
    ALLOCATE (profil_froz_hydro_ns(kjpindex, nslm, nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable profil_froz_hydro_ns','','')
    profil_froz_hydro_ns(:,:,:) = zero
    
    ALLOCATE (soilmoist(kjpindex,nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable soilmoist','','')

    ALLOCATE (soilmoist_liquid(kjpindex,nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable soilmoist_liquid','','')

    ALLOCATE (soil_wet_ns(kjpindex,nslm,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable soil_wet_ns','','')

    ALLOCATE (soil_wet_litter(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable soil_wet_litter','','')

    ALLOCATE (qflux(kjpindex,nslm,nstm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable qflux','','')

    ALLOCATE (tmat(kjpindex,nslm,3),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tmat','','')

    ALLOCATE (stmat(kjpindex,nslm,3),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable stmat','','')

    ALLOCATE (nroot(kjpindex,nvm, nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable nroot','','')

    ALLOCATE (kfact_root(kjpindex, nslm, nstm), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable kfact_root','','')

    ALLOCATE (kfact(nslm, nscm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable kfact','','')

    ALLOCATE (zz(nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable zz','','')

    ALLOCATE (dz(nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable dz','','')
    
    ALLOCATE (dh(nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable dh','','')

    ALLOCATE (mc_lin(kjpindex, imin:imax),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mc_lin','','')

    ALLOCATE (k_lin(kjpindex, nslm, imin:imax),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable k_lin','','')

    ALLOCATE (d_lin(kjpindex, nslm, imin:imax),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable d_lin','','')

    ALLOCATE (a_lin(kjpindex, nslm, imin:imax),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable a_lin','','')

    ALLOCATE (b_lin(kjpindex, nslm, imin:imax),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable b_lin','','')

!!! xuhui: this has to be defined even if not ok_freeze_cwrr

    IF (ok_freeze_cwrr) THEN
       ALLOCATE (profil_froz_hydro(kjpindex, nslm),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable profil_froz_hydrol','','')
       profil_froz_hydro(:,:) = zero

       ALLOCATE (temp_hydro(kjpindex, nslm),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable temp_hydro','','')
       temp_hydro(:,:) = 280.
    ENDIF

    ALLOCATE (mcl(kjpindex, nslm, nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mcl','','')

    ALLOCATE (undermcr(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable undermcr','','')

    ALLOCATE (tot_watveg_beg(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tot_watveg_beg','','')
    
    ALLOCATE (tot_watveg_end(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tot_watvag_end','','')
    
    ALLOCATE (tot_watsoil_beg(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tot_watsoil_beg','','')
    
    ALLOCATE (tot_watsoil_end(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tot_watsoil_end','','')
    
    ALLOCATE (delsoilmoist(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable delsoilmoist','','')
    
    ALLOCATE (delintercept(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable delintercept','','')
    
    ALLOCATE (delswe(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable delswe','','')
    
    ALLOCATE (snow_beg(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable snow_beg','','')
    
    ALLOCATE (snow_end(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable snow_end','','')
    
    !! 4 Open restart input file and read data for HYDROLOGIC process
       IF (printlev>=3) WRITE (numout,*) ' we have to read a restart file for HYDROLOGIC variables'

       CALL ioconf_setatt_p('UNITS', '-')
       CALL restget_p (rest_id, 'moistc', nbp_glo, nslm , nstm, kjit, .TRUE., mc, "gather", nbp_glo, index_g)
       !
       CALL ioconf_setatt_p('UNITS', '-')
       CALL restget_p (rest_id, 'moistcl', nbp_glo, nslm , nstm, kjit, .TRUE., mcl, "gather", nbp_glo, index_g)

       CALL ioconf_setatt_p('UNITS', '-')
       CALL ioconf_setatt_p('LONG_NAME','us')
       CALL restget_p (rest_id, 'us', nbp_glo, nvm, nstm, nslm, kjit, .TRUE., us, "gather", nbp_glo, index_g)
       !
       var_name= 'free_drain_coef'
       CALL ioconf_setatt_p('UNITS', '-')
       CALL ioconf_setatt_p('LONG_NAME','Coefficient for free drainage at bottom of soil')
       CALL restget_p (rest_id, var_name, nbp_glo, nstm, 1, kjit, .TRUE., free_drain_coef, "gather", nbp_glo, index_g)
       !
       var_name= 'zwt_force'
       CALL ioconf_setatt_p('UNITS', 'm')
       CALL ioconf_setatt_p('LONG_NAME','Prescribed water table depth')
       CALL restget_p (rest_id, var_name, nbp_glo, nstm, 1, kjit, .TRUE., zwt_force, "gather", nbp_glo, index_g)
       !
       var_name= 'water2infilt'
       CALL ioconf_setatt_p('UNITS', '-')
       CALL ioconf_setatt_p('LONG_NAME','Remaining water to be infiltrated on top of the soil')
       CALL restget_p (rest_id, var_name, nbp_glo, nstm, 1, kjit, .TRUE., water2infilt, "gather", nbp_glo, index_g)
       !
       var_name= 'ae_ns'
       CALL ioconf_setatt_p('UNITS', 'kg/m^2')
       CALL ioconf_setatt_p('LONG_NAME','Bare soil evap on each soil type')
       CALL restget_p (rest_id, var_name, nbp_glo, nstm, 1, kjit, .TRUE., ae_ns, "gather", nbp_glo, index_g)
       !
       var_name= 'snow'        
       CALL ioconf_setatt_p('UNITS', 'kg/m^2')
       CALL ioconf_setatt_p('LONG_NAME','Snow mass')
       CALL restget_p (rest_id, var_name, nbp_glo, 1  , 1, kjit, .TRUE., snow, "gather", nbp_glo, index_g)
       !
       var_name= 'snow_age'
          CALL ioconf_setatt_p('UNITS', 'd')
          CALL ioconf_setatt_p('LONG_NAME','Snow age')
       CALL restget_p (rest_id, var_name, nbp_glo, 1  , 1, kjit, .TRUE., snow_age, "gather", nbp_glo, index_g)
       !
       var_name= 'snow_nobio'
          CALL ioconf_setatt_p('UNITS', 'kg/m^2')
          CALL ioconf_setatt_p('LONG_NAME','Snow on other surface types')
       CALL restget_p (rest_id, var_name, nbp_glo, nnobio  , 1, kjit, .TRUE., snow_nobio, "gather", nbp_glo, index_g)
       !
       var_name= 'snow_nobio_age'
          CALL ioconf_setatt_p('UNITS', 'd')
          CALL ioconf_setatt_p('LONG_NAME','Snow age on other surface types')
       CALL restget_p (rest_id, var_name, nbp_glo, nnobio  , 1, kjit, .TRUE., snow_nobio_age, "gather", nbp_glo, index_g)
       !
       var_name= 'qsintveg'
          CALL ioconf_setatt_p('UNITS', 'kg/m^2')
          CALL ioconf_setatt_p('LONG_NAME','Intercepted moisture')
       CALL restget_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, .TRUE., qsintveg, "gather", nbp_glo, index_g)

       var_name= 'evap_bare_lim_ns'
          CALL ioconf_setatt_p('UNITS', '?')
          CALL ioconf_setatt_p('LONG_NAME','?')
       CALL restget_p (rest_id, var_name, nbp_glo, nstm, 1, kjit, .TRUE., evap_bare_lim_ns, "gather", nbp_glo, index_g)
       CALL setvar_p (evap_bare_lim_ns, val_exp, 'NO_KEYWORD', 0.0)
!       DO jv = 1,nvm
!           evap_bare_lim_pft(:,jv) = evap_bare_lim_ns(:,pref_soil_veg(jv))
!       ENDDO

       var_name= 'resdist'
          CALL ioconf_setatt_p('UNITS', '-')
          CALL ioconf_setatt_p('LONG_NAME','soiltile values from previous time-step')
       CALL restget_p (rest_id, var_name, nbp_glo, nstm, 1, kjit, .TRUE., resdist, "gather", nbp_glo, index_g)

       var_name= 'vegtot_old'
          CALL ioconf_setatt_p('UNITS', '-')
          CALL ioconf_setatt_p('LONG_NAME','vegtot from previous time-step')
       CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., vegtot_old, "gather", nbp_glo, index_g)       
       
       ! Read drysoil_frac. It will be initalized later in hydrol_var_init if the varaible is not find in restart file.
          CALL ioconf_setatt_p('UNITS', '')
          CALL ioconf_setatt_p('LONG_NAME','Function of litter wetness')
       CALL restget_p (rest_id, 'drysoil_frac', nbp_glo, 1  , 1, kjit, .TRUE., drysoil_frac, "gather", nbp_glo, index_g)


    !! 5 get restart values if none were found in the restart file
       !
       !Config Key   = HYDROL_MOISTURE_CONTENT
       !Config Desc  = Soil moisture on each soil tile and levels
       !Config If    = HYDROL_CWRR       
       !Config Def   = 0.3
       !Config Help  = The initial value of mc if its value is not found
       !Config         in the restart file. This should only be used if the model is 
       !Config         started without a restart file.
       !Config Units = [m3/m3]
       !
       CALL setvar_p (mc, val_exp, 'HYDROL_MOISTURE_CONTENT', 0.3_r_std)

       ! Initialize mcl as mc if it is not found in the restart file
       IF ( ALL(mcl(:,:,:)==val_exp) ) THEN
          mcl(:,:,:) = mc(:,:,:)
       END IF

       
       !Config Key   = US_INIT
       !Config Desc  = US_NVM_NSTM_NSLM
       !Config If    = HYDROL_CWRR       
       !Config Def   = 0.0
       !Config Help  = The initial value of us (relative moisture) if its value is not found
       !Config         in the restart file. This should only be used if the model is 
       !Config         started without a restart file.
       !Config Units = [-]
       !
       DO jsl=1,nslm
          CALL setvar_p (us(:,:,:,jsl), val_exp, 'US_INIT', zero)
       ENDDO
       !
       !Config Key   = ZWT_FORCE
       !Config Desc  = Prescribed water depth, dimension nstm
       !Config If    = HYDROL_CWRR       
       !Config Def   = undef undef undef
       !Config Help  = The initial value of zwt_force if its value is not found
       !Config         in the restart file. undef corresponds to a case whith no forced WT. 
       !Config         This should only be used if the model is started without a restart file.
       !Config Units = [m]
       
       ALLOCATE (zwt_default(nstm),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable zwt_default','','')
       zwt_default(:) = undef_sechiba
       CALL setvar_p (zwt_force, val_exp, 'ZWT_FORCE', zwt_default )

       zforce = .FALSE.
       DO jst=1,nstm
          IF (zwt_force(1,jst) <= zmaxh) zforce = .TRUE. ! AD16*** check if OK with vertical_soil
       ENDDO
       !
       !Config Key   = FREE_DRAIN_COEF
       !Config Desc  = Coefficient for free drainage at bottom, dimension nstm
       !Config If    = HYDROL_CWRR       
       !Config Def   = 1.0 1.0 1.0
       !Config Help  = The initial value of free drainage coefficient if its value is not found
       !Config         in the restart file. This should only be used if the model is 
       !Config         started without a restart file.
       !Config Units = [-]
              
       ALLOCATE (free_drain_max(nstm),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable free_drain_max','','')
       free_drain_max(:)=1.0

       CALL setvar_p (free_drain_coef, val_exp, 'FREE_DRAIN_COEF', free_drain_max)
       DEALLOCATE(free_drain_max)

       !
       !Config Key   = WATER_TO_INFILT
       !Config Desc  = Water to be infiltrated on top of the soil
       !Config If    = HYDROL_CWRR    
       !Config Def   = 0.0
       !Config Help  = The initial value of free drainage if its value is not found
       !Config         in the restart file. This should only be used if the model is 
       !Config         started without a restart file.
       !Config Units = [mm]
       !
       CALL setvar_p (water2infilt, val_exp, 'WATER_TO_INFILT', zero)
       !
       !Config Key   = EVAPNU_SOIL
       !Config Desc  = Bare soil evap on each soil if not found in restart
       !Config If    = HYDROL_CWRR  
       !Config Def   = 0.0
       !Config Help  = The initial value of bare soils evap if its value is not found
       !Config         in the restart file. This should only be used if the model is 
       !Config         started without a restart file.
       !Config Units = [mm]
       !
       CALL setvar_p (ae_ns, val_exp, 'EVAPNU_SOIL', zero)
       !
       !Config Key  = HYDROL_SNOW
       !Config Desc  = Initial snow mass if not found in restart
       !Config If    = OK_SECHIBA
       !Config Def   = 0.0
       !Config Help  = The initial value of snow mass if its value is not found
       !Config         in the restart file. This should only be used if the model is 
       !Config         started without a restart file.
       !Config Units =
       !
       CALL setvar_p (snow, val_exp, 'HYDROL_SNOW', zero)
       !
       !Config Key   = HYDROL_SNOWAGE
       !Config Desc  = Initial snow age if not found in restart
       !Config If    = OK_SECHIBA
       !Config Def   = 0.0
       !Config Help  = The initial value of snow age if its value is not found
       !Config         in the restart file. This should only be used if the model is 
       !Config         started without a restart file.
       !Config Units = ***
       !
       CALL setvar_p (snow_age, val_exp, 'HYDROL_SNOWAGE', zero)
       !
       !Config Key   = HYDROL_SNOW_NOBIO
       !Config Desc  = Initial snow amount on ice, lakes, etc. if not found in restart
       !Config If    = OK_SECHIBA
       !Config Def   = 0.0
       !Config Help  = The initial value of snow if its value is not found
       !Config         in the restart file. This should only be used if the model is 
       !Config         started without a restart file.
       !Config Units = [mm]
       !
       CALL setvar_p (snow_nobio, val_exp, 'HYDROL_SNOW_NOBIO', zero)
       !
       !Config Key   = HYDROL_SNOW_NOBIO_AGE
       !Config Desc  = Initial snow age on ice, lakes, etc. if not found in restart
       !Config If    = OK_SECHIBA
       !Config Def   = 0.0
       !Config Help  = The initial value of snow age if its value is not found
       !Config         in the restart file. This should only be used if the model is 
       !Config         started without a restart file.
       !Config Units = ***
       !
       CALL setvar_p (snow_nobio_age, val_exp, 'HYDROL_SNOW_NOBIO_AGE', zero)
       !
       !Config Key   = HYDROL_QSV
       !Config Desc  = Initial water on canopy if not found in restart
       !Config If    = OK_SECHIBA
       !Config Def   = 0.0
       !Config Help  = The initial value of moisture on canopy if its value 
       !Config         is not found in the restart file. This should only be used if
       !Config         the model is started without a restart file. 
       !Config Units = [mm]
       !
       CALL setvar_p (qsintveg, val_exp, 'HYDROL_QSV', zero)

       IF (ok_freeze_cwrr) THEN  
          CALL setvar_p (profil_froz_hydro, val_exp, 'NO_KEYWORD', zero)
          CALL setvar_p (profil_froz_hydro_ns, val_exp, 'NO_KEYWORD', zero)
          CALL setvar_p (kk, val_exp, 'NO_KEYWORD', 276.48)
          CALL setvar_p (kk_moy, val_exp, 'NO_KEYWORD', 276.48)
          CALL setvar_p (temp_hydro, val_exp, 'NO_KEYWORD', 280.)
       ENDIF
       
    !! 6 Vegetation array      
       !
       ! If resdist is not in restart file, initialize with soiltile
       IF ( MINVAL(resdist) .EQ.  MAXVAL(resdist) .AND. MINVAL(resdist) .EQ. val_exp) THEN
          resdist(:,:) = soiltile(:,:)
       ENDIF
       
       !
       !  Remember that it is only frac_nobio + SUM(veget_max(,:)) that is equal to 1. Thus we need vegtot
       !
       IF ( ALL(vegtot_old(:) == val_exp) ) THEN
          ! vegtot_old was not found in restart file
          DO ji = 1, kjpindex
             vegtot_old(ji) = SUM(veget_max(ji,:))
          ENDDO
       ENDIF
       
       ! In the initialization phase, vegtot must take the value from previous time-step. 
       ! This is because hydrol_main is done before veget_max is updated in the end of the time step. 
       vegtot(:) = vegtot_old(:)
       
       !
       !
       ! compute the masks for veget

       mask_veget(:,:) = 0
       mask_soiltile(:,:) = 0

       DO jst=1,nstm
          DO ji = 1, kjpindex
             IF(soiltile(ji,jst) .GT. min_sechiba) THEN
                mask_soiltile(ji,jst) = 1
             ENDIF
          END DO
       ENDDO
          
       DO jv = 1, nvm
          DO ji = 1, kjpindex
             IF(veget_max(ji,jv) .GT. min_sechiba) THEN
                mask_veget(ji,jv) = 1
             ENDIF
          END DO
       END DO

       humrelv(:,:,:) = SUM(us,dim=4)

          
       !! 7a. Set vegstress
     
       var_name= 'vegstress'
          CALL ioconf_setatt_p('UNITS', '-')
          CALL ioconf_setatt_p('LONG_NAME','Vegetation growth moisture stress')
       CALL restget_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, .TRUE., vegstress, "gather", nbp_glo, index_g)

       vegstressv(:,:,:) = humrelv(:,:,:)
       ! Calculate vegstress if it is not found in restart file
       IF (ALL(vegstress(:,:)==val_exp)) THEN
          DO jv=1,nvm
             DO ji=1,kjpindex
                vegstress(ji,jv)=vegstress(ji,jv) + vegstressv(ji,jv,pref_soil_veg(jv))
             END DO
          END DO
       END IF
       !! 7b. Set humrel   
       ! Read humrel from restart file
       var_name= 'humrel'
          CALL ioconf_setatt_p('UNITS', '')
          CALL ioconf_setatt_p('LONG_NAME','Relative humidity')
       CALL restget_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, .TRUE., humrel, "gather", nbp_glo, index_g)

       ! Calculate humrel if it is not found in restart file
       IF (ALL(humrel(:,:)==val_exp)) THEN
          ! set humrel from humrelv, assuming equi-repartition for the first time step
          humrel(:,:) = zero
          DO jv=1,nvm
             DO ji=1,kjpindex
                humrel(ji,jv)=humrel(ji,jv) + humrelv(ji,jv,pref_soil_veg(jv))      
             END DO
          END DO
       END IF

       ! Read evap_bare_lim from restart file
       var_name= 'evap_bare_lim'
          CALL ioconf_setatt_p('UNITS', '')
          CALL ioconf_setatt_p('LONG_NAME','Limitation factor for bare soil evaporation')
       CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., evap_bare_lim, "gather", nbp_glo, index_g)

       ! Calculate evap_bare_lim if it was not found in the restart file.
       IF ( ALL(evap_bare_lim(:) == val_exp) ) THEN
          DO ji = 1, kjpindex
             evap_bare_lim(ji) =  SUM(evap_bare_lim_ns(ji,:)*vegtot(ji)*soiltile(ji,:))
          ENDDO
       END IF


    ! Read from restart file       
    ! The variables tot_watsoil_beg, tot_watsoil_beg and snwo_beg will be initialized in the end of 
    ! hydrol_initialize if they were not found in the restart file.
       
    var_name= 'tot_watveg_beg'
       CALL ioconf_setatt_p('UNITS', '?')
       CALL ioconf_setatt_p('LONG_NAME','?')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., tot_watveg_beg, "gather", nbp_glo, index_g)
    
    var_name= 'tot_watsoil_beg'
       CALL ioconf_setatt_p('UNITS', '?')
       CALL ioconf_setatt_p('LONG_NAME','?')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., tot_watsoil_beg, "gather", nbp_glo, index_g)
    
    var_name= 'snow_beg'
       CALL ioconf_setatt_p('UNITS', '?')
       CALL ioconf_setatt_p('LONG_NAME','?')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., snow_beg, "gather", nbp_glo, index_g)
       
 
    ! Initialize variables for explictsnow module by reading restart file
    IF (ok_explicitsnow) THEN
       CALL explicitsnow_initialize( kjit,     kjpindex, rest_id,    snowrho,   &
                                     snowtemp, snowdz,   snowheat,   snowgrain)
    END IF

    
    IF (printlev>=3) WRITE (numout,*) ' hydrol_init done '
    
  END SUBROUTINE hydrol_init


!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_clear
!!
!>\BRIEF        Deallocate arrays 
!!
!_ ================================================================================================================================
!_ hydrol_clear

  SUBROUTINE hydrol_clear()

    ! Allocation for soiltile related parameters
    IF ( ALLOCATED (nvan)) DEALLOCATE (nvan)
    IF ( ALLOCATED (avan)) DEALLOCATE (avan)
    IF ( ALLOCATED (mcr)) DEALLOCATE (mcr)
    IF ( ALLOCATED (mcs_mineral)) DEALLOCATE (mcs_mineral)
    IF ( ALLOCATED (ks)) DEALLOCATE (ks)
    IF ( ALLOCATED (pcent)) DEALLOCATE (pcent)
    IF ( ALLOCATED (mcfc_mineral)) DEALLOCATE (mcfc_mineral)
    IF ( ALLOCATED (mcw_mineral)) DEALLOCATE (mcw_mineral)
    IF ( ALLOCATED (mcs)) DEALLOCATE (mcs)
    IF ( ALLOCATED (mcfc)) DEALLOCATE (mcfc)
    IF ( ALLOCATED (mcw)) DEALLOCATE (mcw)
    IF ( ALLOCATED (VG_m)) DEALLOCATE (VG_m)
    IF ( ALLOCATED (VG_n)) DEALLOCATE (VG_n)
    IF ( ALLOCATED (VG_alpha)) DEALLOCATE (VG_alpha)
    IF ( ALLOCATED (VG_psi_fc)) DEALLOCATE (VG_psi_fc)
    IF ( ALLOCATED (VG_psi_wp)) DEALLOCATE (VG_psi_wp)
    IF ( ALLOCATED (mc_awet)) DEALLOCATE (mc_awet)
    IF ( ALLOCATED (mc_adry)) DEALLOCATE (mc_adry)
!!pss+
!    IF ( ALLOCATED (mcs_grid)) DEALLOCATE (mcs_grid)
!    IF ( ALLOCATED (mcw_grid)) DEALLOCATE (mcw_grid)
!!pss-
    ! Other arrays
    IF (ALLOCATED (mask_veget)) DEALLOCATE (mask_veget)
    IF (ALLOCATED (irrig_fin)) DEALLOCATE (irrig_fin)
    IF (ALLOCATED (mask_soiltile)) DEALLOCATE (mask_soiltile)
    IF (ALLOCATED (humrelv)) DEALLOCATE (humrelv)
    IF (ALLOCATED (vegstressv)) DEALLOCATE (vegstressv)
    IF (ALLOCATED (us)) DEALLOCATE (us)
    IF (ALLOCATED  (precisol)) DEALLOCATE (precisol)
    IF (ALLOCATED  (throughfall)) DEALLOCATE (throughfall)
    IF (ALLOCATED  (precisol_ns)) DEALLOCATE (precisol_ns)
    IF (ALLOCATED  (free_drain_coef)) DEALLOCATE (free_drain_coef)
    IF (ALLOCATED  (frac_bare_ns)) DEALLOCATE (frac_bare_ns)
    IF (ALLOCATED  (water2infilt)) DEALLOCATE (water2infilt)
    IF (ALLOCATED  (ae_ns)) DEALLOCATE (ae_ns)
    IF (ALLOCATED  (evap_bare_lim_ns)) DEALLOCATE (evap_bare_lim_ns)
    IF (ALLOCATED  (rootsink)) DEALLOCATE (rootsink)
    IF (ALLOCATED  (subsnowveg)) DEALLOCATE (subsnowveg)
    IF (ALLOCATED  (subsnownobio)) DEALLOCATE (subsnownobio)
    IF (ALLOCATED  (snowmelt)) DEALLOCATE (snowmelt)
    IF (ALLOCATED  (icemelt)) DEALLOCATE (icemelt)
    IF (ALLOCATED  (subsinksoil)) DEALLOCATE (subsinksoil)
    IF (ALLOCATED  (mx_eau_var)) DEALLOCATE (mx_eau_var)
    IF (ALLOCATED  (vegtot)) DEALLOCATE (vegtot)
    IF (ALLOCATED  (vegtot_old)) DEALLOCATE (vegtot_old)
    IF (ALLOCATED  (resdist)) DEALLOCATE (resdist)
    IF (ALLOCATED  (tot_watveg_beg)) DEALLOCATE (tot_watveg_beg)
    IF (ALLOCATED  (tot_watveg_end)) DEALLOCATE (tot_watveg_end)
    IF (ALLOCATED  (tot_watsoil_beg)) DEALLOCATE (tot_watsoil_beg)
    IF (ALLOCATED  (tot_watsoil_end)) DEALLOCATE (tot_watsoil_end)
    IF (ALLOCATED  (delsoilmoist)) DEALLOCATE (delsoilmoist)
    IF (ALLOCATED  (delintercept)) DEALLOCATE (delintercept)
    IF (ALLOCATED  (snow_beg)) DEALLOCATE (snow_beg)
    IF (ALLOCATED  (snow_end)) DEALLOCATE (snow_end)
    IF (ALLOCATED  (delswe)) DEALLOCATE (delswe)
    IF (ALLOCATED  (undermcr)) DEALLOCATE (undermcr)
    IF (ALLOCATED  (v1)) DEALLOCATE (v1)
    IF (ALLOCATED  (humtot)) DEALLOCATE (humtot)
    IF (ALLOCATED  (resolv)) DEALLOCATE (resolv)
    IF (ALLOCATED  (k)) DEALLOCATE (k)
    IF (ALLOCATED  (kk)) DEALLOCATE (kk)
    IF (ALLOCATED  (kk_moy)) DEALLOCATE (kk_moy)
    IF (ALLOCATED  (avan_mod_tab)) DEALLOCATE (avan_mod_tab)
    IF (ALLOCATED  (nvan_mod_tab)) DEALLOCATE (nvan_mod_tab)
    IF (ALLOCATED  (a)) DEALLOCATE (a)
    IF (ALLOCATED  (b)) DEALLOCATE (b)
    IF (ALLOCATED  (d)) DEALLOCATE (d)
    IF (ALLOCATED  (e)) DEALLOCATE (e)
    IF (ALLOCATED  (f)) DEALLOCATE (f)
    IF (ALLOCATED  (g1)) DEALLOCATE (g1)
    IF (ALLOCATED  (ep)) DEALLOCATE (ep)
    IF (ALLOCATED  (fp)) DEALLOCATE (fp)
    IF (ALLOCATED  (gp)) DEALLOCATE (gp)
    IF (ALLOCATED  (rhs)) DEALLOCATE (rhs)
    IF (ALLOCATED  (srhs)) DEALLOCATE (srhs)
    IF (ALLOCATED  (tmc)) DEALLOCATE (tmc)
    IF (ALLOCATED  (tmcs)) DEALLOCATE (tmcs)
    IF (ALLOCATED  (tmcr)) DEALLOCATE (tmcr)
    IF (ALLOCATED  (tmcfc)) DEALLOCATE (tmcfc)
    IF (ALLOCATED  (tmcw)) DEALLOCATE (tmcw)
    IF (ALLOCATED  (tmc_litter)) DEALLOCATE (tmc_litter)
!gmjc top 5 layer mc for grazing
    IF (ALLOCATED  (tmc_trampling)) DEALLOCATE (tmc_trampling)
!end gmjc
    IF (ALLOCATED  (tmc_litt_mea)) DEALLOCATE (tmc_litt_mea)
    IF (ALLOCATED  (tmc_litter_res)) DEALLOCATE (tmc_litter_res)
    IF (ALLOCATED  (tmc_litter_wilt)) DEALLOCATE (tmc_litter_wilt)
    IF (ALLOCATED  (tmc_litter_field)) DEALLOCATE (tmc_litter_field)
    IF (ALLOCATED  (tmc_litter_sat)) DEALLOCATE (tmc_litter_sat)
    IF (ALLOCATED  (tmc_litter_awet)) DEALLOCATE (tmc_litter_awet)
    IF (ALLOCATED  (tmc_litter_adry)) DEALLOCATE (tmc_litter_adry)
    IF (ALLOCATED  (tmc_litt_wet_mea)) DEALLOCATE (tmc_litt_wet_mea)
    IF (ALLOCATED  (tmc_litt_dry_mea)) DEALLOCATE (tmc_litt_dry_mea)
    IF (ALLOCATED  (ru_ns)) DEALLOCATE (ru_ns)
    IF (ALLOCATED  (dr_ns)) DEALLOCATE (dr_ns)
    IF (ALLOCATED  (tr_ns)) DEALLOCATE (tr_ns)
    IF (ALLOCATED  (vegetmax_soil)) DEALLOCATE (vegetmax_soil)
    IF (ALLOCATED  (mc)) DEALLOCATE (mc)
    IF (ALLOCATED  (soilmoist)) DEALLOCATE (soilmoist)
    IF (ALLOCATED  (soilmoist_liquid)) DEALLOCATE (soilmoist_liquid)
    IF (ALLOCATED  (soil_wet_ns)) DEALLOCATE (soil_wet_ns)
    IF (ALLOCATED  (soil_wet_litter)) DEALLOCATE (soil_wet_litter)
    IF (ALLOCATED  (qflux)) DEALLOCATE (qflux)
    IF (ALLOCATED  (tmat)) DEALLOCATE (tmat)
    IF (ALLOCATED  (stmat)) DEALLOCATE (stmat)
    IF (ALLOCATED  (nroot)) DEALLOCATE (nroot)
    IF (ALLOCATED  (kfact_root)) DEALLOCATE (kfact_root)
    IF (ALLOCATED  (kfact)) DEALLOCATE (kfact)
    IF (ALLOCATED  (zz)) DEALLOCATE (zz)
    IF (ALLOCATED  (dz)) DEALLOCATE (dz)
    IF (ALLOCATED  (dh)) DEALLOCATE (dh)
    IF (ALLOCATED  (mc_lin)) DEALLOCATE (mc_lin)
    IF (ALLOCATED  (k_lin)) DEALLOCATE (k_lin)
    IF (ALLOCATED  (d_lin)) DEALLOCATE (d_lin)
    IF (ALLOCATED  (a_lin)) DEALLOCATE (a_lin)
    IF (ALLOCATED  (b_lin)) DEALLOCATE (b_lin)
    IF ( ALLOCATED (refSOC_1d)) DEALLOCATE (refSOC_1d)

    IF (TOPM_calcul) CALL topmodel_clear()

  END SUBROUTINE hydrol_clear

!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_tmc_update
!!
!>\BRIEF        This routine updates the soil moisture profiles when the vegetation fraction have changed. 
!!
!! DESCRIPTION  :
!! 
!!    This routine update tmc and mc with variation of veget_max (LAND_USE or DGVM activated)
!! 
!!
!!
!!
!! RECENT CHANGE(S) : Adaptation to excluding nobio from soiltile(1)
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_tmc_update
  SUBROUTINE hydrol_tmc_update ( kjpindex, veget_max, soiltile, qsintveg, drain_upd, runoff_upd)

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                            :: kjpindex      !! domain size
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)     :: veget_max     !! max fraction of vegetation type
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT (in)   :: soiltile      !! Fraction of each soil tile (0-1, unitless)

    !! 0.2 Output variables
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)         :: drain_upd        !! Change in drainage due to decrease in vegtot
                                                                              !! on mc [kg/m2/dt]
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)         :: runoff_upd       !! Change in runoff due to decrease in vegtot
                                                                              !! on water2infilt[kg/m2/dt]
    
    !! 0.3 Modified variables
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout)  :: qsintveg   !! Amount of water in the canopy interception 

    !! 0.4 Local variables
    INTEGER(i_std)                           :: ji, jv, jst,jsl
    LOGICAL                                  :: soil_upd        !! True if soiltile changed since last time step
    LOGICAL                                  :: vegtot_upd      !! True if vegtot changed since last time step
    LOGICAL                                  :: error=.FALSE.   !! If true, exit in the end of subroutine
    REAL(r_std), DIMENSION(kjpindex,nstm)    :: vmr             !! Change in soiltile (within vegtot)
    REAL(r_std), DIMENSION(kjpindex)         :: vmr_sum
    REAL(r_std), DIMENSION(kjpindex)         :: delvegtot    
    REAL(r_std), DIMENSION(kjpindex,nslm)    :: mc_dilu         !! Total loss of moisture content
    REAL(r_std), DIMENSION(kjpindex)         :: infil_dilu      !! Total loss for water2infilt
    REAL(r_std), DIMENSION(kjpindex,nstm)    :: tmc_old         !! tmc before calculations
    REAL(r_std), DIMENSION(kjpindex,nstm)    :: water2infilt_old!! water2infilt before calculations
    REAL(r_std), DIMENSION (kjpindex,nvm)    :: qsintveg_old    !! qsintveg before calculations
    REAL(r_std), DIMENSION(kjpindex)         :: test
    REAL(r_std), DIMENSION(kjpindex,nslm,nstm) :: mcaux        !! serves to hold the chnage in mc when vegtot decreases

    !! 0. For checks

    IF (check_cwrr) THEN
       ! Save soil moisture for later use
       tmc_old(:,:) = tmc(:,:) 
       water2infilt_old(:,:) = water2infilt(:,:)
       qsintveg_old(:,:) = qsintveg(:,:)
    ENDIF
    
    !! 1. If a PFT has disapperead as result from a veget_max change, 
    !!    then add canopy water to surface water.
    !     Other adaptations of qsintveg are delt by the normal functioning of hydrol_canop

    DO ji=1,kjpindex
       IF (vegtot_old(ji) .GT.min_sechiba) THEN
          DO jv=1,nvm
             IF ((veget_max(ji,jv).LT.min_sechiba).AND.(qsintveg(ji,jv).GT.0.)) THEN
                jst=pref_soil_veg(jv) ! soil tile index
                water2infilt(ji,jst) = water2infilt(ji,jst) + qsintveg(ji,jv)/(resdist(ji,jst)*vegtot_old(ji))
                qsintveg(ji,jv) = zero
             ENDIF
          ENDDO
       ENDIF
    ENDDO
   
    !! 2. We now deal with the changes of soiltile and corresponding soil moistures
    !!    Because sum(soiltile)=1 whatever vegtot, we need to distinguish two cases:
    !!    - when vegtot changes (meaning that the nobio fraction changes too),
    !!    - and when vegtot does not changes (a priori the most frequent case)

    vegtot_upd = SUM(ABS((vegtot(:)-vegtot_old(:)))) .GT. zero ! True if at least one land point with a vegtot change
    runoff_upd(:) = zero
    drain_upd(:) = zero
    IF (vegtot_upd) THEN
       ! We find here the processing specific to the chnages of nobio fraction and vegtot

       delvegtot(:) = vegtot(:) - vegtot_old(:)

       DO jst=1,nstm
          DO ji=1,kjpindex

             IF (delvegtot(ji) .GT. min_sechiba) THEN

                !! 2.1. If vegtot increases (nobio decreases), then the mc in each soiltile is decreased
                !!      assuming the same proportions for each soiltile, and each soil layer
                
                mc(ji,:,jst) = mc(ji,:,jst) * vegtot_old(ji)/vegtot(ji) ! vegtot cannot be zero as > vegtot_old
                water2infilt(ji,jst) = water2infilt(ji,jst) * vegtot_old(ji)/vegtot(ji)

             ELSE

                !! 2.2 If vegtot decreases (nobio increases), then the mc in each soiltile should increase,
                !!     but should not exceed mcs
                !!     For simplicity, we choose to send the corresponding water volume to drainage
                !!     We do the same for water2infilt but send the excess to surface runoff

                IF (vegtot(ji) .GT.min_sechiba) THEN
                   mcaux(ji,:,jst) =  mc(ji,:,jst) * (vegtot_old(ji)-vegtot(ji))/vegtot(ji) ! mcaux is the delta mc
                ELSE ! we just have nobio in the grid-cell
                   mcaux(ji,:,jst) =  mc(ji,:,jst)
                ENDIF
                
                drain_upd(ji) = drain_upd(ji) + dz(2) * ( trois*mcaux(ji,1,jst) + mcaux(ji,2,jst) )/huit
                DO jsl = 2,nslm-1
                   drain_upd(ji) = drain_upd(ji) + dz(jsl) * (trois*mcaux(ji,jsl,jst)+mcaux(ji,jsl-1,jst))/huit &
                        + dz(jsl+1) * (trois*mcaux(ji,jsl,jst)+mcaux(ji,jsl+1,jst))/huit
                ENDDO
                drain_upd(ji) = drain_upd(ji) + dz(nslm) * (trois*mcaux(ji,nslm,jst) + mcaux(ji,nslm-1,jst))/huit

                IF (vegtot(ji) .GT.min_sechiba) THEN
                   runoff_upd(ji) = runoff_upd(ji) + water2infilt(ji,jst) * (vegtot_old(ji)-vegtot(ji))/vegtot(ji)
                ELSE ! we just have nobio in the grid-cell
                   runoff_upd(ji) = runoff_upd(ji) + water2infilt(ji,jst)
                ENDIF

             ENDIF
             
          ENDDO
       ENDDO
       
    ENDIF
    
    !! 3. At the end of step 2, we are back to a case where vegtot changes are treated, so we can use soiltile
    !!    as a fraction of vegtot to process the mc transfers between soil tiles due to the changes of vegetation map
   
    !! 3.1 Check if soiltiles changed since last time step
    soil_upd=SUM(ABS(soiltile(:,:)-resdist(:,:))) .GT. zero
    IF (printlev>=3) WRITE (numout,*) 'soil_upd ', soil_upd
        
    IF (soil_upd) THEN
     
       !! 3.2 Define the change in soiltile
       vmr(:,:) = soiltile(:,:) - resdist(:,:)  ! resdist is the previous values of soiltiles, previous timestep, so before new map

       ! Total area loss by the three soil tiles
       DO ji=1,kjpindex
          vmr_sum(ji)=SUM(vmr(ji,:),MASK=vmr(ji,:).LT.zero)
       ENDDO

       !! 3.3 Shrinking soil tiles
       !! 3.3.1 Total loss of moisture content from the shrinking soil tiles, expressed by soil layer
       mc_dilu(:,:)=zero
       DO jst=1,nstm
          DO jsl = 1, nslm
             DO ji=1,kjpindex
                IF ( vmr(ji,jst) < -min_sechiba ) THEN
                   mc_dilu(ji,jsl) = mc_dilu(ji,jsl) + mc(ji,jsl,jst) * vmr(ji,jst) / vmr_sum(ji)
                ENDIF
             ENDDO
          ENDDO
       ENDDO

       !! 3.3.2 Total loss of water2inft from the shrinking soil tiles
       infil_dilu(:)=zero
       DO jst=1,nstm
          DO ji=1,kjpindex
             IF ( vmr(ji,jst) < -min_sechiba ) THEN
                infil_dilu(ji) = infil_dilu(ji) + water2infilt(ji,jst) * vmr(ji,jst) / vmr_sum(ji)
             ENDIF
          ENDDO
       ENDDO

       !! 3.4 Each gaining soil tile gets moisture proportionally to both the total loss and its areal increase 

       ! As the original mc from each soil tile are in [mcr,mcs] and we do weighted avrage, the new mc are in [mcr,mcs]
       ! The case where the soiltile is created (soiltile_old=0) works as the other cases

       ! 3.4.1 Update mc(kjpindex,nslm,nstm) !m3/m3
       DO jst=1,nstm
          DO jsl = 1, nslm
             DO ji=1,kjpindex
                IF ( vmr(ji,jst) > min_sechiba ) THEN
                   mc(ji,jsl,jst) = ( mc(ji,jsl,jst) * resdist(ji,jst) + mc_dilu(ji,jsl) * vmr(ji,jst) ) / soiltile(ji,jst)
                   ! NB : soiltile can not be zero for case vmr > zero, see slowproc_veget
                ENDIF
             ENDDO
          ENDDO
       ENDDO
!
!       DO jst=1,nstm
!          IF ( vmr(1,jst) > zero ) THEN
!             WRITE(numout,*) 'zdcheck2 jst=',jst,'soiltile,resdist,vmr',soiltile(1,jst),resdist(1,jst),vmr(1,jst)
!          ENDIF
!       ENDDO
       
       ! 3.4.2 Update water2inft
       DO jst=1,nstm
          DO ji=1,kjpindex
             IF ( vmr(ji,jst) > min_sechiba ) THEN !donc soiltile>0     
                water2infilt(ji,jst) = ( water2infilt(ji,jst) * resdist(ji,jst) + infil_dilu(ji) * vmr(ji,jst) ) / soiltile(ji,jst)
             ENDIF !donc resdist>0
          ENDDO
       ENDDO

       ! 3.4.3 Case where soiltile < min_sechiba 
       DO jst=1,nstm
          DO ji=1,kjpindex
             IF ( soiltile(ji,jst) .LT. min_sechiba ) THEN
                water2infilt(ji,jst) = zero
                mc(ji,:,jst) = zero
             ENDIF
          ENDDO
       ENDDO

    ENDIF ! soil_upd

    !! 4. Update tmc and humtot
    
    DO jst=1,nstm
       DO ji=1,kjpindex
             tmc(ji,jst) = dz(2) * ( trois*mc(ji,1,jst) + mc(ji,2,jst) )/huit
             DO jsl = 2,nslm-1
                tmc(ji,jst) = tmc(ji,jst) + dz(jsl) * (trois*mc(ji,jsl,jst)+mc(ji,jsl-1,jst))/huit &
                     + dz(jsl+1) * (trois*mc(ji,jsl,jst)+mc(ji,jsl+1,jst))/huit
             ENDDO
             tmc(ji,jst) = tmc(ji,jst) + dz(nslm) * (trois*mc(ji,nslm,jst) + mc(ji,nslm-1,jst))/huit
             tmc(ji,jst) = tmc(ji,jst) + water2infilt(ji,jst)
             ! WARNING tmc is increased by includes water2infilt(ji,jst)
       ENDDO
    ENDDO

    humtot(:) = zero
    DO jst=1,nstm
       DO ji=1,kjpindex
          humtot(ji) = humtot(ji) + vegtot(ji) * soiltile(ji,jst) * tmc(ji,jst) ! average over grid-cell (i.e. total land)
       ENDDO
    ENDDO

    !! 5. Check
    IF (check_cwrr) THEN
       DO ji=1,kjpindex
          test(ji) = SUM(tmc(ji,:)*soiltile(ji,:)*vegtot(ji)) - SUM(tmc_old(ji,:)*resdist(ji,:)*vegtot_old(ji)) + &
               SUM(qsintveg(ji,:)) - SUM(qsintveg_old(ji,:)) + (drain_upd(ji) + runoff_upd(ji))    
          IF ( ABS(test(ji)) .GT.  10.*allowed_err ) THEN
             WRITE(numout,*) 'tmc update WRONG: ji',ji
             WRITE(numout,*) 'tot water avant:',SUM(tmc_old(ji,:)*resdist(ji,:)*vegtot_old(ji)) + SUM(qsintveg_old(ji,:))
             WRITE(numout,*) 'tot water apres:',SUM(tmc(ji,:)*soiltile(ji,:)*vegtot(ji)) + SUM(qsintveg(ji,:))
             WRITE(numout,*) 'err:',test(ji)
             WRITE(numout,*) 'allowed_err:',allowed_err
             WRITE(numout,*) 'tmc:',tmc(ji,:)
             WRITE(numout,*) 'tmc_old:',tmc_old(ji,:)
             WRITE(numout,*) 'qsintveg:',qsintveg(ji,:)
             WRITE(numout,*) 'qsintveg_old:',qsintveg_old(ji,:)
             WRITE(numout,*) 'SUMqsintveg:',SUM(qsintveg(ji,:))
             WRITE(numout,*) 'SUMqsintveg_old:',SUM(qsintveg_old(ji,:))
             WRITE(numout,*) 'veget_max:',veget_max(ji,:)
             WRITE(numout,*) 'soiltile:',soiltile(ji,:)
             WRITE(numout,*) 'resdist:',resdist(ji,:)
             WRITE(numout,*) 'vegtot:',vegtot(ji)
             WRITE(numout,*) 'vegtot_old:',vegtot_old(ji)
             WRITE(numout,*) 'drain_upd:',drain_upd(ji)
             WRITE(numout,*) 'runoff_upd:',runoff_upd(ji)
             WRITE(numout,*) 'vmr:',vmr(ji,:)
             WRITE(numout,*) 'vmr_sum:',vmr_sum(ji)
             DO jst=1,nstm
                WRITE(numout,*) 'mc(',jst,'):',mc(ji,:,jst)
             ENDDO
             WRITE(numout,*) 'water2infilt:',water2infilt(ji,:)
             WRITE(numout,*) 'water2infilt_old:',water2infilt_old(ji,:)
             WRITE(numout,*) 'infil_dilu:',infil_dilu(ji)
             WRITE(numout,*) 'mc_dilu:',mc_dilu(ji,:)

             error=.TRUE.
             CALL ipslerr_p(2, 'hydrol_tmc_update', 'Error in water balance', 'We STOP in the end of this subroutine','')
          ENDIF
       ENDDO
    ENDIF

    !! Now that the work is done, update resdist
    resdist(:,:) = soiltile(:,:)

    !
    !!  Exit if error was found previously in this subroutine
    !
    IF ( error ) THEN
       WRITE(numout,*) 'One or more errors have been detected in hydrol_tmc_update. Model stops.'
       CALL ipslerr_p(3, 'hydrol_tmc_update', 'We will STOP now.',&
                  & 'One or several fatal errors were found previously.','')
    END IF

    IF (printlev>=3) WRITE (numout,*) ' hydrol_tmc_update done '

  END SUBROUTINE hydrol_tmc_update


  SUBROUTINE hydrol_rotation_update( ip, kjpindex, rot_matrix, old_veget_max, veget_max, soiltile, qsintveg )

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                          :: ip, kjpindex      !! domain size
    REAL(r_std),DIMENSION (nvm), INTENT (in)            :: old_veget_max     !! max fraction of vegetation type
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)   :: veget_max     !! max fraction of vegetation type
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT (in) :: soiltile      !! Fraction of each soil tile (0-1, unitless)
    REAL(r_std), DIMENSION (nvm, nvm), INTENT(in)       :: rot_matrix    !! rotation matrix

    !! 0.3 Modified variables
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout)   :: qsintveg   !! Amount of water in the canopy interception 
!    REAL(r_std), DIMENSION (kjpindex, nstm), INTENT(inout) :: resdist    !! Soiltile from previous time-step
!  resdist is defined in MODULE hydrol, no need to be put in the argument list

    !! 0.4 Local variables
    INTEGER(i_std)                           :: ji, jv, jst, jst1, jst2, jsl, jsrc, jtar
    LOGICAL                                  :: soil_upd        !! True if soiltile changed since last time step
    LOGICAL                                  :: error=.FALSE.   !! If true, exit in the end of subroutine
!!    REAL(r_std), DIMENSION(kjpindex,nstm)    :: vmr             !! Change in soiltile
!!    REAL(r_std), DIMENSION(kjpindex)         :: vmr_sum
!!    REAL(r_std), DIMENSION(kjpindex,nslm)    :: mc_dilu         !! Total loss of moisture content
!!    REAL(r_std), DIMENSION(kjpindex)         :: infil_dilu      !! Total loss for water2infilt
!!    REAL(r_std), DIMENSION(kjpindex,nstm)    :: tmc_old         !! tmc before calculations
!!    REAL(r_std), DIMENSION(kjpindex,nstm)    :: water2infilt_old!! water2infilt before calculations
!!    REAL(r_std), DIMENSION (kjpindex,nvm)    :: qsintveg_old    !! qsintveg before calculations
!!    REAL(r_std), DIMENSION(kjpindex)         :: test

    REAL(r_std), DIMENSION(nslm,nstm)    :: mc_dilu         !! Total loss of moisture content
    REAL(r_std), DIMENSION(nslm,nstm)    :: mc_old          !! temporary file to store mc
    REAL(r_std), DIMENSION(nstm)         :: infil_dilu      !! Total loss for water2infilt
    REAL(r_std), DIMENSION(nstm)    :: tmc_old         !! tmc before calculations
    REAL(r_std), DIMENSION(nstm)    :: water2infilt_old!! water2infilt before calculations
    REAL(r_std), DIMENSION(nvm)     :: qsintveg_old    !! qsintveg before calculations
    REAL(r_std), DIMENSION(nvm)     :: maxfrac, maxfrac_new
    REAL(r_std), DIMENSION(nstm,nstm) :: rot_matrix_tile !! relative portion of soil tile to be transferred
    REAL(r_std)                     :: test

    !! 0. Check if soiltiles changed since last time step
!    soil_upd=SUM(ABS(soiltile(:,:)-resdist(:,:))) .GT. zero
    maxfrac = old_veget_max(:)
    maxfrac_new = old_veget_max(:)
    rot_matrix_tile(:,:) = 0.0
    DO jsrc = 1,nvm
        DO jtar = 1,nvm
            IF (rot_matrix(jsrc,jtar) .GT. 0.0) THEN
                maxfrac_new(jtar) = maxfrac_new(jtar) + maxfrac(jsrc) * rot_matrix(jsrc,jtar)
                maxfrac_new(jsrc) = maxfrac_new(jsrc) - maxfrac(jsrc) * rot_matrix(jsrc,jtar)
                jst1 = pref_soil_veg(jsrc)
                jst2 = pref_soil_veg(jtar)
                rot_matrix_tile(jst1,jst2) = rot_matrix_tile(jst1,jst2) + &
                    rot_matrix(jsrc,jtar) * maxfrac(jsrc) / resdist(ip, jst1 )
                IF ( (rot_matrix_tile(jst1,jst2) .GT. 1.0) .OR. &
                     (rot_matrix_tile(jst1,jst2) .LE. 0.0) ) THEN
                    !!!! make sure the fraction is in (0,1)
                    WRITE(numout,*) 'pref_soil_veg',pref_soil_veg
                    WRITE(numout,*) 'jsrc, jtar,', jsrc, jtar
                    WRITE(numout,*) 'jst1, jst2,', jst1, jst2
                    WRITE(numout,*) 'maxfrac(jsrc), rot_matrix(jsrc,jtar)', maxfrac(jsrc), rot_matrix(jsrc,jtar)
                    WRITE(numout,*) 'pref_soil_veg(jsrc), resdist(ip,pref_soil_veg(jsrc))', pref_soil_veg(jsrc), resdist(ip,pref_soil_veg(jsrc))
                    STOP 'soiltile error in hydrol_rotation'
                ENDIF
            ENDIF
        ENDDO
    ENDDO

    !!!! check if the vegetation conversion is successful
    IF ( SUM(ABS(maxfrac_new - veget_max(ip,:))) .GT. min_sechiba ) THEN
        WRITE(numout,*) 'maxfrac',maxfrac
        WRITE(numout,*) 'maxfrac_new', maxfrac_new
        WRITE(numout,*) 'veget_max(ip,:)', veget_max(ip,:)
        STOP 'hydrol_rotation: fraction conversion error'
    ENDIF
    IF (printlev>=4) THEN 
        WRITE(numout,*) 'xuhui: hydrol_rotation' 
        WRITE(numout,*) 'resdist(ip,:)', resdist(ip,:)
        WRITE(numout,*) 'soiltile(ip,:)', soiltile(ip,:)
        DO jsrc = 1,nstm 
            WRITE(numout,*) 'jsrc, rot_matrix_tile(jsrc,:)', jsrc, rot_matrix_tile(jsrc,:)
        ENDDO
    ENDIF

    IF (SUM(rot_matrix) .GT. 0) THEN
        soil_upd = .TRUE.
    ENDIF
    IF (printlev>=3) WRITE (numout,*) 'soil_upd ', soil_upd

    IF (check_cwrr) THEN
       ! Save soil moisture for later use
       tmc_old(:) = tmc(ip,:) 
       water2infilt_old(:) = water2infilt(ip,:)
       qsintveg_old(:) = qsintveg(ip,:)
    ENDIF

    !! 1. If a PFT has disapperead as result from a veget_max change, 
    !!    then add canopy water to surface water.
     DO jv=1,nvm
        IF ( (maxfrac_new(jv) .LT. min_sechiba) .AND.  (qsintveg(ip, jv) .GT. 0.0) )  THEN
            jst = pref_soil_veg(jv)
            water2infilt(ip,jst) = water2infilt(ip,jst) + qsintveg(ip,jv)/maxfrac(jv)
            qsintveg(ip,jv) = zero
        ENDIF
     ENDDO
    
    mc_old = mc(ip,:,:)
    water2infilt_old(:) = water2infilt(ip,:)
    !! 2. Compute new soil moisture if soiltile changed
    IF (soil_upd) THEN
        DO jtar = 1,nstm
          mc_dilu(:,:) = zero
          infil_dilu(:) = zero
          IF ( SUM(rot_matrix_tile(:,jtar)) .GT. min_sechiba ) THEN
            DO jsrc = 1,nstm
              IF ( rot_matrix_tile(jsrc,jtar) .GT. min_sechiba ) THEN
                mc_dilu(:,jsrc) = mc_old(:,jsrc)
                infil_dilu(jsrc) = water2infilt_old(jsrc)
              ENDIF ! rot_matrix_tile(jsrc,jtar) > 0
            ENDDO
            !!! actually do the rotation
            mc(ip,:,jtar) = mc_old(:,jtar) * resdist(ip,jtar) * (1.0 - SUM(rot_matrix_tile(jtar,:)))
            water2infilt(ip,jtar) = water2infilt_old(jtar) * resdist(ip,jtar) * (1.0 - SUM(rot_matrix_tile(jtar,:)))
            DO jsrc = 1,nstm
                mc(ip,:,jtar) = mc(ip,:,jtar) + resdist(ip,jsrc) * rot_matrix_tile(jsrc,jtar) * mc_dilu(:,jsrc)
                water2infilt(ip,jtar) = water2infilt(ip,jtar) + resdist(ip,jsrc) * rot_matrix_tile(jsrc,jtar) * infil_dilu(jsrc)
            ENDDO
            IF ( soiltile(ip,jtar) .LE. 0. ) THEN
                WRITE(numout,*) 'jtar, soiltile(ip,jtar)',jtar, soiltile(ip,jtar)
                STOP 'hydrol_rotation_update: target tile has no proportion'
            ENDIF
            mc(ip,:,jtar) = mc(ip,:,jtar) / soiltile(ip,jtar)
            water2infilt(ip,jtar) = water2infilt(ip,jtar) / soiltile(ip,jtar)
          ENDIF ! SUM(rot_matrix_tile(:,jtar)) > 0
        ENDDO

       ! 2.3.3 Case where soiltile < min_sechiba 
       DO jst=1,nstm
          IF ( soiltile(ip,jst) .LT. min_sechiba ) THEN
             water2infilt(ip,jst) = zero
             mc(ip,:,jst) = zero
          ENDIF
       ENDDO

       IF (printlev>=4) THEN
            WRITE(numout,*) 'mc_old(1,:)',mc_old(1,:)
            WRITE(numout,*) 'mc(ip,1,:)',mc(ip,1,:)
            WRITE(numout,*) 'water2infilt_old(:)', water2infilt_old(:)
            WRITE(numout,*) 'water2infilt(ip,:)', water2infilt(ip,:)
       ENDIF


    ENDIF ! soil_upd 


    !2.3.3 we compute tmc(kjpindex,nstm) and humtot!
    DO jst=1,nstm
         tmc(ip,jst) = dz(2) * ( trois*mc(ip,1,jst) + mc(ip,2,jst) )/huit
         DO jsl = 2,nslm-1
            tmc(ip,jst) = tmc(ip,jst) + dz(jsl) * (trois*mc(ip,jsl,jst)+mc(ip,jsl-1,jst))/huit &
                 + dz(jsl+1) * (trois*mc(ip,jsl,jst)+mc(ip,jsl+1,jst))/huit
         ENDDO
         tmc(ip,jst) = tmc(ip,jst) + dz(nslm) * (trois*mc(ip,nslm,jst) + mc(ip,nslm-1,jst))/huit
         tmc(ip,jst) = tmc(ip,jst) + water2infilt(ip,jst)
         ! WARNING tmc is increased by water2infilt(ip,jst), but mc is not modified !
    ENDDO

    humtot(ip) = zero
    DO jst=1,nstm
        humtot(ip) = humtot(ip) + soiltile(ip,jst) * tmc(ip,jst)
    ENDDO

    !! 4 check
    IF (check_cwrr) THEN
!       DO ji=1,kjpindex
        ji = ip
          test = ABS(SUM(tmc(ji,:)*soiltile(ji,:)) - SUM(tmc_old(:)*resdist(ji,:)) + &
               SUM(qsintveg(ji,:)) - SUM(qsintveg_old(:))) ! sum(soiltile)=1
          IF ( test .GT.  allowed_err ) THEN
             WRITE(numout,*) 'hydrol_rotation_update WRONG: ji',ji
             WRITE(numout,*) 'tot water before:',SUM(tmc_old(:)*resdist(ji,:)) + SUM(qsintveg_old(:))
             WRITE(numout,*) 'tot water after:',SUM(tmc(ji,:)*soiltile(ji,:)) + SUM(qsintveg(ji,:))
             WRITE(numout,*) 'err:',test
             WRITE(numout,*) 'allowed_err:',allowed_err
             WRITE(numout,*) 'tmc:',tmc(ji,:)
             WRITE(numout,*) 'tmc_old:',tmc_old(:)
             WRITE(numout,*) 'qsintveg:',qsintveg(ji,:)
             WRITE(numout,*) 'qsintveg_old:',qsintveg_old(:)
             WRITE(numout,*) 'SUMqsintveg:',SUM(qsintveg(ji,:))
             WRITE(numout,*) 'SUMqsintveg_old:',SUM(qsintveg_old(:))
             WRITE(numout,*) 'veget_max:',veget_max(ji,:)
             WRITE(numout,*) 'soiltile:',soiltile(ji,:)
             WRITE(numout,*) 'resdist:',resdist(ji,:)
             DO jst=1,nstm
                WRITE(numout,*) 'mc(',jst,'):',mc(ji,:,jst)
             ENDDO
             WRITE(numout,*) 'water2infilt:',water2infilt(ji,:)
             WRITE(numout,*) 'water2infilt_old:',water2infilt_old(:)

             error=.TRUE.
             CALL ipslerr_p(2, 'hydrol_rotation_update', 'Error in water balance', 'We STOP in the end of this subroutine','')
          ENDIF
!       ENDDO
    ENDIF

    !! Now that the work is done, update resdist
    resdist(:,:) = soiltile(:,:)

    !
    !!  Exit if error was found previously in this subroutine
    !
    IF ( error ) THEN
       WRITE(numout,*) 'One or more errors have been detected in hydrol_tmc_update. Model stops.'
       CALL ipslerr_p(3, 'hydrol_tmc_update', 'We will STOP now.',&
                  & 'One or several fatal errors were found previously.','')
    END IF

    IF (printlev>=3) WRITE (numout,*) ' hydrol_rotation_update done '

  END SUBROUTINE hydrol_rotation_update

!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_var_init
!!
!>\BRIEF        This routine initializes hydrologic parameters to define K and D, and diagnostic hydrologic variables.  
!!
!! DESCRIPTION  :
!! - 1 compute the depths
!! - 2 compute the profile for roots
!! - 3 compute the profile for a and n Van Genuchten parameter
!! - 4 compute the linearized values of k, a, b and d for the resolution of Fokker Planck equation
!! - 5 water reservoirs initialisation
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_var_init

  SUBROUTINE hydrol_var_init (kjpindex, veget, veget_max, soiltile, njsc, &
       mx_eau_var, shumdiag_perma, &
       drysoil_frac, qsintveg, mc_layh, mcl_layh, mc_layh_s, mcl_layh_s, &
!gmjc
       tmc_topgrass, humcste_use, altmax, swdown) !!xchen add swdown
!end gmjc 

    ! interface description

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    ! input scalar 
    INTEGER(i_std), INTENT(in)                          :: kjpindex      !! Domain size (number of grid cells) (1)
    ! input fields
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)   :: veget_max     !! PFT fractions within grid-cells (1; 1)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)   :: veget         !! Effective fraction of vegetation by PFT (1; 1)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)    :: njsc          !! Index of the dominant soil textural class 
                                                                         !! in the grid cell (1-nscm, unitless) 
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT (in) :: soiltile      !! Fraction of each soil tile within vegtot (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: swdown           !!xchen Down-welling surface short-wave flux  


    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: mx_eau_var    !! Maximum water content of the soil 
                                                                         !! @tex $(kg m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (out) :: shumdiag_perma!! Percent of porosity filled with water (mc/mcs)
                                                                         !! used for the thermal computations
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)    :: drysoil_frac  !! function of litter humidity
    REAL(r_std), DIMENSION (kjpindex,nslm), INTENT (out):: mc_layh       !! Volumetric soil moisture content for each layer in hydrol(liquid+ice) [m3/m3]
    REAL(r_std), DIMENSION (kjpindex,nslm,nstm), INTENT (out):: mc_layh_s   !! Volumetric soil moisture content for each layer in hydrol(liquid+ice) [m3/m3]
    REAL(r_std), DIMENSION (kjpindex,nslm), INTENT (out):: mcl_layh      !! Volumetric soil moisture content for each layer in hydrol(liquid) [m3/m3]
    REAL(r_std), DIMENSION (kjpindex,nslm,nstm), INTENT (out):: mcl_layh_s  !! Volumetric soil moisture content for each layer in hydrol(liquid) [m3/m3]

    !! 0.3 Modified variables
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout)  :: qsintveg    !! Water on vegetation due to interception
                                                                         !! @tex $(kg m^{-2})$ @endtex  
!gmjc top 5 layer grassland soil moisture for grazing
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)       :: tmc_topgrass
!end gmjc
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)  :: humcste_use
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)   :: altmax
    !! 0.4 Local variables

    INTEGER(i_std)                                      :: ji, jv, jp    !! Grid-cell and PFT indices (1)
    INTEGER(i_std)                                      :: jst, jsc, jsl !! Soiltile, Soil Texture, and Soil layer indices (1)
    INTEGER(i_std)                                      :: i             !! Index (1)
    REAL(r_std)                                         :: m             !! m=1-1/n (unitless)
    REAL(r_std)                                         :: frac          !! Relative linearized VWC (unitless)
    REAL(r_std)                                         :: avan_mod      !! VG parameter a modified from  exponantial profile
                                                                         !! @tex $(mm^{-1})$ @endtex
    REAL(r_std)                                         :: nvan_mod      !! VG parameter n  modified from  exponantial profile
                                                                         !! (unitless)
    REAL(r_std), DIMENSION(nslm,nscm)                   :: afact, nfact  !! Multiplicative factor for decay of a and n with depth
                                                                         !! (unitless)
    ! parameters for "soil densification" with depth
    REAL(r_std)                                         :: dp_comp       !! Depth at which the 'compacted' value of ksat
                                                                         !! is reached (m)
    REAL(r_std)                                         :: f_ks          !! Exponential factor for decay of ksat with depth 
                                                                         !! @tex $(m^{-1})$ @endtex 
    ! Fixed parameters from fitted relationships
    REAL(r_std)                                         :: n0            !! fitted value for relation log((n-n0)/(n_ref-n0)) = 
                                                                         !! nk_rel * log(k/k_ref) 
                                                                         !! (unitless)
    REAL(r_std)                                         :: nk_rel        !! fitted value for relation log((n-n0)/(n_ref-n0)) = 
                                                                         !! nk_rel * log(k/k_ref) 
                                                                         !! (unitless)
    REAL(r_std)                                         :: a0            !! fitted value for relation log((a-a0)/(a_ref-a0)) = 
                                                                         !! ak_rel * log(k/k_ref) 
                                                                         !! @tex $(mm^{-1})$ @endtex
    REAL(r_std)                                         :: ak_rel        !! fitted value for relation log((a-a0)/(a_ref-a0)) = 
                                                                         !! ak_rel * log(k/k_ref)
                                                                         !! (unitless) 
    REAL(r_std)                                         :: kfact_max     !! Maximum factor for Ks decay with depth (unitless)
    REAL(r_std)                                         :: k_tmp, tmc_litter_ratio
    INTEGER(i_std), PARAMETER                           :: error_level = 3 !! Error level for consistency check
                                                                           !! Switch to 2 tu turn fatal errors into warnings
    REAL(r_std), DIMENSION (kjpindex,nslm)              :: alphavg         !! VG param a modified with depth at each node
                                                                           !! @tex $(mm^{-1})$ @endtexe
    REAL(r_std), DIMENSION (kjpindex,nslm)              :: nvg             !! VG param n modified with depth at each node
                                                                           !! (unitless)
    INTEGER(i_std)                                      :: jiref           !! To identify the mc_lins where k_lin and d_lin 
                                                                           !! need special treatment
    REAL(r_std)                                         :: nroot_tmp
    REAL(r_std)               :: zx1

!_ ================================================================================================================================

!!??Aurelien: Les 3 parametres qui suivent pourait peut-être mis dans hydrol_init?
    !
    !
    !Config Key   = CWRR_NKS_N0 
    !Config Desc  = fitted value for relation log((n-n0)/(n_ref-n0)) = nk_rel * log(k/k_ref)
    !Config Def   = 0.0
    !Config If    = HYDROL_CWRR 
    !Config Help  =
    !Config Units = [-]
    n0 = 0.0
    CALL getin_p("CWRR_NKS_N0",n0)

    !! Check parameter value (correct range)
    IF ( n0 < zero ) THEN
       CALL ipslerr_p(error_level, "hydrol_var_init.", &
            &     "Wrong parameter value for CWRR_NKS_N0.", &
            &     "This parameter should be non-negative. ", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !Config Key   = CWRR_NKS_POWER
    !Config Desc  = fitted value for relation log((n-n0)/(n_ref-n0)) = nk_rel * log(k/k_ref)
    !Config Def   = 0.0
    !Config If    = HYDROL_CWRR 
    !Config Help  =
    !Config Units = [-]
    nk_rel = 0.0
    CALL getin_p("CWRR_NKS_POWER",nk_rel)

    !! Check parameter value (correct range)
    IF ( nk_rel < zero ) THEN
       CALL ipslerr_p(error_level, "hydrol_var_init.", &
            &     "Wrong parameter value for CWRR_NKS_POWER.", &
            &     "This parameter should be non-negative. ", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !Config Key   = CWRR_AKS_A0 
    !Config Desc  = fitted value for relation log((a-a0)/(a_ref-a0)) = ak_rel * log(k/k_ref)
    !Config Def   = 0.0
    !Config If    = HYDROL_CWRR 
    !Config Help  =
    !Config Units = [1/mm]
    a0 = 0.0
    CALL getin_p("CWRR_AKS_A0",a0)

    !! Check parameter value (correct range)
    IF ( a0 < zero ) THEN
       CALL ipslerr_p(error_level, "hydrol_var_init.", &
            &     "Wrong parameter value for CWRR_AKS_A0.", &
            &     "This parameter should be non-negative. ", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !Config Key   = CWRR_AKS_POWER
    !Config Desc  = fitted value for relation log((a-a0)/(a_ref-a0)) = ak_rel * log(k/k_ref)
    !Config Def   = 0.0
    !Config If    = HYDROL_CWRR 
    !Config Help  =
    !Config Units = [-]
    ak_rel = 0.0
    CALL getin_p("CWRR_AKS_POWER",ak_rel)

    !! Check parameter value (correct range)
    IF ( nk_rel < zero ) THEN
       CALL ipslerr_p(error_level, "hydrol_var_init.", &
            &     "Wrong parameter value for CWRR_AKS_POWER.", &
            &     "This parameter should be non-negative. ", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !Config Key   = KFACT_DECAY_RATE
    !Config Desc  = Factor for Ks decay with depth
    !Config Def   = 2.0
    !Config If    = HYDROL_CWRR 
    !Config Help  =  
    !Config Units = [1/m]
    f_ks = 2.0
    CALL getin_p ("KFACT_DECAY_RATE", f_ks)

    !! Check parameter value (correct range)
    IF ( f_ks < zero ) THEN
       CALL ipslerr_p(error_level, "hydrol_var_init.", &
            &     "Wrong parameter value for KFACT_DECAY_RATE.", &
            &     "This parameter should be positive. ", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !Config Key   = KFACT_STARTING_DEPTH
    !Config Desc  = Depth for compacted value of Ks 
    !Config Def   = 0.3
    !Config If    = HYDROL_CWRR 
    !Config Help  =  
    !Config Units = [m]
    dp_comp = 0.3
    CALL getin_p ("KFACT_STARTING_DEPTH", dp_comp)

    !! Check parameter value (correct range)
    IF ( dp_comp <= zero ) THEN
       CALL ipslerr_p(error_level, "hydrol_var_init.", &
            &     "Wrong parameter value for KFACT_STARTING_DEPTH.", &
            &     "This parameter should be positive. ", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !Config Key   = KFACT_MAX
    !Config Desc  = Maximum Factor for Ks increase due to vegetation
    !Config Def   = 10.0
    !Config If    = HYDROL_CWRR 
    !Config Help  =
    !Config Units = [-]
    kfact_max = 10.0
    CALL getin_p ("KFACT_MAX", kfact_max)

    !! Check parameter value (correct range)
    IF ( kfact_max < 10. ) THEN
       CALL ipslerr_p(error_level, "hydrol_var_init.", &
            &     "Wrong parameter value for KFACT_MAX.", &
            &     "This parameter should be greater than 10. ", &
            &     "Please, check parameter value in run.def. ")
    END IF

    
    !-
    !! 1 Create local variables in mm for the vertical depths
    !!   Vertical depth variables (znh, dnh, dlh) are stored in module vertical_soil_var in m.
    DO jsl=1,nslm
       zz(jsl) = znh(jsl)*mille
       dz(jsl) = dnh(jsl)*mille
       dh(jsl) = dlh(jsl)*mille
    ENDDO

    !!xchen                                                                     
    !!DO ji = 1, kjpindex      
    !!   DO jsl = 1,nslm                                          
    !!      IF (swdown(ji) .GT. 180) THEN                                                          
    !!         mcw_mineral(njsc(ji)) = 0.5*mcw_mineral(njsc(ji))                        
    !!      ENDIF
    !!   ENDDO
    !!ENDDO
    !! consider impacts of SOC on mcs,mcw,mcf
    mcs(:)=mcs_mineral(njsc(:))
    mcw(:)=mcw_mineral(njsc(:))
    mcfc(:)=mcfc_mineral(njsc(:))

    IF (use_refSOC_hydrol) THEN
     IF (nscm .NE. 12) THEN
          CALL ipslerr_p(3, 'hydrol_var_init', 'Error: use_refSOC_hydrol=true only apply to USDA','','')
     ELSE
       zx1 = 0.
       DO ji = 1, kjpindex
          zx1 = MIN((refSOC_1d(ji)/soilc_max),1.) 
          mcs(ji) = zx1 * poros_org + (1.-zx1 )  * mcs_mineral(njsc(ji))

          !! use Van Genuchten equation to calculate mcw and mcfc
          mcw(ji) = mcr(njsc(ji)) + (mcs(ji)-mcr(njsc(ji)))/(1.+(VG_alpha(njsc(ji))*VG_psi_wp(njsc(ji)))**VG_n(njsc(ji)))**VG_m(njsc(ji))
          mcfc(ji) = mcr(njsc(ji)) + (mcs(ji)-mcr(njsc(ji)))/(1.+(VG_alpha(njsc(ji))*VG_psi_fc(njsc(ji)))**VG_n(njsc(ji)))**VG_m(njsc(ji))
       ENDDO
     ENDIF
    ENDIF

    !-
    !! 2 Compute the root density profile 
    !! Note: if ok_dynroot, nroot is modified at each time step in hydrol_soil
       DO ji=1, kjpindex
          !-
          !! The three following equations concerning nroot computation are derived from the integrals 
          !! of equations C9 to C11 of De Rosnay's (1999) PhD thesis (page 158).
          !! The occasional absence of minus sign before humcste parameter is correct.
          DO jv = 1,nvm
             nroot(ji,jv,1) = zero
             !!xchen
             !!IF (swdown(ji) .GT. 180) THEN                       
             !!   humcste_use(ji,jv)=20*humcste(jv)
             !!ELSE
             !!   humcste_use(ji,jv)=0.5*humcste(jv)
             !!ENDIF
             humcste_use(ji,jv)=humcste(jv)

             DO jsl = 2, nslm-1
                nroot(ji,jv,jsl) = (EXP(-humcste_use(ji,jv)*zz(jsl)/mille)) * &
                        & (EXP(humcste_use(ji,jv)*dz(jsl)/mille/deux) - &
                        & EXP(-humcste_use(ji,jv)*dz(jsl+1)/mille/deux))/ &
                        & (EXP(-humcste_use(ji,jv)*dz(2)/mille/deux) &
                        & -EXP(-humcste_use(ji,jv)*zz(nslm)/mille))
             ENDDO ! jsl = 2, nslm-1

             nroot(ji,jv,nslm) = (EXP(humcste_use(ji,jv)*dz(nslm)/mille/deux) -un) * &
                  & EXP(-humcste_use(ji,jv)*zz(nslm)/mille) / &
                  & (EXP(-humcste_use(ji,jv)*dz(2)/mille/deux) &
                  & -EXP(-humcste_use(ji,jv)*zz(nslm)/mille))

          !! if ok_pc: nroot is set zero below ALT, and then re-normalized to 1 along the depth
            IF ( ok_pc ) THEN
             nroot_tmp=zero
              DO jsl = 1, nslm
               IF (znh(jsl) .LT. altmax(ji,jv)) THEN
                 nroot_tmp =nroot_tmp+nroot(ji,jv,jsl)
               ELSEIF (altmax(ji,jv) .GT. zero) THEN
                 nroot(ji,jv,jsl)=zero
               ENDIF
              ENDDO ! jsl = 1, nslm
              IF (nroot_tmp .GT. zero) THEN
                 nroot(ji,jv,:)=nroot(ji,jv,:)/nroot_tmp
              ENDIF
           !ENDDO ! jsl = 1, nslm
           !IF (nroot_tmp .GT. zero) nroot(ji,jv,:)=nroot(ji,jv,:)/nroot_tmp
            ENDIF ! ok_pc

      ENDDO ! jv = 1,nvm
    ENDDO ! DO ji=1, kjpindex

    !! for check
    IF ( ANY(ABS(SUM(nroot(:,:,:),DIM=3)-un) > min_sechiba) ) THEN 
      WRITE(numout,*) 'WARNING in hydrol_var_init: nroot: vertical summation does not equal to 1'
    ENDIF

    !-
    !! 3 Compute the profile for a and n
    !-

    ! For every soil texture
    DO jsc = 1, nscm 
       DO jsl=1,nslm
          ! PhD thesis of d'Orgeval, 2006, p81, Eq. 4.38; d'Orgeval et al. 2008, Eq. 2
          ! Calibrated against Hapex-Sahel measurements
          kfact(jsl,jsc) = MIN(MAX(EXP(- f_ks * (zz(jsl)/mille - dp_comp)), un/kfact_max),un)
          ! PhD thesis of d'Orgeval, 2006, p81, Eqs. 4.39; 4.42, and Fig 4.14 
          
          nfact(jsl,jsc) = ( kfact(jsl,jsc) )**nk_rel
          afact(jsl,jsc) = ( kfact(jsl,jsc) )**ak_rel
       ENDDO
    ENDDO

    ! For every pixel
    DO jp = 1, kjpindex
       !-
       !! 4 Compute the linearized values of k, a, b and d
       !!   The effect of kfact_root on ks thus on k, a, n and d, is taken into account further in the code,
       !!   in hydrol_soil_coef.
       !-
       ! Calculate the matrix coef for Dublin model (de Rosnay, 1999; p149)
       ! piece-wise linearised hydraulic conductivity k_lin=alin * mc_lin + b_lin
       ! and diffusivity d_lin in each interval of mc, called mc_lin,
       ! between imin, for residual mcr, and imax for saturation mcs.

       ! We define 51 bounds for 50 bins of mc between mcr and mcs
       mc_lin(jp,imin)=mcr(njsc(jp))
       mc_lin(jp,imax)=mcs(jp)
       DO ji= imin+1, imax-1 ! ji=2,50
          mc_lin(jp,ji) = mcr(njsc(jp)) + (ji-imin)*(mcs(jp)-mcr(njsc(jp)))/(imax-imin)
       ENDDO

       DO jsl = 1, nslm
          ! From PhD thesis of d'Orgeval, 2006, p81, Eq. 4.42
          nvan_mod = n0 + (nvan(njsc(jp))-n0) * nfact(jsl,njsc(jp))
          avan_mod = a0 + (avan(njsc(jp))-a0) * afact(jsl,njsc(jp))
          m = un - un / nvan_mod
          ! Creation of arrays for SP-MIP output by landpoint
          nvan_mod_tab(jsl,njsc(jp)) = nvan_mod
          avan_mod_tab(jsl,njsc(jp)) = avan_mod
          ! We apply Van Genuchten equation for K(theta) based on Ks(z)=ks(jsc) * kfact(jsl,jsc)
          DO ji = imax,imin,-1 
             frac=MIN(un,(mc_lin(jp,ji)-mcr(njsc(jp)))/(mcs(jp)-mcr(njsc(jp))))
             k_lin(jp,jsl,ji) = ks(njsc(jp)) * kfact(jsl,njsc(jp)) * (frac**0.5) * ( un - ( un - frac ** (un/m)) ** m )**2
          ENDDO

          ! k_lin should not be zero, nor too small
          ! We track jiref, the bin under which mc is too small and we may get zero k_lin     
          ji=imax-1
          DO WHILE ((k_lin(jp,jsl,ji) > 1.e-32) .and. (ji>0))
             jiref=ji
             ji=ji-1
          ENDDO
          DO ji=jiref-1,imin,-1
             k_lin(jp,jsl,ji)=k_lin(jp,jsl,ji+1)/10.
          ENDDO
         
          DO ji = imin,imax-1 ! ji=1,50
             ! We deduce a_lin and b_lin based on continuity between segments k_lin = a_lin*mc-lin+b_lin
             a_lin(jp,jsl,ji) = (k_lin(jp,jsl,ji+1)-k_lin(jp,jsl,ji)) / (mc_lin(jp,ji+1)-mc_lin(jp,ji))
             b_lin(jp,jsl,ji)  = k_lin(jp,jsl,ji) - a_lin(jp,jsl,ji)*mc_lin(jp,ji)

             ! We calculate the d_lin for each mc bin, from Van Genuchten equation for D(theta)
             ! d_lin is constant and taken as the arithmetic mean between the values at the bounds of each bin 
             IF (ji.NE.imin .AND. ji.NE.imax-1) THEN
                frac=MIN(un,(mc_lin(jp,ji)-mcr(njsc(jp)))/(mcs(jp)-mcr(njsc(jp))))
                d_lin(jp,jsl,ji) =(k_lin(jp,jsl,ji) / (avan_mod*m*nvan_mod)) *  &
                     ( (frac**(-un/m))/(mc_lin(jp,ji)-mcr(njsc(jp))) ) * &
                     (  frac**(-un/m) -un ) ** (-m)
                frac=MIN(un,(mc_lin(jp,ji+1)-mcr(njsc(jp)))/(mcs(jp)-mcr(njsc(jp))))
                d_lin(jp,jsl,ji+1) =(k_lin(jp,jsl,ji+1) / (avan_mod*m*nvan_mod))*&
                     ( (frac**(-un/m))/(mc_lin(jp,ji+1)-mcr(njsc(jp))) ) * &
                     (  frac**(-un/m) -un ) ** (-m)
                d_lin(jp,jsl,ji) = undemi * (d_lin(jp,jsl,ji)+d_lin(jp,jsl,ji+1))
             ELSE IF(ji.EQ.imax-1) THEN
                d_lin(jp,jsl,ji) =(k_lin(jp,jsl,ji) / (avan_mod*m*nvan_mod)) * &
                     ( (frac**(-un/m))/(mc_lin(jp,ji)-mcr(njsc(jp))) ) *  &
                     (  frac**(-un/m) -un ) ** (-m)
             ENDIF
          ENDDO

          ! Special case for ji=imin
          d_lin(jp,jsl,imin) = d_lin(jp,jsl,imin+1)/1000.

          ! We adjust d_lin where k_lin was previously adjusted otherwise we might get non-monotonous variations
          ! We don't want d_lin = zero
          DO ji=jiref-1,imin,-1
             d_lin(jp,jsl,ji)=d_lin(jp,jsl,ji+1)/10.
          ENDDO

       ENDDO
    ENDDO
    

    ! Output of alphavg and nvg at each node for SP-MIP
    DO jsl = 1, nslm
       alphavg(:,jsl) = avan_mod_tab(jsl,njsc(:))*1000. ! from mm-1 to m-1
       nvg(:,jsl) = nvan_mod_tab(jsl,njsc(:))
    ENDDO
    CALL xios_orchidee_send_field("alphavg",alphavg) ! in m-1
    CALL xios_orchidee_send_field("nvg",nvg) ! unitless

    !! 5 Water reservoir initialisation
    !
!!$    DO jst = 1,nstm
!!$       DO ji = 1, kjpindex
!!$          mx_eau_var(ji) = mx_eau_var(ji) + soiltile(ji,jst)*&
!!$               &   zmaxh*mille*mcs(ji)
!!$       END DO
!!$    END DO
!!$    IF (check_CWRR) THEN
!!$       IF ( ANY ( ABS( mx_eau_var(:) - zmaxh*mille*mcs(:) ) > min_sechiba ) ) THEN
!!$          ji=MAXLOC ( ABS( mx_eau_var(:) - zmaxh*mille*mcs(:) ) , 1)
!!$          WRITE(numout, *) "Erreur formule simplifiÃ©e mx_eau_var ! ", mx_eau_var(ji), zmaxh*mille*mcs(ji)
!!$          WRITE(numout, *) "err = ",ABS(mx_eau_var(ji) - zmaxh*mille*mcs(ji))
!!$          STOP 1
!!$       ENDIF
!!$    ENDIF

    mx_eau_var(:) = zero
    mx_eau_var(:) = zmaxh*mille*mcs(:) 

    DO ji = 1,kjpindex 
       IF (vegtot(ji) .LE. zero) THEN
          mx_eau_var(ji) = mx_eau_nobio*zmaxh
          ! Aurelien: what does vegtot=0 mean? is it like frac_nobio=1? But if 0<frac_nobio<1 ???
       ENDIF

    END DO

    ! Compute the litter humidity, shumdiag and fry
    shumdiag_perma(:,:) = zero
    humtot(:) = zero
    tmc(:,:) = zero
!gmjc top 5 layer grassland soil moisture for grazing
    tmc_topgrass(:) = zero
!end gmjc

    ! Loop on soiltiles to compute the variables (ji,jst)
    DO jst=1,nstm 
       DO ji = 1, kjpindex
          tmcs(ji,jst)=zmaxh* mille*mcs(ji)
          tmcr(ji,jst)=zmaxh* mille*mcr(njsc(ji))
          tmcfc(ji,jst)=zmaxh* mille*mcfc(ji)
          tmcw(ji,jst)=zmaxh* mille*mcw(ji)
       ENDDO
    ENDDO
       
    ! The total soil moisture for each soiltile:
    DO jst=1,nstm
       DO ji=1,kjpindex
          tmc(ji,jst)= dz(2) * ( trois*mc(ji,1,jst)+ mc(ji,2,jst))/huit
       END DO
    ENDDO

    DO jst=1,nstm 
       DO jsl=2,nslm-1
          DO ji=1,kjpindex
             tmc(ji,jst) = tmc(ji,jst) + dz(jsl) * ( trois*mc(ji,jsl,jst) + mc(ji,jsl-1,jst))/huit &
                  & + dz(jsl+1)*(trois*mc(ji,jsl,jst) + mc(ji,jsl+1,jst))/huit
          END DO
       END DO
    ENDDO

    DO jst=1,nstm 
       DO ji=1,kjpindex
          tmc(ji,jst) = tmc(ji,jst) +  dz(nslm) * (trois * mc(ji,nslm,jst) + mc(ji,nslm-1,jst))/huit
          tmc(ji,jst) = tmc(ji,jst) + water2infilt(ji,jst)
       ENDDO
    END DO

!JG: hydrol_tmc_update should not be called in the initialization phase. Call of hydrol_tmc_update makes the model restart differenlty.    
!    ! If veget has been updated before restart (with LAND USE or DGVM),
!    ! tmc and mc must be modified with respect to humtot conservation.
!   CALL hydrol_tmc_update ( kjpindex, veget_max, soiltile, qsintveg)

    ! The litter variables:
    ! level 1
    DO jst=1,nstm 
       DO ji=1,kjpindex
          tmc_litter(ji,jst) = dz(2) * (trois*mcl(ji,1,jst)+mcl(ji,2,jst))/huit
!gmjc top 5 layer mc for grazing
          tmc_trampling(ji,jst) = dz(2) *(trois*mc(ji,1,jst)+mc(ji,2,jst))/huit
!end gmjc
          tmc_litter_wilt(ji,jst) = dz(2) * mcw(ji) / deux
          tmc_litter_res(ji,jst) = dz(2) * mcr(njsc(ji)) / deux
          tmc_litter_field(ji,jst) = dz(2) * mcfc(ji) / deux
          tmc_litter_sat(ji,jst) = dz(2) * mcs(ji) / deux
          tmc_litter_awet(ji,jst) = dz(2) * mc_awet(njsc(ji)) / deux
          tmc_litter_adry(ji,jst) = dz(2) * mc_adry(njsc(ji)) / deux
       ENDDO
    END DO
!gmjc top 5 layer mc for grazing
    ! sum from level 2 to 5
    DO jst=1,nstm
       DO jsl=2,6
          DO ji=1,kjpindex
             tmc_trampling(ji,jst) = tmc_trampling(ji,jst) + dz(jsl) * &
                  & ( trois*mc(ji,jsl,jst) + mc(ji,jsl-1,jst))/huit &
                  & + dz(jsl+1)*(trois*mc(ji,jsl,jst) + mc(ji,jsl+1,jst))/huit
          END DO
       END DO
    END DO
!end gmjc
    ! sum from level 2 to 4
    DO jst=1,nstm 
       DO jsl=2,4
          DO ji=1,kjpindex
             tmc_litter(ji,jst) = tmc_litter(ji,jst) + dz(jsl) * & 
                  & ( trois*mcl(ji,jsl,jst) + mcl(ji,jsl-1,jst))/huit &
                  & + dz(jsl+1)*(trois*mcl(ji,jsl,jst) + mcl(ji,jsl+1,jst))/huit
             tmc_litter_wilt(ji,jst) = tmc_litter_wilt(ji,jst) + &
                  &(dz(jsl)+ dz(jsl+1))*& 
                  & mcw(ji)/deux
             tmc_litter_res(ji,jst) = tmc_litter_res(ji,jst) + &
                  &(dz(jsl)+ dz(jsl+1))*& 
                  & mcr(njsc(ji))/deux
             tmc_litter_sat(ji,jst) = tmc_litter_sat(ji,jst) + &
                  &(dz(jsl)+ dz(jsl+1))* & 
                  & mcs(ji)/deux
             tmc_litter_field(ji,jst) = tmc_litter_field(ji,jst) + &
                  & (dz(jsl)+ dz(jsl+1))* & 
                  & mcfc(ji)/deux
             tmc_litter_awet(ji,jst) = tmc_litter_awet(ji,jst) + &
                  &(dz(jsl)+ dz(jsl+1))* & 
                  & mc_awet(njsc(ji))/deux
             tmc_litter_adry(ji,jst) = tmc_litter_adry(ji,jst) + &
                  & (dz(jsl)+ dz(jsl+1))* & 
                  & mc_adry(njsc(ji))/deux
          END DO
       END DO
    END DO


    DO jst=1,nstm 
       DO ji=1,kjpindex
          ! here we set that humrelv=0 in PFT1
          humrelv(ji,1,jst) = zero
       ENDDO
    END DO

!gmjc top 5 layer grassland soil moisture for grazing
    tmc_topgrass(:) = tmc_trampling(:,3)/(SUM(dz(1:6))+dz(7)/2)
!WRITE (numout,*) 'sechiba inittmc_topgrass',tmc_topgrass
!end gmjc

    ! Calculate shumdiag_perma for thermosoil
    ! Use resdist instead of soiltile because we here need to have 
    ! shumdiag_perma at the value from previous time step.
    ! Here, soilmoist is only used as a temporary variable to calculate shumdiag_perma
    ! (based on resdist=soiltile from previous timestep, but normally equal to soiltile)
    ! For consistency with hydrol_soil, we want to calculate a grid-cell average
    soilmoist(:,:) = zero
    DO jst=1,nstm
       DO ji=1,kjpindex
          soilmoist(ji,1) = soilmoist(ji,1) + resdist(ji,jst) * &
               dz(2) * ( trois*mc(ji,1,jst) + mc(ji,2,jst) )/huit
          DO jsl = 2,nslm-1
             soilmoist(ji,jsl) = soilmoist(ji,jsl) + resdist(ji,jst) * &
                  ( dz(jsl) * (trois*mc(ji,jsl,jst)+mc(ji,jsl-1,jst))/huit &
                  + dz(jsl+1) * (trois*mc(ji,jsl,jst)+mc(ji,jsl+1,jst))/huit )
          END DO
          soilmoist(ji,nslm) = soilmoist(ji,nslm) + resdist(ji,jst) * &
               dz(nslm) * (trois*mc(ji,nslm,jst) + mc(ji,nslm-1,jst))/huit
       ENDDO
    ENDDO
    DO ji=1,kjpindex
        soilmoist(ji,:) = soilmoist(ji,:) * vegtot_old(ji) ! grid cell average 
    ENDDO
    
    ! -- shumdiag_perma for restart
    !  For consistency with hydrol_soil, we want to calculate a grid-cell average
    DO jsl = 1, nslm
       DO ji=1,kjpindex        
          shumdiag_perma(ji,jsl) = soilmoist(ji,jsl) / dh(jsl)/mcs(ji)*mcs_mineral(njsc(ji))
          shumdiag_perma(ji,jsl) = MAX(MIN(shumdiag_perma(ji,jsl), un), zero) 
       ENDDO
    ENDDO
               
    ! Calculate drysoil_frac if it was not found in the restart file
    ! For simplicity, we set drysoil_frac to 0.5 in this case
    IF (ALL(drysoil_frac(:) == val_exp)) THEN
       DO ji=1,kjpindex
          drysoil_frac(ji) = 0.5
       END DO
    END IF

    profil_froz_hydro_ns(:,:,:) = 0.0
    IF (ok_freeze_cwrr) THEN
       profil_froz_hydro(:,:) = 0.0
       temp_hydro(:,:) = 0.0
    ENDIF
    
    !! Calculate the volumetric soil moisture content (mc_layh and mcl_layh) needed in 
    !! thermosoil for the thermal conductivity. 
    ! These values are only used in thermosoil_init in absence of a restart file
    mc_layh(:,:) = zero
    mcl_layh(:,:) = zero
    
    mc_layh_s = MAX(mc, min_sechiba)
    mcl_layh_s = MAX(mc, min_sechiba)
    DO jst=1,nstm
       DO jsl=1,nslm
          DO ji=1,kjpindex
            mc_layh(ji,jsl) = mc_layh(ji,jsl) + mc(ji,jsl,jst) * resdist(ji,jst)  * vegtot_old(ji)
            mcl_layh(ji,jsl) = mcl_layh(ji,jsl) + mcl(ji,jsl,jst) * resdist(ji,jst) * vegtot_old(ji)
         ENDDO
      END DO
    END DO

    IF (printlev>=3) WRITE (numout,*) ' hydrol_var_init done '

  END SUBROUTINE hydrol_var_init


!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_snow
!!
!>\BRIEF        This routine computes snow processes. 
!!
!! DESCRIPTION  :
!! - 0 initialisation
!! - 1 On vegetation
!! - 1.1 Compute snow masse
!! - 1.2 Sublimation 
!! - 1.2.1 Check that sublimation on the vegetated fraction is possible.
!! - 1.3. snow melt only if temperature positive
!! - 1.3.1 enough snow for melting or not
!! - 1.3.2 not enough snow
!! - 1.3.3 negative snow - now snow melt
!! - 1.4 Snow melts only on weight glaciers
!! - 2 On Land ice
!! - 2.1 Compute snow
!! - 2.2 Sublimation 
!! - 2.3 Snow melt only for continental ice fraction
!! - 2.3.1 If there is snow on the ice-fraction it can melt
!! - 2.4 Snow melts only on weight glaciers 
!! - 3 On other surface types - not done yet
!! - 4 computes total melt (snow and ice)
!! - 5 computes snow age on veg and ice (for albedo)
!! - 5.1 Snow age on vegetation
!! - 5.2 Snow age on ice
!! - 6 Diagnose the depth of the snow layer
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_snow

  SUBROUTINE hydrol_snow (kjpindex, precip_rain, precip_snow , temp_sol_new, soilcap,&
       & frac_nobio, totfrac_nobio, vevapsno, snow, snow_age, snow_nobio, snow_nobio_age, &
       & tot_melt, snowdepth,snowmelt)

    ! 
    ! interface description

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    ! input scalar 
    INTEGER(i_std), INTENT(in)                               :: kjpindex      !! Domain size
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: precip_rain   !! Rainfall
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: precip_snow   !! Snow precipitation
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: temp_sol_new  !! New soil temperature
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: soilcap       !! Soil capacity
    REAL(r_std), DIMENSION (kjpindex,nnobio), INTENT(in)     :: frac_nobio    !! Fraction of continental ice, lakes, ...
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: totfrac_nobio !! Total fraction of continental ice+lakes+ ...

    !! 0.2 Output variables

    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: tot_melt      !! Total melt from snow and ice  
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: snowmelt      !! Snow melt
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: snowdepth     !! Snow depth

    !! 0.3 Modified variables

    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)         :: vevapsno      !! Snow evaporation
    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)         :: snow          !! Snow mass [Kg/m^2]
    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)         :: snow_age      !! Snow age
    REAL(r_std), DIMENSION (kjpindex,nnobio), INTENT(inout)  :: snow_nobio    !! Ice water balance
    REAL(r_std), DIMENSION (kjpindex,nnobio), INTENT(inout)  :: snow_nobio_age!! Snow age on ice, lakes, ...

    !! 0.4 Local variables

    INTEGER(i_std)                               :: ji, jv
    REAL(r_std), DIMENSION (kjpindex)             :: d_age  !! Snow age change
    REAL(r_std), DIMENSION (kjpindex)             :: xx     !! temporary
    REAL(r_std)                                   :: snowmelt_tmp !! The name says it all !
    REAL(r_std)                                   :: snow_d1k !! The amount of snow that corresponds to a 1K cooling

!_ ================================================================================================================================

    !
    ! for continental points
    !

    !
    !!_0 initialisation
    !
    DO jv = 1, nnobio
       DO ji=1,kjpindex
          subsnownobio(ji,jv) = zero
       ENDDO
    ENDDO
    DO ji=1,kjpindex
       subsnowveg(ji) = zero
       snowmelt(ji) = zero
       icemelt(ji) = zero
       subsinksoil(ji) = zero
       tot_melt(ji) = zero
    ENDDO
    !
    !! 1 On vegetation
    !
    DO ji=1,kjpindex
       !
    !! 1.1 Compute snow masse
       !
       snow(ji) = snow(ji) + (un - totfrac_nobio(ji))*precip_snow(ji)
       !
       !
    !! 1.2 Sublimation 
       !      Separate between vegetated and no-veget fractions 
       !      Care has to be taken as we might have sublimation from the
       !      the frac_nobio while there is no snow on the rest of the grid.
       !
       IF ( snow(ji) > snowcri ) THEN
          subsnownobio(ji,iice) = frac_nobio(ji,iice)*vevapsno(ji)
          subsnowveg(ji) = vevapsno(ji) - subsnownobio(ji,iice)
       ELSE
          ! Correction Nathalie - Juillet 2006.
          ! On doit d'abord tester s'il existe un frac_nobio!
          ! Pour le moment je ne regarde que le iice
          IF ( frac_nobio(ji,iice) .GT. min_sechiba) THEN
             subsnownobio(ji,iice) = vevapsno(ji)
             subsnowveg(ji) = zero
          ELSE 
             subsnownobio(ji,iice) = zero
             subsnowveg(ji) = vevapsno(ji)
          ENDIF
       ENDIF
       ! here vevapsno bas been separated into a bio and nobio fractions, without changing the total
       !
       !
    !! 1.2.1 Check that sublimation on the vegetated fraction is possible.
       !
       IF (subsnowveg(ji) .GT. snow(ji)) THEN
          ! What could not be sublimated goes into subsinksoil
          IF( (un - totfrac_nobio(ji)).GT.min_sechiba) THEN
             subsinksoil (ji) = (subsnowveg(ji) - snow(ji))/ (un - totfrac_nobio(ji))
          END IF
          ! Sublimation is thus limited to what is available
          ! Then, evavpsnow is reduced, of subsinksoil
          subsnowveg(ji) = snow(ji)
          snow(ji) = zero
          vevapsno(ji) = subsnowveg(ji) + subsnownobio(ji,iice)
       ELSE
          snow(ji) = snow(ji) - subsnowveg(ji)
       ENDIF
       !
    !! 1.3. snow melt only if temperature positive
       !
       IF (temp_sol_new(ji).GT.tp_00) THEN
          !
          IF (snow(ji).GT.sneige) THEN
             !
             snowmelt(ji) = (un - frac_nobio(ji,iice))*(temp_sol_new(ji) - tp_00) * soilcap(ji) / chalfu0
             !
    !! 1.3.1 enough snow for melting or not
             !
             IF (snowmelt(ji).LT.snow(ji)) THEN
                snow(ji) = snow(ji) - snowmelt(ji)
             ELSE
                snowmelt(ji) = snow(ji)
                snow(ji) = zero
             END IF
             !
          ELSEIF (snow(ji).GE.zero) THEN
             !
    !! 1.3.2 not enough snow
             !
             snowmelt(ji) = snow(ji)
             snow(ji) = zero
          ELSE
             !
    !! 1.3.3 negative snow - now snow melt
             !
             snow(ji) = zero
             snowmelt(ji) = zero
             WRITE(numout,*) 'hydrol_snow: WARNING! snow was negative and was reset to zero. '
             !
          END IF

       ENDIF
    !! 1.4 Snow melts above a threshold
       ! Ice melt only if there is more than a given mass : maxmass_snow,
       ! But the snow cannot melt more in one time step to what corresponds to
       ! a 1K cooling. This will lead to a progressive melting of snow above
       ! maxmass_snow but it is needed as a too strong cooling can destabilise the model.
       IF ( snow(ji) .GT. maxmass_snow ) THEN
          snow_d1k = un * soilcap(ji) / chalfu0
          snowmelt(ji) = snowmelt(ji) + MIN((snow(ji) - maxmass_snow),snow_d1k)
          snow(ji) = snow(ji) - snowmelt(ji)
          IF ( printlev >= 3 ) WRITE (numout,*) "Snow was above maxmass_snow (", maxmass_snow,") and we melted ", snowmelt(ji)
       ENDIF
       
    END DO
    !
    !! 2 On Land ice
    !
    DO ji=1,kjpindex
       !
    !! 2.1 Compute snow
       !
       !!??Aurelien: pkoi mettre precip_rain en dessous? We considere liquid precipitations becomes instantly snow?  
       snow_nobio(ji,iice) = snow_nobio(ji,iice) + frac_nobio(ji,iice)*precip_snow(ji) + &
            & frac_nobio(ji,iice)*precip_rain(ji)
       !
    !! 2.2 Sublimation 
       !      Was calculated before it can give us negative snow_nobio but that is OK
       !      Once it goes below a certain values (-maxmass_snow for instance) we should kill
       !      the frac_nobio(ji,iice) !
       !
       snow_nobio(ji,iice) = snow_nobio(ji,iice) - subsnownobio(ji,iice)
       !
    !! 2.3 Snow melt only for continental ice fraction
       !
       snowmelt_tmp = zero
       IF (temp_sol_new(ji) .GT. tp_00) THEN
          !
    !! 2.3.1 If there is snow on the ice-fraction it can melt
          !
          snowmelt_tmp = frac_nobio(ji,iice)*(temp_sol_new(ji) - tp_00) * soilcap(ji) / chalfu0
          !
          IF ( snowmelt_tmp .GT. snow_nobio(ji,iice) ) THEN
             snowmelt_tmp = MAX( zero, snow_nobio(ji,iice))
          ENDIF
          snowmelt(ji) = snowmelt(ji) + snowmelt_tmp
          snow_nobio(ji,iice) = snow_nobio(ji,iice) - snowmelt_tmp
          !
       ENDIF
       !
    !! 2.4 Snow melts over a threshold
       !   Ice melt only if there is more than a given mass : maxmass_snow, 
       !   But the snow cannot melt more in one time step to what corresponds to
       !   a 1K cooling. This will lead to a progressive melting of snow above
       !   maxmass_snow but it is needed as a too strong cooling can destabilise the model.
       !
       IF ( snow_nobio(ji,iice) .GT. maxmass_snow ) THEN
          snow_d1k = un * soilcap(ji) / chalfu0
          icemelt(ji) = MIN((snow_nobio(ji,iice) - maxmass_snow),snow_d1k)
          snow_nobio(ji,iice) = snow_nobio(ji,iice) - icemelt(ji)

          IF ( printlev >= 3 ) WRITE (numout,*) "Snow was above maxmass_snow ON ICE (", maxmass_snow,") and we melted ", icemelt(ji)
       ENDIF

    END DO

    !
    !! 3 On other surface types - not done yet
    !
    IF ( nnobio .GT. 1 ) THEN
       WRITE(numout,*) 'WE HAVE',nnobio-1,' SURFACE TYPES I DO NOT KNOW'
       WRITE(numout,*) 'CANNOT TREAT SNOW ON THESE SURFACE TYPES'
       CALL ipslerr_p(3,'hydrol_snow','nnobio > 1 not allowded','Cannot treat snow on these surface types.','')
    ENDIF

    !
    !! 4 computes total melt (snow and ice)
    !
    DO ji = 1, kjpindex
       tot_melt(ji) = icemelt(ji) + snowmelt(ji)
    ENDDO

    !
    !! 5 computes snow age on veg and ice (for albedo)
    !
    DO ji = 1, kjpindex
       !
    !! 5.1 Snow age on vegetation
       !
       IF (snow(ji) .LE. zero) THEN
          snow_age(ji) = zero
       ELSE
          snow_age(ji) =(snow_age(ji) + (un - snow_age(ji)/max_snow_age) * dt_sechiba/one_day) &
               & * EXP(-precip_snow(ji) / snow_trans)
       ENDIF
       !
    !! 5.2 Snow age on ice
       !
       ! age of snow on ice: a little bit different because in cold regions, we really
       ! cannot negect the effect of cold temperatures on snow metamorphism any more.
       !
       IF (snow_nobio(ji,iice) .LE. zero) THEN
          snow_nobio_age(ji,iice) = zero
       ELSE
          !
          d_age(ji) = ( snow_nobio_age(ji,iice) + &
               &  (un - snow_nobio_age(ji,iice)/max_snow_age) * dt_sechiba/one_day ) * &
               &  EXP(-precip_snow(ji) / snow_trans) - snow_nobio_age(ji,iice)
          IF (d_age(ji) .GT. min_sechiba ) THEN
             xx(ji) = MAX( tp_00 - temp_sol_new(ji), zero )
             xx(ji) = ( xx(ji) / 7._r_std ) ** 4._r_std
             d_age(ji) = d_age(ji) / (un+xx(ji))
          ENDIF
          snow_nobio_age(ji,iice) = MAX( snow_nobio_age(ji,iice) + d_age(ji), zero )
          !
       ENDIF

    ENDDO

    !
    !! 6 Diagnose the depth of the snow layer
    !

    DO ji = 1, kjpindex
       snowdepth(ji) = snow(ji) /sn_dens
    ENDDO

    IF (printlev>=3) WRITE (numout,*) ' hydrol_snow done '

  END SUBROUTINE hydrol_snow

   
!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_canop
!!
!>\BRIEF        This routine computes canopy processes.
!!
!! DESCRIPTION  :
!! - 1 evaporation off the continents
!! - 1.1 The interception loss is take off the canopy. 
!! - 1.2 precip_rain is shared for each vegetation type
!! - 1.3 Limits the effect and sum what receives soil
!! - 1.4 swap qsintveg to the new value
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_canop

  SUBROUTINE hydrol_canop (kjpindex, precip_rain, vevapwet, veget_max, veget, qsintmax, &
       & qsintveg,precisol,tot_melt)

    ! 
    ! interface description
    !

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                               :: kjpindex    !! Domain size
    ! input fields
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: precip_rain !! Rain precipitation
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(in)        :: vevapwet    !! Interception loss
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(in)        :: veget_max   !! max fraction of vegetation type
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(in)        :: veget       !! Fraction of vegetation type 
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(in)        :: qsintmax    !! Maximum water on vegetation for interception
    REAL(r_std), DIMENSION  (kjpindex), INTENT (in)          :: tot_melt    !! Total melt

    !! 0.2 Output variables

    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(out)       :: precisol    !! Water fallen onto the ground (throughfall+Totmelt)

    !! 0.3 Modified variables

    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(inout)     :: qsintveg    !! Water on vegetation due to interception

    !! 0.4 Local variables

    INTEGER(i_std)                                           :: ji, jv
    REAL(r_std), DIMENSION (kjpindex,nvm)                    :: zqsintvegnew

!_ ================================================================================================================================

    ! boucle sur les points continentaux
    ! calcul de qsintveg au pas de temps suivant
    ! par ajout du flux interception loss
    ! calcule par enerbil en fonction
    ! des calculs faits dans diffuco
    ! calcul de ce qui tombe sur le sol
    ! avec accumulation dans precisol
    ! essayer d'harmoniser le traitement du sol nu
    ! avec celui des differents types de vegetation
    ! fait si on impose qsintmax ( ,1) = 0.0
    !
    ! loop for continental subdomain
    !
    !
    !! 1 evaporation off the continents
    !
    !! 1.1 The interception loss is take off the canopy. 
    DO jv=2,nvm
       qsintveg(:,jv) = qsintveg(:,jv) - vevapwet(:,jv)
    END DO

    !     It is raining :
    !! 1.2 precip_rain is shared for each vegetation type
    !
    qsintveg(:,1) = zero
    DO jv=2,nvm
       qsintveg(:,jv) = qsintveg(:,jv) + veget(:,jv) * ((1-throughfall_by_pft(jv))*precip_rain(:))
          ! should consider the veget_max(:,jv) to intercept rain
          ! xuhui 20151216
          !DO ji = 1,kjpindex
          !    IF (veget_max(ji,jv) .GT. zero) THEN ! otherwise, there will be no interception
          !        qsintveg(ji,jv) = qsintveg(ji,jv) + veget(ji,jv) / veget_max(ji,jv) * (1 - throughfall_by_pft(jv)) * precip_rain(ji)
          !    ENDIF
          !ENDDO
    END DO

    !
    !! 1.3 Limits the effect and sum what receives soil
    !
    precisol(:,1)=veget_max(:,1)*precip_rain(:)
    DO jv=2,nvm
       DO ji = 1, kjpindex
          zqsintvegnew(ji,jv) = MIN (qsintveg(ji,jv),qsintmax(ji,jv)) 
          precisol(ji,jv) = (veget(ji,jv)*throughfall_by_pft(jv)*precip_rain(ji)) + &
               qsintveg(ji,jv) - zqsintvegnew (ji,jv) + &
               (veget_max(ji,jv) - veget(ji,jv))*precip_rain(ji)
       ENDDO
    END DO
        
    ! Precisol is currently the same as throughfall, save it for diagnostics
    throughfall(:,:) = precisol(:,:)

    DO jv=1,nvm
       DO ji = 1, kjpindex
          IF (vegtot(ji).GT.min_sechiba) THEN
             precisol(ji,jv) = precisol(ji,jv)+tot_melt(ji)*veget_max(ji,jv)/vegtot(ji)
          ENDIF
       ENDDO
    END DO
    !   
    !
    !! 1.4 swap qsintveg to the new value
    !
    DO jv=2,nvm
       qsintveg(:,jv) = zqsintvegnew (:,jv)
    END DO

    IF (printlev>=3) WRITE (numout,*) ' hydrol_canop done '

  END SUBROUTINE hydrol_canop


!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_vegupd
!!
!>\BRIEF        Vegetation update   
!!
!! DESCRIPTION  :
!!   The vegetation cover has changed and we need to adapt the reservoir distribution 
!!   and the distribution of plants on different soil types.
!!   You may note that this occurs after evaporation and so on have been computed. It is
!!   not a problem as a new vegetation fraction will start with humrel=0 and thus will have no
!!   evaporation. If this is not the case it should have been caught above.
!!
!! - 1 Update of vegetation is it needed?
!! - 2 calculate water mass that we have to redistribute
!! - 3 put it into reservoir of plant whose surface area has grown
!! - 4 Soil tile gestion
!! - 5 update the corresponding masks
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_vegupd

  SUBROUTINE hydrol_vegupd(kjpindex, veget, veget_max, soiltile, qsintveg, frac_bare, drain_upd, runoff_upd)


    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    ! input scalar 
    INTEGER(i_std), INTENT(in)                            :: kjpindex 
    ! input fields
    REAL(r_std), DIMENSION (kjpindex, nvm), INTENT(in)    :: veget            !! New vegetation map
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)     :: veget_max        !! Max. fraction of vegetation type
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT (in)   :: soiltile         !! Fraction of each soil tile within vegtot (0-1, unitless)

    !! 0.2 Output variables
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(out)    :: frac_bare        !! Fraction(of veget_max) of bare soil
                                                                              !! in each vegetation type
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)         :: drain_upd        !! Change in drainage due to decrease in vegtot
                                                                              !! on mc [kg/m2/dt]
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)         :: runoff_upd       !! Change in runoff due to decrease in vegtot
                                                                              !! on water2infilt[kg/m2/dt]
    

    !! 0.3 Modified variables

    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout)  :: qsintveg         !! Water on old vegetation 

    !! 0.4 Local variables

    INTEGER(i_std)                                 :: ji,jv,jst

!_ ================================================================================================================================

    !! 1 If veget has been updated at last time step (with LAND USE or DGVM),
    !! tmc and mc must be modified with respect to humtot conservation.
    CALL hydrol_tmc_update ( kjpindex, veget_max, soiltile, qsintveg, drain_upd, runoff_upd)


    ! Compute the masks for veget
    
    mask_veget(:,:) = 0
    mask_soiltile(:,:) = 0
    
    DO jst=1,nstm
       DO ji = 1, kjpindex
          IF(soiltile(ji,jst) .GT. min_sechiba) THEN
             mask_soiltile(ji,jst) = 1
          ENDIF
       END DO
    ENDDO
          
    DO jv = 1, nvm
       DO ji = 1, kjpindex
          IF(veget_max(ji,jv) .GT. min_sechiba) THEN
             mask_veget(ji,jv) = 1
          ENDIF
       END DO
    END DO

    ! Compute vegetmax_soil 
    vegetmax_soil(:,:,:) = zero
    DO jv = 1, nvm
       jst = pref_soil_veg(jv)
       DO ji=1,kjpindex
          ! for veget distribution used in sechiba via humrel
          IF (mask_soiltile(ji,jst).GT.0 .AND. vegtot(ji) > min_sechiba) THEN
             vegetmax_soil(ji,jv,jst)=veget_max(ji,jv)/soiltile(ji,jst)
          ENDIF
       ENDDO
    ENDDO

    ! Calculate frac_bare (previosly done in slowproc_veget)
    DO ji =1, kjpindex
       IF( veget_max(ji,1) .GT. min_sechiba ) THEN
          frac_bare(ji,1) = un
       ELSE
          frac_bare(ji,1) = zero
       ENDIF
    ENDDO
    DO jv = 2, nvm
       DO ji =1, kjpindex
          IF( veget_max(ji,jv) .GT. min_sechiba ) THEN
             frac_bare(ji,jv) = un - veget(ji,jv)/veget_max(ji,jv)
          ELSE
             frac_bare(ji,jv) = zero
          ENDIF
       ENDDO
    ENDDO

    ! Tout dans cette routine est maintenant certainement obsolete (veget_max etant constant) en dehors des lignes 
    ! suivantes et le calcul de frac_bare:
    frac_bare_ns(:,:) = zero
    DO jst = 1, nstm
       DO jv = 1, nvm
          DO ji =1, kjpindex
             IF(vegtot(ji) .GT. min_sechiba) THEN
                frac_bare_ns(ji,jst) = frac_bare_ns(ji,jst) + vegetmax_soil(ji,jv,jst) * frac_bare(ji,jv) / vegtot(ji)
             ENDIF
          END DO
       ENDDO
    END DO
    
    IF (printlev>=3) WRITE (numout,*) ' hydrol_vegupd done '

  END SUBROUTINE hydrol_vegupd


!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_flood
!!
!>\BRIEF        This routine computes the evolution of the surface reservoir (floodplain).  
!!
!! DESCRIPTION  :
!! - 1 Take out vevapflo from the reservoir and transfer the remaining to subsinksoil
!! - 2 Compute the total flux from floodplain floodout (transfered to routing)
!! - 3 Discriminate between precip over land and over floodplain
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_flood

  SUBROUTINE hydrol_flood (kjpindex, vevapflo, flood_frac, flood_res, floodout)

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    ! input scalar 
    INTEGER(i_std), INTENT(in)                               :: kjpindex         !!
    ! input fields
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: flood_frac       !! Fraction of floodplains in grid box

    !! 0.2 Output variables

    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: floodout         !! Flux to take out from floodplains

    !! 0.3 Modified variables

    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)         :: flood_res        !! Floodplains reservoir estimate
    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)         :: vevapflo         !! Evaporation over floodplains

    !! 0.4 Local variables

    INTEGER(i_std)                                           :: ji, jv           !! Indices
    REAL(r_std), DIMENSION (kjpindex)                        :: temp             !! 

!_ ================================================================================================================================
    !- 
    !! 1 Take out vevapflo from the reservoir and transfer the remaining to subsinksoil 
    !-
    DO ji = 1,kjpindex
       temp(ji) = MIN(flood_res(ji), vevapflo(ji))
    ENDDO
    DO ji = 1,kjpindex
       flood_res(ji) = flood_res(ji) - temp(ji)
       subsinksoil(ji) = subsinksoil(ji) + vevapflo(ji) - temp(ji)
       vevapflo(ji) = temp(ji)
    ENDDO

    !- 
    !! 2 Compute the total flux from floodplain floodout (transfered to routing) 
    !-
    DO ji = 1,kjpindex
       floodout(ji) = vevapflo(ji) - flood_frac(ji) * SUM(precisol(ji,:))
    ENDDO

    !-
    !! 3 Discriminate between precip over land and over floodplain
    !-
    DO jv=1, nvm
       DO ji = 1,kjpindex
          precisol(ji,jv) = precisol(ji,jv) * (1 - flood_frac(ji))
       ENDDO
    ENDDO 

    IF (printlev>=3) WRITE (numout,*) ' hydrol_flood done'

  END SUBROUTINE hydrol_flood


!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_soil
!!
!>\BRIEF        This routine computes soil processes with CWRR scheme (Richards equation solved by finite differences).
!! Note that the water fluxes are in kg/m2/dt_sechiba.
!!
!! DESCRIPTION  :
!! 0. Initialisation, and split 2d variables to 3d variables, per soil tile
!! -- START MAIN LOOP (prognostic loop to update mc and mcl) OVER SOILTILES
!! 1. FIRSTLY, WE CHANGE MC BASED ON EXTERNAL FLUXES, ALL APPLIED AT THE SOIL SURFACE
!! 1.1 Reduces water2infilt and water2extract to their difference 
!! 1.2 To remove water2extract (including bare soilevaporation) from top layer
!! 1.3 Infiltration
!! 1.4 Reinfiltration of surface runoff : compute temporary surface water and extract from runoff
!! 2. SECONDLY, WE UPDATE MC FROM DIFFUSION, INCLUDING DRAINAGE AND ROOTSINK
!!    This will act on mcl (liquid water content) only
!! 2.1 K and D are recomputed after infiltration
!! 2.2 Set the tridiagonal matrix coefficients for the diffusion/redistribution scheme
!! 2.3 We define mcl (liquid water content) based on mc and profil_froz_hydro_ns
!! 2.4 We calculate the total SM at the beginning of the routine tridiag for water conservation check
!! 2.5 Defining where diffusion is solved : everywhere
!! 2.6 We define the system of linear equations for mcl redistribution
!! 2.7 Solves diffusion equations
!! 2.8 Computes drainage = bottom boundary condition, consistent with rhs(ji,jsl=nslm)
!! 2.9 For water conservation check during redistribution, we calculate the total liquid SM
!!     at the end of the routine tridiag, and we compare the difference with the flux...
!! 3. AFTER DIFFUSION/REDISTRIBUTION
!! 3.1 Updating mc, as all the following checks against saturation will compare mc to mcs
!! 3.2 Correct here the possible over-saturation values (subroutine hydrol_soil_smooth_over_mcs2 acts on mc)
!!     Here hydrol_soil_smooth_over_mcs2 discard all excess as ru_corr_ns, oriented to either ru_ns or dr_ns
!! 3.3 Negative runoff is reported to drainage
!! 3.4 Optional block to force saturation below zwt_force
!! 3.5 Diagnosing the effective water table depth
!! 3.6 Diagnose under_mcr to adapt water stress calculation below
!! 4. At the end of the prognostic calculations, we recompute important moisture variables
!! 4.1 Total soil moisture content (water2infilt added below)
!! 4.2 mcl is a module variable; we update it here for calculating bare soil evaporation,
!! 5. Optional check of the water balance of soil column (if check_cwrr)
!! 5.1 Computation of the vertical water fluxes
!! 5.2 Total mc conservation
!! 5.3 Total mc should not reach zero, or the tridiag solver will have problems
!! 6. SM DIAGNOSTICS FOR OTHER ROUTINES, MODULES, OR NEXT STEP
!! 6.1 Total soil moisture, soil moisture at litter levels, soil wetness, us, humrelv, vesgtressv
!! 6.2 We need to turn off evaporation when is_under_mcr
!! 6.3 Calculate the volumetric soil moisture content (mc_layh and mcl_layh) needed in thermosoil
!! 6.4 The hydraulic conductivities exported here are the ones used in the diffusion/redistribution
!! -- ENDING THE MAIN LOOP ON SOILTILES
!! 7. Summing 3d variables into 2d variables
!! 8. XIOS export of local variables, including water conservation checks
!! 9. COMPUTING EVAP_BARE_LIM_NS FOR NEXT TIME STEP, WITH A LOOP ON SOILTILES
!!    The principle is to run a dummy integration of the water redistribution scheme
!!    to check if the SM profile can sustain a potential evaporation.
!!    If not, the dummy integration is redone from the SM profile of the end of the normal integration,
!!    with a boundary condition leading to a very severe water limitation: mc(1)=mcr
!! 10. evap_bar_lim is the grid-cell scale beta
!! 11. Exit if error was found previously in this subroutine
!!
!! RECENT CHANGE(S) : 2016 by A. Ducharne
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) :
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_soil

  SUBROUTINE hydrol_soil (kjpindex, veget_max, soiltile, njsc, reinf_slope, &
       & transpir, vevapnu, vevapnu_pft, evapot, evapot_penm, runoff, drainage, &
       & returnflow, reinfiltration, irrigation, irrig_demand_ratio, &
       & tot_melt, evap_bare_lim,  shumdiag, shumdiag_perma,&
!       & tot_melt, evap_bare_lim, evap_bare_lim_pft, shumdiag, shumdiag_perma,&
       & k_litt, litterhumdiag, humrel,vegstress, drysoil_frac, irrig_fin, &
       & is_crop_soil, stempdiag,snow, &
       & snowdz, tot_bare_soil, u, v, tq_cdrag, &
       & mc_layh, mcl_layh, &
       & mc_layh_s, mcl_layh_s, &
!gmjc
       & tmc_topgrass, swdown) !!xchen add swdown
!end gmjc
    ! 
    ! interface description

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                               :: kjpindex 
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (in)       :: veget_max        !! Map of max vegetation types [-]
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)         :: njsc             !! Index of the dominant soil textural class 
                                                                                 !!   in the grid cell (1-nscm, unitless)
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT (in)      :: soiltile         !! Fraction of each soil tile within vegtot (0-1, unitless)
    LOGICAL, DIMENSION (nstm), INTENT (in)                   :: is_crop_soil     !! Whether the tile is under cropland
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(in)        :: transpir         !! Transpiration  
                                                                                 !!  @tex $(kg m^{-2} dt\_sechiba^{-1})$ @endtex
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)           :: reinf_slope      !! Fraction of surface runoff that reinfiltrates
                                                                                 !!  (unitless, [0-1])
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: returnflow       !! Water returning to the soil from the bottom
                                                                                 !!  @tex $(kg m^{-2} dt\_sechiba^{-1})$ @endtex
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: reinfiltration   !! Water returning to the top of the soil
                                                                                 !!  @tex $(kg m^{-2} dt\_sechiba^{-1})$ @endtex
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: irrigation       !! Irrigation
                                                                                 !!  @tex $(kg m^{-2} dt\_sechiba^{-1})$ @endtex
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(in)        :: irrig_demand_ratio  !! ratio of irrigation water [unitless,0-1]
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: evapot           !! Potential evaporation
                                                                                 !!  @tex $(kg m^{-2} dt\_sechiba^{-1})$ @endtex
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: evapot_penm      !! Potential evaporation "Penman" (Milly's correction)
                                                                                 !!  @tex $(kg m^{-2} dt\_sechiba^{-1})$ @endtex
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: tot_melt         !! Total melt from snow and ice
                                                                                 !!  @tex $(kg m^{-2} dt\_sechiba^{-1})$ @endtex
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in)       :: stempdiag        !! Diagnostic temp profile from thermosoil
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: snow             !! Snow mass
                                                                                 !!  @tex $(kg m^{-2})$ @endtex
    REAL(r_std), DIMENSION (kjpindex,nsnow),INTENT(in)       :: snowdz           !! Snow depth (m)
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: tot_bare_soil    !! Total evaporating bare soil fraction 
                                                                                 !!  (unitless, [0-1])
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)             :: u,v              !! Horizontal wind speed
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)             :: tq_cdrag         !! Surface drag coefficient
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: swdown           !!xchen Down-welling surface short-wave flux  

    !! 0.2 Output variables

    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: runoff           !! Surface runoff
                                                                                 !!  @tex $(kg m^{-2} dt\_sechiba^{-1})$ @endtex
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: drainage         !! Drainage
                                                                                 !!  @tex $(kg m^{-2} dt\_sechiba^{-1})$ @endtex
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: evap_bare_lim    !! Limitation factor (beta) for bare soil evaporation  
                                                                                 !! on each soil column (unitless, [0-1])
    REAL(r_std), DIMENSION (kjpindex,nslm), INTENT (out)     :: shumdiag         !! Relative soil moisture in each diag soil layer 
                                                                                 !! with respect to (mcfc-mcw) (unitless, [0-1])
    REAL(r_std), DIMENSION (kjpindex,nslm), INTENT (out)     :: shumdiag_perma   !! Percent of porosity filled with water (mc/mcs)
                                                                                 !! in each diag soil layer (for the thermal computations)
                                                                                 !! (unitless, [0-1])
    REAL(r_std), DIMENSION (kjpindex), INTENT (out)          :: k_litt           !! Litter approximated hydraulic conductivity
                                                                                 !!  @tex $(mm d^{-1})$ @endtex
    REAL(r_std), DIMENSION (kjpindex), INTENT (out)          :: litterhumdiag    !! Mean of soil_wet_litter across soil tiles
                                                                                 !! (unitless, [0-1])
    REAL(r_std), DIMENSION (kjpindex, nvm), INTENT(out)      :: vegstress        !! Veg. moisture stress (only for vegetation 
                                                                                 !! growth) (unitless, [0-1])
    REAL(r_std), DIMENSION (kjpindex), INTENT (out)          :: drysoil_frac     !! Function of the litter humidity
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(out)       :: irrig_fin        !! final application of irrigation [mm]
    REAL(r_std), DIMENSION (kjpindex,nslm), INTENT (out)     :: mc_layh          !! Volumetric water content (liquid + ice) for each soil layer
                                                                                 !! averaged over the mesh (for thermosoil)
                                                                                 !!  @tex $(m^{3} m^{-3})$ @endtex
    REAL(r_std), DIMENSION (kjpindex,nslm), INTENT (out)     :: mcl_layh         !! Volumetric liquid water content for each soil layer 
                                                                                 !! averaged over the mesh (for thermosoil)
                                                                                 !!  @tex $(m^{3} m^{-3})$ @endtex 
    REAL(r_std), DIMENSION (kjpindex,nslm,nstm), INTENT (out)  :: mc_layh_s          !! Volumetric soil moisture content for each layer in hydrol(liquid + ice) [m3/m3]
    REAL(r_std), DIMENSION (kjpindex,nslm,nstm), INTENT (out)  :: mcl_layh_s         !! Volumetric soil moisture content for each layer in hydrol(liquid) [m3/m3]
!gmjc
   REAL(r_std), DIMENSION (kjpindex), INTENT (out)          :: tmc_topgrass
!end gmjc
    !! 0.3 Modified variables

    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)         :: vevapnu          !! Bare soil evaporation
                                                                                 !!  @tex $(kg m^{-2} dt\_sechiba^{-1})$ @endtex
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(inout)     :: vevapnu_pft          !! Bare soil evaporation
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (inout)    :: humrel           !! Relative humidity (0-1, dimensionless)

    !! 0.4 Local variables

    INTEGER(i_std)                                 :: ji, jv, jsl, jst           !! Indices
    REAL(r_std), PARAMETER                         :: frac_mcs = 0.66            !! Temporary depth
    REAL(r_std)                                    :: tot_irrig_frac            !! temporary sum of irrigated fraction of vegetation
    REAL(r_std), DIMENSION(kjpindex)               :: temp                       !! Temporary value for fluxes
    REAL(r_std), DIMENSION(kjpindex)               :: tmcold                     !! Total SM at beginning of hydrol_soil (kg/m2)
    REAL(r_std), DIMENSION(kjpindex)               :: tmcint                     !! Ancillary total SM (kg/m2)
    REAL(r_std), DIMENSION(kjpindex,nslm)          :: mcint                      !! To save mc values for future use
    REAL(r_std), DIMENSION(kjpindex,nslm)          :: mclint                     !! To save mcl values for future use
    LOGICAL, DIMENSION(kjpindex,nstm)              :: is_under_mcr               !! Identifies under residual soil moisture points
    LOGICAL, DIMENSION(kjpindex)                   :: is_over_mcs                !! Identifies over saturated soil moisture points 
    REAL(r_std), DIMENSION(kjpindex)               :: deltahum,diff              !!
    LOGICAL(r_std), DIMENSION(kjpindex)            :: test                       !!
    REAL(r_std), DIMENSION(kjpindex)               :: water2extract              !! Water flux to be extracted at the soil surface
                                                                                 !!  @tex $(kg m^{-2} dt\_sechiba^{-1})$ @endtex
    REAL(r_std), DIMENSION(kjpindex)               :: reinfiltration_soil        !! Water from the routing back to the top of the 
                                                                                 !! soil @tex $(kg m^{-2} dt\_sechiba^{-1})$ @endtex
    REAL(r_std), DIMENSION(kjpindex, nstm)         :: irrigation_soil            !! Water from irrigation returning to soil moisture 
                                                                                 !!  @tex $(kg m^{-2} dt\_sechiba^{-1})$ @endtex
    REAL(r_std), DIMENSION(kjpindex)               :: flux_infilt                !! Water to infiltrate
                                                                                 !!  @tex $(kg m^{-2})$ @endtex
    REAL(r_std), DIMENSION(kjpindex)               :: flux_bottom                !! Flux at bottom of the soil column
                                                                                 !!  @tex $(kg m^{-2})$ @endtex
    REAL(r_std), DIMENSION(kjpindex)               :: flux_top                   !! Flux at top of the soil column (for bare soil evap)
                                                                                 !!  @tex $(kg m^{-2})$ @endtex
    REAL(r_std), DIMENSION (kjpindex,nstm)         :: qinfilt_ns                 !! Effective infiltration flux per soil tile
                                                                                 !!  @tex $(kg m^{-2})$ @endtex   
    REAL(r_std), DIMENSION (kjpindex)              :: qinfilt                    !! Effective infiltration flux  
                                                                                 !!  @tex $(kg m^{-2})$ @endtex 
    REAL(r_std), DIMENSION (kjpindex,nstm)         :: ru_infilt_ns               !! Surface runoff from hydrol_soil_infilt per soil tile
                                                                                 !!  @tex $(kg m^{-2})$ @endtex   
    REAL(r_std), DIMENSION (kjpindex)              :: ru_infilt                  !! Surface runoff from hydrol_soil_infilt 
                                                                                 !!  @tex $(kg m^{-2})$ @endtex 
    REAL(r_std), DIMENSION (kjpindex,nstm)         :: ru_corr_ns                 !! Surface runoff produced to correct excess per soil tile
                                                                                 !!  @tex $(kg m^{-2})$ @endtex
    REAL(r_std), DIMENSION (kjpindex)              :: ru_corr                    !! Surface runoff produced to correct excess 
                                                                                 !!  @tex $(kg m^{-2} dt\_sechiba^{-1})$ @endtex  
    REAL(r_std), DIMENSION (kjpindex,nstm)         :: ru_corr2_ns                !! Correction of negative surface runoff per soil tile
                                                                                 !!  @tex $(kg m^{-2})$ @endtex
    REAL(r_std), DIMENSION (kjpindex)              :: ru_corr2                   !! Correction of negative surface runoff
                                                                                 !!  @tex $(kg m^{-2})$ @endtex 
    REAL(r_std), DIMENSION (kjpindex,nstm)         :: dr_corr_ns                 !! Drainage produced to correct excess
                                                                                 !!  @tex $(kg m^{-2})$ @endtex
    REAL(r_std), DIMENSION (kjpindex,nstm)         :: dr_corrnum_ns              !! Drainage produced to correct numerical errors in tridiag
                                                                                 !!  @tex $(kg m^{-2})$ @endtex    
    REAL(r_std), DIMENSION (kjpindex)              :: dr_corr                    !! Drainage produced to correct excess 
                                                                                 !!  @tex $(kg m^{-2} dt\_sechiba^{-1})$ @endtex
    REAL(r_std), DIMENSION (kjpindex)              :: dr_corrnum                 !! Drainage produced to correct numerical errors in tridiag
                                                                                 !!  @tex $(kg m^{-2} dt\_sechiba^{-1})$ @endtex
    REAL(r_std), DIMENSION (kjpindex,nslm)         :: dmc                        !! Delta mc when forcing saturation (zwt_force) 
                                                                                 !!  @tex $(m^{3} m^{-3})$ @endtex 
    REAL(r_std), DIMENSION (kjpindex,nstm)         :: dr_force_ns                !! Delta drainage when forcing saturation (zwt_force) 
                                                                                 !!  per soil tile  @tex $(kg m^{-2})$ @endtex
    REAL(r_std), DIMENSION (kjpindex)              :: dr_force                   !! Delta drainage when forcing saturation (zwt_force)
                                                                                 !!  @tex $(kg m^{-2})$ @endtex  
    REAL(r_std), DIMENSION (kjpindex,nstm)         :: wtd_ns                     !! Effective water table depth (m)
    REAL(r_std), DIMENSION (kjpindex)              :: wtd                        !! Mean water table depth in the grid-cell (m)
    LOGICAL                                        :: error=.FALSE.              !! If true, exit in the end of subroutine

    ! For the calculation of soil_wet_ns and us/humrel/vegstress 
    REAL(r_std), DIMENSION (kjpindex,nslm)         :: sm                         !! Soil moisture of each layer (liquid phase)
                                                                                 !!  @tex $(kg m^{-2})$ @endtex
    REAL(r_std), DIMENSION (kjpindex,nslm)         :: smt                        !! Soil moisture of each layer (liquid+solid phase)
                                                                                 !!  @tex $(kg m^{-2})$ @endtex
    REAL(r_std), DIMENSION (kjpindex,nslm)         :: smw                        !! Soil moisture of each layer at wilting point
                                                                                 !!  @tex $(kg m^{-2})$ @endtex 
    REAL(r_std), DIMENSION (kjpindex,nslm)         :: smf                        !! Soil moisture of each layer at field capacity
                                                                                 !!  @tex $(kg m^{-2})$ @endtex   
    REAL(r_std), DIMENSION (kjpindex,nslm)         :: sms                        !! Soil moisture of each layer at saturation
                                                                                 !!  @tex $(kg m^{-2})$ @endtex 
    REAL(r_std), DIMENSION (kjpindex,nslm)         :: sms_tmp
    REAL(r_std), DIMENSION (kjpindex,nslm)         :: smf_tmp
    REAL(r_std), DIMENSION (kjpindex,nslm)         :: smw_tmp
    REAL(r_std), DIMENSION (kjpindex,nslm)         :: sm_nostress                !! Soil moisture of each layer at which us reaches 1
                                                                                 !!  @tex $(kg m^{-2})$ @endtex 
    ! For water conservation checks (in mm/dtstep unless otherwise mentioned)
    REAL(r_std), DIMENSION (kjpindex,nstm)         :: check_infilt_ns             !! Water conservation diagnostic at routine scale
    REAL(r_std), DIMENSION (kjpindex,nstm)         :: check1_ns                   !! Water conservation diagnostic at routine scale
    REAL(r_std), DIMENSION (kjpindex,nstm)         :: check_tr_ns                 !! Water conservation diagnostic at routine scale
    REAL(r_std), DIMENSION (kjpindex,nstm)         :: check_over_ns               !! Water conservation diagnostic at routine scale
    REAL(r_std), DIMENSION (kjpindex,nstm)         :: check_under_ns              !! Water conservation diagnostic at routine scale
    REAL(r_std), DIMENSION(kjpindex)               :: tmci                        !! Total soil moisture at beginning of routine (kg/m2)
    REAL(r_std), DIMENSION(kjpindex)               :: tmcf                        !! Total soil moisture at end of routine (kg/m2)
    REAL(r_std), DIMENSION(kjpindex)               :: diag_tr                     !! Transpiration flux
    REAL(r_std), DIMENSION (kjpindex)              :: check_infilt                !! Water conservation diagnostic at routine scale 
    REAL(r_std), DIMENSION (kjpindex)              :: check1                      !! Water conservation diagnostic at routine scale 
    REAL(r_std), DIMENSION (kjpindex)              :: check_tr                    !! Water conservation diagnostic at routine scale 
    REAL(r_std), DIMENSION (kjpindex)              :: check_over                  !! Water conservation diagnostic at routine scale 
    REAL(r_std), DIMENSION (kjpindex)              :: check_under                 !! Water conservation diagnostic at routine scale

    ! Variables for calculation of a soil resistance, option do_rsoil (following the formulation of Sellers et al 1992, implemented in Oleson et al. 2008)
    REAL(r_std)                     		   :: speed			 !! magnitude of wind speed required for Aerodynamic resistance 
    REAL(r_std)                                    :: ra			 !! diagnosed aerodynamic resistance
    REAL(r_std), DIMENSION(kjpindex)               :: mc_rel                     !! first layer relative soil moisture, required for rsoil
    REAL(r_std), DIMENSION(kjpindex)               :: evap_soil                  !! soil evaporation from Oleson et al 2008
    REAL(r_std), DIMENSION(kjpindex,nstm)          :: r_soil_ns	                 !! soil resistance from Oleson et al 2008
    REAL(r_std), DIMENSION(kjpindex)               :: r_soil	                 !! soil resistance from Oleson et al 2008
    REAL(r_std), DIMENSION(kjpindex)               :: tmcs_litter                !! Saturated soil moisture in the 4 "litter" soil layers
    REAL(r_std), DIMENSION(nslm)                   :: nroot_tmp                  !! Temporary variable to calculate the nroot

    ! For CMIP6 and SP-MIP : ksat and matric pressure head psi(theta)
    REAL(r_std)                                    :: mc_ratio, mvg, avg
    REAL(r_std)                                    :: psi                        !! Matric head (per soil layer and soil tile) [mm=kg/m2]
    REAL(r_std), DIMENSION (kjpindex,nslm)         :: psi_moy                    !! Mean matric head per soil layer [mm=kg/m2]  
    REAL(r_std), DIMENSION (kjpindex,nslm)         :: ksat                       !! Saturated hydraulic conductivity at each node (mm/d)  

    ! Temp vars
    REAL(r_std), DIMENSION(kjpindex)               :: pcent_gc                   !! Rearrange from pcent(njsc(kjpindex)) to pcent_gc(kjpindex)

!_ ================================================================================================================================

    !! 0.1 Arrays with DIMENSION(kjpindex)
    
    irrigation_soil(:, :) = zero
    qflux(:,:,:) = zero
    mc_layh(:,:) = zero ! for thermosoil
    mcl_layh(:,:) = zero ! for thermosoil
    kk(:,:,:) = zero
    kk_moy(:,:) = zero
    undermcr(:) = zero ! needs to be initialized outside from jst loop
    ksat(:,:) = zero
    psi_moy(:,:) = zero
    irrig_fin(:,:) = zero
    r_soil_ns = zero

    IF (ok_freeze_cwrr) THEN
       
       ! 0.1 Calculate the temperature and fozen fraction at the hydrological levels
       
       ! AD16*** This subroutine could probably be simplified massively given
       ! that hydro and T share the same vertical discretization 
       ! Here stempdiag is in from thermosoil and temp_hydro is out
       CALL hydrol_calculate_temp_hydro(kjpindex, stempdiag, snow,snowdz)
       
       ! Calculates profil_froz_hydro_ns as a function of temp_hydro, and mc if ok_thermodynamical_freezing
       ! These values will be kept till the end of the prognostic loop
       DO jst=1,nstm
          CALL hydrol_soil_froz(kjpindex,jst,njsc)
       ENDDO

    ELSE
  
       profil_froz_hydro_ns(:,:,:) = zero
             
    ENDIF
    
    !! 0.2 Split 2d variables to 3d variables, per soil tile
    !  Here, the evaporative fluxes are distributed over the soiltiles as a function of the 
    !    corresponding control factors; they are normalized to vegtot
    !  At step 7, the reverse transformation is used for the fluxes produced in hydrol_soil
    !    flux_cell(ji)=sum(flux_ns(ji,jst)*soiltile(ji,jst)*vegtot(ji))
    CALL hydrol_split_soil (kjpindex, veget_max, soiltile, vevapnu, vevapnu_pft, transpir, humrel, evap_bare_lim, tot_bare_soil)

    !! 0.3 Common variables related to routing, with all return flow applied to the soil surface
    ! The fluxes coming from the routing are uniformly splitted into the soiltiles,
    !    but are normalized to vegtot like the above fluxes: 
    !            flux_ns(ji,jst)=flux_cell(ji)/vegtot(ji)
    ! It is the case for : irrigation_soil(ji) and reinfiltration_soil(ji) cf below
    ! It is also the case for subsinksoil(ji), which is divided by (1-tot_frac_nobio) at creation in hydrol_snow
    ! AD16*** The transformation in 0.2 and 0.3 is likely to induce conservation problems
    !         when tot_frac_nobio NE 0, since sum(soiltile) NE vegtot in this case
    
    DO ji=1,kjpindex
       IF(vegtot(ji).GT.min_sechiba) THEN
          reinfiltration_soil(ji) = (returnflow(ji) + reinfiltration(ji))/vegtot(ji)
!          irrigation_soil(ji) = irrigation(ji)/vegtot(ji)
          ! only crop soil(nstm>=4) will be irrigated...
          tot_irrig_frac = zero
          DO jv=2,nvm
              IF ( ok_LAIdev(jv) .AND. ( irrig_demand_ratio(ji,jv) .GT. zero )) THEN
!                  IF (pref_soil_veg(jv) .LT. 4) THEN 
                  IF (.NOT. is_crop_soil(pref_soil_veg(jv)) ) THEN 
                      ! the demanded irrigation is not in the crop soil column
                      WRITE(numout,*) "pft ", jv, " pref_soil_veg(jv) ", pref_soil_veg(jv), "is not crop soil"
                      STOP "hydrol irrig"
                  ENDIF
                  tot_irrig_frac = tot_irrig_frac + veget_max(ji,jv)
              ENDIF
          ENDDO          
!          IF (tot_irrig_frac .GT. zero) THEN
!              irrigation_soil(ji) = irrigation(ji)/tot_irrig_frac
!          ELSE
!              irrigation_soil(ji) = zero
!          ENDIF
          irrig_fin(:,1) = zero
          DO jv = 2,nvm
              ! old: since we only have one crop soil, we pour all irrigation water into
              ! this soil column
              !!!IF ( irrig_demand_ratio(ji,jv) .GT. zero ) THEN
              !!!    irrig_fin(ji,jv) = irrigation_soil(ji)
              !!!ENDIF
              ! now we have several crop soil tile
              ! we irrigate separately considering soil tiles
              ! assuming each crop has its own tile respectively
              IF ((irrig_demand_ratio(ji,jv) .GT. zero) .AND. (pref_soil_veg(jv) .GE. 4)) THEN
                  irrig_fin(ji,jv) = irrigation(ji) * irrig_demand_ratio(ji,jv) / veget_max(ji,jv)
                  irrigation_soil(ji,pref_soil_veg(jv)) = irrigation_soil(ji,pref_soil_veg(jv)) + &
                             & irrigation(ji) * irrig_demand_ratio(ji,jv) / soiltile(ji,pref_soil_veg(jv))
              ENDIF
          ENDDO
       ELSE
          reinfiltration_soil(ji) = zero
       ENDIF
    ENDDO     

!!    WRITE(numout,*) "irrig xuhui: irrigation_soil(1,:): ", irrigation_soil(1,:) 
    
    !! -- START MAIN LOOP (prognostic loop to update mc and mcl) OVER SOILTILES
    !!    The called subroutines work on arrays with DIMENSION(kjpindex),
    !!    recursively used for each soiltile jst
    
    DO jst = 1,nstm

       is_under_mcr(:,jst) = .FALSE.
       is_over_mcs(:) = .FALSE.
       
       !! 0.4. Keep initial values for future check-up
       
       ! Total moisture content (including water2infilt) is saved for balance checks at the end
       ! In hydrol_tmc_update, tmc is increased by water2infilt(ji,jst), but mc is not modified !
       tmcold(:) = tmc(:,jst)
       
       ! The value of mc is kept in mcint (nstm dimension removed), in case needed for water balance checks 
       DO jsl = 1, nslm
          DO ji = 1, kjpindex
             mcint(ji,jsl) = mask_soiltile(ji,jst) * mc(ji,jsl,jst)
          ENDDO
       ENDDO
       !
       ! Initial total moisture content : tmcint does not include water2infilt, contrarily to tmcold
       DO ji = 1, kjpindex
          tmcint(ji) = dz(2) * ( trois*mcint(ji,1) + mcint(ji,2) )/huit 
       ENDDO
       DO jsl = 2,nslm-1
          DO ji = 1, kjpindex
             tmcint(ji) = tmcint(ji) + dz(jsl) &
                  & * (trois*mcint(ji,jsl)+mcint(ji,jsl-1))/huit &
                  & + dz(jsl+1) * (trois*mcint(ji,jsl)+mcint(ji,jsl+1))/huit
          ENDDO
       ENDDO
       DO ji = 1, kjpindex
          tmcint(ji) = tmcint(ji) + dz(nslm) &
               & * (trois * mcint(ji,nslm) + mcint(ji,nslm-1))/huit
       ENDDO

       !! 1. FIRSTLY, WE CHANGE MC BASED ON EXTERNAL FLUXES, ALL APPLIED AT THE SOIL SURFACE
       !!   Input = water2infilt(ji,jst) + irrigation_soil(ji) + reinfiltration_soil(ji) + precisol_ns(ji,jst)
       !!      - negative evaporation fluxes (MIN(ae_ns(ji,jst),zero)+ MIN(subsinksoil(ji),zero))
       !!   Output = MAX(ae_ns(ji,jst),zero) + subsinksoil(ji) = positive evaporation flux = water2extract
       ! In practice, negative subsinksoil(ji) is not possible 

       !! 1.1 Reduces water2infilt and water2extract to their difference 

       ! Compares water2infilt and water2extract to keep only difference
       ! Here, temp is used as a temporary variable to store the min of water to infiltrate vs evaporate
       DO ji = 1, kjpindex
          IF ( is_crop_soil(jst) ) THEN ! crop soil
              temp(ji) = MIN(water2infilt(ji,jst) + irrigation_soil(ji, jst) + reinfiltration_soil(ji) &
                         - MIN(ae_ns(ji,jst),zero) - MIN(subsinksoil(ji),zero) + precisol_ns(ji,jst), &
                           MAX(ae_ns(ji,jst),zero) + MAX(subsinksoil(ji),zero) )
          ELSE
              temp(ji) = MIN(water2infilt(ji,jst) + reinfiltration_soil(ji) &
                         - MIN(ae_ns(ji,jst),zero) - MIN(subsinksoil(ji),zero) + precisol_ns(ji,jst), &
                           MAX(ae_ns(ji,jst),zero) + MAX(subsinksoil(ji),zero) )
          ENDIF
       ENDDO

       ! The water to infiltrate at the soil surface is either 0, or the difference to what has to be evaporated
       !   - the initial water2infilt (right hand side) results from qsintveg changes with vegetation updates
       !   - irrigation_soil is the input flux to the soil surface from irrigation 
       !   - reinfiltration_soil is the input flux to the soil surface from routing 'including returnflow)
       !   - eventually, water2infilt holds all fluxes to the soil surface except precisol (reduced by water2extract)
       DO ji = 1, kjpindex
          IF ( is_crop_soil(jst) ) THEN ! crop soil
               water2infilt(ji,jst) = water2infilt(ji,jst) + irrigation_soil(ji,jst) + reinfiltration_soil(ji) &
                    - MIN(ae_ns(ji,jst),zero) - MIN(subsinksoil(ji),zero) + precisol_ns(ji,jst) &
                    - temp(ji) 
          ELSE
                water2infilt(ji,jst) = water2infilt(ji,jst) + reinfiltration_soil(ji) &
                    - MIN(ae_ns(ji,jst),zero) - MIN(subsinksoil(ji),zero) + precisol_ns(ji,jst) &
                    - temp(ji) 
          ENDIF
       ENDDO       
!!       WRITE(numout,*) "irrig xuhui: water2infilt(1,4): ", water2infilt(1,4)
             
       ! The water to evaporate from the sol surface is either the difference to what has to be infiltrated, or 0
       !   - subsinksoil is the residual from sublimation is the snowpack is not sufficient
       !   - how are the negative values of ae_ns taken into account ??? 
       DO ji = 1, kjpindex
          water2extract(ji) =  MAX(ae_ns(ji,jst),zero) + MAX(subsinksoil(ji),zero) - temp(ji) 
       ENDDO

       ! Here we acknowledge that subsinksoil is part of ae_ns, but ae_ns is not used further
       ae_ns(:,jst) = ae_ns(:,jst) + subsinksoil(:)  

       !! 1.2 To remove water2extract (including bare soil) from top layer
       flux_top(:) = water2extract(:)

       !! 1.3 Infiltration

       !! Definition of flux_infilt
       DO ji = 1, kjpindex
          ! Initialise the flux to be infiltrated  
          flux_infilt(ji) = water2infilt(ji,jst) 
       ENDDO
       
       !! K and D are computed for the profile of mc before infiltration
       !! They depend on the fraction of soil ice, given by profil_froz_hydro_ns
       CALL hydrol_soil_coef(kjpindex,jst,njsc)


       !! Infiltration and surface runoff are computed
       !! Infiltration stems from comparing liquid water2infilt to initial total mc (liquid+ice)
       !! The conductivity comes from hydrol_soil_coef and relates to the liquid phase only
       !  This seems consistent with ok_freeze
       CALL hydrol_soil_infilt(kjpindex, jst, njsc, flux_infilt, qinfilt_ns, ru_infilt_ns, &
            check_infilt_ns)
       ru_ns(:,jst) = ru_infilt_ns(:,jst) 


       !! 1.4 Reinfiltration of surface runoff : compute temporary surface water and extract from runoff
       ! Evrything here is liquid
       ! RK: water2infilt is both a volume for future reinfiltration (in mm) and a correction term for surface runoff (in mm/dt_sechiba) 
       IF ( .NOT. doponds ) THEN ! this is the general case...
          DO ji = 1, kjpindex
             water2infilt(ji,jst) = reinf_slope(ji) * ru_ns(ji,jst)
          ENDDO
       ELSE
          DO ji = 1, kjpindex           
             water2infilt(ji,jst) = zero
          ENDDO
       ENDIF
       !
       DO ji = 1, kjpindex           
          ru_ns(ji,jst) = ru_ns(ji,jst) - water2infilt(ji,jst)
       END DO

       !! 2. SECONDLY, WE UPDATE MC FROM DIFFUSION, INCLUDING DRAINAGE AND ROOTSINK
       !!    This will act on mcl only
       
       !! 2.1 K and D are recomputed after infiltration
       !! They depend on the fraction of soil ice, still given by profil_froz_hydro_ns 
       CALL hydrol_soil_coef(kjpindex,jst,njsc)
 
       !! 2.2 Set the tridiagonal matrix coefficients for the diffusion/redistribution scheme
       !! This process will further act on mcl only, based on a, b, d from hydrol_soil_coef
       CALL hydrol_soil_setup(kjpindex,jst)

       !! 2.3 We define mcl (liquid water content) based on mc and profil_froz_hydro_ns 
       DO jsl = 1, nslm
          DO ji =1, kjpindex
             mcl(ji,jsl,jst)= MIN( mc(ji,jsl,jst), mcr(njsc(ji)) + &
                  (un-profil_froz_hydro_ns(ji,jsl,jst))*(mc(ji,jsl,jst)-mcr(njsc(ji))) )
             ! we always have mcl<=mc
             ! if mc>mcr, then mcl>mcr; if mc=mcr,mcl=mcr; if mc<mcr, then mcl<mcr
             ! if profil_froz_hydro_ns=0 (including NOT ok_freeze_cwrr) we keep mcl=mc
          ENDDO
       ENDDO

       ! The value of mcl is kept in mclint (nstm dimension removed), used in the flux computation after diffusion
       DO jsl = 1, nslm
          DO ji = 1, kjpindex
             mclint(ji,jsl) = mask_soiltile(ji,jst) * mcl(ji,jsl,jst)
          ENDDO
       ENDDO

       !! 2.3bis Diagnostic of the matric potential used for redistribution by Richards/tridiag (in m)
       !  We use VG relationship giving psi as a function of mc (mcl in our case)
       !  With patches against numerical pbs when (mc_ratio - un) becomes very slightly negative (gives NaN) 
       !  or if psi become too strongly negative (pbs with xios output)
       DO jsl=1, nslm
          DO ji = 1, kjpindex
             IF (soiltile(ji,jst) .GT. zero) THEN
                mvg = un - un / nvan_mod_tab(jsl,njsc(ji))
                avg = avan_mod_tab(jsl,njsc(ji))*1000. ! to convert in m-1
                mc_ratio = MAX( 10.**(-14*mvg), (mcl(ji,jsl,jst) - mcr(njsc(ji)))/(mcs(ji) - mcr(njsc(ji))) )**(-un/mvg)
                psi = - MAX(zero,(mc_ratio - un))**(un/nvan_mod_tab(jsl,njsc(ji))) / avg ! in m 
                psi_moy(ji,jsl) = psi_moy(ji,jsl) + soiltile(ji,jst) * psi ! average across soil tiles
             ENDIF
          ENDDO
       ENDDO

       !! 2.4 We calculate the total SM at the beginning of the routine tridiag for water conservation check
       !  (on mcl only, since the diffusion only modifies mcl)
       tmci(:) = dz(2) * ( trois*mcl(:,1,jst) + mcl(:,2,jst) )/huit
       DO jsl = 2,nslm-1
          tmci(:) = tmci(:) + dz(jsl) * (trois*mcl(:,jsl,jst)+mcl(:,jsl-1,jst))/huit &
               + dz(jsl+1) * (trois*mcl(:,jsl,jst)+mcl(:,jsl+1,jst))/huit
       ENDDO
       tmci(:) = tmci(:) + dz(nslm) * (trois*mcl(:,nslm,jst) + mcl(:,nslm-1,jst))/huit

       IF (ok_freeze_cwrr) THEN
          CALL hydrol_soil_coef(kjpindex,jst,njsc)
          DO ji =1, kjpindex
             DO jsl = 1, nslm
                mcl(ji,jsl,jst)= MIN(mc(ji,jsl,jst),mcr(njsc(ji))+(1-profil_froz_hydro_ns(ji,jsl,jst))*(mc(ji,jsl,jst)-mcr(njsc(ji))))
             ENDDO
          ENDDO
       ELSE
          mcl(:,:,jst)=mc(:,:,jst)
       ENDIF
       !! 2.5 Defining where diffusion is solved : everywhere
       !! Since mc>mcs is not possible after infiltration, and we accept that mc<mcr
       !! (corrected later by shutting off all evaporative fluxes in this case)
       !  Nothing is done if resolv=F
       resolv(:) = (mask_soiltile(:,jst) .GT. 0)

       !! 2.6 We define the system of linear equations for mcl redistribution,
       !! based on the matrix coefficients from hydrol_soil_setup
       !! following the PhD thesis of de Rosnay (1999), p155-157
       !! The bare soil evaporation (subtracted from infiltration) is used directly as flux_top
       ! rhs for right-hand side term; fp for f'; gp for g'; ep for e'; with flux=0 !
       
       !- First layer
       DO ji = 1, kjpindex
          tmat(ji,1,1) = zero
          tmat(ji,1,2) = f(ji,1)
          tmat(ji,1,3) = g1(ji,1)
          rhs(ji,1)    = fp(ji,1) * mcl(ji,1,jst) + gp(ji,1)*mcl(ji,2,jst) &
               &  - flux_top(ji) - (b(ji,1)+b(ji,2))/deux *(dt_sechiba/one_day) - rootsink(ji,1,jst)
       ENDDO
       !- soil body
       DO jsl=2, nslm-1
          DO ji = 1, kjpindex
             tmat(ji,jsl,1) = e(ji,jsl)
             tmat(ji,jsl,2) = f(ji,jsl)
             tmat(ji,jsl,3) = g1(ji,jsl)
             rhs(ji,jsl) = ep(ji,jsl)*mcl(ji,jsl-1,jst) + fp(ji,jsl)*mcl(ji,jsl,jst) &
                  & +  gp(ji,jsl) * mcl(ji,jsl+1,jst) & 
                  & + (b(ji,jsl-1) - b(ji,jsl+1)) * (dt_sechiba/one_day) / deux & 
                  & - rootsink(ji,jsl,jst) 
          ENDDO
       ENDDO       
       !- Last layer, including drainage
       DO ji = 1, kjpindex
          jsl=nslm
          tmat(ji,jsl,1) = e(ji,jsl)
          tmat(ji,jsl,2) = f(ji,jsl)
          tmat(ji,jsl,3) = zero
          rhs(ji,jsl) = ep(ji,jsl)*mcl(ji,jsl-1,jst) + fp(ji,jsl)*mcl(ji,jsl,jst) &
               & + (b(ji,jsl-1) + b(ji,jsl)*(un-deux*free_drain_coef(ji,jst))) * (dt_sechiba/one_day) / deux &
               & - rootsink(ji,jsl,jst)
       ENDDO
       !- Store the equations in case needed again
       DO jsl=1,nslm
          DO ji = 1, kjpindex
             srhs(ji,jsl) = rhs(ji,jsl)
             stmat(ji,jsl,1) = tmat(ji,jsl,1)
             stmat(ji,jsl,2) = tmat(ji,jsl,2)
             stmat(ji,jsl,3) = tmat(ji,jsl,3) 
          ENDDO
       ENDDO

!       WRITE(numout,*) 'mcl step2 done'
       
       !! 2.7 Solves diffusion equations, but only in grid-cells where resolv is true, i.e. everywhere (cf 2.2)
       !!     The result is an updated mcl profile

       CALL hydrol_soil_tridiag(kjpindex,jst)

       !! 2.8 Computes drainage = bottom boundary condition, consistent with rhs(ji,jsl=nslm)
       ! dr_ns in mm/dt_sechiba, from k in mm/d
       ! This should be done where resolv=T, like tridiag (drainage is part of the linear system !)
       DO ji = 1, kjpindex
          IF (resolv(ji)) THEN
             dr_ns(ji,jst) = mask_soiltile(ji,jst)*k(ji,nslm)*free_drain_coef(ji,jst) * (dt_sechiba/one_day)
          ELSE
             dr_ns(ji,jst) = zero
          ENDIF
       ENDDO

       !! 2.9 For water conservation check during redistribution AND CORRECTION, 
       !!     we calculate the total liquid SM at the end of the routine tridiag
       tmcf(:) = dz(2) * ( trois*mcl(:,1,jst) + mcl(:,2,jst) )/huit
       DO jsl = 2,nslm-1
          tmcf(:) = tmcf(:) + dz(jsl) * (trois*mcl(:,jsl,jst)+mcl(:,jsl-1,jst))/huit &
               + dz(jsl+1) * (trois*mcl(:,jsl,jst)+mcl(:,jsl+1,jst))/huit
       ENDDO
       tmcf(:) = tmcf(:) + dz(nslm) * (trois*mcl(:,nslm,jst) + mcl(:,nslm-1,jst))/huit
          
       !! And we compare the difference with the flux...
       ! Normally, tcmf=tmci-flux_top(ji)-transpir-dr_ns
       DO ji=1,kjpindex
          diag_tr(ji)=SUM(rootsink(ji,:,jst))
       ENDDO
       ! Here, check_tr_ns holds the inaccuracy during the redistribution phase
       check_tr_ns(:,jst) = tmcf(:)-(tmci(:)-flux_top(:)-dr_ns(:,jst)-diag_tr(:))

       !! We solve here the numerical errors that happen when the soil is close to saturation
       !! and drainage very high, and which lead to negative check_tr_ns: the soil dries more 
       !! than what is demanded by the fluxes, so we need to increase the fluxes.
       !! This is done by increasing the drainage.
       !! There are also instances of positive check_tr_ns, larger when the drainage is high
       !! They are similarly corrected by a decrease of dr_ns, in the limit of keeping a positive drainage.
       DO ji=1,kjpindex
          IF ( check_tr_ns(ji,jst) .LT. zero ) THEN
              dr_corrnum_ns(ji,jst) = -check_tr_ns(ji,jst)
          ELSE
              dr_corrnum_ns(ji,jst) = -MIN(dr_ns(ji,jst),check_tr_ns(ji,jst))             
          ENDIF
          dr_ns(ji,jst) = dr_ns(ji,jst) + dr_corrnum_ns(ji,jst) ! dr_ns increases/decrease if check_tr negative/positive
       ENDDO
       !! For water conservation check during redistribution
       IF (check_cwrr2) THEN          
          check_tr_ns(:,jst) = tmcf(:)-(tmci(:)-flux_top(:)-dr_ns(:,jst)-diag_tr(:)) 
       ENDIF

       !! 3. AFTER DIFFUSION/REDISTRIBUTION

       !! 3.1 Updating mc, as all the following checks against saturation will compare mc to mcs
       !      The frozen fraction is constant, so that any water flux to/from a layer changes
       !      both mcl and the ice amount. The assumption behind this is that water entering/leaving
       !      a soil layer immediately freezes/melts with the proportion profil_froz_hydro_ns/(1-profil_...)
       DO jsl = 1, nslm
          DO ji =1, kjpindex
             mc(ji,jsl,jst) = MAX( mcl(ji,jsl,jst), mcl(ji,jsl,jst) + &
                  profil_froz_hydro_ns(ji,jsl,jst)*(mc(ji,jsl,jst)-mcr(njsc(ji))) )
             ! if profil_froz_hydro_ns=0 (including NOT ok_freeze_cwrr) we get mc=mcl
          ENDDO
       ENDDO

       !! 3.2 Correct here the possible over-saturation values (subroutine hydrol_soil_smooth_over_mcs2 acts on mc)
       !    Oversaturation results from numerical inaccuracies and can be frequent if free_drain_coef=0
       !    Here hydrol_soil_smooth_over_mcs2 discard all excess as ru_corr_ns, oriented to either ru_ns or dr_ns
       !    The former routine hydrol_soil_smooth_over_mcs, which keeps most of the excess in the soiltile
       !    after smoothing, first downward then upward, is kept in the module but not used here 
       dr_corr_ns(:,jst) = zero
       ru_corr_ns(:,jst) = zero
       call hydrol_soil_smooth_over_mcs2(kjpindex, jst, njsc, is_over_mcs, ru_corr_ns, check_over_ns)
       
       ! In absence of freezing, if F is large enough, the correction of oversaturation is sent to drainage       
       DO ji = 1, kjpindex
          IF ((free_drain_coef(ji,jst) .GE. 0.5) .AND. (.NOT. ok_freeze_cwrr) ) THEN
             dr_corr_ns(ji,jst) = ru_corr_ns(ji,jst) 
             ru_corr_ns(ji,jst) = zero
          ENDIF
       ENDDO
       dr_ns(:,jst) = dr_ns(:,jst) + dr_corr_ns(:,jst)
       ru_ns(:,jst) = ru_ns(:,jst) + ru_corr_ns(:,jst)

       !! 3.3 Negative runoff is reported to drainage
       !  Since we computed ru_ns directly from hydrol_soil_infilt, ru_ns should not be negative
             
       ru_corr2_ns(:,jst) = zero
       DO ji = 1, kjpindex
          IF (ru_ns(ji,jst) .LT. zero) THEN
             IF (printlev>=3)  WRITE (numout,*) 'NEGATIVE RU_NS: runoff and drainage before correction',&
                  ru_ns(ji,jst),dr_ns(ji,jst)
             dr_ns(ji,jst)=dr_ns(ji,jst)+ru_ns(ji,jst)
             ru_corr2_ns(ji,jst) = -ru_ns(ji,jst)
             ru_ns(ji,jst)= 0.
          END IF          
       ENDDO

       !! 3.4 Optional block to force saturation below zwt_force
       ! This block is not compatible with freezing; in this case, mcl must be corrected too
       ! We test if zwt_force(1,jst) <= zmaxh, to avoid steps 1 and 2 if unnecessary
       
       IF (zwt_force(1,jst) <= zmaxh) THEN

          !! We force the nodes below zwt_force to be saturated
          !  As above, we compare mc to mcs 
          DO jsl = 1,nslm
             DO ji = 1, kjpindex
                dmc(ji,jsl) = zero
                IF ( ( zz(jsl) >= zwt_force(ji,jst)*mille ) ) THEN
                   dmc(ji,jsl) = mcs(ji) - mc(ji,jsl,jst) ! addition to reach mcs (m3/m3) = positive value
                   mc(ji,jsl,jst) = mcs(ji)
                ENDIF
             ENDDO
          ENDDO
          
          !! To ensure conservation, this needs to be balanced by a negative change in drainage (in kg/m2/dt)
          DO ji = 1, kjpindex
             dr_force_ns(ji,jst) = dz(2) * ( trois*dmc(ji,1) + dmc(ji,2) )/huit ! top layer = initialization
          ENDDO
          DO jsl = 2,nslm-1 ! intermediate layers
             DO ji = 1, kjpindex
                dr_force_ns(ji,jst) = dr_force_ns(ji,jst) + dz(jsl) &
                     & * (trois*dmc(ji,jsl)+dmc(ji,jsl-1))/huit &
                     & + dz(jsl+1) * (trois*dmc(ji,jsl)+dmc(ji,jsl+1))/huit
             ENDDO
          ENDDO
          DO ji = 1, kjpindex
             dr_force_ns(ji,jst) = dr_force_ns(ji,jst) + dz(nslm) & ! bottom layer
                  & * (trois * dmc(ji,nslm) + dmc(ji,nslm-1))/huit
             dr_ns(ji,jst) = dr_ns(ji,jst) - dr_force_ns(ji,jst) ! dr_force_ns is positive and dr_ns must be reduced
          END DO

       ELSE          

          dr_force_ns(:,jst) = zero 

       ENDIF

       !! 3.5 Diagnosing the effective water table depth:
       !!     Defined as as the smallest jsl value when mc(jsl) is no more at saturation (mcs), starting from the bottom
       !      If there is a part of the soil which is saturated but underlain with unsaturated nodes,
       !      this is not considered as a water table
       DO ji = 1, kjpindex
          wtd_ns(ji,jst) = undef_sechiba ! in meters
          jsl=nslm
          DO WHILE ( (mc(ji,jsl,jst) .EQ. mcs(ji)) .AND. (jsl > 1) )
             wtd_ns(ji,jst) = zz(jsl)/mille ! in meters
             jsl=jsl-1   
          ENDDO
       ENDDO

       !! 3.6 Diagnose under_mcr to adapt water stress calculation below
       !      This routine does not change tmc but decides where we should turn off ET to prevent further mc decrease
       !      Like above, the tests are made on total mc, compared to mcr
       CALL hydrol_soil_smooth_under_mcr(kjpindex, jst, njsc, is_under_mcr, check_under_ns)
 
       !! 4. At the end of the prognostic calculations, we recompute important moisture variables

       !! 4.1 Total soil moisture content (water2infilt added below)
       DO ji = 1, kjpindex
          tmc(ji,jst) = dz(2) * ( trois*mc(ji,1,jst) + mc(ji,2,jst) )/huit
       ENDDO
       DO jsl = 2,nslm-1
          DO ji = 1, kjpindex
             tmc(ji,jst) = tmc(ji,jst) + dz(jsl) &
                  & * (trois*mc(ji,jsl,jst)+mc(ji,jsl-1,jst))/huit &
                  & + dz(jsl+1) * (trois*mc(ji,jsl,jst)+mc(ji,jsl+1,jst))/huit
          ENDDO
       ENDDO
       DO ji = 1, kjpindex
          tmc(ji,jst) = tmc(ji,jst) + dz(nslm) &
               & * (trois * mc(ji,nslm,jst) + mc(ji,nslm-1,jst))/huit
       END DO

       !! 4.2 mcl is a module variable; we update it here for calculating bare soil evaporation,
       !!     and in case we would like to export it (xios)
       DO jsl = 1, nslm
          DO ji =1, kjpindex
             mcl(ji,jsl,jst)= MIN( mc(ji,jsl,jst), mcr(njsc(ji)) + &
                  (un-profil_froz_hydro_ns(ji,jsl,jst))*(mc(ji,jsl,jst)-mcr(njsc(ji))) )
             ! if profil_froz_hydro_ns=0 (including NOT ok_freeze_cwrr) we keep mcl=mc
          ENDDO
       ENDDO
       
       !! 5. Optional check of the water balance of soil column (if check_cwrr)

       IF (check_cwrr) THEN

          !! 5.1 Computation of the vertical water fluxes
          CALL hydrol_soil_flux(kjpindex,jst,mclint,flux_top)
          
          !! 5.2 Total mc conservation 
          DO ji = 1,kjpindex   
             deltahum(ji) = (tmc(ji,jst) - tmcold(ji))
             diff(ji) = flux_infilt(ji) - flux_top(ji) - SUM(rootsink(ji,:,jst)) &
                   -ru_ns(ji,jst) - dr_ns(ji,jst)
             test(ji) = (ABS(deltahum(ji)-diff(ji))*mask_soiltile(ji,jst) .GT. allowed_err)
 
             IF (test(ji)) THEN              
                WRITE (numout,*)'CWRR water conservation pb:',ji,jst,njsc(ji),deltahum(ji)-diff(ji)
                WRITE (numout,*)'tmc,tmcold,diff',tmc(ji,jst),tmcold(ji),deltahum(ji)
                WRITE(numout,*) 'evapot,evapot_penm,ae_ns,flux_top',evapot(ji),evapot_penm(ji),&
                     ae_ns(ji,jst),flux_top(ji)
                WRITE (numout,*)'ru_ns,dr_ns,SUM(rootsink)',ru_ns(ji,jst),dr_ns(ji,jst), &
                     SUM(rootsink(ji,:,jst))
                WRITE (numout,*)'precisol, flux_infilt',precisol_ns(ji,jst)
                WRITE (numout,*)'irrigation, reinfiltration', &
                      irrigation_soil(ji,jst),reinfiltration_soil(ji)
                WRITE (numout,*)'mc',mc(ji,:,jst) ! along jsl
                WRITE (numout,*)'qflux',qflux(ji,:,jst) ! along jsl
                WRITE (numout,*)'k', k(ji,:) ! along jsl
                WRITE (numout,*)'soiltile',soiltile(ji,jst)
                WRITE (numout,*)'veget_max', veget_max(ji,:)
                
                error=.TRUE.
                CALL ipslerr_p(2, 'hydrol_soil', 'We will STOP in the end of this subroutine.',&
                     & 'CWRR water balance check','')
             ENDIF
          ENDDO

          !! 5.3 Total mc should not reach zero, or the tridiag solver will have problems
          DO ji = 1,kjpindex
             IF(MINVAL(mc(ji,:,jst)).LT. min_sechiba) THEN
                WRITE (numout,*)'CWRR MC NEGATIVE', &
                     ji,lalo(ji,:),MINLOC(mc(ji,:,jst)),jst,mc(ji,:,jst)
                WRITE(numout,*) 'evapot,evapot_penm,ae_ns,flux_top',evapot(ji),evapot_penm(ji),&
                     ae_ns(ji,jst),flux_top(ji)
                WRITE (numout,*)'ru_ns,dr_ns,SUM(rootsink)',ru_ns(ji,jst),dr_ns(ji,jst), &
                     SUM(rootsink(ji,:,jst))
                WRITE (numout,*)'precisol, flux_infilt',precisol_ns(ji,jst)
                WRITE (numout,*)'irrigation, reinfiltration', &
                      irrigation_soil(ji,jst),reinfiltration_soil(ji)
                WRITE (numout,*)'mc',mc(ji,:,jst) ! along jsl
                WRITE (numout,*)'qflux',qflux(ji,:,jst) ! along jsl
                WRITE (numout,*)'k', k(ji,:) ! along jsl
                WRITE (numout,*)'soiltile',soiltile(ji,jst)
                WRITE (numout,*)'veget_max', veget_max(ji,:)              

                error=.TRUE.
                CALL ipslerr_p(2, 'hydrol_soil', 'We will STOP in the end of this subroutine.',&
                     & 'CWRR MC NEGATIVE','')
             ENDIF
          END DO

       ENDIF

       !! 6. SM DIAGNOSTICS FOR OTHER ROUTINES, MODULES, OR NEXT STEP
       !    Starting here, mc and mcl should not change anymore 
       
       !! 6.1 Total soil moisture, soil moisture at litter levels, soil wetness, us, humrelv, vesgtressv
       !!     (based on mc)

       !! In output, tmc includes water2infilt(ji,jst)
       DO ji=1,kjpindex
          tmc(ji,jst) = tmc(ji,jst) + water2infilt(ji,jst)
       END DO
!gmjc top 5 layer mc for grazing
       ! The trampling depth is the 5 top levels of the soil
       ! Compute various field of soil moisture for the litter (used for stomate
       ! and for albedo)

       DO ji=1,kjpindex
          tmc_trampling(ji,jst) = dz(2) * (trois*mc(ji,1,jst)+mc(ji,2,jst))/huit
          tmc_trampling(ji,jst) = tmc_trampling(ji,jst)
       END DO

       ! sum from level 1 to 5

       DO jsl=2,6
          DO ji=1,kjpindex
             tmc_trampling(ji,jst) = tmc_trampling(ji,jst) + dz(jsl) * &
                  & ( trois*mc(ji,jsl,jst) + mc(ji,jsl-1,jst))/huit &
                  & + dz(jsl+1)*(trois*mc(ji,jsl,jst) + mc(ji,jsl+1,jst))/huit
          END DO
       END DO

!    tmc_topgrass(:) = tmc_trampling(:,3)
!WRITE (numout,*) 'sechiba tmc_trampling',tmc_trampling(:,jst),tmc(:,jst)
!WRITE (numout,*) 'sechiba mc',mc(:,1,jst)
!end gmjc

       ! The litter is the 4 top levels of the soil
       ! Compute various field of soil moisture for the litter (used for stomate and for albedo)
       ! We exclude the frozen water from the calculation
       DO ji=1,kjpindex
          tmc_litter(ji,jst) = dz(2) * ( trois*mcl(ji,1,jst)+ mcl(ji,2,jst))/huit
       END DO
       ! sum from level 1 to 4
       DO jsl=2,4
          DO ji=1,kjpindex
             tmc_litter(ji,jst) = tmc_litter(ji,jst) + dz(jsl) * & 
                  & ( trois*mcl(ji,jsl,jst) + mcl(ji,jsl-1,jst))/huit &
                  & + dz(jsl+1)*(trois*mcl(ji,jsl,jst) + mcl(ji,jsl+1,jst))/huit
          END DO
       END DO

       ! Subsequent calculation of soil_wet_litter (tmc-tmcw)/(tmcfc-tmcw)
       ! Based on liquid water content
       DO ji=1,kjpindex
          soil_wet_litter(ji,jst) = MIN(un, MAX(zero,&
               & (tmc_litter(ji,jst)-tmc_litter_wilt(ji,jst)) / &
               & (tmc_litter_field(ji,jst)-tmc_litter_wilt(ji,jst)) ))
       END DO
    !!xchen                                                                     
    !!DO ji = 1, kjpindex                                                         
    !!   DO jsl = 1,nslm                                                          
    !!      IF (swdown(ji) .GT. 180) THEN                                         
    !!         !!mcw_mineral(njsc(ji)) = 1.5*mcw_mineral(njsc(ji))                  
    !!         dz(jsl) = 1.5*dz(jsl)
    !!      ENDIF                                                                 
    !!   ENDDO                                                                    
    !!ENDDO    
       ! Preliminary calculation of various soil moistures (for each layer, in kg/m2)
       sm(:,1)  = dz(2) * (trois*mcl(:,1,jst) + mcl(:,2,jst))/huit
       smt(:,1) = dz(2) * (trois*mc(:,1,jst) + mc(:,2,jst))/huit
       smw(:,1) = dz(2) * (quatre*mcw(:))/huit
       smf(:,1) = dz(2) * (quatre*mcfc(:))/huit
       sms(:,1) = dz(2) * (quatre*mcs(:))/huit
       smw_tmp(:,1) = dz(2) * (quatre*mcw_mineral(njsc(:)))/huit
       smf_tmp(:,1) = dz(2) * (quatre*mcfc_mineral(njsc(:)))/huit
       sms_tmp(:,1) = dz(2) * (quatre*mcs_mineral(njsc(:)))/huit
       DO jsl = 2,nslm-1
          sm(:,jsl)  = dz(jsl) * (trois*mcl(:,jsl,jst)+mcl(:,jsl-1,jst))/huit &
               + dz(jsl+1) * (trois*mcl(:,jsl,jst)+mcl(:,jsl+1,jst))/huit
          smt(:,jsl) = dz(jsl) * (trois*mc(:,jsl,jst)+mc(:,jsl-1,jst))/huit &
               + dz(jsl+1) * (trois*mc(:,jsl,jst)+mc(:,jsl+1,jst))/huit
          smw(:,jsl) = dz(jsl) * ( quatre*mcw(:) )/huit &
               + dz(jsl+1) * ( quatre*mcw(:) )/huit
          smf(:,jsl) = dz(jsl) * ( quatre*mcfc(:) )/huit &
               + dz(jsl+1) * ( quatre*mcfc(:) )/huit
          sms(:,jsl) = dz(jsl) * ( quatre*mcs(:) )/huit &
               + dz(jsl+1) * ( quatre*mcs(:) )/huit
          smw_tmp(:,jsl) = dz(jsl) * ( quatre*mcw_mineral(njsc(:)) )/huit + dz(jsl+1) * ( quatre*mcw_mineral(njsc(:)) )/huit
          smf_tmp(:,jsl) = dz(jsl) * ( quatre*mcfc_mineral(njsc(:)) )/huit + dz(jsl+1) * ( quatre*mcfc_mineral(njsc(:)) )/huit
          sms_tmp(:,jsl) = dz(jsl) * ( quatre*mcs_mineral(njsc(:)) )/huit + dz(jsl+1) * ( quatre*mcs_mineral(njsc(:)) )/huit
       ENDDO
       sm(:,nslm)  = dz(nslm) * (trois*mcl(:,nslm,jst) + mcl(:,nslm-1,jst))/huit  
       smt(:,nslm) = dz(nslm) * (trois*mc(:,nslm,jst) + mc(:,nslm-1,jst))/huit      
       smw(:,nslm) = dz(nslm) * (quatre*mcw(:))/huit
       smf(:,nslm) = dz(nslm) * (quatre*mcfc(:))/huit
       sms(:,nslm) = dz(nslm) * (quatre*mcs(:))/huit
       smw_tmp(:,nslm) = dz(nslm) * (quatre*mcw_mineral(njsc(:)))/huit
       smf_tmp(:,nslm) = dz(nslm) * (quatre*mcfc_mineral(njsc(:)))/huit
       sms_tmp(:,nslm) = dz(nslm) * (quatre*mcs_mineral(njsc(:)))/huit
       ! sm_nostress = soil moisture of each layer at which us reaches 1, here at the middle of [smw,smf]
       DO jsl = 1,nslm
          sm_nostress(:,jsl) = smw(:,jsl) + pcent(njsc(:)) * (smf(:,jsl)-smw(:,jsl))
       END DO

       ! Saturated litter soil moisture for rsoil
       tmcs_litter(:) = zero
       DO jsl = 1,4
          tmcs_litter(:) = tmcs_litter(:) + sms(:,jsl)
       END DO
             
       ! Soil wetness profiles (W-Ww)/(Ws-Ww)
       ! soil_wet_ns is the ratio of available soil moisture to max available soil moisture
       ! (ie soil moisture at saturation minus soil moisture at wilting point).
       ! soil wet is a water stress for stomate, to control C decomposition
       ! Based on liquid water content
       DO jsl=1,nslm
          DO ji=1,kjpindex
             soil_wet_ns(ji,jsl,jst) = MIN(un, MAX(zero, &
                  (sm(ji,jsl)-smw_tmp(ji,jsl))*(sms(ji,jsl)-smw(ji,jsl))/(sms_tmp(ji,jsl)-smw_tmp(ji,jsl))**2 ))
          END DO
       END DO

       ! Compute us and the new humrelv to use in sechiba (with loops on the vegetation types)
       ! This is the water stress for transpiration (diffuco) and photosynthesis (diffuco)
       ! humrel is never used in stomate
       ! Based on liquid water content

       ! -- PFT1
       humrelv(:,1,jst) = zero        
       ! -- Top layer
       DO jv = 2,nvm
          DO ji=1,kjpindex
             !- Here we make the assumption that roots do not take water from the 1st layer. 
             us(ji,jv,jst,1) = zero
             humrelv(ji,jv,jst) = zero ! initialisation of the sum
          END DO
       ENDDO

       IF (ok_freeze_cwrr) THEN
           CALL hydrol_soil_coef(kjpindex,jst,njsc)
       ENDIF

       !! Dynamic nroot to optimize water use: the root profile used to weight the water stress function 
       !! of each soil layer is updated at each time step in order to match the soil water profile 
       !! (the soil water content of each layer available for transpiration)
       IF (ok_dynroot) THEN
          pcent_gc(:) = pcent(njsc(:))

          CALL hydrol_update_nroot(kjpindex, nslm, nvm, sm, smw, smf, sms, smf_tmp, smw_tmp, pcent_gc, nroot)
       ENDIF ! (ok_dynroot) THEN

       ! -- Intermediate and bottom layers
       DO jsl = 2,nslm
          DO jv = 2, nvm
             DO ji=1,kjpindex
                ! AD16*** Although plants can only withdraw liquid water, we compute here the water stress
                ! based on mc and the corresponding thresholds mcs, pcent, or potentially mcw and mcfc
                ! This is consistent with assuming that ice is uniformly distributed within the poral space
                ! In such a case, freezing makes mcl and the "liquid" porosity smaller than the "total" values
                ! And it is the same for all the moisture thresholds, which are proportional to porosity.
                ! Since the stress is based on relative moisture, it could thus independent from the porosity
                ! at first order, thus independent from freezing.   
                ! 26-07-2017: us and humrel now based on liquid soil moisture, so the stress is stronger
                IF(new_watstress) THEN
                   IF((sm(ji,jsl)-smw(ji,jsl)) .GT. min_sechiba) THEN
                      us(ji,jv,jst,jsl) = MIN(un, MAX(zero, &
                           (EXP(- alpha_watstress * &
                           ( (smf(ji,jsl) - smw(ji,jsl)) / ( sm_nostress(ji,jsl) - smw(ji,jsl)) ) * &
                           ( (sm_nostress(ji,jsl) - sm(ji,jsl)) / ( sm(ji,jsl) - smw(ji,jsl)) ) ) ) ))&
                           * nroot(ji,jv,jsl)
                   ELSE
                      us(ji,jv,jst,jsl) = 0.
                   ENDIF
                ELSE
                   us(ji,jv,jst,jsl) = MIN(un, MAX(zero, &
                        (sm(ji,jsl)-smw_tmp(ji,jsl))/(pcent(njsc(ji))*(smf_tmp(ji,jsl)-smw_tmp(ji,jsl)))*(smf(ji,jsl)-smw(ji,jsl))/(smf_tmp(ji,jsl)-smw_tmp(ji,jsl)) )) * nroot(ji,jv,jsl)
                ENDIF
                humrelv(ji,jv,jst) = humrelv(ji,jv,jst) + us(ji,jv,jst,jsl)
             END DO
          END DO
       ENDDO

       !! vegstressv is the water stress for phenology in stomate
       !! It varies linearly from zero at wilting point to 1 at field capacity
       vegstressv(:,:,jst) = zero
       DO jv = 2, nvm
          DO ji=1,kjpindex
             DO jsl=1,nslm
                vegstressv(ji,jv,jst) = vegstressv(ji,jv,jst) + &
                     MIN(un, MAX(zero, (sm(ji,jsl)-smw_tmp(ji,jsl))*(smf(ji,jsl)-smw(ji,jsl))/(smf_tmp(ji,jsl)-smw_tmp(ji,jsl))**2 ) ) &
                     * nroot(ji,jv,jsl)
             END DO
          END DO
       END DO


       ! -- If the PFT is absent, the corresponding humrelv and vegstressv = 0
       DO jv = 2, nvm
          DO ji = 1, kjpindex
             IF (vegetmax_soil(ji,jv,jst) .LT. min_sechiba) THEN
                humrelv(ji,jv,jst) = zero
                vegstressv(ji,jv,jst) = zero
                us(ji,jv,jst,:) = zero
             ENDIF
          END DO
       END DO

       !! 6.2 We need to turn off evaporation when is_under_mcr
       !!     We set us, humrelv and vegstressv to zero in this case
       !!     WARNING: It's different from having locally us=0 in the soil layers(s) where mc<mcr
       !!              This part is crucial to preserve water conservation
       DO jsl = 1,nslm
          DO jv = 2, nvm
             WHERE (is_under_mcr(:,jst))
                us(:,jv,jst,jsl) = zero
             ENDWHERE
          ENDDO
       ENDDO
       DO jv = 2, nvm
          WHERE (is_under_mcr(:,jst))
             humrelv(:,jv,jst) = zero
          ENDWHERE
       ENDDO
       
       ! For consistency in stomate, we also set moderwilt and soil_wet_ns to zero in this case. 
       ! They are used later for shumdiag and shumdiag_perma
       DO jsl = 1,nslm
          WHERE (is_under_mcr(:,jst))
             soil_wet_ns(:,jsl,jst) = zero
          ENDWHERE
       ENDDO

       ! Counting the nb of under_mcr occurences in each grid-cell 
       WHERE (is_under_mcr(:,jst))
          undermcr = undermcr + un
       ENDWHERE

       !! 6.3 Calculate the volumetric soil moisture content (mc_layh and mcl_layh) needed in 
       !!     thermosoil for the thermal conductivity. 
       !! The multiplication by vegtot creates grid-cell average values
       ! *** To be checked for consistency with the use of nobio properties in thermosoil
       mc_layh_s = MAX(mc, min_sechiba)
       mcl_layh_s = MAX(mc, min_sechiba)
       DO jsl=1,nslm
          DO ji=1,kjpindex
             mc_layh(ji,jsl) = mc_layh(ji,jsl) + mc(ji,jsl,jst) * soiltile(ji,jst) * vegtot(ji) 
             mcl_layh(ji,jsl) = mcl_layh(ji,jsl) + mcl(ji,jsl,jst) * soiltile(ji,jst) * vegtot(ji)
          ENDDO
       END DO

       !! 6.4 The hydraulic conductivities exported here are the ones used in the diffusion/redistribution
       ! (no call of hydrol_soil_coef since 2.1)
       ! We average the values of each soiltile and keep the specific value (no multiplication by vegtot)
       DO ji = 1, kjpindex
          kk_moy(ji,:) = kk_moy(ji,:) + soiltile(ji,jst) * k(ji,:) 
          kk(ji,:,jst) = k(ji,:)
       ENDDO
       
       !! 6.5 We also want to export ksat at each node for CMIP6
       !  (In the output, done only once according to field_def_orchidee.xml; same averaging as for kk)
       DO jsl = 1, nslm
          ksat(:,jsl) = ksat(:,jsl) + soiltile(:,jst) * &
               ( ks(njsc(:)) * kfact(jsl,njsc(:)) * kfact_root(:,jsl,jst) ) 
       ENDDO
              
      IF (printlev>=3) WRITE (numout,*) ' prognostic/diagnostic part of hydrol_soil done for jst =', jst          

    END DO  ! end of loop on soiltile 

!gmjc top 5 layer grassland soil moisture for grazing
    ! should be calculated after loop soiltile
    ! tmc_trampling unit mm water
    ! for soil moisture, it should be divided by 5 layer soil depth
    tmc_topgrass(:) = tmc_trampling(:,3)/(SUM(dz(1:6))+dz(7)/2)
!WRITE (numout,*) 'sechiba tmc',tmc(:,jst),tmc_topgrass(:)
!end gmjc

    !! -- ENDING THE MAIN LOOP ON SOILTILES

    !! 7. Summing 3d variables into 2d variables 
    CALL hydrol_diag_soil (kjpindex, veget_max, soiltile, njsc, runoff, drainage, &
         & evapot, vevapnu, returnflow, reinfiltration, irrigation, &
         & shumdiag,shumdiag_perma, k_litt, litterhumdiag, humrel, vegstress, drysoil_frac,tot_melt)
    

    ! Means of wtd, runoff and drainage corrections, across soiltiles    
    wtd(:) = zero 
    ru_corr(:) = zero
    ru_corr2(:) = zero
    dr_corr(:) = zero
    dr_corrnum(:) = zero
    dr_force(:) = zero
    DO jst = 1, nstm
       DO ji = 1, kjpindex 
          wtd(ji) = wtd(ji) + soiltile(ji,jst) * wtd_ns(ji,jst) ! average over vegtot only
          IF (vegtot(ji) .GT. min_sechiba) THEN ! to mimic hydrol_diag_soil
             ! We average the values of each soiltile and multiply by vegtot to transform to a grid-cell mean
             ru_corr(ji) = ru_corr(ji) + vegtot(ji) * soiltile(ji,jst) * ru_corr_ns(ji,jst) 
             ru_corr2(ji) = ru_corr2(ji) + vegtot(ji) * soiltile(ji,jst) * ru_corr2_ns(ji,jst) 
             dr_corr(ji) = dr_corr(ji) + vegtot(ji) * soiltile(ji,jst) * dr_corr_ns(ji,jst) 
             dr_corrnum(ji) = dr_corrnum(ji) + vegtot(ji) * soiltile(ji,jst) * dr_corrnum_ns(ji,jst)
             dr_force(ji) = dr_force(ji) - vegtot(ji) * soiltile(ji,jst) * dr_force_ns(ji,jst)
                                       ! the sign is OK to get a negative drainage flux
          ENDIF
       ENDDO
    ENDDO

    ! Means local variables, including water conservation checks
    ru_infilt(:)=0.
    qinfilt(:)=0.
    check_infilt(:)=0.
    check_tr(:)=0.
    check_over(:)=0.
    check_under(:)=0.
    DO jst = 1, nstm
       DO ji = 1, kjpindex 
          IF (vegtot(ji) .GT. min_sechiba) THEN ! to mimic hydrol_diag_soil
             ! We average the values of each soiltile and multiply by vegtot to transform to a grid-cell mean
             ru_infilt(ji) = ru_infilt(ji) + vegtot(ji) * soiltile(ji,jst) * ru_infilt_ns(ji,jst)
             qinfilt(ji) = qinfilt(ji) + vegtot(ji) * soiltile(ji,jst) * qinfilt_ns(ji,jst)
          ENDIF
       ENDDO
    ENDDO
 
    IF (check_cwrr2) THEN
       DO jst = 1, nstm
          DO ji = 1, kjpindex 
             IF (vegtot(ji) .GT. min_sechiba) THEN ! to mimic hydrol_diag_soil
                ! We average the values of each soiltile and multiply by vegtot to transform to a grid-cell mean
                check_infilt(ji) = check_infilt(ji) + vegtot(ji) * soiltile(ji,jst) * check_infilt_ns(ji,jst)
                check_tr(ji) = check_tr(ji) + vegtot(ji) * soiltile(ji,jst) * check_tr_ns(ji,jst)
                check_over(ji) = check_over(ji) + vegtot(ji) * soiltile(ji,jst) * check_over_ns(ji,jst)
                check_under(ji) =  check_under(ji) + vegtot(ji) * soiltile(ji,jst) * check_under_ns(ji,jst)
             ENDIF
          ENDDO
       ENDDO
    END IF

    !! 8. COMPUTING EVAP_BARE_LIM_NS FOR NEXT TIME STEP, WITH A LOOP ON SOILTILES
    !!    The principle is to run a dummy integration of the water redistribution scheme
    !!    to check if the SM profile can sustain a potential evaporation.
    !!    If not, the dummy integration is redone from the SM profile of the end of the normal integration,
    !!    with a boundary condition leading to a very severe water limitation: mc(1)=mcr

    ! evap_bare_lim = beta factor for bare soil evaporation
    evap_bare_lim(:) = zero
    evap_bare_lim_ns(:,:) = zero
!    evap_bare_lim_pft(:,:) = zero

    ! Loop on soil tiles  
    DO jst = 1,nstm

       !! 8.1 Save actual mc, mcl, and tmc for restoring at the end of the time step
       !!      and calculate tmcint corresponding to mc without water2infilt
       DO jsl = 1, nslm
          DO ji = 1, kjpindex
             mcint(ji,jsl) = mask_soiltile(ji,jst) * mc(ji,jsl,jst)
             mclint(ji,jsl) = mask_soiltile(ji,jst) * mcl(ji,jsl,jst)
          ENDDO
       ENDDO

       DO ji = 1, kjpindex
          temp(ji) = tmc(ji,jst)
          tmcint(ji) = temp(ji) - water2infilt(ji,jst) ! to estimate bare soil evap based on water budget
       ENDDO

       !! 8.2 Since we estimate bare soile evap for the next time step, we update profil_froz_hydro and mcl
       !     (effect of mc only, the change in temp_hydro is neglected)
       IF ( ok_freeze_cwrr ) CALL hydrol_soil_froz(kjpindex,jst,njsc)
        DO jsl = 1, nslm
          DO ji =1, kjpindex
             mcl(ji,jsl,jst)= MIN( mc(ji,jsl,jst), mcr(njsc(ji)) + &
                  (un-profil_froz_hydro_ns(ji,jsl,jst))*(mc(ji,jsl,jst)-mcr(njsc(ji))) )
             ! if profil_froz_hydro_ns=0 (including NOT ok_freeze_cwrr) we keep mcl=mc
          ENDDO
       ENDDO          

       !! 8.3 K and D are recomputed for the updated profile of mc/mcl
       CALL hydrol_soil_coef(kjpindex,jst,njsc)

       !! 8.4 Set the tridiagonal matrix coefficients for the diffusion/redistribution scheme
       CALL hydrol_soil_setup(kjpindex,jst)
       resolv(:) = (mask_soiltile(:,jst) .GT. 0) 

       !! 8.5 We define the system of linear equations, based on matrix coefficients, 

       !- Impose potential evaporation as flux_top in mm/step, assuming the water is available
       ! Note that this should lead to never have evapnu>evapot_penm(ji)

       DO ji = 1, kjpindex
          
          IF (vegtot(ji).GT.min_sechiba) THEN
             
             ! We calculate a reduced demand, by means of a soil resistance (Sellers et al., 1992)
             ! It is based on the liquid SM only, like for us and humrel
             IF (do_rsoil) THEN
                mc_rel(ji) = tmc_litter(ji,jst)/tmcs_litter(ji) ! tmc_litter based on mcl
                ! based on SM in the top 4 soil layers (litter) to smooth variability
                r_soil_ns(ji,jst) = exp(8.206 - 4.255 * mc_rel(ji))
             ELSE
                r_soil_ns(ji,jst) = zero
             ENDIF

             ! Aerodynamic resistance
             speed = MAX(min_wind, SQRT (u(ji)*u(ji) + v(ji)*v(ji)))
             IF (speed * tq_cdrag(ji) .GT. min_sechiba) THEN
                ra = un / (speed * tq_cdrag(ji))
                evap_soil(ji) = evapot_penm(ji) / (un + r_soil_ns(ji,jst)/ra)
             ELSE
                evap_soil(ji) = evapot_penm(ji)
             ENDIF
             
       ! AD16*** et si evap_bare_lim_ns<0 ?? car on suppose que tmcint > tmc(new)
       ! (water2inflit permet de propager de la ponded water d'un pas de temps a l'autre:
       ! peut-on s'en servir pour creer des cas d'evapnu potentielle negative ? a gerer dans diffuco ?)
             
             flux_top(ji) = evap_soil(ji) * &
                  AINT(frac_bare_ns(ji,jst)+un-min_sechiba)
          ELSE
             
             flux_top(ji) = zero
             
          ENDIF
       ENDDO

       IF (ok_freeze_cwrr) THEN
          CALL hydrol_soil_coef(kjpindex,jst,njsc)
          DO ji =1, kjpindex
             DO jsl = 1, nslm
                mcl(ji,jsl,jst)= MIN(mc(ji,jsl,jst),mcr(njsc(ji))+(1-profil_froz_hydro_ns(ji,jsl,jst))*(mc(ji,jsl,jst)-mcr(njsc(ji))))
             ENDDO
          ENDDO
       ELSE
          mcl(:,:,jst)=mc(:,:,jst)
       ENDIF

       ! We don't use rootsinks, because we assume there is no transpiration in the bare soil fraction (??)
       !- First layer
       DO ji = 1, kjpindex
          tmat(ji,1,1) = zero
          tmat(ji,1,2) = f(ji,1)
          tmat(ji,1,3) = g1(ji,1)
          rhs(ji,1)    = fp(ji,1) * mcl(ji,1,jst) + gp(ji,1)*mcl(ji,2,jst) &
               - flux_top(ji) - (b(ji,1)+b(ji,2))/deux *(dt_sechiba/one_day)
       ENDDO
       !- soil body
       DO jsl=2, nslm-1
          DO ji = 1, kjpindex
             tmat(ji,jsl,1) = e(ji,jsl)
             tmat(ji,jsl,2) = f(ji,jsl)
             tmat(ji,jsl,3) = g1(ji,jsl)
             rhs(ji,jsl) = ep(ji,jsl)*mcl(ji,jsl-1,jst) + fp(ji,jsl)*mcl(ji,jsl,jst) &
                  +  gp(ji,jsl) * mcl(ji,jsl+1,jst) &
                  + (b(ji,jsl-1) - b(ji,jsl+1)) * (dt_sechiba/one_day) / deux
          ENDDO
       ENDDO
       !- Last layer
       DO ji = 1, kjpindex
          jsl=nslm
          tmat(ji,jsl,1) = e(ji,jsl)
          tmat(ji,jsl,2) = f(ji,jsl)
          tmat(ji,jsl,3) = zero
          rhs(ji,jsl) = ep(ji,jsl)*mcl(ji,jsl-1,jst) + fp(ji,jsl)*mcl(ji,jsl,jst) &
               + (b(ji,jsl-1) + b(ji,jsl)*(un-deux*free_drain_coef(ji,jst))) * (dt_sechiba/one_day) / deux
       ENDDO
       !- Store the equations for later use (9.6)
       DO jsl=1,nslm
          DO ji = 1, kjpindex
             srhs(ji,jsl) = rhs(ji,jsl)
             stmat(ji,jsl,1) = tmat(ji,jsl,1)
             stmat(ji,jsl,2) = tmat(ji,jsl,2)
             stmat(ji,jsl,3) = tmat(ji,jsl,3)
          ENDDO
       ENDDO

       !! 8.6 Solve the diffusion equation, assuming that flux_top=evapot_penm (updates mcl)
       CALL hydrol_soil_tridiag(kjpindex,jst)

       !! 9.7 Alternative solution with mc(1)=mcr in points where the above solution leads to mcl<mcr 
       ! hydrol_soil_tridiag calculates mc recursively from the top as a fonction of rhs and tmat 
       ! We re-use these the above values, but for mc(1)=mcr and the related tmat
       
       DO ji = 1, kjpindex
          ! by construction, mc and mcl are always on the same side of mcr, so we can use mcl here
          resolv(ji) = (mcl(ji,1,jst).LT.(mcr(njsc(ji))).AND.flux_top(ji).GT.min_sechiba)
       ENDDO
       !! Reset the coefficient for diffusion (tridiag is only solved if resolv(ji) = .TRUE.)O
       DO jsl=1,nslm
          !- The new condition is to put the upper layer at residual soil moisture
          DO ji = 1, kjpindex
             rhs(ji,jsl) = srhs(ji,jsl)
             tmat(ji,jsl,1) = stmat(ji,jsl,1)
             tmat(ji,jsl,2) = stmat(ji,jsl,2)
             tmat(ji,jsl,3) = stmat(ji,jsl,3)
          END DO
       END DO
       
       DO ji = 1, kjpindex
          tmat(ji,1,2) = un
          tmat(ji,1,3) = zero
          rhs(ji,1) = mcr(njsc(ji))
       ENDDO
       
       ! Solves the diffusion equation with new surface bc where resolv=T 
       CALL hydrol_soil_tridiag(kjpindex,jst)

       ! Calculation of total soil moisture content (liquid + frozen)
       IF (ok_freeze_cwrr) THEN
          CALL hydrol_soil_coef(kjpindex,jst,njsc)           
          DO ji =1, kjpindex
             DO jsl = 1, nslm
                mc(ji,jsl,jst)=MAX(mcl(ji,jsl,jst), mcl(ji,jsl,jst)+profil_froz_hydro_ns(ji,jsl,jst)*(mc(ji,jsl,jst)-mcr(njsc(ji))))
             ENDDO
          ENDDO
       ELSE
          mc(:,:,jst)=mcl(:,:,jst)
       ENDIF


       !! Correct bad moisture content due to numerical errors before water balance
       !! 8.8 In both case, we have drainage to be consistent with rhs
       DO ji = 1, kjpindex
          flux_bottom(ji) = mask_soiltile(ji,jst)*k(ji,nslm)*free_drain_coef(ji,jst) * (dt_sechiba/one_day)
       ENDDO
       
       !! 8.9 Water budget to assess the top flux = soil evaporation
       !      Where resolv=F at the 2nd step (9.6), it should simply be the potential evaporation

       ! Total soil moisture content for water budget

       DO jsl = 1, nslm
          DO ji =1, kjpindex
             mc(ji,jsl,jst) = MAX( mcl(ji,jsl,jst), mcl(ji,jsl,jst) + &
                  profil_froz_hydro_ns(ji,jsl,jst)*(mc(ji,jsl,jst)-mcr(njsc(ji))) )
             ! if profil_froz_hydro_ns=0 (including NOT ok_freeze_cwrr) we get mc=mcl
          ENDDO
       ENDDO
       
       DO ji = 1, kjpindex
          tmc(ji,jst) = dz(2) * ( trois*mc(ji,1,jst) + mc(ji,2,jst) )/huit
       ENDDO       
       DO jsl = 2,nslm-1
          DO ji = 1, kjpindex
             tmc(ji,jst) = tmc(ji,jst) + dz(jsl) &
                  * (trois*mc(ji,jsl,jst)+mc(ji,jsl-1,jst))/huit &
                  + dz(jsl+1) * (trois*mc(ji,jsl,jst)+mc(ji,jsl+1,jst))/huit
          ENDDO
       ENDDO
       DO ji = 1, kjpindex
          tmc(ji,jst) = tmc(ji,jst) + dz(nslm) &
               * (trois * mc(ji,nslm,jst) + mc(ji,nslm-1,jst))/huit
       END DO
    
       ! Deduce upper flux from soil moisture variation and bottom flux
       ! TMCi-D-BSE=TMC (BSE=bare soil evap=TMCi-TMC-D)
       ! The numerical errors of tridiag close to saturation cannot be simply solved here,
       ! we can only hope they are not too large because we don't add water at this stage... 
       DO ji = 1, kjpindex
          evap_bare_lim_ns(ji,jst) = mask_soiltile(ji,jst) * &
               (tmcint(ji)-tmc(ji,jst)-flux_bottom(ji))
       END DO

       !! 8.10 evap_bare_lim_ns is turned from an evaporation rate to a beta
       DO ji = 1, kjpindex
          ! Here we weight evap_bare_lim_ns by the fraction of bare evaporating soil. 
          ! This is given by frac_bare_ns, taking into account bare soil under vegetation
          IF(vegtot(ji) .GT. min_sechiba) THEN
             evap_bare_lim_ns(ji,jst) = evap_bare_lim_ns(ji,jst) * frac_bare_ns(ji,jst)
          ELSE
             evap_bare_lim_ns(ji,jst) = 0.
          ENDIF
       END DO

       ! We divide by evapot, which is consistent with diffuco (evap_bare_lim_ns < evapot_penm/evapot)
       ! Further decrease if tmc_litter is below the wilting point

       IF (do_rsoil) THEN
          DO ji=1,kjpindex
             IF (evapot(ji).GT.min_sechiba) THEN
                evap_bare_lim_ns(ji,jst) = evap_bare_lim_ns(ji,jst) / evapot(ji)
             ELSE
                evap_bare_lim_ns(ji,jst) = zero ! not redundant with the is_under_mcr case below
                                                ! but not necessarily useful
             END IF
             evap_bare_lim_ns(ji,jst)=MAX(MIN(evap_bare_lim_ns(ji,jst),1.),0.)
          END DO
       ELSE
          DO ji=1,kjpindex
             IF ((evapot(ji).GT.min_sechiba) .AND. &
                  (tmc_litter(ji,jst).GT.(tmc_litter_wilt(ji,jst)))) THEN
                evap_bare_lim_ns(ji,jst) = evap_bare_lim_ns(ji,jst) / evapot(ji)
             ELSEIF((evapot(ji).GT.min_sechiba).AND. &
                  (tmc_litter(ji,jst).GT.(tmc_litter_res(ji,jst)))) THEN
                evap_bare_lim_ns(ji,jst) =  (un/deux) * evap_bare_lim_ns(ji,jst) / evapot(ji)
                ! This is very arbitrary, with no justification from the literature
             ELSE
                evap_bare_lim_ns(ji,jst) = zero
             END IF
             evap_bare_lim_ns(ji,jst)=MAX(MIN(evap_bare_lim_ns(ji,jst),1.),0.)
          END DO
       ENDIF

       !! 8.11 Set evap_bare_lim_ns to zero if is_under_mcr at the end of the prognostic loop
       !!      (cf us, humrelv, vegstressv in 5.2)
       WHERE (is_under_mcr(:,jst))
          evap_bare_lim_ns(:,jst) = zero
       ENDWHERE

       !! 8.12 Restores mc, mcl, and tmc, to erase the effect of the dummy integrations
       !!      on these prognostic variables 
       DO jsl = 1, nslm
          DO ji = 1, kjpindex
             mc(ji,jsl,jst) = mask_soiltile(ji,jst) * mcint(ji,jsl)
             mcl(ji,jsl,jst) = mask_soiltile(ji,jst) * mclint(ji,jsl)
          ENDDO
       ENDDO
       DO ji = 1, kjpindex
          tmc(ji,jst) = temp(ji)
       ENDDO

    ENDDO !end loop on tiles for dummy integration

    !! 9. evap_bar_lim is the grid-cell scale beta
    DO ji = 1, kjpindex
       evap_bare_lim(ji) =  SUM(evap_bare_lim_ns(ji,:)*vegtot(ji)*soiltile(ji,:))
       r_soil(ji) =  SUM(r_soil_ns(ji,:)*vegtot(ji)*soiltile(ji,:))
    ENDDO

    !! 10. XIOS export of local variables, including water conservation checks
    
    CALL xios_orchidee_send_field("ksat",ksat) ! mm/d (for CMIP6, once)
    CALL xios_orchidee_send_field("psi_moy",psi_moy) ! mm (for SP-MIP)
    CALL xios_orchidee_send_field("wtd",wtd) ! in m
    CALL xios_orchidee_send_field("ru_corr",ru_corr/dt_sechiba)   ! adjustment flux added to surface runoff (included in runoff)
    CALL xios_orchidee_send_field("ru_corr2",ru_corr2/dt_sechiba)
    CALL xios_orchidee_send_field("dr_corr",dr_corr/dt_sechiba)   ! adjustment flux added to drainage (included in drainage)
    CALL xios_orchidee_send_field("dr_corrnum",dr_corrnum/dt_sechiba) 
    CALL xios_orchidee_send_field("dr_force",dr_force/dt_sechiba) ! adjustement flux added to drainage to sustain a forced wtd
    CALL xios_orchidee_send_field("qinfilt",qinfilt/dt_sechiba)
    CALL xios_orchidee_send_field("ru_infilt",ru_infilt/dt_sechiba)
    CALL xios_orchidee_send_field("r_soil",r_soil) ! s/m

    IF (check_cwrr2) THEN
       CALL xios_orchidee_send_field("check_infilt",check_infilt/dt_sechiba)
       CALL xios_orchidee_send_field("check_tr",check_tr/dt_sechiba)
       CALL xios_orchidee_send_field("check_over",check_over/dt_sechiba)
       CALL xios_orchidee_send_field("check_under",check_under/dt_sechiba)    
    END IF

    !! 11. Exit if error was found previously in this subroutine
    
    IF ( error ) THEN
       WRITE(numout,*) 'One or more errors have been detected in hydrol_soil. Model stops.'
       CALL ipslerr_p(3, 'hydrol_soil', 'We will STOP now.',&
                  & 'One or several fatal errors were found previously.','')
    END IF

    IF (printlev>=3) WRITE(numout,*) 'hydrol_soil done'

  END SUBROUTINE hydrol_soil


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_soil_infilt
!!
!>\BRIEF        Infiltration
!!
!! DESCRIPTION  :
!! 1. We calculate the total SM at the beginning of the routine
!! 2. Infiltration process
!! 2.1 Initialization of time counter and infiltration rate
!! 2.2 Infiltration layer by layer, accounting for an exponential law for subgrid variability
!! 2.3 Resulting infiltration and surface runoff
!! 3. For water conservation check, we calculate the total SM at the beginning of the routine,
!!    and export the difference with the flux
!! 5. Local verification 
!!
!! RECENT CHANGE(S) : 2016 by A. Ducharne
!! Adding checks and interactions variables with hydrol_soil, but the processes are unchanged
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_soil_infilt

  SUBROUTINE hydrol_soil_infilt(kjpindex, ins, njsc, flux_infilt, qinfilt_ns, ru_infilt, check)

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    ! GLOBAL (in or inout)
    INTEGER(i_std), INTENT(in)                        :: kjpindex        !! Domain size
    INTEGER(i_std), INTENT(in)                        :: ins
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)  :: njsc            !! Index of the dominant soil textural class in the grid cell
                                                                         !!  (1-nscm, unitless)
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)    :: flux_infilt     !! Water to infiltrate
                                                                         !!  @tex $(kg m^{-2})$ @endtex

    !! 0.2 Output variables
    REAL(r_std), DIMENSION(kjpindex,nstm), INTENT(out) :: check       !! delta SM - flux (mm/dt_sechiba)
    REAL(r_std), DIMENSION(kjpindex,nstm), INTENT(out) :: ru_infilt   !! Surface runoff from soil_infilt (mm/dt_sechiba)
    REAL(r_std), DIMENSION(kjpindex,nstm), INTENT(out) :: qinfilt_ns  !! Effective infiltration flux (mm/dt_sechiba)

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                :: ji, jsl      !! Indices
    REAL(r_std), DIMENSION (kjpindex)             :: wat_inf_pot  !! infiltrable water in the layer
    REAL(r_std), DIMENSION (kjpindex)             :: wat_inf      !! infiltrated water in the layer
    REAL(r_std), DIMENSION (kjpindex)             :: dt_tmp       !! time remaining before the end of the time step
    REAL(r_std), DIMENSION (kjpindex)             :: dt_inf       !! the time it takes to complete the infiltration in the 
                                                                  !! layer 
    REAL(r_std)                                   :: k_m          !! the mean conductivity used for the saturated front
    REAL(r_std), DIMENSION (kjpindex)             :: infilt_tmp   !! infiltration rate for the considered layer 
    REAL(r_std), DIMENSION (kjpindex)             :: infilt_tot   !! total infiltration 
    REAL(r_std), DIMENSION (kjpindex)             :: flux_tmp     !! rate at which precip hits the ground

    REAL(r_std), DIMENSION(kjpindex)              :: tmci         !! total SM at beginning of routine (kg/m2)
    REAL(r_std), DIMENSION(kjpindex)              :: tmcf         !! total SM at end of routine (kg/m2)
    

!_ ================================================================================================================================

    ! If data (or coupling with GCM) was available, a parameterization for subgrid rainfall could be performed

    !! 1. We calculate the total SM at the beginning of the routine
    IF (check_cwrr2) THEN
       tmci(:) = dz(2) * ( trois*mc(:,1,ins) + mc(:,2,ins) )/huit
       DO jsl = 2,nslm-1
          tmci(:) = tmci(:) + dz(jsl) * (trois*mc(:,jsl,ins)+mc(:,jsl-1,ins))/huit &
               + dz(jsl+1) * (trois*mc(:,jsl,ins)+mc(:,jsl+1,ins))/huit
       ENDDO
       tmci(:) = tmci(:) + dz(nslm) * (trois*mc(:,nslm,ins) + mc(:,nslm-1,ins))/huit
    ENDIF

    !! 2. Infiltration process

    !! 2.1 Initialization

    DO ji = 1, kjpindex
       !! First we fill up the first layer (about 1mm) without any resistance and quasi-immediately
       wat_inf_pot(ji) = MAX((mcs(ji)-mc(ji,1,ins)) * dz(2) / deux, zero)
       wat_inf(ji) = MIN(wat_inf_pot(ji), flux_infilt(ji))
       mc(ji,1,ins) = mc(ji,1,ins) + wat_inf(ji) * deux / dz(2)
       !
    ENDDO

    !! Initialize a countdown for infiltration during the time-step and the value of potential runoff
    dt_tmp(:) = dt_sechiba / one_day
    infilt_tot(:) = wat_inf(:)
    !! Compute the rate at which water will try to infiltrate each layer
    ! flux_temp is converted here to the same unit as k_m
    flux_tmp(:) = (flux_infilt(:)-wat_inf(:)) / dt_tmp(:)

    !! 2.2 Infiltration layer by layer 
    DO jsl = 2, nslm-1
       DO ji = 1, kjpindex
          !! Infiltrability of each layer if under a saturated one
          ! This is computed by an simple arithmetic average because 
          ! the time step (30min) is not appropriate for a geometric average (advised by Haverkamp and Vauclin)
          k_m = (k(ji,jsl) + ks(njsc(ji))*kfact(jsl-1,njsc(ji))*kfact_root(ji,jsl,ins)) / deux 

          IF (ok_freeze_cwrr) THEN
             IF (temp_hydro(ji, jsl) .LT. ZeroCelsius) THEN
                k_m = k(ji,jsl)
             ENDIF
          ENDIF

          !! We compute the mean rate at which water actually infiltrate:
          ! Subgrid: Exponential distribution of k around k_m, but average p directly used 
          ! See d'Orgeval 2006, p 78, but it's not fully clear to me (AD16***)
          infilt_tmp(ji) = k_m * (un - EXP(- flux_tmp(ji) / k_m)) 

          !! From which we deduce the time it takes to fill up the layer or to end the time step...
          wat_inf_pot(ji) =  MAX((mcs(ji)-mc(ji,jsl,ins)) * (dz(jsl) + dz(jsl+1)) / deux, zero)
          IF ( infilt_tmp(ji) > min_sechiba) THEN
             dt_inf(ji) =  MIN(wat_inf_pot(ji)/infilt_tmp(ji), dt_tmp(ji))
             ! The water infiltration TIME has to limited by what is still available for infiltration.
             IF ( dt_inf(ji) * infilt_tmp(ji) > flux_infilt(ji)-infilt_tot(ji) ) THEN
                dt_inf(ji) = MAX(flux_infilt(ji)-infilt_tot(ji),zero)/infilt_tmp(ji)
             ENDIF
          ELSE
             dt_inf(ji) = dt_tmp(ji)
          ENDIF

          !! The water enters in the layer
          wat_inf(ji) = dt_inf(ji) * infilt_tmp(ji)
          ! bviously the moisture content
          mc(ji,jsl,ins) = mc(ji,jsl,ins) + &
               & wat_inf(ji) * deux / (dz(jsl) + dz(jsl+1))
          ! the time remaining before the next time step
          dt_tmp(ji) = dt_tmp(ji) - dt_inf(ji)
          ! and finally the infilt_tot (which is just used to check if there is a problem, below) 
          infilt_tot(ji) = infilt_tot(ji) + infilt_tmp(ji) * dt_inf(ji)
       ENDDO
    ENDDO

    !! 2.3 Resulting infiltration and surface runoff
    ru_infilt(:,ins) = flux_infilt(:) - infilt_tot(:)
    qinfilt_ns(:,ins) = infilt_tot(:)

    !! 3. For water conservation check: we calculate the total SM at the beginning of the routine
    !!    and export the difference with the flux
    IF (check_cwrr2) THEN
       tmcf(:) = dz(2) * ( trois*mc(:,1,ins) + mc(:,2,ins) )/huit
       DO jsl = 2,nslm-1
          tmcf(:) = tmcf(:) + dz(jsl) * (trois*mc(:,jsl,ins)+mc(:,jsl-1,ins))/huit &
               + dz(jsl+1) * (trois*mc(:,jsl,ins)+mc(:,jsl+1,ins))/huit
       ENDDO
       tmcf(:) = tmcf(:) + dz(nslm) * (trois*mc(:,nslm,ins) + mc(:,nslm-1,ins))/huit
       ! Normally, tcmf=tmci+infilt_tot
       check(:,ins) = tmcf(:)-(tmci(:)+infilt_tot(:))
    ENDIF
    
    !! 5. Local verification
    DO ji = 1, kjpindex
       IF (infilt_tot(ji) .LT. -min_sechiba .OR. infilt_tot(ji) .GT. flux_infilt(ji) + min_sechiba) THEN
          WRITE (numout,*) 'Error in the calculation of infilt tot', infilt_tot(ji)
          WRITE (numout,*) 'k, ji, jst, mc', k(ji,1:2), ji, ins, mc(ji,1,ins)
          CALL ipslerr_p(3, 'hydrol_soil_infilt', 'We will STOP now.','Error in calculation of infilt tot','')
       ENDIF
    ENDDO

  END SUBROUTINE hydrol_soil_infilt


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_soil_smooth_under_mcr
!!
!>\BRIEF        : Modifies the soil moisture profile to avoid under-residual values, 
!!                then diagnoses the points where such "excess" values remain. 
!!
!! DESCRIPTION  :
!! The "excesses" under-residual are corrected from top to bottom, by transfer of excesses 
!! to the lower layers. The reverse transfer is performed to smooth any remaining "excess" in the bottom layer.
!! If some "excess" remain afterwards, the entire soil profile is at the threshold value (mcs or mcr), 
!! and the remaining "excess" is necessarily concentrated in the top layer. 
!! This allowing diagnosing the flag is_under_mcr. 
!! Eventually, the remaining "excess" is split over the entire profile
!! 1. We calculate the total SM at the beginning of the routine
!! 2. Smoothes the profile to avoid negative values of punctual soil moisture
!! Note that we check that mc > min_sechiba in hydrol_soil 
!! 3. For water conservation check, We calculate the total SM at the beginning of the routine,
!!    and export the difference with the flux
!!
!! RECENT CHANGE(S) : 2016 by A. Ducharne
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_soil_smooth_under_mcr

  SUBROUTINE hydrol_soil_smooth_under_mcr(kjpindex, ins, njsc, is_under_mcr, check)

    !- arguments

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                         :: kjpindex     !! Domain size
    INTEGER(i_std), INTENT(in)                         :: ins          !! Soiltile index (1-nstm, unitless)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)   :: njsc         !! Index of the dominant soil textural class in grid cell 
                                                                       !! (1-nscm, unitless)    
    
    !! 0.2 Output variables

    LOGICAL, DIMENSION(kjpindex,nstm), INTENT(out)     :: is_under_mcr !! Flag diagnosing under residual soil moisture
    REAL(r_std), DIMENSION(kjpindex,nstm), INTENT(out) :: check        !! delta SM - flux

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                       :: ji,jsl
    REAL(r_std)                          :: excess
    REAL(r_std), DIMENSION(kjpindex)     :: excessji
    REAL(r_std), DIMENSION(kjpindex)     :: tmci      !! total SM at beginning of routine
    REAL(r_std), DIMENSION(kjpindex)     :: tmcf      !! total SM at end of routine

!_ ================================================================================================================================       

    !! 1. We calculate the total SM at the beginning of the routine
    IF (check_cwrr2) THEN
       tmci(:) = dz(2) * ( trois*mc(:,1,ins) + mc(:,2,ins) )/huit
       DO jsl = 2,nslm-1
          tmci(:) = tmci(:) + dz(jsl) * (trois*mc(:,jsl,ins)+mc(:,jsl-1,ins))/huit &
               + dz(jsl+1) * (trois*mc(:,jsl,ins)+mc(:,jsl+1,ins))/huit
       ENDDO
       tmci(:) = tmci(:) + dz(nslm) * (trois*mc(:,nslm,ins) + mc(:,nslm-1,ins))/huit
    ENDIF

    !! 2. Smoothes the profile to avoid negative values of punctual soil moisture

    ! 2.1 smoothing from top to bottom
    DO jsl = 1,nslm-2
       DO ji=1, kjpindex
          excess = MAX(mcr(njsc(ji))-mc(ji,jsl,ins),zero)
          mc(ji,jsl,ins) = mc(ji,jsl,ins) + excess
          mc(ji,jsl+1,ins) = mc(ji,jsl+1,ins) - excess * &
               &  (dz(jsl)+dz(jsl+1))/(dz(jsl+1)+dz(jsl+2))
       ENDDO
    ENDDO

    jsl = nslm-1
    DO ji=1, kjpindex
       excess = MAX(mcr(njsc(ji))-mc(ji,jsl,ins),zero)
       mc(ji,jsl,ins) = mc(ji,jsl,ins) + excess
       mc(ji,jsl+1,ins) = mc(ji,jsl+1,ins) - excess * &
            &  (dz(jsl)+dz(jsl+1))/dz(jsl+1)
    ENDDO

    jsl = nslm
    DO ji=1, kjpindex
       excess = MAX(mcr(njsc(ji))-mc(ji,jsl,ins),zero)
       mc(ji,jsl,ins) = mc(ji,jsl,ins) + excess
       mc(ji,jsl-1,ins) = mc(ji,jsl-1,ins) - excess * &
            &  dz(jsl)/(dz(jsl-1)+dz(jsl))
    ENDDO

    ! 2.2 smoothing from bottom to top
    DO jsl = nslm-1,2,-1
       DO ji=1, kjpindex
          excess = MAX(mcr(njsc(ji))-mc(ji,jsl,ins),zero)
          mc(ji,jsl,ins) = mc(ji,jsl,ins) + excess
          mc(ji,jsl-1,ins) = mc(ji,jsl-1,ins) - excess * &
               &  (dz(jsl)+dz(jsl+1))/(dz(jsl-1)+dz(jsl))
       ENDDO
    ENDDO

    ! 2.3 diagnoses is_under_mcr(ji), and updates the entire profile 
    ! excess > 0
    DO ji=1, kjpindex
       excessji(ji) = mask_soiltile(ji,ins) * MAX(mcr(njsc(ji))-mc(ji,1,ins),zero)
    ENDDO
    DO ji=1, kjpindex
       mc(ji,1,ins) = mc(ji,1,ins) + excessji(ji) ! then mc(1)=mcr
       is_under_mcr(ji,ins) = (excessji(ji) .GT. min_sechiba)
    ENDDO

    ! 2.4 The amount of water corresponding to excess in the top soil layer is redistributed in all soil layers
      ! -excess(ji) * dz(2) / deux donne le deficit total, negatif, en mm
      ! diviser par la profondeur totale en mm donne des delta_mc identiques en chaque couche, en mm
      ! retransformes en delta_mm par couche selon les bonnes eqs (eqs_hydrol.pdf, Eqs 13-15), puis sommes
      ! retourne bien le deficit total en mm
    DO jsl = 1, nslm
       DO ji=1, kjpindex
          mc(ji,jsl,ins) = mc(ji,jsl,ins) - excessji(ji) * dz(2) / (deux * zmaxh*mille)
       ENDDO
    ENDDO
    ! This can lead to mc(jsl) < mcr depending on the value of excess,
    ! but this is no major pb for the diffusion
    ! Yet, we need to prevent evaporation if is_under_mcr
    
    !! Note that we check that mc > min_sechiba in hydrol_soil 

    ! We just make sure that mc remains at 0 where soiltile=0
    DO jsl = 1, nslm
       DO ji=1, kjpindex
          mc(ji,jsl,ins) = mask_soiltile(ji,ins) * mc(ji,jsl,ins)
       ENDDO
    ENDDO

    !! 3. For water conservation check, We calculate the total SM at the beginning of the routine,
    !!    and export the difference with the flux
    IF (check_cwrr2) THEN
       tmcf(:) = dz(2) * ( trois*mc(:,1,ins) + mc(:,2,ins) )/huit
       DO jsl = 2,nslm-1
          tmcf(:) = tmcf(:) + dz(jsl) * (trois*mc(:,jsl,ins)+mc(:,jsl-1,ins))/huit &
               + dz(jsl+1) * (trois*mc(:,jsl,ins)+mc(:,jsl+1,ins))/huit
       ENDDO
       tmcf(:) = tmcf(:) + dz(nslm) * (trois*mc(:,nslm,ins) + mc(:,nslm-1,ins))/huit
       ! Normally, tcmf=tmci since we just redistribute the deficit 
       check(:,ins) = tmcf(:)-tmci(:)
    ENDIF
       
  END SUBROUTINE hydrol_soil_smooth_under_mcr


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_soil_smooth_over_mcs
!!
!>\BRIEF        : Modifies the soil moisture profile to avoid over-saturation values, 
!!                by putting the excess in ru_ns
!!                Thus, no point remain where such "excess" values remain (is_over_mcs becomes useless) 
!!
!! DESCRIPTION  :
!! The "excesses" over-saturation are corrected from top to bottom, by transfer of excesses 
!! to the lower layers. The reverse transfer is performed to smooth any remaining "excess" in the bottom layer.
!! If some "excess" remain afterwards, the entire soil profile is at the threshold value (mcs or mcr), 
!! and the remaining "excess" is necessarily concentrated in the top layer. 
!! Eventually, the remaining "excess" creates rudr_corr, to be added to ru_ns or dr_ns
!! 1. We calculate the total SM at the beginning of the routine
!! 2. In case of over-saturation we put the water where it is possible by smoothing
!! 3. The excess is transformed into surface runoff, with conversion from m3/m3 to kg/m2
!! 4. For water conservation checks, we calculate the total SM at the beginning of the routine,
!!    and export the difference with the flux
!!
!! RECENT CHANGE(S) : 2016 by A. Ducharne
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_soil_smooth_over_mcs

  SUBROUTINE hydrol_soil_smooth_over_mcs(kjpindex, ins, njsc, is_over_mcs, rudr_corr, check)

    !- arguments

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                           :: kjpindex        !! Domain size
    INTEGER(i_std), INTENT(in)                           :: ins             !! Soiltile index (1-nstm, unitless)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)     :: njsc            !! Index of the dominant soil textural class in grid cell 
                                                                            !! (1-nscm, unitless)
    
    !! 0.2 Output variables

    LOGICAL, DIMENSION(kjpindex), INTENT(out)            :: is_over_mcs     !! Flag diagnosing over saturated soil moisture  
    REAL(r_std), DIMENSION(kjpindex,nstm), INTENT(out)   :: check           !! delta SM - flux
    
    !! 0.3 Modified variables   
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT(inout):: rudr_corr         !! Surface runoff produced to correct excess (mm/dtstep)

    !! 0.4 Local variables

    INTEGER(i_std)                        :: ji,jsl
    REAL(r_std)                           :: excess
    REAL(r_std), DIMENSION(kjpindex)      :: tmci    !! total SM at beginning of routine
    REAL(r_std), DIMENSION(kjpindex)      :: tmcf    !! total SM at end of routine

    !_ ================================================================================================================================

    !! 1. We calculate the total SM at the beginning of the routine
    IF (check_cwrr2) THEN
       tmci(:) = dz(2) * ( trois*mc(:,1,ins) + mc(:,2,ins) )/huit
       DO jsl = 2,nslm-1
          tmci(:) = tmci(:) + dz(jsl) * (trois*mc(:,jsl,ins)+mc(:,jsl-1,ins))/huit &
               + dz(jsl+1) * (trois*mc(:,jsl,ins)+mc(:,jsl+1,ins))/huit
       ENDDO
       tmci(:) = tmci(:) + dz(nslm) * (trois*mc(:,nslm,ins) + mc(:,nslm-1,ins))/huit
    ENDIF

    !! 2. In case of over-saturation we put the water where it is possible by smoothing

    ! 2.1 smoothing from top to bottom
    DO jsl = 1, nslm-2
       DO ji=1, kjpindex
          excess = MAX(mc(ji,jsl,ins)-mcs(ji),zero)
          mc(ji,jsl,ins) = mc(ji,jsl,ins) - excess
          mc(ji,jsl+1,ins) = mc(ji,jsl+1,ins) + excess * &
               &  (dz(jsl)+dz(jsl+1))/(dz(jsl+1)+dz(jsl+2))
       ENDDO
    ENDDO

    jsl = nslm-1
    DO ji=1, kjpindex
       excess = MAX(mc(ji,jsl,ins)-mcs(ji),zero)
       mc(ji,jsl,ins) = mc(ji,jsl,ins) - excess
       mc(ji,jsl+1,ins) = mc(ji,jsl+1,ins) + excess * &
            &  (dz(jsl)+dz(jsl+1))/dz(jsl+1)
    ENDDO

    jsl = nslm
    DO ji=1, kjpindex
       excess = MAX(mc(ji,jsl,ins)-mcs(ji),zero)
       mc(ji,jsl,ins) = mc(ji,jsl,ins) - excess
       mc(ji,jsl-1,ins) = mc(ji,jsl-1,ins) + excess * &
            &  dz(jsl)/(dz(jsl-1)+dz(jsl))
    ENDDO

    ! 2.2 smoothing from bottom to top, leading  to keep most of the excess in the soil column
    DO jsl = nslm-1,2,-1
       DO ji=1, kjpindex
          excess = MAX(mc(ji,jsl,ins)-mcs(ji),zero)
          mc(ji,jsl,ins) = mc(ji,jsl,ins) - excess
          mc(ji,jsl-1,ins) = mc(ji,jsl-1,ins) + excess * &
               &  (dz(jsl)+dz(jsl+1))/(dz(jsl-1)+dz(jsl))
       ENDDO
    ENDDO

    !! 3. The excess is transformed into surface runoff, with conversion from m3/m3 to kg/m2

    DO ji=1, kjpindex
       excess = mask_soiltile(ji,ins) * MAX(mc(ji,1,ins)-mcs(ji),zero)
       mc(ji,1,ins) = mc(ji,1,ins) - excess ! then mc(1)=mcs
       rudr_corr(ji,ins) = rudr_corr(ji,ins) + excess * dz(2) / deux 
       is_over_mcs(ji) = .FALSE.
    ENDDO

    !! 4. For water conservation checks, we calculate the total SM at the beginning of the routine,
    !!    and export the difference with the flux

    IF (check_cwrr2) THEN
       tmcf(:) = dz(2) * ( trois*mc(:,1,ins) + mc(:,2,ins) )/huit
       DO jsl = 2,nslm-1
          tmcf(:) = tmcf(:) + dz(jsl) * (trois*mc(:,jsl,ins)+mc(:,jsl-1,ins))/huit &
               + dz(jsl+1) * (trois*mc(:,jsl,ins)+mc(:,jsl+1,ins))/huit
       ENDDO
       tmcf(:) = tmcf(:) + dz(nslm) * (trois*mc(:,nslm,ins) + mc(:,nslm-1,ins))/huit
       ! Normally, tcmf=tmci-rudr_corr
       check(:,ins) = tmcf(:)-(tmci(:)-rudr_corr(:,ins))
    ENDIF
    
  END SUBROUTINE hydrol_soil_smooth_over_mcs

 !! ================================================================================================================================
!! SUBROUTINE   : hydrol_soil_smooth_over_mcs2
!!
!>\BRIEF        : Modifies the soil moisture profile to avoid over-saturation values, 
!!                by putting the excess in ru_ns
!!                Thus, no point remain where such "excess" values remain (is_over_mcs becomes useless) 
!!
!! DESCRIPTION  :
!! The "excesses" over-saturation are corrected, by directly discarding the excess as rudr_corr,
!! to be added to ru_ns or dr_nsrunoff (via rudr_corr).
!! Therefore, there is no more smoothing, and this helps preventing the saturation of too many layers,
!! which leads to numerical errors with tridiag.
!! 1. We calculate the total SM at the beginning of the routine
!! 2. In case of over-saturation, we directly eliminate the excess via rudr_corr
!!    The calculation of the adjustement flux needs to account for nodes n-1 and n+1.
!! 3. For water conservation checks, we calculate the total SM at the beginning of the routine,
!!    and export the difference with the flux   
!!
!! RECENT CHANGE(S) : 2016 by A. Ducharne
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_soil_smooth_over_mcs2

  SUBROUTINE hydrol_soil_smooth_over_mcs2(kjpindex, ins, njsc, is_over_mcs, rudr_corr, check)

    !- arguments

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                           :: kjpindex        !! Domain size
    INTEGER(i_std), INTENT(in)                           :: ins             !! Soiltile index (1-nstm, unitless)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)     :: njsc            !! Index of the dominant soil textural class in grid cell 
                                                                            !! (1-nscm, unitless)
    
    !! 0.2 Output variables

    LOGICAL, DIMENSION(kjpindex), INTENT(out)            :: is_over_mcs     !! Flag diagnosing over saturated soil moisture  
    REAL(r_std), DIMENSION(kjpindex,nstm), INTENT(out)   :: check           !! delta SM - flux
    
    !! 0.3 Modified variables   
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT(inout):: rudr_corr         !! Surface runoff produced to correct excess (mm/dtstep)

    !! 0.4 Local variables

    INTEGER(i_std)                        :: ji,jsl
    REAL(r_std), DIMENSION(kjpindex,nslm) :: excess
    REAL(r_std), DIMENSION(kjpindex)      :: tmci    !! total SM at beginning of routine
    REAL(r_std), DIMENSION(kjpindex)      :: tmcf    !! total SM at end of routine

!_ ================================================================================================================================       
    !-

    !! 1. We calculate the total SM at the beginning of the routine
    IF (check_cwrr2) THEN
       tmci(:) = dz(2) * ( trois*mc(:,1,ins) + mc(:,2,ins) )/huit
       DO jsl = 2,nslm-1
          tmci(:) = tmci(:) + dz(jsl) * (trois*mc(:,jsl,ins)+mc(:,jsl-1,ins))/huit &
               + dz(jsl+1) * (trois*mc(:,jsl,ins)+mc(:,jsl+1,ins))/huit
       ENDDO
       tmci(:) = tmci(:) + dz(nslm) * (trois*mc(:,nslm,ins) + mc(:,nslm-1,ins))/huit
    ENDIF  

    !! 2. In case of over-saturation, we don't do any smoothing,
    !! but directly eliminate the excess as runoff (via rudr_corr)
    !    we correct the calculation of the adjustement flux, which needs to account for nodes n-1 and n+1  
    !    for the calculation to remain simple and accurate, we directly drain all the oversaturated mc,
    !    without transfering to lower layers        

    !! 2.1 thresholding from top to bottom, with excess defined along jsl
    DO jsl = 1, nslm
       DO ji=1, kjpindex
          excess(ji,jsl) = MAX(mc(ji,jsl,ins)-mcs(ji),zero) ! >=0
          mc(ji,jsl,ins) = mc(ji,jsl,ins) - excess(ji,jsl) ! here mc either does not change or decreases
       ENDDO
    ENDDO

    !! 2.2 To ensure conservation, this needs to be balanced by additional drainage (in kg/m2/dt)                        
    DO ji = 1, kjpindex
       rudr_corr(ji,ins) = dz(2) * ( trois*excess(ji,1) + excess(ji,2) )/huit ! top layer = initialisation  
    ENDDO
    DO jsl = 2,nslm-1 ! intermediate layers      
       DO ji = 1, kjpindex
          rudr_corr(ji,ins) = rudr_corr(ji,ins) + dz(jsl) &
               & * (trois*excess(ji,jsl)+excess(ji,jsl-1))/huit &
               & + dz(jsl+1) * (trois*excess(ji,jsl)+excess(ji,jsl+1))/huit
       ENDDO
    ENDDO
    DO ji = 1, kjpindex
       rudr_corr(ji,ins) = rudr_corr(ji,ins) + dz(nslm) &    ! bottom layer 
            & * (trois * excess(ji,nslm) + excess(ji,nslm-1))/huit
       is_over_mcs(ji) = .FALSE. 
    END DO

    !! 3. For water conservation checks, we calculate the total SM at the beginning of the routine,
    !!    and export the difference with the flux

    IF (check_cwrr2) THEN
       tmcf(:) = dz(2) * ( trois*mc(:,1,ins) + mc(:,2,ins) )/huit
       DO jsl = 2,nslm-1
          tmcf(:) = tmcf(:) + dz(jsl) * (trois*mc(:,jsl,ins)+mc(:,jsl-1,ins))/huit &
               + dz(jsl+1) * (trois*mc(:,jsl,ins)+mc(:,jsl+1,ins))/huit
       ENDDO
       tmcf(:) = tmcf(:) + dz(nslm) * (trois*mc(:,nslm,ins) + mc(:,nslm-1,ins))/huit
       ! Normally, tcmf=tmci-rudr_corr
       check(:,ins) = tmcf(:)-(tmci(:)-rudr_corr(:,ins))
    ENDIF
    
  END SUBROUTINE hydrol_soil_smooth_over_mcs2


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_soil_flux
!!
!>\BRIEF        : This subroutine diagnoses the vertical liquid water fluxes between the 
!!                different soil layers, based on each layer water budget. It also checks the
!!                corresponding water conservation (during redistribution).
!!
!! DESCRIPTION  :
!! 1. Initialize qflux from the bottom, with dr_ns
!! 2. Between layer nslm and nslm-1, by means of water budget knowing mc changes and flux at the lowest interface
!! 3. We go up, and deduct qflux(1:nslm-2), still by means of water budget
!! 4. Water balance verification: pursuing upward water budget, the flux at the surface should equal -flux_top  
!!
!! RECENT CHANGE(S) : 2016 by A. Ducharne to fit hydrol_soil
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_soil_flux

  SUBROUTINE hydrol_soil_flux(kjpindex,ins,mclint,flux_top)
    !
    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                         :: kjpindex        !! Domain size
    INTEGER(i_std), INTENT(in)                         :: ins             !! index of soil type
    REAL(r_std), DIMENSION (kjpindex,nslm), INTENT(in) :: mclint          !! mc values at the beginning of the time step
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)      :: flux_top        !! Exfiltration (bare soil evaporation minus infiltration)
    
    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                     :: jsl,ji
    REAL(r_std), DIMENSION(kjpindex)                   :: temp

    !_ ================================================================================================================================

    !- Compute the diffusion flux at every level from bottom to top (using mcl,mclint, and sink values)
    DO ji = 1, kjpindex

       !! 1. Initialize qflux from the bottom, with dr_ns
       jsl = nslm
       qflux(ji,jsl,ins) = dr_ns(ji,ins)
       !! 2. Between layer nslm and nslm-1, by means of water budget knowing mc changes and flux at the lowest interface
       !     qflux is downward
       jsl = nslm-1
       qflux(ji,jsl,ins) = qflux(ji,jsl+1,ins) & 
            &  + (mcl(ji,jsl,ins)-mclint(ji,jsl) &
            &  + trois*mcl(ji,jsl+1,ins) - trois*mclint(ji,jsl+1)) &
            &  * (dz(jsl+1)/huit) &
            &  + rootsink(ji,jsl+1,ins) 
    ENDDO

    !! 3. We go up, and deduct qflux(1:nslm-2), still by means of water budget
    ! Here, qflux(ji,1,ins) is the downward flux between the top soil layer and the 2nd one
    DO jsl = nslm-2,1,-1
       DO ji = 1, kjpindex
          qflux(ji,jsl,ins) = qflux(ji,jsl+1,ins) & 
               &  + (mcl(ji,jsl,ins)-mclint(ji,jsl) &
               &  + trois*mcl(ji,jsl+1,ins) - trois*mclint(ji,jsl+1)) &
               &  * (dz(jsl+1)/huit) &
               &  + rootsink(ji,jsl+1,ins) &
               &  + (dz(jsl+2)/huit) &
               &  * (trois*mcl(ji,jsl+1,ins) - trois*mclint(ji,jsl+1) &
               &  + mcl(ji,jsl+2,ins)-mclint(ji,jsl+2)) 
       END DO
    ENDDO
    
    !! 4. Water balance verification: pursuing upward water budget, the flux at the surface (temp) should equal -flux_top
    DO ji = 1, kjpindex
       temp(ji) =  qflux(ji,1,ins) + (dz(2)/huit) &
            &  * (trois* (mcl(ji,1,ins)-mclint(ji,1)) + (mcl(ji,2,ins)-mclint(ji,2))) &
            &  + rootsink(ji,1,ins)
    ENDDO

    ! flux_top is positive when upward, while temp is positive when downward
    DO ji = 1, kjpindex
       IF (ABS(flux_top(ji)+temp(ji)).GT. deux*min_sechiba) THEN
          WRITE(numout,*) 'Problem in the water balance, qflux computation', flux_top(ji),temp(ji)
          WRITE(numout,*) 'ji', ji, 'jsl',jsl,'ins',ins
          WRITE(numout,*) 'mclint', mclint(ji,:)
          WRITE(numout,*) 'mcl', mcl(ji,:,ins)
          WRITE (numout,*) 'rootsink', rootsink(ji,1,ins)
          CALL ipslerr_p(3, 'hydrol_soil_flux', 'We will STOP now.',&
               & 'Problem in the water balance, qflux computation','')
       ENDIF
    ENDDO

  END SUBROUTINE hydrol_soil_flux


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_soil_tridiag
!!
!>\BRIEF        This subroutine solves a set of linear equations which has a tridiagonal coefficient matrix. 
!!
!! DESCRIPTION  : It is only applied in the grid-cells where resolv(ji)=TRUE
!! 
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : mcl (global module variable)
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_soil_tridiag 

  SUBROUTINE hydrol_soil_tridiag(kjpindex,ins)

    !- arguments

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                         :: kjpindex        !! Domain size
    INTEGER(i_std), INTENT(in)                         :: ins             !! number of soil type

    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                     :: ji,jsl
    REAL(r_std), DIMENSION(kjpindex)                   :: bet
    REAL(r_std), DIMENSION(kjpindex,nslm)              :: gam

!_ ================================================================================================================================
    DO ji = 1, kjpindex

       IF (resolv(ji)) THEN
          bet(ji) = tmat(ji,1,2)
          mcl(ji,1,ins) = rhs(ji,1)/bet(ji)
       ENDIF
    ENDDO

    DO jsl = 2,nslm
       DO ji = 1, kjpindex
          
          IF (resolv(ji)) THEN

             gam(ji,jsl) = tmat(ji,jsl-1,3)/bet(ji)
             bet(ji) = tmat(ji,jsl,2) - tmat(ji,jsl,1)*gam(ji,jsl)
             mcl(ji,jsl,ins) = (rhs(ji,jsl)-tmat(ji,jsl,1)*mcl(ji,jsl-1,ins))/bet(ji)
          ENDIF

       ENDDO
    ENDDO

    DO ji = 1, kjpindex
       IF (resolv(ji)) THEN
          DO jsl = nslm-1,1,-1
             mcl(ji,jsl,ins) = mcl(ji,jsl,ins) - gam(ji,jsl+1)*mcl(ji,jsl+1,ins)
          ENDDO
       ENDIF
    ENDDO

  END SUBROUTINE hydrol_soil_tridiag


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_soil_coef
!!
!>\BRIEF        Computes coef for the linearised hydraulic conductivity 
!! k_lin=a_lin mc_lin+b_lin and the linearised diffusivity d_lin. 
!!
!! DESCRIPTION  :
!! First, we identify the interval i in which the current value of mc is located.
!! Then, we give the values of the linearized parameters to compute 
!! conductivity and diffusivity as K=a*mc+b and d.
!!
!! RECENT CHANGE(S) : Addition of the dependence to profil_froz_hydro_ns 
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_soil_coef
 
  SUBROUTINE hydrol_soil_coef(kjpindex,ins,njsc)

    IMPLICIT NONE
    !
    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                        :: kjpindex         !! Domain size
    INTEGER(i_std), INTENT(in)                        :: ins              !! Index of soil type
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)  :: njsc             !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)

    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                    :: jsl,ji,i
    REAL(r_std)                                       :: mc_ratio
    REAL(r_std)                                       :: mc_used    !! Used liquid water content
    REAL(r_std)                                       :: x,m
    
!_ ================================================================================================================================

    IF (ok_freeze_cwrr) THEN
    
       ! Calculation of liquid and frozen saturation degrees with respect to residual
       ! x=liquid saturation degree/residual=(mcl-mcr)/(mcs-mcr)
       ! 1-x=frozen saturation degree/residual=(mcfc-mcr)/(mcs-mcr) (=profil_froz_hydro)
       
       DO jsl=1,nslm
          DO ji=1,kjpindex
             
             x = 1._r_std - profil_froz_hydro_ns(ji,jsl,ins)
             
             ! mc_used is used in the calculation of hydrological properties
             ! It corresponds to a liquid mc, but the expression is different from mcl in hydrol_soil,
             ! to ensure that we get the a, b, d of the first bin when mcl<mcr
             mc_used = mcr(njsc(ji))+x*MAX((mc(ji,jsl,ins)-mcr(njsc(ji))),zero) 
             !
             ! calcul de k based on mc_liq
             !
             i= MAX(imin, MIN(imax-1, INT(imin +(imax-imin)*(mc_used-mcr(njsc(ji)))/(mcs(ji)-mcr(njsc(ji))))))
             a(ji,jsl) = a_lin(ji,jsl,i) * kfact_root(ji,jsl,ins) ! in mm/d
             b(ji,jsl) = b_lin(ji,jsl,i) * kfact_root(ji,jsl,ins) ! in mm/d
             d(ji,jsl) = d_lin(ji,jsl,i) * kfact_root(ji,jsl,ins) ! in mm^2/d
             k(ji,jsl) = kfact_root(ji,jsl,ins) * MAX(k_lin(ji,jsl,imin+1), &
                  a_lin(ji,jsl,i) * mc_used + b_lin(ji,jsl,i)) ! in mm/d
          ENDDO ! loop on grid
       ENDDO
             
    ELSE
       ! .NOT. ok_freeze_cwrr
       DO jsl=1,nslm
          DO ji=1,kjpindex 
             
             ! it is impossible to consider a mc<mcr for the binning
             mc_ratio = MAX(mc(ji,jsl,ins)-mcr(njsc(ji)), zero)/(mcs(ji)-mcr(njsc(ji)))
             
             i= MAX(MIN(INT((imax-imin)*mc_ratio)+imin , imax-1), imin)
             a(ji,jsl) = a_lin(ji,jsl,i) * kfact_root(ji,jsl,ins) ! in mm/d
             b(ji,jsl) = b_lin(ji,jsl,i) * kfact_root(ji,jsl,ins) ! in mm/d
             d(ji,jsl) = d_lin(ji,jsl,i) * kfact_root(ji,jsl,ins) ! in mm^2/d
             k(ji,jsl) = kfact_root(ji,jsl,ins) * MAX(k_lin(ji,jsl,imin+1), &
                  a_lin(ji,jsl,i) * mc(ji,jsl,ins) + b_lin(ji,jsl,i))  ! in mm/d
          END DO 
       END DO
    ENDIF
    
  END SUBROUTINE hydrol_soil_coef

!! ================================================================================================================================
!! SUBROUTINE   : hydrol_soil_froz
!!
!>\BRIEF        Computes profil_froz_hydro_ns, the fraction of frozen water in the soil layers. 
!!
!! DESCRIPTION  :
!!
!! RECENT CHANGE(S) : Created by A. Ducharne in 2016.
!!
!! MAIN OUTPUT VARIABLE(S) : profil_froz_hydro_ns
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_soil_froz
 
  SUBROUTINE hydrol_soil_froz(kjpindex,ins,njsc)

    IMPLICIT NONE
    !
    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                        :: kjpindex         !! Domain size
    INTEGER(i_std), INTENT(in)                        :: ins              !! Index of soil type
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)  :: njsc             !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)

    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                    :: jsl,ji,i
    REAL(r_std)                                       :: x,m
    REAL(r_std)                                       :: denom
    REAL(r_std),DIMENSION (kjpindex)                  :: froz_frac_moy
    REAL(r_std),DIMENSION (kjpindex)                  :: smtot_moy
    REAL(r_std),DIMENSION (kjpindex,nslm)             :: mc_ns
    
!_ ================================================================================================================================

!    ONLY FOR THE (ok_freeze_cwrr) CASE
    
       ! Calculation of liquid and frozen saturation degrees above residual moisture
       !   x=liquid saturation degree/residual=(mcl-mcr)/(mcs-mcr)
       !   1-x=frozen saturation degree/residual=(mcfc-mcr)/(mcs-mcr) (=profil_froz_hydro)
       ! It's important for the good work of the water diffusion scheme (tridiag) that the total
       ! liquid water also includes mcr, so mcl > 0 even when x=0
       
       DO jsl=1,nslm
          DO ji=1,kjpindex
             ! Van Genuchten parameter for thermodynamical calculation
             m = 1. -1./nvan(njsc(ji))
           
             IF ((.NOT. ok_thermodynamical_freezing).OR.(mc(ji,jsl, ins).LT.(mcr(njsc(ji))+min_sechiba))) THEN
                ! Linear soil freezing or soil moisture below residual
                IF (temp_hydro(ji, jsl).GE.(ZeroCelsius+fr_dT/2.)) THEN
                   x=1._r_std
                ELSE IF ( (temp_hydro(ji,jsl) .GE. (ZeroCelsius-fr_dT/2.)) .AND. &
                     (temp_hydro(ji,jsl) .LT. (ZeroCelsius+fr_dT/2.)) ) THEN 
                   x=(temp_hydro(ji, jsl)-(ZeroCelsius-fr_dT/2.))/fr_dT
                ELSE 
                   x=0._r_std
                ENDIF
             ELSE IF (ok_thermodynamical_freezing) THEN
                ! Thermodynamical soil freezing
                IF (temp_hydro(ji, jsl).GE.(ZeroCelsius+fr_dT/2.)) THEN
                   x=1._r_std
                ELSE IF ( (temp_hydro(ji,jsl) .GE. (ZeroCelsius-fr_dT/2.)) .AND. &
                     (temp_hydro(ji,jsl) .LT. (ZeroCelsius+fr_dT/2.)) ) THEN
                   ! Factor 2.2 from the PhD of Isabelle Gouttevin
                   x=MIN(((mcs(ji)-mcr(njsc(ji))) &
                        *((2.2*1000.*avan(njsc(ji))*(ZeroCelsius+fr_dT/2.-temp_hydro(ji, jsl)) &
                        *lhf/ZeroCelsius/10.)**nvan(njsc(ji))+1.)**(-m)) / &
                        (mc(ji,jsl, ins)-mcr(njsc(ji))),1._r_std)                
                ELSE
                   x=0._r_std 
                ENDIF
             ENDIF
             
             profil_froz_hydro_ns(ji,jsl,ins) = 1._r_std-x
             
             mc_ns(ji,jsl)=mc(ji,jsl,ins)/mcs(ji)

          ENDDO ! loop on grid
       ENDDO
    
       ! Applay correction on the frozen fraction
       ! Depends on two external parameters: froz_frac_corr and smtot_corr
       froz_frac_moy(:)=zero
       denom=zero
       DO jsl=1,nslm
          froz_frac_moy(:)=froz_frac_moy(:)+dh(jsl)*profil_froz_hydro_ns(:,jsl,ins)
          denom=denom+dh(jsl)
       ENDDO
       froz_frac_moy(:)=froz_frac_moy(:)/denom

       smtot_moy(:)=zero
       denom=zero
       DO jsl=1,nslm-1
          smtot_moy(:)=smtot_moy(:)+dh(jsl)*mc_ns(:,jsl)
          denom=denom+dh(jsl)
       ENDDO
       smtot_moy(:)=smtot_moy(:)/denom

       DO jsl=1,nslm
          profil_froz_hydro_ns(:,jsl,ins)=MIN(profil_froz_hydro_ns(:,jsl,ins)* &
                                              (froz_frac_moy(:)**froz_frac_corr)*(smtot_moy(:)**smtot_corr), max_froz_hydro)
       ENDDO

     END SUBROUTINE hydrol_soil_froz
     

!! ================================================================================================================================
!! SUBROUTINE   : hydrol_soil_setup
!!
!>\BRIEF        This subroutine computes the matrix coef.  
!!
!! DESCRIPTION  : None 
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : matrix coef
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE hydrol_soil_setup(kjpindex,ins)


    IMPLICIT NONE
    !
    !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                        :: kjpindex          !! Domain size
    INTEGER(i_std), INTENT(in)                        :: ins               !! index of soil type

    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std) :: jsl,ji
    REAL(r_std)                        :: temp3, temp4

!_ ================================================================================================================================
    !-we compute tridiag matrix coefficients (LEFT and RIGHT) 
    ! of the system to solve [LEFT]*mc_{t+1}=[RIGHT]*mc{t}+[add terms]: 
    ! e(nslm),f(nslm),g1(nslm) for the [left] vector
    ! and ep(nslm),fp(nslm),gp(nslm) for the [right] vector

    ! w_time=1 (in constantes_soil) indicates implicit computation for diffusion 
    temp3 = w_time*(dt_sechiba/one_day)/deux
    temp4 = (un-w_time)*(dt_sechiba/one_day)/deux

    ! Passage to arithmetic means for layer averages also in this subroutine : Aurelien 11/05/10

    !- coefficient for first layer
    DO ji = 1, kjpindex
       e(ji,1) = zero
       f(ji,1) = trois * dz(2)/huit  + temp3 &
            & * ((d(ji,1)+d(ji,2))/(dz(2))+a(ji,1))
       g1(ji,1) = dz(2)/(huit)       - temp3 &
            & * ((d(ji,1)+d(ji,2))/(dz(2))-a(ji,2))
       ep(ji,1) = zero
       fp(ji,1) = trois * dz(2)/huit - temp4 &
            & * ((d(ji,1)+d(ji,2))/(dz(2))+a(ji,1))
       gp(ji,1) = dz(2)/(huit)       + temp4 &
            & * ((d(ji,1)+d(ji,2))/(dz(2))-a(ji,2))
    ENDDO

    !- coefficient for medium layers

    DO jsl = 2, nslm-1
       DO ji = 1, kjpindex
          e(ji,jsl) = dz(jsl)/(huit)                        - temp3 &
               & * ((d(ji,jsl)+d(ji,jsl-1))/(dz(jsl))+a(ji,jsl-1))

          f(ji,jsl) = trois * (dz(jsl)+dz(jsl+1))/huit  + temp3 &
               & * ((d(ji,jsl)+d(ji,jsl-1))/(dz(jsl)) + &
               & (d(ji,jsl)+d(ji,jsl+1))/(dz(jsl+1)) )

          g1(ji,jsl) = dz(jsl+1)/(huit)                     - temp3 &
               & * ((d(ji,jsl)+d(ji,jsl+1))/(dz(jsl+1))-a(ji,jsl+1))

          ep(ji,jsl) = dz(jsl)/(huit)                       + temp4 &
               & * ((d(ji,jsl)+d(ji,jsl-1))/(dz(jsl))+a(ji,jsl-1))

          fp(ji,jsl) = trois * (dz(jsl)+dz(jsl+1))/huit - temp4 &
               & * ( (d(ji,jsl)+d(ji,jsl-1))/(dz(jsl)) + &
               & (d(ji,jsl)+d(ji,jsl+1))/(dz(jsl+1)) )

          gp(ji,jsl) = dz(jsl+1)/(huit)                     + temp4 &
               & *((d(ji,jsl)+d(ji,jsl+1))/(dz(jsl+1))-a(ji,jsl+1))
       ENDDO
    ENDDO

    !- coefficient for last layer
    DO ji = 1, kjpindex
       e(ji,nslm) = dz(nslm)/(huit)        - temp3 &
            & * ((d(ji,nslm)+d(ji,nslm-1)) /(dz(nslm))+a(ji,nslm-1))
       f(ji,nslm) = trois * dz(nslm)/huit  + temp3 &
            & * ((d(ji,nslm)+d(ji,nslm-1)) / (dz(nslm)) &
            & -a(ji,nslm)*(un-deux*free_drain_coef(ji,ins)))
       g1(ji,nslm) = zero
       ep(ji,nslm) = dz(nslm)/(huit)       + temp4 &
            & * ((d(ji,nslm)+d(ji,nslm-1)) /(dz(nslm))+a(ji,nslm-1))
       fp(ji,nslm) = trois * dz(nslm)/huit - temp4 &
            & * ((d(ji,nslm)+d(ji,nslm-1)) /(dz(nslm)) &
            & -a(ji,nslm)*(un-deux*free_drain_coef(ji,ins)))
       gp(ji,nslm) = zero
    ENDDO

  END SUBROUTINE hydrol_soil_setup

  
!! ================================================================================================================================
!! SUBROUTINE   : hydrol_split_soil
!!
!>\BRIEF        Splits 2d variables into 3d variables, per soiltile (_ns suffix), at the beginning of hydrol
!!              At this stage, the forcing fluxes to hydrol are transformed from grid-cell averages 
!!              to mean fluxes over vegtot=sum(soiltile)  
!!
!! DESCRIPTION  :
!! 1. Split 2d variables into 3d variables, per soiltile
!! 1.1 Throughfall
!! 1.2 Bare soil evaporation
!! 1.2.1 vevapnu_old
!! 1.2.2 ae_ns new
!! 1.3 transpiration
!! 1.4 root sink
!! 2. Verification: Check if the deconvolution is correct and conserves the fluxes
!! 2.1 precisol 
!! 2.2 ae_ns and evapnu
!! 2.3 transpiration
!! 2.4 root sink
!!
!! RECENT CHANGE(S) : 2016 by A. Ducharne to match the simplification of hydrol_soil
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_split_soil

  SUBROUTINE hydrol_split_soil (kjpindex, veget_max, soiltile, vevapnu, vevapnu_pft, transpir, humrel, evap_bare_lim, tot_bare_soil)
    ! 
    ! interface description

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                               :: kjpindex
    REAL(r_std), DIMENSION (kjpindex, nvm), INTENT(in)       :: veget_max        !! max Vegetation map 
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT (in)      :: soiltile         !! Fraction of each soiltile within vegtot (0-1, unitless)
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)           :: vevapnu          !! Bare soil evaporation
    REAL(r_std), DIMENSION (kjpindex, nvm), INTENT (in)      :: vevapnu_pft      !! Bare soil evaporation
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (in)       :: transpir         !! Transpiration
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (in)       :: humrel           !! Relative humidity
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: evap_bare_lim    !!   
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: tot_bare_soil    !! Total evaporating bare soil fraction 

    !! 0.4 Local variables

    INTEGER(i_std)                                :: ji, jv, jsl, jst
    REAL(r_std), DIMENSION (kjpindex)             :: vevapnu_old
    REAL(r_std), DIMENSION (kjpindex,nstm)        :: vevapnu_ns      !! Bare soil evaporation
    REAL(r_std), DIMENSION (kjpindex)             :: tmp_check1
    REAL(r_std), DIMENSION (kjpindex)             :: tmp_check2
    REAL(r_std), DIMENSION (kjpindex,nstm)        :: tmp_check3
    LOGICAL                                       :: error=.FALSE. !! If true, exit in the end of subroutine

!_ ================================================================================================================================
    
    !! 1. Split 2d variables into 3d variables, per soiltile
    
    ! Reminders:
    !  corr_veg_soil(:,nvm,nstm) = PFT fraction per soiltile in each grid-cell
    !      corr_veg_soil(ji,jv,jst)=veget_max(ji,jv)/soiltile(ji,jst) 
    !  soiltile(:,nstm) = fraction of vegtot covered by each soiltile (0-1, unitless) 
    !  vegtot(:) = total fraction of grid-cell covered by PFTs (fraction with bare soil + vegetation)
    !  veget_max(:,nvm) = PFT fractions of vegtot+frac_nobio 
    !  veget(:,nvm) =  fractions (of vegtot+frac_nobio) covered by vegetation in each PFT 
    !       BUT veget(:,1)=veget_max(:,1) 
    !  frac_bare(:,nvm) = fraction (of veget_max) with bare soil in each PFT
    !  tot_bare_soil(:) = fraction of grid mesh covered by all bare soil (=SUM(frac_bare*veget_max))
    !  frac_bare_ns(:,nstm) = evaporating bare soil fraction (of vegtot) per soiltile (defined in hydrol_vegupd)
    
    !! 1.1 Throughfall
    ! Transformation from precisol (flux from PFT jv in m2 of grid-mesh)
    ! to  precisol_ns (flux from contributing PFTs with another unit, in m2 of soiltile)
    precisol_ns(:,:)=zero
    DO jv=1,nvm
       DO ji=1,kjpindex
          jst=pref_soil_veg(jv)
          IF((veget_max(ji,jv).GT.min_sechiba) .AND. ((soiltile(ji,jst)*vegtot(ji)) .GT. min_sechiba)) THEN
             precisol_ns(ji,jst) = precisol_ns(ji,jst) + &
                     precisol(ji,jv) / (soiltile(ji,jst)*vegtot(ji))                
          ENDIF
       END DO
    END DO
    
    !! 1.2 Bare soil evaporation

    vevapnu_ns(:,:) = zero
    DO jv = 1,nvm
        DO jst = 1,nstm
            DO ji=1,kjpindex
                IF (veget_max(ji,jv) .GT. min_sechiba) THEN 
                    vevapnu_ns(ji,jst) = vevapnu_ns(ji,jst) + vevapnu_pft(ji,jv)* &
                            & vegetmax_soil(ji,jv,jst) / vegtot(ji) / veget_max(ji,jv)
                ENDIF
            ENDDO
        ENDDO
    ENDDO
    
    !! 1.2.1 vevapnu_old
! AD16*** vevapnu_old ne sert que pour le split suivant de vevapnu (issu de enerbil) en ae_ns pour hydrol_soil
!           mais il ne semble y avoir aucune bonne raison de contraindre ae_ns en fonction de vevapnu_old
    vevapnu_old(:)=zero
    DO jst=1,nstm
       DO ji=1,kjpindex
          IF ( vegtot(ji) .GT. min_sechiba) THEN
             vevapnu_old(ji)=vevapnu_old(ji)+ &
                  & ae_ns(ji,jst)*soiltile(ji,jst)*vegtot(ji)
          ENDIF
       END DO
    END DO
    
    !! 1.2.2 ae_ns new
! AD16*** les lignes ci-dessous sont excessivement compliquees et ne garantissent pas que ae_ns = 0 si evap_bare_lim=0 
!           c'est notamment le cas pour les 3emes et 6emes conditions
    DO jst=1,nstm
       DO ji=1,kjpindex
          IF (vevapnu_old(ji).GT.min_sechiba) THEN   
             IF(evap_bare_lim(ji).GT.min_sechiba) THEN       
                ae_ns(ji,jst) = vevapnu(ji) * evap_bare_lim_ns(ji,jst)/evap_bare_lim(ji) 
             ELSE
                IF(vevapnu_old(ji).GT.min_sechiba) THEN  
                   ae_ns(ji,jst)=ae_ns(ji,jst) * vevapnu(ji)/vevapnu_old(ji) ! 3ème condition
                ELSE
                   ae_ns(ji,jst)=zero
                ENDIF
             ENDIF
          ELSEIF(frac_bare_ns(ji,jst).GT.min_sechiba) THEN
             IF(evap_bare_lim(ji).GT.min_sechiba) THEN  
                ae_ns(ji,jst) = vevapnu(ji) * evap_bare_lim_ns(ji,jst)/evap_bare_lim(ji)
             ELSE
                IF(tot_bare_soil(ji).GT.min_sechiba) THEN  
                   ae_ns(ji,jst) = vevapnu(ji) * frac_bare_ns(ji,jst)/tot_bare_soil(ji) ! 6ème condition
                ELSE
                   ae_ns(ji,jst) = zero
                ENDIF
             ENDIF
          ENDIF
       END DO
    END DO
! ADNV27072016: we believe the following block should be used (tests needed before committ, since AD16*** had pb with it)    
!!$    ! given the definition of evap_bare_lim, it leads to sum(ae_ns(ji,jst)*soiltile(ji,jst)*vegtot(ji))=vevapnu(ji)
!!$    ae_ns(:,:)=zero
!!$    DO jst=1,nstm
!!$       DO ji=1,kjpindex
!!$          IF(evap_bare_lim(ji).GT.min_sechiba) THEN       
!!$             ae_ns(ji,jst) = vevapnu(ji) * evap_bare_lim_ns(ji,jst)/evap_bare_lim(ji)
!            ELSE
!               ae_ns(ji,jst) = zero
!!$          ENDIF
!!$       ENDDO
!!$    ENDDO
    
    !! 1.3 transpiration
    ! Transformation from transpir (flux from PFT jv in m2 of grid-mesh)
    ! to tr_ns (flux from contributing PFTs with another unit, in m2 of soiltile)
    ! To do next: simplify the use of humrelv(ji,jv,jst) /humrel(ji,jv), since both are equal
    tr_ns(:,:)=zero
    DO jv=1,nvm
       jst=pref_soil_veg(jv)
       DO ji=1,kjpindex
          IF ((humrel(ji,jv).GT.min_sechiba) .AND. ((soiltile(ji,jst)*vegtot(ji)) .GT.min_sechiba))THEN 
             tr_ns(ji,jst)= tr_ns(ji,jst) &
                  + transpir(ji,jv) * (humrelv(ji,jv,jst) / humrel(ji,jv)) &
                  / (soiltile(ji,jst)*vegtot(ji))
                     
                ! xuhui 20151217
                ! tr_ns(ji,jst)=tr_ns(ji,jst)+ vegetmax_soil(ji,jv,jst)*humrelv(ji,jv,jst)* & 
                !       & transpir(ji,jv) / humrel(ji,jv) / veget_max(ji,jv)
             ENDIF
       END DO
    END DO

    !! 1.4 root sink
    ! Transformation from transpir (flux from PFT jv in m2 of grid-mesh)
    ! to root_sink (flux from contributing PFTs and soil layer with another unit, in m2 of soiltile)
    rootsink(:,:,:)=zero
    DO jv=1,nvm
       jst=pref_soil_veg(jv)
       DO jsl=1,nslm
          DO ji=1,kjpindex
             IF ((humrel(ji,jv).GT.min_sechiba) .AND. ((soiltile(ji,jst)*vegtot(ji)) .GT.min_sechiba)) THEN 
                rootsink(ji,jsl,jst) = rootsink(ji,jsl,jst) &
                        + transpir(ji,jv) * (us(ji,jv,jst,jsl) / humrel(ji,jv)) &
                        / (soiltile(ji,jst)*vegtot(ji))                     
                   ! rootsink(ji,1,jst)=0 as us(ji,jv,jst,1)=0
             END IF
          END DO
       END DO
    END DO


    !!! ADNV270716 *** we are here

    !! 2. Verification: Check if the deconvolution is correct and conserves the fluxes (grid-cell average)

    IF (check_cwrr) THEN

       !! 2.1 precisol 

       tmp_check1(:)=zero
       DO jst=1,nstm
          DO ji=1,kjpindex
             tmp_check1(ji)=tmp_check1(ji) + precisol_ns(ji,jst)*soiltile(ji,jst)*vegtot(ji)
          END DO
       END DO
       
       tmp_check2(:)=zero  
       DO jv=1,nvm
          DO ji=1,kjpindex
             tmp_check2(ji)=tmp_check2(ji) + precisol(ji,jv)
          END DO
       END DO

       DO ji=1,kjpindex   
          IF(ABS(tmp_check1(ji) - tmp_check2(ji)).GT.allowed_err) THEN
             WRITE(numout,*) 'PRECISOL SPLIT FALSE:ji=',ji,tmp_check1(ji),tmp_check2(ji)
             WRITE(numout,*) 'err',ABS(tmp_check1(ji)- tmp_check2(ji))
             WRITE(numout,*) 'vegtot',vegtot(ji)
             DO jv=1,nvm
                WRITE(numout,'(a,i2.2,"|",F13.4,"|",F13.4,"|",3(F9.6))') &
                     'jv,veget_max, precisol, vegetmax_soil ', &
                     jv,veget_max(ji,jv),precisol(ji,jv),vegetmax_soil(ji,jv,:)
             END DO
             DO jst=1,nstm
                WRITE(numout,*) 'jst,precisol_ns',jst,precisol_ns(ji,jst)
                WRITE(numout,*) 'soiltile', soiltile(ji,jst)
             END DO
             error=.TRUE.
             CALL ipslerr_p(2, 'hydrol_split_soil', 'We will STOP in the end of this subroutine.',&
                  & 'check_CWRR','PRECISOL SPLIT FALSE')
          ENDIF
       END DO
       
       !! 2.2 ae_ns and evapnu

       tmp_check1(:)=zero
       DO jst=1,nstm
          DO ji=1,kjpindex
             tmp_check1(ji)=tmp_check1(ji) + ae_ns(ji,jst)*soiltile(ji,jst)*vegtot(ji)
          END DO
       END DO

       DO ji=1,kjpindex   

          IF(ABS(tmp_check1(ji) - vevapnu(ji)).GT.allowed_err) THEN
             WRITE(numout,*) 'VEVAPNU SPLIT FALSE:ji, Sum(ae_ns), vevapnu =',ji,tmp_check1(ji),vevapnu(ji)
             WRITE(numout,*) 'err',ABS(tmp_check1(ji)- vevapnu(ji))
             WRITE(numout,*) 'ae_ns',ae_ns(ji,:)
             WRITE(numout,*) 'vegtot',vegtot(ji)
             WRITE(numout,*) 'evap_bare_lim, evap_bare_lim_ns',evap_bare_lim(ji), evap_bare_lim_ns(ji,:)
             WRITE(numout,*) 'tot_bare_soil,frac_bare_ns',tot_bare_soil(ji),frac_bare_ns(ji,:)
             WRITE(numout,*) 'vevapnu_old',vevapnu_old(ji)
             DO jst=1,nstm
                WRITE(numout,*) 'jst,ae_ns',jst,ae_ns(ji,jst)
                WRITE(numout,*) 'soiltile', soiltile(ji,jst)
                WRITE(numout,*) 'veget_max/vegtot/soiltile', veget_max(ji,:)/vegtot(ji)/soiltile(ji,jst)
                WRITE(numout,*) "vegetmax_soil",vegetmax_soil(ji,:,jst)
             END DO
             error=.TRUE.
             CALL ipslerr_p(2, 'hydrol_split_soil', 'We will STOP in the end of this subroutine.',&
                  & 'check_CWRR','VEVAPNU SPLIT FALSE')
          ENDIF
       ENDDO

    !! 2.3 transpiration

       tmp_check1(:)=zero
       DO jst=1,nstm
          DO ji=1,kjpindex
             tmp_check1(ji)=tmp_check1(ji) + tr_ns(ji,jst)*soiltile(ji,jst)*vegtot(ji)
          END DO
       END DO
       
       tmp_check2(:)=zero  
       DO jv=1,nvm
          DO ji=1,kjpindex
             tmp_check2(ji)=tmp_check2(ji) + transpir(ji,jv)
          END DO
       END DO

       DO ji=1,kjpindex   
          IF(ABS(tmp_check1(ji)- tmp_check2(ji)).GT.allowed_err) THEN
             WRITE(numout,*) 'TRANSPIR SPLIT FALSE:ji=',ji,tmp_check1(ji),tmp_check2(ji)
             WRITE(numout,*) 'err',ABS(tmp_check1(ji)- tmp_check2(ji))
             WRITE(numout,*) 'vegtot',vegtot(ji)
             DO jv=1,nvm
                WRITE(numout,*) 'jv,veget_max, transpir',jv,veget_max(ji,jv),transpir(ji,jv)
                DO jst=1,nstm
                   WRITE(numout,*) 'vegetmax_soil:ji,jv,jst',ji,jv,jst,vegetmax_soil(ji,jv,jst)
                END DO
             END DO
             DO jst=1,nstm
                WRITE(numout,*) 'jst,tr_ns',jst,tr_ns(ji,jst)
                WRITE(numout,*) 'soiltile', soiltile(ji,jst)
             END DO
             error=.TRUE.
             CALL ipslerr_p(2, 'hydrol_split_soil', 'We will STOP in the end of this subroutine.',&
                  & 'check_CWRR','TRANSPIR SPLIT FALSE')
          ENDIF

       END DO

    !! 2.4 root sink

       tmp_check3(:,:)=zero
       DO jst=1,nstm
          DO jsl=1,nslm
             DO ji=1,kjpindex
                tmp_check3(ji,jst)=tmp_check3(ji,jst) + rootsink(ji,jsl,jst)
             END DO
          END DO
       ENDDO

       DO jst=1,nstm
          DO ji=1,kjpindex
             IF(ABS(tmp_check3(ji,jst) - tr_ns(ji,jst)).GT.allowed_err) THEN
                WRITE(numout,*) 'ROOTSINK SPLIT FALSE:ji,jst=', ji,jst,&
                     & tmp_check3(ji,jst),tr_ns(ji,jst)
                WRITE(numout,*) 'err',ABS(tmp_check3(ji,jst)- tr_ns(ji,jst))
                WRITE(numout,*) 'HUMREL(jv=1:13)',humrel(ji,:)
                WRITE(numout,*) 'TRANSPIR',transpir(ji,:)
                DO jv=1,nvm 
                   WRITE(numout,*) 'jv=',jv,'us=',us(ji,jv,jst,:)
                ENDDO
                error=.TRUE.
                CALL ipslerr_p(2, 'hydrol_split_soil', 'We will STOP in the end of this subroutine.',&
                  & 'check_CWRR','ROOTSINK SPLIT FALSE')
             ENDIF
          END DO
       END DO

    ENDIF ! end of check_cwrr

!! Exit if error was found previously in this subroutine
    IF ( error ) THEN
       WRITE(numout,*) 'One or more errors have been detected in hydrol_split_soil. Model stops.'
       CALL ipslerr_p(3, 'hydrol_split_soil', 'We will STOP now.',&
                  & 'One or several fatal errors were found previously.','')
    END IF

  END SUBROUTINE hydrol_split_soil
  

!! ================================================================================================================================
!! SUBROUTINE   : hydrol_diag_soil
!!
!>\BRIEF        Calculates diagnostic variables at the grid-cell scale
!!
!! DESCRIPTION  :
!! - 1. Apply mask_soiltile
!! - 2. Sum 3d variables in 2d variables with fraction of vegetation per soil type
!!
!! RECENT CHANGE(S) : 2016 by A. Ducharne for the claculation of shumdiag_perma
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_diag_soil

  SUBROUTINE hydrol_diag_soil (kjpindex, veget_max, soiltile, njsc, runoff, drainage, &
       & evapot, vevapnu, returnflow, reinfiltration, irrigation, &
       & shumdiag,shumdiag_perma, k_litt, litterhumdiag, humrel, vegstress, drysoil_frac, tot_melt) !pss:+

    ! 
    ! interface description

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    ! input scalar 
    INTEGER(i_std), INTENT(in)                               :: kjpindex 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)        :: veget_max       !! Max. vegetation type
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)         :: njsc            !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT (in)      :: soiltile        !! Fraction of each soil tile within vegtot (0-1, unitless)
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: evapot          !! 
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: returnflow      !! Water returning to the deep reservoir
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: reinfiltration  !! Water returning to the top of the soil
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: irrigation      !! Water from irrigation
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: tot_melt        !!

    !! 0.2 Output variables

    REAL(r_std), DIMENSION (kjpindex), INTENT (out)          :: drysoil_frac    !! Function of litter wetness
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: runoff          !! complete runoff
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: drainage        !! Drainage
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (out)      :: shumdiag        !! relative soil moisture
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (out)      :: shumdiag_perma  !! Percent of porosity filled with water (mc/mcs) used for the thermal computations
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: k_litt          !! litter cond.
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: litterhumdiag   !! litter humidity
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)       :: humrel          !! Relative humidity
    REAL(r_std), DIMENSION (kjpindex, nvm), INTENT(out)      :: vegstress       !! Veg. moisture stress (only for vegetation growth)

    !! 0.3 Modified variables 

    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)         :: vevapnu         !!

    !! 0.4 Local variables

    INTEGER(i_std)                                           :: ji, jv, jsl, jst, i
    REAL(r_std), DIMENSION (kjpindex)                        :: mask_vegtot
    REAL(r_std)                                              :: k_tmp, tmc_litter_ratio

!_ ================================================================================================================================
    !
    ! Put the prognostics variables of soil to zero if soiltype is zero

    !! 1. Apply mask_soiltile
    
    DO jst=1,nstm 
       IF (ok_freeze_cwrr) THEN
           CALL hydrol_soil_coef(kjpindex,jst,njsc)
       ENDIF
       DO ji=1,kjpindex

             ae_ns(ji,jst) = ae_ns(ji,jst) * mask_soiltile(ji,jst)
             dr_ns(ji,jst) = dr_ns(ji,jst) * mask_soiltile(ji,jst)
             ru_ns(ji,jst) = ru_ns(ji,jst) * mask_soiltile(ji,jst)
             tmc(ji,jst) =  tmc(ji,jst) * mask_soiltile(ji,jst)

             DO jv=1,nvm
                humrelv(ji,jv,jst) = humrelv(ji,jv,jst) * mask_soiltile(ji,jst)
                DO jsl=1,nslm
                   us(ji,jv,jst,jsl) = us(ji,jv,jst,jsl)  * mask_soiltile(ji,jst)
                END DO
             END DO

             DO jsl=1,nslm          
                mc(ji,jsl,jst) = mc(ji,jsl,jst)  * mask_soiltile(ji,jst)
             END DO

       END DO
    END DO

    runoff(:) = zero
    drainage(:) = zero
    humtot(:) = zero
    shumdiag(:,:)= zero
    shumdiag_perma(:,:)=zero
    k_litt(:) = zero
    litterhumdiag(:) = zero
    tmc_litt_dry_mea(:) = zero
    tmc_litt_wet_mea(:) = zero
    tmc_litt_mea(:) = zero
    humrel(:,:) = zero
    vegstress(:,:) = zero
    IF (ok_freeze_cwrr) THEN
       profil_froz_hydro(:,:)=zero ! initialisation for the mean of profil_froz_hydro_ns
    ENDIF
    
    !! 2. Sum 3d variables in 2d variables with fraction of vegetation per soil type

    DO ji = 1, kjpindex
       mask_vegtot(ji) = 0
       IF(vegtot(ji) .GT. min_sechiba) THEN
          mask_vegtot(ji) = 1
       ENDIF
    END DO
    
    DO ji = 1, kjpindex 
       ! Here we weight ae_ns by the fraction of bare evaporating soil. 
       ! This is given by frac_bare_ns, taking into account bare soil under vegetation
       ae_ns(ji,:) = mask_vegtot(ji) * ae_ns(ji,:) * frac_bare_ns(ji,:)
    END DO

    ! We average the values of each soiltile and multiply by vegtot to transform to a grid-cell mean
    DO jst = 1, nstm
       DO ji = 1, kjpindex 
          drainage(ji) = mask_vegtot(ji) * (drainage(ji) + vegtot(ji)*soiltile(ji,jst) * dr_ns(ji,jst))
          runoff(ji) = mask_vegtot(ji) *  (runoff(ji) +   vegtot(ji)*soiltile(ji,jst) * ru_ns(ji,jst)) &
               &   + (1 - mask_vegtot(ji)) * (tot_melt(ji) + irrigation(ji) + returnflow(ji) + reinfiltration(ji))
          humtot(ji) = mask_vegtot(ji) * (humtot(ji) + vegtot(ji)*soiltile(ji,jst) * tmc(ji,jst)) 
          IF (ok_freeze_cwrr) THEN 
             !  profil_froz_hydro_ns comes from hydrol_soil, to remain the same as in the prognotic loop
             profil_froz_hydro(ji,:)=mask_vegtot(ji) * &
                  (profil_froz_hydro(ji,:) + vegtot(ji)*soiltile(ji,jst) * profil_froz_hydro_ns(ji,:, jst))
          ENDIF
       END DO
    END DO

    ! we add the excess of snow sublimation to vevapnu
    ! - because vevapsno is modified in hydrol_snow if subsinksoil
    ! - it is multiplied by vegtot because it is devided by 1-tot_frac_nobio at creation in hydrol_snow

    DO ji = 1,kjpindex
       vevapnu(ji) = vevapnu (ji) + subsinksoil(ji)*vegtot(ji)
    END DO

    DO jst=1,nstm
       DO jv=1,nvm
          DO ji=1,kjpindex
             IF(veget_max(ji,jv).GT.min_sechiba) THEN
                vegstress(ji,jv)=vegstress(ji,jv)+vegstressv(ji,jv,jst)
                vegstress(ji,jv)= MAX(vegstress(ji,jv),zero)
             ENDIF
          END DO
       END DO
    END DO

    DO jst=1,nstm
       DO jv=1,nvm
          DO ji=1,kjpindex
             humrel(ji,jv)=humrel(ji,jv)+humrelv(ji,jv,jst)
             humrel(ji,jv)=MAX(humrel(ji,jv),zero)
          END DO
       END DO
    END DO

    !! Litter... the goal is to calculate drysoil_frac, to calculate the albedo in condveg
    ! In condveg, drysoil_frac serve to calculate the albedo of drysoil, excluding the nobio contribution which is further added
    ! In conclusion, we calculate drysoil_frac based on moisture averages restricted to the soiltile (no multiplication by vegtot)
    ! BUT THIS IS NOT USED ANYMORE WITH THE NEW BACKGROUNG ALBEDO
    !! k_litt is calculated here as a grid-cell average (for consistency with drainage)   
    !! litterhumdiag, like shumdiag, is averaged over the soiltiles for transmission to stomate
    DO jst=1,nstm       
       DO ji=1,kjpindex
          ! We compute here a mean k for the 'litter' used for reinfiltration from floodplains of ponds        
          IF ( tmc_litter(ji,jst) < tmc_litter_res(ji,jst)) THEN
             i = imin
          ELSE
             tmc_litter_ratio = (tmc_litter(ji,jst)-tmc_litter_res(ji,jst)) / &
                  & (tmc_litter_sat(ji,jst)-tmc_litter_res(ji,jst))
             i= MAX(MIN(INT((imax-imin)*tmc_litter_ratio)+imin, imax-1), imin)
          ENDIF       
          k_tmp = MAX(k_lin(ji,1,i)*ks(njsc(ji)), zero)
          k_litt(ji) = k_litt(ji) + vegtot(ji)*soiltile(ji,jst) * SQRT(k_tmp) ! grid-cell average
       ENDDO      
       DO ji=1,kjpindex
          litterhumdiag(ji) = litterhumdiag(ji) + &
               & soil_wet_litter(ji,jst) * soiltile(ji,jst)

          tmc_litt_wet_mea(ji) =  tmc_litt_wet_mea(ji) + & 
               & tmc_litter_awet(ji,jst)* soiltile(ji,jst)

          tmc_litt_dry_mea(ji) = tmc_litt_dry_mea(ji) + &
               & tmc_litter_adry(ji,jst) * soiltile(ji,jst) 

          tmc_litt_mea(ji) = tmc_litt_mea(ji) + &
               & tmc_litter(ji,jst) * soiltile(ji,jst) 
       ENDDO
    ENDDO
    
    DO ji=1,kjpindex
       IF ( tmc_litt_wet_mea(ji) - tmc_litt_dry_mea(ji) > zero ) THEN
          drysoil_frac(ji) = un + MAX( MIN( (tmc_litt_dry_mea(ji) - tmc_litt_mea(ji)) / &
               & (tmc_litt_wet_mea(ji) - tmc_litt_dry_mea(ji)), zero), - un)
       ELSE
          drysoil_frac(ji) = zero
       ENDIF
    END DO
    
    ! Calculate soilmoist, as a function of total water content (mc)
    ! We average the values of each soiltile and multiply by vegtot to transform to a grid-cell mean
    soilmoist(:,:) = zero
    DO jst=1,nstm
       DO ji=1,kjpindex
             soilmoist(ji,1) = soilmoist(ji,1) + soiltile(ji,jst) * &
                  dz(2) * ( trois*mc(ji,1,jst) + mc(ji,2,jst) )/huit
             DO jsl = 2,nslm-1
                soilmoist(ji,jsl) = soilmoist(ji,jsl) + soiltile(ji,jst) * &
                     ( dz(jsl) * (trois*mc(ji,jsl,jst)+mc(ji,jsl-1,jst))/huit &
                     + dz(jsl+1) * (trois*mc(ji,jsl,jst)+mc(ji,jsl+1,jst))/huit )
             END DO
             soilmoist(ji,nslm) = soilmoist(ji,nslm) + soiltile(ji,jst) * &
                  dz(nslm) * (trois*mc(ji,nslm,jst) + mc(ji,nslm-1,jst))/huit
       END DO
    END DO
    DO ji=1,kjpindex
       soilmoist(ji,:) = soilmoist(ji,:) * vegtot(ji) ! conversion to grid-cell average
    ENDDO

    soilmoist_liquid(:,:) = zero
    DO jst=1,nstm
       DO ji=1,kjpindex
          soilmoist_liquid(ji,1) = soilmoist_liquid(ji,1) + soiltile(ji,jst) * &
               dz(2) * ( trois*mcl(ji,1,jst) + mcl(ji,2,jst) )/huit
          DO jsl = 2,nslm-1
             soilmoist_liquid(ji,jsl) = soilmoist_liquid(ji,jsl) + soiltile(ji,jst) * &
                  ( dz(jsl) * (trois*mcl(ji,jsl,jst)+mcl(ji,jsl-1,jst))/huit &
                  + dz(jsl+1) * (trois*mcl(ji,jsl,jst)+mcl(ji,jsl+1,jst))/huit )
          END DO
          soilmoist_liquid(ji,nslm) = soilmoist_liquid(ji,nslm) + soiltile(ji,jst) * &
               dz(nslm) * (trois*mcl(ji,nslm,jst) + mcl(ji,nslm-1,jst))/huit
       ENDDO
    ENDDO
    DO ji=1,kjpindex
        soilmoist_liquid(ji,:) = soilmoist_liquid(ji,:) * vegtot_old(ji) ! grid cell average 
    ENDDO
    
    
    ! Shumdiag: we start from soil_wet_ns, change the range over which the relative moisture is calculated,
    ! then do a spatial average, excluding the nobio fraction on which stomate doesn't act
    DO jst=1,nstm      
       DO jsl=1,nslm
          DO ji=1,kjpindex
             shumdiag(ji,jsl) = shumdiag(ji,jsl) + soil_wet_ns(ji,jsl,jst) * soiltile(ji,jst) * &
                               ((mcs(ji)-mcw(ji))/(mcfc(ji)-mcw(ji)))
             shumdiag(ji,jsl) = MAX(MIN(shumdiag(ji,jsl), un), zero) 
          ENDDO
       ENDDO
    ENDDO
    
    ! Shumdiag_perma is based on soilmoist / moisture at saturation in the layer
    ! Her we start from grid averages by hydrol soil layer and transform it to the diag levels
    ! We keep a grid-cell average, like for all variables transmitted to ok_freeze

    ! in MICT, shumdiag_perma is used to calculate the "moisture modifier" in soil carbon decomposition rate. Global tests show that
    ! the previous calculation, namely mc/mcs, yielded too small soil carbon stock. So we modify it. (To be discussed in the future)
    DO jsl=1,nslm              
       DO ji=1,kjpindex
          ! note: if use_refSOC_hydrol is not activated, mcs(ji)=mcs_mineral(njsc(ji)) 
          shumdiag_perma(ji,jsl) = soilmoist(ji,jsl) / dh(jsl)/mcs(ji)*mcs_mineral(njsc(ji))
          shumdiag_perma(ji,jsl) = MAX(MIN(shumdiag_perma(ji,jsl), un), zero) 
       ENDDO
    ENDDO
    
  END SUBROUTINE hydrol_diag_soil  


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_alma 
!!
!>\BRIEF        This routine computes the changes in soil moisture and interception storage for the ALMA outputs.  
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_alma

  SUBROUTINE hydrol_alma (kjpindex, index, lstep_init, qsintveg, snow, snow_nobio, soilwet)
    !
    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT (in)                        :: kjpindex     !! Domain size
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)   :: index        !! Indeces of the points on the map
    LOGICAL, INTENT (in)                               :: lstep_init   !! At which time is this routine called ?
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: qsintveg     !! Water on vegetation due to interception
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: snow         !! Snow water equivalent
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (in) :: snow_nobio !! Water balance on ice, lakes, .. [Kg/m^2]

    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: soilwet     !! Soil wetness

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std) :: ji
    REAL(r_std) :: watveg

!_ ================================================================================================================================
    !
    !
    IF ( lstep_init ) THEN
       ! Initialize variables if they were not found in the restart file

       DO ji = 1, kjpindex
          watveg = SUM(qsintveg(ji,:))
          tot_watveg_beg(ji) = watveg
          tot_watsoil_beg(ji) = humtot(ji)
          snow_beg(ji)        = snow(ji) + SUM(snow_nobio(ji,:))
       ENDDO

       RETURN

    ENDIF
    !
    ! Calculate the values for the end of the time step
    !
    DO ji = 1, kjpindex
       watveg = SUM(qsintveg(ji,:)) ! average within the mesh
       tot_watveg_end(ji) = watveg
       tot_watsoil_end(ji) = humtot(ji) ! average within the mesh
       snow_end(ji) = snow(ji)+ SUM(snow_nobio(ji,:)) ! average within the mesh

       delintercept(ji) = tot_watveg_end(ji) - tot_watveg_beg(ji) ! average within the mesh
       delsoilmoist(ji) = tot_watsoil_end(ji) - tot_watsoil_beg(ji)
       delswe(ji)       = snow_end(ji) - snow_beg(ji) ! average within the mesh
    ENDDO
    !
    !
    ! Transfer the total water amount at the end of the current timestep top the begining of the next one.
    !
    tot_watveg_beg = tot_watveg_end
    tot_watsoil_beg = tot_watsoil_end
    snow_beg(:) = snow_end(:)
    !
    DO ji = 1,kjpindex
       IF ( mx_eau_var(ji) > 0 ) THEN
          soilwet(ji) = tot_watsoil_end(ji) / mx_eau_var(ji)
       ELSE
          soilwet(ji) = zero
       ENDIF
    ENDDO
    !
  END SUBROUTINE hydrol_alma
  !


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_calculate_temp_hydro
!!
!>\BRIEF         Calculate the temperature at hydrological levels  
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================


  SUBROUTINE hydrol_calculate_temp_hydro(kjpindex, stempdiag, snow,snowdz)

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                             :: kjpindex 
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in)     :: stempdiag
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)          :: snow
    REAL(r_std),DIMENSION (kjpindex,nsnow), INTENT (in)    :: snowdz


    !! 0.2 Local variables
    
    INTEGER jh, jsl, ji
    REAL(r_std) :: snow_h
    REAL(r_std)  :: lev_diag, prev_diag, lev_prog, prev_prog
    REAL(r_std), DIMENSION(nslm,nslm) :: intfactt
    
    
    DO ji=1,kjpindex
       IF (ok_explicitsnow) THEN 
          !The snow pack is above the surface soil in the new snow model.
          snow_h=0
       ELSE  
          snow_h=snow(ji)/sn_dens
       ENDIF
       
       intfactt(:,:)=0.
       prev_diag = snow_h
       DO jh = 1, nslm
          IF (jh.EQ.1) THEN
             lev_diag = zz(2)/1000./2.+snow_h
          ELSEIF (jh.EQ.nslm) THEN
             lev_diag = zz(nslm)/1000.+snow_h
             
          ELSE
             lev_diag = zz(jh)/1000. &
                  & +(zz(jh+1)-zz(jh))/1000./2.+snow_h
             
          ENDIF
          prev_prog = 0.0
          DO jsl = 1, nslm
             lev_prog = diaglev(jsl)
             IF ((lev_diag.GT.diaglev(nslm).AND. &
                  & prev_diag.LT.diaglev(nslm)-min_sechiba)) THEN
                lev_diag=diaglev(nslm)          
             ENDIF
             intfactt(jh,jsl) = MAX(MIN(lev_diag,lev_prog)-MAX(prev_diag, prev_prog),&
                  & 0.0)/(lev_diag-prev_diag)
             prev_prog = lev_prog
          ENDDO
          IF (lev_diag.GT.diaglev(nslm).AND. &
               & prev_diag.GE.diaglev(nslm)-min_sechiba) intfactt(jh,nslm)=1.
          prev_diag = lev_diag
       ENDDO
    ENDDO
    
    temp_hydro(:,:)=0.
    DO jsl= 1, nslm
       DO jh= 1, nslm
          DO ji = 1, kjpindex
             temp_hydro(ji,jh) = temp_hydro(ji,jh) + stempdiag(ji,jsl)*intfactt(jh,jsl)
          ENDDO
       ENDDO
    ENDDO
    
  END SUBROUTINE hydrol_calculate_temp_hydro


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_nudge
!!
!>\BRIEF         Applay nudging of soil moisture and/or snow variables
!!
!! DESCRIPTION  : Nudging of soil moisture and/or snow variables is done if OK_NUDGE_MC=y and/or OK_NUDGE_SNOW=y in run.def
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN IN-OUTPUT VARIABLE(S) : mc, snowdz, snowrho, snowtemp
!!
!! REFERENCE(S) : 
!!
!! \n
!_ ================================================================================================================================

  SUBROUTINE hydrol_nudge(kjit,   kjpindex, &
                          mc_loc, snowdz, snowrho, snowtemp, soiltile)

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                         :: kjit        !! Timestep number
    INTEGER(i_std), INTENT(in)                         :: kjpindex    !! Domain size
    REAL(r_std), DIMENSION(kjpindex,nstm), INTENT (in) :: soiltile    !! Fraction of each soil tile within vegtot (0-1, unitless)

    !! 0.2 Modified variables
    REAL(r_std), DIMENSION(kjpindex,nslm,nstm), INTENT(inout) :: mc_loc      !! Soil moisture
    REAL(r_std), DIMENSION(kjpindex,nsnow), INTENT(inout)     :: snowdz      !! Snow layer thickness
    REAL(r_std), DIMENSION(kjpindex,nsnow), INTENT(inout)     :: snowrho     !! Snow density
    REAL(r_std), DIMENSION(kjpindex,nsnow), INTENT(inout)     :: snowtemp    !! Snow temperature



    !! 0.3 Locals variables
    REAL(r_std)                                :: tau                   !! Position between to values in nudge mc file
    REAL(r_std), DIMENSION(kjpindex,nslm,nstm) :: mc_read_current       !! mc from file interpolated to current timestep
    REAL(r_std), DIMENSION(kjpindex,nsnow)     :: snowdz_read_current   !! snowdz from file interpolated to current timestep
    REAL(r_std), DIMENSION(kjpindex,nsnow)     :: snowrho_read_current  !! snowrho from file interpolated to current timestep
    REAL(r_std), DIMENSION(kjpindex,nsnow)     :: snowtemp_read_current !! snowtemp from file interpolated to current timestep
    REAL(r_std), DIMENSION(kjpindex)           :: nudgincsm             !! Nudging increment of water in soil moisture
    REAL(r_std), DIMENSION(kjpindex)           :: nudgincswe            !! Nudging increment of water in snow
    REAL(r_std), DIMENSION(kjpindex,nslm,nstm) :: mc_aux                !! Temorary variable for calculation of nudgincsm
    REAL(r_std), DIMENSION(kjpindex,nstm)      :: tmc_aux               !! Temorary variable for calculation of nudgincsm
    REAL(r_std), DIMENSION(iim_g,jjm_g,nslm,1) :: mc_read_glo2D_1       !! mc from file at global 2D(lat,lon) grid per soiltile
    REAL(r_std), DIMENSION(iim_g,jjm_g,nslm,1) :: mc_read_glo2D_2       !! mc from file at global 2D(lat,lon) grid per soiltile
    REAL(r_std), DIMENSION(iim_g,jjm_g,nslm,1) :: mc_read_glo2D_3       !! mc from file at global 2D(lat,lon) grid per soiltile
    REAL(r_std), DIMENSION(iim_g,jjm_g,nsnow,1):: snowdz_read_glo2D     !! snowdz from file at global 2D(lat,lon) grid
    REAL(r_std), DIMENSION(iim_g,jjm_g,nsnow,1):: snowrho_read_glo2D    !! snowrho from file at global 2D(lat,lon) grid
    REAL(r_std), DIMENSION(iim_g,jjm_g,nsnow,1):: snowtemp_read_glo2D   !! snowrho from file at global 2D(lat,lon) grid
    REAL(r_std), DIMENSION(nbp_glo,nslm,nstm)  :: mc_read_glo1D         !! mc_read_glo2D on land-only vector form, in global
    REAL(r_std), DIMENSION(nbp_glo,nsnow)      :: snowdz_read_glo1D     !! snowdz_read_glo2D on land-only vector form, in global
    REAL(r_std), DIMENSION(nbp_glo,nsnow)      :: snowrho_read_glo1D    !! snowdz_read_glo2D on land-only vector form, in global
    REAL(r_std), DIMENSION(nbp_glo,nsnow)      :: snowtemp_read_glo1D   !! snowdz_read_glo2D on land-only vector form, in global
    INTEGER(i_std), SAVE                       :: istart_mc, istart_snow!! start index to read from input file
    INTEGER(i_std)                             :: iend                  !! end index to read from input file
    INTEGER(i_std)                             :: i, j, ji, jg, jst, jsl!! loop index
    INTEGER(i_std)                             :: iim_file, jjm_file, llm_file !! Dimensions in input file
    INTEGER(i_std), SAVE                       :: ttm_mc, ttm_snow      !! Time dimensions in input file
    INTEGER(i_std), SAVE                       :: mc_id, snow_id        !! index for netcdf files
    LOGICAL, SAVE                              :: firsttime_mc=.TRUE.
    LOGICAL, SAVE                              :: firsttime_snow=.TRUE.

 
    !! 1. Nudging of soil moisture
    IF (ok_nudge_mc) THEN

       !! 1.2 Read mc from file, once a day only
       !!     The forcing file must contain daily frequency variable for the full year of the simulation
       IF (MOD(kjit,INT(one_day/dt_sechiba)) == 1) THEN 
          ! Save mc read from file from previous day
          mc_read_prev = mc_read_next

          IF (nudge_interpol_with_xios) THEN
             ! Read mc from input file. XIOS interpolates it to the model grid before it is received here.
             CALL xios_orchidee_recv_field("moistc_interp", mc_read_next)

             ! Read and interpolation the mask for variable mc from input file. 
             ! This is only done to be able to output the mask it later for validation purpose.
             ! The mask corresponds to the fraction of the input source file which was underlaying the model grid cell. 
             ! If the msask is 0 for a model grid cell, then the default value 0.2 set in field_def_orchidee.xml, is used for that grid cell. 
             CALL xios_orchidee_recv_field("mask_moistc_interp", mask_mc_interp)

          ELSE

             ! Only read fields from the file. We here suppose that no interpolation is needed.
             IF (is_root_prc) THEN 
                IF (firsttime_mc) THEN
                   ! Open and read dimenions in file
                   CALL flininfo('nudge_moistc.nc',  iim_file, jjm_file, llm_file, ttm_mc, mc_id)
                   
                   ! Coherence test between dimension in the file and in the model run
                   IF ((iim_file /= iim_g) .OR. (jjm_file /= jjm_g)) THEN
                      WRITE(numout,*) 'hydrol_nudge: iim_file, jjm_file, llm_file, ttm_mc=', &
                           iim_file, jjm_file, llm_file, ttm_mc
                      WRITE(numout,*) 'hydrol_nudge: iim_g, jjm_g=', iim_g, jjm_g
                      CALL ipslerr_p(2,'hydrol_nudge','Problem in coherence between dimensions in nudge_moistc.nc file and model',&
                           'iim_file should be equal to iim_g','jjm_file should be equal to jjm_g')
                   END IF
                   
                   firsttime_mc=.FALSE.
                   istart_mc=julian_diff-1 ! initialize time counter to read
                   IF (printlev>=2) WRITE(numout,*) "Start read nudge_moistc.nc file at time step: ", istart_mc+1
                END IF

                istart_mc=istart_mc+1  ! read next time step in the file
                iend=istart_mc         ! only read 1 time step
                
                ! Read mc from file, one variable per soiltile
                IF (printlev>=3) WRITE(numout,*) &
                     "Read variables moistc_1, moistc_2 and moistc_3 from nudge_moistc.nc at time step: ", istart_mc
                CALL flinget (mc_id, 'moistc_1', iim_g, jjm_g, nslm, ttm_mc, istart_mc, iend, mc_read_glo2D_1)
                CALL flinget (mc_id, 'moistc_2', iim_g, jjm_g, nslm, ttm_mc, istart_mc, iend, mc_read_glo2D_2)
                CALL flinget (mc_id, 'moistc_3', iim_g, jjm_g, nslm, ttm_mc, istart_mc, iend, mc_read_glo2D_3)

                ! Transform from global 2D(iim_g, jjm_g) into into land-only global 1D(nbp_glo)
                ! Put the variables on the 3 soiltiles in the same file
                DO ji = 1, nbp_glo
                   j = ((index_g(ji)-1)/iim_g) + 1
                   i = (index_g(ji) - (j-1)*iim_g)
                   mc_read_glo1D(ji,:,1) = mc_read_glo2D_1(i,j,:,1)
                   mc_read_glo1D(ji,:,2) = mc_read_glo2D_2(i,j,:,1)
                   mc_read_glo1D(ji,:,3) = mc_read_glo2D_3(i,j,:,1)
                END DO
             END IF

             ! Distribute the fields on all processors
             CALL scatter(mc_read_glo1D, mc_read_next)

             ! No interpolation is done, set the mask to 1
             mask_mc_interp(:,:,:) = 1

          END IF ! nudge_interpol_with_xios
       END IF ! MOD(kjit,INT(one_day/dt_sechiba)) == 1
       
      
       !! 1.3 Linear time interpolation between daily fields to the current time step
       tau   = (kjit-1)*dt_sechiba/one_day - AINT((kjit-1)*dt_sechiba/one_day)
       mc_read_current(:,:,:) = (1.-tau)*mc_read_prev(:,:,:) + tau*mc_read_next(:,:,:)

       !! 1.4 Output daily fields and time interpolated fields only for debugging and validation purpose
       CALL xios_orchidee_send_field("mc_read_next", mc_read_next)
       CALL xios_orchidee_send_field("mc_read_current", mc_read_current)
       CALL xios_orchidee_send_field("mc_read_prev", mc_read_prev)
       CALL xios_orchidee_send_field("mask_mc_interp_out", mask_mc_interp)

       
       !! 1.5 Applay nudging of soil moisture using alpha_nudge_mc at each model sechiba time step.
       !!     alpha_mc_nudge calculated using the parameter for relaxation time NUDGE_TAU_MC set in module constantes.
       !!     alpha_nudge_mc is between 0-1
       !!     If alpha_nudge_mc=1, the new mc will be replaced by the one read from file 
       mc_loc(:,:,:) = (1-alpha_nudge_mc)*mc_loc(:,:,:) + alpha_nudge_mc * mc_read_current(:,:,:)
    

       !! 1.6 Calculate diagnostic for nudging increment of water in soil moisture
       mc_aux(:,:,:)  = alpha_nudge_mc * ( mc_read_current(:,:,:) - mc_loc(:,:,:))
       DO jst=1,nstm
          DO ji=1,kjpindex
             tmc_aux(ji,jst) = dz(2) * ( trois*mc_aux(ji,1,jst) + mc_aux(ji,2,jst) )/huit
             DO jsl = 2,nslm-1
                tmc_aux(ji,jst) = tmc_aux(ji,jst) + dz(jsl) *  (trois*mc_aux(ji,jsl,jst)+mc_aux(ji,jsl-1,jst))/huit &
                     + dz(jsl+1) * (trois*mc_aux(ji,jsl,jst)+mc_aux(ji,jsl+1,jst))/huit
             ENDDO
             tmc_aux(ji,jst) = tmc_aux(ji,jst) + dz(nslm) * (trois*mc_aux(ji,nslm,jst) + mc_aux(ji,nslm-1,jst))/huit
          ENDDO
       ENDDO
       
       ! Average over grid-cell
       nudgincsm(:) = zero
       DO jst=1,nstm
          DO ji=1,kjpindex
             nudgincsm(ji) = nudgincsm(ji) + vegtot(ji) * soiltile(ji,jst) * tmc_aux(ji,jst)
          ENDDO
       ENDDO
       
       CALL xios_orchidee_send_field("nudgincsm", nudgincsm)
       
       
    END IF ! IF (ok_nudge_mc)


    !! 2. Nudging of snow variables
    IF (ok_nudge_snow) THEN

       !! 2.1 Read snow variables from file, once a day only
       !!     The forcing file must contain daily frequency values for the full year of the simulation
       IF (MOD(kjit,INT(one_day/dt_sechiba)) == 1) THEN 
          ! Save variables from previous day
          snowdz_read_prev   = snowdz_read_next
          snowrho_read_prev  = snowrho_read_next
          snowtemp_read_prev = snowtemp_read_next
          
          IF (nudge_interpol_with_xios) THEN
             ! Read and interpolation snow variables and the mask from input file
             CALL xios_orchidee_recv_field("snowdz_interp", snowdz_read_next)
             CALL xios_orchidee_recv_field("snowrho_interp", snowrho_read_next)
             CALL xios_orchidee_recv_field("snowtemp_interp", snowtemp_read_next)
             CALL xios_orchidee_recv_field("mask_snow_interp", mask_snow_interp)

          ELSE
             ! Only read fields from the file. We here suppose that no interpolation is needed.
             IF (is_root_prc) THEN 
                IF (firsttime_snow) THEN
                   ! Open and read dimenions in file
                   CALL flininfo('nudge_snow.nc',  iim_file, jjm_file, llm_file, ttm_snow, snow_id)
                   
                   ! Coherence test between dimension in the file and in the model run
                   IF ((iim_file /= iim_g) .OR. (jjm_file /= jjm_g)) THEN
                      WRITE(numout,*) 'hydrol_nudge: iim_file, jjm_file, llm_file, ttm_snow=', &
                           iim_file, jjm_file, llm_file, ttm_snow
                      WRITE(numout,*) 'hydrol_nudge: iim_g, jjm_g=', iim_g, jjm_g
                      CALL ipslerr_p(3,'hydrol_nudge','Problem in coherence between dimensions in nudge_snow.nc file and model',&
                           'iim_file should be equal to iim_g','jjm_file should be equal to jjm_g')
                   END IF
                                         
                   firsttime_snow=.FALSE.
                   istart_snow=julian_diff-1  ! initialize time counter to read
                   IF (printlev>=2) WRITE(numout,*) "Start read nudge_snow.nc file at time step: ", istart_snow+1
                END IF

                istart_snow=istart_snow+1  ! read next time step in the file
                iend=istart_snow      ! only read 1 time step
                
                ! Read snowdz, snowrho and snowtemp from file
                IF (printlev>=3) WRITE(numout,*) &
                     "Read variables snowdz, snowrho and snowtemp from nudge_snow.nc at time step: ", istart_snow
                CALL flinget (snow_id, 'snowdz', iim_g, jjm_g, nsnow, ttm_snow, istart_snow, iend, snowdz_read_glo2D)
                CALL flinget (snow_id, 'snowrho', iim_g, jjm_g, nsnow, ttm_snow, istart_snow, iend, snowrho_read_glo2D)
                CALL flinget (snow_id, 'snowtemp', iim_g, jjm_g, nsnow, ttm_snow, istart_snow, iend, snowtemp_read_glo2D)


                ! Transform from global 2D(iim_g, jjm_g) variables into into land-only global 1D variables (nbp_glo)
                DO ji = 1, nbp_glo
                   j = ((index_g(ji)-1)/iim_g) + 1
                   i = (index_g(ji) - (j-1)*iim_g)
                   snowdz_read_glo1D(ji,:) = snowdz_read_glo2D(i,j,:,1)
                   snowrho_read_glo1D(ji,:) = snowrho_read_glo2D(i,j,:,1)
                   snowtemp_read_glo1D(ji,:) = snowtemp_read_glo2D(i,j,:,1)
                END DO
             END IF

             ! Distribute the fields on all processors
             CALL scatter(snowdz_read_glo1D, snowdz_read_next)
             CALL scatter(snowrho_read_glo1D, snowrho_read_next)
             CALL scatter(snowtemp_read_glo1D, snowtemp_read_next)

             ! No interpolation is done, set the mask to 1
             mask_snow_interp=1

          END IF ! nudge_interpol_with_xios
       END IF ! MOD(kjit,INT(one_day/dt_sechiba)) == 1
       
      
       !! 2.2 Linear time interpolation between daily fields for current time step
       tau   = (kjit-1)*dt_sechiba/one_day - AINT((kjit-1)*dt_sechiba/one_day)
       snowdz_read_current(:,:) = (1.-tau)*snowdz_read_prev(:,:) + tau*snowdz_read_next(:,:)
       snowrho_read_current(:,:) = (1.-tau)*snowrho_read_prev(:,:) + tau*snowrho_read_next(:,:)
       snowtemp_read_current(:,:) = (1.-tau)*snowtemp_read_prev(:,:) + tau*snowtemp_read_next(:,:)

       !! 2.3 Output daily fields and time interpolated fields only for debugging and validation purpose
       CALL xios_orchidee_send_field("snowdz_read_next", snowdz_read_next)
       CALL xios_orchidee_send_field("snowdz_read_current", snowdz_read_current)
       CALL xios_orchidee_send_field("snowdz_read_prev", snowdz_read_prev)
       CALL xios_orchidee_send_field("snowrho_read_next", snowrho_read_next)
       CALL xios_orchidee_send_field("snowrho_read_current", snowrho_read_current)
       CALL xios_orchidee_send_field("snowrho_read_prev", snowrho_read_prev)
       CALL xios_orchidee_send_field("snowtemp_read_next", snowtemp_read_next)
       CALL xios_orchidee_send_field("snowtemp_read_current", snowtemp_read_current)
       CALL xios_orchidee_send_field("snowtemp_read_prev", snowtemp_read_prev)
       CALL xios_orchidee_send_field("mask_snow_interp_out", mask_snow_interp)

       !! 2.4 Applay nudging of snow variables using alpha_nudge_snow at each model sechiba time step.
       !!     alpha_snow_nudge calculated using the parameter for relaxation time NUDGE_TAU_SNOW set in module constantes.
       !!     alpha_nudge_snow is between 0-1
       !!     If alpha_nudge_snow=1, the new snow variables will be replaced by the ones read from file.
       snowdz(:,:) = (1-alpha_nudge_snow)*snowdz(:,:) + alpha_nudge_snow * snowdz_read_current(:,:)
       snowrho(:,:) = (1-alpha_nudge_snow)*snowrho(:,:) + alpha_nudge_snow * snowrho_read_current(:,:)
       snowtemp(:,:) = (1-alpha_nudge_snow)*snowtemp(:,:) + alpha_nudge_snow * snowtemp_read_current(:,:)

       !! 2.5 Calculate diagnostic for the nudging increment of water in snow
       nudgincswe=0.
       DO jg = 1, nsnow 
          nudgincswe(:) = nudgincswe(:) +  &
               alpha_nudge_snow*(snowdz_read_current(:,jg)*snowrho_read_current(:,jg)-snowdz(:,jg)*snowrho(:,jg))
       END DO
       CALL xios_orchidee_send_field("nudgincswe", nudgincswe)
       
    END IF


  END SUBROUTINE hydrol_nudge

!!
!================================================================================================================================
!! SUBROUTINE   : read_refSOC_1dfile
!!
!>\BRIEF          
!!
!! DESCRIPTION  : Read file of soil organic carbon to be used in thermix
!! (insulating effect)
!!                
!!
!! RECENT CHANGE(S) : None
!! 
!! MAIN OUTPUT VARIABLE(S): refSOC : soil organic carbon from data
!!                          
!! REFERENCE(S) :
!!
!! FLOWCHART    : None 
!! \n 
!_
!================================================================================================================================
  SUBROUTINE read_refSOCfile_1dfile(kjpindex, filename_in, varname_in, lalo, data_out)
    
    USE interpweight

    IMPLICIT NONE

    !! 0. Variables and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                        :: kjpindex
    CHARACTER(LEN=*), INTENT(in)            :: filename_in, varname_in
    REAL(r_std), DIMENSION(:,:), INTENT(in)    :: lalo

    !! 0.2 Output variables
    REAL(r_std), DIMENSION(:,:), INTENT(out) :: data_out 

    !! 0.3 Local variables
    INTEGER(i_std) :: ib
    CHARACTER(LEN=80) :: filename
!    REAL(r_std),DIMENSION(kjpindex,1) :: data_out_file                          !! Horizontal temperature field interpolated from file [C]
    INTEGER(i_std),DIMENSION(kjpindex,8) :: neighbours
    REAL(r_std)                                          :: vmin, vmax       !! min/max values to use for the
                                                                             !!   renormalization
    REAL(r_std), DIMENSION(kjpindex)                     :: adata_out         !! Availability of data for  the interpolation
    CHARACTER(LEN=80)                                    :: lonname, latname !! lon, lat names in input file
    REAL(r_std), DIMENSION(:), ALLOCATABLE               :: variabletypevals !! Values for all the types of the variable
                                                                             !!   (variabletypevals(1) = -un, not used)
    CHARACTER(LEN=50)                                    :: fractype         !! method of calculation of fraction
                                                                             !!   'XYKindTime': Input values are kinds 
                                                                             !!     of something with a temporal 
                                                                             !!     evolution on the dx*dy matrix'
    LOGICAL                                              :: nonegative       !! whether negative values should be removed
    CHARACTER(LEN=50)                                    :: maskingtype      !! Type of masking
                                                                             !!   'nomask': no-mask is applied
                                                                             !!   'mbelow': take values below maskvals(1)
                                                                             !!   'mabove': take values above maskvals(1)
                                                                             !!   'msumrange': take values within 2 ranges;
                                                                             !!      maskvals(2) <= SUM(vals(k)) <= maskvals(1)
                                                                             !!      maskvals(1) < SUM(vals(k)) <= maskvals(3)
                                                                             !!       (normalized by maskvals(3))
                                                                             !!   'var': mask values are taken from a 
                                                                             !!     variable inside the file (>0)
    REAL(r_std), DIMENSION(3)                            :: maskvals         !! values to use to mask (according to 
                                                                             !!   `maskingtype') 
    CHARACTER(LEN=250)                                   :: namemaskvar      !! name of the variable to use to mask 
    REAL(r_std)                                          :: data_out_norefinf
    REAL(r_std)             :: data_out_default  !! Default value
   
    INTEGER(i_std) :: ier


    IF (printlev >= 1) WRITE(numout,*) "read_refSOCfile_1dfile start" &
         // TRIM(filename_in) //" for variable " //TRIM(varname_in)


    ! For this case there are not types/categories. We have 'only' a continuos field
    ! Assigning values to vmin, vmax

    vmin = 0. ! unused for this type of interpolation
    vmax = 9999.

!   For this file we do not need neightbours!
    neighbours = 0

    !! Variables for interpweight
    ! Type of calculation of cell fractions
    fractype = 'default'
    ! Name of the longitude and latitude in the input file
    lonname = 'longitude'
    latname = 'latitude'
    ! Default value when no value is get from input file
    data_out_default = 0. 
    ! Reference value when no value is get from input file
    data_out_norefinf = 1.
    ! Should negative values be set to zero from input file?
    nonegative = .FALSE.
    ! Type of mask to apply to the input data (see header for more details)
    maskingtype = 'var' !- 'mabove'
    ! Values to use for the masking (here not used)
    maskvals = (/ undef_sechiba, undef_sechiba, undef_sechiba /)
    ! Name of the variable with the values for the mask in the input file (only if maskkingtype='var') (here not used)
    namemaskvar = 'mask'

    CALL interpweight_2Dcont_tstep(kjpindex, lalo, resolution, neighbours,                            &
      contfrac, filename_in, varname_in, lonname, latname, vmin, vmax, nonegative, maskingtype,        &
      maskvals, namemaskvar, -1, fractype, data_out_default, data_out_norefinf,                         &
      data_out, adata_out)

    IF (printlev >= 5) WRITE(numout,*)'read_refSOCfile_1dfile end'

  END SUBROUTINE read_refSOCfile_1dfile

!!
!================================================================================================================================
!! SUBROUTINE   : hydrol_update_nroot
!!
!>\BRIEF          Calculate nroot profile
!!
!! DESCRIPTION  : Dynamic nroot to optimize water use: the root profile used to weight the water stress function 
!!              of each soil layer is updated at each time step in order to match the soil water profile 
!!              (the soil water content of each layer available for transpiration)
!!
!! RECENT CHANGE(S) : None
!! 
!! MAIN OUTPUT VARIABLE(S): nroot 
!!                          
!! REFERENCE(S) :
!!
!! FLOWCHART    : None 
!! \n 
!_
!================================================================================================================================
  SUBROUTINE hydrol_update_nroot(kjpindex, nslm, nvm, sm, smw, smf, sms, smf_tmp, smw_tmp, pcent, nroot_var)

    ! For the calculation of soil_wet_ns and us/humrel/vegstress 
    INTEGER(i_std)                                     :: kjpindex, nslm, nvm !! Array Dimensions
    REAL(r_std), DIMENSION (kjpindex,nslm), INTENT(in) :: sm                !! Soil moisture of each layer (liquid phase)
                                                                            !!  @tex $(kg m^{-2})$ @endtex
    REAL(r_std), DIMENSION (kjpindex,nslm), INTENT(in) :: smw               !! Soil moisture of each layer at wilting point
                                                                            !!  @tex $(kg m^{-2})$ @endtex 
    REAL(r_std), DIMENSION (kjpindex,nslm), INTENT(in) :: smf               !! Soil moisture of each layer at field capacity
                                                                            !!  @tex $(kg m^{-2})$ @endtex   
    REAL(r_std), DIMENSION (kjpindex,nslm), INTENT(in) :: sms               !! Soil moisture of each layer at saturation
                                                                            !!  @tex $(kg m^{-2})$ @endtex 
    REAL(r_std), DIMENSION (kjpindex,nslm), INTENT(in) :: smf_tmp
    REAL(r_std), DIMENSION (kjpindex,nslm), INTENT(in) :: smw_tmp
    REAL(r_std), DIMENSION (kjpindex),   INTENT(in)    :: pcent             !! Fraction of saturated volumetric soil moisture above 
                                                                            !! which transpir is max (0-1, unitless)

    REAL(r_std), DIMENSION (kjpindex,nvm,nslm), INTENT(inout)  :: nroot_var !! Normalized root length fraction in each soil layer 
                                                                            !! (0-1, unitless)
                                                                            !! DIM = kjpindex * nvm * nslm
    ! Local Variables
    REAL(r_std), DIMENSION(kjpindex,nvm,nslm) :: nroot_tmp_fast             !! Temporary variable to calculate the nroot
    REAL(r_std), DIMENSION(kjpindex,nvm)      :: sum_nroot_tmp, nroot_all   !! Temp var to calculate sum of nroot tmp 
    REAL(r_std), DIMENSION(kjpindex,nslm)     :: first, second, third, fourth  !! Temporary variable for partial calculations

    INTEGER(i_std) :: ji, jv, jsl, ier

  !===========================================================================

    nroot_tmp_fast(:,:,:) = zero
    sum_nroot_tmp(:,:) = zero
    nroot_all = 0

    ! precalculate nroot 
    DO jsl = 2,nslm
       !DIR$ VECTOR ALIGNED
       DO ji = 1,kjpindex
          first(ji,jsl) = (sm(ji,jsl)-smw_tmp(ji,jsl))
          second(ji,jsl) = (pcent(ji)*(smf_tmp(ji,jsl)-smw_tmp(ji,jsl)))
          third(ji,jsl) = (smf(ji,jsl)-smw(ji,jsl))
          fourth(ji,jsl) = (smf_tmp(ji,jsl)-smw_tmp(ji,jsl))
          nroot_all(ji,jsl) = MIN(un,MAX(zero,first(ji,jsl) / second(ji,jsl) * third(ji,jsl) / fourth(ji,jsl)))
       ENDDO ! ji
    ENDDO ! jsl

    ! calculate nroot
    DO jsl = 2,nslm
       DO jv = 1, nvm
          !DIR$ VECTOR ALIGNED
          DO ji = 1,kjpindex
             IF (nroot_var(ji,jv,jsl) .NE. zero) THEN
                IF ( is_tree(jv) ) THEN
                   nroot_tmp_fast(ji,jv,jsl) = nroot_all(ji,jsl) 
                   sum_nroot_tmp(ji,jv) = sum_nroot_tmp(ji,jv) + nroot_tmp_fast(ji,jv,jsl)
                        
                ELSE IF ( .NOT. is_tree(jv) .AND. znt(jsl) .LT. un) THEN
                   nroot_tmp_fast(ji,jv,jsl) = nroot_all(ji,jsl) 
                   sum_nroot_tmp(ji,jv) = sum_nroot_tmp(ji,jv) + nroot_tmp_fast(ji,jv,jsl)
                ENDIF ! 
             ENDIF ! nroot_var(ji,jv,jsl) .NE. zero
          ENDDO ! ji
       ENDDO ! jv
    ENDDO ! jsl

    ! arithmetic mean
    DO jsl = 2,nslm
       DO jv = 1, nvm
          !DIR$ VECTOR ALIGNED
          DO ji = 1,kjpindex
             IF (sum_nroot_tmp(ji,jv) .GT. zero) THEN
                IF (.NOT. is_tree(jv) .AND. znt(jsl) .GE. un ) THEN
                   nroot_var(ji,jv,jsl)=zero
                ELSE
                   nroot_var(ji,jv,jsl)=nroot_tmp_fast(ji,jv,jsl)/sum_nroot_tmp(ji,jv)
                ENDIF
             ENDIF
          ENDDO ! ji
       ENDDO ! jv
    ENDDO ! jsl

  END SUBROUTINE hydrol_update_nroot

  
END MODULE hydrol
