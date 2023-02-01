! =================================================================================================================================
! MODULE 	: stomate_vmax
!
! CONTACT	: orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      	: IPSL (2006). This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        calculates the leaf efficiency.
!!	
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! SVN		:
!! $HeadURL$ 
!! $Date$
!! $Revision$
!! \n
!_ =================================================================================================================================

MODULE stomate_vmax

  ! modules used:

  USE ioipsl_para
  USE stomate_data
  USE constantes
  USE pft_parameters

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC vmax, vmax_clear

  ! first call
  LOGICAL, SAVE                                              :: firstcall_vmax = .TRUE.
!$OMP THREADPRIVATE(firstcall_vmax)
!gmjc
  LOGICAL, SAVE                                           :: ok_Nlim = .FALSE.
!end gmjc
CONTAINS

!! ================================================================================================================================
!! SUBROUTINE	: vmax_clear
!!
!>\BRIEF	  Flag setting 
!!
!!\n DESCRIPTION: This subroutine sets flags ::firstcall_vmax, to .TRUE., and therefore activates   
!!		  section 1.1 of the ::vmax subroutine which writes messages to the output. \n
!!		  This subroutine is called at the end of the subroutine ::stomate_clear, in the 
!!		  module ::stomate.
!!
!! RECENT CHANGE(S):None
!!
!! MAIN OUTPUT VARIABLE(S): ::firstcall_vmax
!!
!! REFERENCE(S)  : None 
!!
!! FLOWCHART     : None
!! \n		  
!_ =================================================================================================================================

  SUBROUTINE vmax_clear
    firstcall_vmax=.TRUE.
  END SUBROUTINE vmax_clear



!! ================================================================================================================================
!! SUBROUTINE    : vmax
!!
!>\BRIEF         This subroutine computes vcmax photosynthesis parameters 
!! given optimal vcmax parameter values and a leaf age-related efficiency.
!!
!! DESCRIPTION (functional, design, flags): 
!! Leaf age classes are introduced to take into account the fact that photosynthetic activity depends on leaf age
!! (Ishida et al., 1999). There are \f$nleafages\f$ classes (constant defined in stomate_constants.f90).
!! This subroutine first calculates the new age of each leaf age-class based on fraction of leaf 
!! that goes from one to another class.                                              
!! Then calculation of the new fraction of leaf in each class is performed.      
!! Last, leaf efficiency is calculated for each PFT and for each leaf age class.
!! vcmax is defined as vcmax25 and vjmax_opt weighted by a mean leaf
!! efficiency. vcmax25 is PFT-dependent constants defined in constants_mtc.f90.
!!
!! This routine is called once at the beginning by stomate_var_init and then at each stomate time step by stomateLpj.
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): vcmax
!!
!! REFERENCE(S)	: 
!! - Ishida, A., A. Uemura, N. Koike, Y. Matsumoto, and A. Lai Hoe (1999),
!! Interactive effects of leaf age and self-shading on leaf structure, photosynthetic
!! capacity and chlorophyll fluorescence in the rain forest tree,
!! dryobalanops aromatica, Tree Physiol., 19, 741-747
!!
!! FLOWCHART    : None
!!
!! REVISION(S)	: None
!! \n
!_ ================================================================================================================================

  SUBROUTINE vmax (npts, dt, &
       leaf_age, leaf_frac, &
       phytomer_age, bm_phytomer, bm_FFB, PHYbm, FFBbm, biomass, & !! yidi
       vcmax, vcmax_cl1,vcmax_cl2,vcmax_cl3,vcmax_cl4, VPD_week, ave_VPD, & !!xchen
!gmjc
       N_limfert)
!end gmjc
    !
    !! 0. Variable and parameter declaration
    !

    !
    !! 0.1 Input variables
    !
    INTEGER(i_std), INTENT(in)                                 :: npts                    !! Domain size (unitless)
    REAL(r_std), INTENT(in)                                    :: dt                      !! time step of stomate (days)

    !
    !! 0.2 Output variables 
    !
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)              :: vcmax                   !! Maximum rate of carboxylation 
                                                                                          !! @tex ($\mu mol m^{-2} s^{-1}$) @endtex

    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)              :: vcmax_cl1               !! @xchen  Maximum rate of carboxylation
                                                                                          !! for leaf age 1
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)              :: vcmax_cl2               !! leaf age 2 
                                                                                          !! 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)              :: vcmax_cl3               !! leaf age 3 
                                                                                          !! 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)              :: vcmax_cl4               !! leaf age 4
                                                                                          !! unit ($\mu mol m^{-2} s^{-1}$) 

    REAL(r_std), DIMENSION(npts), INTENT(in)         :: VPD_week              !! xchen weekly VPD
    REAL(r_std), DIMENSION(npts), INTENT(inout)       :: ave_VPD              !! Vapor Pressure Deficit (kPa)
    REAL(r_std), DIMENSION(npts)       :: ave_VPD_old              !! Vapor Pressure Deficit (kPa)

    !
    !! 0.3 Modified variables
    !
    REAL(r_std), DIMENSION(npts,nvm,nleafages), INTENT(inout)  :: leaf_age                !! Leaf age (days)
    REAL(r_std), DIMENSION(npts,nvm,nleafages), INTENT(inout)  :: leaf_frac               !! fraction of leaves in leaf age 
                                                                                          !! classes 
                                                                                          !! (unitless)
!gmjc

!! yidi
    REAL(r_std), DIMENSION(npts,nvm,nphs), INTENT(inout)       :: phytomer_age            !! age of the phytomers (days)
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: biomass           !! biomass @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm,nphs),INTENT(inout)        :: bm_phytomer             !! Each PHYTOMER mass, from sapabove 
                                                                                          !! @tex $(gC.m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nphs),INTENT(inout)        :: bm_FFB                  !! FRUIT mass for each PHYTOMER, from sapabove 
                                                                                          !! @tex $(gC.m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm),INTENT(inout)             :: PHYbm                   !! PHYTOMER mass, from sapabove
                                                                                          !! @tex $(gC.m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm),INTENT(inout)             :: FFBbm                   !! Fruit mass for all phytomers
                                                                                          !! @tex $(gC.m^{-2})$ @endtex
   ! REAL(r_std), DIMENSION(npts),INTENT(inout)                   :: FFBharvest              !! harvest fruit biomass
                                                                                          !! @tex $(gC.m^{-2})$ @endtex
!! yidi
    ! N fertilization limitation factor for grassland Vcmax and SLA
    REAL(r_std), DIMENSION(:,:), INTENT(in)                    :: N_limfert
!end gmjc
    !
    !! 0.4 Local variables
    !
    REAL(r_std), DIMENSION(npts)                               :: leaf_efficiency         !! leaf efficiency (vcmax/vcmax25)
                                                                                          !! (unitless)
    REAL(r_std), DIMENSION(npts,nvm,nleafages)                 :: d_leaf_frac             !! turnover between age classes
                                                                                          !! (unitless)
    REAL(r_std), DIMENSION(npts,nleafages)                     :: leaf_age_new            !! new leaf age (days)
    REAL(r_std), DIMENSION(npts)                               :: sumfrac                 !! sum of leaf age fractions, 
                                                                                          !! for normalization
                                                                                          !! (unitless)
    REAL(r_std), DIMENSION(npts)                               :: rel_age                 !! relative leaf age (age/critical age)
                                                                                          !! (unitless)
    INTEGER(i_std)                                             :: i,j,m,p                   !! indices (unitless)

!! yidi
    REAL(r_std), DIMENSION(npts,nvm)                           :: ffbharvest_age_crit     !! critical ffbharvest age (days)
    REAL(r_std), DIMENSION(npts,nvm)                           :: phytomer_age_crit       !! critical phytomer age (days)
    REAL(r_std), DIMENSION(npts,nvm)                           :: phytomer_timecst_crit   !! critical phytomer update (days)
!! yidi
!_ ================================================================================================================================

    IF (printlev>=3) WRITE(numout,*) 'Entering vmax'

    !
    !! 1 Initialization
    !

    !
    !! 1.1 first call: info about flags and parameters.
    !

    IF ( firstcall_vmax ) THEN
!gmjc
      ok_Nlim=.false.
      CALL GETIN ('GRM_N_LIMITATION',ok_Nlim)
      IF (printlev >= 2) THEN
          WRITE(numout,*) 'GRM_N_LIMITATION',ok_Nlim
!end gmjc
          WRITE(numout,*) '   > offset (minimum vcmax/vmax_opt):' , vmax_offset
          WRITE(numout,*) '   > relative leaf age at which vmax reaches vcmax_opt:', leafage_firstmax 
          WRITE(numout,*) '   > relative leaf age at which vmax falls below vcmax_opt:', leafage_lastmax
          WRITE(numout,*) '   > relative leaf age at which vmax reaches its minimum:', leafage_old
      END IF
      firstcall_vmax = .FALSE.

    ENDIF

    !
    !! 1.2 initialize output
    !

    vcmax(:,:) = zero
!!xchen
    ave_VPD_old = MAX(2.0, VPD_week)                                           
    ave_VPD_old = LOG(ave_VPD_old)*LOG(ave_VPD_old)*LOG(ave_VPD_old) 

    !
    !! 2 leaf age: general increase and turnover between age classes.
    !

    !
    !! 2.1 increase leaf age
    !
!
!! The age of the leaves in each leaf-age-class increases by 1 time step.
    DO m = 1, nleafages	! Loop over # leaf age classes 
       DO j = 2,nvm	! Loop over # PFTs
          WHERE ( leaf_frac(:,j,m) .GT. min_stomate )

             leaf_age(:,j,m) = leaf_age(:,j,m) + dt
             
          ENDWHERE
       ENDDO	! Loop over # PFTs 

    ENDDO	! Loop over # leaf age classes

!! yidi 
     !! 1. transfer between phytomer classes 2.update the phytomer_Age
     IF (ok_oilpalm) THEN
        DO j = 2,nvm
           IF (is_oilpalm(j)) THEN
             ! DO p = 2, nphs
             !    IF ( phytomer_age(:,j,p) .GE. ffbharvestagecrit(j) ) THEN
             !       FFBharvest(:) = FFBharvest(:) +  bm_FFB(:,j,p)
             !       FFBbm(:,j) = FFBbm(:,j) - bm_FFB(:,j,p)
             !       bm_FFB(:,j,p) = zero
             !       biomass(:,j,isapabove,nelements) = biomass(:,j,isapabove,nelements) - bm_FFB(:,j,p)
             !    ENDIF
             ! ENDDO
              phytomer_timecst_crit(:,j)=phytomer_timecst(j)
              WRITE(numout,*) 'yd: vmax1,nvm=',j,'phytomer_timecst(j)',phytomer_timecst(j)
   !              WHERE ( phytomer_age(:,j,1) .GE. phytomer_timecst_crit(:,j) )  !phytomeragecrit(j)/nphs ) THEN
   !                 PHYbm(:,j) = PHYbm(:,j) - bm_phytomer(:,j,nphs)
   !                 FFBbm(:,j) = FFBbm(:,j) - bm_FFB(:,j,nphs)
   !              ENDWHERE
   !              DO p = nphs,2,-1
   !                 WHERE ( phytomer_age(:,j,1) .GE. phytomer_timecst_crit(:,j))
   !                    phytomer_age(:,j,p) = phytomer_age(:,j,p-1)
   !                    bm_phytomer(:,j,p) = bm_phytomer(:,j,p-1)
   !                    bm_FFB(:,j,p) = bm_FFB(:,j,p-1)
   !                 ENDWHERE
   !              !   phytomer_leaffrac(:,j,p,:) = phytomer_leaffrac(:,j,p-1,:)
   !              ENDDO
   !              WHERE ( phytomer_age(:,j,1) .GE. phytomer_timecst_crit(:,j) )  !phytomeragecrit(j)/nphs ) THEN
   !                 phytomer_age(:,j,1) = zero
   !                 bm_phytomer(:,j,1) = zero
   !                 bm_FFB(:,j,1) = zero
   !              ! phytomer_leaffrac(:,j,1,:) = zero
   !              ENDWHERE
              DO i = 1, npts
                 WRITE(numout,*) 'yd: vmax1,before turn ',j,'nvm,','phytomer_age(i,j,:)',phytomer_age(:,j,:)
                 IF ( phytomer_age(i,j,1) .GE. phytomer_timecst(j) ) THEN  !phytomeragecrit(j)/nphs ) THEN
                    WRITE(numout,*) 'yd: vmax1,before update nph nvm=',j,'nphs=40,bm_phytomer',bm_phytomer(:,j,nphs)
                    WRITE(numout,*) 'yd: vmax1,before update nph nvm=',j,'PHYbm',PHYbm(:,j)
                    WRITE(numout,*) 'yd: vmax1,before update nph nvm=',j,'nphs=40,bm_FFB',bm_FFB(:,j,nphs)
                    WRITE(numout,*) 'yd: vmax1,before update nph nvm=',j,'FFBbm',FFBbm(:,j)
                    WRITE(numout,*) 'yd: vmax1,before update nph nvm=',j,'biomass(isapabove)',biomass(:,j,isapabove,icarbon)
                    PHYbm(i,j) = PHYbm(i,j) - bm_phytomer(i,j,nphs)
                    FFBbm(i,j) = FFBbm(i,j) - bm_FFB(i,j,nphs)
                    biomass(i,j,isapabove,icarbon)=biomass(i,j,isapabove,icarbon)-bm_phytomer(i,j,nphs)-bm_FFB(i,j,nphs)
                    WRITE(numout,*) 'yd: vmax2,after update nph nvm=',j,'nphs=40,bm_phytomer',bm_phytomer(:,j,nphs)
                    WRITE(numout,*) 'yd: vmax2,after update nph nvm=',j,'PHYbm',PHYbm(:,j)
                    WRITE(numout,*) 'yd: vmax2,after update nph nvm=',j,'nphs=40,bm_FFB',bm_FFB(:,j,nphs)
                    WRITE(numout,*) 'yd: vmax2,after update nph nvm=',j,'FFBbm',FFBbm(:,j)
                    WRITE(numout,*) 'yd: vmax2,after update nph nvm=',j,'biomass(isapabove)',biomass(:,j,isapabove,icarbon)
                    DO p = nphs,2,-1
                       phytomer_age(i,j,p) = phytomer_age(i,j,p-1)
                       bm_phytomer(i,j,p) = bm_phytomer(i,j,p-1)
                       bm_FFB(i,j,p) = bm_FFB(i,j,p-1)
                    ENDDO
                    phytomer_age(i,j,1) = zero
                    bm_phytomer(i,j,1) = zero
                    bm_FFB(i,j,1) = zero
                    WRITE(numout,*) 'yd: vmax2,after update nph nvm=',j,'bm_FFB',bm_FFB(:,j,:)
                    WRITE(numout,*) 'yd: vmax2,after update nph nvm=',j,'bm_phytomer',bm_phytomer(:,j,:)
                    WRITE(numout,*) 'yd: vmax2,after update nph nvm=',j,'phytomer_age',phytomer_age(:,j,:)
                 ENDIF
              ENDDO

              phytomer_age(:,j,1) = phytomer_age(:,j,1) + dt
              DO p = 2,nphs
                 WHERE ( phytomer_age(:,j,p) .GT. min_stomate) 
                         ! phytomer_frac(:,j,p) .GT. min_stomate .AND. &
                         ! phytomer_age(:,j,1) .GE. phytomeragecrit(j)/nphs &
                         ! phytomer_age(:,j,p) .GT. min_stomate)
                    phytomer_age(:,j,p) = phytomer_age(:,j,p) + dt
                 ENDWHERE
              ENDDO ! Loop over # phytomers
              WRITE(numout,*) 'yd: vmax3,after age update nvm=',j,'phytomer_age',phytomer_age(:,j,:)
           ENDIF ! is_oilpalm?
        ENDDO ! Loop over # PFTs
     ENDIF ! ok_oilpalm
!! yidi
    !
    !! 2.2 turnover between leaf age classes
    !     d_leaf_frac(:,:,m) = what leaves m-1 and goes into m
    !


    DO j = 2,nvm	! Loop over # PFTs

       !! 2.2.1 fluxes

       !! nothing goes into first age class
       d_leaf_frac(:,j,1) = zero

       !! for others age classes (what goes from m-1 to m)
       DO m = 2, nleafages 
!! leaf_timecst is defined in stomate_constants.f90 as the quotient of the critical leaf age per the number of age classes.
!! The critical leaf age is a PFT-dependent constant defined in stomate_constants.f90, that represents the leaf life span.
!! This time constant (leaf_timecst) determines the turnover between the nleafages different leaf age classes
!! (see section [118] in Krinner et al. (2005)).
          d_leaf_frac(:,j,m) = leaf_frac(:,j,m-1) * dt/leaf_timecst(j)
          !!xchen
            !!WHERE ((VPD_week(:) .LT. 6.0/1.7) .AND. (leaf_frac(:,j,m-1) .GT. 20*d_leaf_frac(:,j,m)))                                   
            !!    d_leaf_frac(:,j,m) = 20*d_leaf_frac(:,j,m)
            !!ENDWHERE 
            !!WHERE (ave_VPD_old(:)*dt/leaf_timecst(j) .GT. 5.0)
            !!    d_leaf_frac(:,j,m) = 5.0/ave_VPD_old(:)*d_leaf_frac(:,j,m)         
            !!ENDWHERE
 
       ENDDO

       !! 2.2.2 new leaf age in class
       !!       new age = ( old age * (old fraction - fractional loss) + fractional increase * age of the source class ) / new fraction
       !!       The leaf age of the youngest class (m=1) is updated into stomate_alloc          
       leaf_age_new(:,:) = zero

       DO m = 2, nleafages-1       ! Loop over age classes
	!! For all age classes except first and last 
          WHERE ( d_leaf_frac(:,j,m) .GT. min_stomate )

             leaf_age_new(:,m) = ( ( (leaf_frac(:,j,m)- d_leaf_frac(:,j,m+1)) * leaf_age(:,j,m) )  + &
                  ( d_leaf_frac(:,j,m) * leaf_age(:,j,m-1) ) ) / &
                  ( leaf_frac(:,j,m) + d_leaf_frac(:,j,m)- d_leaf_frac(:,j,m+1) )

          ENDWHERE

       ENDDO       ! Loop over age classes

	!! For last age class, there is no leaf fraction leaving the class. 

       WHERE ( d_leaf_frac(:,j,nleafages) .GT. min_stomate )

          leaf_age_new(:,nleafages) = ( ( leaf_frac(:,j,nleafages) * leaf_age(:,j,nleafages) )  + &
               ( d_leaf_frac(:,j,nleafages) * leaf_age(:,j,nleafages-1) ) ) / &
               ( leaf_frac(:,j,nleafages) + d_leaf_frac(:,j,nleafages) )

       ENDWHERE

       DO m = 2, nleafages       ! Loop over age classes

          WHERE ( d_leaf_frac(:,j,m) .GT. min_stomate )

             leaf_age(:,j,m) = leaf_age_new(:,m)

          ENDWHERE

       ENDDO       ! Loop over age classes

       !! 2.2.3 calculate new fraction

!!  yidi
       IF (is_oilpalm(j)) THEN
          WRITE(numout,*) 'yd: vmax4,before update  leaf_frac', leaf_frac(:,j,1:4),'nvm=',j
       ENDIF
!!  yidi
       DO m = 2, nleafages       ! Loop over age classes

          ! where the change comes from
          leaf_frac(:,j,m-1) = leaf_frac(:,j,m-1) - d_leaf_frac(:,j,m)

          ! where it goes to
          leaf_frac(:,j,m) = leaf_frac(:,j,m) + d_leaf_frac(:,j,m)

       ENDDO       ! Loop over age classes

!!  yidi
       IF (is_oilpalm(j)) THEN
          WRITE(numout,*) 'yd: vmax4,after update  leaf_frac', leaf_frac(:,j,1:4),'nvm=',j
       ENDIF
!!  yidi
       !! 2.2.4 renormalize fractions in order to prevent accumulation 
       !       of numerical errors

       ! correct small negative values

       DO m = 1, nleafages
          leaf_frac(:,j,m) = MAX( zero, leaf_frac(:,j,m) )
       ENDDO

       ! total of fractions, should be very close to one where there is leaf mass

       sumfrac(:) = zero

       DO m = 1, nleafages       ! Loop over age classes

          sumfrac(:) = sumfrac(:) + leaf_frac(:,j,m)

       ENDDO       ! Loop over age classes

       ! normalize

       DO m = 1, nleafages       ! Loop over age classes

          WHERE ( sumfrac(:) .GT. min_stomate )

             leaf_frac(:,j,m) = leaf_frac(:,j,m) / sumfrac(:) 

          ELSEWHERE

             leaf_frac(:,j,m) = zero

          ENDWHERE

       ENDDO       ! Loop over age classes

    ENDDO         ! Loop over PFTs

    !
    !! 3 calculate vmax as a function of the age
    !

    DO j = 2,nvm
       !!@xchen
       vcmax(:,j) = zero
       vcmax_cl1(:,j) = zero                                                        
       vcmax_cl2(:,j) = zero                                                        
       vcmax_cl3(:,j) = zero                                                        
       vcmax_cl4(:,j) = zero                                                        

       ! sum up over the different age classes
       IF (ok_dgvm .AND. pheno_type(j)==1 .AND. leaf_tab(j)==2) THEN
          ! pheno_typ=evergreen and leaf_tab=needleleaf
          vcmax(:,j) = Vcmax25(j)
          vcmax_cl1(:,j) = Vcmax25(j)                                                    
          vcmax_cl2(:,j) = Vcmax25(j)                                                    
          vcmax_cl3(:,j) = Vcmax25(j)                                                    
          vcmax_cl4(:,j) = Vcmax25(j)

       ELSE 
          ! for deciduous tree
          DO m = 1, nleafages       ! Loop over age classes

             !
             !! 3.1 efficiency in each of the age classes
             !!     it varies from vmax_offset to 1 
             !!     linearly increases from vmax_offset to 1 for 0 < rel_age < leafage_firstmax
             !!     is 1 when leafage_firstmax < rel_age < leafage_lastmax
             !!     linearly decreases from 1 to vmax_offset for leafage_lastmax < rel_age < leafage_firstmax
             !!     vmax_offset for rel_age >= leafage_old
             !!     (Ishida et al., 1999)
             rel_age(:) = leaf_age(:,j,m) / leafagecrit(j)

             leaf_efficiency(:) = MAX( vmax_offset(1), MIN( un, &
                  vmax_offset(1) + (un - vmax_offset(1)) * rel_age(:) / leafage_firstmax, &
                  un - (un - vmax_offset(2)) * ( rel_age(:) - leafage_lastmax ) / &
                  ( leafage_old - leafage_lastmax ) ) )

             !
             !! 3.2 add to mean vmax
             !             
             IF (ok_Nlim .OR. ok_LAIdev(j)) THEN
                vcmax(:,j) = vcmax(:,j) + Vcmax25(j) * N_limfert(:,j)*leaf_efficiency(:) * leaf_frac(:,j,m)
                !!@xchen calculate vcmax for different leaf age classes
                IF (m .EQ. 1) THEN
                  vcmax_cl1(:,j) = vcmax_cl1(:,j) + Vcmax25(j) * N_limfert(:,j)*leaf_efficiency(:)
                ELSE If (m .EQ. 2) THEN
                  vcmax_cl2(:,j) = vcmax_cl2(:,j) + Vcmax25(j) * N_limfert(:,j)*leaf_efficiency(:)
                ELSE If (m .EQ. 3) THEN                                             
                  vcmax_cl3(:,j) = vcmax_cl3(:,j) + Vcmax25(j) * N_limfert(:,j)*leaf_efficiency(:)
                ELSE                                              
                  vcmax_cl4(:,j) = vcmax_cl4(:,j) + Vcmax25(j) * N_limfert(:,j)*leaf_efficiency(:)
                ENDIF

             ELSE
                vcmax(:,j) = vcmax(:,j) + Vcmax25(j) * leaf_efficiency(:) * leaf_frac(:,j,m)
                IF (m .EQ. 1) THEN                                                   
                  vcmax_cl1(:,j) = vcmax_cl1(:,j) + Vcmax25(j) *leaf_efficiency(:)
                ELSE If (m .EQ. 2) THEN                                              
                  vcmax_cl2(:,j) = vcmax_cl2(:,j) + Vcmax25(j) *leaf_efficiency(:)
                ELSE If (m .EQ. 3) THEN                                             
                  vcmax_cl3(:,j) = vcmax_cl3(:,j) + Vcmax25(j) *leaf_efficiency(:)
                ELSE                                                            
                  vcmax_cl4(:,j) = vcmax_cl4(:,j) + Vcmax25(j) *leaf_efficiency(:)
                ENDIF
             ENDIF
          ENDDO ! loop over age classes
       ENDIF

    ENDDO       ! loop over PFTs

    IF (printlev>=4) WRITE(numout,*) 'Leaving vmax'

  END SUBROUTINE vmax

END MODULE stomate_vmax
