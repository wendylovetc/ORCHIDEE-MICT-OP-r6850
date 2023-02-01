! INITIALIZATION OF VARIABLES REGARDING THE CROP CYCLE
! 27/06/2013---xcw 

subroutine Stics_init(&
               kjpindex        ,&
               !nvm             ,&
               f_crop_init     ,&   
               f_crop_recycle      ,&   
               in_cycle                ,&   
               f_sen_lai                ,&   
               onarretesomcourdrp      ,&   
!               nlevobs                 ,&
!               namfobs                 ,&
!               nfloobs                 ,&
!               nlanobs                 ,&
!               nlaxobs                 ,&
!               nmatobs                 ,&
!               nrecobs                 ,&
!               nsenobs                 ,&
!               ndrpobs                 ,&
               nsendltams              ,&
               nsendltai               ,&
               nsenpfeuilverte         ,&
               nsendurvie              ,&
               nsenndurvie             ,&
               densiteequiv            ,&
               nplt                    ,&
               tursla                  ,&
               ssla                     ,&
               pfeuilverte             ,&
               bsenlai                 ,&
               zrac                    ,&
               nrec                    ,& 
               nlan                    ,&
               tcult                   ,&
               udevair                 ,&
               udevcult                ,&
               ndrp                    ,&
               rfvi                    ,&
               nlev                    ,&
               nger                    ,&
               etatvernal              ,&
               caljvc                  ,&
               rfpi                    ,&
               upvt                    ,&
               utp                     ,&
               somcour                 ,&
               somcourdrp              ,&
               somcourutp              ,&
               tdevelop                ,&
               somtemp                 ,&
               somcourfauche           ,&
               stpltger                ,&
               R_stamflax              ,&
               R_stlaxsen              ,&
               R_stsenlan              ,&
               stlevflo                ,&
               nflo                    ,&
               R_stlevdrp              ,&
               R_stflodrp              ,&
               R_stdrpmat              ,&
               nmat                    ,&
               nlax                    ,&
               nrecbutoir              ,&
               group                   ,&
               ndebdes                 ,&
               R_stdrpdes              ,&
               densite                 ,&
               densitelev              ,&
               coeflev                 ,&
               densiteger              ,&
               somelong                 ,&
               somger                  ,&
               humectation             ,&
               nbjhumec                ,&
               somtemphumec            ,&
               stpltlev                ,&
               namf                    ,&
               stmatrec                ,&
               tustress                ,&
               slai                     ,&
               somfeuille              ,&
               pdlai                   ,&
               nbfeuille               ,&
               reajust                 ,&
               ulai                    ,&
               pdulai                  ,&
               efdensite               ,&
               tempeff                 ,&
               nstopfeuille            ,&
               deltai                  ,&
               svmax                    ,&
               nsen                    ,&
               laisen                  ,&
               pdlaisen                ,&
               dltaisenat              ,&
               nsencour                ,&
               dltamsen                ,&
               dltaisen                ,&
               fgellev                 ,&
               gelee                   ,&
               fstressgel              ,&
               R_stlevamf              ,&
               dernier_n               ,&
               durvieI                 ,&
               durvie                  ,&
               ndebsen                 ,&
               somsenreste             ,&
               shumrel                  ,&
               swfac                   ,&
               turfac                  ,&
               senfac                  ,&  
               mafeuiljaune            ,&
               msneojaune              ,&              
               v_dltams                ,&
               fgelflo                 ,&
               pdircarb                ,&
               ircarb                  ,&
               nbgrains                ,&
               pgrain                  ,&
               vitmoy                  ,&
               nbgraingel              ,&
               pgraingel               ,&
               dltags                  ,&
               ftempremp               ,&
               magrain                 ,&
               pdmagrain               ,&
               nbj0remp                ,&
               pdsfruittot             ,&
               repracmax               ,&
               repracmin               ,&
               kreprac                 ,&
               somtemprac              ,&
               urac                    ,&
               reprac                  ,&
               nstoprac                ,&
               c_reserve               ,&
               c_leafb                 ,&
               !biomass                 ,&
               deltgrain               ,&
               gslen                   ,&
               drylen, &
               histgrowthset, hist_sencourset, hist_latestset, doyhiststset, &
               nboxmax,box_ulai, box_ndays, box_lai, box_lairem, box_tdev, box_biom, box_biomrem, box_durage, box_somsenbase)                


  !USE stomate_io
  USE pft_parameters
  USE constantes


  implicit none

  ! DECLARATION
  integer,intent(in)                             ::         kjpindex
  !integer,intent(in)                             ::         nvm
  integer,intent(in)                             ::         nboxmax
  logical, intent(in)::         f_crop_init        
  logical,dimension(kjpindex, nvm), intent(inout)::         f_crop_recycle         
  logical,dimension(kjpindex, nvm), intent(out)::         in_cycle         
  logical,dimension(kjpindex, nvm), intent(out)::         f_sen_lai         
  logical,dimension(kjpindex, nvm), intent(out)::         onarretesomcourdrp         
!  integer,dimension(kjpindex, nvm), intent(out)::         nlevobs                 
!  integer,dimension(kjpindex, nvm), intent(out)::         namfobs                 
!  integer,dimension(kjpindex, nvm), intent(out)::         nfloobs                 
!  integer,dimension(kjpindex, nvm), intent(out)::         nlanobs                 
!  integer,dimension(kjpindex, nvm), intent(out)::         nlaxobs                 
!  integer,dimension(kjpindex, nvm), intent(out)::         nmatobs                 
!  integer,dimension(kjpindex, nvm), intent(out)::         nrecobs                 
!  integer,dimension(kjpindex, nvm), intent(out)::         nsenobs                 
!  integer,dimension(kjpindex, nvm), intent(out)::         ndrpobs                 
  !
  real,dimension(kjpindex, nvm), intent(out)::         nsendltams              
  real,dimension(kjpindex, nvm), intent(out)::         nsendltai               
  real,dimension(kjpindex, nvm), intent(out)::         nsenpfeuilverte         
  real,dimension(kjpindex, nvm), intent(out)::         nsendurvie              
  real,dimension(kjpindex, nvm), intent(out)::         nsenndurvie             
  real,dimension(kjpindex, nvm), intent(out)::         densiteequiv            
  integer,dimension(kjpindex, nvm), intent(out)::         nplt                    
  real,dimension(kjpindex, nvm), intent(out)::         tursla                  
  real,dimension(kjpindex, nvm), intent(out)::         ssla                     
  real,dimension(kjpindex, nvm), intent(out)::         pfeuilverte             
  real,dimension(kjpindex, nvm), intent(out)::         bsenlai                 
  
  ! STICS::LAIdev::DEVELOPMENT
  real,dimension(kjpindex, nvm), intent(out)::         zrac                    
  integer,dimension(kjpindex, nvm), intent(out)::         nrec                     
  integer,dimension(kjpindex, nvm), intent(out)::         nlan                    
  real,dimension(kjpindex, nvm), intent(out)::         tcult                   
  real,dimension(kjpindex, nvm), intent(out)::         udevair                 
  real,dimension(kjpindex, nvm), intent(out)::         udevcult                
  integer,dimension(kjpindex, nvm), intent(out)::         ndrp                    
  real,dimension(kjpindex, nvm), intent(out)::         rfvi                    
  integer,dimension(kjpindex, nvm), intent(out)::         nlev                    
  integer,dimension(kjpindex, nvm), intent(out)::         nger                    
  logical,dimension(kjpindex, nvm), intent(out)::         etatvernal              
  real,dimension(kjpindex, nvm), intent(out)::         caljvc                  
  real,dimension(kjpindex, nvm), intent(out)::         rfpi                    
  real,dimension(kjpindex, nvm), intent(out)::         upvt                    
  real,dimension(kjpindex, nvm), intent(out)::         utp                     
  real,dimension(kjpindex, nvm), intent(out)::         somcour                 
  real,dimension(kjpindex, nvm), intent(out)::         somcourdrp              
  real,dimension(kjpindex, nvm), intent(out)::         somcourutp              
  real,dimension(kjpindex, nvm), intent(out)::         tdevelop                
  real,dimension(kjpindex, nvm), intent(out)::         somtemp                 
  real,dimension(kjpindex, nvm), intent(out)::         somcourfauche           
  real,dimension(kjpindex, nvm), intent(out)::         stpltger                
  real,dimension(kjpindex, nvm), intent(out)::         R_stamflax              
  real,dimension(kjpindex, nvm), intent(out)::         R_stlaxsen              
  real,dimension(kjpindex, nvm), intent(out)::         R_stsenlan              
  real,dimension(kjpindex, nvm), intent(out)::         stlevflo                
  integer,dimension(kjpindex, nvm), intent(out)::         nflo                    
  real,dimension(kjpindex, nvm), intent(out)::         R_stlevdrp              
  real,dimension(kjpindex, nvm), intent(out)::         R_stflodrp              
  real,dimension(kjpindex, nvm), intent(out)::         R_stdrpmat              
  integer,dimension(kjpindex, nvm), intent(out)::         nmat                    
  integer,dimension(kjpindex, nvm), intent(out)::         nlax                    
  integer,dimension(kjpindex, nvm), intent(out)::         nrecbutoir              
  real,dimension(kjpindex, nvm), intent(out)::         group                   
  integer,dimension(kjpindex, nvm), intent(out)::         ndebdes                 
  real,dimension(kjpindex, nvm), intent(out)::         R_stdrpdes              
  real,dimension(kjpindex, nvm), intent(out)::         densite                 
  real,dimension(kjpindex, nvm), intent(out)::         densitelev
  real,dimension(kjpindex, nvm), intent(out)::         coeflev
              
  real,dimension(kjpindex, nvm), intent(out)::         densiteger              
  real,dimension(kjpindex, nvm), intent(out)::         somelong                 
  real,dimension(kjpindex, nvm), intent(out)::         somger                  
  logical,dimension(kjpindex, nvm), intent(out)::         humectation             
  integer,dimension(kjpindex, nvm), intent(out)::         nbjhumec                
  real,dimension(kjpindex, nvm), intent(out)::         somtemphumec            
  real,dimension(kjpindex, nvm), intent(out)::         stpltlev                
  integer,dimension(kjpindex, nvm), intent(out)::         namf                    
  real,dimension(kjpindex, nvm), intent(out)::         stmatrec                
  ! STICS::LAIdev:: LAI calculation
  real,dimension(kjpindex, nvm), intent(out)::         tustress                
  real,dimension(kjpindex, nvm), intent(out)::         slai                     
  real,dimension(kjpindex, nvm), intent(out)::         somfeuille              
  real,dimension(kjpindex, nvm), intent(out)::         pdlai                   
  integer,dimension(kjpindex, nvm), intent(out)::         nbfeuille               
  real,dimension(kjpindex, nvm), intent(out)::         reajust                 
  real,dimension(kjpindex, nvm), intent(out)::         ulai                    
  real,dimension(kjpindex, nvm), intent(out)::         pdulai                  
  real,dimension(kjpindex, nvm), intent(out)::         efdensite               
  real,dimension(kjpindex, nvm), intent(out)::         tempeff                 
  integer,dimension(kjpindex, nvm), intent(out)::         nstopfeuille            
  real,dimension(kjpindex, nvm), intent(out)::         deltai                  
  real,dimension(kjpindex, nvm), intent(out)::         svmax                    
  integer,dimension(kjpindex, nvm), intent(out)::         nsen                    
  real,dimension(kjpindex, nvm), intent(out)::         laisen                  
  real,dimension(kjpindex, nvm), intent(out)::         pdlaisen                
  real,dimension(kjpindex, nvm), intent(out)::         dltaisenat              
  ! STICS:: LAIdev:: LAI senescence 
  integer,dimension(kjpindex, nvm), intent(out)::         nsencour                
  real,dimension(kjpindex, nvm), intent(out)::         dltamsen                
  real,dimension(kjpindex, nvm), intent(out)::         dltaisen                
  real,dimension(kjpindex, nvm), intent(out)::         fgellev                 
  logical,dimension(kjpindex, nvm), intent(out)::         gelee                   
  real,dimension(kjpindex, nvm), intent(out)::         fstressgel              
  real,dimension(kjpindex, nvm), intent(out)::         R_stlevamf              
  integer,dimension(kjpindex, nvm), intent(out)::         dernier_n               
  real,dimension(kjpindex, nvm), intent(out)::         durvieI                 
  real,dimension(kjpindex, nvm), intent(out)::         durvie                  
  integer,dimension(kjpindex, nvm), intent(out)::         ndebsen                 
  real,dimension(kjpindex, nvm), intent(out)::         somsenreste
  ! STICS:: LAIdev:: STRESS             
  real,dimension(kjpindex, nvm), intent(out)::         shumrel                  
  real,dimension(kjpindex, nvm), intent(out)::         swfac                   
  real,dimension(kjpindex, nvm), intent(out)::         turfac                  
  real,dimension(kjpindex, nvm), intent(out)::         senfac                
 
  real,dimension(kjpindex, nvm), intent(out)::         mafeuiljaune                
  real,dimension(kjpindex, nvm), intent(out)::         msneojaune                
  ! STICS:: CARBON ALLOCATION

  real,dimension(kjpindex, nvm, 60), intent(out)::         v_dltams                
  real,dimension(kjpindex, nvm), intent(out)::         fgelflo                
  real,dimension(kjpindex, nvm), intent(out)::         pdircarb                
  real,dimension(kjpindex, nvm), intent(out)::         ircarb                
  real,dimension(kjpindex, nvm), intent(out)::         nbgrains                
  real,dimension(kjpindex, nvm), intent(out)::         pgrain              
  real,dimension(kjpindex, nvm), intent(out)::         vitmoy                
  real,dimension(kjpindex, nvm), intent(out)::         nbgraingel                
  real,dimension(kjpindex, nvm), intent(out)::         pgraingel                
  real,dimension(kjpindex, nvm), intent(out)::         dltags                
  real,dimension(kjpindex, nvm), intent(out)::         ftempremp                
  real,dimension(kjpindex, nvm), intent(out)::         magrain                
  real,dimension(kjpindex, nvm), intent(out)::         pdmagrain               
  integer,dimension(kjpindex, nvm), intent(out)::         nbj0remp                 
  real,dimension(kjpindex, nvm), intent(out)::         pdsfruittot               
  real,dimension(kjpindex, nvm), intent(out)::         repracmax               
  real,dimension(kjpindex, nvm), intent(out)::         repracmin               
  real,dimension(kjpindex, nvm), intent(out)::         kreprac               
  real,dimension(kjpindex, nvm), intent(out)::         somtemprac               
  real,dimension(kjpindex, nvm), intent(out)::         urac              
  real,dimension(kjpindex, nvm), intent(out)::         reprac               

  integer,dimension(kjpindex, nvm), intent(out)::         nstoprac                 
  real,dimension(kjpindex, nvm), intent(out)::         c_reserve               
  real,dimension(kjpindex, nvm), intent(out)::         c_leafb            
  !real,dimension(kjpindex, nvm, nparts ), intent(out)::         biomass           
  real,dimension(kjpindex, nvm), intent(out)::         deltgrain          
  integer,dimension(kjpindex, nvm), intent(out)::         gslen                 
  integer,dimension(kjpindex, nvm), intent(out)::         drylen                 
  real,dimension(kjpindex, nvm, 300, 5), intent(out) :: histgrowthset
  integer, dimension(kjpindex,nvm), intent(out) :: hist_sencourset
  integer, dimension(kjpindex,nvm), intent(out) :: hist_latestset
  integer, dimension(kjpindex,nvm), intent(out) :: doyhiststset

  integer, dimension(kjpindex,nvm,nboxmax), intent(out) :: box_ndays
  real, dimension(kjpindex,nvm,nboxmax), intent(out)    :: box_lai
  real, dimension(kjpindex,nvm,nboxmax), intent(out)    :: box_lairem
  real, dimension(kjpindex,nvm,nboxmax), intent(out)    :: box_tdev
  real, dimension(kjpindex,nvm,nboxmax), intent(out)    :: box_biom
  real, dimension(kjpindex,nvm,nboxmax), intent(out)    :: box_biomrem
  real, dimension(kjpindex,nvm,nboxmax), intent(out)    :: box_durage
  real, dimension(kjpindex,nvm,nboxmax), intent(out)    :: box_somsenbase
  real, dimension(nvm,nboxmax), intent(out)             :: box_ulai
  ! local variables
  integer    :: j, ip, k, mid
  real       :: uinflex

 
  ! 
  DO j= 1, nvm
     
     DO ip = 1, kjpindex 
       IF (f_crop_init .OR. f_crop_recycle(ip, j)) THEN  
          

          ! LAIdev :: FORCED RUN VARIABLES
          in_cycle(ip, j) = .FALSE.
          f_sen_lai(ip, j) = .TRUE.
          onarretesomcourdrp(ip, j) = .TRUE.
!          nlevobs(ip, j) = 999
!          namfobs(ip, j) = 999   
!          nfloobs(ip, j) = 999
!          nlanobs(ip, j) = 999
!          nlaxobs(ip, j) = 999
!          nmatobs(ip, j) = 999
!          nrecobs(ip, j) = 999
!          nsenobs(ip, j) = 999
!          ndrpobs(ip, j) = 999
        
          ! LAIdev ::  SPECIFIC
          nsendltams(ip, j) = 0.
          nsendltai(ip, j) = 0.
          nsenpfeuilverte(ip, j) = 0.
          nsendurvie(ip, j) = 0.
          nsenndurvie(ip, j) = 0.
          densiteequiv(ip, j) = 0.
          nplt(ip, j) = 0
          tursla(ip, j) = 1.
          ssla(ip, j) = 0.
          pfeuilverte(ip, j) = 0.
          bsenlai(ip, j) = 0.

          ! LAIdev ::  DEVELOPMENT 
          zrac(ip, j) = 0.
          nrec(ip, j) = 0
          nlan(ip, j) = 0
          tcult(ip, j) = 0.
          udevair(ip, j) = 0.
          udevcult(ip, j) = 0.
          ndrp(ip, j) = 0
          rfvi(ip, j) = 0.
          nlev(ip, j) = 0
          nger(ip, j) = 0
          etatvernal(ip, j) = .FALSE.
          caljvc(ip, j) = 0.
          rfpi(ip, j) = 0.
          upvt(ip, j) = 0.
          utp(ip, j) = 0.
          somcour(ip, j) = 0.
          somcourdrp(ip, j) = 0.
          somcourutp(ip, j) = 0.
          tdevelop(ip, j) = 0.
          somtemp(ip, j) = 0.
          somcourfauche(ip, j) = 0.
          stpltger(ip, j) = SP_stpltger(j)
          R_stamflax(ip, j) = SP_stamflax(j) 
          R_stlaxsen(ip, j) = SP_stlaxsen(j)
          R_stsenlan(ip, j) = SP_stsenlan(j)
          stlevflo(ip, j) = SP_stlevdrp(j) - SP_stflodrp(j)
          nflo(ip, j) = 0
          R_stlevdrp(ip, j) = SP_stlevdrp(j)
          R_stflodrp(ip, j) = SP_stflodrp(j)
          R_stdrpmat(ip, j) = SP_stdrpmat(j)
          nmat(ip, j) = 0
          nlax(ip, j) = 0
          nrecbutoir(ip, j) = 999
          group(ip, j) = 0.
          ndebdes(ip, j) = 0
          R_stdrpdes(ip, j) = SP_stdrpdes(j)
          densite(ip, j) = 0.
          densitelev(ip, j) = 0.
          coeflev(ip, j) = 1.
          densiteger(ip, j) = 0.
          somelong(ip, j) = 0.
          somger(ip, j) = 0.
          humectation(ip, j) = .FALSE.
          nbjhumec(ip, j) = 0
          somtemphumec(ip, j) = 0.
          stpltlev(ip, j) = 0.
          namf(ip, j) = 0
          stmatrec(ip, j) = 0.
          
          ! LAIdev ::  LAI calculation
          tustress(ip, j) = 1.
          slai(ip, j) = 0.
          somfeuille(ip, j) = 0.
          pdlai(ip, j) = 0.
          nbfeuille(ip, j) = 0
          reajust(ip, j) = 0.
          ulai(ip, j) = 0.
          pdulai(ip, j) = 0.
          efdensite(ip, j) = 0.
          tempeff(ip, j) = 0.
          nstopfeuille(ip, j) = 0
          deltai(ip, j) = 0.
          svmax(ip, j) = 0.
          nsen(ip, j) = 0
          laisen(ip, j) = 0.
          pdlaisen(ip, j) = 0.
          dltaisenat(ip, j) = 0.
          
          ! LAIdev :: LAI senescence
          nsencour(ip, j) = 0
          dltamsen(ip, j) = 0.
          dltaisen(ip, j) = 0.
          fgellev(ip, j) = 0.
          gelee(ip, j) = .FALSE.
          fstressgel(ip, j) = 0.
          R_stlevamf(ip, j) = SP_stlevamf(j)
          dernier_n(ip, j) = 0
          durvieI(ip, j) = 0.
          durvie(ip, j) = 0.
          ndebsen(ip, j) = 0
          somsenreste(ip, j) = 0.
          
          IF (any(ok_LAIdev(:)) .AND. any(SP_codlainet(:) == 3)) THEN
!              write(*,*) 'xuhui: box module to be initialize ', j
              box_ndays(ip,j,:) = 0
              box_lai(ip,j,:) = 0.
              box_lairem(ip,j,:) = 0.
              box_tdev(ip,j,:) = 0.
              box_biom(ip,j,:) = 0.
              box_biomrem(ip,j,:) = 0.
              box_durage(ip,j,:) = 0.
              box_somsenbase(ip,j,:) = 0.
!              write(*,*) 'xuhui: box module initialized ', j
          ENDIF
          ! LAIdev :: STRESS
           
          shumrel(ip, j) = 0.
          swfac(ip, j) = 1.
          turfac(ip, j) = 1.
          senfac(ip, j) = 1.

          mafeuiljaune(ip, j) = 0.
          msneojaune(ip, j) = 0.
          ! STICS: CARBON ALLOCATION 

          v_dltams(ip, j, :) = 0.
          fgelflo(ip, j) = 1.
          pdircarb(ip, j) = 0.
          ircarb(ip, j) = 0.
          nbgrains(ip, j) = 0.
          pgrain(ip, j) = 0.
          vitmoy(ip, j) = 0.
          nbgraingel(ip, j) = 0.
          pgraingel(ip, j) = 0.
          dltags(ip, j) = 0.
          ftempremp(ip, j) = 0.
          magrain(ip, j) = 0.
          pdmagrain(ip, j) = 0.
          nbj0remp(ip, j) = 0
          pdsfruittot(ip, j) = 0.
          repracmax(ip, j) = 0.
          repracmin(ip, j) = 0.
          kreprac(ip, j) = 0.
          somtemprac(ip, j) = 0.
          urac(ip, j) = 0.
          reprac(ip, j) = 0.
          nstoprac(ip, j) = 0
          c_reserve(ip, j) = 0. 
          c_leafb(ip, j)= 0. 
          deltgrain(ip, j) = 0.
          gslen(ip, j) = 0 
          drylen(ip, j) = 0 
          ! FINISH THE INITIALIZATION AND RESET THE VALUE
          IF (any(ok_LAIdev(:)) .AND. any(SP_codlainet(:) == 2)) THEN
              histgrowthset(ip,j,:,:) = 0.
              hist_sencourset(ip,j) = 0
              hist_latestset(ip,j) = 0
              doyhiststset(ip,j) = 0
          ENDIF
          f_crop_recycle(ip, j) = .FALSE.

        ENDIF
     ENDDO ! ip
  ENDDO ! j
    
!    write(*,*) 'xuhui: initialized boxulai'
    IF (any(ok_LAIdev(:)) .AND. any(SP_codlainet(:)==3)) THEN
        box_ulai(:,:) = 0.
        DO j=1,nvm
    !        IF (f_crop_init) THEN
!            IF (f_crop_init .OR. any(f_crop_recycle(:, j))) THEN
                if (SP_nbox(j) < 2) then
!                    write(*,*) 'xuhui: box_ulai not allocated for ', j, 'with nbox ', SP_nbox(j)
                else    
                    uinflex = SP_vlaimax(j)
                    mid = floor(real(SP_nbox(j))/2)
                    if (mid<1) then
                        mid=1
                    elseif (mid .eq. SP_nbox(j)) then
                        mid=SP_nbox(j)-1
                    endif
                    do k=1,mid
                        box_ulai(j,k) = 1+(k-1)*((uinflex-1)/real(mid))
                    enddo
                    do k=mid+1,SP_nbox(j)
                        box_ulai(j,k) = uinflex + (k-mid-1)*(3-uinflex)/real(SP_nbox(j)-mid)
                    enddo
!                    write(*,*) 'xuhui: box_ulai for pft ', j, ': ',box_ulai(j,:)
                endif
                
 !           ENDIF
        ENDDO
    ENDIF

end subroutine Stics_init

