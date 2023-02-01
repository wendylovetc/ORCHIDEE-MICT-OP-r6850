! =================================================================================================================================
! MODULE       : stomate_stics
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Groups the subroutines which interface with or are related with sticslai
!
!!
!!\n DESCRIPTION : None
!!
!! RECENT CHANGE(S) : None
!!
!! REFERENCE(S)	: None
!!
!! SVN :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_stomate/stomate.f90 $
!! $Date: 2017-07-17 16:22:36 +0200 (Mon, 17 Jul 2017) $
!! $Revision: 4509 $
!! \n
!_ ================================================================================================================================

MODULE stomate_stics

  ! Modules used:
  USE netcdf
  USE defprec
  USE grid
  USE constantes
  USE constantes_soil
  USE pft_parameters
  USE ioipsl
  USE ioipsl_para 
  USE mod_orchidee_para
  USE interpol_help  ! necessary for management map input
  USE time, ONLY: year_length_in_days 

  IMPLICIT NONE

  ! Private & public routines
  PRIVATE
  PUBLIC stomate_stics_ManageInput, stomate_stics_read_cycle, stomate_stics_read_rotation, stomate_stics_read_plantdate, stomate_stics_get_cmd
  PUBLIC stomate_stics_NfertInput,  stomat_stics_calc_N_limfert,     stomate_stics_rotation

  ! Public variables (not recommended)
  !PUBLIC 

CONTAINS
! 
!! ================================================================================================================================
!! SUBROUTINE 	: stomate_forcing_CropVar
!!
!>\BRIEF        Interpolate (extract) Crop Variety information for each crop type
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
  SUBROUTINE slowproc_CropVar(nbpt, lalo, neighbours, resolution, contfrac,&  ! In
                              plantcyc, cropvar ) ! Out
    ! 
    ! 
    ! INTEGER(i_std)            :: iv
    ! INTEGER(i_std),DIMENSION(nvm) :: plantdate_default = (/(0,iv=1,nvm)/)
    ! !It should be defined as an array for different crops
    ! !The value should be given outside this program and subject to crop type
    ! !if being defined outside, please comment the above sentence
    
    ! 0.1 INPUT
    ! 
    INTEGER(i_std), INTENT(in)  :: nbpt         ! Number of points for which the data needs to be interpolated (extracted)
    REAL(r_std), INTENT(in) :: lalo(nbpt,2)     ! Vector of latitude and longtitude
    INTEGER(i_std), INTENT(in)  :: neighbours(nbpt,8)   ! Vectors of neighbours for each grid point
    ! (1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW)
    !REAL(r_std),INTENT(in)  :: resolution(nbpt,2)   ! The size in km of each grid box in lat and lon
    REAL(r_std)             :: resolution(nbpt,2)   ! The size in km of each grid box in lat and lon
    REAL(r_std),INTENT(in)  :: contfrac(nbpt)   ! The fraction of land in each grid box
    ! INTEGER(I_std),INTENT(in) :: plantcyc         ! The time of the rotation
    ! in this year (0 means this info unknown)
    ! INTEGER(i_std), INTENT(in)    :: croptype     ! The type of crop being
    ! simulated
    ! not necessary as we read planting date data for all crop types
    ! 
    ! 0.2 OUTPUT
    ! 
    REAL(r_std),ALLOCATABLE, INTENT(out)    :: cropvar(:,:,:)   ! The crop variety matrix
    INTEGER(i_std), INTENT(out) :: plantcyc     ! times of rotation that year
    ! nvm is the number of PFTs, there may not be planting date for all the PFTs
    ! 
    ! 0.3 LOCAL
    ! 
    INTEGER(i_std)      :: nbvmax       ! a parameter for interpolation
    CHARACTER(LEN=80)       :: filename
    INTEGER(i_std)      :: iml, jml, lml, tml, fid, fid1
    INTEGER(i_std)      :: ip, jp, ib, ilf, fopt ! for-loop variable
    INTEGER(i_std)      :: nbexp
    REAL(r_std)         :: lev(1), date, dt
    REAL(r_std)         :: missing_val
    INTEGER(i_std)      :: itau(1)

    INTEGER(i_std)      :: nb_dim
    INTEGER,DIMENSION(flio_max_var_dims) :: l_d_w
    LOGICAL         :: l_ex
    
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)    :: lat_rel, lon_rel
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:,:,:)  :: cropvar_mat ! LON LAT VEGET PLANTCYC TIME_COUNTER
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:,:)    :: cropvar_mat_1 ! LON LAT VEGET, PLANTCYC it is used when input file does not have time counter
    ! INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:)   :: plntdt_mat
    ! Because flinget does not support reading integer matrix, plntdt_mat has to
    ! be real
    ! or we will have to rewrite IOIPSL
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:) :: mask
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)    :: temp_var
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)    :: sub_area
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:,:)   :: sub_index

    REAL(r_std) :: sgn, sum_float
    INTEGER(i_std) :: ivgt, icyc, pltcyc, it
    CHARACTER(LEN=30) :: callsign
    LOGICAL :: ok_interpol
    INTEGER :: ALLOC_ERR
    INTEGER(i_std) :: cps,maxcp,jlf,maxps
    INTEGER(i_std),ALLOCATABLE,DIMENSION(:) :: hashlst1
    REAL(r_std),ALLOCATABLE,DIMENSION(:) :: hashlst2
    REAL(r_std) :: areamax
    LOGICAL ok_var  


    ! croptype = TRIM(croptype) !if croptype is a string
    ! else a switch expression is needed
    filename = "/work/cont003/p529tan/WXH/CropVar_1990.nc" ! default input file
    ! String operation needed
    CALL getin_p('CROPVAR_FILE',filename)

    IF (is_root_prc) THEN
    ! ? what does is_root_prc mean?
        CALL flininfo(filename, iml, jml, lml, tml, fid)
        ! CALL flinclo(fid)
    ENDIF
    CALL bcast(iml)
    CALL bcast(jml)
    ! CALL bcast(lml)
    ! CALL bcast(tml)
    
    ! Printing information for debugging
    WRITE(numout, *) "Xuhui's debug info for slowproc_CropVar #1:"
    WRITE(numout, *) "filename is: ", filename
    WRITE(numout, *) "Dimension 1, lon, iml:", iml
    WRITE(numout, *) "Dimension 2, lat, jml:", jml
    WRITE(numout, *) "Dimension 3, veget, lml:", lml
    WRITE(numout, *) "Dimension 4, time, tml:", tml
    ! apparently, flinget function is not designed to take veget but levels to
    ! be the 
    ! 3rd dimension, modification to lml is needed
    CALL flioopfd(filename,fid1)
    CALL flioinqv(fid1,v_n="CROPVAR", l_ex = l_ex, nb_dims = nb_dim, len_dims =l_d_w)
    IF (lml == 0) THEN
        ! CALL
        ! flioinqv(fid1,v_n="PLNTDT",l_ex=l_ex,nb_dims=nb_dim,len_dims=l_d_w)
        lml=l_d_w(3)
        WRITE(numout, *) "len_dims: ", l_d_w
        WRITE(numout, *) "lml AFTER revision"
        WRITE(numout, *) "lml: ", lml
    ENDIF
    IF (nb_dim.EQ.5) THEN
        plantcyc = l_d_w(4)
        tml = l_d_w(5)
        WRITE(numout,*) 'No. of plantcycle in the input: ', plantcyc
        WRITE(numout,*) 'tml changed to: ', tml
    ELSE
        IF (nb_dim.EQ.4) THEN
            plantcyc = 0
        ELSE
            WRITE(numout,*) "Error: More axis in PlantDate than expected"
            STOP
        ENDIF
    ENDIF
    IF (plantcyc>3) THEN
        WRITE(numout,*) "Rotation cycle more than accepted: ",plantcyc
    ENDIF
    WRITE(numout,*) "plantcyc: ",plantcyc
    CALL flioclo(fid1)
    CALL bcast(lml)
    CALL bcast(tml)
    ! CALL bcast(plantcyc)
    
    ALLOC_ERR=-1
    ALLOCATE(cropvar(nbpt,lml,tml),STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
        WRITE(numout,*) "ERROR IN ALLOCATION OF cropvar: ", ALLOC_ERR
    ENDIF
    WRITE(numout,*) "cropvar ALLOCATED"
    ! 
    !! assumption: the crop variety file is a 0.5 degree planting date (julian
    !date) file
    ! 
    ALLOC_ERR=-1
    ALLOCATE(lat_rel(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR/=0) THEN
        WRITE(numout,*) "ERROR IN ALLOCATION of lat_rel : ",ALLOC_ERR
        STOP 
    ENDIF
    ALLOC_ERR=-1
    ALLOCATE(lon_rel(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
        WRITE(numout,*) "ERROR IN ALLOCATION of lon_rel : ",ALLOC_ERR
        STOP 
    ENDIF
    ALLOC_ERR=-1
    ALLOCATE(mask(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
        WRITE(numout,*) "ERROR IN ALLOCATION of mask : ",ALLOC_ERR
        STOP 
    ENDIF
    
    IF (plantcyc == 0) THEN
        ALLOC_ERR=-1
        ALLOCATE(cropvar_mat_1(iml,jml,lml,tml), STAT=ALLOC_ERR) 
        ! !lml is supposed to be nvm (number of PFTs), if not ,change it
        IF (ALLOC_ERR/=0) THEN
            WRITE(numout,*) "ERROR IN ALLOCATION of cropvar_mat_1 : ",ALLOC_ERR
            STOP 
        ENDIF
        ALLOC_ERR=-1
        ALLOCATE(cropvar_mat(iml,jml,lml,1,tml), STAT=ALLOC_ERR)
        IF (ALLOC_ERR/=0) THEN
            WRITE(numout,*) "ERROR IN ALLOCATION of cropvar_mat: ", ALLOC_ERR
            STOP
        ENDIF
    ELSE
        ALLOC_ERR=-1
        ALLOCATE(cropvar_mat(iml,jml,lml,plantcyc,tml), STAT=ALLOC_ERR)
        IF (ALLOC_ERR/=0) THEN
            WRITE(numout,*) "ERROR IN ALLOCATION of cropvar_mat: ", ALLOC_ERR
            STOP
        ENDIF
    ENDIF
    
    ! input of some attributes
    ! IF (is_root_prc) CALL flinopen(filename, .FALSE., iml, jml, lml, lon_rel,
    ! lat_rel, lev, tml, itau, date, dt, fid)
    ! CALL bcast(lon_rel)
    ! CALL bcast(lat_rel)
    ! CALL bcast(itau)
    ! CALL bcast(date)
    ! CALL bcast(dt)
    IF (is_root_prc) THEN
        CALL flinget(fid, 'LON', iml, jml, lml, tml, 1, 1, lon_rel)
        CALL flinget(fid, 'LAT', iml, jml, lml, tml, 1, 1, lat_rel)
        CALL bcast(lon_rel)
        CALL bcast(lat_rel)
        ! WRITE (numout,*) 'lon_rel size: ', SIZE(lon_rel)
        ! WRITE (numout,*) 'lat_rel size: ', SIZE(lat_rel)
    ENDIF

    ! input of the matrix
    IF (is_root_prc) THEN 
        ! CALL flinget(fid, 'CROPVAR', iml, jml, lml, tml, 1, 1, crop_mat) 
        ! time_counter has to be 1, or it will not match the size of plntdt_mat
        CALL flioopfd(filename,fid1)
        IF (plantcyc>3) THEN
            WRITE(numout,*) "slowproc_CropVar : Planting cycle more than model accepted"
            STOP
        ELSE
            IF (plantcyc > 0) THEN
                CALL fliogetv(fid1,'CROPVAR',cropvar_mat,start=(/1,1,1,1,1/),count=(/iml,jml,lml,plantcyc,tml/))
            ELSEIF (plantcyc == 0) THEN
                CALL fliogetv(fid1,'CROPVAR',cropvar_mat_1,start=(/1,1,1,1/),count=(/iml,jml,lml,tml/))
            ELSE
                WRITE(numout,*) "invalid plantcyc axis : ", plantcyc
            ENDIF
        ENDIF
        ! get missing_val
        CALL fliogeta(fid1,'CROPVAR','missing_value',missing_val)
        CALL flioclo(fid1)
    ENDIF
    IF (plantcyc == 0) THEN
        cropvar_mat(:,:,:,1,:) = cropvar_mat_1(:,:,:,:)
        ! CALL bcast(cropvar_mat)
        DEALLOCATE(cropvar_mat_1)
    ELSE
        ! CALL bcast(cropvar_mat)
    ENDIF
    ! WRITE(numout,*) 'cropvar_mat size: ',SIZE(cropvar_mat)
    ! WRITE(numout,*) 'missing value: ', missing_val
    ! WRITE(numout,*) 'lat(361,284): ',lat_rel(361,284)
    ! WRITE(numout,*) 'lon(361,284): ',lon_rel(361,284)
    ! WRITE(numout,*) 'cropvar(361,284,1,1): ',cropvar_mat(361,284,1,1)
    
    IF (is_root_prc) CALL flinclo(fid)
    
    cropvar(:,:,:) = zero ! nbpt veget time
    IF (plantcyc == 0) THEN
        pltcyc = 1
    ELSE
        pltcyc = plantcyc
    ENDIF
    
    DO it = 1,tml
        DO icyc = 1,1 ! pltcyc ! not dealing with plant cycles for now
            DO ivgt = 1,lml ! ? We can suppose PFTs less than 10 are natural veg without planting date, but not now
                IF (natural(ivgt)) THEN
                    WRITE(numout,*) "veget, plantcyc, time: ", ivgt,icyc,it
                    nbexp = 0
                    ! the number of exceptions
                    
                    ! mask of available value
                    !
                    ! mask is commented because it is of no use right now.
                    mask(:,:) = zero;  ! Defined in constante.f90
                    DO ip = 1,iml
                        DO jp = 1,jml
                            IF ((cropvar_mat(ip,jp,ivgt,icyc,it) .GT. min_sechiba) .AND.  &
                            (cropvar_mat(ip,jp,ivgt,icyc,it) /= missing_val)) THEN
                                mask(ip,jp) = un;  ! Defined in constante.f90
                                ! here we assumed that for each plant cycle at each 
                                ! there might be missing data at different grid
                                ! in this case, mask has to be generated each plant
                                ! cycle each time step
                            ENDIF
                        ENDDO
                    ENDDO
                    
                    ! Interpolation started
                    nbvmax = 200
                    ! the maximum amount of fine grids that one coarse grid may have
                    
                    callsign = "Crop Variety"
                    
                    ok_interpol = .FALSE.
                    
                    DO WHILE ( .NOT. ok_interpol )
                        WRITE(numout,*) "Pojection arrays for ", callsign, ":"
                        WRITE(numout,*) "nbvmax = ", nbvmax
                        
                        ALLOC_ERR = -1
                        ALLOCATE(temp_var(nbvmax,lml), STAT=ALLOC_ERR)
                        IF (ALLOC_ERR /=0) THEN
                            WRITE(numout,*) "ERROR IN ALLOCATION OF temp_var :",ALLOC_ERR
                            STOP
                        ENDIF
                        ALLOC_ERR = -1
                        ALLOCATE(sub_index(nbpt,nbvmax,2), STAT=ALLOC_ERR)
                        IF (ALLOC_ERR /=0) THEN
                            WRITE(numout,*) "ERROR IN ALLOCATION OF sub_index :", ALLOC_ERR
                            STOP
                        ENDIF
                        sub_index(:,:,:) = zero
                        ALLOC_ERR = -1
                        ALLOCATE(sub_area(nbpt, nbvmax), STAT=ALLOC_ERR)
                        IF (ALLOC_ERR /=0) THEN
                            WRITE(numout,*) "ERROR IN ALLOCATION OF sub_area :", ALLOC_ERR
                            STOP
                        ENDIF
                        sub_area(:,:) = zero
                        resolution(:,:) = resolution(:,:)/1000 ! m -> km
                        write(*,*) "resolution updated: ", resolution(1,:), "km"
                        
                        CALL aggregate_p(nbpt, lalo, neighbours, resolution,contfrac, &
                        &                iml, jml, lon_rel, lat_rel, mask, callsign, &
                        &                nbvmax, sub_index, sub_area, ok_interpol)
                        
                        IF ( .NOT. ok_interpol ) THEN
                            DEALLOCATE(temp_var)
                            DEALLOCATE(sub_index)
                            DEALLOCATE(sub_area)
                            nbvmax = nbvmax * 2
                        ENDIF
                    ENDDO
                    
                    ! assign the values to plantdate
                    ! values should be given to all PFTs
                    DO ib = 1, nbpt
                        ! cropvar(ib,:,icyc) = zero
                        
                        ! examing all sub_point we found
                        fopt = COUNT(sub_area(ib,:)>zero)
                        
                        ! confirm that we found some points
                        IF ( fopt .EQ. 0) THEN
                            nbexp = nbexp + 1
                            cropvar(ib,ivgt,it) = val_exp ! cropvar_default(ivgt)
                        ELSE
                            DO ilf = 1,fopt
                                ! !Not to get lat and lon in wrong order
                                temp_var(ilf,ivgt) = cropvar_mat(sub_index(ib,ilf,1),sub_index(ib,ilf,2),ivgt,icyc,it)
                            ENDDO
                            
                            
                            ! sum_float = zero
                            ok_var = .FALSE.
                            maxcp = nbvmax/2;
                            DO WHILE (.NOT. ok_var)
                                
                                ALLOCATE(hashlst1(maxcp), STAT=ALLOC_ERR) ! no. 
                                IF (ALLOC_ERR /=0) THEN
                                    WRITE(numout,*) "ERROR IN ALLOCATION OF hashlst1:", ALLOC_ERR
                                    STOP
                                ENDIF
                                ALLOCATE(hashlst2(maxcp), STAT=ALLOC_ERR) ! area
                                IF (ALLOC_ERR /=0) THEN
                                    WRITE(numout,*) "ERROR IN ALLOCATION OF hashlst2:", ALLOC_ERR
                                    STOP
                                ENDIF
                                hashlst1(:) = 0
                                hashlst2(:) = zero
                                sgn = zero
                                cps = 1;
                                DO ilf = 1,fopt
                                    sgn = sgn + sub_area(ib,ilf)
                                    IF (cps==1) THEN
                                        hashlst1(cps) = temp_var(ilf,ivgt)
                                        hashlst2(cps) = hashlst2(cps) + sub_area(ib,ilf)
                                        cps = cps + 1
                                    ELSE
                                        IF (cps .EQ. maxcp) THEN
                                            EXIT
                                        ELSE
                                            DO jlf = 1,cps
                                                IF (temp_var(ilf,ivgt) .EQ. hashlst1(jlf)) THEN
                                                    hashlst2(jlf) = hashlst2(jlf) +sub_area(ib,ilf)
                                                    EXIT
                                                ELSEIF (jlf .EQ. cps) THEN 
                                                    ! .AND. (temp_var(ilf,ivgt) .NE.
                                                    ! hashlst(jlf,1))
                                                    hashlst1(cps) = temp_var(ilf,ivgt)
                                                    hashlst2(cps) = hashlst2(cps) + sub_area(ib,ilf) 
                                                    ! better to multiply PFT here
                                                    cps = cps + 1;
                                                ENDIF
                                            ENDDO
                                        ENDIF
                                    ENDIF
                                ENDDO
                                
                                WRITE(numout,*) "total variety number: ",cps
                                IF (cps .LT. maxcp) THEN
                                    ok_var = .TRUE.
                                    ! find the variety with maximum area
                                    
                                ELSE
                                    ! ok_var = .FALSE.
                                    maxcp = maxcp * 2
                                    DEALLOCATE(hashlst1)
                                    DEALLOCATE(hashlst2)
                                ENDIF
                            ENDDO
                            
                            ! Normalize the surface
                            IF ( sgn .LT. min_sechiba) THEN
                                nbexp = nbexp + 1
                                cropvar(ib,ivgt,it) = val_exp !plantdate_default(ivgt)
                            ELSE
                                areamax = zero
                                maxps = 0
                                IF (cps .LE. 1) THEN
                                    WRITE(numout,*) "no sub point found, probably an error occured"
                                    STOP
                                ENDIF
                                DO jlf = 1,cps-1
                                    IF (hashlst2(jlf) .GT. areamax) THEN
                                        areamax = hashlst2(jlf)
                                        maxps = jlf
                                    ENDIF
                                ENDDO
                                cropvar(ib,ivgt,it) = hashlst1(maxps)
                                DEALLOCATE(hashlst1)
                                DEALLOCATE(hashlst2)
                            !   ! INMAX
                            !   cropvar(ib,ivgt,icyc) = ANINT(sum_float/sgn)
                            ENDIF
                            
                        ENDIF
                    
                    ENDDO
                    
                    IF ( nbexp .GT. 0) THEN
                        WRITE(numout,*) 'slowproc_PlntDt : default plant date was applied in ', nbexp, 'grid(s)'
                        WRITE(numout,*) 'slowproc_PlntDt : These are either coastal points or having missing data'
                    ENDIF
                    DEALLOCATE (sub_area)
                    DEALLOCATE (sub_index)
                    DEALLOCATE (temp_var)
                    ! WRITE(numout,*) 'Planting Date of Site 1 veget ',ivgt,' :
                    ! ',plantdate(1,ivgt,icyc)
                ENDIF
            ENDDO
            ! End of Veget cycle    
        ENDDO 
        ! End of Rotation cycle
    ENDDO
    ! END of Time Axis cycle
    
    DEALLOCATE (lat_rel)
    DEALLOCATE (lon_rel)
    DEALLOCATE (mask)
    ! DEALLOCATE (sub_area)
    ! DEALLOCATE (sub_index)
    ! DEALLOCATE (temp_date)
    DEALLOCATE (cropvar_mat)
    
    WRITE (numout,*) 'Output Crop Variety:'
    WRITE (numout,*) 'timestep 1:'
    WRITE (numout,*) cropvar(1,:,1)
    IF (plantcyc>1) THEN
        WRITE (numout,*) 'timestep 2:'
        WRITE (numout,*) cropvar(1,:,2)
    ENDIF
    WRITE (numout,*) '***END of DEBUG INFO slowproc_CropVar***'
    ! WRITE (numout,*) 'Manual STOP after slowproc_CropVar'
    ! STOP
    RETURN
    
  END SUBROUTINE slowproc_CropVar
! End of Edition by Xuhui, Nov. 27th 010


!! ================================================================================================================================
!! SUBROUTINE 	: stomate_stics_ManageInput
!!
!>\BRIEF        Interpolate (extract) Planting Date information  for a specific crop typ
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S) : Xuhui Wang, Oct. 18th, 2010
!!
!! MAIN OUTPUT VARIABLE(S): None 
!!
!! REFERENCES	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE stomate_stics_ManageInput(nbpt, lalo, neighbours, resolution, contfrac,strIn, varname, manage, yrlen)
  
!    INTEGER, parameter :: i_std = 4
!    REAL, parameter :: r_std = 8
    ! 
    ! 0.1 INPUT
    ! 
    INTEGER(i_std), INTENT(in)  :: nbpt         ! Number of points for which the data needs to be interpolated (extracted)
    REAL(r_std), INTENT(in) :: lalo(nbpt,2)     ! Vector of latitude and longtitude
    INTEGER(i_std), INTENT(in)  :: neighbours(nbpt,8)   ! Vectors of neighbours for each grid point
    ! (1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW)
    REAL(r_std),INTENT(in)  :: resolution(nbpt,2)   ! The size in km of each grid box in lat and lon
    REAL(r_std),INTENT(in)  :: contfrac(nbpt)   ! The fraction of land in each grid box
    CHARACTER(LEN=30),INTENT(in) :: strIn       ! getin parameter and Call Sign of the management data
    CHARACTER(LEN=30),INTENT(in) :: varname     ! variable name in the nc file
    ! 
    ! 0.2 OUTPUT
    ! 
    REAL(r_std),ALLOCATABLE, INTENT(out)    :: manage(:,:,:)    ! The planting date of the crop: nbpt, veg, year
    ! nvm is the number of PFTs, there may not be planting date for all the PFTs
    INTEGER(i_std), INTENT(out)             :: yrlen            ! year length of the output matrix
    ! 
    ! 0.3 LOCAL
    ! 
    INTEGER(i_std)      :: nbvmax       ! a parameter for interpolation
    REAL(r_std)         :: myres(nbpt,2)
    CHARACTER(LEN=80)       :: filename
    INTEGER(i_std)      :: iml, jml, lml, tml, fid, fid1
    INTEGER(i_std)      :: ip, jp, ib, ilf, fopt, it ! for-loop variable
    INTEGER(i_std)      :: nbexp
    REAL(r_std)         :: lev(1), date, dt
    REAL(r_std)         :: missing_val
    INTEGER(i_std)      :: itau(1)

    INTEGER(i_std)      :: nb_dim
    INTEGER,DIMENSION(flio_max_var_dims) :: l_d_w
    LOGICAL         :: l_ex
    
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)    :: lat_rel, lon_rel
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:,:)    :: manage_mat ! LON LAT VEGET, Time
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:) :: mask
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)    :: temp_data
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)    :: sub_area
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:,:)   :: sub_index

    REAL(r_std) :: sgn, sum_float
    INTEGER(i_std) :: ivgt   ! , icyc, pltcyc
    CHARACTER(LEN=30) :: callsign
    LOGICAL :: ok_interpol
    INTEGER :: ALLOC_ERR
    LOGICAL :: mydebug = .false.

!   ! croptype = TRIM(croptype) !if croptype is a string
!   ! else a switch expression is needed
!   filename = "/work/cont003/p529tan/WXH/plt_date_modif.nc" ! default input
!   file
!   ! String operation needed
    filename = "PlantingDate.nc"
    CALL getin_p(strIn,filename)

    IF (is_root_prc) THEN
    ! ? what does is_root_prc mean?
        CALL flininfo(filename, iml, jml, lml, tml, fid)
    ENDIF
    CALL bcast(iml)
    CALL bcast(jml)
    ! CALL bcast(lml)
    ! CALL bcast(tml)
    
    ! Printing information for debugging
    IF (mydebug) THEN
        WRITE(numout, *) "Xuhui's debug info for stomate_stics_ManageInput #1:"
        WRITE(numout, *) "string in: ", strIn
        WRITE(numout, *) "variable name: ", varname
        WRITE(numout, *) "filename is: ", filename
        WRITE(numout, *) "Dimension 1, lon, iml:", iml
        WRITE(numout, *) "Dimension 2, lat, jml:", jml
        WRITE(numout, *) "Dimension 3, veget, lml:", lml
        WRITE(numout, *) "Dimension 4, time, tml:", tml
    ENDIF
    ! apparently, flinget function is not designed to take veget but levels to
    ! be the 
    ! 3rd dimension, modification to lml is needed

!JG all flio calls must be done by is_root_prc
    IF (is_root_prc) THEN
       CALL flioopfd(filename,fid1)
       CALL flioinqv(fid1,v_n=varname, l_ex = l_ex, nb_dims = nb_dim, len_dims =l_d_w)
       IF (lml == 0) THEN
          ! CALL
          ! flioinqv(fid1,v_n="PLNTDT",l_ex=l_ex,nb_dims=nb_dim,len_dims=l_d_w)
          lml=l_d_w(3)
          IF (mydebug) THEN
              WRITE(numout, *) "len_dims: ", l_d_w
              WRITE(numout, *) "lml AFTER revision"
              WRITE(numout, *) "lml: ", lml
          ENDIF
       ENDIF
       IF (mydebug) THEN
           WRITE(numout,*) "nb_dim: ", nb_dim
           WRITE(numout,*) "resolution: ", resolution(1,:)
       ENDIF
       
       IF (nb_dim .NE. 4) THEN
          WRITE(numout,*) "dimension not supported for ", nb_dim
       ENDIF
       tml = l_d_w(4)
       !yrlen = tml
    END IF
    IF (mydebug) THEN
        WRITE(numout, *) "Now the tml is, :", tml 
        WRITE(numout, *) "Now the lml is:", lml
    ENDIF
     
!JG REMVOVE    CALL flioclo(fid1)
    CALL bcast(lml)
    CALL bcast(tml)
    CALL bcast(nb_dim)
    ! CALL bcast(plantcyc)
    
    ! JG yrlen must not be done after bcast(tml)
    yrlen = tml
    
    ALLOC_ERR=-1
    ALLOCATE(manage(nbpt,lml,tml),STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
        WRITE(numout,*) "ERROR IN ALLOCATION OF manage: ", ALLOC_ERR
    ENDIF
    WRITE(numout,*) "manage ALLOCATED"
    !CALL bcast(manage)
    
    ! 
    ALLOC_ERR=-1
    ALLOCATE(lat_rel(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR/=0) THEN
        WRITE(numout,*) "ERROR IN ALLOCATION of lat_rel : ",ALLOC_ERR
        STOP 
    ENDIF
!    CALL bcast(lat_rel)

    ALLOC_ERR=-1
    ALLOCATE(lon_rel(iml,jml), STAT=ALLOC_ERR)
   IF (ALLOC_ERR/=0) THEN
        WRITE(numout,*) "ERROR IN ALLOCATION of lon_rel : ",ALLOC_ERR
        STOP 
    ENDIF
!    CALL bcast(lon_rel)

    ALLOC_ERR=-1
    ALLOCATE(mask(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
        WRITE(numout,*) "ERROR IN ALLOCATION of mask : ",ALLOC_ERR
        STOP 
    ENDIF
!    CALL bcast(mask)
    

    ALLOC_ERR=-1
    ALLOCATE(manage_mat(iml,jml,lml,tml), STAT=ALLOC_ERR) 
    ! !lml is supposed to be nvm (number of PFTs), if not ,change it
    IF (ALLOC_ERR/=0) THEN
        WRITE(numout,*) "ERROR IN ALLOCATION of manage_mat : ",ALLOC_ERR
        STOP 
    ENDIF
!    CALL bcast(manage_mat)
!    WRITE (numout,*) 'bcast manage_mat'

    ! input of some attributes
    IF (is_root_prc) THEN
! JG with the flioclo, done before this was not ok. Now ok
        CALL flinget(fid, 'LON', iml, jml, lml, tml, 1, 1, lon_rel)
        CALL flinget(fid, 'LAT', iml, jml, lml, tml, 1, 1, lat_rel)
    ENDIF
    CALL bcast(lon_rel)
    CALL bcast(lat_rel)
    WRITE (numout,*) 'lon_rel size: ', SIZE(lon_rel)
    WRITE (numout,*) 'lat_rel size: ', SIZE(lat_rel)
    

    ! input of the matrix
    IF (is_root_prc) THEN 
        ! CALL flinget(fid, 'PLNTDT', iml, jml, lml, tml, 1, 1, plntdt_mat) 
! JG remove CALL flioopfd: already done 
!       CALL flioopfd(filename,fid1)
        CALL fliogetv(fid1,trim(varname),manage_mat,start=(/1,1,1,1/),count=(/iml,jml,lml,tml/))
        ! get missing_val
        CALL fliogeta(fid1,varname,'missing_value',missing_val)
        CALL flioclo(fid1)
    ENDIF
    CALL bcast(manage_mat)
    CALL bcast(missing_val)
    WRITE (numout,*) 'bcast manage_mat'
    

    ! WRITE(numout,*) 'manage_mat size: ',SIZE(manage_mat)
    ! WRITE(numout,*) 'missing value: ', missing_val
    ! WRITE(numout,*) 'lat(361,284): ',lat_rel(361,284)
    ! WRITE(numout,*) 'lon(361,284): ',lon_rel(361,284)
    ! WRITE(numout,*) 'plntdt(361,284,1,1): ',plntdt_mat(361,284,1,1)
    
    IF (is_root_prc) CALL flinclo(fid)
    
    manage(:,:,:) = zero ! nbpt veget year
    
    DO it = 1,tml
        DO ivgt = 1,lml ! ? We can suppose PFTs less than 10 are natural veg without planting date, but not now
            IF (.NOT. natural(ivgt)) THEN
                WRITE(numout,*) "veget, time: ", ivgt,it
                nbexp = 0
                ! the number of exceptions
                
                ! mask of available value
                mask(:,:) = zero;  ! Defined in constante.f90
                DO ip = 1,iml
                    DO jp = 1,jml
                        IF ((manage_mat(ip,jp,ivgt,it) .GT. min_sechiba) .AND. &
                        (manage_mat(ip,jp,ivgt,it) /= missing_val)) THEN
                            mask(ip,jp) = un;  ! Defined in constante.f90
                            ! here we assumed that for each plant cycle at each 
                            ! there might be missing data at different grid
                            ! in this case, mask has to be generated each plant
                            ! cycle each time step
                        ENDIF
                    ENDDO
                ENDDO
                
                ! Interpolation started
                nbvmax = 200
                ! the maximum amount of fine grids that one coarse grid may have
                
                callsign = strIn
                
                ok_interpol = .FALSE.
                
                DO WHILE ( .NOT. ok_interpol )
                    WRITE(numout,*) "Pojection arrays for ", callsign, ":"
                    WRITE(numout,*) "nbvmax = ", nbvmax
                    
                    ALLOC_ERR = -1
                    ALLOCATE(temp_data(nbvmax,lml), STAT=ALLOC_ERR)
                    IF (ALLOC_ERR /=0) THEN
                        WRITE(numout,*) "ERROR IN ALLOCATION OF temp_data :", ALLOC_ERR
                        STOP
                    ENDIF
                    temp_data = zero
                    ALLOC_ERR = -1
                    ALLOCATE(sub_index(nbpt,nbvmax,2), STAT=ALLOC_ERR)
                    IF (ALLOC_ERR /=0) THEN
                        WRITE(numout,*) "ERROR IN ALLOCATION OF sub_index :", ALLOC_ERR
                        STOP
                    ENDIF
                    sub_index(:,:,:) = zero
                    ALLOC_ERR = -1
                    ALLOCATE(sub_area(nbpt, nbvmax), STAT=ALLOC_ERR)
                    IF (ALLOC_ERR /=0) THEN
                        WRITE(numout,*) "ERROR IN ALLOCATION OF sub_area :",ALLOC_ERR
                        STOP
                    ENDIF
                    sub_area(:,:) = zero
                    myres(:,:) = resolution(:,:)/1000  !m -> km
                    write(numout,*) "resolution updated: ", myres(1,:), " km"
                    !CALL bcast(myres)
!                    CALL bcast(myres)
                    
!                    write(*,*) "calling aggregate_p? "
                   CALL aggregate_p(nbpt, lalo, neighbours, myres, contfrac, &
                    &                iml, jml, lon_rel, lat_rel, mask, callsign, &
                    &                nbvmax, sub_index, sub_area, ok_interpol)
!                    write(numout,*) "wu: we finished aggregate_p:) "                   
 
                    IF ( .NOT. ok_interpol ) THEN
                        DEALLOCATE(temp_data)
                        DEALLOCATE(sub_index)
                        DEALLOCATE(sub_area)
                        nbvmax = nbvmax * 2
                    ENDIF
                ENDDO
                
!                WRITE(numout,*) "called aggregate_p"
                ! assign the values to plantdate
                ! values should be given to all PFTs
                DO ib = 1, nbpt
                    ! examing all sub_point we found
                    fopt = COUNT(sub_area(ib,:)>zero)
                    
                    ! confirm that we found some points
                    IF ( fopt .EQ. 0) THEN
                        nbexp = nbexp + 1
                        manage(ib,ivgt,it) = val_exp
                    ELSE
                        DO ilf = 1,fopt
                            ! !Not to get lat and lon in wrong order
                            temp_data(ilf,ivgt) = manage_mat(sub_index(ib,ilf,1),sub_index(ib,ilf,2),ivgt,it)
                        ENDDO
                        
                        sgn = zero
                        sum_float = zero
                        DO ilf = 1,fopt
                            ! average the data weighted by area ! better to multiply
                            ! PFT HERE
                            ! need to add management specific judgem
                                sum_float = sum_float + temp_data(ilf,ivgt)*sub_area(ib,ilf)
                                sgn = sgn + sub_area(ib,ilf)
                        ENDDO
                        
                        ! Normalize the surface
                        ! sgn can be a scaler, however, to prepare it for future
                        ! incorporation of fraction
                        ! I make it a vector with nvm values which are equal to each
                        ! other
                        IF ( sgn .LT. min_sechiba) THEN
                            nbexp = nbexp + 1
                            manage(ib,ivgt,it) = val_exp ! plantdate_default(ivgt)
                        ELSE
                            manage(ib,ivgt,it) = ANINT(sum_float/sgn)
                        ENDIF
                        
                    ENDIF
                
                ENDDO ! ib
                
                IF ( nbexp .GT. 0) THEN
                    WRITE(numout,*) 'stomate_stics_ManageInput : exp_val was applied in', nbexp, 'grid(s)'
                    WRITE(numout,*) 'stomate_stics_ManageInput : These are either coastal points or having missing data'
                ENDIF
                DEALLOCATE (sub_area)
                DEALLOCATE (sub_index)
                DEALLOCATE (temp_data)
                ! WRITE(numout,*) 'Planting Date of Site 1 veget ',ivgt,' :
                ! ',plantdate(1,ivgt,icyc)
            ENDIF
        ENDDO
        ! End of Veget cycle    
    ENDDO
    ! End of Time Axis cycle
    
    DEALLOCATE (lat_rel)
    DEALLOCATE (lon_rel)
    DEALLOCATE (mask)
    DEALLOCATE (manage_mat)
    
    WRITE (numout,*) 'Output Management Date:'
    WRITE (numout,*) 'time_step 1:'
    WRITE (numout,*) manage(1,:,1)
    IF (tml>1) THEN
        WRITE (numout,*) 'time_step 2:'
        WRITE (numout,*) manage(1,:,2)
    ENDIF
    WRITE (numout,*) '***END of DEBUG INFO stomate_stics_ManageInput***'
    RETURN
    
  END SUBROUTINE stomate_stics_ManageInput
! End of Edition by Xuhui, Mar. 16th 2011



! Date: 22.06.2014, Xuhui Wang
! Interpolate (extract) Nitrogen fertilization information 
! for a specific crop type
!! ================================================================================================================================
!! SUBROUTINE 	: stomate_stics_NfertInput 
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
  SUBROUTINE stomate_stics_NfertInput(nbpt, lalo, neighbours, resolution, contfrac,strIn, varname, Nfert, yrlen)
  
!    INTEGER, parameter :: i_std = 4
!    REAL, parameter :: r_std = 8
    ! 
    ! 0.1 INPUT
    ! 
    INTEGER(i_std), INTENT(in)  :: nbpt         ! Number of points for which the data needs to be interpolated (extracted)
    REAL(r_std), INTENT(in) :: lalo(nbpt,2)     ! Vector of latitude and longtitude
    INTEGER(i_std), INTENT(in)  :: neighbours(nbpt,8)   ! Vectors of neighbours for each grid point
    ! (1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW)
    REAL(r_std),INTENT(in)  :: resolution(nbpt,2)   ! The size in km of each grid box in lat and lon
    !REAL(r_std)             :: resolution(nbpt,2)   ! The size in km of each grid box in lat and lon
    REAL(r_std),INTENT(in)  :: contfrac(nbpt)   ! The fraction of land in each grid box
    CHARACTER(LEN=30),INTENT(in) :: strIn       ! getin parameter and Call Sign of the Nfertment data
    CHARACTER(LEN=30),INTENT(in) :: varname     ! variable name in the nc file
    ! 
    ! 0.2 OUTPUT
    ! 
    REAL(r_std),ALLOCATABLE, INTENT(out)    :: Nfert(:,:,:)    ! The planting date of the crop: nbpt, veg, year
    ! nvm is the number of PFTs, there may not be planting date for all the PFTs
    INTEGER(i_std), INTENT(out)             :: yrlen            ! year length of the output matrix
    ! 
    ! 0.3 LOCAL
    ! 
    INTEGER(i_std)      :: nbvmax       ! a parameter for interpolation
    REAL(r_std)         :: myres(nbpt,2)
    CHARACTER(LEN=80)       :: filename
    INTEGER(i_std)      :: iml, jml, lml, tml, fid, fid1
    INTEGER(i_std)      :: ip, jp, ib, ilf, fopt, it ! for-loop variable
    INTEGER(i_std)      :: nbexp
    REAL(r_std)         :: lev(1), date, dt
    REAL(r_std)         :: missing_val
    INTEGER(i_std)      :: itau(1)

    INTEGER(i_std)      :: nb_dim
    INTEGER,DIMENSION(flio_max_var_dims) :: l_d_w
    LOGICAL         :: l_ex
    
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)    :: lat_rel, lon_rel
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:,:)    :: Nfert_mat ! LON LAT VEGET, Time
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:) :: mask
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)    :: temp_data
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)    :: sub_area
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:,:)   :: sub_index

    REAL(r_std) :: sgn, sum_float
    INTEGER(i_std) :: ivgt   ! , icyc, pltcyc
    CHARACTER(LEN=30) :: callsign
    LOGICAL :: ok_interpol
    INTEGER :: ALLOC_ERR

!   ! croptype = TRIM(croptype) !if croptype is a string
!   ! else a switch expression is needed
!   filename = "/work/cont003/p529tan/WXH/plt_date_modif.nc" ! default input
!   file
!   ! String operation needed
    filename = "Nfert.nc"
    CALL getin_p(strIn,filename)

    IF (is_root_prc) THEN
    ! ? what does is_root_prc mean?
        CALL flininfo(filename, iml, jml, lml, tml, fid)
    ENDIF
    CALL bcast(iml)
    CALL bcast(jml)
    ! CALL bcast(lml)
    ! CALL bcast(tml)
    
    ! Printing information for debugging
    WRITE(numout, *) "Xuhui's debug info for stomate_stics_NfertInput #1:"
    WRITE(numout, *) "string in: ", strIn
    WRITE(numout, *) "variable name: ", varname
    WRITE(numout, *) "filename is: ", filename
    WRITE(numout, *) "Dimension 1, lon, iml:", iml
    WRITE(numout, *) "Dimension 2, lat, jml:", jml
    WRITE(numout, *) "Dimension 3, veget, lml:", lml
    WRITE(numout, *) "Dimension 4, time, tml:", tml
    ! apparently, flinget function is not designed to take veget but levels to
    ! be the 
    ! 3rd dimension, modification to lml is needed

!JG all flio calls must be done by is_root_prc
    IF (is_root_prc) THEN
       CALL flioopfd(filename,fid1)
       CALL flioinqv(fid1,v_n=varname, l_ex = l_ex, nb_dims = nb_dim, len_dims =l_d_w)
       IF (lml == 0) THEN
          ! CALL
          ! flioinqv(fid1,v_n="PLNTDT",l_ex=l_ex,nb_dims=nb_dim,len_dims=l_d_w)
          lml=l_d_w(3)
          WRITE(numout, *) "len_dims: ", l_d_w
          WRITE(numout, *) "lml AFTER revision"
          WRITE(numout, *) "lml: ", lml
       ENDIF
       WRITE(numout,*) "nb_dim: ", nb_dim
       WRITE(numout,*) "resolution: ", resolution(1,:)
       
       IF (nb_dim .NE. 4) THEN
          WRITE(numout,*) "dimension not supported for ", nb_dim
       ENDIF
       tml = l_d_w(4)
       !yrlen = tml
    END IF
    
    WRITE(numout, *) "Now the tml is, :", tml 
    WRITE(numout, *) "Now the lml is:", lml
     
!JG REMVOVE    CALL flioclo(fid1)
    CALL bcast(lml)
    CALL bcast(tml)
    CALL bcast(nb_dim)
    ! CALL bcast(plantcyc)
    
    ! JG yrlen must not be done after bcast(tml)
    yrlen = tml
    
    ALLOC_ERR=-1
    ALLOCATE(Nfert(nbpt,lml,tml),STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
        WRITE(numout,*) "ERROR IN ALLOCATION OF Nfert: ", ALLOC_ERR
    ENDIF
    WRITE(numout,*) "Nfert ALLOCATED"
    !CALL bcast(Nfert)
    
    ! 
    ALLOC_ERR=-1
    ALLOCATE(lat_rel(iml,jml), STAT=ALLOC_ERR)
      IF (ALLOC_ERR/=0) THEN
        WRITE(numout,*) "ERROR IN ALLOCATION of lat_rel : ",ALLOC_ERR
        STOP 
    ENDIF
!    CALL bcast(lat_rel)

    ALLOC_ERR=-1
    ALLOCATE(lon_rel(iml,jml), STAT=ALLOC_ERR)
   IF (ALLOC_ERR/=0) THEN
        WRITE(numout,*) "ERROR IN ALLOCATION of lon_rel : ",ALLOC_ERR
        STOP 
    ENDIF
!    CALL bcast(lon_rel)

    ALLOC_ERR=-1
    ALLOCATE(mask(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
        WRITE(numout,*) "ERROR IN ALLOCATION of mask : ",ALLOC_ERR
        STOP 
    ENDIF
!    CALL bcast(mask)
    

    ALLOC_ERR=-1
    ALLOCATE(Nfert_mat(iml,jml,lml,tml), STAT=ALLOC_ERR) 
    ! !lml is supposed to be nvm (number of PFTs), if not ,change it
    IF (ALLOC_ERR/=0) THEN
        WRITE(numout,*) "ERROR IN ALLOCATION of Nfert_mat : ",ALLOC_ERR
        STOP 
    ENDIF
!    CALL bcast(Nfert_mat)
!    WRITE (numout,*) 'bcast Nfert_mat'

    ! input of some attributes
    IF (is_root_prc) THEN
! JG with the flioclo, done before this was not ok. Now ok
        CALL flinget(fid, 'LON', iml, jml, lml, tml, 1, 1, lon_rel)
        CALL flinget(fid, 'LAT', iml, jml, lml, tml, 1, 1, lat_rel)
    ENDIF
    CALL bcast(lon_rel)
    CALL bcast(lat_rel)
    WRITE (numout,*) 'lon_rel size: ', SIZE(lon_rel)
    WRITE (numout,*) 'lat_rel size: ', SIZE(lat_rel)
    

    ! input of the matrix
    IF (is_root_prc) THEN 
        ! CALL flinget(fid, 'PLNTDT', iml, jml, lml, tml, 1, 1, plntdt_mat) 
! JG remove CALL flioopfd: already done 
!       CALL flioopfd(filename,fid1)
        CALL fliogetv(fid1,trim(varname),Nfert_mat,start=(/1,1,1,1/),count=(/iml,jml,lml,tml/))
        ! get missing_val
        CALL fliogeta(fid1,varname,'missing_value',missing_val)
        CALL flioclo(fid1)
    ENDIF
    CALL bcast(Nfert_mat)
    
    IF (is_root_prc) CALL flinclo(fid)
    
    Nfert(:,:,:) = zero ! nbpt veget year
    
    DO it = 1,tml
        DO ivgt = 1,lml ! ? We can suppose PFTs less than 10 are natural veg without planting date, but not now
            IF (.NOT. natural(ivgt)) THEN
                WRITE(numout,*) "veget, time: ", ivgt,it
                nbexp = 0
                ! the number of exceptions
                
                ! mask of available value
                mask(:,:) = zero;  ! Defined in constante.f90
                DO ip = 1,iml
                    DO jp = 1,jml
                        IF ((Nfert_mat(ip,jp,ivgt,it) .GE. zero) .AND. &
                        (Nfert_mat(ip,jp,ivgt,it) /= missing_val)) THEN
                            mask(ip,jp) = un;  ! Defined in constante.f90
                            ! here we assumed that for each plant cycle at each 
                            ! there might be missing data at different grid
                            ! in this case, mask has to be generated each plant
                            ! cycle each time step
                        ENDIF
                    ENDDO
                ENDDO
                
                ! Interpolation started
                nbvmax = 200
                ! the maximum amount of fine grids that one coarse grid may have
                
                callsign = strIn
                
                ok_interpol = .FALSE.
                
                DO WHILE ( .NOT. ok_interpol )
                    WRITE(numout,*) "Pojection arrays for ", callsign, ":"
                    WRITE(numout,*) "nbvmax = ", nbvmax
                    
                    ALLOC_ERR = -1
                    ALLOCATE(temp_data(nbvmax,lml), STAT=ALLOC_ERR)
                    IF (ALLOC_ERR /=0) THEN
                        WRITE(numout,*) "ERROR IN ALLOCATION OF temp_data :", ALLOC_ERR
                        STOP
                    ENDIF
                    ALLOC_ERR = -1
                    ALLOCATE(sub_index(nbpt,nbvmax,2), STAT=ALLOC_ERR)
                    IF (ALLOC_ERR /=0) THEN
                        WRITE(numout,*) "ERROR IN ALLOCATION OF sub_index :", ALLOC_ERR
                        STOP
                    ENDIF
                    sub_index(:,:,:) = zero
                    ALLOC_ERR = -1
                    ALLOCATE(sub_area(nbpt, nbvmax), STAT=ALLOC_ERR)
                    IF (ALLOC_ERR /=0) THEN
                        WRITE(numout,*) "ERROR IN ALLOCATION OF sub_area :",ALLOC_ERR
                        STOP
                    ENDIF
                    sub_area(:,:) = zero
                    myres(:,:) = resolution(:,:)/1000  !m -> km
                    write(numout,*) "resolution updated: ", myres(1,:), " km"
                    !CALL bcast(myres)
                    
!                    write(numout,*) "wu: stepping in aggregate_p? "
                    CALL aggregate_p(nbpt, lalo, neighbours, myres, contfrac, &
                    &                iml, jml, lon_rel, lat_rel, mask, callsign, &
                    &                nbvmax, sub_index, sub_area, ok_interpol)
!                    write(numout,*) "wu: we finished aggregate_p:) "                   
 
                    IF ( .NOT. ok_interpol ) THEN
                        DEALLOCATE(temp_data)
                        DEALLOCATE(sub_index)
                        DEALLOCATE(sub_area)
                        nbvmax = nbvmax * 2
                    ENDIF
                ENDDO
                
!                WRITE(numout,*) "called aggregate_p"
                ! assign the values to plantdate
                ! values should be given to all PFTs
                DO ib = 1, nbpt
                    ! examing all sub_point we found
                    fopt = COUNT(sub_area(ib,:)>zero)
                    
                    ! confirm that we found some points
                    IF ( fopt .EQ. 0) THEN
                        nbexp = nbexp + 1
!                        Nfert(ib,ivgt,it) = val_exp
                        Nfert(ib,ivgt,it) = 0
                    ELSE
                        DO ilf = 1,fopt
                            ! !Not to get lat and lon in wrong order
                            temp_data(ilf,ivgt) = Nfert_mat(sub_index(ib,ilf,1),sub_index(ib,ilf,2),ivgt,it)
                        ENDDO
                        
                        sgn = zero
                        sum_float = zero
                        DO ilf = 1,fopt
                            ! average the data weighted by area ! better to multiply
                            ! PFT HERE
                            ! need to add Nfertment specific judgem
                                sum_float = sum_float + temp_data(ilf,ivgt)*sub_area(ib,ilf)
                                sgn = sgn + sub_area(ib,ilf)
                        ENDDO
                        
                        ! Normalize the surface
                        ! sgn can be a scaler, however, to prepare it for future
                        ! incorporation of fraction
                        ! I make it a vector with nvm values which are equal to each
                        ! other
                        IF ( sgn .LT. min_sechiba) THEN
                            nbexp = nbexp + 1
!                            Nfert(ib,ivgt,it) = val_exp ! plantdate_default(ivgt)
                            Nfert(ib,ivgt,it) = 0 ! plantdate_default(ivgt)
                        ELSE
                            !Nfert(ib,ivgt,it) = ANINT(sum_float/sgn)
                            Nfert(ib,ivgt,it) = sum_float/sgn
                        ENDIF
                        
                    ENDIF
                
                ENDDO ! ib
                
                IF ( nbexp .GT. 0) THEN
                    WRITE(numout,*) 'stomate_stics_NfertInput : 0 was applied in', nbexp, 'grid(s)'
                    WRITE(numout,*) 'stomate_stics_NfertInput : These are either coastal points or having missing data'
                ENDIF
                DEALLOCATE (sub_area)
                DEALLOCATE (sub_index)
                DEALLOCATE (temp_data)
                ! WRITE(numout,*) 'Planting Date of Site 1 veget ',ivgt,' :
                ! ',plantdate(1,ivgt,icyc)
            ENDIF
        ENDDO
        ! End of Veget cycle    
    ENDDO
    ! End of Time Axis cycle
    
    DEALLOCATE (lat_rel)
    DEALLOCATE (lon_rel)
    DEALLOCATE (mask)
    DEALLOCATE (Nfert_mat)
    
  END SUBROUTINE stomate_stics_NfertInput
! End of Edition by Xuhui, Mar. 16th 2011
 
!! ================================================================================================================================
!! SUBROUTINE 	: stomate_stics_calc_N_limfert
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
  SUBROUTINE stomat_stics_calc_N_limfert(npts,N_fert_total,N_limfert)
 
  INTEGER (i_std)                             , INTENT(in)  :: npts
  !REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(in) :: nfertamm
  !REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(in) :: nfertnit
  !REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(in) :: nliquidmanure
  !REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(in) :: nslurry
  !REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(in) :: nsolidmanure
  !REAL(r_std), DIMENSION(npts,nvm)          , INTENT(in) :: legume_fraction
  !REAL(r_std), DIMENSION(npts,nvm)          , INTENT(in) :: soil_fertility
  REAL(r_std), DIMENSION(npts,nvm)          , INTENT(in)  :: N_fert_total
  REAL(r_std), DIMENSION(:,:)               , INTENT(out) :: N_limfert
 
 
  INTEGER(i_std) :: k,j,ip
 
    !N_fert_total(:,:) = 0.0
    !DO k=1,nstocking
    !  N_fert_total(:,:) = (N_fert_total(:,:) + nfertamm(:,:,k) + &
    !                      nfertnit(:,:,k) + nliquidmanure(:,:,k) + &
    !                      nslurry(:,:,k) + nsolidmanure(:,:,k))
    !ENDDO
    !N_fert_total(:,:) = N_fert_total(:,:) * 10000
    !DO j=2,nvm
    !  IF ((management_intensity(j) .EQ. 2).AND. is_c3(j)) THEN
    !    N_fert_total(:,mcut_C3)=N_fert_total(:,j)
    !    N_fert_total(:,mgraze_C3)=N_fert_total(:,j)
    !  ENDIF
    !  IF ((management_intensity(j) .EQ. 2).AND. (.NOT.is_c3(j))) THEN
    !    N_fert_total(:,mcut_C4)=N_fert_total(:,j)
    !    N_fert_total(:,mgraze_C4)=N_fert_total(:,j)
    !  ENDIF
    !  
    !ENDDO
   
 
!    N_limfert(:,:) = 1. + N_effect - N_effect * (0.75 ** (N_fert_total(:,:)/30))
    N_limfert(:,:) = 0
    DO j=2, nvm
        IF ( (SP_neffmax(j) .lt. zero) .AND. ( ok_LAIdev(j) ) ) THEN
            CALL ipslerr_p(3, 'calc_limfert', 'negative N effect', '', '')
        ENDIF
        IF ( ok_LAIdev(j) ) THEN
            N_limfert(:,j) = 1. + SP_neffmax(j) - SP_neffmax(j)*(SP_nsatrat(j) ** (N_fert_total(:,j)/10))
        ELSE
            N_limfert(:,j) = 1.
        ENDIF
    ENDDO

!    WHERE (N_limfert(:,:) .LT. 1.0) ! never reached if N_effect>0
!      N_limfert(:,:) = 1.0
!    ELSEWHERE (N_limfert(:,:) .GT. 2.0) ! will yield mistake if N_effect>1
!      N_limfert(:,:) = 1.+N_effect
!    ENDWHERE
    
  END SUBROUTINE stomat_stics_calc_N_limfert
  ! end of the nitrogen fertilization module


!! ================================================================================================================================
!! SUBROUTINE 	: stomate_stics_read_cycle 
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
  SUBROUTINE stomate_stics_read_cycle(nbpt,lalo,neighbours,resolution,contfrac, & ! in
                        cyc_num, cyc_num_tot ) ! out
    ! 0.1 INPUT
    !    
    INTEGER(i_std), INTENT(in)  :: nbpt         ! Number of points for which the data needs to be interpolated (extracted)
    REAL(r_std), INTENT(in) :: lalo(nbpt,2)     ! Vector of latitude and longtitude
    INTEGER(i_std), INTENT(in)  :: neighbours(nbpt,8)   ! Vectors of neighbours for each grid point
    ! (1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW)
    REAL(r_std),INTENT(in)  :: resolution(nbpt,2)   ! The size in km of each grid box in lat and lon
    !REAL(r_std)             :: resolution(nbpt,2)   ! The size in km of each     !grid box in lat and lon
    REAL(r_std),INTENT(in)  :: contfrac(nbpt)   ! The fraction of land in each grid box

    INTEGER(i_std),INTENT(out), DIMENSION(nbpt,nvm) :: cyc_num ! flag for starting the rotation for kjpindex, nvm
    INTEGER(i_std),INTENT(out), DIMENSION(nbpt) :: cyc_num_tot ! number of current rotation cycle

    ! 0.3 local variables
    INTEGER(i_std)                                        :: yrlen
    CHARACTER(LEN=30)                                     :: strManage
    CHARACTER(LEN=30)                                     :: strVar
    REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE      :: manage
    INTEGER(i_std)                                :: ip

!!!!-------------------------------------------------------------------------------------------

            cyc_num(:,:) = 1 

            CALL getin_p('IMPOSE_ROT',rot_1d)
            ! by default, rot_1d is true
            IF (nbpt==1) THEN
                cyc_num_tot(1) = cyc_rot_max
            ELSE    ! multiple point simulations
                IF (rot_1d) THEN
                    cyc_num_tot(:) = cyc_rot_max
                ELSE
                    strManage = 'NUMROTATE_FILE'
                    strVar = 'CycleNum'
                    CALL stomate_stics_ManageInput(nbpt,lalo,neighbours,resolution,contfrac,strManage,strVar,manage,yrlen)
                    DO ip = 1,nbpt
                        cyc_num_tot(ip) = INT(MAXVAL(manage(ip,:,1)),i_std)
                        IF (cyc_num_tot(ip) .LT. 1) THEN
                            cyc_num_tot(ip) = 1 ! at least one cycle in rotation
                            WRITE(numout,*) 'WARNING xuhui: ip, cyc_num_tot(ip) ', ip, MAXVAL(manage(ip,:,1))
                        ENDIF
                    ENDDO
                    IF (MAXVAL(cyc_num_tot(:))>10) THEN
                        WRITE(numout,*) 'xuhui: rotation cycles longer than 10, likely a bug'
                    ENDIF
                    IF (MAXVAL(cyc_num_tot(:))>cyc_rot_max) THEN
                        WRITE(numout,*) 'xuhui: No. of rotation updated'
                        WRITE(numout,*) 'MAXVAL(cyc_num_tot(:)), cyc_rot_max', MAXVAL(cyc_num_tot(:)), cyc_rot_max
                        cyc_rot_max = MAXVAL(cyc_num_tot(:))
                    ELSEIF (MAXVAL(cyc_num_tot(:))<cyc_rot_max) THEN
                        WRITE(numout,*) 'xuhui: MAXVAL(cyc_num_tot(:))<cyc_rot_max',MAXVAL(cyc_num_tot(:)), cyc_rot_max
                        WRITE(numout,*) 'bad setting in cyc_rot_max'
                        WRITE(numout,*) 'xuhui: No. of rotation updated'
                        cyc_rot_max = MAXVAL(cyc_num_tot(:))
                    ENDIF
                ENDIF
            ENDIF
            IF (ALLOCATED(manage)) DEALLOCATE(manage) ! clear manage for other input purpose

  END SUBROUTINE stomate_stics_read_cycle

!! ================================================================================================================================
!! SUBROUTINE 	: stomate_stics_read_rotation 
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
  SUBROUTINE stomate_stics_read_rotation(nbpt,lalo,neighbours,resolution,contfrac, & ! in
                                   rot_cmd_store )  ! out
    ! 0.1 INPUT
    !    
    INTEGER(i_std), INTENT(in)  :: nbpt         ! Number of points for which the data needs to be interpolated (extracted)
    REAL(r_std), INTENT(in) :: lalo(nbpt,2)     ! Vector of latitude and longtitude
    INTEGER(i_std), INTENT(in)  :: neighbours(nbpt,8)   ! Vectors of neighbours for each grid point
    ! (1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW)
    REAL(r_std),INTENT(in)  :: resolution(nbpt,2)   ! The size in km of each grid box in lat and lon
    !REAL(r_std)             :: resolution(nbpt,2)   ! The size in km of each     !grid box in lat and lon
    REAL(r_std),INTENT(in)  :: contfrac(nbpt)   ! The fraction of land in each grid box

    INTEGER(i_std), INTENT(out), ALLOCATABLE, DIMENSION(:,:,:) :: rot_cmd_store

    ! 0.3 local variables
    INTEGER(i_std)                                        :: yrlen
    CHARACTER(LEN=30)                                     :: strManage
    CHARACTER(LEN=30)                                     :: strVar
    CHARACTER(LEN=30)                                     :: temp_varname !! temporary variable string for rot_1d only
    REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE              :: manage
    INTEGER(i_std)                                        :: ip, kc, ier
    INTEGER(i_std),DIMENSION(rot_cmd_max)                 :: temp_cmd  !! temporary command variable for rot_1d only

!!!!-------------------------------------------------------------------------------------------

            !!!! rotation command
            IF (rot_1d) THEN
                DO kc = 1,cyc_rot_max
                    WRITE(temp_varname,"('CMDROTATE_',i1)"),kc
                    CALL getin_p(temp_varname,temp_cmd)
                    !!!! confirm how getin_p perform when the input command is less than rot_cmd_max
                    WRITE(numout,*) 'xuhui: kc, ', kc, 'temp_cmd, ', temp_cmd
                    DO ip = 1,nbpt
                        rot_cmd_store(ip,:,kc) = temp_cmd
                    ENDDO
                ENDDO
            ELSE ! reading rotation map
                strManage = 'CMDROTATE_FILE'
                strVar = 'Command'
                CALL stomate_stics_ManageInput(nbpt,lalo,neighbours,resolution,contfrac,strManage,strVar,manage,yrlen)
                !!! dimension of manage: kjpindex, cm_max, cyc_rot_max
                IF ( SIZE(manage,3,i_std) .NE. cyc_rot_max ) THEN
                    WRITE(numout,*) 'SIZE(manage,3), cyc_rot_max', SIZE(manage,3,i_std), cyc_rot_max
                    STOP 'input rotation cycle did not match cyc_rot_max'
                ENDIF
                IF ( SIZE(manage,2,i_std) .GT. rot_cmd_max ) THEN
                    IF (SIZE(manage,2,i_std) .GT. 10 ) THEN
                        WRITE(numout,*) 'WARNING: more than 10 commands for one rotation'
                        WRITE(numout,*) 'xuhui: memory issue may occur. Are you sure about so many rotation commands?'
                    ENDIF
                    WRITE(numout,*) 'xuhui: update rot_cmd_max'
                    WRITE(numout,*) ' SIZE(manage,2), rot_cmd_max', SIZE(manage,2,i_std), rot_cmd_max
                    rot_cmd_max = SIZE(manage,2,i_std)
                    DEALLOCATE(rot_cmd_store)
                    ALLOCATE(rot_cmd_store(nbpt,rot_cmd_max,cyc_rot_max),stat=ier)
                    IF (ier/=0) THEN
                        WRITE(numout,*) "ERROR IN RE-ALLOCATION OF rot_cmd_store: ",ier
                        STOP 'stomate_init rot_cmd_store reallocate'
                    ENDIF
                ENDIF
                IF (SIZE(manage,2,i_std) .LT. rot_cmd_max) THEN
                    rot_cmd_store(:,1:SIZE(manage,2,i_std),:) = INT(manage,i_std)
                    rot_cmd_store(:,SIZE(manage,2,i_std)+1:rot_cmd_max,:) = 0
                ELSE  !!! size equivalent
                    rot_cmd_store = INT(manage,i_std)
                ENDIF
            ENDIF
            IF (ALLOCATED(manage)) DEALLOCATE(manage) ! clear manage for other input purpose

  END SUBROUTINE stomate_stics_read_rotation

!! ================================================================================================================================
!! SUBROUTINE 	: stomate_stics_read_plantdate 
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
  SUBROUTINE stomate_stics_read_plantdate(nbpt,lalo,neighbours,resolution,contfrac, & ! in
                                    cyc_num, &
                                    plantdate, plantdate_now) ! out
    ! 0.1 INPUT
    !    
    INTEGER(i_std), INTENT(in)  :: nbpt         ! Number of points for which the data needs to be interpolated (extracted)
    REAL(r_std), INTENT(in) :: lalo(nbpt,2)     ! Vector of latitude and longtitude
    INTEGER(i_std), INTENT(in)  :: neighbours(nbpt,8)   ! Vectors of neighbours for each grid point
    ! (1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW)
    REAL(r_std),INTENT(in)  :: resolution(nbpt,2)   ! The size in km of each grid box in lat and lon
    !REAL(r_std)             :: resolution(nbpt,2)   ! The size in km of each     !grid box in lat and lon
    REAL(r_std),INTENT(in)  :: contfrac(nbpt)   ! The fraction of land in each grid box
    INTEGER(i_std),INTENT(in), DIMENSION(nbpt,nvm) :: cyc_num ! flag for starting the rotation for kjpindex, nvm

    INTEGER(i_std), INTENT(out), DIMENSION (nbpt,nvm,cyc_rot_max)  :: plantdate !! kjpindex, nvm, cyc_rot_max
    INTEGER(i_std), INTENT(out), DIMENSION (nbpt,nvm)  :: plantdate_now !! kjpindex, nvm, plantdate of current cycle

    ! 0.3 local variables
    INTEGER(i_std)                                        :: yrlen
    CHARACTER(LEN=30)                                     :: strManage
    CHARACTER(LEN=30)                                     :: strVar
    REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE              :: manage
    INTEGER(i_std)                                        :: ip, kc, ier, j

!!!!-------------------------------------------------------------------------------------------

        !iplt_1d = .TRUE.
        CALL getin_p('IMPOSE_IPLT',iplt_1d)
        IF (.not. iplt_1d) THEN
            strManage = "IPLT_FILE"
            strVar = "PLNTDT"
            CALL stomate_stics_ManageInput(nbpt,lalo,neighbours,resolution,contfrac,strManage,strVar,manage,yrlen)
            !manage(manage=val_exp) = 1 !when there is exceptional value, we applied 1
            IF ( SIZE(manage,3,i_std) .NE. cyc_rot_max ) THEN
                WRITE(numout,*) 'cyc_rot_max, SIZE(manage,3)', cyc_rot_max, SIZE(manage,3,i_std)
                STOP 'input and specified maximum rotation cycle not matched'
            ENDIF
    
            !!!! deal with exceptional values
            DO ip = 1,nbpt
                DO j = 2,nvm
                    IF ( ok_LAIdev(j) ) THEN
                        DO kc = 1,cyc_rot_max
                            IF ( (manage(ip,j,kc)==val_exp) .OR. (manage(ip,j,kc) .LT. zero) ) THEN
                                manage(ip,j,kc) = 1 + year_length_in_days*(kc-1)
                            ENDIF
                        ENDDO
                    ENDIF
                ENDDO
            ENDDO
    
            DO j = 2,nvm
                IF (nvm_plnt) THEN
                    IF (.NOT. (nvm .EQ. SIZE(manage,2,i_std)) ) THEN
                        WRITE(numout,*) 'nvm_plnt: ', nvm_plnt
                        WRITE(numout,*) 'nvm, SIZE(manage,2)', nvm, SIZE(manage,2,i_std)
                        STOP 'PFT more than input planting date'
                    ENDIF
                    WHERE(manage(:,j,:) >= val_exp)
                        manage(:,j,:) = val_exp
                    ENDWHERE
                    plantdate(:,j,:) = INT(manage(:,j,:),i_std)
                ELSE    
                    IF (MAXVAL(pft_to_mtc(:)) > SIZE(manage,2,i_std)) THEN
                        WRITE(numout,*) 'nvm_plnt: ', nvm_plnt
                        WRITE(numout,*) 'MAXVAL(pft_to_mtc(:)), SIZE(manage,2)', MAXVAL(pft_to_mtc(:)), SIZE(manage,2,i_std)
                        STOP 'PFT more than input planting date'
                    ENDIF
                    WHERE(manage(:,j,:) >= val_exp)
                        manage(:,j,:) = val_exp
                    ENDWHERE
                    plantdate(:,j,:) = INT(manage(:,pft_to_mtc(j),:),i_std) 
                ENDIF
            ENDDO
        ELSE ! iplt_1d
            IF (nvm_plnt) THEN
                IF (nvm .GT. SIZE(SP_iplt0,1,i_std)) THEN
                    WRITE(numout,*) 'nvm_plnt: ', nvm_plnt
                    WRITE(numout,*) 'nvm, SIZE(SP_iplt0,1)', nvm, SIZE(SP_iplt0,1,i_std)
                    STOP 'PFT more than input planting date'
                ENDIF
                DO j = 2,nvm
                    plantdate(:,j,1) = SP_iplt0(j)
                    IF ( cyc_rot_max .GT. 1) THEN
                        plantdate(:,j,2) = SP_iplt1(j)
                    ENDIF
                    IF ( cyc_rot_max .GT. 2) THEN
                        plantdate(:,j,3) = SP_iplt2(j)
                    ENDIF
                ENDDO
            ELSE !!!! using pft_to_mtc
                WRITE(numout,*) 'SP_iplt0, ', SP_iplt0
                IF (MAXVAL(pft_to_mtc(:)) .GT.  SIZE(SP_iplt0,1,i_std)) THEN
                    WRITE(numout,*) 'nvm_plnt: ', nvm_plnt
                    WRITE(numout,*) 'MAXVAL(pft_to_mtc(:)), SIZE(SP_iplt0,1)', MAXVAL(pft_to_mtc(:)), SIZE(SP_iplt0,1,i_std)
                    STOP 'PFT more than input planting date'
                ENDIF
                DO j = 2,nvm
                    plantdate(:,j,1) = SP_iplt0(pft_to_mtc(j))
                    IF ( cyc_rot_max .GT. 1) THEN
                        plantdate(:,j,2) = SP_iplt1(pft_to_mtc(j))
                    ENDIF
                    IF ( cyc_rot_max .GT. 2) THEN
                        plantdate(:,j,3) = SP_iplt2(pft_to_mtc(j))
                    ENDIF
                ENDDO
            ENDIF
        ENDIF ! iplt_1d

        IF (ok_rotate) THEN
            DO ip = 1,nbpt
                DO j = 2,nvm
                    IF (ok_LAIdev(j)) THEN
                        plantdate_now(ip,j) = plantdate(ip,j,cyc_num(ip,j))
                    ENDIF
                ENDDO
            ENDDO
        ELSE
            plantdate_now = plantdate(:,:,1)
        ENDIF
  END SUBROUTINE stomate_stics_read_plantdate

!! ================================================================================================================================
!! SUBROUTINE 	: stomate_stics_get_cmd
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
  SUBROUTINE stomate_stics_get_cmd(cmdin, src_rot, tgt_rot, prc_rot)
  ! 0.1 Input
  INTEGER(i_std),INTENT(in)  :: cmdin
  ! 0.2 Output
  INTEGER(i_std),INTENT(out) :: src_rot, tgt_rot
  REAL(r_std),INTENT(out)    :: prc_rot

  ! 0.3 Local variables

  !!!! --------------------------------------------------------------
    IF (cmdin > 1010000 .OR. cmdin < 0 ) THEN
        WRITE(numout,*) 'cmdin, ',cmdin
        STOP 'cmd error in stomate_stics_get_cmd'
    ENDIF
    IF (cmdin == 0) THEN
        tgt_rot = 0
        src_rot = 0
        
        prc_rot = 0.0
    ELSE
        tgt_rot = MOD(cmdin, 100)
        src_rot = MOD(FLOOR(FLOAT(cmdin)/100), 100)

        prc_rot = FLOAT(FLOOR(FLOAT(cmdin)/10000))/100.0
        IF (printlev >=4) THEN
            WRITE(numout,*) 'xuhui: cmdin, tgt_rot, src_rot, prc_rot', cmdin, src_rot, tgt_rot, prc_rot
        ENDIF
        IF (prc_rot .GT. 1.0 .AND. prc_rot .LT. 1.0+0.01) THEN ! resolve potential  precision issues
            prc_rot = 1.0
        ENDIF
        !!! consistency check
        IF (prc_rot .GT. 1.0) THEN
            WRITE(numout,*) 'percent change larger than 1..., prc_rot',prc_rot
            STOP 'incorrect percent rotation, stomate_stics_get_cmd'
        ENDIF
        IF ( (tgt_rot .GT. nvm) .OR. ( .NOT. (ok_LAIdev(tgt_rot) .OR. (tgt_rot .EQ. 1)) ) ) THEN
            WRITE(numout,*) 'rotation target error: tgt_rot ', tgt_rot
            WRITE(numout,*) 'nvm, ok_LAIdev', nvm, ok_LAIdev
            STOP 'incorrect rotation target, stomate_stics_get_cmd'
        ENDIF
        IF ( (src_rot .GT. nvm) .OR. ( .NOT. (ok_LAIdev(src_rot) .OR. (src_rot .EQ. 1)) ) ) THEN
            WRITE(numout,*) 'rotation target error: src_rot ', src_rot
            WRITE(numout,*) 'nvm, ok_LAIdev', nvm, ok_LAIdev
            STOP 'incorrect rotation source, stomate_stics_get_cmd'
        ENDIF
    ENDIF
  END SUBROUTINE stomate_stics_get_cmd


!! ================================================================================================================================
!! SUBROUTINE 	: stomate_stics_rotation
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
  SUBROUTINE stomate_stics_rotation(kjpindex, ip, matrix_prc, maxfrac,  & ! in
                              in_cycle, cyc_num_tot, plantdate, nrec, nlev, &
                              plantdate_now, &   ! inout
                              soilc_mict, cyc_num, litter, turnover_daily, & 
                              deepC_a, deepC_s, deepC_p, carbon, &
                              age, ind, PFTpresent, senescence, &
                              when_growthinit, everywhere, leaf_frac  )
  !!! This subroutine transfer litter and soil carbon when the cropland is rotated
  !!! update plantdate_now
  !!! veget_max and soil water and thermo properties will be transferred in sechiba_main
  ! 0.1 Input
  INTEGER(i_std),INTENT(in)                     :: kjpindex ! number of land points
  INTEGER(i_std),INTENT(in)                     :: ip ! target point
  REAL(r_std),DIMENSION(nvm,nvm),INTENT(in)     :: matrix_prc ! rotation matrix in % of veget_max
  REAL(r_std),DIMENSION(nvm),INTENT(in)         :: maxfrac ! veget_max before rotation
  LOGICAL, DIMENSION(kjpindex, nvm), INTENT(in)          :: in_cycle 
  INTEGER(i_std), DIMENSION(kjpindex), INTENT(in)   :: cyc_num_tot ! number of current rotation cycle
  INTEGER(i_std), DIMENSION(kjpindex,nvm,cyc_rot_max), INTENT(in)  :: plantdate !! kjpindex, nvm, cyc_rot_max
  INTEGER(i_std), DIMENSION(kjpindex, nvm), INTENT(in)        :: nrec
  INTEGER(i_std), DIMENSION(kjpindex, nvm), INTENT(in)        :: nlev

  ! 0.2 Output
  INTEGER(i_std), DIMENSION (kjpindex,nvm), INTENT(inout)  :: plantdate_now !! kjpindex, nvm, plantdate of current cycle
  REAL(r_std), DIMENSION(ndeep,nvm), INTENT(out)  :: soilc_mict
  INTEGER(i_std),INTENT(inout), DIMENSION(kjpindex,nvm) :: cyc_num ! flag for starting the rotation for kjpindex, nvm
  REAL(r_std), INTENT(inout), DIMENSION(:,:,:,:,:)  :: litter ! Above and below ground metabolic and structural litter per ground area 
  REAL(r_std), INTENT(inout), DIMENSION(:,:,:,:)    :: turnover_daily ! Senescence-driven turnover (better: mortality) of leaves and roots  
  REAL(r_std), INTENT(inout), DIMENSION(:,:,:)      :: deepC_a        ! deep active carbon profile
  REAL(r_std), INTENT(inout), DIMENSION(:,:,:)      :: deepC_s        ! deep slow carbon profile
  REAL(r_std), INTENT(inout), DIMENSION(:,:,:)      :: deepC_p        ! deep passive carbon profile
  REAL(r_std), INTENT(inout), DIMENSION(:,:,:)      :: carbon          ! Soil carbon pools per ground area: active, slow, or passive 
  REAL(r_std), INTENT(inout), DIMENSION(:,:)        :: age                  !! Age of PFT it normalized by biomass - can increase and decrease - (years)
  REAL(r_std), INTENT(inout), DIMENSION(:,:)        :: ind                  !! Vegetation density, number of individuals per unit 
  LOGICAL, INTENT(inout), DIMENSION(:,:)            :: PFTpresent           !! PFT exists (equivalent to veget > 0 for natural PFTs)
  LOGICAL, INTENT(inout), DIMENSION(:,:)            :: senescence           !! The PFT is senescent
  REAL(r_std), INTENT(inout), DIMENSION(:,:)        :: when_growthinit      !! Days since beginning of growing season (days)
  REAL(r_std), INTENT(inout), DIMENSION(:,:)    :: everywhere           !! Is the PFT everywhere in the grid box or very localized 
  REAL(r_std), INTENT(inout), DIMENSION(:,:,:)  :: leaf_frac            !! PFT fraction of leaf mass in leaf age class (0-1, 


  ! 0.3 Local Variables
  REAL(r_std),DIMENSION(nvm,nlitt,nlevs,nelements)         :: dilu_lit !! Litter dilution
  REAL(r_std),DIMENSION(nvm,nparts,nelements)              :: dilu_turn !! Turnover dilution
  REAL(r_std),DIMENSION(nvm,ncarb)                         :: dilu_soil_carbon
  REAL(r_std),DIMENSION(nvm,ndeep,ncarb)                   :: dilu_soil_carbon_vertres !!vertically-resolved Soil C arbondilution (gC/m²)

  REAL(r_std),DIMENSION(ncarb,nvm)                 :: carbon_old
  REAL(r_std),DIMENSION(nlitt,nvm,nlevs,nelements) :: litter_old
  REAL(r_std),DIMENSION(nvm,nparts,nelements)      :: turnover_old
  REAL(r_std), DIMENSION(ndeep,nvm)                :: deepC_a_old
  REAL(r_std), DIMENSION(ndeep,nvm)                :: deepC_s_old 
  REAL(r_std), DIMENSION(ndeep,nvm)                :: deepC_p_old 

  REAL(r_std),DIMENSION(nvm)                    :: maxfrac_new ! veget_max after rotation
  INTEGER(i_std)                                :: i,jtar,jsrc,tempcyc
  INTEGER(i_std),DIMENSION(nvm)        :: cyc_num_old
  INTEGER(i_std)                                :: rngTol = 14
  LOGICAL                                       :: f_convert

  !!!! -------------------------------------------------------

  !1.0 update cyc_num plantdate_now
!  temp_cyc = cyc_num
  cyc_num_old = cyc_num(ip,:)
  !!! source crops
  DO jsrc = 2,nvm
    IF (SUM(matrix_prc(jsrc,:)) .GT. min_stomate) THEN
    !!!! rotation from the jsrc, clean the rotation cycle of this PFT
        cyc_num(ip,jsrc) = 1
    ENDIF
  ENDDO

  maxfrac_new = maxfrac
  litter_old = litter(ip,:,:,:,:)
  turnover_old = turnover_daily(ip,:,:,:)
  IF (ok_pc) THEN
    deepC_a_old = deepC_a(ip,:,:)
    deepC_s_old = deepC_s(ip,:,:)
    deepC_p_old = deepC_p(ip,:,:)
  ELSE
    carbon_old = carbon(ip,:,:)
  ENDIF
  !!! target crops
  DO jtar = 1, nvm
    dilu_lit(:,:,:,:) = zero
    dilu_turn(:,:,:) = zero
    dilu_soil_carbon(:,:) = zero
    dilu_soil_carbon_vertres(:,:,:) = zero
    tempcyc = cyc_num_old(jtar)
    IF (SUM(matrix_prc(:,jtar)) .GT. min_stomate) THEN
        DO jsrc = 2, nvm
            IF (matrix_prc(jsrc,jtar) .GT. min_stomate) THEN
                f_convert = .FALSE.
                IF (in_cycle(ip,jtar)) THEN 
                    ! a complex case when the target PFT is still growing
                    ! this is a bad case we do not allow
                    WRITE(numout,*) 'rotate to a growing PFT: jsrc, jtar',jsrc, jtar
                    WRITE(numout,*) 'in_cycle(ip,jtar) ', in_cycle(ip,jtar)
                    STOP 'non-permitted rotation case occurred'
                ELSE
                    IF ( cyc_num(ip,jtar) == 1 ) THEN 
                        !! simple case when target PFT is awating start
                        IF (cyc_num_old(jsrc) .LT. cyc_num_tot(ip)) THEN
                            cyc_num(ip,jtar) = cyc_num_old(jsrc) + 1
                        ELSE ! complete of one rotation cycle, back to 1
                            cyc_num(ip,jtar) = 1
                        ENDIF
                        plantdate_now(ip,jtar) = plantdate(ip,jtar,cyc_num(ip,jtar))
                        IF ( (nrec(ip,jsrc) < nlev(ip,jsrc)) .AND. (nrec(ip,jsrc) .GT. 0) ) THEN 
                            !harvest occured in the next-year of planting, adjusting planting date next cycle
                            IF (plantdate_now(ip,jtar) > year_length_in_days) plantdate_now(ip,jtar) = plantdate_now(ip,jtar) - year_length_in_days
                        ENDIF
                        !!! when next planting date passed but within tolerance range
                        IF ( (nrec(ip,jsrc) .GE. plantdate_now(ip,jtar)) .AND. (nrec(ip,jsrc) .GT. 0) .AND.  &
                             (nrec(ip,jsrc) .LT. rngTol+plantdate_now(ip,jtar)) ) THEN
                            plantdate_now(ip,jtar) = nrec(ip,jsrc)+2  ! we still grow it the day after tomorrow
                        ENDIF
                        maxfrac_new(jtar) = maxfrac_new(jtar) + matrix_prc(jsrc,jtar) * maxfrac(jsrc)
                        maxfrac_new(jsrc) = maxfrac_new(jsrc) - matrix_prc(jsrc,jtar) * maxfrac(jsrc)
                        f_convert = .TRUE.
                    ELSEIF ( (tempcyc-1 .EQ. cyc_num_old(jsrc)) .OR. (tempcyc==1 .AND. cyc_num_old(jsrc)==cyc_num_tot(ip)) ) THEN
                        ! one cycle ahead waiting to start
                        ! no need to change cyc_num & planting date
                        maxfrac_new(jtar) = maxfrac_new(jtar) + matrix_prc(jsrc,jtar) * maxfrac(jsrc)
                        maxfrac_new(jsrc) = maxfrac_new(jsrc) - matrix_prc(jsrc,jtar) * maxfrac(jsrc)
                        f_convert = .TRUE.
                    ELSE
                        WRITE(numout,*) 'rotation situation unexpected:'
                        WRITE(numout,*) 'cyc_num_old(ip,jsrc), cyc_num_old(ip,jtar)',cyc_num_old(jsrc), cyc_num_old(jtar)
                        WRITE(numout,*) 'cyc_num(ip,jtar)',cyc_num(ip,jtar)
                        WRITE(numout,*) 'maxfrac',maxfrac
                        STOP 'stomate_stics_rotation error with unexpected command'
                    ENDIF ! tempcyc == 1
                ENDIF ! in_cycle(jtar)

                IF (f_convert) THEN
                !!!! transfering soil and litter/turnover carbon pools
                    !!! this is actually the better place to modify f_rot_sech & rot_cmd
                    
                    !!! first is source crops
                    dilu_lit(jsrc,:,:,:) = litter_old(:,jsrc,:,:)
                    dilu_turn(jsrc,:,:) = turnover_old(jsrc,:,:)
                    IF (ok_pc) THEN
                        dilu_soil_carbon_vertres(jsrc,:,iactive) = deepC_a_old(:,jsrc)
                        dilu_soil_carbon_vertres(jsrc,:,islow) = deepC_s_old(:,jsrc)
                        dilu_soil_carbon_vertres(jsrc,:,ipassive) = deepC_p_old(:,jsrc)
                    ELSE    
                        dilu_soil_carbon(jsrc,:) = carbon_old(:,jsrc)
                    ENDIF
                ENDIF ! we indeed converting jsrc to jtar
            ENDIF !jsrc -> jtar
        ENDDO !jsrc
        !!! then is the conversion of litter and carbon to target crop
        ! by design, SUM(matrix_prc(jtar,:)), should be either 0 or 1, this
        ! is simply in case for more complex situations...
        litter(ip,:,jtar,:,:) = litter_old(:,jtar,:,:) *  maxfrac(jtar) * (1.0 - SUM(matrix_prc(jtar,:))) 
        turnover_daily(ip,jtar,:,:) = turnover_old(jtar,:,:) * maxfrac(jtar) * (1.0 - SUM(matrix_prc(jtar,:)))
        IF (ok_pc) THEN
            deepC_a(ip,:,jtar) = deepC_a_old(:,jtar) * maxfrac(jtar) * (1.0 - SUM(matrix_prc(jtar,:)))
            deepC_s(ip,:,jtar) = deepC_s_old(:,jtar) * maxfrac(jtar) * (1.0 - SUM(matrix_prc(jtar,:)))
            deepC_p(ip,:,jtar) = deepC_p_old(:,jtar) * maxfrac(jtar) * (1.0 - SUM(matrix_prc(jtar,:)))
        ELSE
            carbon(ip,:,jtar) = carbon_old(:,jtar) * maxfrac(jtar) * (1.0 - SUM(matrix_prc(jtar,:)))
        ENDIF
        DO jsrc = 1,nvm
            litter(ip,:,jtar,:,:) = litter(ip,:,jtar,:,:) + (maxfrac(jsrc) * matrix_prc(jsrc,jtar) * dilu_lit(jsrc,:,:,:)) 
            turnover_daily(ip,jtar,:,:) = turnover_daily(ip,jtar,:,:) + (maxfrac(jsrc) * matrix_prc(jsrc,jtar) * dilu_turn(jsrc,:,:)) 
            IF (ok_pc) THEN
                deepC_a(ip,:,jtar) = deepC_a(ip,:,jtar) + (maxfrac(jsrc) * matrix_prc(jsrc,jtar) * dilu_soil_carbon_vertres(jsrc,:,iactive))
                deepC_s(ip,:,jtar) = deepC_s(ip,:,jtar) + (maxfrac(jsrc) * matrix_prc(jsrc,jtar) * dilu_soil_carbon_vertres(jsrc,:,islow))
                deepC_p(ip,:,jtar) = deepC_p(ip,:,jtar) + (maxfrac(jsrc) * matrix_prc(jsrc,jtar) * dilu_soil_carbon_vertres(jsrc,:,ipassive))
            ELSE
                carbon(ip,:,jtar) = carbon(ip,:,jtar) + (maxfrac(jsrc) * matrix_prc(jsrc,jtar) * dilu_soil_carbon(jsrc,:))
            ENDIF
        ENDDO
        litter(ip,:,jtar,:,:) = litter(ip,:,jtar,:,:) /  maxfrac_new(jtar)
        turnover_daily(ip,jtar,:,:) = turnover_daily(ip,jtar,:,:) /  maxfrac_new(jtar)
        IF (ok_pc) THEN
            deepC_a(ip,:,jtar) = deepC_a(ip,:,jtar) / maxfrac_new(jtar)
            deepC_s(ip,:,jtar) = deepC_s(ip,:,jtar) / maxfrac_new(jtar)
            deepC_p(ip,:,jtar) = deepC_p(ip,:,jtar) / maxfrac_new(jtar)
        ELSE
            carbon(ip,:,jtar) = carbon(ip,:,jtar) / maxfrac_new(jtar)
        ENDIF
        age(ip,jtar) = zero

        IF ( (maxfrac(jtar) .LT. min_stomate) .AND. (maxfrac_new(jtar) .GT. min_stomate) ) THEN !initialize a previously non-existed PFT
            ind(ip,jtar) = maxfrac_new(jtar)
            PFTpresent(ip,jtar) = .TRUE.
            !everywhere(ip,jtar) = 1.
            senescence(ip,jtar) = .TRUE.
            age(ip,jtar) = zero
            when_growthinit(ip,jtar) = large_value
            leaf_frac(ip,jtar,1) = 1.0
        ENDIF

    ENDIF ! to jtar > 0
  ENDDO !jtar
  IF (printlev>=4) THEN
    WRITE(numout,*) 'maxfrac', maxfrac
    WRITE(numout,*) 'maxfrac_new', maxfrac_new
    DO i = 1,3
        WRITE(numout,*) 'i, carbon_old(i,:)', i, carbon_old(i,:)
        WRITE(numout,*) 'i, carbon(ip,i,:)', i, carbon(ip,i,:)
    ENDDO
    WRITE(numout,*) 'SUM(litter_old(:,:,:,icarbon))',SUM(litter_old(:,:,:,icarbon))
    WRITE(numout,*) 'SUM(litter(ip,:,:,:,icarbon))',SUM(litter(ip,:,:,:,icarbon))
  ENDIF

  DO jsrc = 1,nvm
    IF ( maxfrac_new(jsrc)  .LT.  min_stomate ) THEN
    ! clean the litter and soil carbon pool of the  PFT no longer existed
        cyc_num(ip,jsrc) = 1
        litter(ip,:,jsrc,:,:) = zero
        turnover_daily(ip,jsrc,:,:) = zero
        IF (ok_pc) THEN
            deepC_a(ip,:,jsrc) = zero
            deepC_s(ip,:,jsrc) = zero
            deepC_p(ip,:,jsrc) = zero
            
        ELSE
            carbon(ip,:,jsrc) = zero
        ENDIF
        ind(ip,jsrc) = zero
        PFTpresent(ip,jsrc) = .FALSE.
        senescence(ip,jsrc) = .FALSE.
        age(ip,jsrc) = zero
        when_growthinit(ip,jsrc) = large_value
        everywhere(ip,jsrc) = zero
    ENDIF
  ENDDO

  soilc_mict(:,:) = deepC_a(ip,:,:) + deepC_s(ip,:,:) + deepC_p(ip,:,:)

  !!! now check the maxfrac_new for consistency
  IF ( SUM(maxfrac_new(:)) .GT. 1.0 ) THEN
    WRITE(numout,*) 'vegetation fraction greater than 1.0 after rotation'
    WRITE(numout,*) 'ip, maxfrac_new,', ip, maxfrac_new
    STOP 'stomate_stics_rotation: wrong vegetation fraction'
  ENDIF

  END SUBROUTINE stomate_stics_rotation

END MODULE stomate_stics
