! =================================================================================================================================
! MODULE       : stomate_gluc_common
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       This module contains common fuctions and subroutines used by
!              gross land use change and forestry harvest modules.
!!
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): 
!!
!! REFERENCE(S)	: None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/perso/albert.jornet/ORCHIDEE-MICT/src_stomate/stomate_lcchange.f90 $
!! $Date: 2015-07-30 15:38:45 +0200 (Thu, 30 Jul 2015) $
!! $Revision: 2847 $
!! \n
!_ ================================================================================================================================


MODULE stomate_gluc_constants

  ! modules used:
  
  USE pft_parameters
  USE constantes
  USE constantes_soil_var
  
  IMPLICIT NONE
  
  PUBLIC
  
  INTEGER, ALLOCATABLE, SAVE                  :: indall_tree(:)       !! Indices for all tree PFTs
  INTEGER, ALLOCATABLE, SAVE                  :: indold_tree(:)       !! Indices for old tree cohort only
  INTEGER, ALLOCATABLE, SAVE                  :: indagec_tree(:,:)    !! Indices for secondary tree cohorts, 
                                                                      !! note the sequence is old->young.
  INTEGER, ALLOCATABLE, SAVE                  :: indall_grass(:)      !! Indices for all grass PFTs
  INTEGER, ALLOCATABLE, SAVE                  :: indold_grass(:)      !! Indices for old grasses only
  INTEGER, ALLOCATABLE, SAVE                  :: indagec_grass(:,:)   !! Indices for secondary grass cohorts
                                                                      !! note the sequence is old->young.
  INTEGER, ALLOCATABLE, SAVE                  :: indall_pasture(:)    !! Indices for all pasture PFTs
  INTEGER, ALLOCATABLE, SAVE                  :: indold_pasture(:)    !! Indices for old pasture only
  INTEGER, ALLOCATABLE, SAVE                  :: indagec_pasture(:,:) !! Indices for secondary pasture cohorts
                                                                      !! note the sequence is old->young.
  INTEGER, ALLOCATABLE, SAVE                  :: indall_crop(:)       !! Indices for all crop PFTs
  INTEGER, ALLOCATABLE, SAVE                  :: indold_crop(:)       !! Indices for old crops only
  INTEGER, ALLOCATABLE, SAVE                  :: indagec_crop(:,:)    !! Indices for secondary crop cohorts

  INTEGER, ALLOCATABLE, SAVE                  :: indall_bioe1(:)       !! Indices for all bioe1 PFTs
  INTEGER, ALLOCATABLE, SAVE                  :: indold_bioe1(:)       !! Indices for old bioe1 only
  INTEGER, ALLOCATABLE, SAVE                  :: indagec_bioe1(:,:)    !! Indices for secondary bioe1 cohorts

  INTEGER, SAVE :: num_tree_mulagec,num_grass_mulagec,     &
                   num_pasture_mulagec,num_crop_mulagec,   &
                   num_bioe1_mulagec
  
CONTAINS

! ================================================================================================================================
!! SUBROUTINE   : stomate_gluc_constants_init
!!
!>\BRIEF        Calculate coverage fraction for different age classes of forest,
!!              grass, pasture and crops and also for each metaclass. Note baresoil is excluded. 
!!              
!! DESCRIPTION :
!! Note:
!! 1. "calc_cover" subroutine does not depend on how many age classes
!! there are in each MTC.
!! 2. Fraction of baresoil is excluded here. This means transformation
!! of baresoil to a vegetated PFT is excluded in gross land cover change.
!!  
!!
!! MAIN OUTPUT VARIABLE(S) :  
!!
!! \n
!_ ================================================================================================================================

  SUBROUTINE stomate_gluc_constants_init()

    INTEGER(i_std) :: itree,itree2,igrass,igrass2,ipasture,ipasture2,icrop,icrop2, &
                      ibioe1, ibioe12
    INTEGER(i_std) :: i,j,ivma,staind,endind,ivm
    INTEGER(i_std) :: ier               !! Check errors in netcdf call


    !! 1. We first build all different indices that we are going to use
    !!    in handling the PFT exchanges, three types of indices are built:
    !!     - for all age classes
    !!     - include only oldest age classes
    !!     - include all age classes excpet the oldest ones
    ! We have to build these indices because we would like to extract from
    ! donating PFTs in the sequnce of old->young age classes or the revserse, 
    ! and add in the receving PFTs only in the youngest-age-class PFTs. These 
    ! indicies allow us to know where the different age classes are.

    ! calculate the total number of MTCs for each vegetation type.
    num_tree_mulagec=0          
    num_grass_mulagec=0
    num_pasture_mulagec=0
    num_crop_mulagec=0
    num_bioe1_mulagec=0
    
    !! 1.1 Calculate the number of PFTs for different MTCs and allocate
    !! the old and all indices arrays.

    ! [Note here the sequence to identify tree,pasture,grass,crop] is
    ! critical. The similar sequence is used in the subroutine "calc_cover".
    ! Do not forget to change the sequence there if you modify here.
    DO ivma =2,nvmap
      staind=start_index(ivma)
      IF (nagec_pft(ivma)==1) THEN
        WRITE(numout,*) "Error: metaclass has only a single age group: ",ivma
        WRITE(numout,*) "stomate_gluc_constants_init"
        STOP
      ELSE
        IF (is_tree(staind)) THEN
          num_tree_mulagec = num_tree_mulagec+1
        ELSE IF (is_grassland_manag(staind)) THEN
          num_pasture_mulagec = num_pasture_mulagec+1
        ELSE IF (is_bioe1(staind)) THEN
          num_bioe1_mulagec = num_bioe1_mulagec+1
        ELSE IF (natural(staind)) THEN
          num_grass_mulagec = num_grass_mulagec+1
        ELSE
          num_crop_mulagec = num_crop_mulagec+1
        ENDIF
      ENDIF
    ENDDO
    
    !! Allocate index array
    ! allocate all index

    ALLOCATE(indall_tree(num_tree_mulagec*nagec_tree),stat=ier)
    IF (ier .NE. 0) THEN
       WRITE(numout,*) 'Memory allocation error for indall_tree'
       STOP 'stomate_gluc_constants_init'
    ENDIF

    ALLOCATE(indall_grass(num_grass_mulagec*nagec_herb),stat=ier)
    IF (ier .NE. 0) THEN
       WRITE(numout,*) 'Memory allocation error for indall_grass'
       STOP 'stomate_gluc_constants_init'
    ENDIF

    ALLOCATE(indall_pasture(num_pasture_mulagec*nagec_herb),stat=ier)
    IF (ier .NE. 0) THEN
       WRITE(numout,*) 'Memory allocation error for indall_pasture'
       STOP 'stomate_gluc_constants_init'
    ENDIF

    ALLOCATE(indall_crop(num_crop_mulagec*nagec_herb),stat=ier)
    IF (ier .NE. 0) THEN
       WRITE(numout,*) 'Memory allocation error for indall_crop'
       STOP 'stomate_gluc_constants_init'
    ENDIF

    ! allocate old-ageclass index
    ALLOCATE(indold_tree(num_tree_mulagec),stat=ier)
    IF (ier .NE. 0) THEN
       WRITE(numout,*) 'Memory allocation error for indold_tree'
       STOP 'stomate_gluc_constants_init'
    ENDIF

    ALLOCATE(indold_grass(num_grass_mulagec),stat=ier)
    IF (ier .NE. 0) THEN
       WRITE(numout,*) 'Memory allocation error for indold_grass'
       STOP 'stomate_gluc_constants_init'
    ENDIF

    ALLOCATE(indold_pasture(num_pasture_mulagec),stat=ier)
    IF (ier .NE. 0) THEN
       WRITE(numout,*) 'Memory allocation error for indold_pasture'
       STOP 'stomate_gluc_constants_init'
    ENDIF

    ALLOCATE(indold_crop(num_crop_mulagec),stat=ier)
    IF (ier .NE. 0) THEN
       WRITE(numout,*) 'Memory allocation error for indold_crop'
       STOP 'stomate_gluc_constants_init'
    ENDIF

    ! To to be able to go through the following code, when there are no bio1
    ! meta-calsses, we set num_bioe1_mulagec as 1.
    IF (num_bioe1_mulagec .EQ. 0) num_bioe1_mulagec = 1

    ALLOCATE(indall_bioe1(num_bioe1_mulagec*nagec_bioe1),stat=ier)
    IF (ier .NE. 0) THEN
       WRITE(numout,*) 'Memory allocation error for indall_crop'
       STOP 'stomate_gluc_constants_init'
    ENDIF
    indall_bioe1(:) = 0

    ALLOCATE(indold_bioe1(num_bioe1_mulagec),stat=ier)
    IF (ier .NE. 0) THEN
       WRITE(numout,*) 'Memory allocation error for indold_crop'
       STOP 'stomate_gluc_constants_init'
    ENDIF
    indold_bioe1(:) = 0

    ALLOCATE(indagec_bioe1(num_bioe1_mulagec,nagec_bioe1-1))
    IF (ier .NE. 0) THEN
       WRITE(numout,*) 'Memory allocation error for indagec_crop'
       STOP 'stomate_gluc_constants_init'
    ENDIF
    indagec_bioe1(:,:) = 0

    !! 1.2 Fill the oldest-age-class and all index arrays
    itree=0
    igrass=0
    ipasture=0
    icrop=0
    ibioe1=0
    itree2=1
    igrass2=1
    ipasture2=1
    icrop2=1
    ibioe12=1
    DO ivma =2,nvmap
      staind=start_index(ivma)
      IF (is_tree(staind)) THEN
        itree=itree+1
        indold_tree(itree) = staind+nagec_pft(ivma)-1
        DO j = 0,nagec_pft(ivma)-1
          indall_tree(itree2+j) = staind+j
        ENDDO
        itree2=itree2+nagec_pft(ivma)
      ELSE IF (is_bioe1(staind)) THEN
        ibioe1 = ibioe1+1
        indold_bioe1(ipasture) = staind+nagec_pft(ivma)-1
        DO j = 0,nagec_pft(ivma)-1
          indall_bioe1(ibioe12+j) = staind+j
        ENDDO
        ibioe12=ibioe12+nagec_pft(ivma)
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

    ! allocate old-ageclass index
    ALLOCATE(indagec_tree(num_tree_mulagec,nagec_tree-1),stat=ier)     
    IF (ier .NE. 0) THEN
       WRITE(numout,*) 'Memory allocation error for indagec_tree'
       STOP 'stomate_gluc_constants_init'
    ENDIF

    ALLOCATE(indagec_grass(num_grass_mulagec,nagec_herb-1),stat=ier)     
    IF (ier .NE. 0) THEN
       WRITE(numout,*) 'Memory allocation error for indagec_grass'
       STOP 'stomate_gluc_constants_init'
    ENDIF

    ALLOCATE(indagec_pasture(num_pasture_mulagec,nagec_herb-1),stat=ier)
    IF (ier .NE. 0) THEN
       WRITE(numout,*) 'Memory allocation error for indagec_pasture'
       STOP 'stomate_gluc_constants_init'
    ENDIF

    ALLOCATE(indagec_crop(num_crop_mulagec,nagec_herb-1),stat=ier)
    IF (ier .NE. 0) THEN
       WRITE(numout,*) 'Memory allocation error for indagec_crop'
       STOP 'stomate_gluc_constants_init'
    ENDIF

    ! fill the non-oldest age class index arrays when number of age classes
    ! is more than 1.
    itree=0
    igrass=0
    ipasture=0
    icrop=0
    ibioe1=0
    DO ivma = 2,nvmap
      staind=start_index(ivma)
      IF (nagec_pft(ivma) > 1) THEN
        IF (is_tree(staind)) THEN
          itree=itree+1
          DO j = 1,nagec_tree-1
            indagec_tree(itree,j) = staind+nagec_tree-j-1
          ENDDO
        ELSE IF (is_bioe1(staind)) THEN
          ibioe1=ibioe1+1
          DO j = 1,nagec_bioe1-1
            indagec_bioe1(ibioe1,j) = staind+nagec_bioe1-j-1
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

    write (numout,*) "indices calculated"
  END SUBROUTINE stomate_gluc_constants_init

  SUBROUTINE stomate_gluc_constants_init_clear()

    IF (ALLOCATED(indall_tree)) DEALLOCATE(indall_tree)
    IF (ALLOCATED(indall_grass)) DEALLOCATE(indall_grass)
    IF (ALLOCATED(indall_pasture)) DEALLOCATE(indall_pasture)
    IF (ALLOCATED(indall_crop)) DEALLOCATE(indall_crop)

    IF (ALLOCATED(indold_tree)) DEALLOCATE(indold_tree)
    IF (ALLOCATED(indold_grass)) DEALLOCATE(indold_grass)
    IF (ALLOCATED(indold_pasture)) DEALLOCATE(indold_pasture)
    IF (ALLOCATED(indold_crop)) DEALLOCATE(indold_crop)

    IF (ALLOCATED(indagec_tree)) DEALLOCATE(indagec_tree)
    IF (ALLOCATED(indagec_grass)) DEALLOCATE(indagec_grass)
    IF (ALLOCATED(indagec_pasture)) DEALLOCATE(indagec_pasture)
    IF (ALLOCATED(indagec_crop)) DEALLOCATE(indagec_crop)

    IF (ALLOCATED(indall_bioe1)) DEALLOCATE(indall_bioe1)
    IF (ALLOCATED(indold_bioe1)) DEALLOCATE(indold_bioe1)
    IF (ALLOCATED(indagec_bioe1)) DEALLOCATE(indagec_bioe1)

  END SUBROUTINE stomate_gluc_constants_init_clear

END MODULE stomate_gluc_constants
