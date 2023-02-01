MODULE interregxy

  ! Modules used :

  USE constantes
  USE mod_orchidee_para
  USE grid
  USE module_llxy
  USE polygones

  IMPLICIT NONE

  PRIVATE
  PUBLIC interregxy_aggr2d, interregxy_aggrve

CONTAINS

  SUBROUTINE interregxy_aggr2d(nbpt, lalo, neighbours, resolution, contfrac, &
       &                     iml, jml, lon_rel, lat_rel, mask, callsign, &
       &                     incmax, indinc, areaoverlap, ok)
    !
    ! INPUT
    ! 
    INTEGER(i_std), INTENT(in)   :: nbpt                 ! Number of points for which the data needs to be interpolated
    REAL(r_std), INTENT(in)      :: lalo(nbpt,2)         ! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), INTENT(in)   :: neighbours(nbpt,NbNeighb)! Vector of neighbours for each grid point (1=N, 2=E, 3=S, 4=W)
    REAL(r_std), INTENT(in)      :: resolution(nbpt,2)   ! The size in km of each grid-box in X and Y
    REAL(r_std), INTENT(in)      :: contfrac(nbpt)       ! Fraction of land in each grid box.
    INTEGER(i_std), INTENT(in)   :: iml, jml             ! Size of the source grid
    REAL(r_std), INTENT(in)      :: lon_rel(iml, jml)    ! Longitudes for the source grid
    REAL(r_std), INTENT(in)      :: lat_rel(iml, jml)    ! Latitudes for the souce grid
    INTEGER(i_std), INTENT(in)   :: mask(iml, jml)       ! Mask which retains only the significative points
                                                         ! of the source grid.
    CHARACTER(LEN=*), INTENT(in) :: callsign             ! Allows to specify which variable is beeing treated
    INTEGER(i_std), INTENT(in)   :: incmax               ! Maximum point of the source grid we can store.
    !
    ! Output
    !
    INTEGER(i_std), INTENT(out)    :: indinc(nbpt,incmax,2)
    REAL(r_std), INTENT(out)       :: areaoverlap(nbpt,incmax)
    LOGICAL, OPTIONAL, INTENT(out) :: ok                 ! return code
    !
    ! Local Variables
    !
    ! The number of sub-division we wish to have on each side of the rectangle.
    ! sbd=1 means that on each side we take only the middle point.
    !
    INTEGER(i_std), PARAMETER  :: sbd=5
    INTEGER(i_std), PARAMETER  :: idom=1
    !
    INTEGER(i_std)                                :: i, is, js, in
    INTEGER(i_std)                                :: ilow, iup, jlow, jup
    REAL(r_std)                                   :: dle, dlw, dls, dln
    !
    REAL(r_std), DIMENSION((sbd+1)*4+1)           :: lons, lats, ri, rj
    !
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:)   :: nbov
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:,:) :: sourcei, sourcej
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)    :: overlap
    LOGICAL                                       :: checkprint=.FALSE.
    !
    !
    REAL(r_std), DIMENSION(2)   :: BBlon, BBlat
    !
    WRITE(*,*) "interregxy_aggr2d working on ", callsign
    !
    ALLOCATE(sourcei(iim_g, jjm_g, incmax))
    ALLOCATE(sourcej(iim_g, jjm_g, incmax))
    ALLOCATE(nbov(iim_g, jjm_g))
    ALLOCATE(overlap(iim_g, jjm_g, incmax))
    nbov(:,:) = 0
    !
    BBlon(1) = MINVAL(lalo(:,2))
    BBlon(2) = MAXVAL(lalo(:,2))
    BBlat(1) = MINVAL(lalo(:,1))
    BBlat(2) = MAXVAL(lalo(:,1))
    !
    ! Go through the entire source grid and try to place each grid box in a target grid box.
    !
    DO is=1,iml
       DO js=1,jml
          !
          ! 1.0 Get the center of the box and their corresponding grid point 
          !
          lons(1) = lon_rel(is,js)
          lats(1) = lat_rel(is,js)
          !
          ! Assume all grids go from West to Est
          dle=(lon_rel(MIN(is+1,iml),js)-lon_rel(is,js))/2.0
          dlw=(lon_rel(is,js)-lon_rel(MAX(is-1,1),js))/2.0
          ! 
          IF ( lat_rel(1,1) < lat_rel(iml,jml)) THEN
             ! Grid in SN direction
             dls=(lat_rel(is,js)-lat_rel(is,MAX(js-1,1)))/2.0
             dln=(lat_rel(is,MIN(js+1,jml))-lat_rel(is,js))/2.0
          ELSE
             dls=(lat_rel(is,js)-lat_rel(is,MIN(js+1,jml)))/2.0
             dln=(lat_rel(is,MAX(js-1,1))-lat_rel(is,js))/2.0
          ENDIF
          ! Treat the border cases
          IF ( dlw < dle/2. ) dlw = dle
          IF ( dle < dlw/2. ) dle = dlw
          !
          IF ( dls < dln/2. ) dls = dln
          IF ( dln < dls/2. ) dln = dls
          !
          DO i=0,sbd
             ! Put the points in the order : SW, SE, NE, NW
             lons(2+i) = lons(1) - dlw + i*((dlw+dle)/FLOAT(sbd+1))
             lats(2+i) = lats(1) - dls
             !
             lons(2+(sbd+1)+i) = lons(1) + dle
             lats(2+(sbd+1)+i) = lats(1) - dls + i*((dls+dln)/FLOAT(sbd+1))
             !
             lons(2+2*(sbd+1)+i) = lons(1) + dle - i*((dlw+dle)/FLOAT(sbd+1))
             lats(2+2*(sbd+1)+i) = lats(1) + dln
             !
             lons(2+3*(sbd+1)+i) = lons(1) - dlw
             lats(2+3*(sbd+1)+i) = lats(1) + dln - i*((dls+dln)/FLOAT(sbd+1))
          ENDDO
          !
          ! Test that the polygone overlaps with the Bounding box. We might have projections
          ! which are not valid outside of the RCM domain.
          !
          IF ( ((MINVAL(lons) > BBlon(1) .AND. MINVAL(lons) < BBlon(2)) .OR.   &
             &  (MAXVAL(lons) > BBlon(1) .AND. MAXVAL(lons) < BBlon(2))) .AND. & 
             & ((MINVAL(lats) > BBlat(1) .AND. MINVAL(lats) < BBlat(2)) .OR.   &
             &  (MAXVAL(lats) > BBlat(1) .AND. MAXVAL(lats) < BBlat(2)))) THEN
             !
             !
             CALL Grid_Toij (lons, lats, ri, rj)
             !
             !
             ! Find the points of the WRF grid boxes we need to scan to place this box (is, js) 
             ! of the source grid.
             !
             ilow = MAX(FLOOR(MINVAL(ri)),1)
             iup  = MIN(CEILING(MAXVAL(ri)), iim_g)
             jlow = MAX(FLOOR(MINVAL(rj)),1)
             jup  = MIN(CEILING(MAXVAL(rj)), jjm_g)
             !
             !
             IF ( ilow < iup .OR. jlow < jup ) THEN
                !
                ! Find the grid boxes in WRF and compute the overlap area.
                !
                CALL interregxy_findpoints(iim_g, jjm_g, lon_rel(is,js), lon_rel(is,js), is, js, ri, rj, incmax, &
                     &    checkprint, nbov, sourcei, sourcej, overlap)
                
             ENDIF
             !
             !
          ENDIF
       ENDDO
    ENDDO
    !
    ! Gather the land points from the global grid
    !
    areaoverlap(:,:) = zero
    DO in=1,nbpt
       !
       is = ilandindex(in)
       js = jlandindex(in)
       !
       indinc(in,1:nbov(is,js),1) = sourcei(is,js,1:nbov(is,js))
       indinc(in,1:nbov(is,js),2) = sourcej(is,js,1:nbov(is,js))
       areaoverlap(in,1:nbov(is,js)) = overlap(is,js,1:nbov(is,js))
    ENDDO
    !
    ! The model will stop if incmax is not sufficient. So if we are here all is OK.
    !
    ok=.TRUE.
    !
  END SUBROUTINE interregxy_aggr2d
!
! ===================================================================================
!
  SUBROUTINE interregxy_aggrve(nbpt, lalo, neighbours, resolution, contfrac, &
       &                iml, lon_rel, lat_rel, resol_lon, resol_lat, callsign, &
       &                incmax, indinc, areaoverlap, ok)
    !
    ! INPUT
    ! 
    INTEGER(i_std), INTENT(in)   :: nbpt                 ! Number of points for which the data needs to be interpolated
    REAL(r_std), INTENT(in)      :: lalo(nbpt,2)         ! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), INTENT(in)   :: neighbours(nbpt,NbNeighb)   ! Vector of neighbours for each grid point (1=N, 2=E, 3=S, 4=W)
    REAL(r_std), INTENT(in)      :: resolution(nbpt,2)   ! The size in km of each grid-box in X and Y
    REAL(r_std), INTENT(in)      :: contfrac(nbpt)       ! Fraction of land in each grid box.
    INTEGER(i_std), INTENT(in)   :: iml                  ! Size of the source grid
    REAL(r_std), INTENT(in)      :: lon_rel(iml)         ! Longitudes for the source grid
    REAL(r_std), INTENT(in)      :: lat_rel(iml)         ! Latitudes for the source grid
    REAL(r_std), INTENT(in)      :: resol_lon, resol_lat ! Resolution in meters of the source grid
    CHARACTER(LEN=*), INTENT(in) :: callsign             ! Allows to specify which variable is beeing treated
    INTEGER(i_std), INTENT(in)   :: incmax               ! Maximum point of the source grid we can store.
    !
    ! Output
    !
    INTEGER(i_std), INTENT(out)    :: indinc(nbpt,incmax)
    REAL(r_std), INTENT(out)       :: areaoverlap(nbpt,incmax)
    LOGICAL, OPTIONAL, INTENT(out) :: ok                 ! return code
    !
    !
    ! Local Variables
    !
    ! The number of sub-division we wish to have on each side of the rectangle.
    ! sbd=1 means that on each side we take only the middle point.
    !
    INTEGER(i_std), PARAMETER  :: sbd=5
    INTEGER(i_std), PARAMETER  :: idom=1
    !
    INTEGER(i_std)                                :: i, is, js, in
    INTEGER(i_std)                                :: ilow, iup, jlow, jup
    REAL(r_std)                                   :: dle, dlw, dls, dln, coslat
    !
    REAL(r_std), DIMENSION((sbd+1)*4+1)           :: lons, lats, ri, rj
    !
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:)   :: nbov
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:,:) :: sourcei, sourcej
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)    :: overlap
    !
    !
    REAL(r_std), DIMENSION(2)   :: BBlon, BBlat
    REAL(r_std)                 :: half_resol_lon, half_resol_lat, trtmp
    LOGICAL                     :: checkprint=.FALSE.
    !
    WRITE(*,*) "interregxy_aggrve working on ", callsign
    !
    ALLOCATE(sourcei(iim_g, jjm_g, incmax))
    ALLOCATE(sourcej(iim_g, jjm_g, incmax))
    ALLOCATE(nbov(iim_g, jjm_g))
    ALLOCATE(overlap(iim_g, jjm_g, incmax))
    nbov(:,:) = 0
    !
    BBlon(1) = MINVAL(lalo(:,2))
    BBlon(2) = MAXVAL(lalo(:,2))
    BBlat(1) = MINVAL(lalo(:,1))
    BBlat(2) = MAXVAL(lalo(:,1))
    !
    ! Go through the entire source grid and try to place each grid box in a target grid box.
    !
    DO is=1,iml
       !
       ! 1.0 Get the center of the box and their corresponding grid point 
       !
       lons(1) = lon_rel(is)
       lats(1) = lat_rel(is)
       ! Some pre-calculations
       trtmp = 180./(2.*pi)
       half_resol_lon = resol_lon/2.0
       half_resol_lat = resol_lat/2.0
       !
       dln = half_resol_lat*trtmp/R_Earth
       dls = half_resol_lat*trtmp/R_Earth
       !
       DO i=0,sbd
          ! Put the points in the order : SW, SE, NE, NW
          lats(2+i) = lats(1) - dls
          coslat = COS( lats(2+i) * pi/180. )
          dlw = half_resol_lon*trtmp/R_Earth/coslat
          lons(2+i) = lons(1) - dlw + i*((2.*dlw)/FLOAT(sbd+1))
          !
          lats(2+(sbd+1)+i) = lats(1) - dls + i*((dls+dln)/FLOAT(sbd+1))
          coslat = COS( lats(2+(sbd+1)+i) * pi/180. )
          dle = half_resol_lon*trtmp/R_Earth/coslat
          lons(2+(sbd+1)+i) = lons(1) + dle
          !
          lats(2+2*(sbd+1)+i) = lats(1) + dln
          coslat = COS( lats(2+2*(sbd+1)+i) * pi/180. )
          dle = half_resol_lon*trtmp/R_Earth/coslat
          lons(2+2*(sbd+1)+i) = lons(1) + dle - i*((2*dle)/FLOAT(sbd+1))
          !
          lats(2+3*(sbd+1)+i) = lats(1) + dln - i*((dls+dln)/FLOAT(sbd+1))
          coslat = COS( lats(2+3*(sbd+1)+i) * pi/180. )
          dlw = half_resol_lon*trtmp/R_Earth/coslat
          lons(2+3*(sbd+1)+i) = lons(1) - dlw
       ENDDO
       !
       ! Test that the polygone overlaps with the Bounding box. We might have projections
       ! which are not valid outside of the RCM domain.
       !
       IF ( ((MINVAL(lons) > BBlon(1) .AND. MINVAL(lons) < BBlon(2)) .OR.   &
            &  (MAXVAL(lons) > BBlon(1) .AND. MAXVAL(lons) < BBlon(2))) .AND. & 
            & ((MINVAL(lats) > BBlat(1) .AND. MINVAL(lats) < BBlat(2)) .OR.   &
            &  (MAXVAL(lats) > BBlat(1) .AND. MAXVAL(lats) < BBlat(2)))) THEN
          !
          CALL Grid_Toij (lons, lats, ri, rj)
          !
          !
          ! Find the points of the WRF grid boxes we need to scan to place this box (is, js) 
          ! of the source grid.
          !
          ilow = MAX(FLOOR(MINVAL(ri)),1)
          iup  = MIN(CEILING(MAXVAL(ri)), iim_g)
          jlow = MAX(FLOOR(MINVAL(rj)),1)
          jup  = MIN(CEILING(MAXVAL(rj)), jjm_g)
          !
          !
          IF ( ilow < iup .OR. jlow < jup ) THEN
             !
             ! Find the grid boxes in WRF and compute the overlap area.
             !
             CALL interregxy_findpoints(iim_g, jjm_g, lon_rel(is), lon_rel(is), is, 1, ri, rj, incmax, &
                  &    checkprint, nbov, sourcei, sourcej, overlap)
             
          ENDIF
          !
          !
       ENDIF
    ENDDO
    !
    ! Gather the land points from the global grid
    !
    areaoverlap(:,:) = zero
    DO in=1,nbpt
       !
       is = ilandindex(in)
       js = jlandindex(in)
       !
       indinc(in,1:nbov(is,js)) = sourcei(is,js,1:nbov(is,js))
       areaoverlap(in,1:nbov(is,js)) = overlap(is,js,1:nbov(is,js))
    ENDDO
    !
    ! The model will stop if incmax is not sufficient. So if we are here all is OK.
    !
    ok=.TRUE.
    !
  END SUBROUTINE interregxy_aggrve
!
! ===================================================================================
!
  SUBROUTINE interregxy_findpoints(nb_we, nb_sn, lon_cen, lat_cen, is, js, ri, rj, &
       &                           nbpt_max, checkprint, nbpt, sourcei, sourcej, overlap_area)
    !
    !
    !   This Subroutine finds the overlap of grid box ri,rj with the WRF grid.  (ri,rj) are the WRF indexes of the (is,js) point
    !   of the source grid.
    !   Each time we find an overlap of (is,js) with a WRF grid box (ix,jx) we add that in the (sourcei,sourcej) fields.
    !
    !
    ! Arguments
    !
    INTEGER(i_std), INTENT(in)                      :: nb_we, nb_sn
    REAL(r_std), INTENT(in)                         :: lon_cen, lat_cen
    INTEGER(i_std), INTENT(in)                      :: is, js
    REAL(r_std), DIMENSION(:), INTENT(in)           :: ri, rj
    INTEGER(i_std), INTENT(in)                      :: nbpt_max
    LOGICAL, INTENT(in)                             :: checkprint
    INTEGER(i_std), DIMENSION(:,:), INTENT(out)     :: nbpt
    INTEGER(i_std), DIMENSION(:,:,:), INTENT(out)   :: sourcei, sourcej
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)      :: overlap_area
    !
    ! Local
    ! 
    INTEGER(i_std)                                  :: ilow, iup, jlow, jup, nbr
    INTEGER(i_std)                                  :: ix, jx, i

    INTEGER(i_std)                                  :: nbint, nbcross, nbsorted, nbfinal, nbtmpa, nbtmpb
    !
    INTEGER(i_std)                                  :: nbdots=5, dimpoly = 200
    REAL(r_std), DIMENSION(4,2)                     :: poly_wrf
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)        :: poly_in, poly_tmpa, poly_tmpb, poly_int, poly_cross
    !
    REAL(r_std)                                     :: ax, ay, overlaparea
    !
    INTEGER(i_std), PARAMETER                       :: icheck=50, jcheck=46
    !
    !
    !
    nbr = SIZE(ri)
    !
    IF ( dimpoly < 4*(nbr-1) ) THEN
       dimpoly =  4*nbr
    ENDIF
    !
    ALLOCATE(poly_in(nbr-1,2))
    ALLOCATE(poly_tmpa(dimpoly,2))
    ALLOCATE(poly_tmpb(dimpoly,2))
    ALLOCATE(poly_int(dimpoly,2))
    ALLOCATE(poly_cross(dimpoly,2))
    !
    ! Scan all Determine the range of indicis of the WRF grid we need to scan 
    ! within the box of the source grid.
    !
    ilow = MAX(FLOOR(MINVAL(ri)),1)
    iup  = MIN(CEILING(MAXVAL(ri)), nb_we)
    jlow = MAX(FLOOR(MINVAL(rj)),1)
    jup  = MIN(CEILING(MAXVAL(rj)), nb_sn)
    !
    !
    poly_wrf(:,1) = (/ -0.5, 0.5, 0.5, -0.5 /)
    poly_wrf(:,2) = (/ -0.5, -0.5, 0.5, 0.5 /)
    !
    DO ix = ilow,iup
       DO jx = jlow,jup
          !
          ! Bring the points of the source grid into the +-0.5 range of the WRF grid ...
          ! and remove the center point.
          !
          poly_in(1:nbr-1,1) = ri(2:nbr) - ix
          poly_in(1:nbr-1,2) = rj(2:nbr) - jx
          !
          IF ( checkprint .AND. ix == icheck .AND. jx == jcheck ) THEN
             WRITE(*,*) "X = ", poly_in(1:nbr-1,1)
             WRITE(*,*) "Y = ", poly_in(1:nbr-1,2)
          ENDIF
          !
          CALL polygones_extend(nbr-1, poly_in, nbdots, nbtmpa, poly_tmpa)
          CALL polygones_extend(4, poly_wrf, nbdots, nbtmpb, poly_tmpb)
          !
          IF ( checkprint .AND. ix == icheck .AND. jx == jcheck ) THEN
             WRITE(*,*) "X_extend = ", poly_tmpa(1:nbtmpa,1)
             WRITE(*,*) "Y_extend = ", poly_tmpa(1:nbtmpa,2)
          ENDIF
          !
          CALL polygones_intersection(nbtmpa, poly_tmpa, nbtmpb, poly_tmpb, nbint, poly_int) 
          CALL polygones_crossing(nbr-1, poly_in,  4, poly_wrf, nbcross, poly_cross)
          !
          IF ( checkprint .AND. ix == icheck .AND. jx == jcheck ) THEN
              WRITE(*,*) "nbint = nbcross = ", nbint, nbcross
           ENDIF
          !
          IF ( nbint+nbcross > dimpoly ) THEN
             WRITE(*,*) "The intersection of polygones is larger than the memory foressen : ",  nbint+nbcross, dimpoly
             STOP "interregxy_compoverlap"
          ENDIF
          !
          !
          IF ( nbint+nbcross > 2 ) THEN

             poly_tmpa(1:nbint,:) = poly_int(1:nbint,:)
             poly_tmpa(nbint+1:nbint+nbcross,:) =  poly_cross(1:nbcross,:)
             
             IF ( checkprint .AND. ix == icheck .AND. jx == jcheck ) THEN
                WRITE(*,*) "X_overlap = ", poly_tmpa(1:(nbint+nbcross),1)
                WRITE(*,*) "Y_overlap = ", poly_tmpa(1:(nbint+nbcross),2)
             ENDIF

             CALL polygones_cleanup(nbint+nbcross, poly_tmpa, nbfinal, poly_tmpb)
             IF ( nbint+nbcross < 3 ) THEN
                WRITE(*,*) "poly_tmpa X : ", poly_tmpa(1:(nbint+nbcross),1)
                WRITE(*,*) "poly_tmpa Y : ", poly_tmpa(1:(nbint+nbcross),2)
                WRITE(*,*) "Cleaning goes too far : ", nbint+nbcross, " to ", nbfinal
                WRITE(*,*) "poly_tmpb X : ", poly_tmpb(1:nbfinal,1)
                WRITE(*,*) "poly_tmpb Y : ", poly_tmpb(1:nbfinal,2)
             ENDIF

             IF ( checkprint .AND. ix == icheck .AND. jx == jcheck ) THEN
                WRITE(*,*) "X_final = ", poly_tmpb(1:nbfinal,1)
                WRITE(*,*) "Y_final = ", poly_tmpb(1:nbfinal,2)
             ENDIF
             IF ( nbfinal > 2 ) THEN
                CALL polygones_convexhull(nbfinal, poly_tmpb, nbsorted, poly_tmpa)
                IF ( checkprint .AND. ix == icheck .AND. jx == jcheck ) THEN
                   WRITE(*,*) "X_sorted = ", poly_tmpa(1:nbsorted,1)
                   WRITE(*,*) "Y_sorted = ", poly_tmpa(1:nbsorted,2)
                ENDIF
                CALL polygones_area(nbsorted, poly_tmpa, dxWRF(ix,jx), dyWRF(ix,jx), overlaparea)
             ELSE
                overlaparea = 0.0
             ENDIF

             IF ( checkprint .AND. ix == icheck .AND. jx == jcheck ) THEN
                WRITE(*,*) ix, jx, "overlaparea = ", overlaparea
                WRITE(*,*) "============================================================================="
             ENDIF

             IF ( overlaparea > 0.0 ) THEN
                nbpt(ix,jx) = nbpt(ix,jx)+1
                IF ( nbpt(ix,jx) > nbpt_max ) THEN
                   WRITE(*,*) "ERROR in interregxy_findpoints. "
                   WRITE(*,*) "Consider your algorithm to estimate nbpt_max. "
                   WRITE(*,*) "You should change it as it does not give enough points."
                   WRITE(*,*) "nbpt(ix,jx) = ", nbpt(ix,jx)
                   WRITE(*,*) "nbpt_max = ", nbpt_max
                   STOP "interregxy_findpoints"
                ENDIF
                sourcei(ix,jx,nbpt(ix,jx)) = is
                sourcej(ix,jx,nbpt(ix,jx)) = js
                overlap_area(ix,jx,nbpt(ix,jx)) =  overlaparea
                !
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    !
    DEALLOCATE(poly_in, poly_tmpa, poly_tmpb, poly_int, poly_cross)
    !
  END SUBROUTINE interregxy_findpoints

END MODULE interregxy
