!  ==============================================================================================================================\n
!  MODULE 	: This module provides a few basic operations for convex polygons on the sphere and in longitude
!                 latitude coordinates. These operations are based on the haversine equations for great circle arcs.
!
!                 Definition of polygons : the basic assumption is that they describe grid boxes and are thus 
!                 provided with more points than strictly necessary. Each polygon starts at index 1 on a vertex and
!                 alternates mid-points of segments and vertices. The mid-point of segments are kept as they are
!                 useful elements for the operations on the grid boxes.
!
!                 The module provides the following subroutine and functions :
!                 haversine_reglatlontoploy : Lists the polygons and all their attributes (centre point, 
!                                             neighbouring grids, indices in i and j on the original grid) for 
!                                             regular lon/lat grid.
!                 haversine_regxytoploy : As above but for grid boxes on the sphere projected onto a regular
!                                         X/Y grid.
!                 haversine_singlepointpoly : Computes the polygone around a given point with a known area.
!                 haversine_polyheadings : Compute the heading out of the polygon for each vertex and mid-point of
!                                          segment.
!                 haversine_polysort : Sort the polygons so that all points are in clockwise order.
!                 haversine_polyseglen : Compute the length of all segments of the polygon using great circles.
!                 haversine_polyseglen : Compute the length of all segments of the polygon on a regular lat lon grid.
!                 haversine_clockwise : Get the indices which sort a polygon in a clockwise order starting from an
!                                       initial vertex.
!                 haversine_polyarea : Computes the area covered by a polygon on the sphere.
!                 haversine_laloarea : Compute the area for a regular lat lon grid.
!                 haversine_xyarea : Compute the area for the special case where the grid box size in x and y are already 
!                                    given by the projection.
!                 haversine_heading : Initial heading between start point and end point along a great circle.
!                 haversine_distance : Compute the distance between 2 points along the great circle.
!                 haversine_radialdis : Compute the coordinates found in a given heading and distance.
!                 haversine_dtor : Degrees to radians
!                 haversine_rtod : Radians to degrees
!
!  CONTACT      : jan.polcher@lmd.jussieu.fr
!
!  LICENCE      : IPSL (2016)
!  This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        
!! 
!! RECENT CHANGE(S): None 
!! 
!! REFERENCE(S) : None
!! 
!_ ================================================================================================================================
MODULE haversine

  USE defprec
  USE constantes_var
  USE module_llxy

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: haversine_reglatlontoploy, haversine_regxytoploy, haversine_singlepointploy, &
       &    haversine_polyheadings, haversine_polysort, &
       &    haversine_polyseglen, haversine_laloseglen, haversine_clockwise, haversine_polyarea, &
       &    haversine_laloarea, haversine_xyarea

CONTAINS
!!  =============================================================================================================================
!! SUBROUTINE: haversine_reglatlontoploy   
!!
!>\BRIEF       Lists the polygons and all their attributes (centre point, 
!              neighbouring grids, indices in i and j on the original grid) for 
!              regular lon/lat grid.
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================
  SUBROUTINE haversine_reglatlontoploy(iim, jjm, lon, lat, nbpt, index_loc, global, &
       &                              nbseg, lonpoly, latpoly, center, neighb_loc, iorig, jorig)
    !
    ! This subroutine constructs a series of polygons out of the grid boxes which are listed
    ! in the index_loc array. The polygons are ordered according to the indexing space and independently
    ! of the geographical coordinates.
    !
    ! 0 interface
    !
    IMPLICIT NONE
    !
    ! 0.1 input  !
    ! Size of cartesian grid
    INTEGER(i_std), INTENT(in)                                 :: iim, jjm
    ! Longitudes on cartesian grid
    REAL(r_std), DIMENSION(iim,jjm), INTENT(in)                :: lon
    ! Latitudes on cartesian grid
    REAL(r_std), DIMENSION(iim,jjm), INTENT(in)                :: lat
    ! Size of the gathered domain
    INTEGER(i_std), INTENT(in)                                 :: nbpt
    ! Index of land point on 2D map (in local position)
    INTEGER(i_std), DIMENSION(nbpt), INTENT(in)                :: index_loc
    ! Is it a global grid ?
    LOGICAL, INTENT(in)                                        :: global
    ! Number of segments
    INTEGER(i_std), INTENT(in)                                 :: nbseg
    !
    ! 0.2 Ouput
    !
    REAL(r_std), DIMENSION(nbpt,nbseg*2), INTENT(out)          :: lonpoly
    REAL(r_std), DIMENSION(nbpt,nbseg*2), INTENT(out)          :: latpoly
    REAL(r_std), DIMENSION(nbpt,2), INTENT(out)                :: center
    INTEGER(i_std), DIMENSION(nbpt,nbseg*2), INTENT(out)       :: neighb_loc
    INTEGER(i_std), DIMENSION(nbpt), INTENT(out)               :: iorig, jorig
    !
    !
    ! 0.3 Local variables
    !
    INTEGER(i_std), DIMENSION(iim,jjm) :: correspondance
    INTEGER(i_std)                     :: i, ip, jp
    INTEGER(i_std)                     :: ipm1, ipp1, jpm1, jpp1
    REAL(r_std)                        :: dlonm1, dlonp1, dlatm1, dlatp1
    !
    !
    correspondance(:,:) = -1
    DO i = 1, nbpt      
       !
       ! 1 find numbers of the latitude and longitude of each point
       !
       ! index of latitude
       jp = INT( (index_loc(i)-1) /iim ) + 1
       ! index of longitude
       ip = index_loc(i) - ( jp-1 ) * iim
       !
       !correspondance(ip,jp) = kindex(i)
       !
       correspondance(ip,jp) = i
       iorig(i)=ip
       jorig(i)=jp
       !
    ENDDO
    !
    !
    ! Go through all the points and build the polygone which defines the grid area.
    ! This polygone will include the mid-point of each segment so that 
    ! we can later compute the direction of the normal to the segment.
    !
    neighb_loc(:,:) = -1
    !
    DO i = 1, nbpt
       !
       ip = iorig(i)
       jp = jorig(i)
       !
       ipm1 = ip-1
       ipp1 = ip+1
       jpm1 = jp-1
       jpp1 = jp+1
       !
       ! Compute the longitude increments
       !
       IF ( ipp1 <= iim ) THEN
          dlonp1 = (lon(ipp1,jp)-lon(ip,jp))/2.0
       ELSE IF ( ipm1 > 0 ) THEN
          dlonp1 = (lon(ip,jp)-lon(ipm1,jp))/2.0
          IF ( global ) ipp1=1
       ELSE
          dlonp1 = undef_sechiba
       ENDIF
       IF ( ipm1 > 0 ) THEN
          dlonm1 = (lon(ip,jp)-lon(ipm1,jp))/2.0
       ELSE IF ( ipp1 <= iim ) THEN
          dlonm1 = (lon(ipp1,jp)-lon(ip,jp))/2.0
          IF ( global ) ipm1=iim
       ELSE
          dlonm1 = undef_sechiba
       ENDIF
       !
       ! Verify that we have at least one valid longitude increment. Else we do not have enough
       ! points in the grid to estimate the position of the vertices in longitude.
       !
       IF ( dlonp1 >= undef_sechiba-1 ) dlonp1 = dlonm1
       IF ( dlonm1 >= undef_sechiba-1 ) dlonm1 = dlonp1
       IF ( dlonp1 >= undef_sechiba-1 .AND. dlonm1 >= undef_sechiba-1 ) THEN
          CALL ipslerr(3, "haversine_reglatlontoploy", "There are not enogh point in longitude", &
               &       "to estimate the bounds of the polygone of the grid box.",&
               &       "Please choose a larger grid.")
       ENDIF
       !
       ! Compute the latitude increments
       !
       IF ( jpp1 <= jjm ) THEN
          dlatp1 = (lat(ip,jpp1)-lat(ip,jp))/2.0
       ELSE IF ( jpm1 > 0 ) THEN
          dlatp1 = (lat(ip,jp)-lat(ip,jpm1))/2.0
       ELSE
          dlatp1 = undef_sechiba
       ENDIF
       IF ( jpm1 > 0 ) THEN
          dlatm1 = (lat(ip,jp)-lat(ip,jpm1))/2.0
       ELSE IF ( jpp1 <= jjm ) THEN
          dlatm1 = (lat(ip,jpp1)-lat(ip,jp))/2.0
       ELSE
          dlatm1 = undef_sechiba
       ENDIF
       !
       ! Verify that we have at least one valid latitude increment. Else we do not have enough
       ! points in the grid to estimate the position of the vertices in latitude.
       !
       IF ( dlatp1 >= undef_sechiba-1 ) dlatp1 = dlatm1
       IF ( dlatm1 >= undef_sechiba-1 ) dlatm1 = dlatp1
       IF ( dlatp1 >= undef_sechiba-1 .AND. dlatm1 >= undef_sechiba-1 ) THEN
          CALL ipslerr(3, "haversine_reglatlontoploy", "There are not enogh point in latitude", &
               &       "to estimate the bounds of the polygone of the grid box.",&
               &       "Please choose a larger grid.")
       ENDIF
       !
       ! The longitude of all elements of the polygone
       !
       lonpoly(i,1) = lon(ip,jp)-dlonm1
       lonpoly(i,2) = lon(ip,jp)
       lonpoly(i,3) = lon(ip,jp)+dlonp1
       lonpoly(i,4) = lon(ip,jp)+dlonp1
       lonpoly(i,5) = lon(ip,jp)+dlonp1
       lonpoly(i,6) = lon(ip,jp)
       lonpoly(i,7) = lon(ip,jp)-dlonm1
       lonpoly(i,8) = lon(ip,jp)-dlonm1
       !
       ! The longitude of all elements of the polygone
       !
       latpoly(i,1) = lat(ip,jp)-dlatp1
       latpoly(i,2) = lat(ip,jp)-dlatp1
       latpoly(i,3) = lat(ip,jp)-dlatp1
       latpoly(i,4) = lat(ip,jp)
       latpoly(i,5) = lat(ip,jp)+dlatm1
       latpoly(i,6) = lat(ip,jp)+dlatm1
       latpoly(i,7) = lat(ip,jp)+dlatm1
       latpoly(i,8) = lat(ip,jp)
       !
       ! Keep the center of the gridbox
       !
       center(i,1) = lon(ip,jp)
       center(i,2) = lat(ip,jp)
       !
       ! Get the neighbours when they exist in the list of land points
       ! There are no neighbours over the North or South poles.
       !
       IF ( ipm1 > 0 .AND. jpm1 > 0 ) neighb_loc(i,1) = correspondance(ipm1,jpm1)
       IF ( jpm1 > 0 ) neighb_loc(i,2) = correspondance(ip,jpm1)
       IF ( ipp1 <= iim .AND. jpm1 > 0 ) neighb_loc(i,3) = correspondance(ipp1,jpm1)
       IF ( ipp1 <= iim ) neighb_loc(i,4) = correspondance(ipp1,jp)
       IF ( ipp1 <= iim .AND. jpp1 <= jjm ) neighb_loc(i,5) = correspondance(ipp1,jpp1)
       IF ( jpp1 <= jjm ) neighb_loc(i,6) = correspondance(ip,jpp1)
       IF ( ipm1 > 0 .AND. jpp1 <= jjm ) neighb_loc(i,7) = correspondance(ipm1,jpp1)
       IF ( ipm1 > 0 ) neighb_loc(i,8) = correspondance(ipm1,jp)
       !
    ENDDO
    !
  END SUBROUTINE haversine_reglatlontoploy
!!
!!  =============================================================================================================================
!! SUBROUTINE: haversine_regxytoploy    
!!
!>\BRIEF       Same as haversine_reglatlontoploy but for grid boxes on the sphere projected onto a regular
!              X/Y grid.
!!
!! DESCRIPTION: Keep in mind that in these projections the straight line assumed between 2 points is not always
!!              the great circle. Thus the distance computed by the haversine formula might deviate a little.
!!
!! \n
!_ ==============================================================================================================================
  SUBROUTINE haversine_regxytoploy(iim, jjm, lon, lat, nbpt, index_loc, proj, &
       &                              nbseg, lonpoly, latpoly, center, neighb_loc, iorig, jorig)
    !
    ! This subroutine constructs a series of polygons out of the grid boxes which are listed
    ! in the index array. This version will go directly from the indexing space to the coordinate as we know that 
    ! we are dealing with a projection of the sphere to the plane where the regular grid is created.
    ! The polygons are ordered according to the indexing space and independently
    ! of the geographical coordinates.
    !
    ! 0 interface
    !
    IMPLICIT NONE
    !
    ! 0.1 input  !
    ! Size of cartesian grid
    INTEGER(i_std), INTENT(in)                                 :: iim, jjm
    ! Longitudes on cartesian grid
    REAL(r_std), DIMENSION(iim,jjm), INTENT(in)                :: lon
    ! Latitudes on cartesian grid
    REAL(r_std), DIMENSION(iim,jjm), INTENT(in)                :: lat
    ! Size of the gathered domain
    INTEGER(i_std), INTENT(in)                                 :: nbpt
    ! Index of land point on 2D map (in local position)
    INTEGER(i_std), DIMENSION(nbpt), INTENT(in)                :: index_loc
    ! Projection ID
    type (proj_info), DIMENSION(1)                             :: proj
    ! Number of segments
    INTEGER(i_std), INTENT(in)                                 :: nbseg
    !
    ! 0.2 Ouput
    !
    REAL(r_std), DIMENSION(nbpt,nbseg*2), INTENT(out)          :: lonpoly
    REAL(r_std), DIMENSION(nbpt,nbseg*2), INTENT(out)          :: latpoly
    REAL(r_std), DIMENSION(nbpt,2), INTENT(out)                :: center
    INTEGER(i_std), DIMENSION(nbpt,nbseg*2), INTENT(out)       :: neighb_loc
    INTEGER(i_std), DIMENSION(nbpt), INTENT(out)               :: iorig, jorig
    !
    !
    ! 0.3 Local variables
    !
    INTEGER(i_std), DIMENSION(iim,jjm) :: correspondance
    INTEGER(i_std)                     :: i, ip, jp
    INTEGER(i_std)                     :: ipm1, ipp1, jpm1, jpp1
    REAL(r_std)                        :: dlonm1, dlonp1, dlatm1, dlatp1
    !
    !
    !
    correspondance(:,:) = -1
    DO i = 1, nbpt      
       !
       ! 1 find numbers of the latitude and longitude of each point
       !
       ! index of latitude
       jp = INT( (index_loc(i)-1) /iim ) + 1
       ! index of longitude
       ip = index_loc(i) - ( jp-1 ) * iim
       !
       !correspondance(ip,jp) = kindex(i)
       !
       correspondance(ip,jp) = i
       iorig(i)=ip
       jorig(i)=jp
       !
    ENDDO
    !
    ! Go through all the points and build the polygone which defines the grid area.
    ! This polygone will include the mid-point of each segment so that 
    ! we can later compute the direction of the normal to the segment.
    !
    neighb_loc(:,:) = -1
    !
    DO i = 1, nbpt
       ! index of latitude
       jp = INT( (index_loc(i)-1) /iim ) + 1
       ! index of longitude
       ip = index_loc(i) - ( jp-1 ) * iim
       !
       ipm1 = ip-1
       ipp1 = ip+1
       jpm1 = jp-1
       jpp1 = jp+1
       !
       ! Get the longitude and latitude throug the projection
       ! The range of possible values for projection depends on the module
       ! which defines these projections. For module_llxy the range is 0-203.
       !
       IF ( proj(1)%code > 0 .AND. proj(1)%code < 203 ) THEN
          CALL ij_to_latlon(proj(1), ip-0.5, jp-0.5, latpoly(i,1), lonpoly(i,1))
          CALL ij_to_latlon(proj(1), ip+0.0, jp-0.5, latpoly(i,2), lonpoly(i,2))
          CALL ij_to_latlon(proj(1), ip+0.5, jp-0.5, latpoly(i,3), lonpoly(i,3))
          CALL ij_to_latlon(proj(1), ip+0.5, jp+0.0, latpoly(i,4), lonpoly(i,4))
          CALL ij_to_latlon(proj(1), ip+0.5, jp+0.5, latpoly(i,5), lonpoly(i,5))
          CALL ij_to_latlon(proj(1), ip+0.0, jp+0.5, latpoly(i,6), lonpoly(i,6))
          CALL ij_to_latlon(proj(1), ip-0.5, jp+0.5, latpoly(i,7), lonpoly(i,7))
          CALL ij_to_latlon(proj(1), ip-0.5, jp+0.0, latpoly(i,8), lonpoly(i,8))
       ELSE
          CALL ipslerr(3, "haversine_regxytoploy", "Unknown projection code", &
               &       "Check proj(1)%code","")
       ENDIF
       !
       ! Keep the center of the gridbox
       !
       center(i,1) = lon(ip,jp)
       center(i,2) = lat(ip,jp)
       !
       ! Get the neighbours when they exist in the list of land points
       ! There are no neighbours over the North or South poles.
       !
       IF ( ipm1 > 0 .AND. jpm1 > 0 ) neighb_loc(i,1) = correspondance(ipm1,jpm1)
       IF ( jpm1 > 0 ) neighb_loc(i,2) = correspondance(ip,jpm1)
       IF ( ipp1 <= iim .AND. jpm1 > 0 ) neighb_loc(i,3) = correspondance(ipp1,jpm1)
       IF ( ipp1 <= iim ) neighb_loc(i,4) = correspondance(ipp1,jp)
       IF ( ipp1 <= iim .AND. jpp1 <= jjm ) neighb_loc(i,5) = correspondance(ipp1,jpp1)
       IF ( jpp1 <= jjm ) neighb_loc(i,6) = correspondance(ip,jpp1)
       IF ( ipm1 > 0 .AND. jpp1 <= jjm ) neighb_loc(i,7) = correspondance(ipm1,jpp1)
       IF ( ipm1 > 0 ) neighb_loc(i,8) = correspondance(ipm1,jp)
       !
    ENDDO
    !
  END SUBROUTINE haversine_regxytoploy
!!  =============================================================================================================================
!! SUBROUTINE: haversine_singlepointploy
!!
!>\BRIEF       Lists the polygons and all their attributes (centre point, 
!              neighbouring grids, indices in i and j on the original grid) for 
!              a single point. A regular lon/lat grid is assumed.
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================
  SUBROUTINE haversine_singlepointploy(iim, jjm, lon, lat, nbpt, index_loc, global, &
       &                              nbseg, lonpoly, latpoly, center, neighb_loc, iorig, jorig)
    !
    ! This subroutine constructs a series of polygons out of the grid boxe which is provided.
    ! The polygon is ordered according to the indexing space and independently
    ! of the geographical coordinates.
    ! It uses the same interface as haversine_reglatlontoploy but can only used with iim=jjm=1
    !
    ! 0 interface
    !
    IMPLICIT NONE
    !
    ! 0.1 input  !
    ! Size of cartesian grid
    INTEGER(i_std), INTENT(in)                                 :: iim, jjm
    ! Longitudes on cartesian grid
    REAL(r_std), DIMENSION(iim,jjm), INTENT(in)                :: lon
    ! Latitudes on cartesian grid
    REAL(r_std), DIMENSION(iim,jjm), INTENT(in)                :: lat
    ! Size of the gathered domain
    INTEGER(i_std), INTENT(in)                                 :: nbpt
    ! Index of land point on 2D map (in local position)
    INTEGER(i_std), DIMENSION(nbpt), INTENT(in)                :: index_loc
    ! Is it a global grid ?
    LOGICAL, INTENT(in)                                        :: global
    ! Number of segments
    INTEGER(i_std), INTENT(in)                                 :: nbseg
    !
    ! 0.2 Ouput
    !
    REAL(r_std), DIMENSION(nbpt,nbseg*2), INTENT(out)          :: lonpoly
    REAL(r_std), DIMENSION(nbpt,nbseg*2), INTENT(out)          :: latpoly
    REAL(r_std), DIMENSION(nbpt,2), INTENT(out)                :: center
    INTEGER(i_std), DIMENSION(nbpt,nbseg*2), INTENT(out)       :: neighb_loc
    INTEGER(i_std), DIMENSION(nbpt), INTENT(out)               :: iorig, jorig
    !
    !
    ! 0.3 Local variables
    !
    REAL(r_std)                        :: area, rlon, rlat
    REAL(r_std), DIMENSION(2)          :: coord
    !
    IF ( iim .NE. 1 .AND. jjm .NE. 1 ) THEN
       CALL ipslerr(3, "haversine_singlepointploy", "Can only be used if iim=jjm=1", &
            &       "Please ensure this routine is called in the","right conditions")
    ENDIF
    !
    iorig(1)=1
    jorig(1)=1
    !
    area = 111111.0*111111.0
    rlon= SQRT(area)/2.0
    rlat=rlon
    WRITE(*,*) "Area : ", area/1.0e6
    !
    ! Set all the variables defining the polygone of the specified area
    !
    neighb_loc(:,:) = -1
    !
    ! Northern point
    coord = haversine_radialdis(lon(1,1), lat(1,1), 0.0, rlon)
    lonpoly(1,2) = coord(1)
    latpoly(1,2) = coord(2)
    ! Eastern point
    coord = haversine_radialdis(lon(1,1), lat(1,1), 90.0, rlon)
    lonpoly(1,4) = coord(1)
    latpoly(1,4) = coord(2)
    ! Souther point
    coord = haversine_radialdis(lon(1,1), lat(1,1), 180.0, rlon)
    lonpoly(1,6) = coord(1)
    latpoly(1,6) = coord(2)
    ! Souther point
    coord = haversine_radialdis(lon(1,1), lat(1,1), 270.0, rlon)
    lonpoly(1,8) = coord(1)
    latpoly(1,8) = coord(2)
    !
    ! Doing the corners
    !
    ! North West
    lonpoly(1,1) = lonpoly(1,8)
    latpoly(1,1) = latpoly(1,2)
    ! North East
    lonpoly(1,3) = lonpoly(1,4)
    latpoly(1,3) = latpoly(1,2)
    ! South East
    lonpoly(1,5) = lonpoly(1,4)
    latpoly(1,5) = latpoly(1,6)
    ! South West
    lonpoly(1,7) = lonpoly(1,8)
    latpoly(1,7) = latpoly(1,6)
    !
    center(1,1) = lon(1,1)
    center(1,2) = lat(1,1)
    !
  END SUBROUTINE haversine_singlepointploy
!!
!!  =============================================================================================================================
!! SUBROUTINE: haversine_polyheadings   
!!
!>\BRIEF       Compute the heading out of the polygon for each vertex and mid-point of
!              segment.
!!
!! DESCRIPTION: This heading is computed by using the great circle between the centre of the polygon and
!!              the point on the boundary considered. The direction is the one facing outwards from the polygon.
!!
!! \n
!_ ==============================================================================================================================
!!
!!
  SUBROUTINE haversine_polyheadings(nbpt, nbseg, lonpoly, latpoly, center, headings)
    !
    ! 0.1 Input variables 
    ! Size of the gathered domain
    INTEGER(i_std), INTENT(in)                                 :: nbpt
    ! Number of segments
    INTEGER(i_std), INTENT(in)                                 :: nbseg
    !
    REAL(r_std), DIMENSION(nbpt,nbseg*2), INTENT(in)           :: lonpoly
    REAL(r_std), DIMENSION(nbpt,nbseg*2), INTENT(in)           :: latpoly
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in)                 :: center
    !
    ! 0.2 Output variables
    REAL(r_std), DIMENSION(nbpt,nbseg*2), INTENT(out)          :: headings
    !
    ! 0.3 Local variables
    !
    INTEGER(i_std) :: i, ns
    !
    ! We compute for each vertex of our polygon (actual vertices and mid-points) the direction to the
    ! centre. The we add 180. to get the opposite direction.
    !
    DO i=1,nbpt
       DO ns=1,nbseg*2
          headings(i,ns) = MOD(haversine_heading(lonpoly(i,ns), latpoly(i,ns), center(i,1), center(i,2))+180.0, 360.0)
       ENDDO
    ENDDO
    !
  END SUBROUTINE haversine_polyheadings
!!
!!  =============================================================================================================================
!! SUBROUTINE: haversine_polysort     
!!
!>\BRIEF Sort the polygons so that all points are in clockwise order.	  
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================
!!
  SUBROUTINE haversine_polysort(nbpt, nbseg, lonpoly, latpoly, headings, neighb)
    !
    ! This subroutine is foreseen for polygones which start with a vertex and then alternate
    ! with mid-points of the segments. The heading at a vertex is the direction (along the 
    ! great circle) between the center of the polygone and this vertex. At a segment mid-point
    ! the direction is the normal facing outward.
    !
    ! 0.1 Input Variables
    ! Size of the gathered domain
    INTEGER(i_std), INTENT(in)                                 :: nbpt
    ! Number of segments
    INTEGER(i_std), INTENT(in)                                 :: nbseg
    !
    REAL(r_std), DIMENSION(nbpt,nbseg*2), INTENT(inout)        :: lonpoly
    REAL(r_std), DIMENSION(nbpt,nbseg*2), INTENT(inout)        :: latpoly
    REAL(r_std), DIMENSION(nbpt,nbseg*2), INTENT(inout)        :: headings
    INTEGER(i_std), DIMENSION(nbpt,nbseg*2), INTENT(inout)     :: neighb
    !
    ! 0.2 Local variables
    !
    INTEGER(i_std) :: i, ic, startindex
    INTEGER(i_std), DIMENSION(nbseg*2)     :: reindex
    REAL(r_std)                            :: starthead
    REAL(r_std), DIMENSION(nbseg*2)        :: lon_loc
    REAL(r_std), DIMENSION(nbseg*2)        :: lat_loc
    REAL(r_std), DIMENSION(nbseg*2)        :: head_loc
    REAL(r_std), DIMENSION(nbseg*2)        :: negb_loc
    !
    DO i=1,nbpt
       head_loc(:) = headings(i,:)
       lon_loc(:) = lonpoly(i,:)
       lat_loc(:) = latpoly(i,:)
       negb_loc(:) = neighb(i,:)
       !
       ! The first vertice of our polygone needs to be the heading closest to 
       ! North (0 degree). No difference is done between vertices and mid-points
       ! of segments.
       !
       starthead=0.0
       !
       CALL haversine_clockwise(nbseg, head_loc, starthead, reindex)
       !
       !
       headings(i,:) = head_loc(reindex(:))
       lonpoly(i,:) = lon_loc(reindex(:))
       latpoly(i,:) = lat_loc(reindex(:))
       neighb(i,:) = negb_loc(reindex(:))
       !
    ENDDO
    !
  END SUBROUTINE haversine_polysort
!!
!!  =============================================================================================================================
!! SUBROUTINE: haversine_polyseglen   
!!
!>\BRIEF       Compute the length of all segments of the polygon using the great circle.
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================
!!
  SUBROUTINE haversine_polyseglen(nbpt, nbseg, lonpoly, latpoly, seglength)
    !
    ! Computes the segment length for each of the polygones. These are
    ! polygones with middle points given for each segment. This we need
    ! to take only every other point.
    !
    !
    ! 0.1 Input Variables
    ! Size of the gathered domain
    INTEGER(i_std), INTENT(in)                                 :: nbpt
    ! Number of segments
    INTEGER(i_std), INTENT(in)                                 :: nbseg
    !
    REAL(r_std), DIMENSION(nbpt,nbseg*2), INTENT(in)           :: lonpoly
    REAL(r_std), DIMENSION(nbpt,nbseg*2), INTENT(in)           :: latpoly
    REAL(r_std), DIMENSION(nbpt,nbseg), INTENT(out)            :: seglength
    !
    ! 0.2 Local variables
    !
    INTEGER(i_std) :: i, iv, istart, iend, iseg, ioff
    REAL(r_std)    :: slpm1, slpp1
    !
    DO i=1,nbpt
       iseg = 0
       !
       ! Find if the first element is a vertex or mid-point of segment.
       ! Get the 2 headings out of the start point. If these headings are larger than 135 degrees
       ! then probably this point is a segment mid-point.
       !
       slpm1 = haversine_heading(lonpoly(i,1), latpoly(i,1), lonpoly(i,nbseg*2), latpoly(i,nbseg*2))
       slpp1 = haversine_heading(lonpoly(i,1), latpoly(i,1), lonpoly(i,2), latpoly(i,2))
       !
       IF ( ABS(MOD(slpp1-slpm1, 360.0)) > 135.0 ) THEN
          ! The polygon starts with a segment mid-point
          ioff = -1
       ELSE
          ioff = 0
       ENDIF
       !
       DO iv=1,nbseg*2,2
          !
          istart=MODULO((iv-1)+ioff,nbseg*2)+1
          iend=MODULO((iv-1)+ioff+2,nbseg*2)+1
          iseg = iseg + 1
          !
          seglength(i,iseg) = haversine_distance(lonpoly(i,istart), latpoly(i,istart), &
               &                                 lonpoly(i,iend), latpoly(i,iend))
       ENDDO
    ENDDO
    !
  END SUBROUTINE haversine_polyseglen
!!
!!  =============================================================================================================================
!! SUBROUTINE: haversine_laloseglen   
!!
!>\BRIEF       Compute the length of all segments of the polygon when on a regular Lat Lon grid.
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================
!!
  SUBROUTINE haversine_laloseglen(nbpt, nbseg, lonpoly, latpoly, seglength)
    !
    ! Computes the segment length for each of the polygones. These are
    ! polygones with middle points given for each segment. This we need
    ! to take only every other point.
    !
    !
    ! 0.1 Input Variables
    ! Size of the gathered domain
    INTEGER(i_std), INTENT(in)                                 :: nbpt
    ! Number of segments
    INTEGER(i_std), INTENT(in)                                 :: nbseg
    !
    REAL(r_std), DIMENSION(nbpt,nbseg*2), INTENT(in)           :: lonpoly
    REAL(r_std), DIMENSION(nbpt,nbseg*2), INTENT(in)           :: latpoly
    REAL(r_std), DIMENSION(nbpt,nbseg), INTENT(out)            :: seglength
    !
    ! 0.2 Local variables
    !
    INTEGER(i_std) :: i, iv, istart, iend, ioff, iseg
    REAL(r_std)    :: coslat, slpm1, slpp1
    !
    DO i=1,nbpt
       iseg = 0
       !
       ! Find if the first element is a vertex or mid-point of segment.
       ! Get the 2 headings out of the start point. If these headings are larger than 135 degrees
       ! then probably this point is a segment mid-point.
       !
       slpm1 = haversine_heading(lonpoly(i,1), latpoly(i,1), lonpoly(i,nbseg*2), latpoly(i,nbseg*2))
       slpp1 = haversine_heading(lonpoly(i,1), latpoly(i,1), lonpoly(i,2), latpoly(i,2))
       !
       IF ( ABS(MOD(slpp1-slpm1, 360.0)) > 135.0 ) THEN
          ! The polygon starts with a segment mid-point
          ioff = -1
       ELSE
          ioff = 0
       ENDIF
       !
       DO iv=1,nbseg*2,2
          !
          istart=MODULO((iv-1)+ioff,nbseg*2)+1
          iend=MODULO((iv-1)+ioff+2,nbseg*2)+1
          iseg = iseg + 1
          !
          !
          IF ( ABS(lonpoly(i,istart)-lonpoly(i,iend)) < EPSILON(lonpoly) ) THEN
             !
             ! Distance along a meridian
             !
             seglength(i,iseg) =  ABS(latpoly(i,istart) - latpoly(i,iend)) * &
                     pi/180. * R_Earth
             !
          ELSE IF ( ABS(latpoly(i,istart)-latpoly(i,iend)) < EPSILON(latpoly) ) THEN
             !
             ! Distance along a circle of constant latitude
             !
             coslat = MAX(COS(latpoly(i,istart) * pi/180.), mincos)
             seglength(i,iseg) = ABS( lonpoly(i,istart) - lonpoly(i,iend) ) * &
                  pi/180. * R_Earth * coslat
             !
          ELSE
             CALL ipslerr(3, "haversine_laloseglen", "The polygon here does not originate from a regular", &
               &       "latitude longitude grid.","")
          ENDIF
          !
       ENDDO
    ENDDO
    !
  END SUBROUTINE haversine_laloseglen
!!
!!  =============================================================================================================================
!! SUBROUTINE: haversine_clockwise   
!!
!>\BRIEF       Get the indices which sort a polygon in a clockwise order starting from an
!              initial vertex given by start.	  
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================
!!
  SUBROUTINE haversine_clockwise(nbseg, heading, start, sortindex)
    !
    ! Find the order of the polygone vertices which start at "start" and 
    ! follow in a clockwise direction.
    !
    ! 0.1 Input Variables
    ! Number of segments
    INTEGER(i_std), INTENT(in)                                 :: nbseg
    REAL(r_std), DIMENSION(nbseg*2), INTENT(in)                :: heading
    REAL(r_std), INTENT(in)                                    :: start
    ! 
    ! 0.2 Output variable
    !
    INTEGER(i_std), DIMENSION(nbseg*2), INTENT(out)            :: sortindex
    !
    ! 0.3 Local variables
    !
    INTEGER(i_std) :: is, js, imin(1)
    REAL(r_std)    :: delang
    REAL(r_std)    :: undef = 9999999999.99999
    REAL(r_std), DIMENSION(nbseg*2) :: workhead
    !
    delang = 360.0/(nbseg*2)
    !
    workhead(:) = heading(:)
    !
    DO is=1,nbseg*2
       !
       ! Compute the difference of heading to the next target angle
       !
       DO js=1,nbseg*2
          IF ( workhead(js) < undef ) THEN
             workhead(js) = MOD(heading(js)-(start+(is-1)*delang)+360.0, 360.0)
             ! Transfer to -180:180 interval
             IF (workhead(js) > 180.0) workhead(js)=workhead(js)-360.0
          ENDIF
       ENDDO
       !
       ! Locate the vertex closest to that target angle
       !
       imin=MINLOC(ABS(workhead))
       sortindex(is) = imin(1)
       !
       ! Mask this vertex so that it is skipped in the next iteration
       !
       workhead(imin(1)) = undef
       !
    ENDDO
    !
  END SUBROUTINE haversine_clockwise
!!
!!  =============================================================================================================================
!! SUBROUTINE: haversine_polyarea   
!!
!>\BRIEF       Computes the area covered by a polygon on the sphere.	  
!!
!! DESCRIPTION:	Computes the area of each polygon based on Girard's theorem. It
!!              states that the area of a polygon of great circles is R**2 times 
!!              the sum of the angles between the polygons minus (N-2)*pi where N 
!!             x is number of corners.  
!!
!! \n
!_ ==============================================================================================================================
!!
  SUBROUTINE haversine_polyarea(nbpt, nbseg, lonpoly, latpoly, area)
    !
    !
    ! 0.1 Input Variables
    ! Size of the gathered domain
    INTEGER(i_std), INTENT(in)                                 :: nbpt
    ! Number of segments
    INTEGER(i_std), INTENT(in)                                 :: nbseg
    !
    REAL(r_std), DIMENSION(nbpt,nbseg*2), INTENT(in)           :: lonpoly
    REAL(r_std), DIMENSION(nbpt,nbseg*2), INTENT(in)           :: latpoly
    REAL(r_std), DIMENSION(nbpt), INTENT(out)                  :: area
    !
    ! 0.2 Local variables
    !
    INTEGER(i_std) :: i, ia, ib1, ib2, iseg
    REAL(r_std)    :: beta1, beta2
    REAL(r_std), DIMENSION(nbseg) :: angles
    !
    DO i=1,nbpt
       iseg=0
       DO ia=1,nbseg*2,2
          !! Index of next vertex, i.e. ia+1
          ib1=MOD(ia+1, nbseg*2)+1
          !! Index of previous vertex, i.e. ia-2
          ib2=MOD(ia+nbseg*2-3, nbseg*2)+1
          iseg=iseg+1
          !
          beta1=haversine_dtor(haversine_heading(lonpoly(i,ia), latpoly(i,ia), lonpoly(i,ib1), latpoly(i,ib1)))
          beta2=haversine_dtor(haversine_heading(lonpoly(i,ia), latpoly(i,ia), lonpoly(i,ib2), latpoly(i,ib2)))
          !
          angles(iseg)=acos(cos(-beta1)*cos(-beta2) + sin(-beta1)*sin(-beta2))
          !
       ENDDO
       area(i) = (sum(angles) - (nbseg-2)*pi)*R_Earth**2
    ENDDO
    !
  END SUBROUTINE haversine_polyarea
!!
!!  =============================================================================================================================
!! SUBROUTINE: haversine_laloarea   
!!
!>\BRIEF       Computes the area of a regular latitude longitude box for which we already have the
!!             the segment length.	  
!!
!! DESCRIPTION:	Just verify that we have 4 segments. 
!!
!! \n
!_ ==============================================================================================================================
!!
  SUBROUTINE haversine_laloarea(nbpt, nbseg, seglen, area)
    !
    !
    ! 0.1 Input Variables
    ! Size of the gathered domain
    INTEGER(i_std), INTENT(in)                                 :: nbpt
    ! Number of segments
    INTEGER(i_std), INTENT(in)                                 :: nbseg
    !
    REAL(r_std), DIMENSION(nbpt,nbseg), INTENT(in)             :: seglen
    REAL(r_std), DIMENSION(nbpt), INTENT(out)                  :: area
    !
    ! 0.2 Local variables
    !
    INTEGER(i_std) :: i
    !
    IF ( nbseg .NE. 4 ) THEN
          CALL ipslerr(3, "haversine_laloarea", "We need to have 4 segments in the polygons", &
               &       "to use this subroutine","")
    ENDIF
    !
    DO i=1,nbpt
       area(i) = (seglen(i,1)+seglen(i,3))/2.0*(seglen(i,2)+seglen(i,4))/2.0
    ENDDO
    !
  END SUBROUTINE haversine_laloarea
!!
!!  =============================================================================================================================
!! SUBROUTINE: haversine_laloarea   
!!
!>\BRIEF       Computes the area in the special case where the projection has given us the size in X and Y of each
!!             grid box.
!!
!! DESCRIPTION:
!!
!! \n
!_ ==============================================================================================================================
!!
  SUBROUTINE haversine_xyarea(nbpt, nbseg, ilandtoij, jlandtoij, dx, dy, area)
    !
    !
    ! 0.1 Input Variables
    ! Size of the gathered domain
    INTEGER(i_std), INTENT(in)                                 :: nbpt
    ! Number of segments
    INTEGER(i_std), INTENT(in)                                 :: nbseg
    !
    INTEGER(i_std), DIMENSION(nbpt), INTENT(in)                :: ilandtoij, jlandtoij
    REAL(r_std), DIMENSION(:,:), INTENT(in)                    :: dx, dy
    REAL(r_std), DIMENSION(nbpt), INTENT(out)                  :: area
    !
    ! 0.2 Local variables
    !
    INTEGER(i_std) :: il
    !
    DO il=1,nbpt
       area(il) = dx(ilandtoij(il),jlandtoij(il))*dy(ilandtoij(il),jlandtoij(il))
    ENDDO
    !
  END SUBROUTINE haversine_xyarea
!!
!!  =============================================================================================================================
!! FUNCTIONS: haversine_heading,  haversine_distance, haversine_dtor, haversine_rtod 
!!
!>\BRIEF Various functions to help with the calculations in this module.	  
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================
  REAL(r_std) FUNCTION haversine_heading(lon_start, lat_start, lon_end, lat_end)
    !
    ! The heading is an angle in degrees between 0 and 360.
    !
    IMPLICIT NONE
    REAL(r_std), INTENT(IN) :: lon_start, lat_start, lon_end, lat_end
    !
    REAL(r_std) :: dlon, lat1, lat2, y, x
    !
    dlon = haversine_dtor(lon_end-lon_start)
    lat1 = haversine_dtor(lat_start)
    lat2 = haversine_dtor(lat_end)
    y = sin(dlon) * cos(lat2)
    x = cos(lat1)* sin(lat2) - sin(lat1) * cos(lat2)* cos(dLon)
    haversine_heading = MOD(haversine_rtod(atan2(y,x))+360.0, 360.0)
    !
  END FUNCTION haversine_heading
  !!
  !!
  !!
  REAL(r_std) FUNCTION haversine_distance(lon_start, lat_start, lon_end, lat_end)
    IMPLICIT NONE
    REAL(r_std), INTENT(IN) :: lon_start, lat_start, lon_end, lat_end
    !
    REAL(r_std) :: dlon, dlat, lat1, lat2, a, c
    !
    dlat = haversine_dtor(lat_end-lat_start)
    dlon = haversine_dtor(lon_end-lon_start)
    lat1 = haversine_dtor(lat_start)
    lat2 = haversine_dtor(lat_end)
    a = sin(dlat/2) * sin(dlat/2) + sin(dlon/2) * sin(dlon/2) * cos(lat1) * cos(lat2)
    c = 2 * atan2(sqrt(a), sqrt(1-a))
    haversine_distance = c * R_Earth
    !
  END FUNCTION haversine_distance
  !!
  !!
  !!
  FUNCTION haversine_radialdis(lon_start, lat_start, head, dis)
    IMPLICIT NONE
    REAL(r_std), INTENT(IN)   :: lon_start, lat_start, head, dis
    REAL(r_std), DIMENSION(2) :: haversine_radialdis
    !
    REAL(r_std) :: lonout, latout, lat1, lon1, lat, lon, dlon, tc
    !
    lon1 = haversine_dtor(lon_start)
    lat1 = haversine_dtor(lat_start)
    tc = haversine_dtor(head)
    !
    lat =asin(sin(lat1)*cos(dis/R_Earth)+cos(lat1)*sin(dis/R_Earth)*cos(tc))
    dlon = atan2(sin(tc)*sin(dis/R_Earth)*cos(lat1),cos(dis/R_Earth)-sin(lat1)*sin(lat))
    lon = MOD(lon1+dlon+pi, 2*pi)-pi
    !
    latout = haversine_rtod(lat)
    lonout = haversine_rtod(lon)
    !
    haversine_radialdis(1) = lonout
    haversine_radialdis(2) = latout
    !
  END FUNCTION haversine_radialdis
  ! --------------------------------------------------------------------
  ! FUNCTION  haversine_dtor():
  !    This function takes a REAL argument in degree and converts it to
  ! the equivalent radian.
  ! --------------------------------------------------------------------
  REAL(r_std) FUNCTION  haversine_dtor(Degree)
    IMPLICIT  NONE
    REAL(r_std), INTENT(IN) :: Degree
    haversine_dtor = Degree * pi/180.0
  END FUNCTION haversine_dtor
  ! --------------------------------------------------------------------
  ! FUNCTION  haversine_rtod():
  !    This function takes a REAL argument in radian and converts it to
  ! the equivalent degree.
  ! --------------------------------------------------------------------
  REAL(r_std) FUNCTION  haversine_rtod(Radian)
    IMPLICIT  NONE
    REAL(r_std), INTENT(IN) :: Radian

    haversine_rtod = Radian * 180.0/pi
  END FUNCTION haversine_rtod
  !!
  !!
  !!
END MODULE haversine
