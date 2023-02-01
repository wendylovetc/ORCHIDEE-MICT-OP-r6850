
MODULE polygones
  !
  ! This module is there to handle polygones. It is an essential tool in interpolation methods.
  !
  USE defprec

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: polygones_pointinside, polygones_extend, polygones_intersection, polygones_cleanup, &
       &    polygones_area, polygones_crossing, polygones_convexhull
  !
  REAL(r_std),PARAMETER                     :: zero=0.0
  !
CONTAINS
!
!=======================================================================================================
!
  SUBROUTINE polygones_pointinside(nvert_in, poly, point_x, point_y, inside)
    !
    ! Routine to determine of point (point_x, point_y) is inside of polygone poly
    ! Mathode based on docume http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
    !
    ! Arguments
    INTEGER(i_std), INTENT(in)              :: nvert_in
    REAL(r_std), DIMENSION(:,:), INTENT(in) :: poly
    REAL(r_std), INTENT(in)                 :: point_x, point_y
    LOGICAL, INTENT(out)                    :: inside
    !
    ! Local
    INTEGER(i_std)                          :: i, j, nvert
    REAL(r_std)                             :: xonline
    !
    !
    inside = .FALSE.
    !
    nvert=SIZE(poly,DIM=1)
    IF ( nvert < nvert_in) THEN
       WRITE(*,*) "Input polygone is smaller than the size proposed in the interface."
       STOP "polygones_pointinside"
    ENDIF
    IF (SIZE(poly,DIM=2) .NE. 2) THEN
       WRITE(*,*) "This cannot be a polygone :", SIZE(poly)
       STOP "polygones_pointinside"
    ENDIF
    !
    j = nvert_in
    DO i=1,nvert_in
       IF ( (poly(i,2) > point_y) .NEQV. (poly(j,2) > point_y) ) THEN
          xonline = (poly(j,1)-poly(i,1))*(point_y-poly(i,2))/(poly(j,2)-poly(i,2))+poly(i,1)
          IF ( point_x < xonline ) THEN
             inside = (.NOT. inside)
          ENDIF
       ENDIF
       j = i
    ENDDO
    !
  END SUBROUTINE polygones_pointinside
!
!=======================================================================================================
!
  SUBROUTINE polygones_lineintersect(la, lb, intersection, point_x, point_y)
    !
    ! Simple alogorith tho determine the intersections of 2 lines.
    !
    ! ARGUMENTS
    REAL(r_std), DIMENSION(2,2), INTENT(in)    :: la, lb
    LOGICAL, INTENT(out)                       :: intersection
    REAL(r_std), INTENT(out)                   :: point_x, point_y
    !
    ! Local
    REAL(r_std)                                :: den, ua, ub
    REAL(r_std)                                :: x1, x2, x3, x4, y1, y2, y3, y4
    !
    intersection = .FALSE.
    !
     x1 = la(1,1)
     y1 = la(1,2)
     x2 = la(2,1)
     y2 = la(2,2)

     x3 = lb(1,1)
     y3 = lb(1,2)
     x4 = lb(2,1)
     y4 = lb(2,2)

     den = (y4-y3)*(x2-x1)-(x4-x3)*(y2-y1)

     IF ( ABS(den) > EPSILON(den) ) THEN
        ua = ((x4-x3)*(y1-y3)-(y4-y3)*(x1-x3))/den
        ub = ((x2-x1)*(y1-y3)-(y2-y1)*(x1-x3))/den
        IF ( ua >= 0 .AND. ua <= 1 .AND. ub >= 0 .AND. ub <= 1) THEN
           intersection = .TRUE.
           point_x = x1 + ua*(x2-x1)
           point_y = y1 + ua*(y2-y1)
        ENDIF
     ENDIF

  END SUBROUTINE polygones_lineintersect
!
!=======================================================================================================
!
  SUBROUTINE polygones_extend(nvert_in, poly_in, nbdots, nvert_out, poly_out)
    !
    ! This function simply add nbdots onto each side of the polygone.
    !
    ! Arguments
    INTEGER(i_std), INTENT(in)               :: nvert_in
    REAL(r_std), DIMENSION(:,:), INTENT(in)  :: poly_in
    INTEGER(i_std), INTENT(in)               :: nbdots
    INTEGER(i_std), INTENT(out)              :: nvert_out
    REAL(r_std), DIMENSION(:,:), INTENT(out) :: poly_out
    !
    ! Local
    !
    INTEGER(i_std)                        :: nvert, nvert_tmp, i, j, id, ipos
    REAL(r_std)                           :: xs, xe, ys, ye
    !
    nvert=SIZE(poly_in,DIM=1)
    IF ( nvert < nvert_in ) THEN
       WRITE(*,*) "Input polygone is smaller than the size proposed in the interface."
       STOP "polygones_extend"
    ENDIF
    IF (SIZE(poly_in,DIM=2) .NE. 2) THEN
       WRITE(*,*) "Poly_in cannot be a polygone :", SIZE(poly_in)
       STOP "polygones_extend"
    ENDIF
    nvert_tmp=SIZE(poly_out,DIM=1)
    IF (SIZE(poly_out,DIM=2) .NE. 2) THEN
       WRITE(*,*) "Poly_out cannot be a polygone :", SIZE(poly_out)
       STOP "polygones_extend"
    ENDIF
    !
    IF ( nvert_tmp < nvert_in*nbdots) THEN
       WRITE(*,*) "Poly_out is too small for extending it :", SIZE(poly_out), " we would need :", nvert_in*nbdots
       STOP "polygones_extend"
    ENDIF
    nvert_out = nvert_in*nbdots
    !
    j = nvert_in
    ipos = 0
    DO i=1,nvert_in
       xs = poly_in(j,1)
       xe = poly_in(i,1)
       ys = poly_in(j,2)
       ye = poly_in(i,2)
       IF (xs == xe ) THEN 
          DO id=1,nbdots 
             ipos = ipos + 1
             poly_out(ipos,1) = xs
             poly_out(ipos,2) = ys + (id-1)*(ye-ys)/nbdots
          ENDDO
       ELSE IF (ys == ye ) THEN
          DO id=1,nbdots
             ipos = ipos + 1
             poly_out(ipos,1) = xs + (id-1)*(xe-xs)/nbdots
             poly_out(ipos,2) = ys
          ENDDO
       ELSE
          DO id=1,nbdots
             ipos = ipos + 1
             poly_out(ipos,2) = ys + (id-1)*(ye-ys)/nbdots
             poly_out(ipos,1) = (xs-xe)*(poly_out(ipos,2)-ye)/(ys-ye)+xe
          ENDDO
       ENDIF
       j = i
    ENDDO
    !
  END SUBROUTINE polygones_extend
!
!=======================================================================================================
!
  SUBROUTINE polygones_intersection(nvert_a, poly_a, nvert_b, poly_b, nbpts, poly_out)
    !
    ! Find the polygon resulting from the intersection of poly_a and poly_b
    !
    ! Arguments
    INTEGER(i_std), INTENT(in)               :: nvert_a
    REAL(r_std), DIMENSION(:,:), INTENT(in)  :: poly_a
    INTEGER(i_std), INTENT(in)               :: nvert_b
    REAL(r_std), DIMENSION(:,:), INTENT(in)  :: poly_b
    INTEGER(i_std), INTENT(out)              :: nbpts
    REAL(r_std), DIMENSION(:,:), INTENT(out) :: poly_out
    !
    ! Local
    INTEGER(i_std)                            :: nvert_tmpa, nvert_tmpb, nvert_out, i, j
    LOGICAL                                   :: inside
    INTEGER(i_std)                            :: in_a, in_b
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:) :: ptin_a, ptin_b
    !
    nvert_tmpa=SIZE(poly_a,DIM=1)
    IF ( nvert_tmpa < nvert_a) THEN
       WRITE(*,*) "Input polygone (poly_a) is smaller than the size proposed in the interface."
       STOP "polygones_intersection"
    ENDIF
    IF (SIZE(poly_a,DIM=2) .NE. 2) THEN
       WRITE(*,*) "Poly_a cannot be a polygone :", SIZE(poly_a)
       STOP "polygones_intersection"
    ENDIF
    ALLOCATE(ptin_a(nvert_a))
    nvert_tmpb=SIZE(poly_b,DIM=1)
    IF ( nvert_tmpb < nvert_a) THEN
       WRITE(*,*) "Input polygone (poly_b) is smaller than the size proposed in the interface."
       STOP "polygones_intersection"
    ENDIF
    IF (SIZE(poly_b,DIM=2) .NE. 2) THEN
       WRITE(*,*) "Poly_b cannot be a polygone :", SIZE(poly_b)
       STOP "polygones_intersection"
    ENDIF
    ALLOCATE(ptin_b(nvert_b))
    !
    !
    in_a = 0
    DO i=1,nvert_a
       CALL polygones_pointinside(nvert_b, poly_b, poly_a(i,1), poly_a(i,2), inside)
       IF ( inside ) THEN
          in_a = in_a + 1
          ptin_a(in_a) = i
       ENDIF
    ENDDO
    !
    in_b = 0
    DO i=1,nvert_b
       CALL polygones_pointinside(nvert_a, poly_a, poly_b(i,1), poly_b(i,2), inside)
       IF ( inside ) THEN
          in_b = in_b + 1
          ptin_b(in_b) = i
       ENDIF
    ENDDO
    !
    nbpts = in_a + in_b
    !
    nvert_out=SIZE(poly_out,DIM=1)
    IF (SIZE(poly_out,DIM=2) .NE. 2) THEN
       WRITE(*,*) "Poly_out cannot be a polygone :", SIZE(poly_out)
       STOP "polygones_intersection"
    ENDIF
    IF ( nbpts <= nvert_out ) THEN
       i = 0
       DO j=1,in_a
          i = i+1
          poly_out(i,:) = poly_a(ptin_a(j),:)
       ENDDO
       DO j=1,in_b
          i = i+1
          poly_out(i,:) = poly_b(ptin_b(j),:)
       ENDDO
    ELSE
       WRITE(*,*) "The intersection polygone has", nbpts, "points and that is more than available in poly_out", SIZE(poly_out)
       STOP "polygones_intersection"
    ENDIF
    !
    DEALLOCATE(ptin_a, ptin_b)
    !
  END SUBROUTINE polygones_intersection
!
!=======================================================================================================
!
  SUBROUTINE polygones_cleanup(nvert_in, poly_in, nvert_out, poly_out)
    !
    ! Clean_up the polygone by deleting points which are redundant and ordering based 
    ! on the proximity of the points. Beware this does not give a convex hull to the surface
    ! inside the polygone.
    !
    ! Arguments
    INTEGER(i_std), INTENT(in)               :: nvert_in
    REAL(r_std), DIMENSION(:,:), INTENT(in)  :: poly_in
    INTEGER(i_std), INTENT(out)              :: nvert_out
    REAL(r_std), DIMENSION(:,:), INTENT(out) :: poly_out
    !
    ! Local
    !
    INTEGER(i_std)                           :: nvert, nvert_tmpin, nvert_tmp
    INTEGER(i_std)                           :: i, j
    INTEGER(i_std), DIMENSION(1)             :: ismin
    REAL(r_std), ALLOCATABLE, DIMENSION(:)   :: dist
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: poly_tmp
    !
    nvert_tmpin=SIZE(poly_in,DIM=1)
    IF ( nvert_tmpin < nvert_in ) THEN
       WRITE(*,*) "Input polygone is smaller than the size proposed in the interface."
       STOP "polygones_cleanup"
    ENDIF
    IF (SIZE(poly_in,DIM=2) .NE. 2) THEN
       WRITE(*,*) "Poly_in cannot be a polygone :", SIZE(poly_in)
       STOP "polygones_cleanup"
    ENDIF
    nvert_tmp=SIZE(poly_out,DIM=1)
    IF (SIZE(poly_out,DIM=2) .NE. 2) THEN
       WRITE(*,*) "Poly_out cannot be a polygone :", SIZE(poly_out)
       STOP "polygones_cleanup"
    ENDIF
    !
    ALLOCATE(poly_tmp(nvert_in,2))
    ALLOCATE(dist(nvert_in))
    !
    nvert = nvert_in-1
    poly_tmp(1:nvert,:) = poly_in(2:nvert_in,:)
    poly_out(1,:) = poly_in(1,:)
    nvert_out = 1
    !
    DO WHILE ( nvert > 0 ) 
       !
       DO j=1,nvert
          dist(j) = SQRT((poly_out(nvert_out,1)-poly_tmp(j,1))**2+(poly_out(nvert_out,2)-poly_tmp(j,2))**2)
       ENDDO
       !
       ismin = MINLOC(dist(1:nvert))
       !
       IF ( dist(ismin(1)) > EPSILON(dist(1)) ) THEN
          !
          ! The smallest distance between 2 vertices is larger than zero :
          ! Add the vertex to poly_out and delete it from poly_tmp
          !
          nvert_out = nvert_out + 1
          IF ( nvert_out > nvert_tmp) THEN
             WRITE(*,*) "Output polygone too small"
             STOP "polygones_cleanup"
          ENDIF
          poly_out(nvert_out,:) = poly_tmp(ismin(1),:)
          poly_tmp(ismin(1):nvert-1,:) = poly_tmp(ismin(1)+1:nvert,:)
       ELSE
          !
          ! Else the vertex already exists in poly_out and thus we just need
          ! to delete it from poly_tmp
          !
          poly_tmp(ismin(1):nvert-1,:) = poly_tmp(ismin(1)+1:nvert,:)
       ENDIF
       !
       ! One less vertices exists in poly_tmp
       !
       nvert = nvert-1
       !
    ENDDO

    DEALLOCATE(poly_tmp, dist)

  END SUBROUTINE polygones_cleanup
!
!=======================================================================================================
!
  SUBROUTINE polygones_area(nvert_in, poly_in, dx, dy, area)
    !
    ! Find the polygon resulting from the intersection of poly_a and poly_b
    !
    ! Arguments
    INTEGER(i_std), INTENT(in)               :: nvert_in
    REAL(r_std), DIMENSION(:,:), INTENT(in)  :: poly_in
    REAL(r_std), INTENT(in)                  :: dx, dy
    REAL(r_std), INTENT(out)                 :: area
    !
    ! Local
    !
    INTEGER(i_std)     :: nvert, i, j
    !
    nvert = SIZE(poly_in,DIM=1)
    IF ( nvert < nvert_in) THEN
       WRITE(*,*) "Input polygone is smaller than the size proposed in the interface."
       STOP "polygones_area"
    ENDIF
    IF (SIZE(poly_in,DIM=2) .NE. 2) THEN
       WRITE(*,*) "Poly_in cannot be a polygone :", SIZE(poly_in)
       STOP "polygones_area"
    ENDIF
    !
    area = zero
    j = nvert_in
    DO i=1,nvert_in
!!       area = area + (dy*(poly_in(j,2)+poly_in(i,2))/2.0)*(dx*(poly_in(j,1)-poly_in(i,1)))
       area = area + dy*dx/2.0*(poly_in(j,2)+poly_in(i,2))*(poly_in(j,1)-poly_in(i,1))
       j = i
    ENDDO
    !
    area = ABS(area)
    !
  END SUBROUTINE polygones_area
!
!=======================================================================================================
!
  SUBROUTINE polygones_crossing(nvert_a, poly_a, nvert_b, poly_b, nbpts, crossings)
    !
    ! Find the polygon resulting from the intersection of poly_a and poly_b
    !
    ! Arguments
    INTEGER(i_std), INTENT(in)               :: nvert_a
    REAL(r_std), DIMENSION(:,:), INTENT(in)  :: poly_a
    INTEGER(i_std), INTENT(in)               :: nvert_b
    REAL(r_std), DIMENSION(:,:), INTENT(in)  :: poly_b
    INTEGER(i_std), INTENT(out)              :: nbpts
    REAL(r_std), DIMENSION(:,:), INTENT(out) :: crossings
    !
    ! Local
    !
    INTEGER(i_std)                            :: nvert, ia, ja, ib, jb
    REAL(r_std), DIMENSION(2,2)               :: la, lb
    REAL(r_std)                               :: x, y
    LOGICAL                                   :: intersect
    !
    nvert=SIZE(poly_a,DIM=1)
    IF ( nvert < nvert_a) THEN
       WRITE(*,*) "Input polygone (poly_a) is smaller than the size proposed in the interface."
       STOP "polygones_crossing"
    ENDIF
    IF (SIZE(poly_a,DIM=2) .NE. 2) THEN
       WRITE(*,*) "Poly_a cannot be a polygone :", SIZE(poly_a)
       STOP "polygones_crossing"
    ENDIF
    nvert=SIZE(poly_b,DIM=1)
    IF ( nvert < nvert_b) THEN
       WRITE(*,*) "Input polygone (poly_b) is smaller than the size proposed in the interface."
       STOP "polygones_crossing"
    ENDIF
    IF (SIZE(poly_b,DIM=2) .NE. 2) THEN
       WRITE(*,*) "Poly_b cannot be a polygone :", SIZE(poly_b)
       STOP "polygones_crossing"
    ENDIF
    nvert=SIZE(crossings,DIM=1)
    IF (SIZE(crossings,DIM=2) .NE. 2) THEN
       WRITE(*,*) "Poly_a cannot be a polygone :", SIZE(poly_a)
       STOP "polygones_crossing"
    ENDIF
    !
    nbpts = 0
    !
    ja = nvert_a
    DO ia=1,nvert_a
       !
       ! Put one side of the polygone a into la
       !
       la(1,1) = poly_a(ja,1)
       la(1,2) = poly_a(ja,2)
       la(2,1) = poly_a(ia,1)
       la(2,2) = poly_a(ia,2)
       !
       jb = nvert_b
       DO ib=1,nvert_b
          !
          ! Put one side of the polygone b into lb
          !
          lb(1,1) = poly_b(jb,1)
          lb(1,2) = poly_b(jb,2)
          lb(2,1) = poly_b(ib,1)
          lb(2,2) = poly_b(ib,2)
          !
          CALL polygones_lineintersect(la, lb, intersect, x, y)
          !
          IF ( intersect ) THEN
             nbpts = nbpts + 1
             IF ( nbpts <= nvert ) THEN
                crossings(nbpts, 1) = x
                crossings(nbpts, 2) = y
             ELSE
                WRITE(*,*) "Polygone to write the intersection points is too small", nvert, nbpts
                STOP "polygones_crossing"
             ENDIF
          ENDIF
          !
          jb = ib
       ENDDO
       ja = ia
    ENDDO
    !

  END SUBROUTINE polygones_crossing
!
!=======================================================================================================
!
  SUBROUTINE polygones_convexhull(nvert_in, poly_in, nvert_out, poly_out)
    !
    ! Routine orders points in poly_in so that they form a convex shape. This is based on the Graham scan
    ! algorith as implemented by Alan Miller : http://jblevins.org/mirror/amiller/
    !
    ! The output polygone can be smaller than the input one as we are computing the convex envelope.
    !
    ! Arguments
    INTEGER(i_std), INTENT(in)               :: nvert_in
    REAL(r_std), DIMENSION(:,:), INTENT(in)  :: poly_in
    INTEGER(i_std), INTENT(out)              :: nvert_out
    REAL(r_std), DIMENSION(:,:), INTENT(out) :: poly_out
    !
    ! Local
    !
    INTEGER(i_std)                           :: nvert, nvert_tmp
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:):: iwk, next, vertex
    REAL(r_std)                              :: xmin, temp, dy, dx2, dy1, dy2
    REAL(r_std)                              :: x1, x2, y1, y2
    REAL(r_std)                              :: xmax, ymin, ymax
    REAL(r_std)                              :: dx, dx1
    REAL(r_std)                              :: dmax, dmax1, dmax2
    REAL(r_std)                              :: dist, dmin
    INTEGER(i_std)                           :: i, i1, i2, j, jp1, jp2, i2save, i3, i2next
    LOGICAL                                  :: points_todo
    INTEGER(i_std), PARAMETER                :: nextinc=20
    !
    nvert_tmp=SIZE(poly_in,DIM=1)
    IF ( nvert_tmp < nvert_in) THEN
       WRITE(*,*) "Input polygone is smaller than the size proposed in the interface."
       STOP "polygones_convexhull"
    ENDIF
    IF ( nvert_tmp <= 2 ) THEN
       WRITE(*,*) "The input polygone is too small to compute a convex hull."
       STOP "polygones_convexhull"
    ENDIF
    IF (SIZE(poly_in,DIM=2) .NE. 2) THEN
       WRITE(*,*) "Poly_a cannot be a polygone :", SIZE(poly_in)
       STOP "polygones_convexhull"
    ENDIF
    IF (SIZE(poly_out,DIM=2) .NE. 2) THEN
       WRITE(*,*) "Poly_b cannot be a polygone :", SIZE(poly_out)
       STOP "polygones_convexhull"
    ENDIF
    !
    ALLOCATE(vertex(nvert_in))
    ALLOCATE(iwk(nvert_in))
    ALLOCATE(next(nvert_in*nextinc))
    !
    !
    !  Choose the points with smallest & largest x- values as the
    !  first two vertices of the polygon.
    !
    IF (poly_in(1,1) > poly_in(nvert_in,1)) THEN
       vertex(1) = nvert_in
       vertex(2) = 1
       xmin = poly_in(nvert_in,1)
       xmax = poly_in(1,1)
    ELSE
       vertex(1) = 1
       vertex(2) = nvert_in
       xmin = poly_in(1,1)
       xmax = poly_in(nvert_in,1)
    END IF
    
    DO i = 2, nvert_in-1
       temp = poly_in(i,1)
       IF (temp < xmin) THEN
          vertex(1) = i
          xmin = temp
       ELSE IF (temp > xmax) THEN
          vertex(2) = i
          xmax = temp
       END IF
    END DO
    !
    !       Special case, xmax = xmin.
    !
    IF (xmax == xmin) THEN
       IF (poly_in(1,2) > poly_in(nvert_in,2)) THEN
          vertex(1) = nvert_in
          vertex(2) = 1
          ymin = poly_in(nvert_in,2)
          ymax = poly_in(1,2)
       ELSE
          vertex(1) = 1
          vertex(2) = nvert_in
          ymin = poly_in(1,2)
          ymax = poly_in(nvert_in,2)
       END IF
       
       DO i = 2, nvert_in-1
          temp = poly_in(i,2)
          IF (temp < ymin) THEN
             vertex(1) = i
             ymin = temp
          ELSE IF (temp > ymax) THEN
             vertex(2) = i
             ymax = temp
          END IF
       END DO
       
       nvert = 2
       IF (ymax == ymin) nvert = 1
       RETURN
    END IF
    !
    !  Set up two initial lists of points; those points above & those below the
    !  line joining the first two vertices.    next(i) will hold the pointer to the
    !  point furthest from the line joining vertex(i) to vertex(i+1) on the left
    !  hand side.
    !
    i1 = vertex(1)
    i2 = vertex(2)
    iwk(i1) = -1
    iwk(i2) = -1
    dx = xmax - xmin
    y1 = poly_in(i1,2)
    dy = poly_in(i2,2) - y1
    dmax = zero
    dmin = zero
    next(1) = -1
    next(2) = -1
    !
    DO i = 1, nvert_in
       IF (i == vertex(1) .OR. i == vertex(2)) CYCLE
       dist = (poly_in(i,2) - y1)*dx - (poly_in(i,1) - xmin)*dy
       IF (dist > zero) THEN
          iwk(i1) = i
          i1 = i
          IF (dist > dmax) THEN
             next(1) = i
             dmax = dist
          END IF
       ELSE IF (dist < zero) THEN
          iwk(i2) = i
          i2 = i
          IF (dist < dmin) THEN
             next(2) = i
             dmin = dist
          END IF
       END IF
    END DO
    !
    !  Ends of lists are indicated by pointers to -ve positions.
    !
    iwk(i1) = -1
    iwk(i2) = -1
    nvert = 2
    !
    j = 1
    !
    !  Start of main process.
    !
    !  Introduce new vertex between vertices j & j+1, if one has been found.
    !  Otherwise increase j.   Exit if no more vertices.
    !
    !
    points_todo = .TRUE.
    DO WHILE ( points_todo )
       !
       DO WHILE ( points_todo .AND. next(j) < 0 )
          IF (j == nvert) points_todo = .FALSE.
          j = j + 1
       ENDDO
          !
       IF ( points_todo ) THEN
          !
          jp1 = j + 1
          IF ( jp1 >= nvert_in*nextinc) THEN
             STOP "polygones_convexhull : please increase nextinc"
          ENDIF
          !
          DO i = nvert, jp1, -1
             vertex(i+1) = vertex(i)
             next(i+1) = next(i)
          END DO
          jp2 = jp1 + 1
          nvert = nvert + 1
          IF (jp2 > nvert) jp2 = 1
          i1 = vertex(j)
          i2 = next(j)
          i3 = vertex(jp2)
          vertex(jp1) = i2
          !
          !  Process the list of points associated with vertex j.   New list at vertex j
          !  consists of those points to the left of the line joining it to the new
          !  vertex (j+1).   Similarly for the list at the new vertex.
          !  Points on or to the right of these lines are dropped.
          !
          x1 = poly_in(i1,1)
          x2 = poly_in(i2,1)
          y1 = poly_in(i1,2)
          y2 = poly_in(i2,2)
          !
          dx1 = x2 - x1
          dx2 = poly_in(i3,1) - x2
          dy1 = y2 - y1
          dy2 = poly_in(i3,2) - y2
          dmax1 = zero
          dmax2 = zero
          next(j) = -1
          next(jp1) = -1
          i2save = i2
          i2next = iwk(i2)
          i = iwk(i1)
          iwk(i1) = -1
          iwk(i2) = -1
          !
          DO WHILE ( i > 0 )
             IF (i /= i2save) THEN
                dist = (poly_in(i,2) - y1)*dx1 - (poly_in(i,1) - x1)*dy1
                IF (dist > zero) THEN
                   iwk(i1) = i
                   i1 = i
                   IF (dist > DMAX1) THEN
                      next(j) = i
                      dmax1 = dist
                   END IF
                ELSE
                   dist = (poly_in(i,2) - y2)*dx2 - (poly_in(i,1) - x2)*dy2
                   IF (dist > zero) THEN
                      iwk(i2) = i
                      i2 = i
                      IF (dist > dmax2) THEN
                         next(jp1) = i
                         dmax2 = dist
                      END IF
                   END IF
                END IF
                i = iwk(i)
             ELSE
                i = i2next
             END IF
             !
             !  Get next point from old list at vertex j.
             !
          ENDDO
          !
          !  End lists with -ve values.
          !
          iwk(i1) = -1
          iwk(i2) = -1
          !
       ENDIF
    ENDDO
    !
    ! Copy the polygone in he right order to the output variable
    !
    nvert_out = nvert
    nvert_tmp=SIZE(poly_out,DIM=1)
    IF ( nvert_tmp < nvert_out) THEN
       WRITE(*,*) "Ouptput polygone is smaller than the size proposed in the interface."
       STOP "polygones_convexhull"
    ENDIF
    !
    DO i=1,nvert
       poly_out(i,:) = poly_in(vertex(i),:)
    ENDDO
    !
    DEALLOCATE(vertex, iwk, next)
    !
  END SUBROUTINE polygones_convexhull

END MODULE polygones
