SUBROUTINE node_db
  USE parameters
  USE geometry
  USE constants
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Subroutine node_db initializes the database associated with each node in   !
!  the computational domain.  (1) it initializes the table that gives the     !
!  triangles associated with each node, n2t(); (2) it computes the area of    !
!  the control volume associated with each node, cv_area(); (3) it computes   !
!  computes the control volume perimeters, cv_perim(); and (4) it builds the  !
!  table that gives, for each node in the domain, all the nodes connected to  !
!  it by an edge.                                                             !
!                                                                             !
!  F. Simoes, November 2004                                                   !
!  Last updated (mm-dd-yyyy): 08-01-2007 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables:
  INTEGER :: e1,e2,e3,i,ierror,j,k,m,p1,p2,p3,t1,t2
  INTEGER :: counter(n_pts)  ! Counts the number of triangles associated with
                             ! each node.
  REAL (KIND=mp) :: x1,x2,xm,y1,y2,ym
  LOGICAL, EXTERNAL :: node_in_triang,pto_in_edge,node_in_bdry

  WRITE (*,'("  Building node database...",$)')

! Find out the maximum number of triangles associated with each node.
  counter = 0
  DO i = 1,n_elems
    DO j = 1,3
      counter(grid(i)%vertex(j)) = counter(grid(i)%vertex(j)) + 1
    END DO
  END DO
  k = MAXVAL(counter)

! Allocate variables.
  ALLOCATE(cv_area(n_pts),cv_perim(n_pts),n2t(n_pts,k+1),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add((8*2 + 4*(k + 1))*n_pts)
  m = MAX(k,4)  ! Needed later for n2t2 table.

! Set-up the node-to-triangle connectivity.
  n2t = 0
  ! The old, slower way of doing it:
  !DO i = 1,n_pts
  !  k = 0
  !  DO j = 1,n_elems
  !    IF (node_in_triang(i,grid(j))) THEN  ! Add triangle j to node i table.
  !      k = k + 1
  !      n2t(i,k+1) = j
  !    END IF
  !  END DO
  !  n2t(i,1) = k  ! Number of triangles associated with node i.
  !END DO
  ! The newer, faster way of doing it:
  DO i = 1,n_elems
    DO j = 1,3
      p1 = grid(i)%vertex(j)
      k = n2t(p1,1)
      k = k + 1
      n2t(p1,k+1) = i
      n2t(p1,1) = k
    END DO
  END DO

! The above n2t table only contains the immediate neighborhood of a node, i.e.,
! it strictly contains the triangles that share that node.  For least-squares
! reconstruction this neighborhood may not constitute a good enough
! computational molecule, such as nodes in the boundaries shared by only one
! or two triangles, which may happen in corners.  The n2t2 table should be used
! for least square reconstruction: it is built from n2t by adding additional
! triangles when the computational molecute is less than three, i.e., when a
! node is shared by less than three triangles.  The n2t2 computational molecule
! has at least 4 triangles, except in very pathological cases with strange grid
! configurations.
  ALLOCATE(n2t2(n_pts,m+1),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(4*(m + 1)*n_pts)
  DO i = 1,n_pts
    DO j = 1,n2t(i,1)+1
      n2t2(i,j) = n2t(i,j)
    END DO
  END DO

  DO i = 1,n_pts
    k = n2t2(i,1)
    SELECT CASE (k)
    CASE (1)  ! Molecules with only one triangle.
      t1 = n2t2(i,2)
      ! Find the second triangle, which shares and edge with t1.
      DO j = 1,3
        e1 = ABS(grid(t1)%edge(j))
        IF (edges(e1)%e(1) > 0 .AND. edges(e1)%e(1) /= t1) THEN
          t2 = edges(e1)%e(1)
          EXIT
        ELSE IF (edges(e1)%e(2) > 0 .AND. edges(e1)%e(2) /= t1) THEN
          t2 = edges(e1)%e(2)
          EXIT
        END IF
      END DO
      n2t2(i,3) = t2
      k = k + 1
      ! Find the other triangles, which share an edge with t2.
      DO j = 1,3
        e1 = ABS(grid(t2)%edge(j))
        IF (edges(e1)%e(1) > 0 .AND. edges(e1)%e(1) /= t1 .AND. &
          edges(e1)%e(1) /= t2) THEN
          k = k + 1
          n2t2(i,k+1) = edges(e1)%e(1)
        ELSE IF (edges(e1)%e(2) > 0 .AND. edges(e1)%e(2) /= t1 .AND. &
          edges(e1)%e(2) /= t2) THEN
          k = k + 1
          n2t2(i,k+1) = edges(e1)%e(2)
        END IF
      END DO
      n2t2(i,1) = k

    CASE (2)  ! Molecules with two triangles
      t1 = n2t2(i,2)
      t2 = n2t2(i,3)
      ! Find the triangle that shares an edge with t1.
      DO j = 1,3
        e1 = ABS(grid(t1)%edge(j))
        IF (edges(e1)%e(1) > 0 .AND. edges(e1)%e(1) /= t1 .AND. &
          edges(e1)%e(1) /= t2) THEN
          k = k + 1
          n2t2(i,k+1) = edges(e1)%e(1)
          EXIT
        ELSE IF (edges(e1)%e(2) > 0 .AND. edges(e1)%e(2) /= t1 .AND. &
          edges(e1)%e(2) /= t2) THEN
          k = k + 1
          n2t2(i,k+1) = edges(e1)%e(2)
          EXIT
        END IF
      END DO
      ! Find the triangle that shares an edge with t2.
      DO j = 1,3
        e1 = ABS(grid(t2)%edge(j))
        IF (edges(e1)%e(1) > 0 .AND. edges(e1)%e(1) /= t1 .AND. &
          edges(e1)%e(1) /= t2) THEN
          k = k + 1
          n2t2(i,k+1) = edges(e1)%e(1)
          EXIT
        ELSE IF (edges(e1)%e(2) > 0 .AND. edges(e1)%e(2) /= t1 .AND. &
          edges(e1)%e(2) /= t2) THEN
          k = k + 1
          n2t2(i,k+1) = edges(e1)%e(2)
          EXIT
        END IF
      END DO
      n2t2(i,1) = k

    END SELECT
  END DO

! Compute control volume areas.  The total area is the sum of 1/3 of the areas
! of all traingles sharing the node.
  DO i = 1,n_pts
    cv_area(i) = zero
    DO j = 2,n2t(i,1)+1
      cv_area(i) = cv_area(i) + grid(n2t(i,j))%area
    END DO
    cv_area(i) = cv_area(i)/3.0_mp
  END DO

! Compute the control volume perimeters.
  DO i = 1,n_pts
    cv_perim(i) = zero

    ! Sum over all the triangles having point 'i' as a vertex.
    DO j = 2,n2t(i,1)+1
      k = n2t(i,j)  ! Triangle k.
      e1 = ABS(grid(k)%edge(1))  ! Edges in triangle k.
      e2 = ABS(grid(k)%edge(2))
      e3 = ABS(grid(k)%edge(3))

      ! Find the midpoints of the relevant edges.
      IF (pto_in_edge(i,edges(e1))) THEN  ! Point i is in edge 1; otherwise...
        p1 = edges(e1)%p(1)  ! One end of the edge...
        p2 = edges(e1)%p(2)  ! ...and the other end of the edge.
        x1 = half*(nodes(p1)%x + nodes(p2)%x)
        y1 = half*(nodes(p1)%y + nodes(p2)%y)
      ELSE  ! ...point i must be in edge 2.
        p1 = edges(e2)%p(1)  ! One end of the edge...
        p2 = edges(e2)%p(2)  ! ...and the other end of the edge.
        x1 = half*(nodes(p1)%x + nodes(p2)%x)
        y1 = half*(nodes(p1)%y + nodes(p2)%y)
      END IF

      IF (pto_in_edge(i,edges(e3))) THEN  ! Point i is in edge 3; otherwise...
        p1 = edges(e3)%p(1)  ! One end of the edge...
        p2 = edges(e3)%p(2)  ! ...and the other end of the edge.
        x2 = half*(nodes(p1)%x + nodes(p2)%x)
        y2 = half*(nodes(p1)%y + nodes(p2)%y)
      ELSE  ! ...point i must be in edge 2.
        p1 = edges(e2)%p(1)  ! One end of the edge...
        p2 = edges(e2)%p(2)  ! ...and the other end of the edge.
        x2 = half*(nodes(p1)%x + nodes(p2)%x)
        y2 = half*(nodes(p1)%y + nodes(p2)%y)
      END IF

      ! Triangle's mid point.
      p1 = grid(k)%vertex(1)
      p2 = grid(k)%vertex(2)
      p3 = grid(k)%vertex(3)
      xm = (nodes(p1)%x + nodes(p2)%x + nodes(p3)%x)*one_third
      ym = (nodes(p1)%y + nodes(p2)%y + nodes(p3)%y)*one_third

      ! Compute distances and add to cv_perim().
      cv_perim(i) = cv_perim(i) + SQRT((xm - x1)**2 + (ym - y1)**2)
      cv_perim(i) = cv_perim(i) + SQRT((xm - x2)**2 + (ym - y2)**2)
    END DO

    ! Now check if point 'i' is on a solid boundary and add the appropriate
    ! part of the wetted perimeter.
    IF (node_in_bdry(i,j,k)) THEN
      p1 = edges(j)%p(1)  ! Edge j.
      p2 = edges(j)%p(2)
      x1 = half*(nodes(p1)%x + nodes(p2)%x)
      y1 = half*(nodes(p1)%y + nodes(p2)%y)
      p1 = edges(k)%p(1)  ! Edge k.
      p2 = edges(k)%p(2)
      x2 = half*(nodes(p1)%x + nodes(p2)%x)
      y2 = half*(nodes(p1)%y + nodes(p2)%y)
      xm = nodes(i)%x
      ym = nodes(i)%y

      cv_perim(i) = cv_perim(i) + SQRT((xm - x1)**2 + (ym - y1)**2)
      cv_perim(i) = cv_perim(i) + SQRT((xm - x2)**2 + (ym - y2)**2)
    END IF

  END DO

! Find out the maximum number of edges associated with each node.
  counter = 0
  DO i = 1,n_edges
    counter(edges(i)%p(1)) = counter(edges(i)%p(1)) + 1
    counter(edges(i)%p(2)) = counter(edges(i)%p(2)) + 1
  END DO
  k = MAXVAL(counter)

! Allocate table space.
  ALLOCATE(n2n(n_pts,k+1),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(4*(k + 1)*n_pts)
  n2n = 0

! Compute the node-to-node connectivity table.
  DO i = 1,n_edges
    p1 = edges(i)%p(1)
    p2 = edges(i)%p(2)
    k = n2n(p1,1)
    k = k + 1
    n2n(p1,k+1) = p2
    n2n(p1,1) = k
    k = n2n(p2,1)
    k = k + 1
    n2n(p2,k+1) = p1
    n2n(p2,1) = k
  END DO

  PRINT *,'Done.'

END SUBROUTINE node_db
