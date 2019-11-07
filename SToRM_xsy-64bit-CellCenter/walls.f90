SUBROUTINE walls
  USE geometry
  USE dep_vars
  USE options
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Sweep of the entire grid to find which edges are wall edges and which      !
!  ones are in the flow.  Wall edges are all the outside edges, i.e., all of  !
!  those that belong to one triangular element only.                          !
!                                                                             !
!  Francisco Simoes, November 2004                                            !
!  Last updated (mm-dd-yyyy): 02-22-2011 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables.
  INTEGER :: edge_count(n_edges),i,ierror,j,k,l,m,n,p,p1,p2
  INTEGER, ALLOCATABLE, DIMENSION(:) :: work
  INTEGER, EXTERNAL :: edge_in_element

  INTERFACE
    LOGICAL FUNCTION in_array(i,a)
    INTEGER, INTENT(IN) :: i
    INTEGER, DIMENSION(:), INTENT(IN) :: a
    END FUNCTION in_array

    LOGICAL FUNCTION in_array2(i,a)
    INTEGER, INTENT(IN) :: i
    INTEGER, DIMENSION(:,:), INTENT(IN) :: a
    END FUNCTION in_array2
  END INTERFACE

  WRITE (*,'("  Initialize boundary arrays...",$)')

! Count how many times each element edge is used by an element. Wall edges are
! those that are used by only one element and inner edges are those used by two
! elements.
  edge_count = 0
  DO n = 1,n_elems
    DO j = 1,3
      i = ABS(grid(n)%edge(j))
      edge_count(i) = edge_count(i) + 1
    END DO
  END DO

! Create separate counters for wall edges and for inner (flow) edges.
  wall_edges = 0  ! Number of wall edges.
  wall_edges1 = 0  ! Number of wall edges.
  flow_edges = 0  ! Number of in-flow edges.
  flow_edges1 = 0  ! Number of in-flow edges,
  n_bpolygon = 0  ! Number of edges in bounding polygon.
  DO i = 1,n_edges
    IF (edge_count(i) == 1) THEN
      p1 = edges(i)%p(1)
      p2 = edges(i)%p(2)
      n_bpolygon = n_bpolygon + 1
      ! Do not include edges that have an inflow, outflow, or free-flow node.
      IF (.NOT.(in_array2(p1,qin_nodes) .OR. in_array2(p1,hbc_nodes)) .AND. &
          .NOT.(in_array2(p2,qin_nodes) .OR. in_array2(p2,hbc_nodes))) THEN
        wall_edges = wall_edges + 1
      ELSE  ! Edge is in boundary, but the fluxes are computed later...
        flow_edges = flow_edges + 1
      END IF
      ! The edges in array walledg1 may contain one node that is an inflow,
      ! outflow, or free-flow node.
      IF (.NOT.((in_array2(p1,qin_nodes) .OR. in_array2(p1,hbc_nodes)) .AND. &
          (in_array2(p2,qin_nodes) .OR. in_array2(p2,hbc_nodes)))) THEN
        wall_edges1 = wall_edges1 + 1
      END IF
    ELSE IF (edge_count(i) == 2) THEN
      flow_edges = flow_edges + 1
      flow_edges1 = flow_edges1 + 1
    ELSE
      PRINT *,'ERROR: computational grid construction error, in edge'
      PRINT *,'assignement (wall or inner edge). Please check the grid'
      PRINT *,'for compliance with the definitions.'
      CALL byebye('SToRM stopped.')
    END IF
  END DO

  IF (wall_edges + flow_edges /= n_edges) THEN  ! Another error check.
    PRINT *,wall_edges,flow_edges,n_edges
    PRINT *,'ERROR: computational grid construction error, in edge'
    PRINT *,'assignement (wall or inner edge). Please check the grid'
    PRINT *,'for compliance with the definitions.'
    CALL byebye('SToRM stopped.')
  END IF

! Allocate space for global arrays.  Make sure that there is at least one
! element in array to prevent failure of the ALLOCATE command.
  wall_edges = MAX(1,wall_edges)
  wall_edges1 = MAX(1,wall_edges1)
  flow_edges = MAX(1,flow_edges)
  flow_edges1 = MAX(1,flow_edges1)
  ALLOCATE(walledg(wall_edges),walledg1(wall_edges1),flowedg(flow_edges), &
    flowedg1(flow_edges1),bpolygon(n_bpolygon),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(4*wall_edges + 4*wall_edges1 + 4*flow_edges + 4*flow_edges1 + &
    4*n_bpolygon)

! Initialize each array with the corresponding pointers.  Each element of the
! array walledg() contains a pointer to an edge that is a solid wall.  No two
! elements of this array contain duplicate pointers, i.e., each element points
! to a unique edge, therefore there are as many wall edges as there are
! elements in the array.  Similarly, array flowedg() points to the edges that
! are flow edges.  Also, bpolygon() points to all the edges that form the
! bounding polygon to the computational grid.
  j = 0
  k = 0
  l = 0
  m = 0
  n = 0
  DO i = 1,n_edges
    IF (edge_count(i) == 1) THEN
      p1 = edges(i)%p(1)
      p2 = edges(i)%p(2)
      l = l + 1
      bpolygon(l) = i
      ! Skip edges that have an inflow, outflow, or free-flow node.
      IF (.NOT.(in_array2(p1,qin_nodes) .OR. in_array2(p1,hbc_nodes)) .AND. &
          .NOT.(in_array2(p2,qin_nodes) .OR. in_array2(p2,hbc_nodes))) THEN
        j = j + 1
        walledg(j) = i
      ELSE  ! Boundary edges with computable fluxes.
        k = k + 1
        flowedg(k) = i
      END IF
      ! The edges in array walledg1 may contain one node that is an inflow,
      ! outflow, or free-flow node.
      IF (.NOT.((in_array2(p1,qin_nodes) .OR. in_array2(p1,hbc_nodes)) .AND. &
          (in_array2(p2,qin_nodes) .OR. in_array2(p2,hbc_nodes)))) THEN
        m = m + 1
        walledg1(m) = i
      END IF
    ELSE IF (edge_count(i) == 2) THEN
      k = k + 1
      flowedg(k) = i
      n = n + 1
      flowedg1(n) = i
    ELSE
      PRINT *,'ERROR: computational grid construction error, in edge'
      PRINT *,'assignement (wall or inner edge). Please check the grid'
      PRINT *,'for compliance with the definitions.'
      CALL byebye('SToRM stopped.')
    END IF
  END DO
! IMPORTANT NOTE ABOUT THE WALL AND FLOW EDGE ARRAYS.  The difference between
! arrays flowedg() and flowedg1() is that the first contains edges in the
! boundary (i.e., inflow, outflow, and free-flow edges), while the latter only
! contains edges inside the computational domain.  The difference between
! arrays walledg() and walledge1() is that none of the end points in walledge()
! are inflow, outflow, or free-flow nodes, while the edges in walledge1() may
! contain one end point that is an inflow, outflow, or free-flow node (but not
! both nodes).  Arrays flowedg1() and walledg1() were developed for the FVT
! solver, in order to facilitate implementation of the boundary conditions.

! Last error check (redundant).
  IF (j /= wall_edges .OR. k /= flow_edges .OR. l /= n_bpolygon .OR. &
      m /= wall_edges1 .OR. n /= flow_edges1 .OR. &
      wall_edges + flow_edges /= n_edges) THEN
    PRINT *,'ERROR: computational grid construction error, in edge'
    PRINT *,'assignement (wall or inner edge). Please check the grid'
    PRINT *,'for compliance with the definitions.'
    CALL byebye('Program SToRM stopped.')
  END IF

! Now find all the nodes that are located on solid boundaries.
  IF (wall_edges > 0) THEN
    ALLOCATE(work(wall_edges*2),STAT=ierror)  ! Working array.
    IF (ierror /= 0) CALL alloc_err(ierror)
    work = 0
    n_wall = 0  ! Number of points located on solid walls.

    ! Find all the wall nodes, without repetitions.
    DO i = 1,wall_edges1
      DO j = 1,2
        p = edges(walledg1(i))%p(j)
        IF ((.NOT. in_array(p,work)) .AND. &  ! Check if point p is in any of
            (.NOT. in_array2(p,qin_nodes)) .AND. &  ! the boundary conditions
            (.NOT. in_array2(p,hbc_nodes))) THEN    ! arrays.
          n_wall = n_wall + 1
          work(n_wall) = p
        END IF
      END DO
    END DO

    ! Allocate and initialize the appropriate array.
    ALLOCATE(wall_pts(n_wall),STAT=ierror)
    IF (ierror /= 0) CALL alloc_err(ierror)
    CALL mem_add(4*n_wall)
    DO i = 1,n_wall
      wall_pts(i) = work(i)
    END DO
    DEALLOCATE(work)
  END IF

  PRINT *,'Done.'

END SUBROUTINE walls
