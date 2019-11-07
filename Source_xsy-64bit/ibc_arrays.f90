SUBROUTINE ibc_arrays
  USE parameters
  USE constants
  USE geometry
  USE dep_vars
  USE options
  USE vbc_arrays
  USE memory
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Calculation of preliminary arrays needed later for implementing the        !
!  inflow boundary conditions.  Needed by both types of boundary conditions   !
!  (i.e., velocity and discharge).  It also prepares the arrays for other     !
!  boundary conditions, such as the free boundary nodes.                      !
!                                                                             !
!  F. Simoes, November 2005                                                   !
!  Last updated (mm-dd-yyyy): 03-21-2011 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables.
  INTEGER :: counter,d1,d2,e,e1,e2,i,i0,ierror,j,k,t,t1,t2
  INTEGER, EXTERNAL :: edge_in_element,is_edge
  REAL (KIND=mp) :: length,n1x,n1y,n2x,n2y
  LOGICAL, EXTERNAL :: pto_in_edge

  INTERFACE
    LOGICAL FUNCTION in_array(i,a)
    INTEGER, INTENT(IN) :: i
    INTEGER, DIMENSION(:), INTENT(IN) :: a
    END FUNCTION in_array
  END INTERFACE

!-----------------------------------------------------------------------------!
!                                                                             !
!               I N F L O W   E D G E S   A N D   N O R M A L S               !
!                                                                             !
!-----------------------------------------------------------------------------!

  d1 = SIZE(qin_nodes,1)  ! Max number of points per boundary.
  d2 = SIZE(qin_nodes,2)  ! Number of inflow boundaries.
  IF (d2 /= n_inflowbdr) THEN
    PRINT *,"ERROR allocating inflow normals,"
    PRINT *,"Dimension does not match count."
    CALL byebye('Program SToRM stopped.')
  END IF

  ALLOCATE(pnormals(d1,d2),bcedgesinflow(0:d1,d2),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(4*(d1 + 1)*d2 + v_size*d1*d2)

  ALLOCATE(l_inflow(d1,d2),h_inflow(d1),a_inflow(d1),u_inflow(d1), &
    v_inflow(d1),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(8*d1*d2 + 4*8*d1)

  DO i0 = 1,n_inflowbdr
    ! Find the inflow boundary edges and normals.
    ! Number of inflow edges for the present boundary condition.
    bcedgesinflow(0,i0) = n_qin(i0) - 1
    DO i = 1,n_qin(i0)-1

      ! Find the edge defined by the two nodes.
      e = is_edge(qin_nodes(i,i0),qin_nodes(i+1,i0))
      IF (e == 0) THEN
        PRINT *,"ERROR in inflow boundary node assignement,"
        PRINT *,"nodes",qin_nodes(i,i0),qin_nodes(i+1,i0)
        PRINT *,"These nodes are not connected by an edge."
        PRINT *,"Please revise your inflow boundary data."
        CALL byebye('Program SToRM stopped.')
      ELSE IF (e < 0) THEN
        PRINT *,"ERROR in inflow boundary node assignement,"
        PRINT *,"nodes",qin_nodes(i,i0)," or ",qin_nodes(i+1,i0),"(or both)."
        PRINT *,"Nodes out of range.  Please revise your inflow boundary data."
        CALL byebye('Program SToRM stopped.')
      END IF
      bcedgesinflow(i,i0) = e

      ! Find the triangle containing the edge.
      t = edge_in_element(edges(e))
      IF (t < 1) THEN
        PRINT *,"ERROR in inflow boundary node assignement,"
        PRINT *,"nodes",qin_nodes(i,i0),qin_nodes(i+1,i0)
        PRINT *,"Edge condition passed; triangle condition failed."
        PRINT *,"Please revise your inflow boundary data."
        CALL byebye('Program SToRM stopped.')
      END IF

      ! Use the information in the computational grid to determine the normal
      ! to the edge that points into the computational domain.  This normal, of
      ! unit length, will determine the direction of the inflow velocity
      ! vectors.
      DO j = 1,3
        k = grid(t)%edge(j)
        IF (e == ABS(k)) THEN
          IF (k > 0) THEN
            ! The edge proceeds counterclockwise in the triangle, therefore it
            ! points in the outwards direction: reverse signs.
            pnormals(i,i0)%x = -edges(e)%normal(1)  ! x-component.
            pnormals(i,i0)%y = -edges(e)%normal(2)  ! y-component.
          ELSE
            ! The edge proceeds clockwise in the triangle, therefore it points
            ! into the triangle, i.e., into the computational domain.
            pnormals(i,i0)%x = edges(e)%normal(1)  ! x-component.
            pnormals(i,i0)%y = edges(e)%normal(2)  ! y-component.
          END IF
          EXIT
        END IF
      END DO

      ! Finally, initialize edge-length array for convenience of use later.
      l_inflow(i,i0) = edges(e)%length

    END DO

    ! All the edge normals have been computed.  Now, transfer that information
    ! to the nodes by averaging from adjacent edge normals.  The normal to node
    ! i = 1 is identical to the normal of the first edge, and the normal to
    ! node i = n_qin(i0) is identical to the normal of the last edge.
    IF (n_qin(i0) > 1) pnormals(n_qin(i0),i0) = pnormals(n_qin(i0)-1,i0)
    DO i = n_qin(i0)-1,2,-1
      pnormals(i,i0)%x = pnormals(i,i0)%x + pnormals(i-1,i0)%x
      pnormals(i,i0)%y = pnormals(i,i0)%y + pnormals(i-1,i0)%y
    END DO

    ! Finally, normalize the vectors to become the unit normal pointing into the
    ! flow  domain, as desired.
    DO i = 1,n_qin(i0)
      length = SQRT(pnormals(i,i0)%x*pnormals(i,i0)%x + &
               pnormals(i,i0)%y*pnormals(i,i0)%y)
      pnormals(i,i0)%x = pnormals(i,i0)%x/length
      pnormals(i,i0)%y = pnormals(i,i0)%y/length
    END DO

  END DO

!-----------------------------------------------------------------------------!
!                                                                             !
!              O U T F L O W   E D G E S   A N D   N O R M A L S              !
!                                                                             !
!-----------------------------------------------------------------------------!

  ! Compute the outward-pointing normals to all the nodes where the stage
  ! boundary condition is enforced.  This also includes free boundaries.

  d1 = SIZE(hbc_nodes,1)  ! Max number of points per boundary.
  d2 = SIZE(hbc_nodes,2)  ! Number of outflow boundaries.
  IF (d2 /= n_outflowbdr) THEN
    PRINT *,"ERROR allocating outflow normals,"
    PRINT *,"Dimension does not match count."
    CALL byebye('Program SToRM stopped.')
  END IF

  ALLOCATE(hnormals(d1,d2),bcedgesoutflow(0:d1,d2),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(4*(d1 + 1)*d2 + v_size*d1*d2)

! Find the strings of outflow boundary edges.
  DO i0 = 1,n_outflowbdr
    bcedgesoutflow(0,i0) = n_hbc(i0) - 1  ! Number of edges for current outlet.
    DO i = 1,n_hbc(i0)-1

      ! Find the edge defined by the two nodes.
      e = is_edge(hbc_nodes(i,i0),hbc_nodes(i+1,i0))
      IF (e == 0) THEN
        PRINT *,"ERROR in outflow boundary node assignement,"
        PRINT *,"nodes",hbc_nodes(i,i0),hbc_nodes(i+1,i0)
        PRINT *,"These nodes are not connected by an edge."
        PRINT *,"Please revise your outflow boundary data."
        CALL byebye('Program SToRM stopped.')
      ELSE IF (e < 0) THEN
        PRINT *,"ERROR in outflow boundary node assignement,"
        PRINT *,"nodes",hbc_nodes(i,i0)," or ",hbc_nodes(i+1,i0),"(or both)."
        PRINT *,"Nodes out of range.  Please revise your outflow boundary data."
        CALL byebye('Program SToRM stopped.')
      END IF
      bcedgesoutflow(i,i0) = e
    END DO
  END DO

! Find the outflow normals.
  DO i0 = 1,n_outflowbdr
    DO i = 1,n_hbc(i0)
      ! For each node, find the two adjacent boundary edges that share that
      ! node.
      counter = 0
      DO j = 1,n_bpolygon
        IF (pto_in_edge(hbc_nodes(i,i0),edges(bpolygon(j)))) THEN
          e2 = bpolygon(j)
          counter = counter + 1
          IF (counter == 2) THEN
            e2 = bpolygon(j)
            EXIT
          END IF
          e1 = bpolygon(j)
        END IF
      END DO

      ! e1 and e2 are the desired edges. Compute the normal at point i by
      ! averaging the outward-pointing normals of these edges.
      t1 = edge_in_element(edges(e1))  ! Find the triangle containing the edge.
      t2 = edge_in_element(edges(e2))

      ! Use the information in the computational grid to determine the normal
      ! to the edge that points out of the computational domain.
      DO j = 1,3
        k = grid(t1)%edge(j)
        IF (e1 == ABS(k)) THEN
          IF (k > 0) THEN
            ! The edge proceeds counterclockwise in the triangle, therefore it
            ! points in the outwards direction, i.e., out of the computational
            ! domain.
            n1x = edges(e1)%normal(1)  ! x-component.
            n1y = edges(e1)%normal(2)  ! y-component.
          ELSE
            ! The edge proceeds clockwise in the triangle, therefore it points
            ! into the triangle: reverse signs.
            n1x = -edges(e1)%normal(1)  ! x-component.
            n1y = -edges(e1)%normal(2)  ! y-component.
          END IF
          EXIT
        END IF
      END DO
      DO j = 1,3
        k = grid(t2)%edge(j)
        IF (e2 == ABS(k)) THEN
          IF (k > 0) THEN
            ! The edge proceeds counterclockwise in the triangle, therefore it
            ! points in the outwards direction, i.e., out of the computational
            ! domain.
            n2x = edges(e2)%normal(1)  ! x-component.
            n2y = edges(e2)%normal(2)  ! y-component.
          ELSE
            ! The edge proceeds clockwise in the triangle, therefore it points
            ! into the triangle: reverse signs.
            n2x = -edges(e2)%normal(1)  ! x-component.
            n2y = -edges(e2)%normal(2)  ! y-component.
          END IF
          EXIT
        END IF
      END DO

      ! Average and normalize.
      hnormals(i,i0)%x = n1x + n2x
      hnormals(i,i0)%y = n1y + n2y
      length = SQRT(hnormals(i,i0)%x*hnormals(i,i0)%x + &
               hnormals(i,i0)%y*hnormals(i,i0)%y)
      hnormals(i,i0)%x = hnormals(i,i0)%x/length
      hnormals(i,i0)%y = hnormals(i,i0)%y/length

    END DO

  END DO

!-----------------------------------------------------------------------------!
!                                                                             !
!                T A N G E N T S   T O   S O L I D   W A L L S                !
!                                                                             !
!-----------------------------------------------------------------------------!

  ! For solid walls with the free-slip condition, compute the vectors at each
  ! wall node that are tangent to the wall.
  !IF (opt_solver == 1 .AND. btype == 1 .AND. n_wall > 0) THEN
  IF (n_wall > 0) THEN
    ALLOCATE(wtangs(n_wall),STAT=ierror)
    IF (ierror /= 0) CALL alloc_err(ierror)
    CALL mem_add(v_size*n_wall)

    ! The procedure is the following: first compute all the normals at each
    ! node, just like it is done for the other types of boundary nodes.  Then,
    ! use those normals to compute the tangents simply by coordinate inversion.
    DO i = 1,n_wall
      ! For each node, find the two adjacent boundary edges that share that
      ! node.
      counter = 0
      DO j = 1,n_bpolygon
        IF (pto_in_edge(wall_pts(i),edges(bpolygon(j)))) THEN
          e2 = bpolygon(j)
          counter = counter + 1
          IF (counter == 2) THEN
            e2 = bpolygon(j)
            EXIT
          END IF
          e1 = bpolygon(j)
        END IF
      END DO

      ! e1 and e2 are the desired edges. Compute the normal at point i by
      ! averaging the outward-pointing normals of these edges.
      t1 = edge_in_element(edges(e1))  ! Find the triangle containing the edge.
      t2 = edge_in_element(edges(e2))

      ! Use the information in the computational grid to determine the normal
      ! to the edge that points out of the computational domain.
      DO j = 1,3
        k = grid(t1)%edge(j)
        IF (e1 == ABS(k)) THEN
          IF (k > 0) THEN
            ! The edge proceeds counterclockwise in the triangle, therefore it
            ! points in the outwards direction, i.e., out of the computational
            ! domain.
            n1x = edges(e1)%normal(1)  ! x-component.
            n1y = edges(e1)%normal(2)  ! y-component.
          ELSE
            ! The edge proceeds clockwise in the triangle, therefore it points
            ! into the triangle: reverse signs.
            n1x = -edges(e1)%normal(1)  ! x-component.
            n1y = -edges(e1)%normal(2)  ! y-component.
          END IF
          EXIT
        END IF
      END DO
      DO j = 1,3
        k = grid(t2)%edge(j)
        IF (e2 == ABS(k)) THEN
          IF (k > 0) THEN
            ! The edge proceeds counterclockwise in the triangle, therefore it
            ! points in the outwards direction, i.e., out of the computational
            ! domain.
            n2x = edges(e2)%normal(1)  ! x-component.
            n2y = edges(e2)%normal(2)  ! y-component.
          ELSE
            ! The edge proceeds clockwise in the triangle, therefore it points
            ! into the triangle: reverse signs.
            n2x = -edges(e2)%normal(1)  ! x-component.
            n2y = -edges(e2)%normal(2)  ! y-component.
          END IF
          EXIT
        END IF
      END DO

      ! Average and normalize.  These are the normals pointing outwards.
      wtangs(i)%x = n1x + n2x
      wtangs(i)%y = n1y + n2y
      length = SQRT(wtangs(i)%x*wtangs(i)%x + wtangs(i)%y*wtangs(i)%y)
      wtangs(i)%x = wtangs(i)%x/length
      wtangs(i)%y = wtangs(i)%y/length
      ! Now invert the coordinates to get the tangents.  Not that it matters,
      ! but the tangents are pointing in the counterclockwise direction...
      length = wtangs(i)%x
      wtangs(i)%x = -wtangs(i)%y
      wtangs(i)%y = length

    END DO

  END IF

  ! For solid walls with the free-slip condition, compute the vectors at each
  ! wall node that are tangent to the wall.
  !IF (opt_solver == 2 .AND. btype == 1 .AND. wall_edges1 > 0) THEN
  IF (wall_edges1 > 0) THEN
    ALLOCATE(wtangs1(wall_edges1),STAT=ierror)
    IF (ierror /= 0) CALL alloc_err(ierror)
    CALL mem_add(v_size*wall_edges1)

    ! The process is simple: rotate the edge normals by 90 degrees.
    DO i = 1,wall_edges1
      j = walledg1(i)
      wtangs1(i)%x = -edges(j)%normal(2)/edges(j)%length
      wtangs1(i)%y = edges(j)%normal(1)/edges(j)%length
    END DO

  END IF

  ! CLARIFICATION: array wtangs() is defined on wall nodes, while array
  ! wtangs1() is defined on wall edges (which are specified in array
  ! walledg1()).

!-----------------------------------------------------------------------------!
!                                                                             !
! L O G I C A L   A R R A Y S   F O R   B O U N D A R Y   C O N D I T I O N S !
!                                                                             !
!-----------------------------------------------------------------------------!

  n_bcedges = n_bpolygon - wall_edges1
  ALLOCATE(bcedges(n_bcedges),hbar(n_bcedges),qbar(n_bcedges),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(4*3*n_bcedges)
! Do some error checking.
  k = 0
  DO i0 = 1,n_inflowbdr
    k = k + n_qin(i0) - 1
  END DO
  DO i0 = 1,n_outflowbdr
    k = k + n_hbc(i0) - 1
  END DO
  IF (k /= n_bcedges) THEN
    PRINT *,"ERROR: number of boundary condition edges does not"
    PRINT *,"match number of boundary condition chains."
    CALL byebye('Program SToRM stopped.')
  END IF
! hbar(i) is .TRUE. when a stage boundary condition is imposed on edge i and
! .FALSE. otherwise; qbar(i) is .TRUE. when a velocity or discharge boundary
! condition is imposed on edge i and .FALSE. otherwise.
  hbar = .FALSE.
  qbar = .FALSE.

  counter = 0
  DO i = 1,n_bpolygon
    j = bpolygon(i)
    IF (.NOT.in_array(j,walledg1)) THEN
      counter = counter + 1
      bcedges(counter) = j
    END IF
  END DO
  IF (counter /= n_bcedges) THEN
    PRINT *,"ERROR in inflow boundary conditions: wall edges"
    PRINT *,"and inflow/outflow edges do not add-up right."
    PRINT *,"Please revise your inflow/outflow boundary data."
    CALL byebye('Program SToRM stopped.')
  END IF

! Find the inflow boundary edges.
  DO i0 = 1,n_inflowbdr
    DO i = 1,n_qin(i0)-1

      ! Find the edge defined by the two nodes.
      e = is_edge(qin_nodes(i,i0),qin_nodes(i+1,i0))
      IF (e == 0) THEN
        PRINT *,"ERROR in inflow boundary node assignement,"
        PRINT *,"nodes",qin_nodes(i,i0),qin_nodes(i+1,i0)
        PRINT *,"These nodes are not connected by an edge."
        PRINT *,"Please revise your stage boundary data."
        CALL byebye('Program SToRM stopped.')
      ELSE IF (e < 0) THEN
        PRINT *,"ERROR in stage boundary node assignement,"
        PRINT *,"nodes",qin_nodes(i,i0)," or ",qin_nodes(i+1,i0),"(or both)."
        PRINT *,"Nodes out of range.  Please revise your stage boundary data."
        CALL byebye('Program SToRM stopped.')
      END IF
      DO j = 1,n_bcedges
        IF (e == bcedges(j)) THEN
          qbar(j) = .TRUE.
          EXIT
        END IF
      END DO
    END DO
  END DO

! Find the stage boundary edges.
  DO i0 = 1,n_outflowbdr
    DO i = 1,n_hbc(i0)-1

      e = is_edge(hbc_nodes(i,i0),hbc_nodes(i+1,i0))
      IF (e == 0) THEN
        PRINT *,"ERROR in stage boundary node assignement,"
        PRINT *,"nodes",hbc_nodes(i,i0),hbc_nodes(i+1,i0)
        PRINT *,"These nodes are not connected by an edge."
        PRINT *,"Please revise your stage boundary data."
        CALL byebye('Program SToRM stopped.')
      ELSE IF (e < 0) THEN
        PRINT *,"ERROR in stage boundary node assignement,"
        PRINT *,"nodes",hbc_nodes(i,i0)," or ",hbc_nodes(i+1,i0),"(or both)."
        PRINT *,"Nodes out of range.  Please revise your stage boundary data."
        CALL byebye('Program SToRM stopped.')
      END IF
      DO j = 1,n_bcedges
        IF (e == bcedges(j)) THEN
          hbar(j) = .TRUE.
          EXIT
        END IF
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------------!
!                                                                             !
!                  B O U N D A R Y   E D G E   N O R M A L S                  !
!                                                                             !
!-----------------------------------------------------------------------------!

! The boundary edge normals are unit-vectors that are normal to all the edges
! where non-wall boundary conditions are applied, and that point into the
! computational domain.
!  ALLOCATE(abcnormals(n_bcedges),STAT=ierror)
!  IF (ierror /= 0) CALL alloc_err(ierror)
!  CALL mem_add(v_size*n_bcedges)

!  DO i = 1,n_bcedges
!    e = bcedges(i)
!    t = edge_in_element(edges(e))
!    DO j = 1,3
!      k = grid(t)%edge(j)
!      IF (e == ABS(K)) THEN
!        IF (k > 0) THEN
!          ! The edge proceeds counterclockwise in the triangle, therefore it
!          ! points in the outwards direction: reverse signs.
!          abcnormals(i)%x = -edges(e)%normal(1)/edges(e)%length  ! x-component.
!          abcnormals(i)%y = -edges(e)%normal(2)/edges(e)%length  ! y-component.
!        ELSE
!          ! The edge proceeds clockwise in the triangle, therefore it points
!          ! into the triangle, i.e., into the computational domain.
!          abcnormals(i)%x = edges(e)%normal(1)/edges(e)%length  ! x-component.
!          abcnormals(i)%y = edges(e)%normal(2)/edges(e)%length  ! y-component.
!        END IF
!        EXIT
!      END IF
!    END DO
!  END DO

END SUBROUTINE ibc_arrays
