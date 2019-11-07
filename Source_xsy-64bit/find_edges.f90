SUBROUTINE find_edges
  USE parameters
  USE geometry
  USE memory
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This subroutine initializes the edge data structures in the SToRM data     !
!  base.  Variable edges is allocated here.                                   !
!                                                                             !
!  Francisco Simoes, March 2004                                               !
!  Last updated (mm-dd-yyyy): 08-13-2007                                      !
!                                                                             !
!-----------------------------------------------------------------------------!


! Local variables:
  INTEGER :: i,ierror,j,k,l
  INTEGER, ALLOCATABLE, DIMENSION(:) :: counter
  REAL(KIND=mp) :: delx,dely
  TYPE(edge) :: edg
  TYPE(edge), ALLOCATABLE, DIMENSION (:) :: tmp_edges
  LOGICAL :: tdata
  LOGICAL, EXTERNAL :: edge_in_array

  WRITE (*,'("  Building edge database...",$)')

! First build the array with the edge information.  Note that it is possible
! that this information may have already be read in subroutine 'tecplot_read',
! therefore a check must be done.  Here, this is accomplished by checking if
! array 'edges' has been previously allocated.
  tdata = .TRUE.  ! Indicates if 'triangle' data is read from tecplot file.
  IF (.NOT. ALLOCATED(edges)) THEN
    tdata = .FALSE.

    ! Allocate temporary storage, array of maximum possible size for the number
    ! of triangles in the computational grid, plus 1.  The extra element in the
    ! array is convenient because it makes programming easier in edge_in_array
    ! and add_edge_to_array.
    ALLOCATE(tmp_edges(n_elems*3 + 1),STAT=ierror)
    IF (ierror /= 0) CALL alloc_err(ierror)

    ! Sweep the computational mesh looking for unique edges and store them in
    ! array tmp_edges.
    n_edges = 0
    DO i = 1,n_elems
      DO j = 1,3
        edg%p(1) = grid(i)%vertex(j)
        edg%p(2) = grid(i)%vertex(MAX(1,MOD(j+1,4)))
        IF (.NOT. edge_in_array(edg,tmp_edges,n_edges)) &
          CALL add_edge_to_array(edg,tmp_edges,n_edges)
      END DO
    END DO

    ! Allocate storage for global array edges.
    ALLOCATE(edges(n_edges),STAT=ierror)
    IF (ierror /= 0) CALL alloc_err(ierror)
    CALL mem_add(e_size*n_edges)

    ! Initialize array edges.
    DO i = 1,n_edges
      edges(i) = tmp_edges(i)
    END DO

    ! Release the space allocated to tmp_edges.
    DEALLOCATE(tmp_edges,STAT=ierror)
    IF (ierror /= 0) CALL alloc_err(ierror)
  END IF

  IF (tdata) THEN
    ! Edge-to-element connectivity already exists and can be used to find the
    ! element-to-edge connectivity.
    DO i = 1,n_elems
      DO j = 1,3
        grid(i)%edge(j) = 0  ! Initialize to test for error condition later.
      END DO
    END DO

    DO j = 1,n_edges
      DO k = 1,2
        i = edges(j)%e(k)
        IF (i < 0) CYCLE
        DO l = 1,3
          edg%p(1) = grid(i)%vertex(l)
          edg%p(2) = grid(i)%vertex(MAX(1,MOD(l+1,4)))
          IF (edg%p(1) == edges(j)%p(1) .AND. &
              edg%p(2) == edges(j)%p(2)) THEN
            ! Edge proceeds counterclockwise in the triangle, i.e., it points
            ! in the outwards direction:
            grid(i)%edge(l) = j
            EXIT
          ELSE IF (edg%p(1) == edges(j)%p(2) .AND. &
                   edg%p(2) == edges(j)%p(1)) THEN
            ! Edge proceeds clockwise in the triangle, i.e., it points into the
            ! element:
            grid(i)%edge(l) = -j
            EXIT
          END IF
        END DO
      END DO
    END DO

    ! Test for error.
    DO i = 1,n_elems
      DO j = 1,3
        IF (grid(i)%edge(j) == 0) THEN  ! ERROR: edge was not initialized.
          WRITE (*,'("PROGRAM ERROR: edge data structure incomplete.",/, &
               & "This is a fatal programming error.  Please contact the",/, &
               & "SToRM programming and maintenance team at the USGS.")')
          CALL byebye('Program stopped.')
        END IF
      END DO
    END DO

  ELSE
    ! Start from scratch: sweep computational mesh and initialize the grid data
    ! structure with the information contained in array edges.  The sequence in
    ! which edges are set in grid()%edge() is important: each node i must be
    ! followed by edge i so that the edge opposite to node i can easily be
    ! found from MAX(1,MOD(i+1,4)).
    DO i = 1,n_elems
      DO j = 1,3
        edg%p(1) = grid(i)%vertex(j)
        edg%p(2) = grid(i)%vertex(MAX(1,MOD(j+1,4)))
        grid(i)%edge(j) = 0  ! Initialize to test for error condition later.

        ! Find the triangle edge (given in edg) in the edge database.
        DO k = 1,n_edges
          IF (edg%p(1) == edges(k)%p(1) .AND. &
              edg%p(2) == edges(k)%p(2)) THEN
            ! Edge proceeds counterclockwise in the triangle, i.e., it points
            ! in the outwards direction:
            grid(i)%edge(j) = k
            EXIT
          ELSE IF (edg%p(1) == edges(k)%p(2) .AND. &
                   edg%p(2) == edges(k)%p(1)) THEN
            ! Edge proceeds clockwise in the triangle, i.e., it points into the
            ! element:
            grid(i)%edge(j) = -k
            EXIT
          END IF
        END DO

        IF (grid(i)%edge(j) == 0) THEN  ! ERROR: edge was not initialized.
          WRITE (*,'("PROGRAM ERROR: edge data structure incomplete.",/, &
               & "This is a fatal programming error.  Please contact the",/, &
               & "SToRM programming and maintenance team at the USGS.")')
          CALL byebye('Program stopped.')
        END IF
      END DO
    END DO
  END IF

! Compute all the edge normals.
  CALL edge_normals(edges,n_edges,nodes,n_pts)

! Compute all the edge lengths.
  DO k = 1,n_edges
    delx = nodes(edges(k)%p(1))%x - nodes(edges(k)%p(2))%x
    dely = nodes(edges(k)%p(1))%y - nodes(edges(k)%p(2))%y
    edges(k)%length = SQRT(delx*delx + dely*dely)
  END DO

  IF (.NOT. tdata) THEN
    ! Initialize the edge-to-element connectivity tables.
    DO i = 1,n_edges  ! Initialize pointers to the "no element" state.
      edges(i)%e(1) = -1
      edges(i)%e(2) = -1
    END DO
    ALLOCATE(counter(n_edges),STAT=ierror)  ! Working array.
    IF (ierror /= 0) CALL alloc_err(ierror)
    counter = 1
    DO i = 1,n_elems
      DO j = 1,3
        k = ABS(grid(i)%edge(j))  ! Edge being processed.
        IF (counter(k) > 2) THEN  ! Check for connectivity error.
          WRITE (*,'("FATAL PROGRAM ERROR: edge data structure incorrect.",/, &
               & "Please check your mesh connectivity data.")')
          DEALLOCATE(counter)
          CALL byebye('Program stopped.')
        END IF

        edges(k)%e(counter(k)) = i  ! %e(1) is always initialized before %e(2).
        counter(k) = counter(k) + 1
      END DO
    END DO
    DEALLOCATE(counter)
    ! NOTES: counter(i) keeps track of which element pointer is set at any
    ! particular time. I.e., counter(k) = 1 means that edges(k)%e(1) is set to
    ! the element number in question (element i in the "DO i = 1,n_elems"
    ! loop).
  END IF

  PRINT *,'Done.'

END SUBROUTINE find_edges
