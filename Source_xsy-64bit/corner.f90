SUBROUTINE corner(id)
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Subroutine CORNER: reads a triangular mesh in Tecplot format (version 9 &  !
!  later of Tecplot) and processes the corners found such that each triangle  !
!  that has 3 nodes located on boundaries is divided into 2 triangles, thus   !
!  guaranteeing that no triangle has more than 2 nodes over boundaries.       !
!  This code is based on program corner2.                                     !
!                                                                             !
!  INPUT:                                                                     !
!    id  file unit for READs and WRITEs, i.e., the file that contains the     !
!        STORM data.                                                          !
!                                                                             !
!  NOTE: this subroutine does not use any of the MODULEs in SToRM.  All the   !
!  variable TYPEs are defined in the body of the subroutine itself.  This     !
!  avoids variable name conflicts, but care must be taken so that any called  !
!  FUNCTION or SUBROUTINE is equally free of USE statements and conforms      !
!  with all definitions in this file.                                         !
!                                                                             !
!  Francisco Simoes, February 2007                                            !
!  Last updated (mm-dd-yyyy): 02-27-2007 by F. Simões                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Machine precision.
  INTEGER, PARAMETER :: mp = KIND(1.0D0) ! = KIND(1.0) for single precision.

! Points:
  TYPE :: point
    REAL(KIND=mp) :: x,y  ! Point coordinates.
  END TYPE point

! Edges (element sides):
  TYPE :: edge
    INTEGER :: p(2)  ! Points to beginning and end coordinate points.
    INTEGER :: e(2)  ! Elements that contain the edge.  e(1) = 0 when the edge
                     ! is a boundary edge (i.e., belongs to only 1 element).
    INTEGER :: f     ! Frequency of the edge (i.e., how many times it is used).
!   REAL(KIND=mp) :: normal(2)  ! Components of the normal to the edge.
!   REAL(KIND=mp) :: length     ! Length of edge.
  END TYPE

! Computational grids:
  TYPE :: triangle
    INTEGER :: vertex(3)        ! Points to coordinates of vertices.
    INTEGER :: edge(3)          ! Points to edge information (can be < 0).
    INTEGER :: s                ! Sum of the frequency of all the edges.
!   REAL(KIND=mp) :: area       ! Area of element.
  END TYPE

! Dummy arguments.
  INTEGER, INTENT(IN) :: id

! Local variables.
  REAL (KIND=mp), ALLOCATABLE, DIMENSION(:) :: cd,e,h,phiu,phiv,u,v,wse,z
  INTEGER :: i,ierror,j,k,l,m,n,slength
  INTEGER :: npts,nelems,nedges  ! Number of nodes, elements, and edges.
  CHARACTER (LEN=160) :: buffer,cellctr,title(2),varloc
  CHARACTER (LEN=20) :: nelem_char,npts_char
  TYPE(point), ALLOCATABLE, DIMENSION(:) :: nodes
  TYPE(triangle), ALLOCATABLE, DIMENSION(:) :: grid
  TYPE(edge), ALLOCATABLE, DIMENSION(:) :: edges,tmp_edges
  TYPE(edge) :: edg
  LOGICAL :: edge_in_array
  LOGICAL, DIMENSION(11) :: datav
  INTEGER, EXTERNAL :: ccstring
  LOGICAL, EXTERNAL :: cell_centered

  INTERFACE
    INTEGER FUNCTION npoints(varno,datav,npts,nelems)
      INTEGER, INTENT(IN) :: varno,nelems,npts
      LOGICAL, DIMENSION(:), INTENT(IN) :: datav
    END FUNCTION npoints
  END INTERFACE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Read grid dimension and allocate memory space for the needed arrays.  The  !
!  read sequence is identical to that in subroutine tecplot.                  !
!                                                                             !
!-----------------------------------------------------------------------------!

! Read the three header line in the Tecplot file. No Tecplot zones allowed.
  DO i = 1,2
    READ (id,'(A)') title(i)
  END DO
  READ (id,'(A)') buffer

! Find the CELLCENTERED string. The way this is implemented, only one set of []
! is allowed in the Tecplot line defining the number of nodes and elements.
  i = INDEX(buffer,'[')
  j = INDEX(buffer,']')
  IF (i < 1 .OR. j < 1) THEN
    varloc = ''  ! There are no cell-centered values in the input data.
  ELSE
    varloc = buffer(i:j)
  END IF

! Parse buffer: read number of points and number of elements.
! Read number of points.
  i = INDEX(buffer,'N=',.FALSE.)
  j = INDEX(buffer,',',.FALSE.)
  READ (buffer(i+2:j-1),*) npts

! Read number of elements.
  i = INDEX(buffer,'E=',.FALSE.)
  j = INDEX(buffer(i:),',',.FALSE.)
  READ (buffer(i+2:i+j-1),*) nelems

! Allocate space for some of the global variables.  The extra space allocated
! allows the mesh to grow in the corner areas without the danger of exceeding
! the array dimensions.
  ALLOCATE(nodes(npts*2),STAT=ierror)  ! Nodal information.
  IF (ierror /= 0) CALL alloc_err(ierror)
  ALLOCATE(grid(nelems*2),STAT=ierror) ! Element (control volume) connectivity.
  IF (ierror /= 0) CALL alloc_err(ierror)

  n = MAX(nelems,npts)
  ALLOCATE(z(n*2),h(n*2),u(n*2),v(n*2),cd(n*2),wse(n*2),phiu(n*2),phiv(n*2), &
    e(n*2),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  ALLOCATE (edges(nelems*3),STAT=ierror)  ! Edge array.
  IF (ierror /= 0) CALL alloc_err(ierror)

!-----------------------------------------------------------------------------!
!                                                                             !
!  Read data from Tecplot file and initialize data structure.                 !
!                                                                             !
!-----------------------------------------------------------------------------!

  datav = .TRUE.
  DO i = 1,npts  ! Variable #1 in the Tecplot order,
    READ (id,*) nodes(i)%x
  END DO
  DO i = 1,npts  ! variable #2 in the Tecplot order,
    READ (id,*) nodes(i)%y
  END DO
! Variable #3 in the Tecplot order, etc.
  IF (cell_centered('03',varloc)) datav(3) = .FALSE.
  DO i = 1,npoints(3,datav,npts,nelems)
    READ (id,*) z(i)
  END DO
  IF (cell_centered('04',varloc)) datav(4) = .FALSE.
  DO i = 1,npoints(4,datav,npts,nelems)
    READ (id,*) h(i)
  END DO
  IF (cell_centered('05',varloc)) datav(5) = .FALSE.
  DO i = 1,npoints(5,datav,npts,nelems)
    READ (id,*) u(i)
  END DO
  IF (cell_centered('06',varloc)) datav(6) = .FALSE.
  DO i = 1,npoints(6,datav,npts,nelems)
    READ (id,*) v(i)
  END DO
  IF (cell_centered('07',varloc)) datav(7) = .FALSE.
  DO i = 1,npoints(7,datav,npts,nelems)
    READ (id,*) cd(i)
  END DO
  IF (cell_centered('08',varloc)) datav(8) = .FALSE.
  DO i = 1,npoints(8,datav,npts,nelems)
    READ (id,*) wse(i)
  END DO
  IF (cell_centered('09',varloc)) datav(9) = .FALSE.
  DO i = 1,npoints(9,datav,npts,nelems)
    READ (id,*) phiu(i)
  END DO
  IF (cell_centered('10',varloc)) datav(10) = .FALSE.
  DO i = 1,npoints(10,datav,npts,nelems)
    READ (id,*) phiv(i)
  END DO
  IF (cell_centered('11',varloc)) datav(11) = .FALSE.
  DO i = 1,npoints(11,datav,npts,nelems)
    READ (id,*) e(i)
  END DO

! Read element connectivity, counterclockwise.
  DO i = 1,nelems
    READ (id,*) grid(i)%vertex(1),grid(i)%vertex(2),grid(i)%vertex(3)
  END DO

  REWIND id

! Allocate temporary storage, array of maximum possible size for the number of
! triangles in the computational grid, plus 1.  The extra element in the array
! is convenient because it makes programming easier in corner_edge_in_array and
! add_edge_to_array.
  ALLOCATE(tmp_edges(nelems*3 + 1),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)

! Sweep the computational mesh looking for unique edges and store them in array
! tmp_edges.
  nedges = 0
  DO i = 1,nelems
    DO j = 1,3
      edg%p(1) = grid(i)%vertex(j)
      edg%p(2) = grid(i)%vertex(MAX(1,MOD(j+1,4)))
      ! Note: in the original code edge_in_array is a FUNCTION.
      edge_in_array = .FALSE.
      DO k = 1,nedges
        IF ((edg%p(1) == tmp_edges(i)%p(1)  .AND. &
             edg%p(2) == tmp_edges(i)%p(2)) .OR.  &
            (edg%p(2) == tmp_edges(i)%p(1)  .AND. &
             edg%p(1) == tmp_edges(i)%p(2))) THEN
          edge_in_array = .TRUE.
          EXIT
        END IF
      END DO
      IF (.NOT. edge_in_array) THEN
        ! Add edge to edge-array.  This code replaces the call to SUBROUTINE
        ! add_edge_to_array in the original program.
        nedges = nedges + 1
        tmp_edges(nedges) = edg
      END IF
    END DO
  END DO

! Initialize array edges.
  DO i = 1,nedges
    edges(i) = tmp_edges(i)
  END DO

! Release the space allocated to tmp_edges.
  DEALLOCATE(tmp_edges,STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)

! Sweep computational mesh and initialize the grid data structure with the
! information contained in array edges.
  DO i = 1,nelems
    DO j = 1,3
      edg%p(1) = grid(i)%vertex(j)
      edg%p(2) = grid(i)%vertex(MAX(1,MOD(j+1,4)))
      grid(i)%edge(j) = 0  ! Initialize to test for error condition later.

      ! Find the traingle edge (given in edg) in the edge database.
      DO k = 1,nedges
        IF (edg%p(1) == edges(k)%p(1) .AND. &
            edg%p(2) == edges(k)%p(2)) THEN
          ! Edge proceeds counterclockwise in the triangle, i.e., it points in
          ! the outwards direction:
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
             & "SToRM programming and maintenance team at the USGS.",/, &
             & "Program stopped!")')
        STOP
      END IF
    END DO
  END DO

!-----------------------------------------------------------------------------!
!                                                                             !
!  Search for all the wall edges, i.e., all the edges that belong to only     !
!  one triangle.                                                              !
!                                                                             !
!-----------------------------------------------------------------------------!

! Compute the frequency of each edge.
  DO i = 1,nedges
    l = 0
    DO j = 1,nelems
      DO k = 1,3
        IF (ABS(grid(j)%edge(k)) == i) l = l + 1
      END DO
    END DO
    IF (l > 2 .OR. l < 1) THEN
      WRITE (*,'(" ERROR: wrong frequency in edge",I7,I7)') i,l
      WRITE (*,'(" Program CORNER stopped.")')
      STOP
    END IF
    edges(i)%f = l
  END DO

! Compute the sum of all edge frequencies in each triangle.
  DO i = 1,nelems
    grid(i)%s = 0
    DO j = 1,3
      k = ABS(grid(i)%edge(j))
      grid(i)%s = grid(i)%s + edges(k)%f
    END DO
  END DO

! Find the triangles that use each edge.
  DO i = 1,nedges
    edges(i)%e = 0
    DO j = 1,nelems
      DO k = 1,3
        IF (ABS(grid(j)%edge(k)) == i) THEN
          edges(i)%e(1) = j
          EXIT
        END IF
      END DO
      IF (edges(i)%e(1) /= 0) EXIT
    END DO

    IF (edges(i)%f ==2) THEN
      DO j = nelems,1,-1
        DO k = 1,3
          IF (ABS(grid(j)%edge(k)) == i) THEN
            edges(i)%e(2) = j
            EXIT
          END IF
        END DO
        IF (edges(i)%e(2) /= 0) EXIT
      END DO
    END IF

  ! Make sure they are ordered in increasing manner.
    IF (edges(i)%e(1) > edges(i)%e(2)) THEN
      l = edges(i)%e(1)
      edges(i)%e(1) = edges(i)%e(2)
      edges(i)%e(2) = l
    END IF

    IF (edges(i)%e(1) == edges(i)%e(2)) THEN
      WRITE (*,'(" ERROR: wrong count in edge",I7)') i
      WRITE (*,'(" Program CORNER stopped.")')
      STOP
    END IF
  END DO

!-----------------------------------------------------------------------------!
!                                                                             !
!  Find the triangles that contain more than one wall edge and process them:  !
!  two triangles are split into four:                                         !
!                                                                             !
!            +----+           +----+                                          !
!            |\   |           |\  /|                                          !
!            | \  |   ---->   | \/ |                                          !
!            |  \ |           | /\ |                                          !
!            |   \|           |/  \|                                          !
!            +----+           +----+                                          !
!                                                                             !
!-----------------------------------------------------------------------------!

  m = nelems
  DO i = 1,m
    IF (grid(i)%s < 4) THEN
      WRITE (*,'(" ERROR: edge frequency too low in triangle",I7)') i
      WRITE (*,'(" Program CORNER stopped.")')
      STOP
    ELSE IF (grid(i)%s > 6) THEN
      WRITE (*,'(" ERROR: edge frequency too low in triangle",I7)') i
      WRITE (*,'(" Program CORNER stopped.")')
      STOP
    ELSE IF (grid(i)%s /= 4) THEN
      CYCLE
    END IF

    ! Find which edge is not a boundary edge (k).
    DO l = 1,3
      k = ABS(grid(i)%edge(l))
      IF (edges(k)%e(1) /= 0) EXIT
    END DO

    ! Find the adjoining triangle (j).
    j = edges(k)%e(1)
    IF (j == i) j = edges(k)%e(2)

    ! New point for triangle generation, located at the midpoint of the edge.
    npts = npts + 1
    nodes(npts)%x = 0.5e0*(nodes(edges(k)%p(1))%x + nodes(edges(k)%p(2))%x)
    nodes(npts)%y = 0.5e0*(nodes(edges(k)%p(1))%y + nodes(edges(k)%p(2))%y)
    ! Interpolate vertex-based quatities.
    IF (datav(3)) z(npts) = 0.5e0*(z(edges(k)%p(1)) + z(edges(k)%p(2)))
    IF (datav(4)) h(npts) = 0.5e0*(h(edges(k)%p(1)) + h(edges(k)%p(2)))
    IF (datav(5)) u(npts) = 0.5e0*(u(edges(k)%p(1)) + u(edges(k)%p(2)))
    IF (datav(6)) v(npts) = 0.5e0*(v(edges(k)%p(1)) + v(edges(k)%p(2)))
    IF (datav(7)) cd(npts) = 0.5e0*(cd(edges(k)%p(1)) + cd(edges(k)%p(2)))
    IF (datav(8)) wse(npts) = 0.5e0*(wse(edges(k)%p(1)) + wse(edges(k)%p(2)))
    IF (datav(9)) phiu(npts) = 0.5e0*(phiu(edges(k)%p(1)) + &
      phiu(edges(k)%p(2)))
    IF (datav(10)) phiv(npts) = 0.5e0*(phiv(edges(k)%p(1)) + &
      phiv(edges(k)%p(2)))
    IF (datav(11)) e(npts) = 0.5e0*(e(edges(k)%p(1)) + e(edges(k)%p(2)))

    ! Subdivide triangle i and j.
    CALL divide(i,k,grid,nelems,npts)
    ! Copy cell-centered quantities (value of descendents is identical to the
    ! value of the parent cell.
    IF (.NOT.datav(3)) z(nelems) = z(i)
    IF (.NOT.datav(4)) h(nelems) = h(i)
    IF (.NOT.datav(5)) u(nelems) = u(i)
    IF (.NOT.datav(6)) v(nelems) = v(i)
    IF (.NOT.datav(7)) cd(nelems) = cd(i)
    IF (.NOT.datav(8)) wse(nelems) = wse(i)
    IF (.NOT.datav(9)) phiu(nelems) = phiu(i)
    IF (.NOT.datav(10)) phiv(nelems) = phiv(i)
    IF (.NOT.datav(11)) e(nelems) = e(i)
    CALL divide(j,k,grid,nelems,npts)
    IF (.NOT.datav(3)) z(nelems) = z(j)
    IF (.NOT.datav(4)) h(nelems) = h(j)
    IF (.NOT.datav(5)) u(nelems) = u(j)
    IF (.NOT.datav(6)) v(nelems) = v(j)
    IF (.NOT.datav(7)) cd(nelems) = cd(j)
    IF (.NOT.datav(8)) wse(nelems) = wse(j)
    IF (.NOT.datav(9)) phiu(nelems) = phiu(j)
    IF (.NOT.datav(10)) phiv(nelems) = phiv(j)
    IF (.NOT.datav(11)) e(nelems) = e(j)

  END DO

!-----------------------------------------------------------------------------!
!                                                                             !
!  Output the newly processed mesh in Tecplot format.                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Set-up CELLCENTERD string for Tecplot file header.
  slength = ccstring(datav,11,cellctr)

  DO i = 1,2
    WRITE (id,'(A)') TRIM(title(i))
  END DO

  WRITE (npts_char,*) npts
  npts_char = ADJUSTL(npts_char)
  i = LEN_TRIM(npts_char)
  WRITE (nelem_char,*) nelems
  nelem_char = ADJUSTL(nelem_char)
  j = LEN_TRIM(nelem_char)
  IF (slength > 0) THEN
    WRITE (id,'(A)') "ZONE N=" // npts_char(1:i) //", E=" // nelem_char(1:j) &
      // ", DATAPACKING=BLOCK, ZONETYPE=FETRIANGLE" // &
      ", VARLOCATION=([" // TRIM(cellctr) // "]=CELLCENTERED)"
  ELSE
    WRITE (id,'(A)') "ZONE N=" // npts_char(1:i) //", E=" // nelem_char(1:j) &
      // ", DATAPACKING=BLOCK, ZONETYPE=FETRIANGLE"
  END IF

! Write the nodal coordinates and solution information.
  DO i = 1,npoints(1,datav,npts,nelems)  ! X-coordinate.
    WRITE (id,'(ES24.15)') nodes(i)%x
  END DO
  DO i = 1,npoints(2,datav,npts,nelems)  ! Y-coordinate.
    WRITE (id,'(ES24.15)') nodes(i)%y
  END DO
! Bed elevation and other variables, in the proper order.
  DO i = 1,npoints(3,datav,npts,nelems)
    WRITE (id,'(ES24.15)') z(i)
  END DO
  DO i = 1,npoints(4,datav,npts,nelems)
    WRITE (id,'(ES24.15)') h(i)
  END DO
  DO i = 1,npoints(5,datav,npts,nelems)
    WRITE (id,'(1ES24.15)') u(i)
  END DO
  DO i = 1,npoints(6,datav,npts,nelems)
    WRITE (id,'(ES24.15)') v(i)
  END DO
  DO i = 1,npoints(7,datav,npts,nelems)
    WRITE (id,'(ES24.15)') cd(i)
  END DO
  DO i = 1,npoints(8,datav,npts,nelems)
    WRITE (id,'(ES24.15)') wse(i)
  END DO
  DO i = 1,npoints(9,datav,npts,nelems)
    WRITE (id,'(ES24.15)') phiu(i)
  END DO
  DO i = 1,npoints(10,datav,npts,nelems)
    WRITE (id,'(ES24.15)') phiv(i)
  END DO
  DO i = 1,npoints(11,datav,npts,nelems)
    WRITE (id,'(ES24.15)') e(i)
  END DO

! Finally, write the element connectivity table.
  DO i = 1,nelems
    WRITE (id,*) grid(i)%vertex(1),grid(i)%vertex(2),grid(i)%vertex(3)
  END DO

  REWIND id

! Release temporary storage of arrays.
  DEALLOCATE(nodes,grid,z,h,u,v,cd,wse,phiu,phiv,e,edges)

END SUBROUTINE corner
