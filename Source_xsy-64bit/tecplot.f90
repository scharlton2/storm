SUBROUTINE tecplot(id)
  USE parameters
  USE geometry
  USE dep_vars
  USE memory
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Read in nodal coordinates and element connectivity from an external data   !
!  file in Tecplot format.  The bed elevation information is also read from   !
!  this file.  The input data format follows that of chapter 19 of Tecplot    !
!  User's Manual, Version 10, September, 2003.  Triangle grids only, no       !
!  Tecplot zones allowed.  Note that variables nodes() and grid() are         !
!  allocated here.  Temporary arrays are used to store , h(), u(), v(), and   !
!  z().  The position of the variables may be at cell centers or at cell      !
!  vertices.                                                                  !
!                                                                             !
!  Francisco Simoes, March 2004                                               !
!  Last updated (mm-dd-yyyy): 08-14-2007 by F. Simões                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments:
  INTEGER, INTENT(IN) :: id  ! Unit number for the external READ statements.

! Local variables:
  INTEGER :: i,ierror,j,n
  CHARACTER (LEN=160) :: buffer,varloc

! Check flag requesting the processing of corners.
  !READ (id,'(A)') buffer
  !REWIND id
  !IF (INDEX(buffer,'::CORN') /= 0) CALL corner(id)

! Check flag requesting renumbering of nodes and elements.
  !READ (id,'(A)') buffer
  !REWIND id
  !IF (INDEX(buffer,'::REORD') /= 0) CALL reorder(id,table1,table2,...)

! Read the three header line in the Tecplot file. No Tecplot zones allowed.
  DO i = 1,3
    READ (id,'(A)') buffer
  END DO

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
  READ (buffer(i+2:j-1),*) n_pts

! Read number of elements.
  i = INDEX(buffer,'E=',.FALSE.)
  j = INDEX(buffer(i:),',',.FALSE.)
  READ (buffer(i+2:i+j-1),*) n_elems

! Allocate space for some of the global variables.
  ALLOCATE(nodes(n_pts),STAT=ierror)  ! Nodal information.
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(p_size*n_pts)
  ALLOCATE(grid(n_elems + 1),STAT=ierror)  ! Element (ctrl vol) connectivity.
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(t_size*(n_elems + 1))

  n = MAX(n_elems,n_pts)
  ALLOCATE(z_tcp(n),STAT=ierror)  ! Bed elevation.
  IF (ierror /= 0) CALL alloc_err(ierror)
  ALLOCATE(h_tcp(n),u_tcp(n),v_tcp(n),STAT=ierror)  ! Solution variables.
  IF (ierror /= 0) CALL alloc_err(ierror)
  ALLOCATE(cd_tcp(n),STAT=ierror)  ! Bed friction coefficient.
  IF (ierror /= 0) CALL alloc_err(ierror)

  call tecplot_read(id,varloc)

END SUBROUTINE tecplot
