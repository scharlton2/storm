SUBROUTINE tecplot(id)
  USE parameters
  USE geometry
  USE dep_vars
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Read in nodal coordinates and element connectivity from an external data   !
!  file in Tecplot format.  The bed elevation information is also read from   !
!  this file.  The input data format follows that of chapter 19 of Tecplot    !
!  User's Manual, Version 10, September, 2003.  Triangle grids only, no       !
!  Tecplot zones allowed.                                                     !
!  Note that variables nodes(), grid(), h(), u(), v(), and z() are allocated  !
!  here.                                                                      !
!                                                                             !
!  Francisco Simoes, March 2004                                               !
!  Last updated (mm-dd-yyyy): 04-02-2004                                      !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments:
  INTEGER, INTENT(IN) :: id  ! Unit number for the external READ statements.

! Local variables:
  INTEGER :: i,ierror,j
  CHARACTER (LEN=80) :: buffer

! Read the three header line in the Tecplot file. No Tecplot zones allowed.
  DO i = 1,3
    READ (id,'(A)') buffer
  END DO

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
  ALLOCATE(grid(n_elems),STAT=ierror)  ! Element (control volume) connectivity.
  IF (ierror /= 0) CALL alloc_err(ierror)
  ALLOCATE(z(n_pts),STAT=ierror)  ! Bed elevation.
  IF (ierror /= 0) CALL alloc_err(ierror)
  ALLOCATE(h(n_pts),u(n_pts),v(n_pts),STAT=ierror)  ! Solution variables.
  IF (ierror /= 0) CALL alloc_err(ierror)
  ALLOCATE(cd(n_pts),rough(n_pts),STAT=ierror)  ! Bed friction coefficient.
  IF (ierror /= 0) CALL alloc_err(ierror)

! Read point coordinates, bed elevation (x,y,z), and friction coefficients.  An
! initial solution must be included here, too, but it may be set to zero if
! necessary.
  DO i = 1,n_pts
    READ (id,*) nodes(i)%x,nodes(i)%y,z(i),h(i),u(i),v(i),cd(i)
  END DO

! Read element connectivity, counterclockwise.
  DO i = 1,n_elems
    READ (id,*) grid(i)%vertex(1),grid(i)%vertex(2),grid(i)%vertex(3)
  END DO

END SUBROUTINE tecplot
