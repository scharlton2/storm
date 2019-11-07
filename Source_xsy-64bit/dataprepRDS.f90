SUBROUTINE dataprepRDS
  USE geometry
  USE dep_vars
  USE options
  USE memory
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Prepares the data for use by the RDS solver.  It allocates needed array    !
!  space, moves the variables that need it from cell centers to cell          !
!  vertices, and releases memory used for temporary arrays.                   !
!                                                                             !
!  Francisco Simoes, February 2007                                            !
!  Last updated (mm-dd-yyyy): 02-02-2009 by F. Simões                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables.
  INTEGER :: i,ierror

! Allocate space for some of the global variables.
  ALLOCATE(z(n_pts),zbx(n_pts),zby(n_pts),STAT=ierror)  ! Bed elevation.
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(8*3*n_pts)
  ALLOCATE(h(n_pts),u(n_pts),v(n_pts),STAT=ierror)  ! Solution variables.
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(8*3*n_pts)
  ALLOCATE(cd(n_pts),rough(n_pts),STAT=ierror)  ! Bed friction coefficient.
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(8*2*n_pts)
  ALLOCATE(delh(n_pts),delu(n_pts),delv(n_pts),STAT=ierror)  ! Grads.
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(v_size*3*n_pts)

! Allocate the working variables.
  ALLOCATE(residual(n_elems,3),phi(n_pts,3),flux(n_edges,3), &
           source(n_elems,3),e_est(n_elems),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  ALLOCATE(u_avg(n_elems,3),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)

! Set-up arrays for interpolation from cells to nodes using a least squares
! technique.
  CALL c2v_setup(nctr2vtx)

! For the data that needs it, interpolate from cell-centered to cell-vertex
! positions.
  IF (.NOT.datav(3)) THEN
    CALL ctr2vtx(nctr2vtx,z_tcp,z)
  ELSE
    DO i = 1,n_pts
      z(i) = z_tcp(i)
    END DO
  ENDIF
  IF (.NOT.datav(4)) THEN
    CALL ctr2vtx(nctr2vtx,h_tcp,h)
  ELSE
    DO i = 1,n_pts
      h(i) = h_tcp(i)
    END DO
  ENDIF
  IF (.NOT.datav(5)) THEN
    CALL ctr2vtx(nctr2vtx,u_tcp,u)
  ELSE
    DO i = 1,n_pts
      u(i) = u_tcp(i)
    END DO
  ENDIF
  IF (.NOT.datav(6)) THEN
    CALL ctr2vtx(nctr2vtx,v_tcp,v)
  ELSE
    DO i = 1,n_pts
      v(i) = v_tcp(i)
    END DO
  ENDIF
  IF (.NOT.datav(7)) THEN
    CALL ctr2vtx(nctr2vtx,cd_tcp,cd)
  ELSE
    DO i = 1,n_pts
      cd(i) = cd_tcp(i)
    END DO
  ENDIF

! Free the memory used for the Tecplot temporary arrays.
  DEALLOCATE(z_tcp,h_tcp,u_tcp,v_tcp,cd_tcp)

! Reset datav() so that it can be used to write to Tecplot-formatted files.
  datav = .TRUE.
  datav(11) = .FALSE.

! Set-up arrays for computation of gradients by least-squares reconstruction.
  CALL lsq_setup(opt_solver,2)

END SUBROUTINE dataprepRDS
