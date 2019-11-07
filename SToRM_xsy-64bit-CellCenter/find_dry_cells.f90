SUBROUTINE find_dry_cells
  USE constants
  USE geometry
  USE dep_vars
  USE options
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  A number of tasks are performed here, but the main objective of this code  !
!  is to find the dry cells and initialize arra 'drycell'.  In the process,   !
!  the magnitude of the cells' velocity is computed and stored in array       !
!  'u_mag'.                                                                   !
!                                                                             !
!  Francisco Simoes, January 2008                                             !
!  Last updated (mm-dd-yyyy): 04-16-2008 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables.
  INTEGER :: i
  LOGICAL, EXTERNAL :: equals

! Initialize dry cell arrays.
  FORALL (i = 1:n_elems) drycell(i) = h(i) < h_dry
  partdry = .FALSE.

! Compute the magnitude of the cell velocity.  This is done here not only for
! the sake of computational efficiency, but also because u_mag(i) is used to
! determine if a cell is dry or not.
  u_mag = zero
  WHERE (.NOT. drycell) u_mag = SQRT(u*u + v*v)
  !DO i = 1,n_elems
  !  IF (drycell(i) CYCLE
  !  u_mag(i) = SQRT(u(i)*u(i) + v(i)*v(i))
  !END DO

  ! Extra loop to catch the dry cells with stage above h_dry but below h_wet.
  DO i = 1,n_elems
    IF (drycell(i)) CYCLE
    IF ((h(i) < h_wet) .AND. equals(u_mag(i),zero)) drycell(i) = .TRUE.
  END DO

END SUBROUTINE find_dry_cells
