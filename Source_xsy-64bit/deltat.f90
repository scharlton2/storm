FUNCTION deltat(crtype,cfl)
  USE parameters
  USE constants
  USE geometry
  USE dep_vars
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Computes time step size for unsteady computations.                         !
!                                                                             !
!  Francisco Simoes, November 2005                                            !
!  Last updated (mm-dd-yyyy): 02-25-2011 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  INTEGER, INTENT(IN) :: crtype  ! Type of DelT (unused for now).
  REAL (KIND=mp), INTENT(IN) :: cfl  ! Courant number.
  REAL (KIND=mp) :: deltat  ! Time step size.

! Local variables.
  INTEGER :: i
  REAL (KIND=mp) :: delt
  REAL (KIND=mp), PARAMETER :: small = 1.0E-9  ! A very small number.

  deltat = vlarge
  DO i = 1,n_pts
    delt = SQRT(u(i)*u(i) + v(i)*v(i)) + SQRT(g*h(i)) + small
    ! Compute the time step using the CFL criterion from Brufau et al. (2004).
    delt = cfl*cv_area(i)/(delt*cv_perim(i))
    deltat = MIN(delt,deltat)
  END DO

END FUNCTION deltat
