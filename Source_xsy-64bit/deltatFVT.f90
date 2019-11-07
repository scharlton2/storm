FUNCTION deltatFVT(crtype,cfl)
  USE parameters
  USE constants
  USE geometry
  USE dep_vars
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Computes time step size for unsteady computations.                         !
!                                                                             !
!  NOTE: this function is used only in the FVT solver, because it is the      !
!  only solver with the arrays 'u_mag' and 'wavec' properly implemented.      !
!  Eventually, when the other solvers are updated they will also be able to   !
!  use this fuction.                                                          !
!                                                                             !
!  Francisco Simoes, November 2005                                            !
!  Last updated (mm-dd-yyyy): 12-11-2007 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  INTEGER, INTENT(IN) :: crtype  ! Type of DelT (unused for now).
  REAL(KIND=mp), INTENT(IN) :: cfl  ! Courant number.
  REAL(KIND=mp) :: deltatFVT  ! Time step size.

! Local variables.
  INTEGER :: i
  REAL(KIND=mp) :: delt
  REAL (KIND=mp), PARAMETER :: small = 1.0E-9  ! A very small number.

  deltatFVT = vlarge
  DO i = 1,n_elems
    delt = u_mag(i) + wavec(i) + small
    ! Compute the time step using the CFL criterion from Brufau et al. (2004).
    delt = cfl*cv_area(i)/(delt*cv_perim(i))
    deltatFVT = MIN(delt,deltatFVT)
  END DO

END FUNCTION deltatFVT
