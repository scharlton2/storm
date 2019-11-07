SUBROUTINE centered(n,r)
  USE parameters
  USE constants
  USE dep_vars
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Central scheme for residual distribution schemes (beta = 1/3 for all i).   !
!  Use with care.                                                             !
!                                                                             !
!  INPUT:                                                                     !
!    n       element index.                                                   !
!                                                                             !
!  OUTPUT:                                                                    !
!    r(i,j)  residuals for node i of the element n (3 by 3 array).  In the    !
!            calling program, simply add r(i,j) to phi(grid(n)%vertex(i),j)   !
!            to update the nodal residuals.                                   !
!                                                                             !
!  Francisco Simoes, October 2005                                             !
!  Last updated (mm-dd-yyyy): 10-31-2005 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables:
  INTEGER, INTENT(IN) :: n
  REAL(KIND=mp), INTENT(OUT) :: r(3,3)

! Local variables:
  INTEGER :: i,j
  REAL(KIND=mp), PARAMETER :: onethird = one/3.0_mp

  DO i = 1,3
    DO j = 1,3
      r(i,j) = onethird*residual(n,j)
    END DO
  END DO

END SUBROUTINE centered
