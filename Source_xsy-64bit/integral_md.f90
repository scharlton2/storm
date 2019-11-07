FUNCTION integral_md(node,phi,ndim)
  USE parameters
  USE geometry
  USE constants
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This function returns the integral of 'phi' over the median-dual cell of   !
!  node 'node'.  The integrand is assumed to vary linearly within each        !
!  triangle.                                                                  !
!                                                                             !
!  INPUT:                                                                     !
!    node  node of interest (center of control volume);                       !
!    phi   quantity to integrate;                                             !
!    ndim  dimension of array 'phi', used so that the check bounds option of  !
!          the Fortran compiler can work well in the subroutine.              !
!                                                                             !
!  Francisco Simoes, January 2007                                             !
!  Last updated (mm-dd-yyyy): 08-20-2007 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables:
  INTEGER, INTENT(IN) :: ndim,node
  REAL (KIND=mp), INTENT(IN) :: phi(ndim)
  REAL (KIND=mp) :: integral_md
  REAL (KIND=mp), PARAMETER :: c1 = 2.03703703703704, &   ! 11/54
                               c2 = 6.481481481481481E-2  ! 7/108

! Local variables:
  INTEGER :: j,k,l,m,n
  TYPE(triangle) :: t

  integral_md = zero
  n = n2t(node,1)
  DO l = 1,n
    t = grid(n2t(node,l+1))
    DO m = 1,3  ! Find the node of the triangle that matches 'node'.
      IF (t%vertex(m) == node) THEN
        j = MAX(1,MOD(m+1,4))
        k = MAX(1,MOD(j+1,4))
        EXIT
      END IF
    END DO
    j = t%vertex(j)
    k = t%vertex(k)
    integral_md = integral_md + t%area*(c1*phi(node) + c2*(phi(j) + phi(k)))
  END DO

END FUNCTION integral_md
