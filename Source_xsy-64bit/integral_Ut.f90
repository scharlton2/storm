FUNCTION integral_Ut(node,phi,ndim)
  USE parameters
  USE geometry
  USE constants
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This function returns the integral of 'phi' over all the triangles that    !
!  contain node 'node'.  The integrand is assumed to vary linearly within     !
!  each triangle.                                                             !
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
  REAL (KIND=mp) :: integral_Ut

! Local variables:
  INTEGER :: i,j,k,l,n
  TYPE(triangle) :: t

  integral_Ut = zero
  n = n2t(node,1)
  DO l = 1,n
    t = grid(n2t(node,l+1))
    i = t%vertex(1)
    j = t%vertex(2)
    k = t%vertex(3)
    integral_Ut = integral_Ut + t%area*(phi(i) + phi(j) + phi(k))*one_third
  END DO

END FUNCTION integral_Ut
