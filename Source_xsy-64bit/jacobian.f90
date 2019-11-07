SUBROUTINE jacobian(ki,nxi,nyi,h,u,v)
  USE parameters
  USE constants
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Computes the Jacobian matrix Ki based on the cell-face interior normal     !
!  vector--eq. (3.35) of Caraeni (2000) modified for 2D--for the adopted      !
!  conservative linearization.                                                !
!                                                                             !
!  INPUT:                                                                     !
!    nxi      x-component of the inward-pointing normal to edge i;            !
!    nyi      y-component of the inward-pointing normal to edge i;            !
!    h        element-linearized water surface elevation;                     !
!    u        x-component of the element-linearized velocity;                 !
!    v        y-component of the element-linearized velocity.                 !
!                                                                             !
!  OUTPUT:                                                                    !
!    ki       the Jacobian matrix.                                            !
!                                                                             !
!  Francisco Simoes, November 2005                                            !
!  Last updated (mm-dd-yyyy): 11-18-2005 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

  ! Dummy variables:
  REAL(KIND=mp), INTENT(IN) :: nxi,nyi,h,u,v
  REAL(KIND=mp), INTENT(OUT) :: ki(3,3)

  ! Local variables:
  REAL(KIND=mp) :: c2,ni,ni2,uni

  uni = u*nxi + v*nyi
  ni2 = nxi*nxi + nyi*nyi
  ni = SQRT(ni2)
  c2 = g*h

  ! Linearized element matrix:
  ki(1,1) = zero;           ki(1,2) = nxi;         ki(1,3) = nyi
  ki(2,1) = c2*nxi - u*uni; ki(2,2) = u*nxi + uni; ki(2,3) = u*nyi
  ki(3,1) = c2*nyi - v*uni; ki(3,2) = v*nxi;       ki(3,3) = v*nyi + uni
  ki = ki*half

END SUBROUTINE jacobian
