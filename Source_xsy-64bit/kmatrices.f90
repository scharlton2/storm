SUBROUTINE kmatrices(kplus,kminus,nxi,nyi,h,u,v)
  USE parameters
  USE constants
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Computes the residual distribution matrices.  Follows the ideas of         !
!  Caraeni (2000).                                                            !
!                                                                             !
!  INPUT:                                                                     !
!    nxi      x-component of the inward-pointing normal to edge i;            !
!    nyi      y-component of the inward-pointing normal to edge i;            !
!    h        element-linearized water surface elevation;                     !
!    u        x-component of the element-linearized velocity;                 !
!    v        y-component of the element-linearized velocity.                 !
!                                                                             !
!  OUTPUT:                                                                    !
!    kplus    the residual distribution matrix associated with the positive   !
!             eigenvalues of the linearized element matrix;                   !
!    kminus   the residual distribution matrix associated with the negative   !
!             eigenvalues of the linearized element matrix.                   !
!                                                                             !
!  Debugging completed 11-03-2004 using MATLAB.                               !
!                                                                             !
!  Francisco Simoes, October 2004                                             !
!  Last updated (mm-dd-yyyy): 11-18-2005 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

  ! Dummy variables:
  REAL(KIND=mp), INTENT(IN) :: nxi,nyi,h,u,v
  REAL(KIND=mp), INTENT(OUT) :: kminus(3,3),kplus(3,3)

  ! Local variables:
  REAL(KIND=mp) :: c,c2,f1,f2,lambda(3),ni,ni2,uni
  REAL(KIND=mp) :: ki(3,3)

  uni = u*nxi + v*nyi
  ni2 = nxi*nxi + nyi*nyi
  ni = SQRT(ni2)
  c2 = g*h
  c = SQRT(c2)

  ! Eigenvalues:
  lambda(1) = uni + c*ni
  lambda(2) = uni
  lambda(3) = uni - c*ni

  !WRITE (*,'(3(2X,ES12.5))') lambda(1),lambda(2),lambda(3)

  ! Linearized element matrix:
  CALL jacobian(ki,nxi,nyi,h,u,v)

  IF (lambda(3) >= zero) THEN  ! All eigenvalues are positive.
    kplus = ki
    kminus = zero

  ELSE IF (lambda(1) <= zero) THEN  ! All eigenvalues negative.
    kminus = ki
    kplus = zero

  ELSE IF (lambda(1) > zero .AND. lambda(2) <= zero) THEN
    ! One eigenvalue is positive, two are negative.

    kplus(1,1) = c - uni/ni; kplus(1,2) = nxi/ni; kplus(1,3) = nyi/ni
    f1 = (u*ni + c*nxi)/ni2
    f2 = (v*ni + c*nyi)/ni2
    kplus(2,1) = f1*(c*ni - uni)
    kplus(2,2) = f1*nxi
    kplus(2,3) = f1*nyi
    kplus(3,1) = f2*(c*ni - uni)
    kplus(3,2) = f2*nxi
    kplus(3,3) = f2*nyi
    kplus = kplus*lambda(1)/c/4.0_mp
    kminus = ki - kplus

  ELSE IF (lambda(3) < zero .AND. lambda(2) >= zero) THEN
    ! One eigenvalue is negative, the other two are positive.

    kminus(1,1) = c + uni/ni; kminus(1,2)= -nxi/ni; kminus(1,3) = -nyi/ni
    f1 = (c*nxi - u*ni)/ni2
    f2 = (c*nyi - v*ni)/ni2
    kminus(2,1) = -f1*(uni + c*ni)
    kminus(2,2) = f1*nxi
    kminus(2,3) = f1*nyi
    kminus(3,1) = -f2*(uni + c*ni)
    kminus(3,2) = f2*nxi
    kminus(3,3) = f2*nyi
    kminus = kminus*lambda(3)/c/4.0_mp
    kplus = ki - kminus

  ELSE
    ! Error: this is just a safety net and the program should never get here.
    PRINT *,""
    PRINT *,"FATAL ERROR: eigenvalue selection failure in Ki decomposition."
    PRINT *,"Please report error to SToRM's development team."
    CALL byebye('Program SToRM stopped.')

  END IF

END SUBROUTINE kmatrices
