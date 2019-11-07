SUBROUTINE limiterBJ(beta,q,gradq,ndim)
  USE parameters
  USE constants
  USE geometry
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Limits the gradient using the Barth and Jespersen (1989) technique.  The   !
!  implementation follows closely the equations in section 4.1 of Anastasiou  !
!  and Chan (1997).                                                           !
!                                                                             !
!  INPUT:                                                                     !
!    beta    a parameter used in the limiter (1 <= beta <= 2). If beta = 1    !
!            this is the minmod limiter; if beta = 2 it is the Superbee;      !
!    q       the scalar quantity of interest (it can be u, v, z, etc.);       !
!    gradq   the gradient of q;                                               !
!    ndim    dimension of arrays 'q' and 'gradq' are used so that the check   !
!            bounds option of the Fortran compiler can work well in the       !
!            subroutine.                                                      !
!                                                                             !
!  INPUT:                                                                     !
!    gradq   the gradient of q, modified by the appropriate use of the        !
!            limiter.                                                         !
!                                                                             !
!  F. Simoes, March 2007                                                      !
!  Last updated (mm-dd-yyyy): 04-08-2008 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  INTEGER, INTENT(IN) :: ndim
  REAL (KIND=mp), INTENT(IN) :: beta,q(ndim)
  TYPE(vector), INTENT(INOUT) :: gradq(ndim)

! Local variables.
  INTEGER :: i,j,k
  REAL (KIND=mp) :: a,b,delq,phi,q0,qj,qmax,qmin,qneighb,r(3)
  !INTEGER, EXTERNAL :: neighbor
  LOGICAL, EXTERNAL :: equals

  DO i = 1,n_elems  ! Sweep over all elements.
    q0 = q(i)
    DO j = 1,3  ! Sweep over element faces.
      ! Get the neighbor of i sharing face j:
      !k = neighbor(i,ABS(grid(i)%edge(j)))
      k = t2t3(j,i)
      !IF (k < 1) THEN  ! Boundary edge.
      !  r(j) = one
      !  CYCLE
      !END IF
      IF (k < 1) THEN  ! Boundary node.
        qneighb = q0
      ELSE
        qneighb = q(k)
      END IF
      ! Eqs. (19) of Anastasiou and Chan (1997).
      qmin = MIN(q0,qneighb)
      qmax = MAX(q0,qneighb)
      ! Unlimited gradient reconstruction.
      qj = q0 + gradq(i)%x*rc(i,j)%x + gradq(i)%y*rc(i,j)%y
      delq = qj - q0
      ! Eqs. (18) of Anastasiou and Chan (1997).
      IF (equals(delq,zero)) THEN
        r(j) = one
      ELSE IF (delq < 0) THEN
        r(j) = (qmin - q0)/delq
      ELSE
        r(j) = (qmax - q0)/delq
      END IF
      ! Eq. (17) of Anastasiou and Chan (1997).
      a = MIN(beta*r(j),one)
      b = MIN(r(j),beta)
      r(j) = MAX(a,b)
    END DO
    ! Eq. (16) of Anastasiou and Chan (1997).
    phi = MIN(r(1),r(2),r(3))
    ! Gradient modified for eq. (15) of Anastasiou and Chan (1997).
    gradq(i)%x = phi*gradq(i)%x
    gradq(i)%y = phi*gradq(i)%y
  END DO

END SUBROUTINE limiterBJ
