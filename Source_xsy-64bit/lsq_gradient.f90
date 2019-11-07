SUBROUTINE lsq_gradient(solver,node,phi,ndim,phix,phiy)
  USE parameters
  USE geometry
  USE constants
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Carry out the least-squares gradient construction of variable phi,         !
!  following technique of Mavriplis (2003).                                   !
!                                                                             !
!  INPUT:                                                                     !
!    solver     = 1 for RDS and other vertex-oriented solvers;                !
!               = 2 for FVT and other cellcentered solvers;                   !
!    node       global index of the node (or element, if solver = 2) where    !
!               the gradient is going to be computed;                         !
!    phi        array with the quantity of interest.                          !
!    ndim       dimension of array 'phi', used so that the check bounds       !
!               option of the Fortran compiler can work well in the           !
!               subroutine.                                                   !
!                                                                             !
!  OUTPUT:                                                                    !
!    phix,phiy  x- and y-component of the gradient of phi.                    !
!                                                                             !
!  F. Simoes, January 2007                                                    !
!  Last updated (mm-dd-yyyy): 08-20-2007 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments:
  INTEGER, INTENT(IN) :: ndim,node,solver
  REAL (KIND=mp), INTENT(IN) :: phi(ndim)
  REAL (KIND=mp), INTENT(OUT) :: phix,phiy

! Local variables:
  INTEGER :: i,k,n
  REAL (KIND=mp) :: di,dphi,ei

  SELECT CASE (solver)

!-----------------------------------------------------------------------------!
!                                 RDS solver                                  !
!-----------------------------------------------------------------------------!
  CASE (1)
    di = zero; ei = zero
    n = n2n(node,1)
    i = node  ! This makes it easier to follow eqs (31-35) in Mavriplis (2003).
    DO k = 1,n
      dphi = phi(n2n(i,k+1)) - phi(i)
      di = di + lsqweight(i,k)*dphi*dx(i,k)
      ei = ei + lsqweight(i,k)*dphi*dy(i,k)
    END DO

    ! Solve the system using Cramer's rule.
    phix = (ci(i)*di - bi(i)*ei)/det(i)
    phiy = (ai(i)*ei - di*bi(i))/det(i)

!-----------------------------------------------------------------------------!
!                                 FVT solver                                  !
!-----------------------------------------------------------------------------!
  CASE (2)
    di = zero; ei = zero
    n = t2t(node,1)
    i = node  ! This makes it easier to follow eqs (31-35) in Mavriplis (2003).
    DO k = 1,n
      dphi = phi(t2t(i,k+1)) - phi(i)
      di = di + lsqweight(i,k)*dphi*dx(i,k)
      ei = ei + lsqweight(i,k)*dphi*dy(i,k)
    END DO

    ! Solve the system using Cramer's rule.
    phix = (ci(i)*di - bi(i)*ei)/det(i)
    phiy = (ai(i)*ei - di*bi(i))/det(i)

  CASE DEFAULT
    PRINT *,"ERROR: least-squares gradient computation, invalid option."
    CALL byebye('Program STORM stopped.')

  END SELECT

END SUBROUTINE lsq_gradient
