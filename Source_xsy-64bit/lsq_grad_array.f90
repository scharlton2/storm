SUBROUTINE lsq_grad_array(solver,phi,gradphi,ndim)
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
!    solver = 1 for RDS and other vertex-oriented solvers;                    !
!             2 for FVT and other cellcentered solvers;                       !
!    phi        array with the quantity of interest;                          !
!    ndim       dimension of arrays 'phi' and 'gradphi', used so that the     !
!               check bounds option of the Fortran compiler can work well in  !
!               the subroutine.                                               !
!                                                                             !
!  OUTPUT:                                                                    !
!    gradphi    gradient of phi.                                              !
!                                                                             !
!  Francisco Simoes, March 2007                                               !
!  Last updated (mm-dd-yyyy): 04-14-2008 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  INTEGER, INTENT(IN) :: ndim,solver
  REAL (KIND=mp), INTENT(IN) :: phi(ndim)
  TYPE(vector), INTENT(OUT) :: gradphi(ndim)

! Local variables.
  INTEGER :: i,k,n
  REAL (KIND=mp) :: di,dphi,ei

  SELECT CASE (solver)

!-----------------------------------------------------------------------------!
!                                 RDS solver                                  !
!-----------------------------------------------------------------------------!
  CASE (1)
    DO i = 1,n_pts
      di = zero; ei = zero
      n = n2n(i,1)
      DO k = 1,n
        dphi = phi(n2n(i,k+1)) - phi(i)
        di = di + lsqweight(i,k)*dphi*dx(i,k)
        ei = ei + lsqweight(i,k)*dphi*dy(i,k)
      END DO

      ! Solve the system using Cramer's rule.
      gradphi(i)%x = (ci(i)*di - bi(i)*ei)/det(i)
      gradphi(i)%y = (ai(i)*ei - di*bi(i))/det(i)
    END DO

!-----------------------------------------------------------------------------!
!                                 FVT solver                                  !
!-----------------------------------------------------------------------------!
  CASE (2)
    DO i = 1,n_elems
      di = zero; ei = zero
      n = t2t(i,1)
      DO k = 1,n
        dphi = phi(t2t(i,k+1)) - phi(i)
        di = di + lsqweight(i,k)*dphi*dx(i,k)
        ei = ei + lsqweight(i,k)*dphi*dy(i,k)
      END DO

      ! Solve the system using Cramer's rule.
      gradphi(i)%x = (ci(i)*di - bi(i)*ei)/det(i)
      gradphi(i)%y = (ai(i)*ei - di*bi(i))/det(i)
    END DO

  CASE DEFAULT
    PRINT *,"ERROR: least-squares gradient computation, invalid option."
    CALL byebye('Program STORM stopped.')

  END SELECT

END SUBROUTINE lsq_grad_array
