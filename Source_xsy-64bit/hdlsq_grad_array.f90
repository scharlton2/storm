SUBROUTINE hdlsq_grad_array(phi,gradphi,ndim)
  USE parameters
  USE geometry
  USE constants
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Carry out the least-squares gradient construction of variable phi at the,  !
!  center of a triangular control volume, following technique of Mavriplis    !
!  (2003).  Use for cell-centered solvers only.                               !
!                                                                             !
!  INPUT:                                                                     !
!    phi        array with the quantity of interest;                          !
!    ndim       dimension of arrays 'phi' and 'gradphi', used so that the     !
!               check bounds option of the Fortran compiler can work well in  !
!               the subroutine.                                               !
!                                                                             !
!  OUTPUT:                                                                    !
!    gradphi    gradient of phi.                                              !
!                                                                             !
!  F. Simoes, September 2007                                                  !
!  Last updated (mm-dd-yyyy): 10-20-2007 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  INTEGER, INTENT(IN) :: ndim
  REAL (KIND=mp), INTENT(IN) :: phi(ndim)
  TYPE(vector), INTENT(OUT) :: gradphi(ndim)

! Local variables.
  INTEGER :: i,k,n
  REAL (KIND=mp) :: di,dphi,ei

  DO i = 1,n_elems
    di = zero; ei = zero
    n = t2tHD(i,1)
    DO k = 1,n
      dphi = phi(t2tHD(i,k+1)) - phi(i)
      di = di + lsqweight(i,k)*dphi*dx(i,k)
      ei = ei + lsqweight(i,k)*dphi*dy(i,k)
    END DO

    ! Solve the system using Cramer's rule.
    gradphi(i)%x = (ci(i)*di - bi(i)*ei)/det(i)
    gradphi(i)%y = (ai(i)*ei - di*bi(i))/det(i)
  END DO

END SUBROUTINE hdlsq_grad_array
