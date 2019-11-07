SUBROUTINE gaussgrad(phi,gradphi,ndim)
  USE parameters
  USE geometry
  USE constants
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Computation of the gradients using the technique described in Anastasiou   !
!  and Chan (1997) for the inviscid fluxes, eqs. (11-13).  Don't forget to    !
!  CALL subroutine 'gaussgrad_setup' to set-up the coefficients used in       !
!  these calculations.                                                        !
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
!  Francisco Simoes, February 2009                                            !
!  Last updated (mm-dd-yyyy): 02-18-2009 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  INTEGER, INTENT(IN) :: ndim
  REAL (KIND=mp), INTENT(IN) :: phi(ndim)
  TYPE(vector), INTENT(OUT) :: gradphi(ndim)

! Local variables.
  INTEGER :: i

! Compute the primary cell gradient of quantity Phi.
  DO i = 1,ndim
    gwgrad(i)%x = gcoefx(i,1)*phi(gconn(i,1)) + gcoefx(i,2)*phi(gconn(i,2)) + &
                  gcoefx(i,3)*phi(gconn(i,3))
    gwgrad(i)%y = gcoefy(i,1)*phi(gconn(i,1)) + gcoefy(i,2)*phi(gconn(i,2)) + &
                  gcoefy(i,3)*phi(gconn(i,3))
  END DO

! Compute the final, weighted gradient of quantity Phi.
  DO i = 1,ndim
    gradphi(i)%x = half*(gwgrad(i)%x + gweight(i,1)*gwgrad(gconn(i,1))%x + &
                   gweight(i,2)*gwgrad(gconn(i,2))%x + &
                   gweight(i,3)*gwgrad(gconn(i,3))%x)
    gradphi(i)%y = half*(gwgrad(i)%y + gweight(i,1)*gwgrad(gconn(i,1))%y + &
                   gweight(i,2)*gwgrad(gconn(i,2))%y + &
                   gweight(i,3)*gwgrad(gconn(i,3))%y)
  END DO

END SUBROUTINE gaussgrad
