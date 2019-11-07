FUNCTION Shields_tau(d)
  USE parameters
  USE constants
  USE options
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Returns the value of the critical bed shear stress, in the Shields         !
!  diagram, corresponding to diameter "d".  The curve used is the modified    !
!  Soulsby curve.
!                                                                             !
!  INPUT:                                                                     !
!    d              Value of the particle diameter of interest (m).           !
!                                                                             !
!  OUTPUT:                                                                    !
!    Shields_tau    Critical bed shear stress for particles of diameter d     !
!                   (n/m^2).                                                  !
!                                                                             !
!  Francisco Simoes, July 2014                                                !
!  Last updated (mm-dd-yyyy): 07-07-2014 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  REAL (KIND=mp), INTENT(IN) :: d
  REAL (KIND=mp) :: Shields_tau

! Local variables.
  REAL (KIND=mp) :: dstar,s,theta

  s = rhos/rho
  dstar = (g*(s - one)/visc**2)**(one_third)*d
! Modified Soulsby curve using nondimensional parameters.  This is where new
! curves can be programmed in, in the D* - Theta space.
  theta = 0.30_mp/(one + 1.2_mp*dstar) + 0.055_mp*(one - EXP(-0.020_mp*dstar))
! Put dimensions in the final result (N/m^2).
  Shields_tau = g*(rhos-rho)*d*theta

END FUNCTION Shields_tau
