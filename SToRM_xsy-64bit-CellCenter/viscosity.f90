FUNCTION viscosity(t)
  USE parameters
  USE constants
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This function computes the kinematic viscosity of water as a function of   !
!  water temperature, SI units.  For reference, see eq. (1.12) of             !
!  Yang (1996).                                                               !
!                                                                             !
!  F. Simoes, April 2004                                                      !
!  Last updated (mm-dd-yyyy): 04-07-2004                                      !
!                                                                             !
!-----------------------------------------------------------------------------!

! Arguments:
  REAL(KIND=mp), INTENT(IN) :: t  ! Water temperature, in degrees Centigrade.
  REAL(KIND=mp) :: viscosity  ! Water viscosity, in m^2/s.

  viscosity = zero  ! Zero means error.
  IF (t > 0.0_mp .AND. t < 100.0_mp) THEN
    viscosity = REAL(1.792e-6,mp)/(one + 0.0337_mp*t + 0.000221_mp*t*t)
  END IF

END FUNCTION viscosity
