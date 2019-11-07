FUNCTION azimuth2angle(x)
  USE parameters  ! Defines KIND = mp.
  USE constants
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Converts the azimuth to a conventional mathematical angle.                 !
!                                                                             !
!  INPUT:                                                                     !
!    x               Azimuth (angle, in degrees, measured from the Northing   !
!                    in the clockwise direction).  Must be in [0,360].        !
!                                                                             !
!  OUTPUT:                                                                    !
!    azimuth2angle   Conventional mathematical angle, in radians.             !
!                                                                             !
!  Francisco Simoes, May 2013                                                 !
!  Last updated (mm-dd-yyyy): 05-02-2013                                      !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  REAL (KIND=mp), INTENT(IN) :: x
  REAL (KIND=mp) :: azimuth2angle
! Local variables.
  REAL (KIND=mp) :: a

  IF (x < zero .OR. x > 360_mp) THEN
    PRINT *,"Wrong argument in wind direction: value < 0 or > 360 degrees."
    CALL byebye("Program SToRM stopped.")
  END IF

  a = (360_mp - x) + 90_mp
  IF (a < zero) THEN
    a = a + 360_mp
  ELSE IF (a >= 360_mp) THEN
    a = a - 360_mp
  END IF

  azimuth2angle = a*pi/180.0_mp

END FUNCTION azimuth2angle
