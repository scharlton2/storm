SUBROUTINE hequivalent(e,hleft,hright)
  USE parameters
  USE constants
  USE geometry
  USE dep_vars
  USE options
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Computes the equivalent water depth at an edge using eq. (29) of Kuiry et  !
!  al. (2008) and resets the left and right water depth to this value.        !
!                                                                             !
!  INPUT:                                                                     !
!    e               edge where the procedure is applied.                     !
!                                                                             !
!  INPUT:                                                                     !
!    hleft,hright    contain the value of the equivalent water depth for the  !
!                    specified edge.  They are identical and                  !
!                    interchangeable.                                         !
!                                                                             !
!  NOTE: this is an experimental feature.  It seems to be unecessary for      !
!  second order approximation of the water surface profiles, according to     !
!  pg. 235 of Kuiry et al. (2008).  I haven't thought this through, yet...    !
!                                                                             !
!  Francisco Simoes, April 2008                                               !
!  Last updated (mm-dd-yyyy): 02-04-2009 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  REAL (KIND=mp), INTENT(INOUT) :: hleft,hright
  TYPE(edge), INTENT(IN) :: e

! Local variables.
  REAL (KIND=mp) :: h1,h2


  ! Compute the equivalent water depth at the element's edge.
  IF (ABS(hleft - hright)/MAX(hleft,hright) < delta_hequiv) THEN
    h1 = hvtx(e%p(1))
    h2 = hvtx(e%p(2))
    hright = SQRT((h1*h1 + h1*h2 + h2*h2)*one_third)
    hleft = hright
  END IF

END SUBROUTINE hequivalent
