SUBROUTINE adjust_water_depth
  USE parameters
  USE dep_vars
  USE geometry
  USE constants
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This subroutine is needed if the bed elevation adjustment algorithm of     !
!  Brufau and Garcia-Navarro (2003) is used.  It corrects the water depth to  !
!  a value compatible with the unchanged bed elevation.  Note that its use    !
!  may result in loss of water volume in the computational domain.            !
!                                                                             !
!  Francisco Simoes, 3 February 2009                                          !
!  Last updated (mm-dd-yyyy): 03-03-2009 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local avriables.
  REAL (KIND=mp) :: delta
  INTEGER :: i

  DO i = 1,n_elems
    delta = zstore(i) - z(i)
    h(i) = h(i) - delta
  END DO
  WHERE (h < zero) h = zero
  z = zstore
  zvtx = zvtxstore

END SUBROUTINE adjust_water_depth
