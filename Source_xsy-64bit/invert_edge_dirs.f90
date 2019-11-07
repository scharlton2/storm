FUNCTION invert_edge_dirs(e)
  USE parameters
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  A function to reverse the direction of the normal and tangent vectors      !
!  associated with edge 'e'.  The function always returns -1 to be useful     !
!  to subroutine 'invic_fluxes'.                                              !
!                                                                             !
!  Francisco Simoes, August 2007                                              !
!  Last updated (mm-dd-yyyy): 08-05-2007 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables:
  REAL (KIND=mp) :: invert_edge_dirs
  TYPE(edge), INTENT(INOUT) :: e

  e%normal(1) = -e%normal(1)
  e%normal(2) = -e%normal(2)
  e%tang(1) = -e%tang(1)
  e%tang(2) = -e%tang(2)

  invert_edge_dirs = -1.0_mp

END FUNCTION invert_edge_dirs
