LOGICAL FUNCTION pto_in_edge(i,e)
  USE parameters
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This function determines if node i belongs to edge e.                      !
!                                                                             !
!  F. Simoes, October 2005                                                    !
!  Last updated (mm-dd-yyyy): 11-01-2005 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables:
  INTEGER, INTENT(IN) :: i
  TYPE(edge), INTENT(IN) :: e

  pto_in_edge = .TRUE.
  IF (i /= e%p(1) .AND. i /= e%p(2)) pto_in_edge = .FALSE.

END FUNCTION pto_in_edge
