LOGICAL FUNCTION node_in_triang(n,t)
  USE parameters
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This function returns .TRUE. if node n belongs to triangle t, and .FALSE.  !
!  otherwise.                                                                 !
!                                                                             !
!  F. Simoes, November 2004                                                   !
!  Last updated (mm-dd-yyyy): 10-20-2005 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables:
  INTEGER, INTENT(IN) :: n
  TYPE(triangle), INTENT(IN) :: t

  node_in_triang = .FALSE.
  IF (t%vertex(1) == n .OR. t%vertex(2) == n .OR. t%vertex(3) == n) &
    node_in_triang = .TRUE.

END FUNCTION node_in_triang
