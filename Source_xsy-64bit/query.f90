INTEGER FUNCTION query(q)
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  A simple function that returns 1 if q is .FALSE., and 2 if q is .TRUE.,    !
!  to simplify subroutine 'header'.                                           !
!                                                                             !
!  F. Simoes, November 2005                                                   !
!  Last updated (mm-dd-yyyy): 11-19-2005 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variable.
  LOGICAL :: q

  query = 1
  IF (q) query = 2

END FUNCTION query
