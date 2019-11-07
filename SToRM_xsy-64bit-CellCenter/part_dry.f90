LOGICAL FUNCTION part_dry(e,n,i,j)
  USE options
  USE dep_vars
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Determines how many nodes are dry in a partially dry triangle.  It takes   !
!  the triangle 'e' as input and outputs 'n' (the number of dry nodes) and    !
!  the local indeces of those nodes ('i' and 'j').  It returns .TRUE. if one  !
!  or two nodes are dry, and .FALSE. if none or all nodes are dry, in which   !
!  case the BGN algorithm is not applied.  Note that the indeces 'i' and 'j'  !
!  only contain the correct information if n = 1 or 2 upon exit of            !
!  'part_dry', as follows: if n = 1, 'i' contains the index of the dry node   !
!  and 'j' is ignored; if n = 2, 'i' and 'j' contain the indeces of the two   !
!  dry nodes.                                                                 !
!                                                                             !
!  Francisco Simoes, October 2005                                             !
!  Last updated (mm-dd-yyyy): 10-31-2005 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables.
  TYPE(triangle), INTENT(IN) :: e
  INTEGER, INTENT(OUT) :: n,i,j

! Local variables.
  INTEGER :: k

  part_dry = .FALSE.

! Find out how many nodes are dry.
  n = 0
  DO k = 1,3
    IF (h(e%vertex(k)) < h_dry) THEN
      n = n + 1
      j = k
      IF (n == 1) i = k
    END IF
  END DO
  IF (n == 0 .OR. n == 3) RETURN  ! None or all of the nodes are dry.

  part_dry = .TRUE.

END FUNCTION part_dry
