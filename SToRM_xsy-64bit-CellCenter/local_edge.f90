INTEGER FUNCTION local_edge(element,edg)
  USE parameters
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Given a global edge index, this function searches the element edge         !
!  connectivity table to produce the local table index of that edge.  If the  !
!  element does not contain the edge, the function returns -1.                !
!                                                                             !
!  INPUT:                                                                     !
!    element  triangle of interest;                                           !
!    edg      edge to look for.                                               !
!                                                                             !
!  OUTPUT:                                                                    !
!    local_edge  the local index that triangle 'element' uses to point to     !
!                edge 'edg'.  It returns -1 if 'element' does not contain     !
!                'edg' as one of its sides.                                   !
!                                                                             !
!  Francisco Simoes, March 2007                                               !
!  Last updated (mm-dd-yyyy): 03-22-2007 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables.
  TYPE(triangle), INTENT(IN) :: element
  INTEGER, INTENT(IN) :: edg

! Local variables.
  INTEGER :: i

  DO i = 1,3
    IF (ABS(element%edge(i)) == edg) THEN
      local_edge = i
      RETURN
    END IF
  END DO

  local_edge = -1

END FUNCTION local_edge
