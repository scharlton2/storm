INTEGER FUNCTION is_edge(p1,p2)
  USE geometry
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This function finds if the two points p1 and p2 define an edge in the      !
!  geometry database built by SToRM and, if they do, returns the index of     !
!  that edge.                                                                 !
!                                                                             !
!  INPUT:                                                                     !
!    p1,p2      global indexes of two points.                                 !
!                                                                             !
!  OUTPUT:                                                                    !
!    is_edge   -1 if any of the input indexes are out of range;               !
!               0 if the two points do not define an edge;                    !
!               n if the two points define the edge with index n in the       !
!                 SToRM geometry data base.                                   !
!                                                                             !
!  F. Simoes, November 2004                                                   !
!  Last updated (mm-dd-yyyy): 11-07-2005 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

  ! Dummy arguments:
  INTEGER, INTENT(IN) :: p1,p2

  ! Local variables:
  INTEGER :: e,e1,e2

! Check for validity of the dummy argument variables.
  IF (p1 > n_pts .OR. p1 < 1 .OR. p2 > n_pts .OR. p2 < 1) THEN
    is_edge = -1  ! Return an error condition (points out of range).
    RETURN
  END IF

! Sweep all the edges to find the right edge.
  DO e = 1,n_edges
    e1 = edges(e)%p(1)
    e2 = edges(e)%p(2)
    IF ((e1 == p1 .AND. e2 == p2) .OR. (e1 == p2 .AND. e2 == p1)) THEN
      is_edge = e  ! The desired result was found.
      RETURN
    END IF
  END DO

  is_edge = 0  ! Error, no edge found.

END FUNCTION is_edge
