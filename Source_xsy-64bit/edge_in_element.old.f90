INTEGER FUNCTION edge_in_element(e)
  USE geometry
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Determines to which element does edge 'e' belongs to and returns the       !
!  element number.  Note, however, that this function only returns the first  !
!  element that it encounters (in array grid), therefore it is only useful    !
!  for edges on the boundary.                                                 !
!                                                                             !
!  INPUT:                                                                     !
!    e               global edge number, as it is used in array 'edges'.      !
!                                                                             !
!  OUTPUT:                                                                    !
!    edge_in_element 0 if edge 'e' is not found in any triangle in grid();    !
!                    n triangle number in array grid().                       !
!                                                                             !
!  F. Simoes, 2005                                                            !
!  Last updated (mm-dd-yyyy): 11-07-2005 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!
! Dummy variables.
  INTEGER, INTENT(IN) :: e

! Local variables.
  INTEGER :: i,j

  DO i = 1,n_elems
    DO j = 1,3
      IF (e == ABS(grid(i)%edge(j))) THEN
        edge_in_element = i  ! Return the element # where the edge was found.
        RETURN
      END IF
    END DO
  END DO

  edge_in_element = 0  ! Error condition.

END FUNCTION edge_in_element
