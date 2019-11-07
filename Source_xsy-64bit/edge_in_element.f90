INTEGER FUNCTION edge_in_element(e)
  USE parameters
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Returns an element containing edge 'e'.  Note that this function only      !
!  returns the first element in the data structure, therefore it is only      !
!  useful for edges on the boundary.                                          !
!                                                                             !
!  INPUT:                                                                     !
!    e                 edge.                                                  !
!                                                                             !
!  OUTPUT:                                                                    !
!    edge_in_element   triangle number in array grid().                       !
!                                                                             !
!  F. Simoes, 2005                                                            !
!  Last updated (mm-dd-yyyy): 03-20-2007 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables.
  TYPE(edge), INTENT(IN) :: e

  edge_in_element = e%e(1)
  IF (edge_in_element < 1) edge_in_element = e%e(2)

END FUNCTION edge_in_element
