INTEGER FUNCTION neighbor(elemnt,edg)
  USE geometry
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Finds the triangle that is neighbor of element "elemnt" and that shares    !
!  edge "edg" with it.  Note that typical calls to this function are          !
!                                                                             !
!    n = neighbor(i,ABS(grid(i)%edge(j)))  ! j = 1,2,3                        !
!                                                                             !
!  or                                                                         !
!                                                                             !
!    n = neighbor(edges(i)%e(j),i)  ! j = 1,2                                 !
!                                                                             !
!  Note that neighbor(edges(i)%e(1),i) = edges(i)%e(2) and vice-versa.  Note  !
!  also the presence of ABS() in the first call, because in SToRM's database  !
!  grid()%edge() < 0 if the normal points into the element.  Note also that   !
!  this function does not do any error checking, which means that 'edg' must  !
!  really be an edge of 'elemnt' for it to work properly.                     !
!                                                                             !
!  F. Simoes, March 2007                                                      !
!  Last updated (mm-dd-yyyy): 08-23-2007 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  INTEGER, INTENT(IN) :: edg,elemnt

  neighbor = edges(edg)%e(1)
  IF (neighbor == elemnt) neighbor = edges(edg)%e(2)

END FUNCTION neighbor
