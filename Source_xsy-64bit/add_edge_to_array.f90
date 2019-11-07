SUBROUTINE add_edge_to_array(edg,edges_array,n)
  USE parameters
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This subroutine adds edge edg to array edges_array and increments n.  n    !
!  is the number of elements in edges_array (n >= 0).                         !
!                                                                             !
!  Francisco Simoes, March 2004                                               !
!  Last updated (mm-dd-yyyy): 03-11-2004                                      !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments:
  INTEGER, INTENT(INOUT) :: n
  TYPE(edge), INTENT(IN) :: edg
  TYPE(edge), INTENT(INOUT) :: edges_array(n+1)

  n = n + 1
  edges_array(n) = edg

END SUBROUTINE add_edge_to_array
