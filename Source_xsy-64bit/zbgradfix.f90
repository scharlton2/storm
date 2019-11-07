SUBROUTINE zbgradfix
  USE geometry
  USE constants
  USE dep_vars
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!                                                                             !
!  Francisco Simoes, November 2007                                            !
!  Last updated (mm-dd-yyyy): 11-03-2007 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables.
  INTEGER :: i,k,l
  INTEGER, EXTERNAL :: edge_in_element

  DO l = 1,wall_edges1
    k = walledg1(l)
    i = edge_in_element(edges(k))
    delz(i)%y = zero
  END DO

END SUBROUTINE zbgradfix
