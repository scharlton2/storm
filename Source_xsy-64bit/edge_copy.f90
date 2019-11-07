SUBROUTINE edge_copy(in_edge,out_edge)
  USE parameters
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This subroutines copies all the contents of edge 'in_edge' into the        !
!  variable 'out_edge'.  The data TYPE of bothe variables is 'edge'.          !
!                                                                             !
!  INPUT:                                                                     !
!    in_edge    'edge' data TYPE variable with the known contents.            !
!                                                                             !
!  OUTPUT:                                                                    !
!    out_edge   copy of 'in_edge'.                                            !
!                                                                             !
!  Francisco Simoes, August 2007                                              !
!  Last updated (mm-dd-yyyy): 08-05-2007 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables.
  TYPE(edge), INTENT(IN) :: in_edge
  TYPE(edge), INTENT(OUT) :: out_edge

  out_edge%p(1) = in_edge%p(1)
  out_edge%p(2) = in_edge%p(2)
  out_edge%e(1) = in_edge%e(1)
  out_edge%e(2) = in_edge%e(2)
  out_edge%normal(1) = in_edge%normal(1)
  out_edge%normal(2) = in_edge%normal(2)
  out_edge%tang(1) = in_edge%tang(1)
  out_edge%tang(2) = in_edge%tang(2)
  out_edge%length = in_edge%length

END SUBROUTINE edge_copy
