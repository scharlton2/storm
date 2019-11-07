SUBROUTINE edge_normals(side,nside,pto,npto)
  USE parameters
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This subroutine computes the normals to all the edges.  The magnitude of   !
!  the normals is the edge length.  The edges are assumed to be directional   !
!  and point from point 1 to point 2 (Blazek, 2001, pg 79).  If the edge      !
!  appears in the element in a counterclockwise manner, the normal points in  !
!  the outward direction, otherwise it points into the element.  The tangent  !
!  vectors are also computed here, pointing from node 1 to node 2 in a        !
!  counterclockwise manner from the normals.                                  !
!                                                                             !
!  Francisco Simoes, March 2004                                               !
!  Last updated (mm-dd-yyyy): 08-05-2004                                      !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments:
  INTEGER, INTENT(IN) :: npto   ! Number of points in array pto.
  INTEGER, INTENT(IN) :: nside  ! Number of edges in array side.
  TYPE(edge), DIMENSION(nside) :: side  ! Array containing the triangle edges.
  TYPE(point), DIMENSION(npto) :: pto  ! Array containing vertex coordinates.

! Local variables:
  REAL(KIND=mp) :: x1,x2,y1,y2
  INTEGER :: i

  DO i = 1,nside
    ! Edge end points:
    x1 = pto(side(i)%p(1))%x ; y1 = pto(side(i)%p(1))%y
    x2 = pto(side(i)%p(2))%x ; y2 = pto(side(i)%p(2))%y
    ! Compute normal:
    side(i)%normal(1) = y2 - y1  ! x-component of the normal.
    side(i)%normal(2) = x1 - x2  ! y-    "     "   "    "
    side(i)%tang(1) = x2 - x1  ! x-component of the tangent vector.
    side(i)%tang(2) = y2 - y1  ! y-    "     "   "     "      "
  END DO

END SUBROUTINE edge_normals
