SUBROUTINE t_areas(grid,ngrid,pto,npto,edg,nedg)
  USE parameters
  USE constants
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This subroutine computes the areas of all the triangular elements in the   !
!  mesh system.  It uses Gauss' formula to compute each triangle area.  The   !
!  order in which the vertices are defined (i.e., clockwise or counter-       !
!  clockwise) is not important. (Blazek, 2001, pg 134)  It also computes      !
!  other important geometric quantities, such as the coordinates of the       !
!  center of the triangle and its perimeter.                                  !
!                                                                             !
!  Francisco Simoes, March 2004                                               !
!  Last updated (mm-dd-yyyy): 08-14-2007                                      !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments:
  INTEGER, INTENT(IN) :: ngrid  ! Number of triangles in array grid.
  INTEGER, INTENT(IN) :: npto   ! Number of points in array pto.
  INTEGER, INTENT(IN) :: nedg  ! Number of edges in array edge.
  TYPE(triangle), DIMENSION(ngrid) :: grid  ! Array containing the triangles.
  TYPE(point), DIMENSION(npto) :: pto  ! Array containing vertex coordinates.
  TYPE(edge), DIMENSION(nedg)  :: edg  ! Array containing the edge data.

! Local variables:
  TYPE(triangle) :: elem
  REAL(KIND=mp) :: area,x1,x2,x3,y1,y2,y3
  INTEGER :: i

  DO i = 1,ngrid
    elem = grid(i)  ! elem is used for speed.
    ! Triangle vertices:
    x1 = pto(elem%vertex(1))%x ; y1 = pto(elem%vertex(1))%y
    x2 = pto(elem%vertex(2))%x ; y2 = pto(elem%vertex(2))%y
    x3 = pto(elem%vertex(3))%x ; y3 = pto(elem%vertex(3))%y
    ! Compute area:
    area = half*ABS((x1 - x2)*(y1 + y2) + (x2 - x3)*(y2 + y3) + &
           (x3 - x1)*(y3 + y1))
    grid(i)%area = area

    ! Compute also the triangle's center.
    grid(i)%xc = one_third*(x1 + x2 + x3)
    grid(i)%yc = one_third*(y1 + y2 + y3)

    ! Compute also the triangle's perimeter.
    grid(i)%perim = edg(ABS(elem%edge(1)))%length + &
                    edg(ABS(elem%edge(2)))%length + &
                    edg(ABS(elem%edge(3)))%length
  END DO

END SUBROUTINE t_areas
