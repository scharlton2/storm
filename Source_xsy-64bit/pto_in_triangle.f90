LOGICAL FUNCTION pto_in_triangle(e,x,y)
  USE constants
  USE geometry
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Finds if point (x,y) is inside a triangle e.  The funtion does not need    !
!  that the vertices are given in any particular order and determines if      !
!  the point (x,y) is inside the triangle using barycentric coordinates.      !
!                                                                             !
!  INPUT:                                                                     !
!    e                   triangle;                                            !
!    x, y                point to be checked.                                 !
!                                                                             !
!  OUTPUT:                                                                    !
!    pto_in_triangle     .TRUE. if (x,y) is inside of triangle e or on its    !
!                        boundary;                                            !
!                        .FALSE. otherwise.                                   !
!                                                                             !
!  Francisco Simoes, April 2013                                               !
!  Last updated (mm-dd-yyyy): 04-19-2013 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables:
  TYPE(triangle) :: e  ! Triangle of interest.
  REAL (KIND=mp), INTENT(IN) :: x,y  ! Point to be tested.

! Local variables:
  REAL(KIND=mp) :: s,t,tarea,x1,x2,x3,y1,y2,y3

! Area and vertices of triangle.
  pto_in_triangle = .FALSE.
  tarea = e%area
  x1 = nodes(e%vertex(1))%x
  y1 = nodes(e%vertex(1))%y
  x2 = nodes(e%vertex(2))%x
  y2 = nodes(e%vertex(2))%y
  x3 = nodes(e%vertex(3))%x
  y3 = nodes(e%vertex(3))%y

  s = half/tarea*(y1*x3 - x1*y3 + (y3 - y1)*x + (x1 - x3)*y)
  t = half/tarea*(x1*y2 - y1*x2 + (y1 - y2)*x + (x2 - x1)*y)

! The point (x,y,) is inside triangle e only and if only s, t, and 1-s-t are
! all >= 0.
  IF (s < zero) RETURN
  IF (t < zero) RETURN
  IF (one-s-t < zero) RETURN

  pto_in_triangle = .TRUE.

END FUNCTION pto_in_triangle
