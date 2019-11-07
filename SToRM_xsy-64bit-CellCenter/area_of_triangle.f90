FUNCTION area_of_triangle(x1,y1,x2,y2,x3,y3)
  USE parameters
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Computes the area of a triangle.                                           !
!                                                                             !
!  INPUT:                                                                     !
!    x1,y1,x2,y2,x3,y3  coordinates of the vertices of the triangle.          !
!                                                                             !
!  OUTPUT:                                                                    !
!    area_of_triangle   Area of the triangle, > 0 if the vertices are         !
!                       oriented in a counterclockwise manner, < 0 if in a    !
!                       clockwise manner.                                     !
!                                                                             !
!  Francisco Simoes, February 2009                                            !
!  Last updated (mm-dd-yyyy): 02-03-2009 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  REAL(KIND=mp), INTENT(IN) :: x1,x2,x3,y1,y2,y3
  REAL(KIND=mp) :: area_of_triangle

  area_of_triangle = ((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1))/2.0_mp

END FUNCTION triangle
