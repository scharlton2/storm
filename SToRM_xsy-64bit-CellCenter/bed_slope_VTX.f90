SUBROUTINE bed_slope_VTX
  USE parameters
  USE constants
  USE geometry
  USE dep_vars
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Compute the bed slope using vertex z values within a triangle and the      !
!  traditional finite element interpolation.                                  !
!                                                                             !
!  Francisco Simoes, April 2008                                               !
!  Last updated (mm-dd-yyyy): 04-21-2008 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables.
  INTEGER :: i
  REAL (KIND=mp) x1,x2,x3,y1,y2,y3,zz1,zz2,zz3

  DO i = 1,n_elems

    x1 = nodes(grid(i)%vertex(1))%x; y1 = nodes(grid(i)%vertex(1))%y
    x2 = nodes(grid(i)%vertex(2))%x; y2 = nodes(grid(i)%vertex(2))%y
    x3 = nodes(grid(i)%vertex(3))%x; y3 = nodes(grid(i)%vertex(3))%y

    zz1 = zvtx(grid(i)%vertex(1))
    zz2 = zvtx(grid(i)%vertex(2))
    zz3 = zvtx(grid(i)%vertex(3))

    delz(i)%x = half*(zz1*(y2 - y3) + zz2*(y3 - y1) + zz3*(y1 - y2))/grid(i)%area
    delz(i)%y = half*(zz1*(x3 - x2) + zz2*(x1 - x3) + zz3*(x2 - x1))/grid(i)%area

    !delz(i)%x = zero
    !delz(i)%y = zero

  END DO

END SUBROUTINE bed_slope_VTX
