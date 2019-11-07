SUBROUTINE cell_velocity
  USE parameters
  USE geometry
  USE dep_vars
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This subroutine computes the cell velocities by averaging the nodal        !
!  values in each triangle of the mesh.                                       !
!                                                                             !
!  Francisco Simoes, March 2004                                               !
!  Last updated (mm-dd-yyyy): 03-03-2004                                      !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variable:
  INTEGER :: i

  DO i = 1,n_elems
    u_avg(i,1) = (h(grid(i)%vertex(1)) + h(grid(i)%vertex(2)) + &
                  h(grid(i)%vertex(3)))/3.0_mp
    u_avg(i,2) = (u(grid(i)%vertex(1)) + u(grid(i)%vertex(2)) + &
                  u(grid(i)%vertex(3)))/3.0_mp
    u_avg(i,3) = (v(grid(i)%vertex(1)) + v(grid(i)%vertex(2)) + &
                  v(grid(i)%vertex(3)))/3.0_mp
  END DO

END SUBROUTINE cell_velocity
