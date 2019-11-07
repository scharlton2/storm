SUBROUTINE vtx2ctr(varray,carray)
  USE parameters
  USE constants
  USE geometry
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This subroutine computes the values of a desired quantity at the           !
!  geometric centers of the triangles in a STORM grid.                        !
!                                                                             !
!  INPUT:                                                                     !
!    varray    array with the quantity of interest located at the vertices    !
!              of the triangles (dimension = n_pts).                          !
!                                                                             !
!  OUTPUT:                                                                    !
!    carray    array with the same quantity interpolated to the centers of    !
!              the triangles (dimension = n_elems).                           !
!                                                                             !
!  Francisco Simoes, February 2007                                            !
!  Last updated (mm-dd-yyyy): 02-26-2007 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  REAL (KIND=mp), INTENT(IN) :: varray(n_pts)
  REAL (KIND=mp), INTENT(OUT) :: carray(n_elems)

! Local variables.
  INTEGER :: i

  DO i = 1,n_elems
    carray(i) = (varray(grid(i)%vertex(1)) + varray(grid(i)%vertex(2)) + &
      varray(grid(i)%vertex(3)))*one_third
  END DO

END SUBROUTINE vtx2ctr
