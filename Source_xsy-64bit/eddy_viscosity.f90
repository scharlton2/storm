SUBROUTINE eddy_viscosity
  USE parameters
  USE geometry
  USE dep_vars
  USE constants
  USE options
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Parabolic eddy viscosity model.  Computes the eddy viscosity at each       !
!  triangle vertex using the values at the centers of the surrounding         !
!  triangles.  The triangle-to-node interpolation follows the technique of    !
!  Maisano et al. (2006) using the areas of neighboring triangles as          !
!  weights, but only wet triangles are used in the interpolation.             !
!                                                                             !
!  Francisco Simoes, July 2007                                                !
!  Last updated (mm-dd-yyyy): 04-06-2008 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables.
  INTEGER :: i,j,k
  REAL (KIND=mp) :: a,weight,tarea,ustar

  rough = -rough
  DO i = 1,n_pts
    IF (hvtx(i) < h_dry) CYCLE  ! Dry node.
    tarea = zero
    a = zero
    DO j = 2,n2t(i,1)+1
      k = n2t(i,j)
      IF (h(k) > h_dry) THEN
        weight = grid(k)%area
        tarea = tarea + weight
        a = a + weight*ABS(rough(i))*u_mag(i)
      END IF
    END DO
    ustar = SQRT(a/tarea)
    ev(i) = omega*hvtx(i)*ustar
  END DO

END SUBROUTINE eddy_viscosity
