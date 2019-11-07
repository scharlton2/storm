SUBROUTINE ctr2vtxg(varc,gradvarc,varv)
  USE parameters
  USE constants
  USE geometry
  USE dep_vars
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Interpolation from the centers to the vertices of triangles.  This         !
!  subroutine uses the gradient of the variable of interest defined over      !
!  each triangle, and only uses wet triangles.                                !
!                                                                             !
!  INPUT:                                                                     !
!    varc      variable of interest located at the center of the triangles;   !
!    gradvarc  gradient of varc defined over each triangle.                   !
!                                                                             !
!  OUTPUT:                                                                    !
!    varv      values of the variable of interest at the vertices of the      !
!              computational grid.                                            !
!                                                                             !
!  Francisco Simoes, January 2008                                             !
!  Last updated (mm-dd-yyyy): 01-14-2008 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  REAL (KIND=mp), INTENT(IN) :: varc(n_elems)
  TYPE(vector), INTENT(IN) :: gradvarc(n_elems)
  REAL (KIND=mp), INTENT(OUT) :: varv(n_pts)

! Local variables.
  INTEGER :: counter,i,j,k
  REAL (KIND=mp) :: newvar

  DO i = 1,n_pts
    counter = 0  ! Counts the number of wet nodes.
    newvar = zero
    DO j = 2,n2t(i,1) + 1
      k = n2t(i,j)
      IF (h(k) < h_dry) CYCLE  ! Dry triangle.
      counter = counter + 1
      newvar = newvar + varc(k) + (nodes(i)%x - grid(k)%xc)*gradvarc(k)%x + &
                                  (nodes(i)%y - grid(k)%yc)*gradvarc(k)%y
    END DO
    varv(i) = newvar/MAX(one,REAL(counter,mp))  ! Avoids division by zero.
  END DO

END SUBROUTINE ctr2vtxg
