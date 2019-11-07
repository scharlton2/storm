SUBROUTINE fvsources
  USE parameters
  USE geometry
  USE dep_vars
  USE constants
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This subroutine computes the source/sink terms using the finite volume     !
!  formulation on the median dual-mesh.                                       !
!                                                                             !
!  Francisco Simoes, January 2007                                             !
!  Last updated (mm-dd-yyyy): 08-20-2007 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables:
  INTEGER :: i
  REAL (KIND=mp) :: a,b,c,zx,zy
  REAL (KIND=mp), EXTERNAL :: integral_md
  REAL (KIND=mp), EXTERNAL :: integral_Ut

! For mobil bed computations, the bed slope should be computed here for each
! time step. For fixed bed, it can be computed once-and-for-all at the
! beginning of the calculations.
  DO i = 1,n_pts
    CALL lsq_gradient(1,i,z,n_pts,zx,zy)
    zbx(i) = -zx*h(i)
    zby(i) = -zy*h(i)
  END DO

! Compute the source/sink term due to bed slope using the simplest type of
! approximation: constant within the cell.
!  DO i = 1,n_pts
!    source(i,1) = zero
!    source(i,2) = g*zbx(i)*cv_area(i)
!    source(i,3) = g*zby(i)*cv_area(i)
!  END DO

! Compute the source/sink term due to bed slope on the median-dual mesh.
  DO i = 1,n_pts
    source(i,1) = zero
    source(i,2) = g*integral_md(i,zbx,n_pts)
    source(i,3) = g*integral_md(i,zby,n_pts)
  END DO

! Compute the source/sink term due to bed slope on the set of triangles that
! contain point i.
!  DO i = 1,n_pts
!    source(i,1) = zero
!    source(i,2) = g*integral_Ut(i,zbx,n_pts)
!    source(i,3) = g*integral_Ut(i,zby,n_pts)
!  END DO

!  DO i = 1,n_pts
!    a = g*zbx(i)*cv_area(i)
!    b = g*integral_md(i,zbx,n_pts)
!    c = g*integral_Ut(i,zbx,n_pts)
!    WRITE (*,'(I5,A,3ES24.15)') i,'  x',a,b,c
!    a = g*zby(i)*cv_area(i)
!    b = g*integral_md(i,zby,n_pts)
!    c = g*integral_Ut(i,zby,n_pts)
!    WRITE (*,'(I5,A,3ES24.15)') i,'  y',a,b,c
!  END DO
!  STOP

END SUBROUTINE fvsources
