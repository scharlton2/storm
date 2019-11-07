SUBROUTINE windforcing2
  USE parameters
  USE constants
  USE geometry
  USE dep_vars
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Compute parameters needed by the wind forcing term computations.           !
!                                                                             !
!  Francisco Simoes, April 2013                                               !
!  Last updated (mm-dd-yyyy): 05-02-2013                                      !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables.
  INTEGER :: i,ierror
  REAL (KIND=mp) :: a
  REAL (KIND=mp), ALLOCATABLE, DIMENSION(:) :: wcoef0  ! Temporary work array.
  REAL (KIND=mp), EXTERNAL :: azimuth2angle

! Allocate variables.
  ALLOCATE(wcoef0(n_elems),wcoef1(n_elems),wcoef2(n_elems),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(8*n_elems*2)
! Move wind friction coefficient (w_fric), direction (w_dir), and magnitude
! (w_mag) to center of triangles.
  CALL vtx2ctr(w_fric,wcoef0)  ! From vertex to center of triangles.
  CALL vtx2ctr(w_dir,wcoef1)
  CALL vtx2ctr(w_mag,wcoef2)
  DO i = 1,n_elems
    a = azimuth2angle(wcoef1(i))  ! Convert from azimuth to angle (degrees).
    wcoef1(i) = wcoef0(i)*wcoef2(i)*COS(a)*wcoef2(i)/rho
    wcoef2(i) = wcoef0(i)*wcoef2(i)*SIN(a)*wcoef2(i)/rho
  END DO
! Free unneeded memory.
  DEALLOCATE(wcoef0,w_fric,w_mag,w_dir)

END SUBROUTINE windforcing2
