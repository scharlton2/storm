SUBROUTINE positivity_vtx
  USE constants
  USE geometry
  USE dep_vars
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  It was discovered that after some interpolation procedures and boundary    !
!  fixes, some of the water depths at the vertices were negative.  This       !
!  subroutine corrects those values by resetting 'hvtx' to zero and setting   !
!  'zetavtx' to compatible values.                                            !
!                                                                             !
!  Francisco Simoes, April 2008                                               !
!  Last updated (mm-dd-yyyy): 04-21-2008 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables.
  INTEGER :: i,j

  DO i = 1,n_pts
    IF (hvtx(i) < zero) THEN
      hvtx(i) = zero
      zetavtx(i) = zvtx(i)
    END IF
  END DO

END SUBROUTINE positivity_vtx
