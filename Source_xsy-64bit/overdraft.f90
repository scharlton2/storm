SUBROUTINE overdraft
  USE parameters
  USE constants
  USE geometry
  USE dep_vars
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This subroutine computes the overdrafted mass in an element with negative  !
!  water depth and takes that mass from the element designated in its data    !
!  structure.                                                                 !
!                                                                             !
!  Francisco Simoes, August 2007                                              !
!  Last updated (mm-dd-yyyy): 08-14-2007 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables.
  INTEGER :: i,j,k
  REAL (KIND=mp) :: vol

! The algorithm is repeated three times to reduce the chance of hanging cells.
  DO k = 1,3
    DO i = 1,n_elems
      IF (h(i) < zero) THEN
        vol = h(i)*grid(i)%area  ! Overdraft volume.
        j = grid(i)%o_cell
        h(j) = h(j) + vol/grid(j)%area
        h(i) = zero
      END IF
    END DO
  END DO

END SUBROUTINE overdraft
