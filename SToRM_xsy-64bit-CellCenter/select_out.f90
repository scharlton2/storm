LOGICAL FUNCTION select_out(iter)
  USE io
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This function simply determines if it is time to print-out an              !
!  intermediate solution by comparing the iteration value, 'iter', and the    !
!  contents of array iout().  For this subroutine to work, the values of      !
!  iout() must be ordered in increasing value and the dimension of iout()     !
!  must be >= 1.                                                              !
!                                                                             !
!  Francisco Simoes, November 2005                                            !
!  Last updated (mm-dd-yyyy): 11-28-2005 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  INTEGER, INTENT(IN) :: iter

  select_out = .FALSE.
  IF (iter == iout(p_iout)) THEN
    select_out = .TRUE.
    p_iout = p_iout + 1
    p_iout = MIN(p_iout,n_iout)  ! To avoid array-out-of-bounds error.
  END IF

END FUNCTION select_out
