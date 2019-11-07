SUBROUTINE alloc_err(i)
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Print-out of error message when allocation of variables fails.             !
!                                                                             !
!  Francisco Simoes, December 2003                                            !
!  Last updated (mm-dd-yyyy): 07-19-2005                                      !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  INTEGER, INTENT(IN) :: i

  WRITE (*,'(/"PROGRAM FAILURE: allocation of memory unsuccessful.",/, &
           & "This is a computer systems failure mostly likely due to",/, &
           & "insufficient memory or lack of other computing resources.",/, &
           & "Reducing the size of the SToRM run may help.",/, &
           & "System error no.",I5)') i
  CALL byebye("Program SToRM stopped!")

END SUBROUTINE alloc_err
