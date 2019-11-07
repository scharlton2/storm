SUBROUTINE err_opts(string,lineno)
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Prints an error message tailored to the options file and stops the run.    !
!                                                                             !
!  Francisco Simoes, October 2004                                             !
!  Last updated (mm-dd-yyyy): 11-11-2005                                      !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  INTEGER, INTENT(IN) :: lineno
  CHARACTER (LEN=*) :: string

! Local variables.
  CHARACTER (LEN=40) :: buffer

  WRITE (buffer,*) lineno
  buffer = ADJUSTL(buffer)
  PRINT *,''
  PRINT *,'ERROR: illegal value in option ',string
  PRINT *,'in line ',TRIM(buffer),' of the options file.'
  CALL byebye('Program SToRM stopped.')

END SUBROUTINE err_opts
