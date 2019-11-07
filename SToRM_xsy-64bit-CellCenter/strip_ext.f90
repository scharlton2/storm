SUBROUTINE strip_ext(fn)
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Strips the extension (including the dot) from a file name.  The dummy      !
!  argument variable 'fn' contains the file name on input.  On output, 'fn'   !
!  contains the same file name but without the extension, truncated from the  !
!  last '.' in the filename (the '.' is removed, too).  If 'fn' doens't have  !
!  any extension, 'fn' is returned untouched.                                 !
!                                                                             !
!  Francisco Simoes, November 2005                                            !
!  Last updated (mm-dd-yyyy): 11-28-2005 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  CHARACTER (LEN=*), INTENT(INOUT) :: fn

! Local variables.
  INTEGER :: i

  fn = ADJUSTL(fn)
  i = INDEX(fn,'.',.TRUE.)
  IF (i /= 0) fn = fn(1:i-1)

END SUBROUTINE strip_ext
