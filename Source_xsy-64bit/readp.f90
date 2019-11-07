INTEGER FUNCTION readp(string)
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Function readp reads the first integer in string, then removes it from     !
!  string, shifting its contents to the left. It returns the integer read or  !
!  zero, if the reading failed.                                               !
!                                                                             !
!  F. Simoes, November 2004                                                   !
!  Last updated (mm-dd-yyyy): 11-16-2004 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

  ! Dummy variables:
  CHARACTER (LEN=*), INTENT(INOUT) :: string

  ! Local variables:
  INTEGER :: j

  IF (LEN_TRIM(string) == 0) THEN
    readp = 0
  ELSE
    READ(string,*) readp
    j = INDEX(string," ")
    string = string(j:)
    string = ADJUSTL(string)
  END IF

END FUNCTION readp
