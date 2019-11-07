LOGICAL FUNCTION readline(funit,line,lineno)
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This subroutine reads the input file line by line.  The behavior is the    !
!  following:                                                                 !
!                                                                             !
!  - all blank spaces are removed from the beginning (left) of the line;      !
!  - if the first non-blank character is # the line is a comment line,        !
!    therefore it is skipped without further processing and the next line is  !
!    read;                                                                    !
!  - blank lines are skipped and the next line is read;                       !
!  - all other lines are converted to upper case and returned to the calling  !
!    program.                                                                 !
!                                                                             !
!  Note that only the first 40 characters of each line are read-in, which     !
!  means that all data in column 41 and following are automatically           !
!  discarded.                                                                 !
!                                                                             !
!  The function returns .TRUE. if the line was read successfuly, and .FALSE.  !
!  if an error or an End-Of-File was encountered.                             !
!                                                                             !
!  Francisco Simoes, October 2004                                             !
!  Last updated (mm-dd-yyyy): 11-16-2004                                      !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments:
  INTEGER, INTENT(IN) :: funit  ! File unit number for READ statements.
  INTEGER, INTENT(INOUT) :: lineno  ! Line number of line read, for debugging.
  CHARACTER(LEN=40), INTENT(OUT) :: line  ! Line read from file.

! Local variables:
  INTEGER :: i
  CHARACTER, EXTERNAL :: to_upper

  DO WHILE (.TRUE.)

    line = ''
    lineno = lineno + 1
    READ (funit,'(A)',ERR = 1,END = 1) line
    line = ADJUSTL(line)

    ! Comment lines are ignored.
    IF (line(1:1) == '#') CYCLE

    ! Blank lines are skipped.
    IF (LEN_TRIM(line) == 0) CYCLE

    ! All other lines are processed to convert from lower to upper case, and
    ! then passed to the calling program.
    DO i = 1,LEN_TRIM(line)
      line(i:i) = to_upper(line(i:i))
    END DO

    readline = .TRUE.
    RETURN  ! Normal return to calling program.

  END DO

! Error or end-of-file encountered.
1 readline = .FALSE.

END FUNCTION readline
