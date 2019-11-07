LOGICAL FUNCTION fnames(funit,fname)
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This subroutine reads file unit funit line by line until the desired data  !
!  is found.  The behavior is the following:                                  !
!                                                                             !
!  - all blank spaces are removed from the beginning (left) of the line;      !
!  - if the first non-blank character is # the line is a comment line,        !
!    therefore it is skipped without further processing and the next line is  !
!    read;                                                                    !
!  - blank lines are skipped and the next line is read;                       !
!  - all other lines passed to variable fname without further processing.     !
!                                                                             !
!  Note that only the first 80 characters of each line are read-in, which     !
!  means that all data in column 81 and following are automatically           !
!  discarded.                                                                 !
!                                                                             !
!  The function returns .TRUE. if the line was read successfuly, and .FALSE.  !
!  if an error or an End-Of-File was encountered.                             !
!                                                                             !
!  Francisco Simoes, November 2004                                            !
!  Last updated (mm-dd-yyyy): 11-16-2004                                      !
!                                                                             !
!-----------------------------------------------------------------------------!

  ! Dummy arguments:
  INTEGER, INTENT(IN) :: funit
  CHARACTER (LEN=80), INTENT(OUT) :: fname

  DO WHILE (.TRUE.)

    fname = ''
    READ (funit,'(A)',ERR=1,END=1) fname
    fname = ADJUSTL(fname)

    ! Comment lines are ignored.
    IF (fname(1:1) == '#') CYCLE

    ! Blank lines are skipped.
    IF (LEN_TRIM(fname) == 0) CYCLE

    ! Initialize file name variable and return successfuly.
    fnames = .TRUE.
    RETURN

  END DO

! Error or end-of-file encountered.
1 fnames = .FALSE.

END FUNCTION fnames
