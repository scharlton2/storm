PROGRAM strip
  USE dflib ! For functions GETARG and NARGS of Microsoft Fortran Library.
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Program strip  -  version 0.9                                              !
!                                                                             !
!  The objective of this program is to remove all blank and comment lines,    !
!  as well as any other comments mixed in with the code, from Fortran source  !
!  files.  It works by reading a list of file names from file                 !
!  'filenames.txt' and processing those files sequentially.  The result is    !
!  regurgitated into a single output file that contains all the input files   !
!  sequentially concatenated.  The name of the output file can be given in    !
!  the command line, e.g.,                                                    !
!                                                                             !
!     strip outputfile.f90                                                    !
!                                                                             !
!  If not given in the command line, the program asks for it.                 !
!                                                                             !
!  Note that the output file should have the proper Fortran extension (.f90   !
!  or .f95).  Note alse that, in order for the program to compile             !
!  successfully, the module and subroutine files must appear in a propper     !
!  sequence.  That must be accomplished by placing the file names in file     !
!  'filenames.txt' in the required sequence.                                  !
!                                                                             !
!  Program strip will ignore all files (in the 'filenames.txt' file) that do  !
!  not have a Fortran extension (.f90 or .f95) or that have the string        !
!  '.old.' in their file name.  The program is case sensitive and uses only   !
!  lower case characters.                                                     !
!                                                                             !
!  Francisco Simoes, March 2005                                               !
!  Last updated (mm-dd-yyyy): 03-04-2005                                      !
!                                                                             !
!-----------------------------------------------------------------------------!

  INTEGER :: i,ierror
  CHARACTER (LEN=80) :: fname,line,outfile

  ! Parse command line to find the name for the output file.
  IF (NARGS() == 1) THEN
    WRITE (*,'(" Enter the OUTPUT file name: ",$)')  ! The $ edit descriptor is
    READ (*,'(A)') outfile                           ! not Fortran 90 standard.
    outfile = ADJUSTL(outfile)
  ELSE
    CALL GETARG(1,outfile)
    outfile = ADJUSTL(outfile)
  END IF
  OPEN (15,FILE=outfile,ACTION='WRITE')
  REWIND (15)

  ! Open file containing all the filenames of the files to be stripped and
  ! concatenated.
  OPEN (10,FILE='filenames.txt',STATUS='OLD',ACTION='READ',IOSTAT=ierror)
  IF (ierror /= 0) THEN
    WRITE (*,'(" ERROR: unable to open file ",A)') 'filenames.txt'
    WRITE (*,'(" Program strip stopped!")')
    STOP
  END IF

  DO  ! Main DO-loop over all the files in file 'filenames.txt'.

    READ (10,'(A)',IOSTAT=ierror) fname
    IF (ierror < 0) EXIT  ! EOF encountered.
    IF (ierror > 0) THEN  ! Error encountered.
      WRITE (*,'(" ERROR: unable to read file ",A)') 'filenames.txt'
      WRITE (*,'(" Program strip stopped!")')
      STOP
    END IF
    fname = '..\' // ADJUSTL(fname)
    fname = TRIM(fname)

    ! If a file contains the string '.old.' it means that the file is an older
    ! version of an existing file and should not be processed.
    IF (INDEX(fname,'.old.') /= 0) THEN
      PRINT *,'*** Skipping file ' // TRIM(fname)
      CYCLE
    END IF

    ! If a file does not end in the suffix '.f90' or '.f95' it should not be
    ! processed.
    IF (INDEX(fname,'.f90',.TRUE.)+3 /= LEN_TRIM(fname) .AND. &
      INDEX(fname,'.f95',.TRUE.)+3 /= LEN_TRIM(fname)) THEN
      PRINT *,'--- Skipping file ' // TRIM(fname)
      CYCLE
    END IF

    ! Open the next file in the sequence.
    PRINT *,'Processing file ' // TRIM(fname)
    OPEN (12,FILE=fname,STATUS='OLD',ACTION='READ',IOSTAT=ierror)
    IF (ierror /= 0) THEN
      PRINT *,'ERROR: unable to open file ' // TRIM(fname)
      CYCLE
    END IF

    DO  ! Inner loop to read and process each file.

      READ (12,'(A)',IOSTAT=ierror) line
      IF (ierror < 0) EXIT  ! EOF encountered.
      IF (ierror > 0) THEN  ! Error encountered.
        PRINT *,' ERROR: unable to read file ' // TRIM(fname)
        EXIT
      END IF

      line = ADJUSTL(line)
      IF (LEN_TRIM(line) == 0 .OR. line(1:1) == '!') CYCLE

      i = INDEX(line,'!')
      IF ( i > 0) line = line(1:i-1)

      WRITE (15,'(A)') TRIM(line)

    END DO  ! Inner loop, for each file.

    CLOSE (12)

  END DO  ! Outer loop, over all the files.

  CLOSE (10)
  CLOSE (15)

END PROGRAM strip
