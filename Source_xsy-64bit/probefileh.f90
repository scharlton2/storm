INTEGER FUNCTION probefileh()
  USE options
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This function opens a file to write the probe output data.  The header is  !
!  set for filds having 12 characters in width (ES12.5).  It returns the      !
!  unit number used in I/O WRITE statements.                                  !
!                                                                             !
!  Francisco Simoes, October 2008                                             !
!  Last updated (mm-dd-yyyy): 10-23-2008 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables.
  INTEGER :: i,ierror
  CHARACTER (LEN=80) :: scratch
  LOGICAL, EXTERNAL :: get_iounit

! Check for a valid file name.
  IF (LEN_TRIM(probe_file) == 0) CALL byebye("Gage output filename missing.")

! Open file.
  IF (.NOT. get_iounit(i)) CALL byebye('ERROR in get_iounit(p_unit).')
  OPEN (i,FILE=TRIM(probe_file),IOSTAT=ierror)
  IF (ierror /= 0) THEN
    WRITE (*,'(" ERROR: unable to open output file ",A)') TRIM(probe_file)
    CALL byebye('SToRM stopped.')
  END IF
  probefileh = i

! Write header.
  WRITE (scratch,*) nprobe
  scratch = ADJUSTL(scratch)
  scratch = "('    Time    '," // TRIM(scratch) // &
    "(2X,I8,':h',2X,I8,':u',2X,I8,':v'))"
  WRITE (i,scratch)(probepts(i),probepts(i),probepts(i),i=1,nprobe)

! Set-up FORMAT string.
  WRITE (scratch,*) nprobe
  scratch = ADJUSTL(scratch)
  probe_format = "(ES12.4," // TRIM(scratch) // "(3ES12.4))"

END FUNCTION probefileh
