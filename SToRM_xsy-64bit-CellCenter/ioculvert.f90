SUBROUTINE ioculvert
  USE options
  USE dep_vars
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This SUBROUTINE reads all the culvert information from an external file.   !
!  First, it sweeps the file to find the correct array sizes for memory       !
!  allocation; then, it allocates memory and reads each culvert coordinates   !
!  and performance curve (as a table) into the appropriate variables.         !
!                                                                             !
!  F. Simoes, January 2013                                                    !
!  Last updated (mm-dd-yyyy): 04-25-2013 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables:
  INTEGER :: counter,cvunit,i,ierror,j,k,line
  CHARACTER (LEN=40) :: buffer
  LOGICAL :: stats
  LOGICAL, EXTERNAL :: readline,get_iounit

  IF (get_iounit(cvunit)) THEN
    OPEN (cvunit,FILE=cvfile,STATUS='OLD',IOSTAT=ierror)
    IF (ierror /= 0) THEN
      WRITE (*,'(/" ERROR: unable to open culvert file ",A)') &
        TRIM(cvfile)
      CAll byebye('  Program SToRM stopped.')
    END IF
  ELSE
    WRITE (*,'(/" ERROR: unable to open culvert file ",A)') &
      TRIM(cvfile)
    CAll byebye('  Program SToRM stopped.')
  END IF

  line = 0
  counter = 0

  stats = readline(cvunit,buffer,line)  ! Number of culverts.
  IF (.NOT. stats) THEN
    PRINT *
    PRINT *,' ERROR reading number of culverts at line',line
    CALL byebye('  Program SToRM stopped.')
  END IF
  READ(buffer,*) nculvert

! Sweep over all culvert data to count maximum array dimensions for variable
! allocation.
  DO i = 1,nculvert
    stats = readline(cvunit,buffer,line)
    IF (.NOT. stats) THEN
      PRINT *
      PRINT *,' ERROR reading culvert coords at line',line
      CALL byebye('  Program SToRM stopped.')
    END IF
    stats = readline(cvunit,buffer,line)
    IF (.NOT. stats) THEN
      PRINT *
      PRINT *,' ERROR reading culvert coords at line',line
      CALL byebye('  Program SToRM stopped.')
    END IF
    stats = readline(cvunit,buffer,line)
    IF (.NOT. stats) THEN
      PRINT *
      PRINT *,' ERROR reading culvert data at line',line
      CALL byebye('  Program SToRM stopped.')
    END IF
    READ(buffer,*) k
    DO j = 1,k
      stats = readline(cvunit,buffer,line)
      IF (.NOT. stats) THEN
        PRINT *
        PRINT *,' ERROR reading culvert data at line',line
        CALL byebye('  Program SToRM stopped.')
      END IF
    END DO
    counter = MAX(counter,k)  ! Allocation size for culvert performance curve.
  END DO

  REWIND (cvunit)

! Allocate variables.
  ALLOCATE(cvptoin(3,nculvert),cvptoout(2,nculvert), ncvtable(nculvert), &
    cvhead(counter,nculvert),cvq(counter,nculvert),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(4*nculvert + 8*5*nculvert + 8*2*counter*nculvert)
  ALLOCATE(cvtrigin(nculvert),cvtrigout(nculvert),cvrin(2,nculvert), &
    cvsrc(nculvert),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(4*2*nculvert + 8*3*nculvert)

  line = 0
  stats = readline(cvunit,buffer,line)  ! Number of culverts.
  READ(buffer,*,IOSTAT=ierror) nculvert
  IF (ierror /= 0) THEN
    PRINT *,' ERROR IN LINE: ',line
    PRINT *,' Invalid format reading number of culverts in file.'
    CALL byebye('  Program SToRM stopped.')
  END IF
  DO i = 1,nculvert
    ! Read coordinates of inflow (inlet): x,y, and z.
    stats = readline(cvunit,buffer,line)
    READ(buffer,*,IOSTAT=ierror) cvptoin(1,i),cvptoin(2,i),cvptoin(3,i)
    IF (ierror /= 0) THEN
      PRINT *,' ERROR IN LINE: ',line
      PRINT *,' Invalid format reading coordinates of culvert from file.'
      CALL byebye('  Program SToRM stopped.')
    END IF
    ! Read coordinates of outflow (outlet): x and y only.
    stats = readline(cvunit,buffer,line)
    READ(buffer,*,IOSTAT=ierror) cvptoout(1,i),cvptoout(2,i)
    IF (ierror /= 0) THEN
      PRINT *,' ERROR IN LINE: ',line
      PRINT *,' Invalid format reading coordinates of culvert from file.'
      CALL byebye('  Program SToRM stopped.')
    END IF
    ! Read number of entries in performance curve table.
    stats = readline(cvunit,buffer,line)
    READ(buffer,*,IOSTAT=ierror) ncvtable(i)
    IF (ierror /= 0) THEN
      PRINT *,' ERROR IN LINE: ',line
      PRINT *,' Invalid format reading culvert performance curve.'
      CALL byebye('  Program SToRM stopped.')
    END IF
    IF (ncvtable(i) < 1) THEN
      PRINT *,' ERROR IN LINE: ',line
      PRINT *,' Too few entries in culvert performance curve.'
      CALL byebye('  Program SToRM stopped.')
    END IF
    DO j = 1,ncvtable(i)
      stats = readline(cvunit,buffer,line)
      READ(buffer,*,IOSTAT=ierror) cvhead(j,i),cvq(j,i)
      IF (ierror /= 0) THEN
        PRINT *,' ERROR IN LINE: ',line
        PRINT *,' Invalid format reading culvert performance curve.'
        CALL byebye('  Program SToRM stopped.')
      END IF
    END DO
  END DO

  CLOSE (cvunit)

END SUBROUTINE ioculvert
