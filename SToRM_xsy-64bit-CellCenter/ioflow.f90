SUBROUTINE ioflow(funit)
  USE options
  USE dep_vars
  USE vbc_arrays
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Subroutine ioflow reads the flow boundary conditions from an external      !
!  file and prepares the corresponding data structures used by SToRM.         !
!                                                                             !
!  F. Simoes, November 2008                                                   !
!  Last updated (mm-dd-yyyy): 10-01-2011 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables:
  INTEGER, INTENT(IN) :: funit

! Local variables:
  INTEGER :: bdrtype,counter,i,ierror,j,k,kq,kh,kn,line,maxh,maxnin,maxnout, &
             maxnchan,maxqin,p
  INTEGER, EXTERNAL :: is_edge,readp
  CHARACTER (LEN=40) :: buffer
  LOGICAL :: stats
  LOGICAL, EXTERNAL :: readline

!-----------------------------------------------------------------------------!
!                                                                             !
!  First loop: reads the entire file, counting the number of inflow and       !
!  outflow boundaries.  Free boundaries, where no boundary conditions are     !
!  applied, are considered outflow boundaries with htype = 0.                 !
!                                                                             !
!-----------------------------------------------------------------------------!

  line = 0
  n_inflowbdr = 0
  n_outflowbdr = 0
  n_chanbdr = 0
  maxqin = 0  ! Max number of inflow time series entries.
  maxh = 0  ! Max number of outflow time series entries.
  maxnin = 0  ! Max number of inflow nodes.
  maxnout = 0  ! Max number of outflow nodes.
  maxnchan = 0  ! Max number ot type 3 nodes (temporary Buckeye Lake nodes).

  DO
    stats = readline(funit,buffer,line)
    IF (.NOT. stats) EXIT  ! Finished reading all the file.
    READ(buffer,*) bdrtype,i,j
    IF (bdrtype == 1) THEN
      n_inflowbdr = n_inflowbdr + 1
      maxqin = MAX(j,maxqin)
    ELSE IF (bdrtype == 2) THEN
      n_outflowbdr = n_outflowbdr + 1
      maxh = MAX(j,maxh)
    ELSE IF (bdrtype == 3) THEN
      n_chanbdr = n_chanbdr + 1
    ELSE
      PRINT *
      PRINT *,' ERROR: invalid I/O boundary type at line',line
      CALL byebye('  Program SToRM stopped.')
    END IF

    ! Loop over the lines corresponding to the time series.
    ! If boundary is outflow with htype = 0 (free boundary) or = 2 (critical
    ! outflow), there is no time series.
    IF ((bdrtype == 2 .AND. i /= 1) .OR. bdrtype == 3) THEN
      j = 0
    END IF
    DO k = 1,j
      stats = readline(funit,buffer,line)
      IF (.NOT. stats) THEN
        PRINT *
        PRINT *,' ERROR: number of time series entries at line',line
        CALL byebye('  Program SToRM stopped.')
      END IF
    END DO

    ! Read the nodes belonging to the boundary.
    stats = readline(funit,buffer,line)
    IF (.NOT. stats) THEN
      PRINT *
      PRINT *,' ERROR: reading bounary nodes at line',line
      CALL byebye('  Program SToRM stopped.')
    END IF
    READ(buffer,*) k  ! Number of nodes in the boundary.
    IF (k < 2) THEN
      PRINT *
      PRINT *,' ERROR: number of nodes < 2 at line',line
      CALL byebye('  Program SToRM stopped.')
    END IF
    IF (bdrtype == 1) THEN
      maxnin = MAX(k,maxnin)
    ELSE IF (bdrtype == 2) THEN
      maxnout = MAX(k,maxnout)
    ELSE IF (bdrtype == 3) THEN
      maxnchan = MAX(k,maxnchan)
    END IF

    ! The following nested DO WHILE loops read all the nodes that belong to the
    ! boundary.  There may be multiple nodes per line as long as they are
    ! separated by one or more blank spaces.
    counter = 0
    DO WHILE (counter < k)
      stats = readline(funit,buffer,line)  ! Read each line.
      p = readp(buffer)
      DO WHILE (p > 0)  ! Read all the points in the same line.
        counter = counter + 1
        p = readp(buffer)
        IF (counter >= k) EXIT
      END DO
    END DO

  END DO
  REWIND (funit)

!------------------------------------------------------------------------------
!                                                                             !
!  ALLOCATE arrays and prepare to read the data.                              !
!                                                                             !
!------------------------------------------------------------------------------

  ALLOCATE(hbc(n_outflowbdr),htype(n_outflowbdr),n_timeserh(n_outflowbdr), &
    timeloch(maxh,n_outflowbdr),timeserh(maxh,n_outflowbdr),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(2*4*n_outflowbdr + 8*n_outflowbdr + 2*8*maxh*n_outflowbdr)
  ALLOCATE(qin(n_inflowbdr),qin2(n_inflowbdr),vtype(n_inflowbdr), &
    n_timeserq(n_inflowbdr),timelocq(maxqin,n_inflowbdr), &
    timeserq(maxqin,n_inflowbdr),timeserq2(maxqin,n_inflowbdr),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(2*4*n_inflowbdr + 2*8*n_inflowbdr + 3*8*maxqin*n_inflowbdr)
  ALLOCATE(n_qin(n_inflowbdr),n_hbc(n_outflowbdr), &
    qin_nodes(maxnin,n_inflowbdr),hbc_nodes(maxnout,n_outflowbdr),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(4*n_inflowbdr*n_outflowbdr)
  CALL mem_add(8*maxnin*n_inflowbdr + 8*maxnout*n_outflowbdr)
  ALLOCATE(n_channel(n_chanbdr),channel(maxnchan,n_chanbdr),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(4*n_chanbdr + 4*n_chanbdr*maxnchan)

! Initialize some of the arrays to -1.  This may be inportant later, when using
! function in_array2 to check for wall edges.
  DO j = 1,n_inflowbdr
    DO i = 1,maxnin
      qin_nodes(i,j) = -1
    END DO
  END DO
  DO j = 1,n_outflowbdr
    DO i = 1,maxnout
      hbc_nodes(i,j) = -1
    END DO
  END DO
  DO j = 1,n_chanbdr
    DO i = 1,maxnchan
      channel(i,j) = -1
    END DO
  END DO

! Allocate some of the arrays used to keep track of computed inflows and
! outflows during the time stepping of the solution.
  ALLOCATE(cqin(n_inflowbdr),cqin2(n_inflowbdr),cqin3(n_inflowbdr), &
    cqout(n_outflowbdr),cqout2(n_outflowbdr),cqout3(n_outflowbdr),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(3*8*n_inflowbdr + 3*8*n_outflowbdr)

!------------------------------------------------------------------------------
!                                                                             !
!  Read in the data.  For each boundary, read the type, the time table, and   !
!  the string of nodes belonging to the boundary.  Boundary types can be      !
!  mixed, but all the information corresponding to each boundary must be      !
!  grouped.                                                                   !
!                                                                             !
!------------------------------------------------------------------------------

  line = 0
  kq = 0
  kh = 0
  kn = 0
  DO k = 1,n_inflowbdr + n_outflowbdr + n_chanbdr  ! Read all data groups.
    stats = readline(funit,buffer,line)
    IF (.NOT. stats) EXIT
    READ(buffer,*) bdrtype,i

    IF (bdrtype == 1) THEN  ! It's an inflow boundary -------------------------
      READ(buffer,*) bdrtype,i,j
      kq = kq + 1
      vtype(kq) = i
      n_timeserq(kq) = j
      counter = 0

      ! First read the time series.
      ! The time series containing the boundary conditions (both inflow and
      ! outflow) are stored in arrays, with one entry for each time value.  The
      ! times corresponding to each array entry are stored in arrays timeloch()
      ! and timelocq() for the outflow and the inflow boundaries, respectively.
      ! Arrays timeserq() and timeserq2() contain the inflow conditions and
      ! array timeserh() contains the outflow conditions.
      DO WHILE (counter < n_timeserq(kq))
        stats = readline(funit,buffer,line)  ! Read each line.
        IF (.NOT. stats) EXIT
        counter = counter + 1

        IF (vtype(kq) == 3) THEN  ! Both velocity components are specified
          READ(buffer,*,IOSTAT=ierror) timelocq(counter,kq), &
                                     timeserq(counter,kq),timeserq2(counter,kq)
          IF (ierror /= 0) THEN
            PRINT *,' ERROR IN LINE: ',line
            PRINT *,' Invalid format reading inflow boundary condition values.'
            CALL byebye('  Program SToRM stopped.')
          END IF

        ELSE  ! Only one velocity or discharge value used.
          READ(buffer,*,IOSTAT=ierror) timelocq(counter,kq), &
                                       timeserq(counter,kq)
          IF (ierror /= 0) THEN
            PRINT *,' ERROR IN LINE: ',line
            PRINT *,' Invalid format reading inflow boundary condition values.'
            CALL byebye('  Program SToRM stopped.')
          END IF
        END IF
      END DO

      ! Do a quick check of the data.
      DO i = 1,n_timeserq(kq) - 1
        IF (timelocq(i+1,kq) < timelocq(i,kq)) THEN
          PRINT *,'ERROR: boundary conditions time series must be ordered'
          PRINT *,'with increasing time data. Program SToRM stopped.'
          CALL byebye('')
        END IF
      END DO

      ! Now read the nodes belonging to the inflow boundary.
      stats = readline(funit,buffer,line)
      READ(buffer,*) n_qin(kq)  ! Number of nodes in the inlet boundary.
      IF (n_qin(kq) <= 1) CALL byebye('ERROR: n_qin <= 1.')

      ! qin_nodes(,) contains the nodes in the inlet boundary.
      ! The following nested DO WHILE loop reads all the nodes that belong to
      ! the inlet boundary.  There may be multiple nodes per line as long as
      ! they are separated by one or more blank spaces.
      counter = 0
      DO WHILE (counter < n_qin(kq))
        stats = readline(funit,buffer,line)  ! Read each line.
        IF (.NOT. stats) EXIT

        p = readp(buffer)
        DO WHILE (p > 0)  ! Read all the points in the same line.
          counter = counter + 1
          qin_nodes(counter,kq) = p
          p = readp(buffer)
          IF (counter >= n_qin(kq)) EXIT
        END DO
      END DO

    ELSE IF (bdrtype == 2) THEN  ! It's an outflow boundary -------------------
      kh = kh + 1
      htype(kh) = i
      IF (i == 1) THEN
        READ(buffer,*) bdrtype,i,j
      ELSE
        j = 0
      END IF
      n_timeserh(kh) = j
      counter = 0

      ! First read the time series, but only for regular subcritical outflow
      ! boundary, with htype() = 1.
      DO WHILE (counter < n_timeserh(kh))
        stats = readline(funit,buffer,line)  ! Read each line.
        IF (.NOT. stats) EXIT
        counter = counter + 1

        READ(buffer,*,IOSTAT=ierror) timeloch(counter,kh),timeserh(counter,kh)
        IF (ierror /= 0) THEN
          PRINT *,' ERROR IN LINE: ',line
          PRINT *,' Invalid format reading outflow boundary condition values.'
          CALL byebye('  Program SToRM stopped.')
        END IF
      END DO

      ! Do a quick check of the data.
      DO i = 1,n_timeserh(kh) - 1
        IF (timeloch(i+1,kh) < timeloch(i,kh)) THEN
          PRINT *,'ERROR: boundary conditions time series must be ordered'
          PRINT *,'with increasing time data. Program SToRM stopped.'
          CALL byebye('')
        END IF
      END DO

      ! Now read the nodes belonging to the outflow boundary.
      stats = readline(funit,buffer,line)
      READ(buffer,*) n_hbc(kh)  ! Number of nodes in the outlet boundary.
      IF (n_hbc(kh) <= 1) CALL byebye('ERROR: n_hbc <= 1.')

      ! hbc_nodes(,) contains the nodes with specified water depth, as well as
      ! all the free boundary nodes.
      ! The following nested DO WHILE loop reads all the nodes that belong to
      ! the outlet boundary.  There may be multiple nodes per line as long as
      ! they are separated by one or more blank spaces.
      counter = 0
      DO WHILE (counter < n_hbc(kh))
        stats = readline(funit,buffer,line)  ! Read each line.
        IF (.NOT. stats) EXIT

        p = readp(buffer)
        DO WHILE (p > 0)  ! Read all the points in the same line.
          counter = counter + 1
          hbc_nodes(counter,kh) = p
          p = readp(buffer)
          IF (counter >= n_hbc(kh)) EXIT
        END DO
      END DO

    ELSE IF (bdrtype == 3) THEN  ! It's a channel -----------------------------
      ! Read the nodes belonging to a channel.
      kn = kn + 1
      stats = readline(funit,buffer,line)
      READ(buffer,*) n_channel(kn)  ! Number of nodes in the channel.
      IF (n_channel(kn) <= 1) CALL byebye('ERROR: n_channel <= 1.')

      ! channel(,) contains the nodes of each specified channel.  The following
      ! nested DO WHILE loop reads all the nodes that belong to the channel.
      ! There may be multiple nodes per line as long as they are separated by
      ! one or more blank spaces.
      counter = 0
      DO WHILE (counter < n_channel(kn))
        stats = readline(funit,buffer,line)  ! Read each line.
        IF (.NOT. stats) EXIT

        p = readp(buffer)
        DO WHILE (p > 0)  ! Read all the points in the same line.
          counter = counter + 1
          channel(counter,kn) = p
          p = readp(buffer)
          IF (counter >= n_channel(kn)) EXIT
        END DO
      END DO

    END IF
  END DO

END SUBROUTINE ioflow
