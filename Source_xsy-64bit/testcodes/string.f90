PROGRAM string
  IMPLICIT NONE

  ! This program tests readp.f90 and most of ioflow.f90.
  ! F. Simões, 11/16/2004

  REAL :: hbc,qin
  INTEGER :: counter,funit,line,n_h,h_nodes(20),n_qin,p,qin_edges(20)
  INTEGER, EXTERNAL :: readp
  CHARACTER (LEN=40) :: buffer
  LOGICAL :: stats
  LOGICAL, EXTERNAL :: readline

  line = 0
  funit = 12
  OPEN(funit,FILE="string_test.txt")

  ! Read inflow discharge boundary information.
  stats = readline(funit,buffer,line)
  READ(buffer,*) qin  ! Discharge at inlet.
  stats = readline(funit,buffer,line)
  READ(buffer,*) n_qin  ! Number of nodes in the inlet boundary.

  counter = 0
  DO WHILE (counter < n_qin)
    stats = readline(funit,buffer,line)  ! Read each line.
    IF (.NOT.stats) EXIT

    p = readp(buffer)
    DO WHILE (p > 0)  ! Read all the points in the same line.
      counter = counter + 1
      qin_edges(counter) = p
      p = readp(buffer)
      if (counter >= n_qin) EXIT
    END DO
  END DO

  ! Finally, read the fixed stage boundary nodes.
  stats = readline(funit,buffer,line)
  READ(buffer,*) hbc  ! Boundary water depth.
  stats = readline(funit,buffer,line)
  READ(buffer,*) n_h  ! Number of nodes with the specified water depth.

  counter = 0
  DO WHILE (counter < n_h)
    stats = readline(funit,buffer,line)  ! Read each line.
    IF (.NOT.stats) EXIT

    p = readp(buffer)
    DO WHILE (p > 0)  ! Read all the points in the same line.
      counter = counter + 1
      h_nodes(counter) = p
      p = readp(buffer)
      if (counter >= n_h) EXIT
    END DO
  END DO

  PRINT *,qin,hbc
  PRINT *,n_qin,n_h
  PRINT *,qin_edges
  PRINT *,h_nodes

END PROGRAM string

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
