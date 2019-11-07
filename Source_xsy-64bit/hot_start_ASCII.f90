SUBROUTINE hot_start_ASCII
  USE parameters
  USE geometry
  USE dep_vars
  USE io
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  A temporary subroutine to read cell centered values of the dependent       !
!  variables from a Tecplot compatible ASCII file.  This subroutine checks    !
!  if the solution file exists and, if it does, reads it.  Should only be     !
!  used if hot start is used with CGNS I/O, and was developed for iRIC 2.0.   !
!                                                                             !
!  Francisco Simoes, October 2012                                             !
!  Last updated (mm-dd-yyyy): 11-03-2016 by F. Simões                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables.
  INTEGER :: i,ierror,j,k,lunit,pto
  REAL (KIND=mp) :: work
  CHARACTER (LEN=256) :: scratch  ! Same length as hsfile.
  LOGICAL :: f
  LOGICAL, EXTERNAL :: get_iounit

! Use CGNS hot-start file name to find the ASCII file to load.
  i = INDEX(hsfile,'\',.TRUE.)
  IF (i == 0) CALL byebye('Error in hot-start ASCII file name.')
  scratch = hsfile(1:i) // "Tecplot_sol_file.dat"

! Check if file exists and RETURN if it doesn't.
  INQUIRE (FILE=scratch,EXIST=f)
  IF (.NOT. f) RETURN

! Open file and read cell-centered values of the dependent variables.
  IF (.NOT. get_iounit(lunit)) CALL byebye('ERROR in get_iounit(lunit).')
  OPEN (lunit,FILE=scratch,IOSTAT=ierror)
  IF (ierror /= 0) THEN
    WRITE (*,'(" ERROR: unable to open CGNS-Tecplot hot-start file.")')
    CALL byebye('SToRM stopped.')
  END IF

! Read titles.
  DO i = 1,3
    READ (lunit,'(A)') scratch
  END DO

! Skip x, y, and z.
  DO i = 1,n_pts*3
    READ (lunit,*) work
  END DO

! Read h, u, v, and cd.
  DEALLOCATE(u_tcp,v_tcp,h_tcp)
  ALLOCATE(u_tcp(n_elems),v_tcp(n_elems),h_tcp(n_elems),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  DO i = 1,n_elems
    READ (lunit,*) h_tcp(i)
  END DO
  DO i = 1,n_elems
    READ (lunit,*) u_tcp(i)
  END DO
  DO i = 1,n_elems
    READ (lunit,*) v_tcp(i)
  END DO

! Reset datav values from TRUE (cell vertex in CGNS) to FALSE (cell centered).
  datav(4) = .FALSE.  ! h
  datav(5) = .FALSE.  ! u
  datav(6) = .FALSE.  ! v

  CLOSE (lunit)

END SUBROUTINE hot_start_ASCII
