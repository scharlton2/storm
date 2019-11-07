SUBROUTINE rvegdata
  USE options
  USE geometry
  USE dep_vars
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Read the roughness coefficients for the vegetation flow resistance model.  !
!  This is hardcoded to use the vegetation drag model of Tsihrintzis (2001).  !
!  The first data line must be the number of nodes (blank lines and lines     !
!  starting with # are skipped).  Then, each line contains the values of      !
!  gamma and kappa of eq. (18) of Tsihrintzis (2001).                         !
!                                                                             !
!  Francisco Simoes, November 2007                                            !
!  Last updated (mm-dd-yyyy): 02-02-2009                                      !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables.
  INTEGER :: f,i,ierror,lineno,n
  REAL (KIND=mp), ALLOCATABLE, DIMENSION(:) :: tmp1,tmp2
  CHARACTER(LEN=40) :: buffer
  LOGICAL, EXTERNAL :: get_iounit,readline

  IF (.NOT. get_iounit(f)) THEN
    WRITE (*,'(" ERROR: unable to open file ",A)') TRIM(f_fric_facts)
    CALL byebye('SToRM stopped.')
  END IF
  OPEN (f,FILE=f_fric_facts,STATUS='OLD',IOSTAT=ierror)
  IF (ierror /= 0) THEN
    WRITE (*,'(" ERROR: unable to open file ",A)') TRIM(f_fric_facts)
    CALL byebye('SToRM stopped.')
  END IF

  lineno = 0
  IF (readline(f,buffer,lineno)) THEN
    READ (buffer,*) n
  ELSE
    WRITE (*,'(" Error reading file ",A)') TRIM(f_fric_facts)
    WRITE (*,'(" in line",I4)') lineno
    CALL byebye('SToRM stopped.')
  END IF

  ALLOCATE(tmp1(n),tmp2(n),STAT=ierror)  ! Temporary working arrays.
  IF (ierror /= 0) CALL alloc_err(ierror)

! Read the data:
  DO i = 1,n
    IF (readline(f,buffer,lineno)) THEN
      READ (buffer,*) tmp1(i),tmp2(i)
    ELSE
      WRITE (*,'(" Error reading file ",A)') TRIM(f_fric_facts)
      WRITE (*,'(" in line",I4)') lineno
      CALL byebye('SToRM stopped.')
    END IF
  END DO

  SELECT CASE (opt_solver)

!-----------------------------------------------------------------------------!
!                                 RDS solver                                  !
!-----------------------------------------------------------------------------!
  CASE (1)  ! Residual distribution scheme options.
    ALLOCATE(rcoef1(n_pts),rcoef2(n_pts),STAT=ierror)  ! Allocate global vars.
    IF (ierror /= 0) CALL alloc_err(ierror)
    CALL mem_add(8*n_pts*2)

    IF (n == n_elems) THEN
      ! Data is cell-centered: interpolate to vertices.
      CALL ctr2vtx(nctr2vtx,tmp1,rcoef1)
      CALL ctr2vtx(nctr2vtx,tmp2,rcoef2)

    ELSE IF (n == n_pts) THEN
      ! Data is good: simply read the variables.
      rcoef1 = tmp1
      rcoef2 = tmp2

    ELSE
      WRITE (*,'(" Error reading file ",A)') TRIM(f_fric_facts)
      WRITE (*,'(" Invalid value in line",I4)') lineno
      CALL byebye('SToRM stopped.')
    END IF

!-----------------------------------------------------------------------------!
!                                 FVT solver                                  !
!-----------------------------------------------------------------------------!
  CASE (2,3)  ! Centered finite volume solvers.
    ALLOCATE(rcoef1(n_elems),rcoef2(n_elems),STAT=ierror)  ! Global arrays.
    IF (ierror /= 0) CALL alloc_err(ierror)
    CALL mem_add(8*n_elems*2)

    IF (n == n_elems) THEN
      ! Data is good: simply read the variables.
      rcoef1 = tmp1
      rcoef2 = tmp2

    ELSE IF (n == n_pts) THEN
      ! Data is at vertices: interpolate to cell center.
      CALL vtx2ctr(tmp1,rcoef1)
      CALL vtx2ctr(tmp2,rcoef2)

    ELSE
      WRITE (*,'(" Error reading file ",A)') TRIM(f_fric_facts)
      WRITE (*,'(" Invalid value in line",I4)') lineno
      CALL byebye('SToRM stopped.')
    END IF

  END SELECT

! Free the memory used for the temporary arrays and close datafile.
  DEALLOCATE(tmp1,tmp2)
  CLOSE (f)

END SUBROUTINE rvegdata
