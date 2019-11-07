SUBROUTINE windforcing
  USE options
  USE geometry
  USE dep_vars
  USE constants
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Read the necessary parameters to compute the wind forcing terms.           !
!                                                                             !
!  The first entry is the number of lines, which must be equal to n_pts or    !
!  to n_elems (blank lines and lines starting with # are skipped).  The data  !
!  is organized in three columns.  The first column must contain the product  !
!  of Cd with RHOa, where Cd is the usual wind shear coefficient at 10 m and  !
!  RHOa is the density of air.  The second and third columns have the x and   !
!  y components of the wind velocity.                                         !
!                                                                             !
!  Francisco Simoes, November 2007                                            !
!  Last updated (mm-dd-yyyy): 02-02-2009                                      !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables.
  INTEGER :: f,i,ierror,lineno,n
  REAL (KIND=mp) :: windmag
  REAL (KIND=mp), ALLOCATABLE, DIMENSION(:) :: tmp1,tmp2,tmp3,wcoef0
  CHARACTER(LEN=40) :: buffer
  LOGICAL, EXTERNAL :: get_iounit,readline

  IF (.NOT. get_iounit(f)) THEN
    WRITE (*,'(" ERROR: unable to open file ",A)') TRIM(wind_file)
    CALL byebye('SToRM stopped.')
  END IF
  OPEN (f,FILE=wind_file,STATUS='OLD',IOSTAT=ierror)
  IF (ierror /= 0) THEN
    WRITE (*,'(" ERROR: unable to open file ",A)') TRIM(wind_file)
    CALL byebye('SToRM stopped.')
  END IF

  lineno = 0
  IF (readline(f,buffer,lineno)) THEN
    READ (buffer,*) n
  ELSE
    WRITE (*,'(" Error reading file ",A)') TRIM(wind_file)
    WRITE (*,'(" in line",I4)') lineno
    CALL byebye('SToRM stopped.')
  END IF

  ALLOCATE(tmp1(n),tmp2(n),tmp3(n),STAT=ierror)  ! Temporary working arrays.
  IF (ierror /= 0) CALL alloc_err(ierror)

! Read the data:
  DO i = 1,n
    IF (readline(f,buffer,lineno)) THEN
      READ (buffer,*) tmp1(i),tmp2(i),tmp3(i)
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
    ! Allocate global vars.
    ALLOCATE(wcoef0(n_pts),wcoef1(n_pts),wcoef2(n_pts),STAT=ierror)
    IF (ierror /= 0) CALL alloc_err(ierror)
    CALL mem_add(8*n_pts*2)

    IF (n == n_elems) THEN
      ! Data is cell-centered: interpolate to vertices.
      CALL ctr2vtx(nctr2vtx,tmp1,wcoef0)
      CALL ctr2vtx(nctr2vtx,tmp2,wcoef1)
      CALL ctr2vtx(nctr2vtx,tmp3,wcoef2)

    ELSE IF (n == n_pts) THEN
      ! Data is good: simply read the variables.
      wcoef0 = tmp1
      wcoef1 = tmp2
      wcoef2 = tmp3

    ELSE
      WRITE (*,'(" Error reading file ",A)') TRIM(f_fric_facts)
      WRITE (*,'(" Invalid value in line",I4)') lineno
      CALL byebye('SToRM stopped.')
    END IF

    ! Now compute the wind forcing terms.
    DO i = 1,n_pts
      windmag = SQRT(wcoef1(i)*wcoef1(i) + wcoef2(i)*wcoef2(i))
      wcoef1(i) = wcoef0(i)*wcoef1(i)*windmag/rho
      wcoef2(i) = wcoef0(i)*wcoef2(i)*windmag/rho
    END DO

!-----------------------------------------------------------------------------!
!                                 FVT solver                                  !
!-----------------------------------------------------------------------------!
  CASE (2,3)  ! Centered finite volume solvers.
    ! Allocate global vars.
    ALLOCATE(wcoef0(n_elems),wcoef1(n_elems),wcoef2(n_elems),STAT=ierror)
    IF (ierror /= 0) CALL alloc_err(ierror)
    CALL mem_add(8*n_elems*2)

    IF (n == n_elems) THEN
      ! Data is good: simply read the variables.
      wcoef0 = tmp1
      wcoef1 = tmp2
      wcoef2 = tmp3

    ELSE IF (n == n_pts) THEN
      ! Data is at vertices: interpolate to cell center.
      CALL vtx2ctr(tmp1,wcoef0)
      CALL vtx2ctr(tmp2,wcoef1)
      CALL vtx2ctr(tmp3,wcoef2)

    ELSE
      WRITE (*,'(" Error reading file ",A)') TRIM(f_fric_facts)
      WRITE (*,'(" Invalid value in line",I4)') lineno
      CALL byebye('SToRM stopped.')
    END IF

    ! Now compute the wind forcing terms.
    DO i = 1,n_elems
      windmag = SQRT(wcoef1(i)*wcoef1(i) + wcoef2(i)*wcoef2(i))
      wcoef1(i) = wcoef0(i)*wcoef1(i)*windmag/rho
      wcoef2(i) = wcoef0(i)*wcoef2(i)*windmag/rho
    END DO

  END SELECT

! Free the memory used for the temporary arrays and close datafile.
  DEALLOCATE(tmp1,tmp2,tmp3,wcoef0)
  CLOSE (f)

END SUBROUTINE windforcing
