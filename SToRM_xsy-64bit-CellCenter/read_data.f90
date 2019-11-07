SUBROUTINE read_data(filename,title)
  USE io
  USE options
  USE dep_vars
  !USE IFPORT ! For portability functions IARGC and GETARG using IVF.
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This subroutine reads the name of the main input file and chooses the      !
!  appropriate data loader based on the extension:                            !
!                                                                             !
!     Extension    Data loader                                                !
!    ---------------------------------------------------------------------    !
!      .run        Multiple file ASCII (STORM's original data format).        !
!      .cgn        CGNS standard.                                             !
!      .cgns       CGNS standard.                                             !
!                                                                             !
!  Upon exit from the subroutine, the dummy argument variables 'filename'     !
!  and 'title' contain, respectively, the name of the main input file         !
!  (without the extension) and a title to be written to the output file.      !
!                                                                             !
!  Francisco Simoes, December 2003                                            !
!  Last updated (mm-dd-yyyy): 11-04-2016 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  CHARACTER (LEN=*), INTENT(OUT) :: filename,title

! Local variables.
  INTEGER :: ierror,tmp_unit
  CHARACTER (LEN=3) :: fext  ! For file extension query.
  CHARACTER (LEN=80) :: buffer,out_file
  LOGICAL :: f,fnexist
  INTEGER, EXTERNAL :: probefileh
  CHARACTER (LEN=3), EXTERNAL :: get_extension
  LOGICAL, EXTERNAL :: get_iounit

! Get input file from command line.
  IF (IARGC() == 0) THEN
10  WRITE (*,'("  Enter the INPUT file name: ",$)') ! The $ edit descriptor is
    READ (*,'(A)') filename                         ! not Fortran 90 standard.
    filename = ADJUSTL(filename)
    INQUIRE (FILE=filename,EXIST=fnexist)
    IF (.NOT.fnexist) THEN
      WRITE (*,'(2A/)') 'ERROR: file not found, ', TRIM(filename)
      WRITE (*,'("Press <Enter> to continue")'); READ (*,'(A)') buffer
      GO TO 10
    END IF
  ELSE
    CALL GETARG(1,filename)
    filename = ADJUSTL(filename)
  END IF

  fext = get_extension(filename)

! If the file is in CGNS format, use the appropriate data loader.
  IF (fext == "CGN") THEN
    cgns = .TRUE.
    !CALL byebye('CGNS not running in this version of STORM.')
    CALL read_cgns(filename,title)

    ! Write data structures.  Note that the solutions field may be appended to
    ! this file later.
    out_file = TRIM(filename) // ".output"
    IF (.NOT. get_iounit(o_unit)) CALL byebye('ERROR in get_iounit(o_unit).')
    OPEN (o_unit,FILE=out_file,IOSTAT=ierror)
    IF (ierror /= 0) THEN
      WRITE (*,'(" ERROR: unable to open output file ",A)') TRIM(out_file)
      CALL byebye('SToRM stopped.')
    END IF

    TITLE = "CGNS I/O"
    CALL header(o_unit,filename,title)
    !CALL strip_ext(filename)

! The data is in multiple ASCII data files.  Use the Multiple File ASCII (mfa)
! data loader.
  ELSE IF (fext == "RUN") THEN
    cgns = .FALSE.
    CALL read_mfa(filename,title)
    CALL strip_ext(filename)
  ELSE

! Unrecognized extension: issue error message and stop the program.
    WRITE (*,'(A)') 'ERROR: unrecognized file extension.'
    CALL byebye('SToRM stopped.')
  END IF

! Read datafile with culvert information.
  IF (culvert) CALL ioculvert

! Open and write header in probe output file.
  IF (nprobe > 0) p_unit = probefileh()

! Delete CGNS-Tecplot hot start file, if it exists.
  INQUIRE (FILE="Tecplot_sol_file.dat",EXIST=f)
  IF (f) THEN
    IF (.NOT. get_iounit(tmp_unit)) &
      CALL byebye('ERROR in get_iounit(tmp_unit).')
      OPEN (tmp_unit,FILE="Tecplot_sol_file.dat",IOSTAT=ierror)
      IF (ierror /= 0) THEN
        WRITE (*,'(" ERROR: unable to open CGNS-Tecplot hot-start file.")')
        CALL byebye('SToRM stopped.')
      END IF
    CLOSE (tmp_unit,STATUS='DELETE')
  END IF

END SUBROUTINE read_data
