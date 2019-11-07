SUBROUTINE read_mfa(filename,title)
  USE parameters
  USE geometry
  USE io
  USE options
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This is the data loader when the data is in multiple file ASCII (mfa)      !
!  format, which is the original format developed for and used in STORM.      !
!  Upon exit from this subroutine, the dummy argument variable 'title'        !
!  contains a title to be written to the output file.  The dummy argument     !
!  variable 'filename' contains the name of the main input file (with         !
!  extension '.run') that contains all the names of the files where the data  !
!  is stored.                                                                 !
!                                                                             !
!  Francisco Simoes, June 2006                                                !
!  Last updated (mm-dd-yyyy): 03-22-2011 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  CHARACTER (LEN=*), INTENT(IN) :: filename
  CHARACTER (LEN=*), INTENT(OUT) :: title

! Local variables.
  INTEGER :: bcs_unit,iofile,ierror,opsfile,tecpfile
  CHARACTER (LEN=80) :: bcs_file,grid_file,options_file,out_file, &
                        sol_file
  LOGICAL :: fnexist,optionals,rconn,stats
  LOGICAL, EXTERNAL :: fnames,get_iounit

!-----------------------------------------------------------------------------!
!                                                                             !
!  Read the input file name.  This file must contain the names of the files   !
!  with the data needed by the SToRM's run of interest  One file name per     !
!  line, with a maximum of 80 characters per file name.  Lines starting with  !
!  the character # are comment lines and are discarded; blank lines are       !
!  skipped.  The following files are needed (by order in which they must      !
!  appear in the file:                                                        !
!                                                                             !
!  -  geometry data file name (in Tecplot format)                             !
!  -  boundary conditions file name                                           !
!  -  output file name where some of the data goes to (mainly from            !
!     subroutine 'header').                                                   !
!  -  main solutions file (in Tecplot format)                                 !
!  -  options file name (may be absent)                                       !
!                                                                             !
!  See the comments below about the appropriate data format required for      !
!  each of the input files.                                                   !
!                                                                             !
!-----------------------------------------------------------------------------!

  IF (.NOT. get_iounit(iofile)) THEN
    WRITE (*,'(" ERROR: unable to open file ",A)') TRIM(filename)
    CALL byebye('SToRM stopped.')
  END IF
  OPEN (iofile,FILE=filename,STATUS='OLD',IOSTAT=ierror)
  IF (ierror /= 0) THEN
    WRITE (*,'(" ERROR: unable to open file ",A)') TRIM(filename)
    CALL byebye('SToRM stopped.')
  END IF

! Read file names.
! Read grid file name.
  stats = fnames(iofile,grid_file)
! Read the boundary conditions file name.
  stats = fnames(iofile,bcs_file)
! Read output file name (where headers and such go to).
  stats = fnames(iofile,out_file)
! Read solution output file (in Tecplot format).
  stats = fnames(iofile,sol_file)
! Read file name with the run's options.
  optionals = fnames(iofile,options_file)

  CLOSE (iofile)

! Read-in the data in Tecplot format.
  INQUIRE (file=grid_file,EXIST=fnexist)
  IF (.NOT. fnexist) THEN
    WRITE (*,'(2A)') 'ERROR: file not found, ', TRIM(grid_file)
    CALL byebye('SToRM stopped.')
  END IF
  IF (.NOT. get_iounit(tecpfile)) THEN
    WRITE (*,'(" ERROR: unable to open grid file ",A)') TRIM(grid_file)
    CALL byebye('SToRM stopped.')
  END IF
  OPEN (tecpfile,FILE=grid_file,STATUS='OLD',IOSTAT=ierror)
  IF (ierror /= 0) THEN
    WRITE (*,'(" ERROR: unable to open grid file ",A)') TRIM(grid_file)
    CALL byebye('SToRM stopped.')
  END IF

  CALL tecplot(tecpfile)
  CLOSE (tecpfile)

! Set-up default on all options.
  CALL initialize_options

! Read options from external file.
  IF (optionals) THEN
    IF (get_iounit(opsfile)) THEN
      OPEN (opsfile,FILE=options_file,STATUS='OLD',IOSTAT=ierror)
      IF (ierror /= 0) THEN
        WRITE (*,'(/" ERROR: unable to open options file ",A)') &
          TRIM(options_file)
        WRITE (*,'(" SToRM run will proceed with all default values."/)')
      ELSE
        CALL read_options(opsfile)
        CLOSE (opsfile)
      END IF
    ELSE
      WRITE (*,'(/" ERROR: unable to open options file ",A)') &
        TRIM(options_file)
      WRITE (*,'(" SToRM run will proceed with all default values."/)')
    END IF
  END IF

  ! Open boundary conditions file.
  INQUIRE (file=bcs_file,EXIST=fnexist)
  IF (.NOT. fnexist) THEN
    WRITE (*,'(2A)') 'ERROR: file not found, ', TRIM(bcs_file)
    CALL byebye('SToRM stopped.')
  END IF
  IF (get_iounit(bcs_unit)) THEN
    OPEN (bcs_unit,FILE=bcs_file,STATUS='OLD',IOSTAT=ierror)
    IF (ierror /= 0) THEN
      WRITE (*,'(" ERROR: unable to open boundary conditions file ",A)') &
        TRIM(bcs_file)
      CALL byebye('SToRM stopped.')
    END IF
  ELSE
    WRITE (*,'(" ERROR: unable to open boundary conditions file ",A)') &
      TRIM(bcs_file)
    CALL byebye('SToRM stopped.')
  END IF

  IF (read_conn) INQUIRE(FILE=f_read_conn,EXIST=rconn)

  IF (read_conn .AND. rconn) THEN  ! Read all the mesh connectivity data from a
                                   ! binary external file.---------------------

    ! Read-in mesh connectivity arrays.
    CALL mesh_in(f_read_conn)

    ! Read-in the boundary conditions file, in case discharge and/or stage
    ! values changed.  Note that stage and discharge values may change, and
    ! even the type of these boundary conditions, but not the nodes at which
    ! they are enforced.
    CALL ioflow(bcs_unit)

  ELSE  ! Build mesh connectivity arrays.--------------------------------------

    IF (read_conn) THEN
      WRITE (*,'(/2A)') 'File not found: ',TRIM(f_read_conn)
      WRITE (*,'(A/)') 'STORM building the mesh connectivity arrays...'
    END IF

    ! Prepare edge information and data structure.
    CALL find_edges

    ! Compute triangle areas.
    CALL t_areas(grid,n_elems,nodes,n_pts,edges,n_edges)

    ! Read boundary conditions.
    CALL bcs(bcs_unit)

    ! Set-up node data structure.
    CALL node_db

    ! Set-up element connectivity data structure.
    CALL elemnt_db

  END IF

  CLOSE (bcs_unit)

  ! Write data structures.  Note that the solutions field may be appended to
  ! this file later.
  IF (.NOT. get_iounit(o_unit)) THEN
    PRINT *,'ERROR in get_iounit(o_unit).'
    CALL byebye('Program SToRM stopped.')
  END IF
  OPEN (o_unit,FILE=out_file,IOSTAT=ierror)
  IF (ierror /= 0) THEN
    WRITE (*,'(" ERROR: unable to open output file ",A)') TRIM(out_file)
    CALL byebye('SToRM stopped.')
  END IF

  CALL header(o_unit,filename,title)

! Write grid file.
  IF (write_conn) CALL mesh_out(f_write_conn)

! Write solution field in Tecplot format.
  IF (.NOT. get_iounit(s_unit)) CALL byebye('ERROR in get_iounit(s_unit).')
  OPEN (s_unit,FILE=sol_file,IOSTAT=ierror)
  IF (ierror /= 0) THEN
    WRITE (*,'(" ERROR: unable to open output file ",A)') TRIM(sol_file)
    CALL byebye('SToRM stopped.')
  END IF

END SUBROUTINE read_mfa
