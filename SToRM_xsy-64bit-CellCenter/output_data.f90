SUBROUTINE output_data(itermax,ttime,resmax,filename,title)
  USE options
  USE geometry
  USE dep_vars
  USE constants
  USE io
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Main output subroutine for solution variables.                             !
!                                                                             !
!  Francisco Simoes, April 2004                                               !
!  Last updated (mm-dd-yyyy): 06-02-2017 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables:
  INTEGER, INTENT(IN) :: itermax  ! Number of iterations achieved.
  REAL(KIND=mp), INTENT(IN) :: resmax  ! Maximum value of the residuals.
  REAL(KIND=mp), INTENT(IN) :: ttime !rmcd mod
  CHARACTER (LEN=*) :: title, filename !rmcd mod 'filename'

! Local variables:
  INTEGER :: i,ierror,n,slength,value(8)
  INTEGER, EXTERNAL :: ccstring
  ! Note: the size of cellctr must be >= n*3, where n is the dimension of
  ! datav().
  CHARACTER :: cellctr*80
  CHARACTER :: hour*10,today*8,zone*5  ! Used in intrinsic time subroutine.
  CHARACTER (LEN=3) :: month(12)
  CHARACTER (LEN=20) :: nelem_char,npts_char
  LOGICAL, EXTERNAL :: get_iounit

  month = (/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct', &
            'Nov','Dec'/)

! If using CGNS, write the solution in the CGNS file.
  IF (cgns) THEN
    CALL write_cgns(ttime)
    CALL write_cgns_cell(ttime)
    IF (.NOT. tecplot_print) RETURN
    ! Note: this RETURN statement prevents output to screen when running from
    ! the MD_SWMS interface.  The drawback is that no ASCII output is done
    ! for the Tecplot file or for the main output file.  This could probably
    ! be redone a little bit better, perhaps similarly to routine spit_out.
  END IF

!-----------------------------------------------------------------------------!
!                                                                             !
! Write the solution in Tecplot format.                                       !
!                                                                             !
!-----------------------------------------------------------------------------!

! Set-up CELLCENTERD string for Tecplot file header.
  slength = ccstring(datav,11,cellctr)

! This code opens the Tecplot output file if tecplot_print = .TRUE. when using
! CGNS I/O.  This implementation only works if the output of data is done only
! upon exit of SToRM.
  IF (cgns) THEN
    IF (.NOT. get_iounit(s_unit)) CALL byebye('ERROR in get_iounit(s_unit).')
    OPEN (s_unit,FILE="Tecplot_sol_file.dat",IOSTAT=ierror)
    IF (ierror /= 0) THEN
      WRITE (*,'(" ERROR: unable to open CGNS-Tecplot output file.")')
      CALL byebye('SToRM stopped.')
    END IF
  END IF

! First, write the Tecplot header.
  WRITE (s_unit,'(A)') 'TITLE = "'//TRIM(title)//'"'
  WRITE (s_unit,'(A)') 'VARIABLES = "X", "Y", "Z", "H", "U", "V", "CD", &
    &"PHI_H", "PHI_U", "PHI_V", "ERROR"'
  WRITE (npts_char,*) n_pts
  npts_char = ADJUSTL(npts_char)
  WRITE (nelem_char,*) n_elems
  nelem_char = ADJUSTL(nelem_char)
  IF (slength > 0) THEN
    WRITE (s_unit,'(A)') "ZONE N=" // TRIM(npts_char) //", E=" // &
      TRIM(nelem_char) // ", DATAPACKING=BLOCK, ZONETYPE=FETRIANGLE" // &
      ", VARLOCATION=([" // TRIM(cellctr) // "]=CELLCENTERED)"
  ELSE
    WRITE (s_unit,'(A)') "ZONE N=" // TRIM(npts_char) //", E=" // &
      TRIM(nelem_char) // ", DATAPACKING=BLOCK, ZONETYPE=FETRIANGLE"
  END IF

  CALL tecplot_write(s_unit)

  CLOSE (s_unit)

!-----------------------------------------------------------------------------!
!                                                                             !
! Complete and close the main output file.                                    !
!                                                                             !
!-----------------------------------------------------------------------------!

! Program's header file.
  IF (output_head) THEN
    WRITE (o_unit,'(//80("=")//)')
    WRITE (o_unit,'("Number of iterations:",I9)') itermax
    WRITE (o_unit,'("Maximum residual:",ES13.5)') resmax

    WRITE (o_unit,'(/24X,A)')"Computed solution"
    WRITE (o_unit,'(8X,"N",7X,"U",13X,"V",13X,"H",13X,"Z",13X,"CD")')
    WRITE (o_unit,'(80("-"))')

    SELECT CASE (opt_solver)
    CASE (1)  ! RDS solver.
      n = n_pts

    CASE (2)  ! Finite volume solver using triangles for control volumes (FVT).
      n = n_elems

    END SELECT

    DO i = 1,n
      WRITE (o_unit,'(I9,5(2X,ES12.5))') i,u(i),v(i),h(i),z(i),cd(i)
    END DO

    CALL DATE_AND_TIME(today,hour,zone,value)
    WRITE (o_unit,'(//80("=")//)')
    WRITE (o_unit,'(A,I2,1X,A,1X,I4)')"Run ended at "//hour(1:2)//":"// &
      & hour(3:4)//":"//hour(5:6)//" on ",value(3),month(value(2)),value(1)
  END IF

  CLOSE (o_unit)

END SUBROUTINE output_data
