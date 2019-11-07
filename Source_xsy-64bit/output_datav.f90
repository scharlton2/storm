SUBROUTINE output_dataV(itermax,resmax,title)
  USE options
  USE geometry
  USE dep_vars
  USE constants
  USE io
  !USE WriteCGNSData
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Main output subroutine for solution variables.  Use instead of subroutine  !
!  output_data() to write vertex data directly to the final Tecplot file.     !
!  Use only for debugging purposes.                                           !
!                                                                             !
!  Francisco Simoes, October 2007                                             !
!  Last updated (mm-dd-yyyy): 04-18-2008 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables:
  INTEGER, INTENT(IN) :: itermax  ! Number of iterations achieved.
  REAL(KIND=mp), INTENT(IN) :: resmax  ! Maximum value of the residuals.
  CHARACTER (LEN=*) :: title

! Local variables:
  INTEGER :: i,j,slength,value(8)
  INTEGER, EXTERNAL :: ccstring
  ! Note: the size of cellctr must be >= n*3, where n is the dimension of
  ! datav().
  CHARACTER :: cellctr*80
  CHARACTER :: hour*10,today*8,zone*5  ! Used in intrinsic time subroutine.
  CHARACTER (LEN=3) :: month(12)
  CHARACTER (LEN=20) :: nelem_char,npts_char

  month = (/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct', &
            'Nov','Dec'/)

! If using CGNS, write the solution in the CGNS file.
  !IF (cgns) CALL write_TimeStep_CGNS(TRIM(filename), 0, 0.)

!-----------------------------------------------------------------------------!
!                                                                             !
! Write the solution in Tecplot format.                                       !
!                                                                             !
!-----------------------------------------------------------------------------!

! Set-up CELLCENTERD string for Tecplot file header.
  slength = ccstring(datav,11,cellctr)

! First, write the Tecplot header.
  WRITE (s_unit,'(A)') 'TITLE = "'//TRIM(title)//'"'
  WRITE (s_unit,'(A)') 'VARIABLES = "X", "Y", "Z", "H", "U", "V", "CD", &
    &"PHI_H", "PHI_U", "PHI_V", "ERROR"'
  WRITE (npts_char,*) n_pts
  npts_char = ADJUSTL(npts_char)
  WRITE (nelem_char,*) n_elems
  nelem_char = ADJUSTL(nelem_char)
  WRITE (s_unit,'(A)') "ZONE N=" // TRIM(npts_char) //", E=" // &
    TRIM(nelem_char) // ", DATAPACKING=BLOCK, ZONETYPE=FETRIANGLE" // &
    ", VARLOCATION=([7,8,9,10,11]=CELLCENTERED)"

  DO i = 1,n_pts
    WRITE (s_unit,'(ES24.15)') nodes(i)%x
  END DO
  DO i = 1,n_pts
    WRITE (s_unit,'(ES24.15)') nodes(i)%y
  END DO
  DO i = 1,n_pts
    WRITE (s_unit,'(ES24.15)') zvtx(i)
  END DO
  DO i = 1,n_pts
    WRITE (s_unit,'(ES24.15)') hvtx(i)
  END DO
  DO i = 1,n_pts
    WRITE (s_unit,'(ES24.15)') uvtx(i)
  END DO
  DO i = 1,n_pts
    WRITE (s_unit,'(ES24.15)') vvtx(i)
  END DO
  DO i = 1,n_elems
    WRITE (s_unit,'(ES24.15)') cd(i)
  END DO
  DO j = 1,3
    DO i = 1,n_elems
      WRITE (s_unit,'(ES24.15)') phi(i,j)
    END DO
  END DO
  DO i = 1,n_elems
    WRITE (s_unit,'(ES24.15)') e_est(i)
  END DO

! Finally, write the element connectivity table.
  DO i = 1,n_elems
    WRITE (s_unit,*) grid(i)%vertex(1),grid(i)%vertex(2),grid(i)%vertex(3)
  END DO

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
    DO i = 1,n_pts
      WRITE (o_unit,'(I9,5(2X,ES12.5))') i,u(i),v(i),h(i),z(i),cd(i)
    END DO

    CALL DATE_AND_TIME(today,hour,zone,value)
    WRITE (o_unit,'(//80("=")//)')
    WRITE (o_unit,'(A,I2,1X,A,1X,I4)')"Run ended at "//hour(1:2)//":"// &
      & hour(3:4)//":"//hour(5:6)//" on ",value(3),month(value(2)),value(1)
  END IF

  CLOSE (o_unit)

END SUBROUTINE output_dataV
