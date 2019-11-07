SUBROUTINE spit_outv(iter,fname,t)
  USE geometry
  USE dep_vars
  USE options
  USE constants
  USE io
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Export the solution to a Tecplot-formatted file.  The output filename is   !
!  constructed from the iteration number 'iter' and the root file name        !
!  'fname' using subroutine ifname (with additional information from array    !
!  iout() ).  The title for the Tecplot file is created from variable 't'.    !
!  This version spits out the vertex values, just replace the usual           !
!  'spit-out' SUBROUTINE with this one.                                       !
!                                                                             !
!  Francisco Simoes, April 2008                                               !
!  Last updated (mm-dd-yyyy): 03-22-2011 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  INTEGER, INTENT(IN) :: iter  ! Iteration number.
  CHARACTER (LEN=*) :: fname,t

! Local variables.
  INTEGER :: f_unit,ierror,slength
  INTEGER, EXTERNAL :: ccstring
  CHARACTER (LEN=80) :: buffer,cellctr,fn,nelem_char,npts_char

  LOGICAL, EXTERNAL :: get_iounit

! Generate output file name.
  CALL ifname(iter,max_iout,fname,fn)

! Write solution.  First open file...
  IF (.NOT. get_iounit(f_unit)) CALL byebye('ERROR in get_iounit().')
  OPEN (f_unit,FILE=fn,IOSTAT=ierror)
  IF (ierror /= 0) THEN
    WRITE (*,'(" ERROR: unable to open output file ",A)') TRIM(fn)
    CALL byebye('SToRM stopped.')
  END IF

! Set-up CELLCENTERD string for Tecplot file header.
  slength = ccstring(datav,11,cellctr)

! Write Tecplot header.
  WRITE (buffer,*) iter
  buffer = ADJUSTL(buffer)
  WRITE (f_unit,'(A)') 'TITLE = "' // TRIM(t) // ' [ITER ' // TRIM(buffer) // &
    &']"'
  WRITE (f_unit,'(A)') 'VARIABLES = "X", "Y", "Z", "H", "U", "V", "CD", &
    &"PHI_H", "PHI_U", "PHI_V", "ERROR"'
  WRITE (npts_char,*) n_pts
  npts_char = ADJUSTL(npts_char)
  WRITE (nelem_char,*) n_elems
  nelem_char = ADJUSTL(nelem_char)
  WRITE (f_unit,'(A)') "ZONE N=" // TRIM(npts_char) //", E=" // &
    TRIM(nelem_char) // ", DATAPACKING=BLOCK, ZONETYPE=FETRIANGLE" // &
    ", VARLOCATION=([7,8,9,10,11]=CELLCENTERED)"

  CALL tecplot_writev(f_unit)

  CLOSE (f_unit)

END SUBROUTINE spit_outv
