INTEGER FUNCTION ccstring(a,n,string)
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This function prepares a string to be used to write a Tecplot-formatted    !
!  ASCII datafile with the STORM variables.  The string is meant to go where  !
!  the VARLOCATION and CELLCENTERED Tecplot qualifiers usually go, and        !
!  between straight parenthesis ([]).                                         !
!                                                                             !
!  INPUT:                                                                     !
!    a         the array containing the location of the variables to be       !
!              written to the Tecplot file (a(i) = .TRUE. if variable i is    !
!              located at the vertices; a(i) = .FALSE. if it is located at    !
!              the cell centers);                                             !
!    n         dimension of array a().                                        !
!                                                                             !
!  OUTPUT:                                                                    !
!    string    a string with the proper values for insertion in the Tecplot   !
!              file header as                                                 !
!                "VARLOCATION=(["//TRIM(string)//"]=CELLCENTERED)";           !
!    ccstring  length of string (0 if string is empty.                        !
!                                                                             !
!  Francisco Simoes, April 2004                                               !
!  Last updated (mm-dd-yyyy): 03-22-2006 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments:
  INTEGER, INTENT(IN) :: n
  LOGICAL, INTENT(IN), DIMENSION(n) :: a
  CHARACTER (LEN=*), INTENT(OUT) :: string

! Local variables:
  INTEGER :: i,slength
  CHARACTER (LEN=2) :: c

  string = ''
! First do the single digit loop.  This adds the string '03,' to string when
! i = 3, etc.
  DO i = 1,9
    WRITE (c,'(I1)') i
    IF (.NOT.a(i)) string = TRIM(string) // '0' // TRIM(ADJUSTL(c)) &
      // ','
  END DO
! Now do the two digit loop.  This adds the string '13,' to cellctr when
! i = 13, etc.
  DO i = 10,n
    WRITE (c,'(I2)') i
    IF (.NOT.a(i)) string = TRIM(string) // c // ','
  END DO
  ! Remove the last comma.
  slength = LEN_TRIM(string)
  IF (slength > 1) string(slength:slength) = ''

  ccstring = MAX(0,slength-1)

END FUNCTION ccstring
