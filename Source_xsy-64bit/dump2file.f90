SUBROUTINE dump2file(array,n,fname)
  USE parameters
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Writes an array to a file.                                                 !
!                                                                             !
!  INPUT:                                                                     !
!    array    array to be printed to file (REAL (KIND=mp);                    !
!    n        dimension of the array (single dimension);                      !
!    fname    name of output file.                                            !
!                                                                             !
!  Francisco Simoes, 07-10-2012                                               !
!  Last updated (mm-dd-yyyy): 07-10-2012 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy argument.
  INTEGER, INTENT(IN) :: n
  REAL (KIND=mp), INTENT(IN) :: array(n)
  !INTEGER, INTENT(IN) :: array(n)
  CHARACTER(*), INTENT(IN) :: fname

! Local variables.
  INTEGER :: i,ierror,ounit
  LOGICAL, EXTERNAL :: get_iounit

  IF (get_iounit(ounit)) THEN
    OPEN (ounit,FILE=TRIM(fname),STATUS='UNKNOWN',IOSTAT=ierror)
    DO i = 1,n
      WRITE (ounit,*) i,array(i)
    END DO
    CLOSE (ounit)
  END IF
END SUBROUTINE dump2file
