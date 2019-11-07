LOGICAL FUNCTION get_iounit(iounit)
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This function returns an unused Fortran UNIT specifier to use in an OPEN   !
!  statement.  The argument "iounit" contains the UNIT number in the          !
!  interval [10,109].  Usually, several UNIT numbers below 10 are reserved    !
!  (e.g., units 0, 5, and 6 are associated with preconnected devices by the   !
!  Compaq Fortran compiler), therefore those are skipped.  The upper limit    !
!  of the interval can go to 2**31-1, but there is no need for it in most     !
!  practical applications.  get_iounit will be .FALSE. if no free Fortran     !
!  unit can be found in the specified interval, otherwise it will return      !
!  .TRUE.                                                                     !
!                                                                             !
!  Francisco Simoes, 09-23-2006                                               !
!  Last updated (mm-dd-yyyy): 07-03-2012 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!


! Dummy argument.
  INTEGER, INTENT(OUT) :: iounit

! Local variables.
  INTEGER :: i,ierr
  LOGICAL :: is_open

  get_iounit = .FALSE.
  iounit = -1

  DO i = 10,109

    INQUIRE(UNIT=i,OPENED=is_open,IOSTAT=ierr)

    IF (ierr == 0 .AND. .NOT.is_open) THEN
      iounit = i
      get_iounit = .TRUE.
      RETURN
    END IF

  END DO

END FUNCTION get_iounit
