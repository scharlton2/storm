LOGICAL FUNCTION findttcheck(strg)
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This function looks for unauthorized characters in the argument string     !
!  (strg).  The only authorized characters in strg are numbers (0-9) and dot  !
!  (.).  But there must be no more than one dot in strg.  This function is    !
!  used in "findtt".                                                          !
!                                                                             !
!  F. Simoes, August 2014                                                     !
!  Last updated (mm-dd-yyyy): 08-21-2014                                      !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables.
  CHARACTER (*), INTENT(IN) :: strg
! Local variables.
  INTEGER :: c,i,np

  findttcheck = .TRUE.
  np = 0
  DO i = 1,LEN_TRIM(strg)
    c = IACHAR(strg(i:i))
    IF ((c < 48 .OR. c > 57) .AND. c /= 46) THEN
      findttcheck = .FALSE.
      RETURN
    END IF
    IF (c == 46) np = np + 1
  END DO
  IF (np > 1) findttcheck = .FALSE.

END FUNCTION findttcheck
