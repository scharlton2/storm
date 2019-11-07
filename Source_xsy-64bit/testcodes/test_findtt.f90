MODULE parameters
  IMPLICIT NONE
  SAVE

! Machine precision:
  INTEGER, PARAMETER :: mp = KIND(1.0D0)  ! = KIND(1.0) for single precision.
END MODULE parameters

PROGRAM test_findtt
  USE parameters
  IMPLICIT NONE
  CHARACTER (LEN=40) :: times
  REAL (KIND=mp) :: dt = 0.1_mp
  INTEGER, EXTERNAL :: findtt

! Get input file from command line.
  IF (NARGS() == 1) THEN
    WRITE (*,'(" Enter data: ",$)')
    READ (*,'(A)') times
  ELSE
    CALL GETARG(1,times)
  END IF

  times = ADJUSTL(times)
  PRINT *,dt
  PRINT *,findtt(times,dt)

END PROGRAM test_findtt

INTEGER FUNCTION findtt(strg,dt)
  USE parameters
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  A function to convert time into a specific number of time steps.  Usage:   !
!                                                                             !
!  INPUT:                                                                     !
!    strg     string with the value to be decoded.  Trailing blanks are ok,   !
!             but all leading blanks must be eliminated before the call to    !
!             findtt (use ADJUSTL);                                           !
!    dt       time step size, in seconds.                                     !
!                                                                             !
!  OUTPUT:                                                                    !
!    findtt   number of time steps for the combination of the given           !
!             function arguments.  Returns -1 if there is an error.           !
!                                                                             !
!  Usage of "strg": it must be a number (a decimal point is ok) or a number   !
!  followed by a letter.  If it is a number, the value is converted to        !
!  integer (truncated) to be used directly as number of time steps.  The      !
!  number can be followed by one (and only one) letter as follows: s or S     !
!  for seconds; m or M for minutes; h or H for hours; and d or D for days.    !
!  The value is converted to seconds and divided by the time step, and the    !
!  result is truncated to integer to yield a number of time steps.            !
!                                                                             !
!  F. Simoes, August 2014                                                     !
!  Last updated (mm-dd-yyyy): 08-21-2014                                      !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables.
  REAL (KIND=mp), INTENT(IN) :: dt
  CHARACTER (*), INTENT(IN) :: strg
! Local variables.
  REAL (KIND=mp) :: s
  INTEGER :: l
  LOGICAL, EXTERNAL :: findttcheck

  l = LEN_TRIM(strg)
  findtt = -1
  SELECT CASE(strg(l:l))  ! Seconds.
    CASE ('S','s')
      IF (.NOT.findttcheck(strg(1:l-1))) RETURN  ! Error in string.
      READ (strg(1:l-1),*) s
      findtt = INT(s/dt)
    CASE ('M','m')
      IF (.NOT.findttcheck(strg(1:l-1))) RETURN  ! Error in string.
      READ (strg(1:l-1),*) s
      s = s*60.0_mp  ! Convert minutes to seconds.
      findtt = INT(s/dt)
    CASE ('H','h')
      IF (.NOT.findttcheck(strg(1:l-1))) RETURN  ! Error in string.
      READ (strg(1:l-1),*) s
      s = s*3600.0_mp  ! Convert hours to seconds.
      findtt = INT(s/dt)
    CASE ('D','d')
      IF (.NOT.findttcheck(strg(1:l-1))) RETURN  ! Error in string.
      READ (strg(1:l-1),*) s
      s = s*86400.0_mp  ! Convert days to seconds.
      findtt = INT(s/dt)
    CASE DEFAULT  ! Time steps.
      IF (.NOT.findttcheck(strg(1:l))) RETURN
      READ (strg(1:l),*) s
      findtt = INT(s)
  END SELECT

END FUNCTION findtt

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
