PROGRAM funcnoarg
  IMPLICIT NONE

! See how a function without argument works...

  INTEGER :: i
  INTEGER, EXTERNAL :: omm

  i = omm()
  PRINT *,i

END PROGRAM funcnoarg

INTEGER FUNCTION omm
  IMPLICIT NONE

  omm = 5
END FUNCTION omm
