MODULE mem
  IMPLICIT NONE
  SAVE

  INTEGER :: counter = 1

END MODULE mem


SUBROUTINE check(i)
  USE mem
  IMPLICIT NONE

  INTEGER :: i

  counter = counter + i
  PRINT *,counter

END SUBROUTINE check


PROGRAM testc
  IMPLICIT NONE

! Note that counter is initialized to 1 in the MODULE block...

  CALL check(0)
  CALL Check(5)
  CALL check(3)
  CALL check(0)

END PROGRAM testc
