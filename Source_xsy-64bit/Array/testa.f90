SUBROUTINE expand(n)
  IMPLICIT NONE

  INTEGER :: i,n
  INTEGER :: a(n)

  DO i = 1,n
    a(i) = i*10
  END DO

  DO i = 1,n
    PRINT *,a(i)
  END DO

END SUBROUTINE expand

PROGRAM testa
  IMPLICIT NONE

  INTEGER :: k

  PRINT *,'Input k:'
  READ (*,*) k
  CALL expand(k)

END PROGRAM testa
