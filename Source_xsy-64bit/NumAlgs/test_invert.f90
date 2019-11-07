PROGRAM test_invert
  USE parameters
  IMPLICIT NONE

  INTEGER :: i
  REAL(KIND=mp) :: a(3,3)
  LOGICAL :: t
  LOGICAL, EXTERNAL :: invert

  a(1,1) = 3.4_mp;   a(1,2) = 5.9_mp;     a(1,3) = 89.1_mp
  a(2,1) = 909.0_mp; a(2,2) = 0_mp;       a(2,3) = 1024.0_mp
  a(3,1) = -43.9_mp; a(3,2) = -1005.1_mp; a(3,3) = 89.1_mp

  t = invert(a)
  PRINT *,t
  DO i = 1,3
    PRINT *,a(i,1),a(i,2),a(i,3)
  END DO

END PROGRAM test_invert
