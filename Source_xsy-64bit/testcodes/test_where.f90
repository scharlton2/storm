PROGRAM test_where
  IMPLICIT NONE

  INTEGER, PARAMETER :: n = 10
  INTEGER :: i
  REAL  :: u(n),v(n),u_mag(n)
  LOGICAL :: drycell(n)

  drycell = .FALSE.
  drycell(2) = .TRUE. ; drycell(3) = .TRUE. ; drycell(7) = .TRUE.

  FORALL (i = 1:n)
    u(i) = REAL(i)
    v(i) = REAL(i)
  END FORALL

  u_mag = 0.0
  WHERE (.NOT. drycell) u_mag = SQRT(u*u + v*v)

  DO i = 1,n
    PRINT *,i,drycell(i),u(i),v(i),u_mag(i),SQRT(u(i)*u(i) + v(i)*v(i))
  END DO

END PROGRAM test_where
