PROGRAM inflow
  IMPLICIT NONE

! Test an idea about implementing inflow boundary conditions in SToRM.

! Machine precision:
  INTEGER, PARAMETER :: mp = KIND(1.0D0)  ! = KIND(1.0) for single precision.
  REAL(mp), PARAMETER :: zero = 0.0, two = 2.0, vsmall = 0.000001
  REAL(mp) :: at,ht,q,qp,qpp
  REAL(mp), ALLOCATABLE, DIMENSION(:) :: a,h,l,u,up,upp,v,vp
  INTEGER :: i,n

! Read the data and allocate memory.
  OPEN (10,FILE="inflow.dat")
  READ (10,*) n,q
  ALLOCATE(a(n),h(n),l(n),u(n),up(n),upp(n),v(n),vp(n))
  DO i = 1,n
    READ (10,*) l(i),h(i)
  END DO

! Initialize variables.
  v = zero
  vp = zero
  u = zero
  up = zero
  upp = zero
  ht = zero
  at = zero
  DO i = 1,n-1
    a(i) = l(i)*(h(i) + h(i+1))/two
    ht = ht + (h(i) + h(i+1))/two
    at = at + a(i)
  END DO
  DO i = 1,n-1
    IF (a(i) > vsmall) THEN
      u(i) = (h(i) + h(i+1))*q/ht/at/two
    ELSE
      u(i) = zero
    END IF
  END DO

! First pass to compute v(i).
  v(1) = zero
  v(n) = zero
  DO i = 2,n-1
    IF (h(i) > vsmall) v(i) = (u(i-1) + u(i))/two
  END DO

! Second pass to compute the new u(i).
  DO i = 1,n-1
    up(i) = (v(i) + v(i+1))/two
  END DO

! Compute the new discharge.
  qp = zero
  DO i = 1,n-1
    qp = qp + up(i)*a(i)
  END DO

! Adjust v(i) accordingly.
  vp(1) = zero
  vp(n) = zero
  DO i = 2,n-1
    IF (h(i) > vsmall) vp(i) = v(i)*q/qp
  END DO

! Compute the adjusted discharge.
  DO i = 1,n-1
    upp(i) = (vp(i) + vp(i+1))/two
  END DO
  qpp = zero
  DO i = 1,n-1
    qpp = qpp + upp(i)*a(i)
  END DO

! Print results.
  PRINT *,q,qp,qpp
  DO i = 1,n
    PRINT *,i,vp(i),upp(i)
  END DO

END PROGRAM inflow
