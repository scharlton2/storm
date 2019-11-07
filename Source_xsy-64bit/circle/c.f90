PROGRAM c
  USE parameters
  IMPLICIT NONE

  TYPE(point) :: ctr,p1,p2,p3
  REAL(KIND=mp) :: r
  LOGICAL, EXTERNAL :: circle

  p1%x =  1.0; p1%y =  7.0
  p2%x =  8.0; p2%y =  6.0
  p3%x =  7.0; p3%y = -2.0

  IF (circle(p1,p2,p3,ctr,r)) THEN
    PRINT *,"RADIUS =",r
    PRINT *,"CENTER =",ctr
  ELSE
    PRINT *,"Program failed."
  END IF

END PROGRAM c
