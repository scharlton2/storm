MODULE parameters
  IMPLICIT NONE
  SAVE

! Machine precision:
  INTEGER, PARAMETER :: mp = KIND(1.0D0)  ! = KIND(1.0) for single precision.

END MODULE parameters

FUNCTION f(x,y,dfdx,dfdy,df2dx,df2dy)
  USE parameters
  IMPLICIT NONE

! Computes the bed elevation values (zb) for a function of x and y and its
! respective slopes.  Also computes the slopes of zb**2.

! Dummy arguments.
  REAL (KIND=mp), INTENT(IN) :: x,y
  REAL (KIND=mp), INTENT(OUT) :: dfdx,dfdy  ! Derivatives of f.
  REAL (KIND=mp), INTENT(OUT) :: df2dx,df2dy  ! Derivatives of f**2.
  REAL (KIND=mp) :: f

! Local variables.
  REAL (KIND=mp) :: a,b,c

! A simple plane.
  a = 1.0_mp
  b = 0.20_mp
  c = 0.50_mp
  f = a + b*x + c*y
  dfdx = b
  dfdy = c
  !fsquared = a**2 + 2*a*b*x + 2*a*c*y + 2*b*c*x*y + b*b*x*x + c*c*y*y
  df2dx = 2.0_mp*(a*b + b*c*y + b*b*x)
  df2dy = 2.0_mp*(a*c + b*c*x + c*c*y)

END FUNCTION f


!-----------------------------------------------------------------------------!
!                                                                             !
!  This program is used to test the discretization of several forms of the    !
!  source term due to bed slope.  In particular, the interest resides in      !
!  determining the accuracy of the methods in the case of partially dry       !
!  triangles.                                                                 !
!                                                                             !
!  Francisco Simoes, February 2008                                            !
!  Last updated (mm-dd-yyyy): 02-22-2008 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!


PROGRAM bed_slope
  USE parameters
  IMPLICIT NONE

! Local variables.
  REAL (KIND=mp), DIMENSION(3) :: eta,etaside,h,hside,nx,ny,x,y,zb,zbside,zbsq
  REAL (KIND=mp) :: area,etaavg,h0,havg,havg2,sc1x,sc1y,sc2x,sc2y,sc3x,sc3y, &
                    sc4x,sc4y,sc5x,sc5y,sc6x,sc6y,sc7x,sc7y,xc,yc,zbavg, &
                    zbsqx,zbsqy,zbx,zby,zc
  INTEGER :: i
  REAL (KIND=mp), EXTERNAL :: f

! Coordinates of vertices of triangle, defined counterclockwise.
  x(1) = 1.0_mp; x(2) = 6.0_mp; x(3) = 6.0_mp
  y(1) = 1.0_mp; y(2) = 1.0_mp; y(3) = 6.0_mp
  x(1) = 1.0_mp; x(2) = 6.0_mp; x(3) = 1.0_mp
  y(1) = 1.0_mp; y(2) = 1.0_mp; y(3) = 6.0_mp
  xc = (x(1) + x(2) + x(3))/3.0_mp
  yc = (y(1) + y(2) + y(3))/3.0_mp

! Area of triangle.
  area = 0.50_mp*ABS((x(1) - x(2))*(y(1) + y(2)) + (x(2) - x(3))*(y(2) + y(3)) &
         + (x(3) - x(1))*(y(3) + y(1)))

! Define bed elevation at the triangle's verices.
  !zb(1) = 10.0_mp; zb(2) = 11.0_mp; zb(3) = 11.0_mp
  zb(1) = f(x(1),y(1),zbx,zby,zbsqx,zbsqy)
  zb(2) = f(x(2),y(2),zbx,zby,zbsqx,zbsqy)
  zb(3) = f(x(3),y(3),zbx,zby,zbsqx,zbsqy)

! Define stage at the triangle's vertices.
  eta(1) = 2.80_mp; eta(2) = 2.910_mp; eta(3) = 3.0_mp

! Define water depth at the triangle's vertices.
  h = eta - zb

! Average quantities.
  zbavg = (zb(1) + zb(2) + zb(3))/3.0_mp      ! Bed elevation.
  etaavg = (eta(1) + eta(2) + eta(3))/3.0_mp  ! Stage.
  havg = etaavg - zbavg                       ! Water depth.
  havg2 = (h(1) + h(2) + h(3))/3.0_mp

! Bed gradient, finite element style.
  !zbx = ((y(2) - y(3))*zb(1) + (y(3) - y(1))*zb(2) + (y(1) - y(2))*zb(3)) &
  !      /2.0_mp/area
  !zby = ((x(3) - x(2))*zb(1) + (x(1) - x(3))*zb(2) + (x(2) - x(1))*zb(3)) &
  !      /2.0_mp/area

! Classic cell-centered approximation.
  sc1x = havg*zbx*area
  sc1y = havg*zby*area

! Using eq. (15) of Bradford and Sanders (2002).
  !zbsq = zb*zb
  !zbsqx = ((y(2) - y(3))*zbsq(1) + (y(3) - y(1))*zbsq(2) + &
  !        (y(1) - y(2))*zbsq(3))/2.0_mp/area
  !zbsqy = ((x(3) - x(2))*zbsq(1) + (x(1) - x(3))*zbsq(2) + &
  !        (x(2) - x(1))*zbsq(3))/2.0_mp/area
  zc = f(xc,yc,zbx,zby,zbsqx,zbsqy)

  sc2x = (etaavg*zbx - zbsqx/2.0_mp)*area
  sc2y = (etaavg*zby - zbsqy/2.0_mp)*area

  sc3x = (etaavg*zbx - zbavg*zbx)*area  ! Perhaps a better way to do it...
  sc3y = (etaavg*zby - zbavg*zby)*area

! Using the divergence form, i.e., eq. (27) of Valiani and Begnudelli (2006).
  nx(1) = y(2) - y(1); ny(1) = x(1) - x(2)
  nx(2) = y(3) - y(2); ny(2) = x(2) - x(3)
  nx(3) = y(1) - y(3); ny(3) = x(3) - x(1)
  xc = (x(1) + x(2))/2.0_mp
  yc = (y(1) + y(2))/2.0_mp
  zbside(1) = f(xc,yc,zbx,zby,zbsqx,zbsqy)
  xc = (x(2) + x(3))/2.0_mp
  yc = (y(2) + y(3))/2.0_mp
  zbside(2) = f(xc,yc,zbx,zby,zbsqx,zbsqy)
  xc = (x(3) + x(1))/2.0_mp
  yc = (y(3) + y(1))/2.0_mp
  zbside(3) = f(xc,yc,zbx,zby,zbsqx,zbsqy)
  sc4x = 0.0_mp; sc4y = 0.0_mp
  DO i = 1,3
    sc4x = sc4x - 0.50_mp*((etaavg - zbside(i))**2)*nx(i)
    sc4y = sc4y - 0.50_mp*((etaavg - zbside(i))**2)*ny(i)
  END DO

! Using Gauss' theorem.
  sc5x = 0.0_mp; sc5y = 0.0_mp
  DO i = 1,3
    sc5x = sc5x + zbside(i)*nx(i)
    sc5y = sc5y + zbside(i)*ny(i)
  END DO
  sc5x = sc5x/area*havg*area
  sc5y = sc5y/area*havg*area

! Using Simoes' expansion.
  etaside(1) = (eta(1) + eta(2))/2.0_mp
  etaside(2) = (eta(2) + eta(3))/2.0_mp
  etaside(3) = (eta(3) + eta(1))/2.0_mp
  sc6x = 0.0_mp; sc6y = 0.0_mp
  DO i = 1,3
    sc6x = sc6x + etaside(i)*zbside(i)*nx(i) - zc*(etaside(i) + zbside(i))*nx(i)
    sc6y = sc6y + etaside(i)*zbside(i)*ny(i) - zc*(etaside(i) + zbside(i))*ny(i)
  END DO

! Using the 2D generalization of Jin (2001).
  hside(1) = (h(1) + h(2))/2.0_mp
  hside(2) = (h(2) + h(3))/2.0_mp
  hside(3) = (h(3) + h(1))/2.0_mp
  sc7x = 0.0_mp; sc7y = 0.0_mp; h0 = 0.0_mp
  DO i = 1,3
    sc7x = sc7x + zbside(i)*nx(i)
    sc7y = sc7y + zbside(i)*ny(i)
    h0 = h0 + hside(i)
  END DO
  sc7x = sc7x*h0/3.0_mp
  sc7y = sc7y*h0/3.0_mp

! Print out the results.
  PRINT *
  !PRINT '(X,A,3(ES12.5,5X))','zb =',zb
  PRINT *,'zb =',zb
  PRINT *,'zb @ c =',zbavg,zc
  PRINT *,'area =',area
  PRINT *,'eta avg =',etaavg
  PRINT *,'eta =',eta
  PRINT *,'h avg =',havg,havg2
  PRINT *,'h =',eta - zb
  PRINT *,'h side =',hside
  PRINT *,'dzb/dx =',zbx
  PRINT *,'dzb/dy =',zby
  PRINT *,'nx =',nx
  PRINT *,'ny =',ny
  PRINT *

  PRINT *,'Std Sc  =',sc1x,sc1y!/area/havg
  PRINT *,'B&S Sc  =',sc2x,sc2y!/area/havg
  PRINT *,'B&S Sc  =',sc3x,sc3y!/area/havg
  PRINT *,'V&B Sc  =',sc4x,sc4y!/area/havg
  PRINT *,'Gauss   =',sc5x,sc5y!/area/havg
  PRINT *,'Simoes  =',sc6x,sc6y
  PRINT *,'Jin     =',sc7x,sc7y

  !PRINT *,0.0_mp
  !PRINT *,0.50_mp
  !PRINT *,1.0_mp/2.0_mp

END PROGRAM bed_slope
