LOGICAL FUNCTION circle(p1,p2,p3,c,r)
  USE parameters
  USE constants
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Function circle computes the center and radius of a circle defined by      !
!  three points following the method of SSchmitt.pdf (http://home.att.net/    !
!  ~srscmitt/circle3pts.html, accessed 5/6/2005).                             !
!                                                                             !
!  INPUT:                                                                     !
!    p1,p2,p3   the three points defining the circle.                         !
!                                                                             !
!  OUTPUT:                                                                    !
!    c          coordinates of the center of the circle;                      !
!    r          radius of the circle;                                         !
!    circle     .TRUE. if the computations yielded a circle;                  !
!               .FALSE. if the three points do not define a circle, e.g.,     !
!               colinear or coincident points.                                !
!                                                                             !
!  F. Simoes, May 2005                                                        !
!  Last updated (mm-dd-yyyy): 05-06-2005                                      !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables:
  TYPE(point), INTENT(IN) :: p1,p2,p3
  TYPE(point), INTENT(OUT) :: c
  REAL(KIND=mp), INTENT(OUT) :: r

! Local variables:
  REAL(KIND=mp) :: m11,m12,m13,m14
  REAL(KIND=mp) :: a(3,3)
  REAL(KIND=mp), EXTERNAL :: det3by3
  LOGICAL, EXTERNAL :: equals

! Find minor 11.
  a(1,1) = p1%x; a(2,1) = p2%x; a(3,1) = p3%x
  a(1,2) = p1%y; a(2,2) = p2%y; a(3,2) = p3%y
  a(1,3) = one;  a(2,3) = one;  a(3,3) = one
  m11 = det3by3(a)

! Find minor 12.
  a(1,1) = p1%x*p1%x + p1%y*p1%y
  a(2,1) = p2%x*p2%x + p2%y*p2%y
  a(3,1) = p3%x*p3%x + p3%y*p3%y
  a(1,2) = p1%y; a(2,2) = p2%y; a(3,2) = p3%y
  a(1,3) = one;  a(2,3) = one;  a(3,3) = one
  m12 = det3by3(a)

! Find minor 13.  a(i,1) is the same as in minor 12.
  a(1,2) = p1%x; a(2,2) = p2%x; a(3,2) = p3%x
  a(1,3) = one;  a(2,3) = one;  a(3,3) = one
  m13 = det3by3(a)

! Find minor 14.  a(i,1) and a(i,2) are the same as in minor 13.
  a(1,3) = p1%y; a(2,3) = p2%y; a(3,3) = p3%y
  m14 = det3by3(a)

  IF (equals(m11,zero)) THEN
    ! The three points do not define a circle.
    circle = .FALSE.
  ELSE
    c%x = half*m12/m11
    c%y = -half*m13/m11
    r = SQRT(c%x*c%x + c%y*c%y + m14/m11)
    circle = .TRUE.
  END IF

END FUNCTION circle
