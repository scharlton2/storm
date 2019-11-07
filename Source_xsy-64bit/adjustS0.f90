SUBROUTINE adjustS0(flag,e,gradzb,gradzb2)
  USE parameters
  USE constants
  USE geometry
  USE dep_vars
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Adjustment of the bed slope for partially wet triangles.  This is done by  !
!  lowering each side that has negative water depth to the bed elevation      !
!  corresponding to the zero depth level.  The gradients are then computed    !
!  using the finite element technique using shape functions.                  !
!                                                                             !
!  INPUT:                                                                     !
!    flag     .TRUE. if the gradient of z^2 is also computed, .FALSE. if      !
!             not;                                                            !
!    e        element index.                                                  !
!                                                                             !
!  OUTPUT:                                                                    !
!    gradzb   the gradient of the bed elevation;                              !
!    gradzb2  the gradient of the bed elevation squared.  This variable is    !
!             optional, and is used when the technique of Bradford and        !
!             Sanders (2002) for the bed source terms is used.                !
!                                                                             !
!  Francisco Simoes, October 2007                                             !
!  Last updated (mm-dd-yyyy): 10-29-2007 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  INTEGER, INTENT(IN) :: e
  LOGICAL, INTENT(IN) :: flag
  TYPE(vector), INTENT(OUT) :: gradzb,gradzb2

! Local variables:
  INTEGER :: i,p1,p2
  REAL (KIND=mp) :: a,b(3),c(3),hmid,x(3),y(3),zmid(3),z2mid(3)

  gradzb%x = zero
  gradzb%y = zero

  DO i = 1,3
    p1 = edges(ABS(grid(e)%edge(i)))%p(1)
    p2 = edges(ABS(grid(e)%edge(i)))%p(2)
    hmid = h(e) + delh(e)%x*rc(e,i)%x + delh(e)%y*rc(e,i)%y
    zmid(i) = half*(zvtx(p1) + zvtx(p2))
    IF (hmid < zero) zmid(i) = zmid(i) + hmid
    x(i) = half*(nodes(p1)%x + nodes(p2)%x)
    y(i) = half*(nodes(p1)%y + nodes(p2)%y)
  END DO
  b(1) = y(2) - y(3) ; c(1) = x(3) - x(2)
  b(2) = y(3) - y(1) ; c(2) = x(1) - x(3)
  b(3) = y(1) - y(2) ; c(3) = x(2) - x(1)
  a = half*grid(e)%area
  gradzb%x = (b(1)*zmid(1) + b(2)*zmid(2) + b(3)*zmid(3))/a
  gradzb%y = (c(1)*zmid(1) + c(2)*zmid(2) + c(3)*zmid(3))/a

  IF (flag) THEN
    z2mid = zmid*zmid
    gradzb2%x = (b(1)*z2mid(1) + b(2)*z2mid(2) + b(3)*z2mid(3))/a
    gradzb2%y = (c(1)*z2mid(1) + c(2)*z2mid(2) + c(3)*z2mid(3))/a
  END IF

END SUBROUTINE adjustS0
