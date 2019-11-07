SUBROUTINE adjustZb(flag,e,gradzb,gradzb2)
  USE parameters
  USE constants
  USE geometry
  USE dep_vars
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Adjustment of the bed slope for partially wet triangles.  This is done by  !
!  lowering each vertex that has negative water depth to the bed elevation    !
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
!  Francisco Simoes, April 2008                                               !
!  Last updated (mm-dd-yyyy): 04-16-2008 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  INTEGER, INTENT(IN) :: e
  LOGICAL, INTENT(IN) :: flag
  TYPE(vector), INTENT(OUT) :: gradzb,gradzb2

! Local variables:
  INTEGER :: i,j
  REAL (KIND=mp) :: a,b(3),c(3),hp(3),x(3),y(3),zp(3)
  LOGICAL :: adjust

  adjust = .FALSE.

  DO i = 1,3
    j = grid(e)%vertex(i)
    ! Water depth at vertex.
    hp(i) = h(e) + delh(e)%x*(nodes(j)%x - grid(e)%xc) + &
                   delh(e)%y*(nodes(j)%y - grid(e)%yc)
    IF (hp(i) < zero) adjust = .TRUE.
  END DO

! 'adjust' triggers the adjustment.
  IF (adjust) THEN
    DO i = 1,3
      j = grid(e)%vertex(i)
      hp(i) = MIN(zero,hp(i))
      zp(i) = zvtx(j) + hp(i)  ! Subtract negative water depth (lower node).
      x(i) = nodes(j)%x
      y(i) = nodes(j)%y
    END DO
    b(1) = y(2) - y(3) ; c(1) = x(3) - x(2)
    b(2) = y(3) - y(1) ; c(2) = x(1) - x(3)
    b(3) = y(1) - y(2) ; c(3) = x(2) - x(1)
    a = half/grid(e)%area
    gradzb%x = a*(b(1)*zp(1) + b(2)*zp(2) + b(3)*zp(3))
    gradzb%y = a*(c(1)*zp(1) + c(2)*zp(2) + c(3)*zp(3))

  ELSE  ! If no adjustment...
    gradzb%x = delz(e)%x
    gradzb%y = delz(e)%y
  END IF

  IF (flag) THEN
    gradzb2%x = 2.0_mp*z(e)*gradzb%x
    gradzb2%y = 2.0_mp*z(e)*gradzb%y
  END IF

END SUBROUTINE adjustZb
