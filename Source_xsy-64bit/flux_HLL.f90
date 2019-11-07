SUBROUTINE flux_HLL(celledge,hleft,uleft,vleft,hright,uright,vright,h0,flux)
  USE parameters
  USE constants
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Evaluates the inviscid numerical fluxes through edge 'celledge' using      !
!  the HLL Riemann solver.  The cell states at the left and right of the      !
!  edge must be given.  The computation follows closely the presentation in   !
!  page 679 of Yoon and Kang (2004).                                          !
!                                                                             !
!  INPUT:                                                                     !
!    celledge              a cell edge defined by the SToRM data structure.   !
!                          It must be pointing from left to right;            !
!    hleft,uleft,vleft     state at the left-side cell sharing edge           !
!                          'celledge';                                        !
!    hright,uright,vright  state at the right-side cell sharing edge          !
!                          'celledge';                                        !
!    h0                    water depth threshold for dry cell.                !
!                                                                             !
!  OUTPUT:                                                                    !
!    flux                  Roe's flux. flux(1) corresponds to the continuity  !
!                          equation; flux(2) to the u-momentum equation; and  !
!                          flux(3) to the v-momentum equation.                !
!                                                                             !
!  Francisco Simoes, April 2008                                               !
!  Last updated (mm-dd-yyyy): 04-19-2008 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  TYPE(edge), INTENT(IN) :: celledge
  REAL (KIND=mp), INTENT(IN) :: h0,hleft,hright,uleft,uright,vleft,vright
  REAL (KIND=mp), INTENT(OUT) :: flux(3)

! Local variables.
  REAL (KIND=mp) :: aux1,aux2,cleft,cright,cstar,nx,ny,sleft,sright,ustar
  REAL (KIND=mp) :: fleft(3),fright(3),gleft(3),gright(3)

  nx = celledge%normal(1)
  ny = celledge%normal(2)
  cleft = SQRT(g*hleft)
  cright = SQRT(g*hright)
  ustar = half*((uleft + uright)*nx + (vleft + vright)*ny) + cleft - cright
  cstar = half*(cleft + cright) + ((uleft - uright)*nx + &
          (vleft - vright)*ny)/4.0_mp

! Compute sleft, eq. (7) of Yoon and Kang (2004).
  IF (hleft < h0) THEN  ! Left side dry.
    sleft = uright*nx + vright*ny - 2.0_mp*cright
  ELSE IF (hright < h0) THEN  ! Right side dry.
    sleft = uleft*nx + vleft*ny - cleft
  ELSE  ! Both sides wet.
    aux1 = uleft*nx + vleft*ny - cleft
    aux2 = ustar - cstar
    sleft = MIN(aux1,aux2)
  END IF

! Compute sleft, eq. (8) of Yoon and Kang (2004).
  IF (hleft < h0) THEN  ! Left side dry.
    sright = uright*nx + vright*ny + cright
  ELSE IF (hright < h0) THEN  ! Right side dry.
    sright = uleft*nx + vleft*ny + 2.0_mp*cleft
  ELSE  ! Both sides wet.
    aux1 = uright*nx + vright*ny - cright
    aux2 = ustar + cstar
    sright = MAX(aux1,aux2)
  END IF

! Compute the left inviscid flux.
  fleft(1) = hleft*uleft*nx
  fleft(2) = (hleft*uleft*uleft + half*g*hleft*hleft)*nx
  fleft(3) = hleft*uleft*vleft*nx
  gleft(1) = hleft*vleft*ny
  gleft(2) = hleft*uleft*vleft*ny
  gleft(3) = (hleft*vleft*vleft + half*g*hleft*hleft)*ny

! Compute the right inviscid flux.
  fright(1) = hright*uright*nx
  fright(2) = (hright*uright*uright + half*g*hright*hright)*nx
  fright(3) = hright*uright*vright*nx
  gright(1) = hright*vright*ny
  gright(2) = hright*uright*vright*ny
  gright(3) = (hright*vright*vright + half*g*hright*hright)*ny

! Now compute the flux, eq. (6) of Yoon and Kang (2004).
  IF (sleft > zero .AND. sright < zero) THEN  ! Temporary debugging code...
    PRINT *,"Problem with HLL flux: debug."
    CALL byebye("Program SToRM stopped.")
  END IF
  aux1 = one
  IF (sleft > zero) THEN
    flux(1) = fleft(1) + gleft(1)
    flux(2) = fleft(2) + gleft(2)
    flux(3) = fleft(3) + gleft(3)
  ELSE IF (sright < zero) THEN
    flux(1) = fright(1) + gright(1)
    flux(2) = fright(2) + gright(2)
    flux(3) = fright(3) + gright(3)
  ELSE  ! sleft <= zero and sright >= zero.
    !aux1 = MAX(sright - sleft,one)
    aux1 = sright - sleft
    flux(1) = sright*(fleft(1) + gleft(1)) - sleft*(fright(1) + gright(1)) + &
           sleft*sright*(hright - hleft)/aux1
    flux(2) = sright*(fleft(2) + gleft(2)) - sleft*(fright(2) + gright(2)) + &
           sleft*sright*(hright*uright - hleft*uleft)/aux1
    flux(2) = sright*(fleft(3) + gleft(3)) - sleft*(fright(3) + gright(3)) + &
           sleft*sright*(hright*vright - hleft*vleft)/aux1
  END IF

END SUBROUTINE flux_HLL
