SUBROUTINE flux_bgnvc(celledge,hleft,uleft,vleft,hright,uright,vright,flux)
  USE parameters
  USE constants
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Evaluates the inviscid numerical fluxes through edge 'celledge' using      !
!  Roe's flux function.  The cell states at the left and right of the edge    !
!  must be given.  The computation follows closely the presentation in        !
!  page 1051 of Brufau et al. (2004).                                         !
!                                                                             !
!  INPUT:                                                                     !
!    celledge              a cell edge defined by the SToRM data structure;   !
!    hleft,uleft,vleft     state at the left-side cell sharing edge           !
!                          'celledge';                                        !
!    hright,uright,vright  state at the right-side cell sharing edge          !
!                          'celledge';                                        !
!                                                                             !
!  OUTPUT:                                                                    !
!    flux                  Roe's flux. flux(1) corresponds to the continuity  !
!                          equation; flux(2) to the u-momentum equation; and  !
!                          flux(3) to the v-momentum equation.                !
!                                                                             !
!  Francisco Simoes, August 2007                                              !
!  Last updated (mm-dd-yyyy): 12-11-2007 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  TYPE(edge), INTENT(IN) :: celledge
  REAL (KIND=mp), INTENT(IN) :: hleft,hright,uleft,uright,vleft,vright
  REAL (KIND=mp), INTENT(OUT) :: flux(3)

! Local variables.
  INTEGER :: i
  REAL (KIND=mp) :: cavg,denom,hleftsqr,hrightsqr,nx,ny,uavg,vavg
  REAL (KIND=mp) :: a(3,3),aq(3),fileft(3),firight(3),gileft(3),giright(3), &
                    lev(3,3),qdiff(3),rev(3,3),temp1(3,3),temp2(3,3)

  nx = celledge%normal(1)
  ny = celledge%normal(2)

! Compute the left viscous flux.
  fileft(1) = hleft*uleft*nx
  fileft(2) = (hleft*uleft*uleft + half*g*hleft*hleft)*nx
  fileft(3) = hleft*uleft*vleft*nx
  gileft(1) = hleft*vleft*ny
  gileft(2) = hleft*uleft*vleft*ny
  gileft(3) = (hleft*vleft*vleft + half*g*hleft*hleft)*ny

! Compute the right viscous flux.
  firight(1) = hright*uright*nx
  firight(2) = (hright*uright*uright + half*g*hright*hright)*nx
  firight(3) = hright*uright*vright*nx
  giright(1) = hright*vright*ny
  giright(2) = hright*uright*vright*ny
  giright(3) = (hright*vright*vright + half*g*hright*hright)*ny

! Now compute Roe's correction based on the flux Jacobian using average state
! variables following Brufau and Garcia-Navarro (2000).
  hleftsqr=SQRT(hleft)
  hrightsqr = SQRT(hright)
  denom = hrightsqr + hleftsqr
  uavg = (uright*hrightsqr + uleft*hleftsqr)/denom
  vavg = (vright*hrightsqr + vleft*hleftsqr)/denom
  cavg = SQRT(half*g*(hleft + hright))

  qdiff(1) = hright - hleft
  qdiff(2) = hright*uright - hleft*uleft
  qdiff(3) = hright*vright - hleft*vleft

! Right eigenvector matrix P:
  rev(1,1) = one; rev(1,2) = zero; rev(1,3) = one
  rev(2,1) = uavg + cavg*nx; rev(2,2) = -cavg*ny; rev(2,3) = uavg - cavg*nx
  rev(3,1) = vavg + cavg*ny; rev(3,2) = cavg*nx; rev(3,3) = vavg - cavg*ny

! Left eigenvector matrix P^-1:
  lev(1,1) = uavg*nx + vavg*ny + cavg
  lev(1,2) = nx
  lev(1,3) = ny
  lev(2,1) = 2.0_mp*(uavg*ny - vavg*nx)
  lev(2,2) = -2.0_mp*ny
  lev(2,3) = 2.0_mp*nx
  lev(3,1) = uavg*nx + vavg*ny + cavg
  lev(3,2) = -nx
  lev(3,3) = vavg - ny
  lev = half*lev/cavg

! Matrix with the absolute value of the eigenvalues:
  a(1,1) = ABS(uavg*nx + vavg*ny + cavg); a(1,2) = zero; a(1,3) = zero
  a(2,1) = zero; a(2,2) = ABS(uavg*nx + vavg*ny); a(2,3) = zero
  a(3,1) = zero; a(3,2) = zero; a(3,3) = ABS(uavg*nx + vavg*ny - cavg)

! Compute right-most term of eq. (5) of Anastasiou & Chan (1997):
  temp1 = MATMUL(rev,a)
  temp2 = MATMUL(temp1,lev)
  aq = MATMUL(temp2,qdiff)

! Finally, add all the terms to obtain the total inviscid flux through the
! edge, as presented in eq. (5) of Anastasiou & Chan (1997).
  DO i = 1,3
    flux(i) = half*(firight(i) + giright(i) + fileft(i) + gileft(i) - aq(i))
  END DO

END SUBROUTINE flux_bgnvc
