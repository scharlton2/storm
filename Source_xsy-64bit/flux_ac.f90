SUBROUTINE flux_ac(celledge,hlefty,uleft,vleft,hrighty,uright,vright,flux)
  USE parameters
  USE constants
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Evaluates the inviscid numerical fluxes through edge 'celledge' using      !
!  Roe's flux function.  The cell states at the left and right of the edge    !
!  must be given.  The computation follows closely the presentation in        !
!  section 3.1 of Anastasiou and Chan (1997).                                 !
!                                                                             !
!  INPUT:                                                                     !
!    celledge              a cell edge defined by the SToRM data structure;   !
!    hleft,uleft,vlefty    state at the left-side cell sharing edge           !
!                          'celledge';                                        !
!    hright,uright,vrighty state at the right-side cell sharing edge          !
!                          'celledge';                                        !
!                                                                             !
!  OUTPUT:                                                                    !
!    flux                  Roe's flux. flux(1) corresponds to the continuity  !
!                          equation; flux(2) to the u-momentum equation; and  !
!                          flux(3) to the v-momentum equation.                !
!                                                                             !
!  NOTE: it uses the entropy fix of Alcrudo et al. (1992) via parameter       !
!  'delta', which must be in the interval [0.1,1.0].  To deactivate, simply   !
!  set 'delta' to zero.                                                       !
!                                                                             !
!  Francisco Simoes, March 2007                                               !
!  Last updated (mm-dd-yyyy): 04-02-2008 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  TYPE(edge), INTENT(IN) :: celledge
  REAL (KIND=mp), INTENT(IN) :: hlefty,hrighty,uleft,uright,vleft,vright
  REAL (KIND=mp), INTENT(OUT) :: flux(3)

! Local variables.
  INTEGER :: i
  REAL (KIND=mp) :: cavg,denom,havg,hleft,hleftsqr,hright,hrightsqr,nx,ny, &
                    uavg,vavg
  REAL (KIND=mp) :: a(3,3),aq(3),eigv(3),fileft(3),firight(3),gileft(3), &
                    giright(3),lev(3,3),qdiff(3),rev(3,3),temp1(3,3),temp2(3,3)
  REAL (KIND=mp) :: delta = 0.10_mp  ! For the entropy fix.
  REAL (KIND=mp), PARAMETER :: small = 1.0e-9  ! For wet-dry fronts.

  nx = celledge%normal(1)
  ny = celledge%normal(2)

  hleft = hlefty
  hright = hrighty

! Compute the left viscous flux.
  fileft(1) = hleft*uleft
  fileft(2) = hleft*uleft*uleft + half*g*hleft*hleft
  fileft(3) = hleft*uleft*vleft
  gileft(1) = hleft*vleft
  gileft(2) = hleft*uleft*vleft
  gileft(3) = hleft*vleft*vleft + half*g*hleft*hleft

! Compute the right viscous flux.
  firight(1) = hright*uright
  firight(2) = hright*uright*uright + half*g*hright*hright
  firight(3) = hright*uright*vright
  giright(1) = hright*vright
  giright(2) = hright*uright*vright
  giright(3) = hright*vright*vright + half*g*hright*hright

! Now compute Roe's correction based on the flux Jacobian using average state
! variables following Brufau and Garcia-Navarro (2000).
  hleft = MAX(hleft,small)
  hleftsqr=SQRT(hleft)
  hright = MAX(hright,small)
  hrightsqr = SQRT(hright)
  denom = hrightsqr + hleftsqr
  uavg = (uright*hrightsqr + uleft*hleftsqr)/denom
  vavg = (vright*hrightsqr + vleft*hleftsqr)/denom
  cavg = SQRT(half*g*(hleft + hright))

! Plain average, as in Anastasiou and Chan (1997).
  !havg = half*(hleft + hright)
  !uavg = half*(hleft*uleft + hright*uright)/havg
  !vavg = half*(hleft*vleft + hright*vright)/havg
  !cavg = SQRT(g*havg)

  qdiff(1) = hright - hleft
  qdiff(2) = hright*uright - hleft*uleft
  qdiff(3) = hright*vright - hleft*vleft

! Right eigenvector matrix, eq. (8) of Anastasiou & Chan (1997):
  rev(1,1) = zero; rev(1,2) = one; rev(1,3) = one
  rev(2,1) = ny; rev(2,2) = uavg - cavg*nx; rev(2,3) = uavg + cavg*nx
  rev(3,1) = -nx; rev(3,2) = vavg - cavg*ny; rev(3,3) = vavg + cavg*ny

! Left eigenvector matrix, eq. (9) of Anastasiou & Chan (1997):
  lev(1,1) = -(uavg*ny - vavg*nx)
  lev(1,2) = ny
  lev(1,3) = -nx
  lev(2,1) = half*(uavg*nx + vavg*ny)/cavg + half
  lev(2,2) = -half*nx/cavg
  lev(2,3) = -half*ny/cavg
  lev(3,1) = -half*(uavg*nx + vavg*ny)/cavg + half
  lev(3,2) = half*nx/cavg
  lev(3,3) = half*ny/cavg

! Diagonal matrix with the absolute value of the eigenvalues:
  eigv(1) = MAX(ABS(uavg*nx + vavg*ny),delta)
  eigv(2) = MAX(ABS(uavg*nx + vavg*ny - cavg),delta)
  eigv(3) = MAX(ABS(uavg*nx + vavg*ny + cavg),delta)
  a = zero
  a(1,1) = eigv(1)
  a(2,2) = eigv(2)
  a(3,3) = eigv(3)

! Compute right-most term of eq. (5) of Anastasiou & Chan (1997):
  temp1 = MATMUL(rev,a)
  temp2 = MATMUL(temp1,lev)
  aq = MATMUL(temp2,qdiff)

! Finally, add all the terms to obtain the total inviscid flux through the
! edge, as presented in eq. (5) of Anastasiou & Chan (1997).
  DO i = 1,3
    flux(i) = half*(firight(i)*nx + giright(i)*ny + fileft(i)*nx + &
              gileft(i)*ny - aq(i))
  END DO

END SUBROUTINE flux_ac
