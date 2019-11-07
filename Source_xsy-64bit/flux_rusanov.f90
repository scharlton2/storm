SUBROUTINE flux_rusanov(celledge,hleft,uleft,vleft,sqrtghleft,hright,uright, &
                        vright,sqrtghright,flux)
  USE parameters
  USE constants
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Evaluates the inviscid numerical fluxes through edge 'celledge' using      !
!  the Rusanov flux function.  The cell states at the left and right of the   !
!  edge must be given.  More details about the Rusanov flux can be found in   !
!  page 184 of Toro (2001), for example.                                      !
!                                                                             !
!  INPUT:                                                                     !
!    celledge               a cell edge defined by the SToRM data structure.  !
!                           Its normal must point from the left-side cell to  !
!                           the right-side cell.                              !
!    hleft,uleft,vleft      state at the left-side cell sharing edge          !
!                           'celledge';                                       !
!    hright,uright,vright   state at the right-side cell sharing edge         !
!                           'celledge';                                       !
!    sqrtghleft,sqrtghright wave celerity, SQRT(g*h), of the left and right   !
!                           cells, respectively;                              !
!                                                                             !
!  OUTPUT:                                                                    !
!    flux                   Rusanov flux. flux(1) corresponds to the          !
!                           continuity equation; flux(2) to the u-momentum    !
!                           equation; and flux(3) to the v-momentum           !
!                           equation.                                         !
!                                                                             !
!  Francisco Simoes, December 2007                                            !
!  Last updated (mm-dd-yyyy): 12-11-2007 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  TYPE(edge), INTENT(IN) :: celledge
  REAL (KIND=mp), INTENT(IN) :: hleft,hright,sqrtghleft,sqrtghright,uleft, &
                                uright,vleft,vright
  REAL (KIND=mp), INTENT(OUT) :: flux(3)

! Local variables.
  INTEGER :: i
  REAL (KIND=mp) :: nx,ny,splus,ul,ur
  REAL (KIND=mp) :: fileft(3),firight(3),gileft(3),giright(3),qdiff(3)

  nx = celledge%normal(1)
  ny = celledge%normal(2)

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

! Wave speed estimates, eq. [10.43] of Toro (2001).
  ul = uleft*nx + vleft*ny
  ur = uright*nx + vright*ny
  splus = MAX(ABS(ul) + sqrtghleft,ABS(ur) + sqrtghright)

  qdiff(1) = hright - hleft
  qdiff(2) = hright*uright - hleft*uleft
  qdiff(3) = hright*vright - hleft*vleft

! Finally, add all the terms to obtain the total inviscid flux through the
! edge, as presented in eq. (5) of Anastasiou & Chan (1997).
  DO i = 1,3
    flux(i) = half*(firight(i)*nx + giright(i)*ny + fileft(i)*nx + &
              gileft(i)*ny - splus*qdiff(i))
  END DO

END SUBROUTINE flux_rusanov
