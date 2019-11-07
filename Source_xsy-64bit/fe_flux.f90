SUBROUTINE fe_flux
  USE parameters
  USE geometry
  USE dep_vars
  USE constants
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This subroutine computes the finite element flux at the midpoints of all   !
!  cell edges. The equations are formulated very similarly to those in eqs.   !
!  (1-2) of Brufau and Garcia-Navarro (2000), for example.  Upon exit, the    !
!  array flux() contains the line integral of the convective flux through     !
!  the edge, for every edge of the computational domain, i.e., the            !
!  right-hand side of eq. (7) in Sleigh et al. (1998).  Calculations are      !
!  made strictly in the finite element sense, in which the dependent          !
!  variables are h, hu, and hv.                                               !
!                                                                             !
!  Francisco Simoes, November 2005                                            !
!  Last updated (mm-dd-yyyy): 11-19-2004 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables:
  INTEGER :: i,j,p1,p2
  REAL(KIND=mp) :: h1,h2,hface,u1,u2,uface,v1,v2,vface
  TYPE(edge) :: side

  DO j = 1,flow_edges
    i = flowedg(j)
    side = edges(i)
    p1 = side%p(1)
    p2 = side%p(2)

    ! Convective flux, continuity equation.  Later, some of these products can
    ! be put in arrays to save CPU time at the expense of added memory
    ! requirements.
    h1 = h(p1) ; u1 = h1*u(p1) ; v1 = h1*v(p1)
    h2 = h(p2) ; u2 = h2*u(p2) ; v2 = h2*v(p2)
    hface = (h1 + h2)*half
    uface = (u1 + u2)*half
    vface = (v1 + v2)*half
    flux(i,1) = uface*side%normal(1) + vface*side%normal(2)

    ! Convective flux, x-momentum.
    flux(i,2) = (uface*uface/hface + half*g*hface*hface)*side%normal(1) + &
                uface*vface/hface*side%normal(2)

    ! Convective flux, y-momentum.
    flux(i,3) = uface*vface/hface*side%normal(1) + &
                (vface*vface/hface + half*g*hface*hface)*side%normal(2)
  END DO

  DO j = 1,wall_edges  ! Solid walls.
    i = walledg(j)
    side = edges(i)
    p1 = side%p(1)
    p2 = side%p(2)
    hface = (h(p1) + h(p2))*half
    hface = half*g*hface*hface

    flux(i,1) = zero
    flux(i,2) = hface*side%normal(1)
    flux(i,3) = hface*side%normal(2)
  END DO

END SUBROUTINE fe_flux
