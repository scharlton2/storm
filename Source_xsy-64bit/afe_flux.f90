SUBROUTINE afe_flux
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
!  right-hand side of eq. (7) in Sleigh et al. (1998).                        !
!                                                                             !
!  In this subroutine, the interpolations are not carried out in the          !
!  strictest finite-element sense.  Instead, first all the major quantities   !
!  are computed at the nodes and then interpolated to the mid-faces, a        !
!  technique that I call "approximate finite element" interpolation.          !
!                                                                             !
!  Francisco Simoes, March 2004                                               !
!  Last updated (mm-dd-yyyy): 11-19-2004 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables:
  INTEGER :: i,j,p1,p2
  REAL(KIND=mp) :: h1,h2,q1,q2,u1,u2,v1,v2
  TYPE(edge) :: side

  DO j = 1,flow_edges
    i = flowedg(j)
    side = edges(i)
    p1 = side%p(1)
    p2 = side%p(2)

    ! Convective flux, continuity equation.  Later, some of these products can
    ! be put in arrays to save CPU time at the expense of added memory
    ! requirements.
    h1 = h(p1) ; u1 = u(p1) ; v1 = v(p1)
    h2 = h(p2) ; u2 = u(p2) ; v2 = v(p2)
    flux(i,1) = half*(h1*u1 + h2*u2)*side%normal(1) + &
                half*(h1*v1 + h2*v2)*side%normal(2)

    ! Convective flux, x-momentum.
    q1 = h1*u1*u1 + half*h1*h1*g
    q2 = h2*u2*u2 + half*h2*h2*g
    flux(i,2) = half*(q1 + q2)*side%normal(1)
    q1 = h1*u1*v1
    q2 = h2*u2*v2
    flux(i,2) = flux(i,2) + half*(q1 + q2)*side%normal(2)

    ! Convective flux, y-momentum.
    flux(i,3) = half*(q1 + q2)*side%normal(1)
    q1 = h1*v1*v1 + half*h1*h1*g
    q2 = h2*v2*v2 + half*h2*h2*g
    flux(i,3) = flux(i,3) + half*(q1 + q2)*side%normal(2)
  END DO

  DO j = 1,wall_edges  ! Solid walls.
    i = walledg(j)
    side = edges(i)
    p1 = side%p(1)
    p2 = side%p(2)
    h1 = h(p1)
    h2 = h(p2)
    q1 = half*(h1*h1 + h2*h2)
    q1 = half*g*q1

    flux(i,1) = zero
    flux(i,2) = q1*side%normal(1)
    flux(i,3) = q1*side%normal(2)
  END DO

END SUBROUTINE afe_flux
