SUBROUTINE outflowbyH2
  USE parameters
  USE geometry
  USE dep_vars
  USE vbc_arrays
  USE options
  USE constants
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Subcritical outflow boundary condition, water depth (or stage) specified.  !
!  Uses the Riemann invariants as shown in Yoon and Kang (2004). The fluxes   !
!  are computed directly.                                                     !
!                                                                             !
!  Francisco Simoes, April 2008                                               !
!  Last updated (mm-dd-yyyy): 02-25-2011 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables.
  INTEGER :: i,icycle,k,l,m
  REAL (KIND=mp) :: dh,dotprod,eps,f(3),hleft,hright,sqrtghleft, &
                    sqrtghright,uedg,uleft,uright,vedg,vleft,vright,zmiddle
  TYPE(edge) :: ework  ! Working edge.
  INTEGER, EXTERNAL :: edge_in_element,local_edge
  REAL (KIND=mp), EXTERNAL :: invert_edge_dirs
  LOGICAL, EXTERNAL :: equals

  DO icycle = 1,n_outflowbdr
    IF (htype(icycle) /= 1) CYCLE

    DO l = 1,bcedgesoutflow(0,icycle)
      k = bcedgesoutflow(l,icycle)
      i = edge_in_element(edges(k))
      CALL edge_copy(edges(k),ework)
      eps = one

      ! Compute the left-hand side state variables from (element i).
      m = local_edge(grid(i),k)
      hleft = h(i) + delh(i)%x*rc(i,m)%x + delh(i)%y*rc(i,m)%y
      IF (drycell(i) .OR. hleft < h_dry) THEN
        uedg = zero
        vedg = zero
        IF ((.NOT. drycell(i)) .AND. hleft < zero) partdry(i) = .TRUE.
        hleft = MAX(hleft,zero)
      ELSE
        uedg = u(i) + delu(i)%x*rc(i,m)%x + delu(i)%y*rc(i,m)%y
        vedg = v(i) + delv(i)%x*rc(i,m)%x + delv(i)%y*rc(i,m)%y
      END IF
      ! This means that the normal points into element i: reverse it.
      IF (grid(i)%edge(m) < 0) eps = invert_edge_dirs(ework)

      ! In this implementation, the prescribed stage at the boundary nodes is
      ! enforced in a gradual manner.  The relaxation parameter is 'alpha' and
      ! 'kappa' is the minimum increment allowed.  Set kappa = [0.001,0.01] for
      ! laboratory channels, and kappa = [0.01,0.1] for larger channels and
      ! natural rivers. alpha = [0,1].
      !zmiddle = half*(zvtx(ework%p(1)) + zvtx(ework%p(2)))
      !hright = MAX(zero,hbc(icycle) - zmiddle)
      !incr = alpha*(hright - hleft)
      ! To avoid infinitesimal steps:
      !IF (ABS(incr) < kappa) incr = hright - hleft
      !hright = hright - incr

      ! In the following code, the prescribed h is enforced instantaneously.
      zmiddle = half*(zvtx(ework%p(1)) + zvtx(ework%p(2)))
      hright = MAX(zero,hbc(icycle) - zmiddle)

      IF (hleft < h_dry .AND. hright < h_dry) CYCLE  ! Dry edge.

      IF (hright < h_dry) THEN  ! Treat as a solid wall.
        uright = zero ; uleft = uedg
        vright = zero ; vleft = vedg
      ELSE
        uleft = uedg*ework%normal(1) + vedg*ework%normal(2)
        vleft = uedg*ework%tang(1) + vedg*ework%tang(2)
        sqrtghleft = SQRT(g*hleft)
        sqrtghright = SQRT(g*hright)
        uright = uleft + 2.0_mp*(sqrtghleft - sqrtghright)
        vright = vleft
        uleft = uedg
        vleft = vedg
        uedg = uright
        vedg = vright
        uright = uedg*ework%normal(1) + vedg*ework%tang(1)
        vright = uedg*ework%normal(2) + vedg*ework%tang(2)
      END IF

      ! Prevent flow reversal into the domain.  This is a numerical valve that
      ! sets the normal velocity to zero whenever the velocity vector points
      ! into the domain.  Here, it is set to act like a solid wall with no slip
      ! when the velocity vector at the boundary points into the computational
      ! domain.
      IF (hvalve) THEN
        dotprod = uright*ework%normal(1) + vright*ework%normal(2)
        IF (dotprod < zero) THEN
          hright = hleft
          uright = zero
          vright = zero
        END IF
      END IF

      ! Fluxes enforced exactly -- eq. (35) of Yoon and Kang (2004).
      f(1) = hright*(uright*ework%normal(1) + vright*ework%normal(2))
      f(2) = hright*uright*(uright*ework%normal(1) + &
             vright*ework%normal(2)) + half*g*hright*hright*ework%normal(1)
      f(3) = hright*vright*(uright*ework%normal(1) + &
             vright*ework%normal(2)) + half*g*hright*hright*ework%normal(2)

      ! Debugging discharges, computed from the left and from the right of the
      ! outflow edges.
      !cqout(icycle) = cqout(icycle) + (uright*ework%normal(1) + &
      !                vright*ework%normal(2))*hright*ework%length
      !cqout2(icycle) = cqout2(icycle) + (uleft*ework%normal(1) + &
      !                 vleft*ework%normal(2))*hleft*ework%length
      cqout3(icycle) = cqout3(icycle) + f(1)*ework%length

      ! Add the edge flux to the residual's array for element i.
      dh = hvtx(ework%p(1)) - hvtx(ework%p(2))
      dh = g*dh*dh/12.0_mp
      dh = zero

      phi(i,1) = phi(i,1) - f(1)*ework%length
      !phi(i,2) = phi(i,2) - (f(2) + dh*rc(i,m)%x)*ework%length
      !phi(i,3) = phi(i,3) - (f(3) + dh*rc(i,m)%y)*ework%length
      phi(i,2) = phi(i,2) - (f(2) + dh*ework%normal(1))*ework%length
      phi(i,3) = phi(i,3) - (f(3) + dh*ework%normal(2))*ework%length
    END DO
  END DO

END SUBROUTINE outflowbyH2
