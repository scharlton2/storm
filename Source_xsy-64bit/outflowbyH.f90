SUBROUTINE outflowbyH
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
!  Uses the Riemann invariants as shown in Anastasiou and Chan (1997).        !
!                                                                             !
!  Francisco Simoes, April 2008                                               !
!  Last updated (mm-dd-yyyy): 02-27-2011 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables.
  INTEGER :: i,icycle,k,l,m
  REAL (KIND=mp) :: dotprod,eps,f(3),hleft,hright,sqrtghleft, &
                    sqrtghright,uedg,uleft,uright,vedg,vleft,vright,zmiddle
  LOGICAL :: set_wall
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
        uright = uleft + (sqrtghleft - sqrtghright)
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
      set_wall = .FALSE.  ! Set 'temporary wall' variable.
      IF (hvalve) THEN
        dotprod = uright*ework%normal(1) + vright*ework%normal(2)
        IF (dotprod < zero) THEN
          set_wall = .TRUE.
          hright = hleft
          uright = zero
          vright = zero
        END IF
      END IF

      !IF (ABS(hleft - hright)/MAX(hleft,hright) > delta_hshock) THEN
        ! For edges with a large depth variation use a method that deals well
        ! with shocks, such as Roe's flux.
        CALL flux_ac(ework,hleft,uleft,vleft,hright,uright,vright,f)
        !CALL flux_bgnvc(ework,hleft,uleft,vleft,hright,uright,vright,f)
      !ELSE
        ! When no shocks are present use a simpler method to compute the
        ! fluxes, such as Rusanov's.  Note that the use of the equivalent
        ! depth invalidates the above hleft/hright computation.
      !  CALL hequivalent(ework,hleft,hright)
      !  CALL flux_rusanov(ework,hleft,uleft,vleft,wavec(i),hright,uright, &
      !                  vright,SQRT(g*hright),f)
        !CALL flux_HLL(k,ework,hleft,uleft,vleft,hright,uright,vright,h_dry,f)
      !END IF

      ! Add the edge flux to the residual's array for element i.
      !IF (hvalve .AND. dotprod < zero) f(1) = zero
      IF (set_wall) f(1) = zero

      ! Debugging discharges, computed from the left and from the right of the
      ! outflow edges.
      !cqout(icycle) = cqout(icycle) + (uright*ework%normal(1) + &
      !                vright*ework%normal(2))*hright*ework%length
      !cqout2(icycle) = cqout2(icycle) + (uleft*ework%normal(1) + &
      !                 vleft*ework%normal(2))*hleft*ework%length
      cqout3(icycle) = cqout3(icycle) + f(1)*ework%length

      phi(i,1) = phi(i,1) - f(1)*ework%length
      phi(i,2) = phi(i,2) - f(2)*ework%length
      phi(i,3) = phi(i,3) - f(3)*ework%length
    END DO
  END DO

END SUBROUTINE outflowbyH
