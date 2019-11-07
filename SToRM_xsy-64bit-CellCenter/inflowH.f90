SUBROUTINE inflowH
  USE parameters
  USE geometry
  USE dep_vars
  USE vbc_arrays
  USE options
  USE constants
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  In this type of inflow boundary, the stage is specified at the boundary    !
!  and the velocity is then computed using the Riemann invariants.            !
!                                                                             !
!  Francisco Simoes, April 2009                                               !
!  Last updated (mm-dd-yyyy): 02-28-2011 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables.
  INTEGER :: counter,i,icycle,k,l,m,n
  REAL (KIND=mp) :: eps,f(3),hleft,hright,sqrtghleft,sqrtghright,uedg, &
                    uleft,uright,vedg,vleft,vright,zmiddle
  TYPE(edge) :: ework  ! Working edge.
  INTEGER, EXTERNAL :: edge_in_element,local_edge
  REAL (KIND=mp), EXTERNAL :: invert_edge_dirs

  DO icycle = 1,n_inflowbdr
    IF (vtype(icycle) /= 4) CYCLE

    DO l = 1,bcedgesinflow(0,icycle)
      k = bcedgesinflow(l,icycle)
      i = edge_in_element(edges(k))
      CALL edge_copy(edges(k),ework)
      eps = one

      ! Compute the left-hand side state variables from (element i).
      m = local_edge(grid(i),k)
      hleft = h(i) + delh(i)%x*rc(i,m)%x + delh(i)%y*rc(i,m)%y
      hleft = MAX(hleft,zero)
      IF (h(i) < h_dry) THEN
        uedg = zero
        vedg = zero
      ELSE
        uedg = u(i) + delu(i)%x*rc(i,m)%x + delu(i)%y*rc(i,m)%y
        vedg = v(i) + delv(i)%x*rc(i,m)%x + delv(i)%y*rc(i,m)%y
      END IF
      ! This means that the normal points into element i: reverse it.
      IF (grid(i)%edge(m) < 0) eps = invert_edge_dirs(ework)

      ! Compute the right-hand side state boundary conditions. 'qin()' is the
      ! stage at the inflow boundary.
      zmiddle = half*(zvtx(ework%p(1)) + zvtx(ework%p(2)))
      hright = MAX(zero,qin(icycle) - zmiddle)
      uleft = uedg*ework%normal(1) + vedg*ework%normal(2)
      sqrtghleft = SQRT(g*hleft)
      sqrtghright = SQRT(g*hright)
      uright = uleft + (sqrtghleft - sqrtghright)
      vright = uright*ework%normal(2)
      uright = uright*ework%normal(1)
      uleft = uedg
      vleft = vedg
      IF (hleft < h_dry .AND. hright < h_dry) CYCLE  ! Dry edge.

      !IF (ABS(hleft - hright)/MAX(hleft,hright) > delta_hshock) THEN
        ! For edges with a large depth variation use a method that deals well
        ! with shocks, such as Roe's flux.
        CALL flux_ac(ework,hleft,uleft,vleft,hright,uright,vright,f)
        !CALL flux_bgnvc(ework,hleft,uleft,vleft,hright,uright,vright,f)
      !ELSE
        ! When no shocks are present use a simpler method to compute the
        ! fluxes, such as Rusanov's.  Note that the use of the equivalent
        ! depth invalidates the above hleft/hright computation.
        !CALL hequivalent(ework,hleft,hright)
      !  CALL flux_rusanov(ework,hleft,uleft,vleft,wavec(i),hright,uright, &
      !                    vright,SQRT(g*hright),f)
        !CALL flux_HLL(ework,hleft,uleft,vleft,hright,uright,vright,h_dry,f)
      !END IF

      ! Debugging discharges, computed from the left and from the right of the
      ! inflow edges.
      !cqin(icycle) = cqin(icycle) + (uright*ework%normal(1) + &
      !               vright*ework%normal(2))*hright*ework%length
      !cqin2(icycle) = cqin2(icycle) + (uleft*ework%normal(1) + &
      !                vleft*ework%normal(2))*hleft*ework%length
      cqin3(icycle) = cqin3(icycle) + f(1)*ework%length

      ! Add the edge flux to the residual's array for element i.
      DO n = 1,3
        phi(i,n) = phi(i,n) - f(n)*ework%length
      END DO
    END DO
  END DO

END SUBROUTINE inflowH
