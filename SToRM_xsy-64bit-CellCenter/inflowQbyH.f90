SUBROUTINE inflowQbyH
  USE parameters
  USE geometry
  USE dep_vars
  USE vbc_arrays
  USE options
  USE constants
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Enforce discharge along the specified chain of grid edges according to     !
!  flow depth.  Uses Riemann invariants and an iterative procedure as per     !
!  the equations in page 682 of Yoon and Kang (2004).                         !
!                                                                             !
!  Francisco Simoes, April 2008                                               !
!  Last updated (mm-dd-yyyy): 04-18-2013 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables.
  INTEGER :: counter,i,icycle,k,l,m
  REAL (KIND=mp) :: atotal,eps,f(3),hleft,hright,newatotal,qfrac, &
                    qtotal,sqrtghleft,uedg,uleft,uright,vedg,vleft,vright
  TYPE(edge) :: ework  ! Working edge.
  LOGICAL :: null_area
  INTEGER, EXTERNAL :: edge_in_element,local_edge
  REAL (KIND=mp), EXTERNAL :: invert_edge_dirs
  LOGICAL, EXTERNAL :: equals

  DO icycle = 1,n_inflowbdr
    IF (vtype(icycle) /= 1) CYCLE
    ! Zero inflow boundaries should be implemented as free slip walls.  For
    ! now, the inflow edges are simply skipped in the calculation of the
    ! inviscid fluxes.
    IF (qin(icycle) < vsmall) CYCLE

    ! First guess the total flow area using the stage computed from the inside
    ! of the computational domain.
    atotal = zero
    DO l = 1,bcedgesinflow(0,icycle)
      k = bcedgesinflow(l,icycle)
      i = edge_in_element(edges(k))
      m = local_edge(grid(i),k)

      !IF (drycell(i)) CYCLE  ! Dry cell.

      ! Use the stage computed from the inside of the computational domain.
      IF (drycell(i)) THEN
        hleft = zero
      ELSE
        hleft = h(i) + delh(i)%x*rc(i,m)%x + delh(i)%y*rc(i,m)%y
        hleft = MAX(hleft,zero)
      END IF
      atotal = atotal + hleft*edges(k)%length
    END DO

    null_area = equals(atotal,zero)

    ! Iterative procedure.
    counter = 0
    DO
      newatotal = zero
      qtotal = zero
      DO l = 1,bcedgesinflow(0,icycle)
        k = bcedgesinflow(l,icycle)
        i = edge_in_element(edges(k))
        CALL edge_copy(edges(k),ework)
        eps = one

        !IF (drycell(i)) CYCLE  ! Dry cell.

        ! Compute the left-hand side state variables from (element i).
        m = local_edge(grid(i),k)
        hleft = h(i) + delh(i)%x*rc(i,m)%x + delh(i)%y*rc(i,m)%y
        IF (null_area) hleft = bch_max  ! For cases where atotal = 0.
        IF (hleft < h_dry) THEN  ! Dry edge.
          IF (hleft < zero) partdry(i) = .TRUE.
          !CYCLE
          hleft = MAX(hleft,zero)
          uedg = zero
          vedg = zero
        ELSE
          uedg = u(i) + delu(i)%x*rc(i,m)%x + delu(i)%y*rc(i,m)%y
          vedg = v(i) + delv(i)%x*rc(i,m)%x + delv(i)%y*rc(i,m)%y
        END IF
        ! This means that the normal points into element i: reverse it.
        IF (grid(i)%edge(m) < 0) eps = invert_edge_dirs(ework)

        ! Compute the right-hand side state boundary conditions.  Note that
        ! here the stage is taken from outside the domain, hright, which is
        ! inconsistent with the computation of atotal.  Nonetheless, this is
        ! an approximation that isn't far from the truth, specially for steady
        ! state computations.  qfrac if the fraction of the discharge through
        ! the edge...
        !qfrac = qin(icycle)*hleft*edges(k)%length/atotal
        !qfrac = qfrac/(hleft*edges(k)%length) !..and now is the flow velocity.
        qfrac = qin(icycle)/atotal
        uright = -qfrac
        uleft = uedg*ework%normal(1) + vedg*ework%normal(2)
        sqrtghleft = SQRT(g*hleft)
        hright = 2.0_mp*sqrtghleft + uleft - uright
        hright = hright*hright/4.0_mp/g  ! This uses the Riemann invariant...
        !hright = hleft  !..and this simply sets a zero depth-gradient.
        uright = -qfrac*ework%normal(1)
        vright = -qfrac*ework%normal(2)
        uleft = uedg
        vleft = vedg
        IF (hright < h_dry) CYCLE  ! Dry edge.

        newatotal = newatotal + hright*edges(k)%length
        qtotal = qtotal + (uright*ework%normal(1) + &
                           vright*ework%normal(2))*hright*ework%length
      END DO

      ! Convergence criterion to exit the DO-loop.
      !IF (ABS((newatotal - atotal)/MAX(newatotal,atotal)) < 0.02 .OR. &
      !  counter > 5) EXIT
      qtotal = ABS(qtotal)
      IF (ABS(qtotal - qin(icycle))/qin(icycle) < 0.005_mp .OR. counter > 50) &
        EXIT

      ! Fixed point iteration: update area and cycle.
      counter = counter + 1
      atotal = newatotal
    END DO  ! End of main iteration cycle.
    !PRINT *,counter,qtotal

    ! This DO-loop finally applies the boundary condition.
    atotal = newatotal
    DO l = 1,bcedgesinflow(0,icycle)
      k = bcedgesinflow(l,icycle)
      i = edge_in_element(edges(k))
      CALL edge_copy(edges(k),ework)
      eps = one

      !IF (drycell(i)) CYCLE  ! Dry cell.

      ! Compute the left-hand side state variables from (element i).
      m = local_edge(grid(i),k)
      hleft = h(i) + delh(i)%x*rc(i,m)%x + delh(i)%y*rc(i,m)%y
      IF (null_area) hleft = bch_max  ! For cases where atotal = 0.
      IF (hleft < h_dry) THEN  ! Dry edge.
        IF (hleft < zero) partdry(i) = .TRUE.
        !CYCLE
        hleft = MAX(hleft,zero)
        uedg = zero
        vedg = zero
      ELSE
        uedg = u(i) + delu(i)%x*rc(i,m)%x + delu(i)%y*rc(i,m)%y
        vedg = v(i) + delv(i)%x*rc(i,m)%x + delv(i)%y*rc(i,m)%y
      END IF
      ! This means that the normal points into element i: reverse it.
      IF (grid(i)%edge(m) < 0) eps = invert_edge_dirs(ework)

      ! Compute the right-hand side state boundary conditions.  Note that
      ! here the stage is taken from outside the domain, hright, which is
      ! inconsistent with the computation of atotal.  Nonetheless, this is
      ! an approximation that isn't far from the truth, specially for steady
      ! state computations.  qfrac if the fraction of the discharge through
      ! the edge...
      !qfrac = qin(icycle)*hleft*edges(k)%length/atotal
      !qfrac = qfrac/(hleft*edges(k)%length) !..and now is the flow velocity.
      qfrac = qin(icycle)/atotal
      uright = -qfrac
      uleft = uedg*ework%normal(1) + vedg*ework%normal(2)
      sqrtghleft = SQRT(g*hleft)
      hright = 2.0_mp*sqrtghleft + uleft - uright
      hright = hright*hright/4.0_mp/g  ! This uses the Riemann invariant...
      !hright = hleft  !..and this simply sets a zero depth-gradient.
      uright = -qfrac*ework%normal(1)
      vright = -qfrac*ework%normal(2)
      uleft = uedg
      vleft = vedg
      IF (hright < h_dry) CYCLE  ! Dry edge.

      ! Fluxes enforced exactly -- eq. (35) of Yoon and Kang (2004).
      f(1) = hright*(uright*ework%normal(1) + vright*ework%normal(2))
      f(2) = hright*uright*(uright*ework%normal(1) + vright*ework%normal(2)) + &
             half*g*hright*hright*ework%normal(1)
      f(3) = hright*vright*(uright*ework%normal(1) + vright*ework%normal(2)) + &
             half*g*hright*hright*ework%normal(2)

      !IF (ABS(hleft - hright)/MAX(hleft,hright) > delta_hshock) THEN
      !  ! For edges with a large depth variation use a method that deals well
      !  ! with shocks, such as Roe's flux.
      !  CALL flux_ac(ework,hleft,uleft,vleft,hright,uright,vright,f)
      !  !CALL flux_bgnvc(ework,hleft,uleft,vleft,hright,uright,vright,f)
      !ELSE
      !  ! When no shocks are present use a simpler method to compute the
      !  ! fluxes, such as Rusanov's.  Note that the use of the equivalent
      !  ! depth invalidates the above hleft/hright computation.
      !  CALL hequivalent(ework,hleft,hright)
      !  CALL flux_rusanov(ework,hleft,uleft,vleft,wavec(i),hright,uright, &
      !                    vright,SQRT(g*hright),f)
      !  !CALL flux_HLL(ework,hleft,uleft,vleft,hright,uright,vright,h_dry,f)
      !END IF

      ! Debugging discharges, computed from the left and from the right of the
      ! inflow edges.
      !cqin(icycle) = cqin(icycle) + (uright*ework%normal(1) + &
      !               vright*ework%normal(2))*hright*ework%length
      !cqin2(icycle) = cqin2(icycle) + (uleft*ework%normal(1) + &
      !                vleft*ework%normal(2))*hleft*ework%length
      cqin3(icycle) = cqin3(icycle) + f(1)*ework%length

      ! Add the edge flux to the residual's array for element i.
      phi(i,1) = phi(i,1) - f(1)*ework%length
      phi(i,2) = phi(i,2) - f(2)*ework%length
      phi(i,3) = phi(i,3) - f(3)*ework%length
    END DO
  END DO

! This DO-loop is identical to the previous DO-loop in every way except that it
! is done for vtype(icycle) = 5, where inflow only occurs where triangles are
! wet.  The difference in the code resides in the three statements that check
! for dry cells:
!
!   IF (drycell(i)) CYCLE  ! Dry cell.
!
  DO icycle = 1,n_inflowbdr
    IF (vtype(icycle) /= 5) CYCLE
    ! Zero inflow boundaries should be implemented as free slip walls.  For
    ! now, the inflow edges are simply skipped in the calculation of the
    ! inviscid fluxes.
    IF (qin(icycle) < vsmall) CYCLE

    ! First guess the total flow area using the stage computed from the inside
    ! of the computational domain.
    atotal = zero
    DO l = 1,bcedgesinflow(0,icycle)
      k = bcedgesinflow(l,icycle)
      i = edge_in_element(edges(k))
      m = local_edge(grid(i),k)

      IF (drycell(i)) CYCLE  ! Dry cell.

      ! Use the stage computed from the inside of the computational domain.
      IF (drycell(i)) THEN
        hleft = zero
      ELSE
        hleft = h(i) + delh(i)%x*rc(i,m)%x + delh(i)%y*rc(i,m)%y
        hleft = MAX(hleft,zero)
      END IF
      atotal = atotal + hleft*edges(k)%length
    END DO

    null_area = equals(atotal,zero)

    ! Iterative procedure.
    counter = 0
    DO
      newatotal = zero
      qtotal = zero
      DO l = 1,bcedgesinflow(0,icycle)
        k = bcedgesinflow(l,icycle)
        i = edge_in_element(edges(k))
        CALL edge_copy(edges(k),ework)
        eps = one

        IF (drycell(i)) CYCLE  ! Dry cell.

        ! Compute the left-hand side state variables from (element i).
        m = local_edge(grid(i),k)
        hleft = h(i) + delh(i)%x*rc(i,m)%x + delh(i)%y*rc(i,m)%y
        IF (null_area) hleft = bch_max  ! For cases where atotal = 0.
        IF (hleft < h_dry) THEN  ! Dry edge.
          IF (hleft < zero) partdry(i) = .TRUE.
          !CYCLE
          hleft = MAX(hleft,zero)
          uedg = zero
          vedg = zero
        ELSE
          uedg = u(i) + delu(i)%x*rc(i,m)%x + delu(i)%y*rc(i,m)%y
          vedg = v(i) + delv(i)%x*rc(i,m)%x + delv(i)%y*rc(i,m)%y
        END IF
        ! This means that the normal points into element i: reverse it.
        IF (grid(i)%edge(m) < 0) eps = invert_edge_dirs(ework)

        ! Compute the right-hand side state boundary conditions.  Note that
        ! here the stage is taken from outside the domain, hright, which is
        ! inconsistent with the computation of atotal.  Nonetheless, this is
        ! an approximation that isn't far from the truth, specially for steady
        ! state computations.  qfrac if the fraction of the discharge through
        ! the edge...
        !qfrac = qin(icycle)*hleft*edges(k)%length/atotal
        !qfrac = qfrac/(hleft*edges(k)%length) !..and now is the flow velocity.
        qfrac = qin(icycle)/atotal
        uright = -qfrac
        uleft = uedg*ework%normal(1) + vedg*ework%normal(2)
        sqrtghleft = SQRT(g*hleft)
        hright = 2.0_mp*sqrtghleft + uleft - uright
        hright = hright*hright/4.0_mp/g  ! This uses the Riemann invariant...
        !hright = hleft  !..and this simply sets a zero depth-gradient.
        uright = -qfrac*ework%normal(1)
        vright = -qfrac*ework%normal(2)
        uleft = uedg
        vleft = vedg
        IF (hright < h_dry) CYCLE  ! Dry edge.

        newatotal = newatotal + hright*edges(k)%length
        qtotal = qtotal + (uright*ework%normal(1) + &
                           vright*ework%normal(2))*hright*ework%length
      END DO

      ! Convergence criterion to exit the DO-loop.
      !IF (ABS((newatotal - atotal)/MAX(newatotal,atotal)) < 0.02 .OR. &
      !  counter > 5) EXIT
      qtotal = ABS(qtotal)
      IF (ABS(qtotal - qin(icycle))/qin(icycle) < 0.005_mp .OR. counter > 50) &
        EXIT

      ! Fixed point iteration: update area and cycle.
      counter = counter + 1
      atotal = newatotal
    END DO  ! End of main iteration cycle.
    !PRINT *,counter,qtotal

    ! This DO-loop finally applies the boundary condition.
    atotal = newatotal
    DO l = 1,bcedgesinflow(0,icycle)
      k = bcedgesinflow(l,icycle)
      i = edge_in_element(edges(k))
      CALL edge_copy(edges(k),ework)
      eps = one

      IF (drycell(i)) CYCLE  ! Dry cell.

      ! Compute the left-hand side state variables from (element i).
      m = local_edge(grid(i),k)
      hleft = h(i) + delh(i)%x*rc(i,m)%x + delh(i)%y*rc(i,m)%y
      IF (null_area) hleft = bch_max  ! For cases where atotal = 0.
      IF (hleft < h_dry) THEN  ! Dry edge.
        IF (hleft < zero) partdry(i) = .TRUE.
        !CYCLE
        hleft = MAX(hleft,zero)
        uedg = zero
        vedg = zero
      ELSE
        uedg = u(i) + delu(i)%x*rc(i,m)%x + delu(i)%y*rc(i,m)%y
        vedg = v(i) + delv(i)%x*rc(i,m)%x + delv(i)%y*rc(i,m)%y
      END IF
      ! This means that the normal points into element i: reverse it.
      IF (grid(i)%edge(m) < 0) eps = invert_edge_dirs(ework)

      ! Compute the right-hand side state boundary conditions.  Note that
      ! here the stage is taken from outside the domain, hright, which is
      ! inconsistent with the computation of atotal.  Nonetheless, this is
      ! an approximation that isn't far from the truth, specially for steady
      ! state computations.  qfrac if the fraction of the discharge through
      ! the edge...
      !qfrac = qin(icycle)*hleft*edges(k)%length/atotal
      !qfrac = qfrac/(hleft*edges(k)%length) !..and now is the flow velocity.
      qfrac = qin(icycle)/atotal
      uright = -qfrac
      uleft = uedg*ework%normal(1) + vedg*ework%normal(2)
      sqrtghleft = SQRT(g*hleft)
      hright = 2.0_mp*sqrtghleft + uleft - uright
      hright = hright*hright/4.0_mp/g  ! This uses the Riemann invariant...
      !hright = hleft  !..and this simply sets a zero depth-gradient.
      uright = -qfrac*ework%normal(1)
      vright = -qfrac*ework%normal(2)
      uleft = uedg
      vleft = vedg
      IF (hright < h_dry) CYCLE  ! Dry edge.

      ! Fluxes enforced exactly -- eq. (35) of Yoon and Kang (2004).
      f(1) = hright*(uright*ework%normal(1) + vright*ework%normal(2))
      f(2) = hright*uright*(uright*ework%normal(1) + vright*ework%normal(2)) + &
             half*g*hright*hright*ework%normal(1)
      f(3) = hright*vright*(uright*ework%normal(1) + vright*ework%normal(2)) + &
             half*g*hright*hright*ework%normal(2)

      !IF (ABS(hleft - hright)/MAX(hleft,hright) > delta_hshock) THEN
      !  ! For edges with a large depth variation use a method that deals well
      !  ! with shocks, such as Roe's flux.
      !  CALL flux_ac(ework,hleft,uleft,vleft,hright,uright,vright,f)
      !  !CALL flux_bgnvc(ework,hleft,uleft,vleft,hright,uright,vright,f)
      !ELSE
      !  ! When no shocks are present use a simpler method to compute the
      !  ! fluxes, such as Rusanov's.  Note that the use of the equivalent
      !  ! depth invalidates the above hleft/hright computation.
      !  CALL hequivalent(ework,hleft,hright)
      !  CALL flux_rusanov(ework,hleft,uleft,vleft,wavec(i),hright,uright, &
      !                    vright,SQRT(g*hright),f)
      !  !CALL flux_HLL(ework,hleft,uleft,vleft,hright,uright,vright,h_dry,f)
      !END IF

      ! Debugging discharges, computed from the left and from the right of the
      ! inflow edges.
      !cqin(icycle) = cqin(icycle) + (uright*ework%normal(1) + &
      !               vright*ework%normal(2))*hright*ework%length
      !cqin2(icycle) = cqin2(icycle) + (uleft*ework%normal(1) + &
      !                vleft*ework%normal(2))*hleft*ework%length
      cqin3(icycle) = cqin3(icycle) + f(1)*ework%length

      ! Add the edge flux to the residual's array for element i.
      phi(i,1) = phi(i,1) - f(1)*ework%length
      phi(i,2) = phi(i,2) - f(2)*ework%length
      phi(i,3) = phi(i,3) - f(3)*ework%length
    END DO
  END DO

END SUBROUTINE inflowQbyH
