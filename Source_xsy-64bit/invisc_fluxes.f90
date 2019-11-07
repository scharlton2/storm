SUBROUTINE invisc_fluxes
  USE parameters
  USE geometry
  USE dep_vars
  USE vbc_arrays
  USE options
  USE constants
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Computes the contributions to the residual due to the inviscid,            !
!  convective fluxes.  It uses the technique of Anastasiou and Chan (1997)    !
!  with modifications and improvements, especially in what concerns the       !
!  cell drying/wetting process.                                               !
!                                                                             !
!  Francisco Simoes, March 2007                                               !
!  Last updated (mm-dd-yyyy): 08-15-2016 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables.
  INTEGER :: i,j,k,l,m,n
  REAL (KIND=mp) :: atotal,dotprod,eps,f(3),h0,h1,h2,hleft,hright,incr, &
                    qfrac,sqrtghleft,sqrtghright,uedg,uleft,uright,vedg, &
                    vleft,vright,z0,zmiddle
  LOGICAL :: set_wall
  TYPE(edge) :: ework  ! Working edge.
  INTEGER, EXTERNAL :: edge_in_element,local_edge
  REAL (KIND=mp), EXTERNAL :: invert_edge_dirs

!-----------------------------------------------------------------------------!
!                                                                             !
!                         I N T E R I O R   E D G E S                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! First, compute the contributions of all the interior edges.
  DO l = 1,flow_edges1
    k = flowedg1(l)  ! Global number of the edge being processed.
    i = edges(k)%e(1)
    j = edges(k)%e(2)
    CALL edge_copy(edges(k),ework)
    eps = one
    set_wall = .FALSE.

    IF (drycell(i) .AND. drycell(j)) CYCLE  ! Both cells dry => edge dry.

    ! Compute the left-hand side state variables (from the side of element i).
    m = local_edge(grid(i),k)
    ! Make sure that the working normal points from element i to element j.
    IF (grid(i)%edge(m) < 0) eps = invert_edge_dirs(ework)
    IF (drycell(i)) THEN
      uleft = zero  ! Dry or partially dry node: use the suggestion in page 376
      vleft = zero  ! of Begnudelli and Sanders (2006) (no slip).
      hleft = zero ! h(i)
    ELSE
      hleft = h(i) + delh(i)%x*rc(i,m)%x + delh(i)%y*rc(i,m)%y
      IF (hleft < zero) THEN  ! Flag cell for bed slope treatment.
        hleft = zero
        partdry(i) = .TRUE.
      END IF
      uleft = u(i) + delu(i)%x*rc(i,m)%x + delu(i)%y*rc(i,m)%y
      vleft = v(i) + delv(i)%x*rc(i,m)%x + delv(i)%y*rc(i,m)%y
    END IF

    ! Compute the right-hand side state variables (from the side of element j).
    n = local_edge(grid(j),k)
    IF (drycell(j)) THEN
      uright = zero  ! Dry or partially dry node: use the suggestion in page
      vright = zero  ! 376 of Begnudelli and Sanders (2006) (no slip).
      hright = zero ! h(j)
    ELSE
      hright = h(j) + delh(j)%x*rc(j,n)%x + delh(j)%y*rc(j,n)%y
      IF (hright < zero) THEN  ! Flag cell for bed slope treatment.
        hright = zero
        partdry(j) = .TRUE.
      END IF
      uright = u(j) + delu(j)%x*rc(j,n)%x + delu(j)%y*rc(j,n)%y
      vright = v(j) + delv(j)%x*rc(j,n)%x + delv(j)%y*rc(j,n)%y
    END IF

    IF (hleft < h_dry .AND. hright < h_dry) CYCLE  ! Dry edge.

    ! ALGORITHM I1: algorithm to deal with the interface between dry and wet
    ! cells.  A similar procedure must be implemented for the viscous flux.
    IF (drycell(i) .OR. drycell(j)) THEN
      h0 = MAX(hright,hleft)
      z0 = half*(zvtx(ework%p(1)) + zvtx(ework%p(2)))
      h0 = h0 + z0  ! Stage at the edge.
      ! If the stage at he edge is lower than the bed elevation at the center
      ! of the dry element, set that edge flux to zero.
      IF (drycell(i)) THEN
        IF (z(i) + h_wet > h0) THEN
          set_wall = .TRUE.
          hleft = hright
          ! Here, the slip wall condition is enforced.  The magnitude of the
          ! projection of the velocity along the boundary edge is temporarily
          ! stored in variable vright.  If the no-slip wall condition is
          ! desired instead, comment out the next three lines.
          !vleft = uright*ework%tang(1) + vright*ework%tang(2)
          !uleft = vleft*ework%tang(1)
          !vleft = vleft*ework%tang(2)
        END IF
      END IF
      IF (drycell(j)) THEN
        IF (z(j) + h_wet > h0) THEN
          set_wall = .TRUE.
          hright = hleft
          ! The same as bove: the slip wall condition is enforced.  The
          ! magnitude of the projection of the velocity along the boundary edge
          ! is temporarily stored in variable vright.  If the no-slip wall
          ! condition is desired instead, comment out the next three lines.
          !vright = uleft*ework%tang(1) + vleft*ework%tang(2)
          !uright = vright*ework%tang(1)
          !vright = vright*ework%tang(2)
        END IF
      END IF
    END IF

    IF (set_wall) THEN
      ! Use physical fluxes.
      f(1) = zero
      ! The next three lines contain the method in eq. (29) of Kuiry et al.
      ! (2008).  Comment out if not desired.
      !h1 = hvtx(ework%p(1))
      !h2 = hvtx(ework%p(2))
      !hleft = SQRT((h1*h1 + h1*h2 + h2*h2)*one_third)
      hleft = MAX(hleft,hright)
      f(2) = half*g*hleft*hleft*ework%normal(1)  ! Note that hleft = hright.
      f(3) = half*g*hleft*hleft*ework%normal(2)
      !IF (set_wall) WRITE (*,'(I8,3ES14.5)') k,f(1),f(2),f(3)
    ELSE
      !IF (ABS(hleft - hright)/MAX(hleft,hright) > delta_hshock) THEN
        ! For edges with a large depth variation use a method that deals well
        ! with shocks, such as Roe's flux.
        CALL flux_ac(ework,hleft,uleft,vleft,hright,uright,vright,f)
        !CALL flux_bgnvc(ework,hleft,uleft,vleft,hright,uright,vright,f)
      !ELSE
        ! When no shocks are present use a simpler method to compute the
        ! fluxes, such as Rusanov's.
        !CALL hequivalent(ework,hleft,hright)
      !  CALL flux_rusanov(ework,hleft,uleft,vleft,wavec(i),hright,uright, &
      !                    vright,wavec(j),f)
        !CALL flux_HLL(ework,hleft,uleft,vleft,hright,uright,vright,h_dry,f)
      !END IF
    END IF

    ! For the hydrostatic correction terms of Bradford and Sanders (2000) see
    ! past implementations of SUBROUTINE invisc_fluxes, such as the one in
    ! Track xsq.

    ! Add the appropriate edge contribs to the corresponding cell residuals.
    ! First, add contributions to cell i...
    phi(i,1) = phi(i,1) - f(1)*ework%length  ! Continuity.

    ! X-momentum.
    phi(i,2) = phi(i,2) - f(2)*ework%length

    ! Y-momentum.
    phi(i,3) = phi(i,3) - f(3)*ework%length

    !... then add contributions to cell j.
    phi(j,1) = phi(j,1) + f(1)*ework%length
    phi(j,2) = phi(j,2) + f(2)*ework%length
    phi(j,3) = phi(j,3) + f(3)*ework%length

    !eps = invert_edge_dirs(ework)
    !tmp1 = uleft
    !tmp2 = vleft
    !tmp3 = hleft
    !uleft = uright
    !vleft = vright
    !hleft = hright
    !uright = tmp1
    !vright = tmp2
    !hright = tmp3
    !IF (set_wall) f(1) = zero
    !CALL flux_ac(ework,hleft,uleft,vleft,hright,uright,vright,f)

    !IF (set_wall) f(1) = zero
    !phi(j,1) = phi(j,1) - f(1)*ework%length
    !phi(j,2) = phi(j,2) - f(2)*ework%length
    !phi(j,3) = phi(j,3) - f(3)*ework%length

  END DO

!-----------------------------------------------------------------------------!
!                                                                             !
!                             W A L L   E D G E S                             !
!                                                                             !
!-----------------------------------------------------------------------------!

! Compute wall fluxes.
  SELECT CASE (btype)
  CASE (0)  ! No-slip walls.
    !CALL wenslp_fluxes(h_dry)
    !CALL wenslp2_fluxes(h_dry)
    CALL wenslp3_fluxes(h_dry)

  CASE (1)  ! Slip walls.
    !CALL weslp_fluxes(h_dry)
    CALL weslp2_fluxes(h_dry)

  END SELECT

!-----------------------------------------------------------------------------!
!                                                                             !
!                   I N F L O W / O U T F L O W   E D G E S                   !
!                                                                             !
!-----------------------------------------------------------------------------!

! Finally, enforce the inflow and outflow boundary conditions.
! Subcritical inflow: velocity specified.

  ! Enforce discharge along the specified chain of grid edges according to
  ! flow depth (vtype = 1 and 5).
  CALL inflowQbyH
  ! Stage is specified at the inlet boundary, instead of discharge (vtype = 4).
  CALL inflowH

!-----------------------------------------------------------------------------!
!                                                                             !
!-----------------------------------------------------------------------------!
! Subcritical outflow: stage specified.

  CALL outflowbyH
  !CALL outflowbyH2
  !CALL outflowbyH3
  CALL crit_outflow  ! Fr = 1.0 at outlet.

!-----------------------------------------------------------------------------!
!                                                                             !
!-----------------------------------------------------------------------------!
! Free flow nodes (e.g., supercritical outflow).  No flow or stage boundary
! conditions are imposed on these edges.

  CALL free_flow

END SUBROUTINE invisc_fluxes
