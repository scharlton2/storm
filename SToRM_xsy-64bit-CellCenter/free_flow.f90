SUBROUTINE free_flow
  USE parameters
  USE geometry
  USE dep_vars
  USE vbc_arrays
  USE options
  USE constants
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Free flow nodes (e.g., supercritical outflow).  No flow or stage boundary  !
!  conditions are imposed on these edges.                                     !
!                                                                             !
!  Francisco Simoes, April 2009                                               !
!  Last updated (mm-dd-yyyy): 02-28-2011 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables.
  INTEGER :: i,icycle,k,l,m,n
  REAL (KIND=mp) :: dotprod,eps,f(3),hleft,hright,sqrtghleft, &
                    sqrtghright,uedg,uleft,uright,vedg,vleft,vright,zmiddle
  LOGICAL :: set_wall
  TYPE(edge) :: ework  ! Working edge.
  INTEGER, EXTERNAL :: edge_in_element,local_edge
  REAL (KIND=mp), EXTERNAL :: invert_edge_dirs
  LOGICAL, EXTERNAL :: equals

  DO icycle = 1,n_outflowbdr
    IF (htype(icycle) /= 0) CYCLE

    DO l = 1,bcedgesoutflow(0,icycle)
      k = bcedgesoutflow(l,icycle)
      i = edge_in_element(edges(k))
      CALL edge_copy(edges(k),ework)
      eps = one

      ! Compute the left-hand side state variables from (element i).
      m = local_edge(grid(i),k)
      hleft = h(i) + delh(i)%x*rc(i,m)%x + delh(i)%y*rc(i,m)%y
      hleft = MAX(hleft,zero)
      IF (hleft < h_dry) CYCLE  ! Dry edge.

      IF (h(i) < h_dry) THEN
        uedg = zero
        vedg = zero
      ELSE
        uedg = u(i) + delu(i)%x*rc(i,m)%x + delu(i)%y*rc(i,m)%y
        vedg = v(i) + delv(i)%x*rc(i,m)%x + delv(i)%y*rc(i,m)%y
      END IF
      ! This means that the normal points into element i: reverse it.
      IF (grid(i)%edge(m) < 0) eps = invert_edge_dirs(ework)
      uleft = uedg ; vleft = vedg

      ! Compute the right-hand side state boundary conditions.
      hright = hleft
      uright = uleft
      vright = vleft

      set_wall = .FALSE.  ! Set 'temporary wall' variable.
      IF (hvalve) THEN
        dotprod = uright*ework%normal(1) + vright*ework%normal(2)
        IF (dotprod < 0) THEN  ! Apply valve to prevent inflow through this
          set_wall = .TRUE.   ! boundary.
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
        !CALL hequivalent(ework,hleft,hright)
      !  CALL flux_rusanov(ework,hleft,uleft,vleft,wavec(i),hright,uright, &
      !                    vright,SQRT(g*hright),f)
        !CALL flux_HLL(ework,hleft,uleft,vleft,hright,uright,vright,h_dry,f)
      !END IF

      ! Add the edge flux to the residual's array for element i.
      IF (set_wall) f(1) = zero
      DO n = 1,3
        phi(i,n) = phi(i,n) - f(n)*ework%length
      END DO

    END DO
  END DO

END SUBROUTINE free_flow
