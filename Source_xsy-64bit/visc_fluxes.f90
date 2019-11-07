SUBROUTINE visc_fluxes
  USE geometry
  USE dep_vars
  USE options
  USE constants
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  The viscous fluxes are computed here on a edge-by-edge basis.  The method  !
!  used follows that of Anastasiou and Chan (1997) but using shape functions  !
!  to calculate the gradients in the triangles of figs. 2(a) and 2(b).  The   !
!  viscous fluxes are calculated from the velocity gradients at the edges     !
!  and using Gauss's theorem to integrate them in the element using the       !
!  mid-point rule on each edge.  More details about the form of the           !
!  governing equation and this type of term integration can be seen in sec.   !
!  2 of Anastasiou and Chan (1997).                                           !
!                                                                             !
!  Francisco Simoes, July 2007                                                !
!  Last updated (mm-dd-yyyy): 02-04-2009 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables.
  INTEGER :: i,j,k,k1,k2,l,m,p1,p2,p3
  REAL (KIND=mp) :: area,b(3),c(3),evisc,eps,fv(3),gv(3),h0,hedge,nx,ny,x3, &
                    xc,y3,yc,z0
  TYPE(vector) :: gradu,gradv
  TYPE(edge) :: ework
  INTEGER, EXTERNAL :: edge_in_element,local_edge
  REAL (KIND=mp), EXTERNAL :: invert_edge_dirs

! First do the boundary edges.
  DO i = 1,n_bpolygon
    j = bpolygon(i)  ! Edge.
    k = edge_in_element(edges(j))  ! Triangle.
    l = local_edge(grid(k),j)
    CALL edge_copy(edges(j),ework)
    IF (grid(k)%edge(l) < 0) THEN  ! Make sure circulation is counterclockwise.
      p2 = edges(j)%p(2)
      p3 = edges(j)%p(1)
      eps = invert_edge_dirs(ework)
    ELSE
      p2 = edges(j)%p(1)
      p3 = edges(j)%p(2)
    END IF

    hedge = half*(hvtx(p2) + hvtx(p3))
    IF (hedge < h_dry) CYCLE  ! Dry edge.

    ! Prepare constants needed in shape functions.
    area = grid(k)%area*one_third*2.0_mp  ! Area to use in shape functions.
    b(1) = nodes(p2)%y - nodes(p3)%y ; c(1) = nodes(p3)%x - nodes(p2)%x
    b(2) = nodes(p3)%y - grid(k)%yc  ; c(2) = grid(k)%xc - nodes(p3)%x
    b(3) = grid(k)%yc - nodes(p2)%y  ; c(3) = nodes(p2)%x - grid(k)%xc
    gradu%x = (b(1)*u(k) + b(2)*uvtx(p2) + b(3)*uvtx(p3))/area
    gradu%y = (c(1)*u(k) + c(2)*uvtx(p2) + c(3)*uvtx(p3))/area
    gradv%x = (b(1)*v(k) + b(2)*vvtx(p2) + b(3)*vvtx(p3))/area
    gradv%y = (c(1)*v(k) + c(2)*vvtx(p2) + c(3)*vvtx(p3))/area

    ! Compute the turbulent eddy viscosity.  This takes the average of the two
    ! end nodes, but only considers those that are wet.
    evisc = zero
    IF (parabev) THEN
      m = 0
      IF (hvtx(p2) > h_dry) THEN
        m = m + 1
        evisc = evisc + ev(p2)
      END IF
      IF (hvtx(p3) > h_dry) THEN
        m = m + 1
        evisc = evisc + ev(p3)
      END IF
      evisc = evisc/REAL(m,mp)
    END IF
    evisc = evisc + visc  ! Add the kinematic viscosity.

    ! Compute the viscous flux for the edge.
    nx = ework%normal(1)*ework%length
    ny = ework%normal(2)*ework%length
    fv(2) = evisc*hedge*gradu%x*nx
    fv(3) = evisc*hedge*gradv%x*nx
    gv(2) = evisc*hedge*gradu%y*ny
    gv(3) = evisc*hedge*gradv%y*ny

    ! Add the local flux to the global array containing the residual.
    DO l = 2,3
      phi(k,l) = phi(k,l) + fv(l) + gv(l)
    END DO

  END DO

! Now do the flow edges.  Each edge has two associated triangles (see fig. 4 of
! Jawahar and Kamath (2000)).  The shape functions are set-up and used in each
! triangle, and then the area-weighted averaged is used to compute the velocity
! gradient at the edge (similar to eqs. (17) and (18) of Jawahar and Kamath
! (2000)).
  DO i = 1,flow_edges1
    j = flowedg1(i)  ! Edge.
    CALL edge_copy(edges(j),ework)

    k1 = edges(j)%e(1)  ! First triangle, with positive contribution.
    l = local_edge(grid(k1),j)
    k2 = edges(j)%e(2)  ! Second triangle, with negative contribution.

    !hedge = half*(hvtx(p1) + hvtx(p2))
    !IF (hedge < h_dry) CYCLE  ! Dry edge.

    ! ALGORITHM I1: algorithm to deal with the interface between dry and wet
    ! cells.  A similar procedure must be implemented for the inviscid flux.
    !IF (drycell(k1) .OR. drycell(k2)) THEN
    !  z0 = half*(zvtx(edges(j)%p(1)) + zvtx(edges(j)%p(2)))
    !  h0 = hedge + z0  ! Stage at the edge.
    !  ! If the stage at he edge is lower than the bed elevation at the center
    !  ! of the dry element, set that edge flux to zero.
    !  IF (drycell(k1)) THEN
    !    IF (z(k1) + h_wet > h0) CYCLE
    !  END IF
    !  IF (drycell(k2)) THEN
    !    IF (z(k2) + h_wet > h0) CYCLE
    !  END IF
    !END IF

    ! ALGORITHM I2: this algorithm ignores the contribution due to the dry cell
    ! when computing the gradient at the edge.
    IF (drycell(k2)) THEN  ! Use gradient from cell k1.
      !Water depth at edge center computed from the wet cell, k1.
      hedge = h(k1) + delh(k1)%x*rc(k1,l)%x + delh(k1)%y*rc(k1,l)%y
      IF (hedge < h_dry) CYCLE  ! Dry edge.
      IF (grid(k1)%edge(l) < 0) THEN  ! Make sure circulation is counterclock.
        p2 = edges(j)%p(2)
        p3 = edges(j)%p(1)
        eps = invert_edge_dirs(ework)
      ELSE
        p2 = edges(j)%p(1)
        p3 = edges(j)%p(2)
      END IF
    ! Compute flow gradients.
      area = grid(k1)%area*one_third*2.0_mp
      b(1) = nodes(p2)%y - nodes(p3)%y ; c(1) = nodes(p3)%x - nodes(p2)%x
      b(2) = nodes(p3)%y - grid(k1)%yc ; c(2) = grid(k1)%xc - nodes(p3)%x
      b(3) = grid(k1)%yc - nodes(p2)%y ; c(3) = nodes(p2)%x - grid(k1)%xc
      gradu%x = (b(1)*u(k1) + b(2)*uvtx(p2) + b(3)*uvtx(p3))/area
      gradu%y = (c(1)*u(k1) + c(2)*uvtx(p2) + c(3)*uvtx(p3))/area
      gradv%x = (b(1)*v(k1) + b(2)*vvtx(p2) + b(3)*vvtx(p3))/area
      gradv%y = (c(1)*v(k1) + c(2)*vvtx(p2) + c(3)*vvtx(p3))/area

    ELSE IF (drycell(k1)) THEN  ! Use gradient from cell k2.
      !Water depth at edge center computed from the wet cell, k2.
      m = local_edge(grid(k2),j)
      hedge = h(k2) + delh(k2)%x*rc(k2,m)%x + delh(k2)%y*rc(k2,m)%y
      IF (hedge < h_dry) CYCLE  ! Dry edge.
      IF (grid(k1)%edge(m) < 0) THEN  ! Make sure circulation is counterclock.
        p2 = edges(j)%p(2)
        p3 = edges(j)%p(1)
      ELSE
        p2 = edges(j)%p(1)
        p3 = edges(j)%p(2)
        eps = invert_edge_dirs(ework)
      END IF
    ! Compute flow gradients.
      area = grid(k2)%area*one_third*2.0_mp
      b(1) = nodes(p2)%y - nodes(p3)%y ; c(1) = nodes(p3)%x - nodes(p2)%x
      b(2) = nodes(p3)%y - grid(k2)%yc ; c(2) = grid(k2)%xc - nodes(p3)%x
      b(3) = grid(k2)%yc - nodes(p2)%y ; c(3) = nodes(p2)%x - grid(k2)%xc
      gradu%x = (b(1)*u(k2) + b(2)*uvtx(p2) + b(3)*uvtx(p3))/area
      gradu%y = (c(1)*u(k2) + c(2)*uvtx(p2) + c(3)*uvtx(p3))/area
      gradv%x = (b(1)*v(k2) + b(2)*vvtx(p2) + b(3)*vvtx(p3))/area
      gradv%y = (c(1)*v(k2) + c(2)*vvtx(p2) + c(3)*vvtx(p3))/area

    ELSE  ! Both cells wet: use area-weighted gradients from both sides.
      IF (grid(k1)%edge(l) < 0) THEN  ! Make sure circulation is counterclock.
        p1 = edges(j)%p(2)
        p2 = edges(j)%p(1)
        eps = invert_edge_dirs(ework)
      ELSE
        p1 = edges(j)%p(1)
        p2 = edges(j)%p(2)
      END IF
      hedge = half*(hvtx(p1) + hvtx(p2))
      ! Compute flow gradients.
      area = (grid(k1)%area + grid(k2)%area)*one_third*2.0_mp
      x3 = grid(k1)%xc  ;  y3 = grid(k1)%yc
      xc = grid(k2)%xc  ;  yc = grid(k2)%yc
      gradu%x = ((yc - y3)*(uvtx(p1) - uvtx(p2)) + (nodes(p1)%y - &
                nodes(p2)%y)*(u(k1) - u(k2)))/area
      gradu%y = ((x3 - xc)*(uvtx(p1) - uvtx(p2)) + (nodes(p2)%x - &
                nodes(p1)%x)*(u(k1) - u(k2)))/area
      gradv%x = ((yc - y3)*(vvtx(p1) - vvtx(p2)) + (nodes(p1)%y - &
                nodes(p2)%y)*(v(k1) - v(k2)))/area
      gradv%y = ((x3 - xc)*(vvtx(p1) - vvtx(p2)) + (nodes(p2)%x - &
                nodes(p1)%x)*(v(k1) - v(k2)))/area
    END IF

    ! Compute the turbulent eddy viscosity.  This takes the average of the two
    ! end nodes, but only considers those that are wet.
    evisc = zero
    IF (parabev) THEN
      m = 0
      IF (hvtx(p1) > h_dry) THEN
        m = m + 1
        evisc = evisc + ev(p1)
      END IF
      IF (hvtx(p2) > h_dry) THEN
        m = m + 1
        evisc = evisc + ev(p2)
      END IF
      evisc = evisc/REAL(m,mp)
    END IF
    evisc = evisc + visc  ! Add the kinematic viscosity.

    ! Compute the viscous flux for the edge.
    nx = ework%normal(1)*ework%length
    ny = ework%normal(2)*ework%length
    fv(2) = evisc*hedge*gradu%x*nx
    fv(3) = evisc*hedge*gradv%x*nx
    gv(2) = evisc*hedge*gradu%y*ny
    gv(3) = evisc*hedge*gradv%y*ny

    ! Add the local flux to the global arrays containing the residual.
    DO l = 2,3
      phi(k1,l) = phi(k1,l) + fv(l) + gv(l)
      phi(k2,l) = phi(k2,l) - fv(l) - gv(l)
    END DO

  END DO

END SUBROUTINE visc_fluxes
