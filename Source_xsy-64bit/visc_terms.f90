SUBROUTINE visc_terms
  USE geometry
  USE dep_vars
  USE options
  USE constants
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This subroutine computes the viscous fluxes by first calculating the       !
!  velocity gradients at the nodes and the using Gauss's theorem to           !
!  integrate them in the element using the mid-point rule on each edge.       !
!  More details about the form of the governing equation and this type of     !
!  term integration can be seen in sec. 2 of Anastasiou and Chan (1997).      !
!                                                                             !
!  Francisco Simoes, November 2006                                            !
!  Last updated (mm-dd-yyyy): 08-20-2007 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables.
  INTEGER :: e,i,j,p1,p2
  REAL (KIND=mp) :: fv(3),gv(3),hedge,nx,ny,uxedge,uyedge,vxedge,vyedge

  CALL lsq_grad_array(opt_solver,u,delu,n_pts)
  CALL lsq_grad_array(opt_solver,v,delv,n_pts)

  DO i = 1,n_elems
    fv = zero
    gv = zero
    DO j = 1,3
      e = grid(i)%edge(j)
      ! Find the normal to the edge.
      IF (e < 0) THEN  ! The edge points into the triangle: invert direction.
        e = ABS(e)
        nx = -edges(e)%normal(1)
        ny = -edges(e)%normal(2)
      ELSE
        nx = edges(e)%normal(1)
        ny = edges(e)%normal(2)
      END IF
      ! Compute the relevant quatities at the edge's mid-point.
      p1 = edges(e)%p(2)
      p2 = edges(e)%p(2)
      hedge = half*(h(p1) + h(p2))
      uxedge = half*(delu(p1)%x + delu(p2)%x)
      uyedge = half*(delu(p1)%y + delu(p2)%y)
      vxedge = half*(delv(p1)%x + delv(p2)%x)
      vyedge = half*(delv(p1)%y + delv(p2)%y)
      ! Compute the viscous fluxes.
      fv(2) = fv(2) + visc*hedge*uxedge*nx
      fv(3) = fv(3) + visc*hedge*vxedge*nx
      gv(2) = gv(2) + visc*hedge*uyedge*ny
      gv(3) = gv(3) + visc*hedge*vyedge*ny
    END DO

    ! Add to the main source terms.
    source(i,1) = source(i,1) + zero           ! Continuity.
    source(i,2) = source(i,2) + fv(2) + gv(2)  ! x-momentum.
    source(i,3) = source(i,3) + fv(3) + gv(3)  ! y-momentum.
  END DO

END SUBROUTINE visc_terms
