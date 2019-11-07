SUBROUTINE weslp_fluxes(hdry)
  USE parameters
  USE constants
  USE geometry
  USE dep_vars
  USE vbc_arrays
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Computes the contributions to the residual due to the inviscid wall        !
!  convective fluxes for a free slip boundary.  It uses the technique of      !
!  Anastasiou and Chan (1997) with modifications and improvements,            !
!  especially in what concerns the cell drying/wetting process.               !
!                                                                             !
!  INPUT:                                                                     !
!    hdry   water depth threshold below which the cell is considered dry.     !
!                                                                             !
!  Francisco Simoes, October 2007                                             !
!  Last updated (mm-dd-yyyy): 02-27-2008 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  REAL (KIND=mp), INTENT(IN) :: hdry

! Local variables.
  INTEGER :: i,k,l,m
  REAL (KIND=mp) :: eps,f(3),hleft,hright,uleft,uright,vleft,vright
  TYPE(edge) :: ework  ! Working edge.
  INTEGER, EXTERNAL :: edge_in_element,local_edge
  REAL (KIND=mp), EXTERNAL :: invert_edge_dirs

  DO l = 1,wall_edges1
    k = walledg1(l)
    i = edge_in_element(edges(k))
    CALL edge_copy(edges(k),ework)
    eps = one

    IF (drycell(i)) CYCLE

    ! Compute the left-hand side state variables from (element i).
    m = local_edge(grid(i),k)
    hleft = h(i) + delh(i)%x*rc(i,m)%x + delh(i)%y*rc(i,m)%y

    IF (hleft < hdry) THEN  ! Do not compute dry cells.
      IF (hleft < zero) partdry(i) = .TRUE.
      CYCLE
    ELSE
      uleft = u(i) + delu(i)%x*rc(i,m)%x + delu(i)%y*rc(i,m)%y
      vleft = v(i) + delv(i)%x*rc(i,m)%x + delv(i)%y*rc(i,m)%y
    END IF
    ! Set-up normal pointing out of domain.
    ! This means that the normal is pointing into element i: reverse it.
    IF (grid(i)%edge(m) < 0) eps = invert_edge_dirs(ework)

    ! Compute the right-hand side state boundary conditions.
    hright = hleft
    ! Temporarily, store the magnitude of the projection of the velocity
    ! along the boundary edge in variable vright.
    vright = uleft*ework%tang(1) + vleft*ework%tang(2)
    uright = vright*ework%tang(1)
    vright = vright*ework%tang(2)

    CALL flux_ac(ework,hleft,uleft,vleft,hright,uright,vright,f)

    ! Add the edge flux to the residual's array for element i.  Note that the
    ! flux of mass through the wall must be zero, therefore f(1) = 0.
    !f(1) = zero
    !phi(i,1) = phi(i,1) - f(1)*ework%length
    phi(i,2) = phi(i,2) - f(2)*ework%length
    phi(i,3) = phi(i,3) - f(3)*ework%length
  END DO

END SUBROUTINE weslp_fluxes
