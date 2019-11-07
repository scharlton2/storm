SUBROUTINE phisubt(phit)
  USE parameters
  USE constants
  USE geometry
  USE dep_vars
  USE options
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Here, the residuals, PHI_t, are computed directly from the Jacobian        !
!  matrices, K_i, linearized for the element according to eq. (3.36) of       !
!  Caraeni (2000).  This is a conservative discretization.  Unfortunately,    !
!  it makes impossible to enforce directly the condition of zero flow across  !
!  solid boundaries unless zero velocities are imposed.  Use with care and    !
!  only when zero velocities are specified at solid walls.                    !
!                                                                             !
!  Francisco Simoes, November 2005                                            !
!  Last updated (mm-dd-yyyy): 01-12-2007 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables.
  REAL(KIND=mp), INTENT(OUT) :: phit(n_elems,3)

! Local variables.
  INTEGER :: i,j,k,l,m
  REAL(KIND=mp) :: hcell,ki(3,3),nxi,nyi,ucell,vcell

  phit = zero  ! Note that this is an implicit DO-loop.

! Main DO-loop over all the elements.
  DO i = 1,n_elems

    ! Cell quantities, linearized as required by the residual distribution
    ! theory.
    hcell = u_avg(i,1)
    ucell = u_avg(i,2)
    vcell = u_avg(i,3)

    !IF (hcell < h_dry) CYCLE  ! Dry cell.
    !IF ((ucell*ucell + vcell*vcell) < u_stag*u_stag) CYCLE  ! Stagnant cell.

    ! Loop over all the element vertices.
    DO j = 1,3
      k = MAX(1,MOD(j+1,4))  ! Element index of edge opposite to node j.
      nxi = edges(ABS(grid(i)%edge(k)))%normal(1)
      nyi = edges(ABS(grid(i)%edge(k)))%normal(2)
      ! We want the normal pointing inwards, therefore a sign change may be
      ! needed:
      IF (grid(i)%edge(k) > 0) THEN
        nxi = -nxi
        nyi = -nyi
      END IF

      CALL jacobian(ki,nxi,nyi,hcell,ucell,vcell)

      ! Solve eq. (3.36) of Caraeni (2000).
      l = grid(i)%vertex(j)
      DO m = 1,3
        phit(i,m) = phit(i,m) + ki(m,1)*h(l) + ki(m,2)*h(l)*u(l) + &
          ki(m,3)*h(l)*v(l)
      END DO

    END DO

  END DO

END SUBROUTINE phisubt
