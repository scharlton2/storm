SUBROUTINE lsq_grad_array4h(solver)
  USE parameters
  USE geometry
  USE constants
  USE options
  USE dep_vars
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Carry out the least-squares gradient construction of variable phi,         !
!  following technique of Mavriplis (2003).  This particular version is for   !
!  the water depth only.  It contains a special treatment of the dry nodes    !
!  that sets the free surface to zero slope.  The computation of the water    !
!  depth gradient is consistent (identical) to the computation of the bed     !
!  elevation gradient and ensures equilibrium at the still water limit.  The  !
!  gradient of the bed must be known for solver = 2 before this subroutine    !
!  can be used.                                                               !
!                                                                             !
!  INPUT:                                                                     !
!    solver = 1 for RDS and other vertex-oriented solvers;                    !
!             2 for FVT and other cell-centered solvers;                      !
!             3 for FVZ and other cell-centered solvers that use the free     !
!               surface elevation as dependent variable instead of water      !
!               depth.                                                        !
!                                                                             !
!  Francisco Simoes, March 2007                                               !
!  Last updated (mm-dd-yyyy): 01-09-2008 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  INTEGER, INTENT(IN) :: solver

! Local variables.
  INTEGER :: i,k,n
  REAL (KIND=mp) :: di,dphi,ei,hnode
  LOGICAL :: etac

  SELECT CASE (solver)

!-----------------------------------------------------------------------------!
!                                 RDS solver                                  !
!-----------------------------------------------------------------------------!

! To do: use the bed slope to compute the water depth gradient for dry nodes.
! See below, in solver 2, how to do this.
  CASE (1)
    DO i = 1,n_pts
      di = zero; ei = zero
      n = n2n(i,1)
      DO k = 1,n
        hnode = h(n2n(i,k+1))  ! Node of the computational molecule.
        IF (hnode < h_dry) hnode = z(i) + h(i) - z(n2n(i,k+1))
        dphi = hnode - h(i)
        di = di + lsqweight(i,k)*dphi*dx(i,k)
        ei = ei + lsqweight(i,k)*dphi*dy(i,k)
      END DO

      ! Solve the system using Cramer's rule.
      delh(i)%x = (ci(i)*di - bi(i)*ei)/det(i)
      delh(i)%y = (ai(i)*ei - di*bi(i))/det(i)
    END DO

!-----------------------------------------------------------------------------!
!                                 FVT solver                                  !
!-----------------------------------------------------------------------------!

  CASE (2)
    DO i = 1,n_elems
      IF (drycell(i)) THEN  ! Dry cell.
        delh(i)%x = zero
        delh(i)%y = zero
        CYCLE
      END IF

      di = zero; ei = zero
      etac = .FALSE.  ! Flag indicates if grad eta is zero.
      n = t2t(i,1)
      DO k = 1,n
        hnode = h(t2t(i,k+1))  ! Node of the computational molecule.
        ! A simpler solution to set a flat free surface...
        !IF (hnode < h_dry .AND. z(t2t(i,k+1)) > h(i) + z(i)) &
        !  hnode = z(i) + h(i) - z(t2t(i,k+1))
        ! ...and a more radical solution.
        IF (drycell(t2t(i,k+1)) .AND. z(t2t(i,k+1)) > h(i) + z(i)) THEN
          etac = .TRUE.
          EXIT
        END IF
        dphi = hnode - h(i)
        di = di + lsqweight(i,k)*dphi*dx(i,k)
        ei = ei + lsqweight(i,k)*dphi*dy(i,k)
      END DO

      IF (etac) THEN
        delh(i)%x = -delz(i)%x  ! Sets water surface elevation to zero slope.
        delh(i)%y = -delz(i)%y
      ELSE
        ! Solve the system using Cramer's rule.
        delh(i)%x = (ci(i)*di - bi(i)*ei)/det(i)
        delh(i)%y = (ai(i)*ei - di*bi(i))/det(i)
      END IF
    END DO

!-----------------------------------------------------------------------------!
!                                 FVZ solver                                  !
!-----------------------------------------------------------------------------!

! In the FVZ solver one works with the free surface elevation instead of the
! water depth.  The gradient 'delh' is now the gradient of the free surface
! elevation.
  CASE (3)
    DO i = 1,n_elems
      IF (drycell(i)) THEN  ! Dry cell.
        delh(i)%x = delz(i)%x
        delh(i)%y = delz(i)%y
        CYCLE
      END IF

      di = zero; ei = zero
      n = t2t(i,1)
      DO k = 1,n
        ! t2t(i,k+1) is the node of the computational molecule.
        IF (drycell(t2t(i,k+1))) THEN
          di = zero  ! Sets the free surface gradient to zero.
          ei = zero
          EXIT
        END IF
        dphi = zeta(t2t(i,k+1)) - zeta(i)
        di = di + lsqweight(i,k)*dphi*dx(i,k)
        ei = ei + lsqweight(i,k)*dphi*dy(i,k)
      END DO

      ! Solve the system using Cramer's rule.
      delh(i)%x = (ci(i)*di - bi(i)*ei)/det(i)
      delh(i)%y = (ai(i)*ei - di*bi(i))/det(i)
    END DO

  CASE DEFAULT
    PRINT *,"ERROR: least-squares gradient computation, invalid option."
    CALL byebye('Program STORM stopped.')

  END SELECT

END SUBROUTINE lsq_grad_array4h

! NOTE: I think that setting the boundary condition gradients here is
! inefficient.  It's best to compute the gradients of h (or of the stage) in
! a post-processing stage, such as it is done in subroutine 'hgradfix'.
