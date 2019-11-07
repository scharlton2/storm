SUBROUTINE hdlsq_grad_array4h
  USE parameters
  USE geometry
  USE constants
  USE options
  USE dep_vars
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Carry out the least-squares gradient construction of variable phi at the,  !
!  center of a triangular control volume, following technique of Mavriplis    !
!  (2003).  This particular version is for the water depth only.  It          !
!  contains a special treatment of the dry nodes that sets the free surface   !
!  to zero slope.  The computation of the water depth gradient is consistent  !
!  (identical) to the computation of the bed elevation gradient and ensures   !
!  equilibrium at the still water limit.  The gradient of the bed must be     !
!  known.  Use for cell-centered solvers only.                                !
!                                                                             !
!  F. Simoes, September 2007                                                  !
!  Last updated (mm-dd-yyyy): 10-25-2007 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables.
  INTEGER :: i,k,n
  REAL (KIND=mp) :: di,dphi,ei,hnode
  LOGICAL :: etac

  DO i = 1,n_elems
    IF (h(i) < h_dry) THEN  ! Dry cell.
      delh(i)%x = zero
      delh(i)%y = zero
      CYCLE
    END IF

    di = zero; ei = zero
    etac = .FALSE.
    n = t2tHD(i,1)
    DO k = 1,n
      hnode = h(t2tHD(i,k+1))  ! Node of the computational molecule.
      ! A simpler solution to set a flat free surface...
      !IF (hnode < h_dry .AND. z(t2tHD(i,k+1))) > h(i) + z(i)) &
      !  hnode = z(i) + h(i) - z(t2tHD(i,k+1))
      ! ...and a more radical solution.
      IF (hnode < h_dry .AND. z(t2tHD(i,k+1)) > h(i) + z(i)) THEN
        etac = .TRUE.
        EXIT
      END IF
      dphi = hnode - h(i)
      di = di + lsqweight(i,k)*dphi*dx(i,k)
      ei = ei + lsqweight(i,k)*dphi*dy(i,k)
    END DO

    IF (etac) THEN
      delh(i)%x = -delz(i)%x  ! Sets water surface to zero slope.
      delh(i)%y = -delz(i)%y
    ELSE
      ! Solve the system using Cramer's rule.
      delh(i)%x = (ci(i)*di - bi(i)*ei)/det(i)
      delh(i)%y = (ai(i)*ei - di*bi(i))/det(i)
    END IF
  END DO

END SUBROUTINE hdlsq_grad_array4h
