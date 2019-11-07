LOGICAL FUNCTION nonlinearb(kplus,kminus,n,r)
  USE parameters
  USE constants
  USE dep_vars
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This subroutine implements the linear nonlinear B scheme (blended), which  !
!  is a high resolution scheme that is both monotone and at least second      !
!  order accurate.  The version implemented here is described in section 2.5  !
!  of Csik et al. (2002).                                                     !
!                                                                             !
!  INPUT:                                                                     !
!    kplus      the residual distribution matrix associated with the          !
!               positive eigenvalues of the linearized element matrix;        !
!    kminus     the residual distribution matrix associated with the          !
!               negative eigenvalues of the linearized element matrix;        !
!    n          element index.                                                !
!                                                                             !
!  OUTPUT:                                                                    !
!    nonlinearb .TRUE. if scheme was successful;                              !
!               .FALSE. if the scheme failed;                                 !
!    r(i,j)     B-scheme residuals for node i of the element n (3 by 3        !
!               array). In the calling program, simply add r(i,j) to          !
!               phi(grid(n)%vertex(i),j) to update the nodal residuals.       !
!                                                                             !
!  Francisco Simoes, October 2004                                             !
!  Last updated (mm-dd-yyyy): 8-28-2006 by F. Simoes                          !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables:
  INTEGER, INTENT(IN) :: n
  REAL(KIND=mp), INTENT(IN) :: kplus(3,3,3),kminus(3,3,3)
  REAL(KIND=mp), INTENT(OUT) :: r(3,3)

! Local variables:
  INTEGER :: i,j
  REAL(KIND=mp) :: coef(3),rlda(3,3),rnarrow(3,3),sum
  LOGICAL :: flag
  LOGICAL, EXTERNAL :: equals,lda,narrow

  r = zero

  flag = narrow(kplus,kminus,n,rnarrow)
  IF (.NOT. flag) THEN  ! If the N scheme fails, use the LDA scheme.
    flag = lda(kplus,n,rlda)
    IF (.NOT. flag) THEN
      nonlinearb = .FALSE.
      RETURN
    END IF
    r = rlda
    nonlinearb = .TRUE.
    RETURN
  END IF

  flag = lda(kplus,n,rlda)
  IF (.NOT. flag) THEN  ! If the LDA scheme fails, use the N scheme.
    flag = narrow(kplus,kminus,n,rnarrow)
    IF (.NOT. flag) THEN
      nonlinearb = .FALSE.
      RETURN
    END IF
    r = rnarrow
    nonlinearb = .TRUE.
    RETURN
  END IF

! Compute the method's weighting coefficients according to eq. (26) of Csik
! (2002).
  DO i = 1,3
    sum = zero
    DO j = 1,3
      sum = sum + ABS(rnarrow(j,i))
    END DO
    IF (equals(sum,zero)) THEN
      nonlinearb = .FALSE.
      RETURN
    END IF
    coef(i) = ABS(residual(n,i))/sum
  END DO

  !PRINT *,coef(1),coef(2),coef(3)

! Finally, compute the B scheme nodal contributions.
  DO i = 1,3
    DO j = 1,3
      r(i,j) = coef(j)*rnarrow(i,j) + (one - coef(j))*rlda(i,j)
    END DO
  END DO

  nonlinearb = .TRUE.

END FUNCTION nonlinearb
