LOGICAL FUNCTION lda(kplus,n,r)
  USE parameters
  USE constants
  USE dep_vars
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This subroutine implements the linear LDA scheme (low diffusion A), which  !
!  is a second-order linear scheme.  Therefore ,the LDA scheme is a           !
!  nonmonotonic scheme.                                                       !
!                                                                             !
!  INPUT:                                                                     !
!    kplus   the residual distribution matrix associated with the positive    !
!            eigenvalues of the linearized element matrix;                    !
!    n       element index.                                                   !
!                                                                             !
!  OUTPUT:                                                                    !
!    lda     .TRUE. if scheme was successful;                                 !
!            .FALSE. if N_hat has no inverse;                                 !
!    r(i,j)  N-scheme residuals for node i of the element n (3 by 3 array).   !
!            In the calling program, simply add r(i,j) to                     !
!            phi(grid(n)%vertex(i),j) to update the nodal residuals.          !
!                                                                             !
!  Francisco Simoes, October 2004                                             !
!  Last updated (mm-dd-yyyy): 11-03-2004 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables:
  INTEGER, INTENT(IN) :: n
  REAL(KIND=mp), INTENT(IN) :: kplus(3,3,3)
  REAL(KIND=mp), INTENT(OUT) :: r(3,3)

! Local variables:
  INTEGER :: i,j,l,m
  REAL(KIND=mp) :: nhat(3,3),beta(3,3),k(3,3)
  LOGICAL :: flag
  LOGICAL, EXTERNAL :: invert

  r = zero

! Compute the N_hat matrix, eq. (20) of Csik et al. (2002).
  nhat = zero
  DO i = 1,3
    DO j = 1,3
      nhat(i,j) = kplus(1,i,j) + kplus(2,i,j) + kplus(3,i,j)
    END DO
  END DO
  flag = invert(nhat)
  IF (.NOT. flag) THEN  ! Matrix N_hat does not exist: abort and return with
    lda = .FALSE.       ! the error flag set.
    RETURN
  END IF

  DO i = 1,3
    ! Matrix k is used for MATMUL.
    DO l = 1,3
      DO m = 1,3
        k(l,m) = kplus(i,l,m)
      END DO
    END DO
    beta = MATMUL(k,nhat)  ! Eq. (24) of Csik et al. (2002).

    DO j = 1,3  ! Eq. (22) of Csik et al (2002).
      r(i,1) = r(i,1) + beta(1,j)*residual(n,j)
      r(i,2) = r(i,2) + beta(2,j)*residual(n,j)
      r(i,3) = r(i,3) + beta(3,j)*residual(n,j)
    END DO
  END DO

  lda = .TRUE.

END FUNCTION lda
