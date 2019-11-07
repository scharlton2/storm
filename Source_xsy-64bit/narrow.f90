LOGICAL FUNCTION narrow(kplus,kminus,n,r)
  USE parameters
  USE constants
  USE dep_vars
  USE geometry
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This subroutine implements the linear N scheme (narrow), which is a first  !
!  order monotone scheme.  The version implemented here is the conservative   !
!  contour-integration-based N scheme described in section 3 of Csik et al.   !
!  (2002).  The source terms are thus included in the N scheme, similarly to  !
!  what is done in eq. (39) of that paper.                                    !
!                                                                             !
!  INPUT:                                                                     !
!    kplus    the residual distribution matrix associated with the positive   !
!             eigenvalues of the linearized element matrix;                   !
!    kminus   the residual distribution matrix associated with the negative   !
!             eigenvalues of the linearized element matrix;                   !
!    n        element index.                                                  !
!                                                                             !
!  OUTPUT:                                                                    !
!    narrow   .TRUE. if scheme was successful;                                !
!             .FALSE. if N_hat has no inverse;                                !
!    r(i,j)   N-scheme residuals for node i of the element n (3 by 3 array).  !
!             In the calling program, simply add r(i,j) to                    !
!             phi(grid(n)%vertex(i),j) to update the nodal residuals.         !
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
  INTEGER :: i,j,k
  REAL(KIND=mp) :: kiwi(3),nhat(3,3),w(3),wi(3),win(3)
  LOGICAL :: flag
  LOGICAL, EXTERNAL :: invert

  r = zero

! Compute the N_hat matrix, eq. (20) of Csik et al. (2002).
  nhat = zero
  DO i = 1,3
    DO j = 1,3
      nhat(i,j) = kminus(1,i,j) + kminus(2,i,j) + kminus(3,i,j)
    END DO
  END DO
  flag = invert(nhat)
  IF (.NOT. flag) THEN  ! Matrix N_hat does not exist: abort and return with
    narrow = .FALSE.    ! the error flag set.
    RETURN
  END IF

! Compute Win as in eq. (27) of Csik et al. (2002).
  kiwi = zero  !                               +
  DO i = 1,3   ! Compute summation of product K W
               !                               i i
    k = grid(n)%vertex(i) ! k is the global node number of node i in element n.
    wi(1) = h(k)
    wi(2) = u(k)
    wi(3) = v(k)
    DO j = 1,3
      kiwi(1) = kiwi(1) + kplus(i,1,j)*wi(j)
      kiwi(2) = kiwi(2) + kplus(i,2,j)*wi(j)
      kiwi(3) = kiwi(3) + kplus(i,3,j)*wi(j)
    END DO
  END DO
  w(1) = kiwi(1) + residual(n,1)
  w(2) = kiwi(2) + residual(n,2)
  w(3) = kiwi(3) + residual(n,3)

  win = zero
  DO j = 1,3
    win(1) = win(1) - nhat(1,j)*w(j)
    win(2) = win(2) - nhat(2,j)*w(j)
    win(3) = win(3) - nhat(3,j)*w(j)
  END DO

! Finally, compute the nodal contribution calculated by the N residual
! distribution scheme for element n and add it to the contributions previously
! calculated from other elements.
  DO i = 1,3
    k = grid(n)%vertex(i)
    w(1) = h(k) - win(1)  ! W = Wi - Win
    w(2) = u(k) - win(2)
    w(3) = v(k) - win(3)
    DO j = 1,3
      r(i,1) = r(i,1) - kplus(i,1,j)*w(j)  ! Minus sign for consistency...
      r(i,2) = r(i,2) - kplus(i,2,j)*w(j)
      r(i,3) = r(i,3) - kplus(i,3,j)*w(j)
    END DO
  END DO

  narrow = .TRUE.

END FUNCTION narrow
