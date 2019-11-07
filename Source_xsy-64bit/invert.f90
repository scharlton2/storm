LOGICAL FUNCTION invert(a)
  USE parameters
  USE constants
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This function inverts matrix A by applying exchange steps (also called     !
!  the method of pivotization), as per the code in PIVOT.F90 of Engeln-       !
!  Mullges and Uhlig (1996).  The original code was customized for StoRM and  !
!  tested for accuracy using MATLAB.                                          !
!                                                                             !
!  INPUT:                                                                     !
!    a       the original 3 by 3 matrix to be inverted (double array).        !
!                                                                             !
!  OUTPUT:                                                                    !
!    a       the inverted matrix (overwritten into the original array);       !
!    invert  .TRUE. if the inversion was successful;                          !
!            .FALSE. if the inversion failed.                                 !
!                                                                             !
!  Francisco Simoes, October 2004                                             !
!  Last updated (mm-dd-yyyy): 10-27-2004 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

  ! Dummy variables:
  REAL(KIND=mp), INTENT(INOUT) :: a(3,3)

  ! Local variables:
  INTEGER :: i,ix,iy,j,k,m,mx(3),my(3),nx,ny
  REAL(KIND=mp) :: b(3,3),dummy,factor,h,pivo,s1,s2

! Store matrix A in array B.
  b = a

! Initialize the pivot vectors MX and MY to zero.
  mx = 0
  my = 0

! Determine the pivot element.
  DO i = 1,3
    pivo = zero

    DO ix = 1,3
      IF (mx(ix) == 0) THEN
        DO iy = 1,3
          IF (my(iy) == 0) THEN
            IF (ABS(b(ix,iy)) > ABS(pivo) ) THEN
              pivo = b(ix,iy)
              nx = ix
              ny = iy
            ENDIF
          ENDIF
        END DO
      ENDIF
    END DO

! If the pivot element is nearly zero, the matrix is numerically singular.
    IF (ABS(pivo) <= 4.0_mp*fmachp) THEN
      invert = .FALSE.
      RETURN
    ENDIF

! Saving the indices of the pivot element.
    mx(nx) = ny
    my(ny) = nx

! Calculation of the matrix elements according to the rules for an exchange
! step.
    dummy = 1.0d0 / pivo
    DO j = 1,3
      IF (j /= nx) THEN
        factor = b(j,ny)*dummy
        DO k = 1,3
          b(j,k) = b(j,k) - b(nx,k)*factor
        END DO
        b(j,ny) = -factor
      ENDIF
    END DO
    DO k = 1,3
      b(nx,k) = b(nx,k)*dummy
    END DO
    b(nx,ny) = dummy

  END DO

! Reverse row and column permutations.
  DO i = 1,2
    DO m = i,3
      IF (mx(m) == i) EXIT
    END DO
    j = m
    IF (j /= i) THEN
      DO k = 1,3
        h = b(i,k)
        b(i,k) = b(j,k)
        b(j,k) = h
      END DO
      mx(j) = mx(i)
      mx(i) = i
    ENDIF
    DO m = i,3
      IF (my(m) == i) EXIT
    END DO
    j = m
    IF (j /= i) THEN
      DO k = 1,3
        h = b(k,i)
        b(k,i) = b(k,j)
        b(k,j) = h
      END DO
      my(j) = my(i)
      my(i) = i
    ENDIF
  END DO

! Forming the difference S = A*B - I, where I is the identity matrix.
! Forming the sum S1 of the absolute values of the diagonal elements of the
! matrix A*B - I and the sum S2 of remaining elements. Theoretically, S1 and S2
! should both equal zero.
  s1 = zero
  s2 = zero
  DO i = 1,3
    DO j = 1,3
      h = zero
      DO k = 1,3
        h = h + a(i,k)*b(k,j)
      END DO
      IF (i == j) THEN
        s1 = s1 + ABS(h - one)
      ELSE
        s2 = s2 + ABS(h)
      ENDIF
    END DO
  END DO

  a = b
  invert = .TRUE.

END FUNCTION invert
