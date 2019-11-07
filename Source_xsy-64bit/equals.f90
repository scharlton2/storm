!PURE FUNCTION equals(a,b) ! This line is Fortran 95.
FUNCTION equals(a,b) ! This line is Fortran 90.
  USE parameters  ! Defines KIND mp.
  USE constants   ! Defines zero.
  IMPLICIT NONE

! A reasonably safe way to compare of two REALs of kind mp...
! F. Simoes, October 2004.

  LOGICAL :: equals
  REAL(KIND=mp), INTENT(IN) :: a,b
  REAL(KIND=mp) :: eps

  eps = ABS(a)*fmachp   ! Scale epsilon.
  IF (eps == zero) THEN
    eps = vsmall        ! If eps underflowed to 0 use a very small positive
                        ! value for epsilon.
  END IF

  IF (ABS(a-b) > eps) THEN
    equals = .FALSE.    ! Not equal if difference > eps.
  ELSE
    equals = .TRUE.     ! Equal otherwise.
  ENDIF

END FUNCTION equals
