PROGRAM heff_sort
  IMPLICIT NONE

! Test of the logic to sort triangle vertices according to height above datum.
! F. Simoes, 02-02-2009.

! Machine precision:
  INTEGER, PARAMETER :: mp = KIND(1.0D0)  ! = KIND(1.0) for single precision.

! Computational grids:
  TYPE :: triangle
    INTEGER :: vertex(3)        ! Points to coordinates of vertices.
  END TYPE

! Array with three INTEGER components.
  TYPE :: vector3
    INTEGER :: p(3)
  END TYPE

! Local variables.
  TYPE(triangle) :: grid(1)
  REAL(KIND=mp) :: zvtx(3)
  TYPE(vector3) :: csortedz(1)
  INTEGER :: aux,i,j,k,n_elems

! Initialization of variables.
  n_elems = 1

! Arbitrary triangle, all vertex elevation different.
  !grid(1)%vertex(1) = 3
  !grid(1)%vertex(2) = 1
  !grid(1)%vertex(3) = 2
  !zvtx(1) = 1.0
  !zvtx(2) = 5.25
  !zvtx(3) = 3.0

! Pathologic case: two vertices equal.
  grid(1)%vertex(1) = 2
  grid(1)%vertex(2) = 3
  grid(1)%vertex(3) = 1
  zvtx(1) = 3.0
  zvtx(2) = 2.1
  zvtx(3) = 3.0

! SToRM code.
  DO i = 1,n_elems
    csortedz(i)%p(1) = grid(i)%vertex(1)
    csortedz(i)%p(2) = grid(i)%vertex(2)
    csortedz(i)%p(3) = grid(i)%vertex(3)
    DO j = 1,2  ! Two passes (worst case scenario).
      DO k = 1,2  ! Bubble-sort a list with 3 elements.
        IF (zvtx(csortedz(i)%p(k)) > zvtx(csortedz(i)%p(k+1))) THEN  ! Swap elements.
          aux = csortedz(i)%p(k)
          csortedz(i)%p(k) = csortedz(i)%p(k+1)
          csortedz(i)%p(k+1) = aux
        END IF
      END DO
    END DO

  END DO

! Output
  PRINT *,csortedz(1)%p(1),csortedz(1)%p(2),csortedz(1)%p(3)

END PROGRAM heff_sort
