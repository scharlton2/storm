FUNCTION det3by3(a)
  USE parameters
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Computes the determinant of a 3 x 3 matrix.                                !
!                                                                             !
!  INPUT:                                                                     !
!    a          3 x 3 matrix.                                                 !
!                                                                             !
!  OUTPUT:                                                                    !
!    det3       the determinant of matrix a.                                  !
!                                                                             !
!  F. Simoes, May 2005                                                        !
!  Last updated (mm-dd-yyyy): 05-06-2005                                      !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables:
  REAL(KIND=mp), INTENT(IN) :: a(3,3)
  REAL(KIND=mp) :: det3by3

  det3by3 = a(1,1)*a(2,2)*a(3,3) + a(1,3)*a(2,1)*a(3,2) + a(1,2)*a(2,3)*a(3,1)&
          - a(1,3)*a(2,2)*a(3,1) - a(2,3)*a(3,2)*a(1,1) - a(1,2)*a(2,1)*a(3,3)

END FUNCTION det3by3
