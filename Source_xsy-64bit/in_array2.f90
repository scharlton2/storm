LOGICAL FUNCTION in_array2(i,a)
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Function in_array determines if any of the elements of a is equal to i.    !
!  This function is designed for the boundary nodes, which have a dimension   !
!  of 2.                                                                      !
!                                                                             !
!  INPUT:                                                                     !
!    i         the value to be found in the array;                            !
!    a         two-dimensional array of integers.                             !
!                                                                             !
!  OUTPUT:                                                                    !
!    in_array  .TRUE. if at least one of the elements of a is equal to i;     !
!              .FALSE. if none of the elements of a is equal to i.            !
!                                                                             !
!  NOTE: this function must be used in an INTERFACE construct in the calling  !
!  program:                                                                   !
!                                                                             !
!    INTERFACE                                                                !
!      LOGICAL FUNCTION in_array2(i,a)                                        !
!      INTEGER, INTENT(IN) :: i                                               !
!      INTEGER, DIMENSION(:,:), INTENT(IN) :: a                               !
!      END FUNCTION in_array2                                                 !
!    END INTERFACE                                                            !
!                                                                             !
!  F. Simoes, November 2010                                                   !
!  Last updated (mm-dd-yyyy): 02-28-2011 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

  ! Dummy variables:
  INTEGER, INTENT(IN) :: i
  INTEGER, DIMENSION(:,:), INTENT(IN) :: a

  ! Local variables:
  INTEGER :: j,k,m,n

  n = SIZE(a,1)
  m = SIZE(a,2)
  DO k = 1,m
    DO j = 1,n
      IF (a(j,k) == i) THEN
        in_array2 = .TRUE.
        RETURN
      END IF
    END DO
  END DO

  in_array2 = .FALSE.

END FUNCTION in_array2
