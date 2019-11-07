LOGICAL FUNCTION in_array(i,a)
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Function in_array determines if any of the elements of a is equal to i.    !
!                                                                             !
!  INPUT:                                                                     !
!    i         the value to be found in the array;                            !
!    a         array of integers.                                             !
!                                                                             !
!  OUTPUT:                                                                    !
!    in_array  .TRUE. if at least one of the elements of a is equal to i;     !
!              .FALSE. if none of the elements of a is equal to i.            !
!                                                                             !
!  NOTE: this function must be used in an INTERFACE construct in the calling  !
!  program:                                                                   !
!                                                                             !
!    INTERFACE                                                                !
!      LOGICAL FUNCTION in_array(i,a)                                         !
!      INTEGER, INTENT(IN) :: i                                               !
!      INTEGER, DIMENSION(:), INTENT(IN) :: a                                 !
!      END FUNCTION in_array                                                  !
!    END INTERFACE                                                            !
!                                                                             !
!  F. Simoes, November 2004                                                   !
!  Last updated (mm-dd-yyyy): 09-19-2006 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

  ! Dummy variables:
  INTEGER, INTENT(IN) :: i
  INTEGER, DIMENSION(:), INTENT(IN) :: a

  ! Local variables:
  INTEGER :: j,n

  n = SIZE(a)
  DO j = 1,n
    IF (a(j) == i) THEN
      in_array = .TRUE.
      RETURN
    END IF
  END DO

  in_array = .FALSE.

END FUNCTION in_array
