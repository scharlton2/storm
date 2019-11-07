SUBROUTINE add_array_ptoH(i,n,iarray)
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Adds an element to the beginning (HEAD) of array "iarray".  The first      !
!  entry of the array is initialized with the value "i" and all other         !
!  elements are shifted up by 1.                                              !
!                                                                             !
!  INPUT:                                                                     !
!    i           value to be added to the head of the array;                  !
!    n           dimension of the array;                                      !
!    iarray      array to be processed (array of integers).                   !
!                                                                             !
!  OUTPUT:                                                                    !
!    n           new dimension of array (= n + 1);                            !
!    iarray      array with one more element located at position 1, i.e.,     !
!                iarray(1) = i.                                               !
!                                                                             !
!  F. Simoes, 1 May 2013                                                      !
!  Last updated (mm-dd-yyyy): 05-08-2013 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  INTEGER, INTENT(IN) :: i
  INTEGER, INTENT(INOUT) :: n
  INTEGER, DIMENSION(n+1), INTENT(INOUT) :: iarray

! Local variables.
  INTEGER :: j

  DO j = n+1,2,-1
    iarray(j) = iarray(j-1)
  END DO
  iarray(1) = i
  n = n + 1

END SUBROUTINE add_array_ptoH
