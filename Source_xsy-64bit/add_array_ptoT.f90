SUBROUTINE add_array_ptoT(i,n,iarray)
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Adds an element to the end (TAIL) of array "iarray".                       !
!                                                                             !
!  INPUT:                                                                     !
!    i           value to be added to the tail of the array;                  !
!    n           dimension of the array;                                      !
!    iarray      array to be processed (array of integers).                   !
!                                                                             !
!  OUTPUT:                                                                    !
!    n           new dimension of array (= n + 1);                            !
!    iarray      array with one more element located at position n, i.e.,     !
!                iarray(n) = i.                                               !
!                                                                             !
!  F. Simoes, 1 May 2013                                                      !
!  Last updated (mm-dd-yyyy): 05-08-2013 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  INTEGER, INTENT(IN) :: i
  INTEGER, INTENT(INOUT) :: n
  INTEGER, DIMENSION(n+1), INTENT(INOUT) :: iarray

  n = n + 1
  iarray(n) = i

END SUBROUTINE add_array_ptoT
