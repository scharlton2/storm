SUBROUTINE rmv_array_pto(i,n,iarray)
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Removes element i from array "iarray".                                     !
!                                                                             !
!  INPUT:                                                                     !
!    i           index of element to be removed (0 < i <= n);                 !
!    n           dimension of the array;                                      !
!    iarray      array to be processed (array of integers).                   !
!                                                                             !
!  OUTPUT:                                                                    !
!    n           new dimension of array (= n - 1);                            !
!    iarray      array with one less element.                                 !
!                the right sequence that forms the contiguous grid edges;     !
!                                                                             !
!  F. Simoes, 1 May 2013                                                      !
!  Last updated (mm-dd-yyyy): 05-08-2013 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  INTEGER, INTENT(IN) :: i
  INTEGER, INTENT(INOUT) :: n
  INTEGER, DIMENSION(n), INTENT(INOUT) :: iarray

! Local variables.
  INTEGER :: j

  DO j = i,n-1
    iarray(j) = iarray(j+1)
  END DO
  n = n - 1

END SUBROUTINE rmv_array_pto
