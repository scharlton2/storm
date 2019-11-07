SUBROUTINE mem_add(m)
  USE memory
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Updates SToRM's memory usage. This subroutine must be called after each    !
!  ALLOCATE statement and it's argument must contain the number of bytes      !
!  ALLOCATEd.                                                                 !
!                                                                             !
!  F. Simoes, July 2007                                                       !
!  Last updated (mm-dd-yyyy): 07-16-2007 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

  INTEGER, INTENT(IN) :: m

  mem_used = mem_used + m

END SUBROUTINE mem_add
