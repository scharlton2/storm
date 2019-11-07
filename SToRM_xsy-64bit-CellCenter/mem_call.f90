INTEGER FUNCTION mem_call()
  USE memory
  IMPLICIT none

!-----------------------------------------------------------------------------!
!                                                                             !
!  Query function that returns the number of bytes ALLOCATEd by SToRM.        !
!                                                                             !
!  F. Simoes, July 2007                                                       !
!  Last updated (mm-dd-yyyy): 07-16-2007 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

  mem_call = mem_used

END FUNCTION mem_call
