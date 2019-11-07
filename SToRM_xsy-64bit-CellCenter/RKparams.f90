!-----------------------------------------------------------------------------!
!                                                                             !
!  This file contains the parameters used in the Runge-Kutta time stepping.   !
!                                                                             !
!  F. Simoes, April 2008                                                      !
!  Last updated (mm-dd-yyyy): 04-11-2008 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

MODULE RKparams
  USE parameters
  IMPLICIT NONE
  SAVE

! Coefficients.
  REAL (KIND=mp), ALLOCATABLE :: rkalpha(:,:),rkbeta(:)

! Temporary variables used in the stepping procedure.
  REAL (KIND=mp), ALLOCATABLE :: rkh(:,:),rku(:,:),rkv(:,:)

END MODULE RKparams
