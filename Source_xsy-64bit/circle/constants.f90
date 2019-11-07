!-----------------------------------------------------------------------------!
!                                                                             !
!  This file contains the physical and numerical constants used in SToRM.     !
!                                                                             !
!  F. Simoes, December 2003                                                   !
!  Last updated (mm-dd-yyyy): 03-03-2004                                      !
!                                                                             !
!-----------------------------------------------------------------------------!

MODULE constants
  USE parameters
  IMPLICIT NONE
  SAVE

! machine precision.
  REAL(KIND=mp), PARAMETER :: fmachp = EPSILON(1.0_mp)
  REAL(KIND=mp), PARAMETER :: vsmall = TINY(1.0_mp)

! Numerical constants.
  REAL(KIND=mp), PARAMETER :: half = 0.50_mp
  REAL(KIND=mp), PARAMETER :: one = 1.0_mp
  REAL(KIND=mp), PARAMETER :: one_third = one/3.0_mp
  REAL(KIND=mp), PARAMETER :: zero = 0.0_mp

! Physical constants.
  REAL(KIND=mp), PARAMETER :: g = 9.81_mp  ! Acceleration due to gravity.

END MODULE constants
