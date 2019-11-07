!-----------------------------------------------------------------------------!
!                                                                             !
!  This file contains the physical and numerical constants used in SToRM.     !
!                                                                             !
!  F. Simoes, December 2003                                                   !
!  Last updated (mm-dd-yyyy): 07-07-2014 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

MODULE constants
  USE parameters
  IMPLICIT NONE
  SAVE

! Machine precision.
  REAL (KIND=mp), PARAMETER :: zero = 0.0_mp
  REAL (KIND=mp), PARAMETER :: fmachp = EPSILON(zero)
  REAL (KIND=mp), PARAMETER :: vsmall = TINY(zero)
  REAL (KIND=mp), PARAMETER :: vlarge = HUGE(zero)

! Numerical constants.
  REAL (KIND=mp), PARAMETER :: half = 0.50_mp
  REAL (KIND=mp), PARAMETER :: one = 1.0_mp
  REAL (KIND=mp), PARAMETER :: one_third = one/3.0_mp
  REAL (KIND=mp), PARAMETER :: pi = 4.0_mp*ATAN(one)

! Physical constants.
  REAL (KIND=mp), PARAMETER :: g = 9.81_mp  ! Acceleration due to gravity.
  REAL (KIND=mp), PARAMETER :: rho = 998.60_mp  ! Density of water @ 18 deg C.
  REAL (KIND=mp), PARAMETER :: rhos = 2650.0_mp  ! Density of quartz.


END MODULE constants
