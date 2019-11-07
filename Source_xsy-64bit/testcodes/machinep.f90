PROGRAM machinep
  IMPLICIT NONE

! Test of machine precision values.
! F. Simoes, 13 August 2007

! Machine precision:
  INTEGER, PARAMETER :: mp = KIND(1.0D0)  ! = KIND(1.0) for single precision.

! Local variables.
  REAL (KIND=mp) :: eps

! Machine precision.
  REAL (KIND=mp), PARAMETER :: zero = 0.0_mp
  REAL (KIND=mp), PARAMETER :: fmachp = EPSILON(zero)
  REAL (KIND=mp), PARAMETER :: vsmall = TINY(zero)
  REAL (KIND=mp), PARAMETER :: vlarge = HUGE(zero)

  PRINT *,'Machine EPSILON:',fmachp
  PRINT *,'Very small number:',vsmall
  PRINT *,'Invert it:',1.0_mp/vsmall
  PRINT *,'Again:',1000000.0_mp/vsmall
  eps = 1000000.0_mp*fmachp  ! Scale a very small number by the numerator.
  IF (eps == 0.0_mp) eps = vsmall  ! Correct underflow to zero.
  PRINT *,'EPSILON:',eps
  PRINT *,'Yet again:',1000000.0_mp/eps  ! This is right!
  PRINT *,'Very large number:',vlarge

! 'eps' is a scaled number that corresponds to a machine precision definition
! of a 'very small number'.  For the example above, 'eps' is a very small
! number that, added to a denominator, will prevent a division by zero.
! 'eps' is based on EPSILON and is scaled by the numerator.  EPSILON makes it
! easy to select a delta for algorithms (such as root locators) that search
! until the calculation is within delta of an estimate.  If delta is too small
! (smaller than the decimal resolution of the data type), the algorithm might
! never halt. Scaling the value returned by EPSILON to the estimate results in
! a delta that ensures search termination.

! Output of this program:
! Machine EPSILON:  2.220446049250313E-016
! Very small number:  2.225073858507201E-308
! Invert it:  4.494232837155790E+307
! Again: Infinity
! EPSILON:  2.220446049250313E-010
! Yet again:  4.503599627370496E+015
! Very large number:  1.797693134862316E+308


END PROGRAM machinep
