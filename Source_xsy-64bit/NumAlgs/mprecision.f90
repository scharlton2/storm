MODULE mprec

! Machine precision:
  INTEGER, PARAMETER :: mp = KIND(1.0D0)  ! = KIND(1.0) for single precision.

  CONTAINS

  INTEGER FUNCTION machpd(x)
    IMPLICIT NONE

    !  This subroutine was adapted from the subroutine with the same name in
    !  Engeln-Mullges and Uhlig (1996).

    REAL(KIND=mp), INTENT(IN) :: x

    machpd = 0
    IF (1.0_mp < x) machpd = 1

  END FUNCTION machpd

END MODULE mprec

PROGRAM mprecision
  USE mprec
  IMPLICIT NONE

  !  This program computes the machine precision, as per Engeln-Mullges and
  !  Uhlig (1996), lines 60-63 of file PIVOT.F90.  It was developed to find out
  !  a better way to implement machpd() in Fortran 90 for use in SToRM.

  REAL(KIND=mp) :: fmachp

  PRINT *,"Program mprecision running..."

  fmachp = 1.0_mp
5 fmachp = 0.5_mp*fmachp
  IF (machpd(1.0_mp + fmachp) == 1) GOTO 5
  fmachp = fmachp*2.0_mp
  PRINT *,"FMACHP     =",fmachp

  ! Now use some of the intrinsic Fortran functions.
  PRINT *,"EPSILON(x) =",EPSILON(1.0_mp)
  PRINT *,"TINY(x)    =",TINY(1.0_mp)

  ! After running this program, it is easy to see that fmachp and EPSILON(x)
  ! are the same.
  PRINT *,"Conclusion: ",fmachp - EPSILON(1.0_mp)

END PROGRAM mprecision
