MODULE parameters
  IMPLICIT NONE
  SAVE

! Machine precision:
  INTEGER, PARAMETER :: mp = KIND(1.0D0)  ! = KIND(1.0) for single precision.

! Numerical constants.
  REAL (KIND=mp), PARAMETER :: zero = 0.0_mp
  REAL (KIND=mp), PARAMETER :: half = 0.50_mp
  REAL (KIND=mp), PARAMETER :: one = 1.0_mp
  REAL (KIND=mp), PARAMETER :: one_third = one/3.0_mp

! Physical constants.
  REAL (KIND=mp), PARAMETER :: g = 9.81_mp  ! Acceleration due to gravity.

END MODULE parameters

PROGRAM kwave
  USE parameters
  IMPLICIT NONE

! To debug SUBROUTINE kinematicwave2d in SToRM.
! F. Simoes, 9/7/07

! Local variables.
  INTEGER :: i,rtype
  REAL (KIND=mp) :: cd,denom,f0,f1,f2,f3,fcoef,h,magu,magzb,s12,x1,x2,x3
  REAL (KIND=mp), PARAMETER :: relerr = 0.005  ! Relative error.
  REAL (KIND=mp), PARAMETER :: abserr = 0.0001  ! Absolute error.
  REAL (KIND=mp), PARAMETER :: small = 1.0E-12  ! A very small number.
  REAL (KIND=mp), PARAMETER :: visc = 1.139E-6  ! Water viscosity.

  REAL (KIND=mp), EXTERNAL :: skin_friction

  WRITE (*,'(A)')"INPUT PARAMETERS"
  WRITE (*,'("Roughness type: ",$)'); READ (*,*) rtype
  WRITE (*,'("Roughness value: ",$)'); READ (*,*) cd
  WRITE (*,'("Value of dzb/dx: ",$)'); READ (*,*) magzb
  WRITE (*,'("Water depth: ",$)'); READ (*,*) h

  IF (rtype == 3) THEN
    ! The computation of the friction factor is done using the bed roughness
    ! height and Colebrook law, which depends on flow velocity.  In this case
    ! the equations are non-linear and require an iterative solution technique.
    ! Here, the Anderson-Bjorck method is used (Anderson and Bjorck, 1973).
    x1 = zero ; x2 = 100.0_mp  ! Initial interval (guessed).
    fcoef = skin_friction(rtype,cd,h,x1,zero,visc)
    denom = fcoef + small
    f1 = x1 - SQRT(g*h*magzb/denom)
    fcoef = skin_friction(rtype,cd,h,x2,zero,visc)
    denom = fcoef + small
    f2 = x2 - SQRT(g*h*magzb/denom)

    DO
      s12 = (f2 - f1)/(x2 - x1)
      x3 = x2 - f2/s12
      fcoef = skin_friction(rtype,cd,h,x3,zero,visc)
      denom = fcoef + small
      f3 = x3 - SQRT(g*h*magzb/denom)

      IF (f3 < small) THEN  ! Frist break-off condition.
        magu = x3
        EXIT
      END IF

      IF (f3*f2 < zero) THEN
        x1 = x2 ; f1 = f2
        x2 = x3 ; f2 = f3
      ELSE
        f0 = one - f3/f2
        IF (f0 <= zero) f0 = half
        x2 = x3
        f1 = f0*f1
        f2 = f3
      END IF

      ! Test taken from page 5 of Engeln-Mullges and Uhlig (1996):
      IF (ABS(x2 - x1) <= ABS(x2)*relerr + abserr) THEN  ! Algorithm converged.
        magu = x1
        IF (ABS(f2) <= ABS(f1)) magu = x2
      END IF
    END DO

  ELSE
    ! No iteration necessary.
    fcoef = skin_friction(rtype,cd,h,zero,zero,visc)
    denom = fcoef + small  ! To avoid divide by zero.
    magu = SQRT(g*h*magzb/denom)
  END IF

  WRITE (*,'(\"Magnitude of velocity: ",E14.5)') magu

END PROGRAM kwave
