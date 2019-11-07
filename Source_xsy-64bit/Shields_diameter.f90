SUBROUTINE Shields_diameter(d,tau,depth)
  USE parameters
  USE constants
  USE geometry
  USE options
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Computation of the Shields diameter: given the value of the bed shear      !
!  stress TAU (N/m^2), find the corresponding value of the grain diameter D   !
!  (mm) in the Shields diagram.                                               !
!                                                                             !
!  INPUT:                                                                     !
!    tau    value of the bed shear stress (N/m^2);                            !
!    depth  water depth.                                                      !
!                                                                             !
!  OUTPUT:                                                                    !
!    d      diameter of particle in incipient motion conditions (mm).         !
!                                                                             !
!  Francisco Simoes, July 2014                                                !
!  Last updated (mm-dd-yyyy): 07-07-2014 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  !REAL (KIND=mp), DIMENSION(n_pts), INTENT(IN) :: tau,depth
  REAL (KIND=mp), DIMENSION(n_pts) :: tau,depth
  REAL (KIND=mp), DIMENSION(n_pts), INTENT(OUT) :: d

! Local variables.
  INTEGER :: i,k
  REAL (KIND=mp) :: e,gg,f1,f2,f3,s12,tmax,tmin,x1,x2,x3
  REAL (KIND=mp), PARAMETER :: dmin = REAL(6.25e-5,mp), dmax = REAL(0.010,mp)
  REAL (KIND=mp), EXTERNAL :: Shields_tau
  LOGICAL, EXTERNAL :: equals

! Bounds of Shields diameters.
  tmin = Shields_tau(dmin)
  tmax = Shields_tau(dmax)

! Relative error tolerance.
  e = 0.01_mp

  d = zero
  DO i = 1,n_pts
    IF (depth(i) < h_dry) CYCLE  ! Dry vertices.
    IF (tau(i) < tmin) CYCLE  ! Immobile bed.
    IF (tau(i) > tmax) THEN  ! Everything moves.
      d(i) = dmax
      CYCLE
    END IF

    ! Anderson-Bjorck algorithm.
    f1 = tmin - tau(i); x1 = dmin
    f2 = tmax - tau(i); x2 = dmax
    k = 0
    DO
      k = k + 1
      s12 = (f2 - f1)/(x2 - x1)
      x3 = x2 - f2/s12
      f3 = Shields_tau(x3) - tau(i)

      IF (equals(f3,zero)) THEN  ! The zero has been found.
        d(i) = x3
        EXIT
      END IF

      IF (f3*f2 < zero) THEN
        x1 = x2
        f1 = f2
        x2 = x3
        f2 = f3
      ELSE
        gg = one - f3/f2
        IF (zero > gg) gg = half
        x2 = x3
        f1 = gg*f1
        f2 = f3
      END IF

      IF (ABS(x2 - x1) < e*ABS(x2)) THEN  ! Check for error tolerance.
        IF (ABS(f1) > ABS(f2)) THEN
          d(i) = x2
        ELSE
          d(i) = x1
        END IF
        EXIT
      END IF

      IF (k > 1000) THEN  ! Safety exit.
        d(i) = REAL(k,mp)
        EXIT
      END IF
    END DO  ! End Anderson-Bjork algorithm.
  END DO

  d = d*1000.0_mp  ! Convert from m to mm.

END SUBROUTINE Shields_diameter
