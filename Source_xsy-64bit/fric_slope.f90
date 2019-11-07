SUBROUTINE fric_slope
  USE parameters
  USE geometry
  USE dep_vars
  USE constants
  USE options
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Compute the source terms due to the bed friction.                          !
!                                                                             !
!  Francisco Simoes, November 2005                                            !
!  Last updated (mm-dd-yyyy): 03-22-2011 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables:
  INTEGER :: i
  REAL(KIND=mp) :: f(3),ln,re
  TYPE(triangle) :: cv

  SELECT CASE (opt_friction)

  CASE (0)  ! Manning's n used (it's also the default).
    DO i = 1,n_pts
      IF (h(i) < h_dry) THEN
        rough(i) = zero
      ELSE
        rough(i) = g*cd(i)*cd(i)/h(i)**one_third
      END IF
    END DO

  CASE (1)  ! Chezy's coefficient.
    DO i = 1,n_pts
      IF (cd(i) < fmachp) THEN
        rough(i) = vlarge
      ELSE
        rough(i) = g/(cd(i)*cd(i))
      END IF
    END DO

  CASE (2)  ! Usual drag coefficient.
    rough = cd

  CASE (3)  ! Bed roughness height plus Colebrook law.  See pg. 294 of Bradford
            ! and Sanders (2002).
    DO i = 1,n_pts
      IF (h(i) < h_dry) THEN
        rough(i) = zero
      ELSE
        re = SQRT(u(i)*u(i) + v(i)*v(i))*h(i)/visc  ! Reynolds number.
        re = MAX(re,1000.0_mp)  ! Limit Re below 1,000.
        re = MIN(re,REAL(2.5e7,mp))  ! Limit Re above 2.5 x 10^7.
        rough(i) = cd(i)
        IF (cd(i)/h(i) > 0.20_mp) rough(i) = 0.20_mp*h(i)
        ln = LOG(1.725_mp/re + (rough(i)/(14.8_mp*h(i)))**1.11)
        rough(i) = 0.204_mp/(ln*ln)
      END IF
    END DO

  CASE (4)  ! Bed roughness height plus Strickler's equation.
    ! First convert ks to Manning's n (ks must be in cm), see dissertation
    ! (chap. 5, pg. 76.).
    rough = ((100.0_mp*cd)**0.16666666666666667_mp)/51.79_mp
    ! Now convert Manning's n to drag coefficient.
    DO i = 1,n_pts
      IF (h(i) < h_dry) THEN
        rough(i) = zero
      ELSE
        rough(i) = g*rough(i)*rough(i)/h(i)**one_third
      END IF
    END DO

  CASE (5)  ! Vegetation drag model of Tsihrintzis (2001).
    DO i = 1,n_pts
      IF (drycell(i)) THEN
        rough(i) = zero
      ELSE
        re = SQRT(u(i)*u(i) + v(i)*v(i))*h(i)/visc  ! Reynolds number.
        rough(i) = half*rcoef1(i)*re**(-rcoef2(i))
      END IF
    END DO

  CASE DEFAULT
    PRINT *,''
    PRINT *,'ERROR: invalid value in opt_friction.'
    CALL byebye('Program SToRM stopped.')
  END SELECT

! Compute the friction term, node by node.
  DO i = 1,n_pts
    rough(i) = -rough(i)*SQRT(u(i)*u(i) + v(i)*v(i))
  END DO

  DO i = 1,n_elems

    ! Skip dry or stagnant cells.
    IF (u_avg(i,1) < h_dry) CYCLE
    IF ((u_avg(i,2)*u_avg(i,2) + u_avg(i,3)*u_avg(i,3)) < u_stag*u_stag) CYCLE

    cv = grid(i)  ! Control volume.

    ! Set-up integration variables.  This part could be made a little bit more
    ! efficient, computationally, by the addition of an array to store some of
    ! these terms.
    source(i,1) = source(i,1) + zero  ! Continuity.
    ! x-momentum:
    f(1) = rough(cv%vertex(1))*u(cv%vertex(1))
    f(2) = rough(cv%vertex(2))*u(cv%vertex(2))
    f(3) = rough(cv%vertex(3))*u(cv%vertex(3))
    source(i,2) = source(i,2) + cv%area*(f(1) + f(2) + f(3))*one_third
    ! y-momentum:
    f(1) = rough(cv%vertex(1))*v(cv%vertex(1))
    f(2) = rough(cv%vertex(2))*v(cv%vertex(2))
    f(3) = rough(cv%vertex(3))*v(cv%vertex(3))
    source(i,3) = source(i,3) + cv%area*(f(1) + f(2) + f(3))*one_third
  END DO

END SUBROUTINE fric_slope
