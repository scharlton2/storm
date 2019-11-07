SUBROUTINE fric_terms
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
!  Francisco Simoes, July 2007                                                !
!  Last updated (mm-dd-yyyy): 03-22-2011 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables:
  INTEGER :: i
  REAL (KIND=mp) :: ln,re

  SELECT CASE (opt_friction)

  CASE (0)  ! Manning's n used (it's also the default).
    DO i = 1,n_elems
      IF (drycell(i)) THEN
        rough(i) = zero
      ELSE
        rough(i) = g*cd(i)*cd(i)/h(i)**one_third
      END IF
    END DO

  CASE (1)  ! Chezy's coefficient.
    DO i = 1,n_elems
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
    DO i = 1,n_elems
      IF (drycell(i)) THEN
        rough(i) = zero
      ELSE
        re = u_mag(i)*h(i)/visc  ! Reynolds number.
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
    DO i = 1,n_elems
      IF (drycell(i)) THEN
        rough(i) = zero
      ELSE
        rough(i) = g*rough(i)*rough(i)/h(i)**one_third
      END IF
    END DO

  CASE (5)  ! Vegetation drag model of Tsihrintzis (2001).
    DO i = 1,n_elems
      IF (drycell(i)) THEN
        rough(i) = zero
      ELSE
        re = u_mag(i)*h(i)/visc  ! Reynolds number.
        rough(i) = half*rcoef1(i)*re**(-rcoef2(i))
      END IF
    END DO

  CASE DEFAULT
    PRINT *,''
    PRINT *,'ERROR: invalid value in opt_friction.'
    CALL byebye('Program SToRM stopped.')
  END SELECT

! Compute the friction term, element by element.  Note the negative sign.
  DO i = 1,n_elems
    rough(i) = -rough(i)*u_mag(i)
  END DO

  DO i = 1,n_elems
    phi(i,2) = phi(i,2) + rough(i)*u(i)*grid(i)%area  ! x-momentum equation.
    phi(i,3) = phi(i,3) + rough(i)*v(i)*grid(i)%area  ! y-momentum equation.
  END DO

END SUBROUTINE fric_terms
