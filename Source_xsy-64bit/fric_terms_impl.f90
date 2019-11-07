SUBROUTINE fric_terms_impl(delta_t)
  USE parameters
  USE geometry
  USE dep_vars
  USE constants
  USE options
  USE RKparams
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Compute the source terms due to the bed friction.  The input argument,     !
!  delta_t, is the physical time step, in seconds.                            !
!                                                                             !
!  NOTE: this code is based on a suggestion by Murillo et al. (2006), but it  !
!  seems not to work very well because, in some cases, spureous eddies        !
!  tend to appear where there should be none.  This seems to be aggravated    !
!  with mesh refinement, i.e., the solution seems good in a given grid, but   !
!  eddies appear at strange locations when the grid is refined.               !
!                                                                             !
!  Francisco Simoes, March 2012                                               !
!  Last updated (mm-dd-yyyy): 07-27-2012 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy argumments:
  REAL (KIND=mp), INTENT(IN) :: delta_t
! Local variables:
  INTEGER :: i
  REAL (KIND=mp) :: ln,re

  rough = zero
  SELECT CASE (opt_friction)

  CASE (0)  ! Manning's n used (it's also the default).
    WHERE (.NOT. drycell) rough = g*cd*cd/h**one_third

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
      IF (.NOT. drycell(i)) THEN
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
    DO i = 1,n_elems
      IF (.NOT. drycell(i)) THEN
        ! First convert ks to Manning's n (ks must be in cm), see dissertation
        ! (chap. 5, pg. 76.).
        rough(i) = ((100.0_mp*cd(i))**0.16666666666666667_mp)/51.79_mp
        ! Now convert Manning's n to drag coefficient.
        rough(i) = g*rough(i)*rough(i)/h(i)**one_third
      END IF
    END DO

  CASE (5)  ! Vegetation drag model of Tsihrintzis (2001).
    DO i = 1,n_elems
      IF (.NOT. drycell(i)) THEN
        re = u_mag(i)*h(i)/visc  ! Reynolds number.
        rough(i) = half*rcoef1(i)*re**(-rcoef2(i))
      END IF
    END DO

  CASE DEFAULT
    PRINT *,''
    PRINT *,'ERROR: invalid value in opt_friction.'
    CALL byebye('Program SToRM stopped.')
  END SELECT

! Update the dependent variables.  This follows what is suggested in, say,
! Murillo et al. (2006) [page 361, section 4, equations after eqs. (12)].
  WHERE (.NOT. drycell) rough = h/(h + delta_t*rough*u_mag)
  u = u*rough
  v = v*rough
  DO i = 1,n_elems
    IF (.NOT. drycell(i)) THEN
      rku(i,rkorder) = h(i)*u(i)
      rkv(i,rkorder) = h(i)*v(i)
    END IF
  END DO

END SUBROUTINE fric_terms_impl
