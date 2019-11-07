FUNCTION skin_friction(ftype,ffactor,depth,uvel,vvel,nu)
  USE parameters
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This function returns the skin friction coefficient for a given set of     !
!  hydraulic conditions.  This is the usual drag coefficient.                 !
!                                                                             !
!  INPUT:                                                                     !
!    ftype          type of friction factor: =0 for Manning's n, =1 for       !
!                   Chezy's, 2 for the usual drag coefficient (in which case  !
!                   the output is equal to the inpur), =3 Colebrook law, =4   !
!                   roughness height based on Strickler's law;                !
!    ffactor        value of the friction factor;                             !
!    depth          water depth;                                              !
!    uvel           flow velocity in the x-axis direction;                    !
!    vvel           flow velocity in the y-axis direction;                    !
!    nu             water viscosity (needed for ftype = 3).                   !
!                                                                             !
!  OUTPUT:                                                                    !
!    skin_friction  ffactor converted to drag coefficient.  Note that         !
!                   tau_w/rho = -skin_friction*SQRT(uvel + vvel).             !
!                                                                             !
!  F. Simoes, 9 August 2005                                                   !
!  Last updated (mm-dd-yyyy): 08-14-2007 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  REAL (KIND=mp), INTENT(IN) :: depth,ffactor,nu,uvel,vvel
  INTEGER, INTENT(IN) :: ftype
  REAL (KIND=mp) :: skin_friction

! Local variables.
  REAL (KIND=mp) :: ln,re,work
  REAL (KIND=mp), PARAMETER :: small = 1.0E-9  ! A very small number.

  SELECT CASE (ftype)

  CASE (0)  ! Manning's n used (it's also the default).
    IF (depth < small) THEN
      skin_friction = one/small
    ELSE
      skin_friction = g*ffactor*ffactor/depth**one_third
    END IF

  CASE (1)  ! Chezy's coefficient.
    IF (ffactor < small) THEN
      skin_friction = one/small
    ELSE
      skin_friction = g/(ffactor*ffactor)
    END IF

  CASE (2)  ! Usual drag coefficient.
    skin_friction = ffactor

  CASE (3)  ! Bed roughness height plus Colebrook law.  See pg. 294 of Bradford
            ! and Sanders (2002).
    re = SQRT(uvel*uvel + vvel*vvel)*depth/nu  ! Reynolds number.
    IF (depth < small) THEN
      skin_friction = one/small
    ELSE
      re = MAX(re,1000.0_mp)  ! Limit Re below 1,000.
      re = MIN(re,REAL(2.5e7,mp))  ! Limit Re above 2.5 x 10^7.
      work = ffactor
      IF (ffactor/depth > 0.20_mp) work = 0.20_mp*depth
      ln = LOG(1.725_mp/re + (work/(14.8_mp*depth))**1.11)
      skin_friction = 0.204_mp/(ln*ln)
    END IF

  CASE (4)  ! Bed roughness height plus Strickler's equation.
    ! First convert ks to Manning's n (ks must be in cm), see dissertation
    ! (chap. 5, pg. 76.).
    work = ((100.0_mp*ffactor)**0.16666666666666667_mp)/51.79_mp
    ! Now convert Manning's n to drag coefficient.
    IF (depth < small) THEN
      skin_friction = one/small
    ELSE
      skin_friction = g*work*work/depth**one_third
    END IF

  CASE DEFAULT
    PRINT *,''
    PRINT *,'ERROR: invalid value in ftype, function skin_friction.'
    STOP
  END SELECT

END FUNCTION skin_friction
