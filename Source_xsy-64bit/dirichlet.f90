SUBROUTINE dirichlet(i,rkstep)
  USE geometry
  USE constants
  USE dep_vars
  USE RKparams
  USE options
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Dirichelet boundary conditions near dry cells, as suggested in Titov and   !
!  Synolakis (1998).                                                          !
!                                                                             !
!  INPUT:                                                                     !
!    i          node number where the boundary condition will be applied;     !
!    rkstep     Runge-Kutta step.                                             !
!                                                                             !
!  Francisco Simoes, February 2009                                            !
!  Last updated (mm-dd-yyyy): 02-18-2009 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  INTEGER, INTENT(IN) :: i,rkstep

! Local variables.
  INTEGER :: k,l,total
  REAL (KIND=mp) :: utotal,vtotal,ztotal

  total = 0
  utotal = zero
  vtotal = zero
  ztotal = zero
  DO k = 1,3
    l = t2t3(k,i)
    IF (l < 1) CYCLE  ! Edge triangle.
    IF (.NOT. drycell(l)) THEN
      total = total + 1
      utotal = utotal + u(l)
      vtotal = vtotal + v(l)
      ztotal = ztotal + zeta(l)
    END IF
  END DO
  IF (total > 0) THEN
    zeta(i) = ztotal/REAL(total,mp)
    h(i) = zeta(i) - z(i)
    IF (h(i) > h_dry) THEN
      u(i) = utotal/REAL(total,mp)
      v(i) = vtotal/REAL(total,mp)
      rku(i,rkstep) = h(i)*u(i)
      rkv(i,rkstep) = h(i)*v(i)
      rkh(i,rkstep) = h(i)
    ELSE
      u(i) = zero
      v(i) = zero
      h(i) = zero
      zeta(i) = z(i)
      rku(i,rkstep) = zero
      rkv(i,rkstep) = zero
      rkh(i,rkstep) = zero
      drycell(i) = .TRUE.
      bdrycell(i) = .FALSE.
    END IF
  ELSE
    u(i) = zero
    v(i) = zero
    h(i) = zero
    zeta(i) = z(i)
    rku(i,rkstep) = zero
    rkv(i,rkstep) = zero
    rkh(i,rkstep) = zero
    drycell(i) = .TRUE.
    bdrycell(i) = .FALSE.
  END IF

END SUBROUTINE dirichlet
