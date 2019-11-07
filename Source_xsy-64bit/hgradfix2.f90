SUBROUTINE hgradfix2(rkstep,solver)
  USE geometry
  USE constants
  USE dep_vars
  USE RKparams
  USE options
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Set a flat free surface on all the partially dry elements.  This           !
!  SUBROUTINE replaces hgradfix, which does something else...                 !
!                                                                             !
!  INPUT:                                                                     !
!    rkstep     Runge-Kutta step, may be needed to enforce boundary           !
!               conditions here.                                              !
!    solver = 2 for FVT and other cell-centered solvers;                      !
!             3 for FVZ and other cell-centered solvers that use the free     !
!               surface elevation as dependent variable instead of water      !
!               depth.  (This option has not been implemented.)               !
!                                                                             !
!  Francisco Simoes, October 2007                                             !
!  Last updated (mm-dd-yyyy): 10-06-2012 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  INTEGER, INTENT(IN) :: rkstep,solver

! Local variables.
  INTEGER :: i,k,l,total
  REAL (KIND=mp) :: htotal,utotal,vtotal,ztotal
  LOGICAL :: setgrad
  INTEGER, EXTERNAL :: edge_in_element

  bdrycell = .FALSE.

  SELECT CASE (solver)

!-----------------------------------------------------------------------------!
!                                 FVT solver                                  !
!-----------------------------------------------------------------------------!

! Set the gradient of the water depth to the symmetric of the bed gradient in
! cells that are partially dry.
  CASE (2)

    DO i = 1,n_elems
      IF (drycell(i)) CYCLE
      setgrad = zeta(i) < zvtx(csortedz(i)%p(3))

      ! Second loop, the override condition: if there are neighboring dry cells
      ! that have elevation below this element's stage, then those cells will
      ! be wetted in the next time step and the zero gradient boundary
      ! condition does not apply.
      DO k = 1,3
        l = t2t3(k,i)
        IF (l < 1) CYCLE  ! Edge triangle.
        IF (drycell(l) .AND. z(l) < zeta(i)) setgrad = .FALSE.
      END DO

      ! Set the free surface gradient to zero by setting the gradient of the
      ! water depth to the symmetric of the bed gradient.
      IF (setgrad) THEN
        delh(i)%x = -delz(i)%x
        delh(i)%y = -delz(i)%y
        bdrycell(i) = .TRUE.

        !delu(i)%x = zero; delu(i)%y = zero
        !delv(i)%x = zero; delv(i)%y = zero

        ! Use Dirichelet boundary conditions near dry cells, as suggested in
        ! Titov and Synolakis (1998).
        !CALL dirichlet(i,rkstep)

        ! Enforce zero velocity and set water depth to the average of nearby
        ! wet cells.
        !total = 0
        !htotal = zero
        !DO k = 1,3
        !  l = t2t3(k,i)
        !  IF (l < 1) CYCLE  ! Edge triangle.
        !  IF (.NOT. drycell(l)) THEN
        !    total = total + 1
        !    htotal = htotal + u(l)
        !  END IF
        !END DO
        !IF (total > 0) h(i) = htotal/REAL(total,mp)

        ! Set velocity, and possibly the gradient of velocity, to zero.
        u(i) = zero
        v(i) = zero
        rku(i,rkstep-1) = zero
        rkv(i,rkstep-1) = zero
        delu(i)%x = zero; delu(i)%y = zero
        delv(i)%x = zero; delv(i)%y = zero
      END IF
    END DO

!-----------------------------------------------------------------------------!
!                                 FVZ solver                                  !
!-----------------------------------------------------------------------------!

! Nothing implemented here.
  CASE (3)
    PRINT *,"ERROR: invalid option in hgradfix2 (solver = 3)."
    CALL byebye('Program STORM stopped.')

  CASE DEFAULT
    PRINT *,"ERROR: invalid option in hgradfix2."
    CALL byebye('Program STORM stopped.')

  END SELECT

END SUBROUTINE hgradfix2
