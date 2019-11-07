SUBROUTINE bed_slope_DFB
  USE parameters
  USE constants
  USE geometry
  USE dep_vars
  USE options
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Compute the source terms due to the bed slope using the "divergence form   !
!  for bed slope source term (DFB)" of Valiani and Begnudelli (2006).  The    !
!  implementation follows eq. (27) modified for a triangle.                   !
!                                                                             !
!  Francisco Simoes, July 2007                                                !
!  Last updated (mm-dd-yyyy): 11-02-2007 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables.
  REAL (KIND=mp) :: eta,hr,ksign,s0(3)
  REAL (KIND=mp), PARAMETER :: three = 3.0_mp
  INTEGER :: i,j,k,p1,p2

  SELECT CASE (opt_solver)

!-----------------------------------------------------------------------------!
!                                 FVT solver                                  !
!-----------------------------------------------------------------------------!

  CASE (2)
    DO i = 1,n_elems
      IF (drycell(i)) CYCLE

      eta = z(i) + h(i)
      s0 = zero
      DO j = 1,3
        k = grid(i)%edge(j)
        ksign = SIGN(one,REAL(k,mp))  ! = -1 if normal points into the element.
        k = ABS(k)
        p1 = edges(k)%p(1)
        p2 = edges(k)%p(2)
        hr = eta - half*(zvtx(p1) + zvtx(p2))  ! Depth at center of edge k.
        IF (hr < h_dry) CYCLE
        !IF (hr < zero) THEN
        !  hr = zero
        !  partdry = .TRUE.
        !  CYCLE
        !END IF
        s0(2) = s0(2) + ksign*hr*hr*edges(k)%normal(1)*edges(k)%length
        s0(3) = s0(3) + ksign*hr*hr*edges(k)%normal(2)*edges(k)%length
      END DO

      phi(i,2) = phi(i,2) + half*g*s0(2)  ! X-momentum equation.
      phi(i,3) = phi(i,3) + half*g*s0(3)  ! Y-momentum equation.
    END DO

!-----------------------------------------------------------------------------!
!                                 FVZ solver                                  !
!-----------------------------------------------------------------------------!

  CASE (3)
    DO i = 1,n_elems
      IF (drycell(i)) CYCLE

      s0 = zero
      DO j = 1,3
        k = grid(i)%edge(j)
        ksign = SIGN(one,REAL(k,mp))  ! = -1 if normal points into the element.
        k = ABS(k)
        p1 = edges(k)%p(1)
        p2 = edges(k)%p(2)
        hr = zeta(i) - half*(zvtx(p1) + zvtx(p2))  ! Depth at center of edge k.
        IF (hr < h_dry) CYCLE
        !IF (hr < zero) THEN
        !  hr = zero
        !  partdry = .TRUE.
        !  CYCLE
        !END IF
        s0(2) = s0(2) + ksign*hr*hr*edges(k)%normal(1)*edges(k)%length
        s0(3) = s0(3) + ksign*hr*hr*edges(k)%normal(2)*edges(k)%length
      END DO

      phi(i,2) = phi(i,2) + half*g*s0(2)  ! X-momentum equation.
      phi(i,3) = phi(i,3) + half*g*s0(3)  ! Y-momentum equation.
    END DO

  CASE DEFAULT
    PRINT *,"ERROR: invalid option in bed_slope_DFB."
    CALL byebye('Program STORM stopped.')

  END SELECT

END SUBROUTINE bed_slope_DFB
