SUBROUTINE error_est(etype)
  USE parameters
  USE constants
  USE geometry
  USE dep_vars
  USE options
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Computes `a posteriori error estimates for the FVT solver.                 !
!                                                                             !
!  INPUT:                                                                     !
!    etype  error estimate chosen.  At the moment there is only one...        !
!                                                                             !
!  Francisco Simoes, July 2007                                                !
!  Last updated (mm-dd-yyyy): 08-01-2007 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  INTEGER :: etype

! Local variables.
  INTEGER :: i,n1,n2,n3
  REAL (KIND=mp) :: gradetax,gradetay,gradzx,gradzy,l

  e_est = zero
  SELECT CASE (etype)

  CASE (1)

    ! This technique uses the root mean square values of the free surface
    ! gradient, as in eq. (30) of Liang et al. (2007).  This estimator is
    ! only good for smooth flows without hydraulic jumps.  I doubt it is good
    ! for anything...
    DO i = 1,n_elems
      IF (h(i) < h_dry) CYCLE

      ! Compute the finite element gradient of the bed slope.
      n1 = grid(i)%vertex(1)
      n2 = grid(i)%vertex(2)
      n3 = grid(i)%vertex(3)
      gradzx = (zvtx(n1)*(nodes(n2)%y - nodes(n3)%y) + &  ! X-component.
               zvtx(n2)*(nodes(n3)%y - nodes(n1)%y) + &
               zvtx(n3)*(nodes(n1)%y - nodes(n2)%y))*half/grid(i)%area
      gradzy = (zvtx(n1)*(nodes(n3)%x - nodes(n2)%x) + &  ! Y-component.
               zvtx(n2)*(nodes(n1)%x - nodes(n3)%x) + &
               zvtx(n3)*(nodes(n2)%x - nodes(n1)%x))*half/grid(i)%area

      ! Find the gradient of the free surface elevation.
      gradetax = gradzx + delh(i)%x
      gradetay = gradzy + delh(i)%y
      e_est(i) = SQRT(gradetax*gradetax + gradetay*gradetay)
    END DO

  CASE (2)

    ! This technique is very similar to the technique above, but uses a form of
    ! scaling inspired by the thoughts in page 32 of Holmes and Connell (1989).
    ! Just like for etype = 1, this technique is good only for smooth flows
    ! without hydraulic jumps.
    DO i = 1,n_elems
      IF (h(i) < h_dry) CYCLE

      ! Compute the finite element gradient of the bed slope.
      n1 = grid(i)%vertex(1)
      n2 = grid(i)%vertex(2)
      n3 = grid(i)%vertex(3)
      gradzx = (zvtx(n1)*(nodes(n2)%y - nodes(n3)%y) + &  ! X-component.
               zvtx(n2)*(nodes(n3)%y - nodes(n1)%y) + &
               zvtx(n3)*(nodes(n1)%y - nodes(n2)%y))*half/grid(i)%area
      gradzy = (zvtx(n1)*(nodes(n3)%x - nodes(n2)%x) + &  ! Y-component.
               zvtx(n2)*(nodes(n1)%x - nodes(n3)%x) + &
               zvtx(n3)*(nodes(n2)%x - nodes(n1)%x))*half/grid(i)%area

      ! Compute element characteristic size.
      l = grid(i)%perim*one_third

      ! Find the gradient of the free surface elevation.
      gradetax = gradzx + delh(i)%x
      gradetay = gradzy + delh(i)%y
      e_est(i) = SQRT(gradetax*gradetax + gradetay*gradetay)*l/(z(i) + h(i))
    END DO

  CASE DEFAULT
    PRINT *,''
    PRINT *,'ERROR: invalid argument in subroutine error_est.'
    CALL byebye('Program SToRM stopped.')

  END SELECT

END SUBROUTINE error_est
