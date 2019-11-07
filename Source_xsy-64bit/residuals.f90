SUBROUTINE residuals
  USE parameters
  USE dep_vars
  USE geometry
  USE constants
  USE options
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Computation of the cell residuals.                                         !
!                                                                             !
!  Francisco Simoes, March 2004                                               !
!  Last updated (mm-dd-yyyy): 01-16-2007 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables:
  INTEGER :: i,j,k

  residual = zero  ! Note that these are implicit DO-loops.
  source = zero

! Compute the fluxes.  Convective fluxes are set-up first.
  SELECT CASE (opt_residual)

  CASE (0)  ! Finite element fluxes.
    CALL fe_flux

    ! Add-up the fluxes.
    DO i = 1,n_elems
      ! The convective residuals are subtracted instead of added to
      ! residual(.,.).  That is because this term is moved here to the
      ! righ-hand side of the governing equations.
      DO j = 1,3
        DO k = 1,3
          ! The flux() must be computed using the normal that points outwardly
          ! in the element.  That is accounted for in the SIGN() term.
          ! grid(i)%edge(j) is > 0 if the normal to edge j points out of
          ! element i, and is < 0 otherwise.
          !PRINT *,i,j,k,flux(ABS(grid(i)%edge(j)),k)
          residual(i,k) = residual(i,k) - flux(ABS(grid(i)%edge(j)),k)* &
                          SIGN(one,REAL(grid(i)%edge(j),mp))
        END DO
      END DO
    END DO

  CASE (1)  ! Approximate finite element fluxes.
    CALL afe_flux

    ! Add-up the fluxes, just like for CASE (0).
    DO i = 1,n_elems
      DO j = 1,3
        DO k = 1,3
          residual(i,k) = residual(i,k) - flux(ABS(grid(i)%edge(j)),k)* &
                          SIGN(one,REAL(grid(i)%edge(j),mp))
        END DO
      END DO
    END DO

  CASE (2)
    CALL phisubt(residual)
    ! Move this term to the right-hand side of the equations.
    residual = -residual

  CASE DEFAULT
    PRINT *,''
    PRINT *,'ERROR: invalid value in opt_residual.'
    CALL byebye('Program SToRM stopped.')
  END SELECT


  IF (opt_src == 0) THEN
    ! Now the viscous fluxes.
    CALL visc_terms

    ! Source terms due to bed slope.
    CALL bed_slope

    ! Source terms due to friction slope.
    CALL fric_slope

    ! Other source terms (e.g., nonlinear waves or Boussinesq terms).
    ! ... to do ...

    ! Now add the contributions due to the various source terms, which are
    ! added to residual(.,.) when they are located in the right-hand side of
    ! the governing equations (and subtracted if they appear on the left-hand
    ! side of the equation).
    DO i = 1,n_elems

      DO k = 1,3
        !PRINT *,i,k,source(i,k)
        residual(i,k) = residual(i,k) + source(i,k)
      END DO
      !WRITE (*,'(I9,3(2X,ES12.5))') i,residual(i,1),residual(i,2),residual(i,3)

    END DO
  END IF

END SUBROUTINE
