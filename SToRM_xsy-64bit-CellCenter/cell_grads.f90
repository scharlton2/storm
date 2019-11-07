SUBROUTINE cell_grads
  USE parameters
  USE geometry
  USE constants
  USE dep_vars
  USE options
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Computation of the gradients of the dependent variables on each cell.      !
!                                                                             !
!  Francisco Simoes, March 2007                                               !
!  Last updated (mm-dd-yyyy): 10-08-2012 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables.
  INTEGER :: i
  REAL (KIND=mp) :: beta

! Zero gradient, meaning constant element (first-order accurate
! reconstruction):
  IF (opt_limiter == 0) THEN
    DO i = 1,n_elems
      delh(i)%x = zero  ! Flat stage.
      delh(i)%y = zero
      delu(i)%x = zero
      delu(i)%y = zero
      delv(i)%x = zero
      delv(i)%y = zero
    END DO
  END IF

! Second-order accurate reconstruction using least-squares gradients.
  IF (opt_limiter == 1 .OR. opt_limiter == 2) THEN
    IF (activateHDLS) THEN  ! Use higher order least squares.
      !CALL hdlsq_grad_array4h
      !CALL hdlsq_grad_array(h,delh,n_elems)
      CALL hdlsq_grad_array(zeta,delh,n_elems)
      CALL hdlsq_grad_array(u,delu,n_elems)
      CALL hdlsq_grad_array(v,delv,n_elems)
    ELSE  ! Use least squares.
      !CALL lsq_grad_array4h(opt_solver)
      !CALL lsq_grad_array(opt_solver,h,delh,n_elems)
      CALL lsq_grad_array(opt_solver,zeta,delh,n_elems)
      CALL lsq_grad_array(opt_solver,u,delu,n_elems)
      CALL lsq_grad_array(opt_solver,v,delv,n_elems)
      !CALL fe_grad_array(zetavtx,n_pts,delh,n_elems)  ! Use fe approx.
      !CALL fe_grad_array(uvtx,n_pts,delu,n_elems)
      !CALL fe_grad_array(vvtx,n_pts,delv,n_elems)
      !CALL gaussgrad(zeta,delh,n_elems)  ! Gauss-based gradients.
      !CALL gaussgrad(u,delu,n_elems)
      !CALL gaussgrad(v,delv,n_elems)
    END IF
  END IF

  IF (opt_limiter == 2) THEN  ! Apply Barth and Jespersen (1989).
    beta = one  ! Minmod.
    !beta = 2.0_mp  ! Superbee.
    !CALL limiterBJ(beta,h,delh,n_elems)
    ! I think that the limiter should be applied to the stage, not the water
    ! depth.  Hopefully, this solves the problem encountered when the limiter
    ! sets delh to zero: delh should be set to -delz instead.
    CALL limiterBJ(beta,zeta,delh,n_elems)
    CALL limiterBJ(beta,u,delu,n_elems)
    CALL limiterBJ(beta,v,delv,n_elems)
  END IF

  ! In the above computations, delh is actuall the gradient of the free surface
  ! elevetain.  Here, it is reset to the gradient of the water depth.
  DO i = 1,n_elems
    delh(i)%x = delh(i)%x - delz(i)%x
    delh(i)%y = delh(i)%y - delz(i)%y
  END DO

END SUBROUTINE cell_grads
