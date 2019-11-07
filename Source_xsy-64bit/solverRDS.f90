SUBROUTINE solverRDS(iter,rmax,fn,title)
  USE constants
  USE geometry
  USE dep_vars
  USE options
  USE io
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Main solution subroutine using the Residual Distrubution Scheme (RDS).     !
!  On output it passes the number of iterations performed (iter) and the      !
!  value of the maximum residual (rmax).                                      !
!                                                                             !
!  INPUT:                                                                     !
!    fn     filename for the output files.  The names of the files that       !
!           contain the intermediate output are built from this name;         !
!    title  title for Tecplot file.                                           !
!                                                                             !
!  OUTPUT:                                                                    !
!    iter   maximum iteration before stopping criteria was reached;           !
!    rmax   maximum residual or other value used in the stopping criteria     !
!           (i.e., quantity that is checked for to determine convergence).    !
!                                                                             !
!  Francisco Simoes, March 2004                                               !
!  Last updated (mm-dd-yyyy): 08-14-2007 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables.
  INTEGER, INTENT(OUT) :: iter
  REAL(KIND=mp), INTENT(OUT) :: rmax
  CHARACTER (LEN=*), INTENT(IN) :: fn,title

! Local variables.
  INTEGER :: counter,i
  REAL(KIND=mp) :: delt,dt,oldh,sumbot,sumtop,time
  LOGICAL :: allow_out,fixed_dt,stop_cr1
  REAL (KIND=mp), PARAMETER :: small = 1.0E-9  ! A very small number.

  INTEGER, EXTERNAL :: distribution
  REAL(KIND=mp), EXTERNAL :: deltat
  LOGICAL, EXTERNAL :: select_out

! Initialize counters, etc.
  iter = 0
  rmax = zero
  allow_out = n_iout > 0
  time = zero
  stop_cr1 = (stop_cr == 1)
  fixed_dt = (dt_op == 0)

  DO

    ! Stop criterion: maximum number of iterations is reached.
    iter = iter + 1
    IF (iter > istop) THEN
      iter = iter - 1
      EXIT
    END IF

    CALL enforce_bcs(h,u,v,z,n_pts)

    ! Element linearization.
    CALL cell_velocity

    CALL residuals

    IF (opt_src == 1) CALL fvsources

    counter = distribution()

    ! Compute time step size.
    IF (.NOT. steady) THEN
      IF (fixed_dt) THEN
        delt = cfl
      ELSE
        delt = deltat(dt_op,cfl)
      END IF
    END IF

    IF (stop_cr1) THEN  ! Used for steady state convergence criterion #1.
      sumtop = zero
      sumbot = zero
    END IF

    ! Update dependent variables.
    DO i = 1,n_pts
      IF (steady) THEN  ! Local time step for each element.
        delt = SQRT(u(i)*u(i) + v(i)*v(i)) + SQRT(g*h(i)) + small
        ! Use CFL criterion from Brufau et al. (2004).
        delt = cfl*cv_area(i)/(delt*cv_perim(i))
      END IF

      ! Set-up the pseudo-time marching scheme.
      dt = delt/cv_area(i)
      oldh = h(i)
      h(i) = oldh + dt*phi(i,1)
      ! Make sure that h > 0 everywhere to avoid math problems.
      h(i) = MAX(h(i),zero)
      IF (h(i) > h_wet) THEN
        IF (opt_src == 0) THEN
          u(i) = (oldh*u(i) + dt*phi(i,2))/h(i)
          v(i) = (oldh*v(i) + dt*phi(i,3))/h(i)
        ELSE
          u(i) = (oldh*u(i) + dt*(phi(i,2) + source(i,2)))/h(i)
          v(i) = (oldh*v(i) + dt*(phi(i,3) + source(i,3)))/h(i)
        END IF
      END IF

      IF (stop_cr1) THEN
        sumtop = sumtop + (h(i) - oldh)*(h(i) - oldh)
        sumbot = sumbot + oldh*oldh
      END IF
    END DO

    ! Check if intermediary output is requested and accomplish it.
    IF (allow_out) THEN
      IF (select_out(iter)) THEN
        CALL errore
        CALL spit_out(iter,fn,title)
      END IF
    END IF

    IF (steady) THEN  ! Steady state stopping criteria.
      SELECT CASE (stop_cr)

      CASE (0)  ! Use the cell residuals and stop if the maximum cell residual
                ! in the is smaller than a defined threshold.
        rmax = zero
        DO i = 1,n_pts
          rmax = MAX(rmax,ABS(phi(i,1)),ABS(phi(i,2)),ABS(phi(i,3)))
        END DO

      CASE (1)  ! Use the relative error method of Lai et al. (2006).
        rmax = SQRT(sumtop/sumbot)

      END SELECT
      PRINT *,iter,counter,rmax
      IF (rmax < rstop) EXIT

    ELSE  ! March in time.
      time = time + delt
      ! Compute and print the maximum cell residual for informative purposes.
      rmax = zero
      DO i = 1,n_pts
        rmax = MAX(rmax,ABS(phi(i,1)),ABS(phi(i,2)),ABS(phi(i,3)))
      END DO
      PRINT *,iter,time,rmax
    END IF

  END DO

! Make solution file pretty.
  CALL enforce_bcs(h,u,v,z,n_pts)
  DO i = 1,n_pts
    IF (h(i) < h_dry) h(i) = zero
  END DO

  ! Compute error estimates.
  CALL errore

END SUBROUTINE solverRDS
