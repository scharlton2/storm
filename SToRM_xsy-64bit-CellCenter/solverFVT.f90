SUBROUTINE solverFVT(iter,ttime,rmax,fn,title)
  USE constants
  USE geometry
  USE dep_vars
  USE options
  USE io
  USE vbc_arrays
  USE RKparams
  IMPLICIT NONE
  INCLUDE "Header Files\cgnslib_f.h"

!-----------------------------------------------------------------------------!
!                                                                             !
!  Main solution subroutine using the finite volume technique (FVT) loosely   !
!  based on Pan and Cheng (1993) (see also Anastasiou and Chan, 1997).  On    !
!  output it passes the number of iterations performed (iter) and the value   !
!  of the maximum residual (rmax).                                            !
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
!  Francisco Simoes, April 2008                                               !
!  Last updated (mm-dd-yyyy): 06-25-2015 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables.
  INTEGER, INTENT(OUT) :: iter
  REAL(KIND=mp), INTENT(OUT) :: rmax, ttime !rmcd mod 'ttime'
  CHARACTER (LEN=*), INTENT(IN) :: fn,title

! Local variables.
  INTEGER :: i,j,k,m,nimin,nomin,rkstep
  REAL (KIND=mp) :: aeff,delt,denom,dotprod,dt,fac1,fac2,h_delta,heff, &
    ksign,stime,sumbot,sumtop,tmass ,umag
  LOGICAL :: allow_out,fixed_dt,kinwave,probeoutp,zchanged
  REAL (KIND=mp), PARAMETER :: small = 1.0E-9  ! A very small number.
  INTEGER, PARAMETER :: nmax = 8  ! Used in controlling print-out information.

  INTEGER, EXTERNAL :: edge_in_element,local_edge
  INTEGER, EXTERNAL :: distribution
  REAL(KIND=mp), EXTERNAL :: deltatFVT,heffective
  LOGICAL, EXTERNAL :: equals,select_out

! Initialize counters, etc.
  iter = 0
  rmax = zero
  allow_out = n_iout > 0
  stime = zero
  ttime = zero
  fixed_dt = (dt_op == 0)
  fac1 = one/(one + epsRSAVG*3.0_mp)       ! Factors used in the residual
  fac2 = epsRSAVG/(one + epsRSAVG*3.0_mp)  ! smoothing operations.
  ev = zero  ! Eddy viscosity.
  DO i = 1,n_elems
    delh(i)%x = zero; delh(i)%y = zero
    delu(i)%x = zero; delu(i)%y = zero
    delv(i)%x = zero; delv(i)%y = zero
  END DO
  probeoutp = nprobe > 0
  IF (.NOT.probeoutp) probeskip = istop + 1
  zchanged = .FALSE.

! Set-up output of inflow and outflow discharges to screen during run time.
! A maximum of nmax discharges are allowed, counting both inflows and outflows.
! Needs to be consistent with the FORMAT statement where the cqin3() and
! cqout3() information are printed.  nimin is the number of inflows printed on
! screen. Priority is given to inflow information, to a maximum of nmax - 1
! inflow discharges.
  nimin = MIN(n_inflowbdr,nmax-1)
! nomin is the number of outflows printed on screen.  A minimum of 1 outflow
! discharge is printed.
  nomin = MIN(n_outflowbdr,nmax-nimin)

! The gradient of the bed is computed here, but it must be moved into the time
! marching DO-loop when computing bed changes.  Unlimited least-squares
! gradient is used here.
  !IF (activateHDLS) THEN
  !  CALL hdlsq_grad_array(z,delz,n_elems)
  !ELSE
  !  CALL lsq_grad_array(opt_solver,z,delz,n_elems)
  !END IF
! Instead of the code above use subroutine 'bed_slope_vtx', which is simpler
! and possibly more accurate when the bed elevation data is defined at the
! vertices.
  CALL bed_slope_VTX

! The following code is only necessary if subroutine 'bed_slope_ETA' is used.
  !z2 = z*z
  !IF (activateHDLS) THEN
  !  CALL hdlsq_grad_array(z2,delz2,n_elems)
  !ELSE
  !  CALL lsq_grad_array(opt_solver,z2,delz2,n_elems)
  !END IF

  CALL get_time_bcs(stime)

! Enforce boundary conditions here to ensure that the vertex variables have the
! proper values.
  CALL enforce_bcs(hvtx,uvtx,vvtx,zvtx,n_pts)
  CALL zeta_bcs

  IF (probeoutp) CALL probefileo(stime,n_pts,hvtx,uvtx,vvtx)

  rough = cd

  IF (cgns) CALL cg_iric_write_sol_time_f(stime,i)

!-----------------------------------------------------------------------------!
!                                                                             !
!                     M A I N   S O L U T I O N   L O O P                     !
!                                                                             !
!-----------------------------------------------------------------------------!

! Set up initial values for RK stepping.
  DO i = 1,n_elems
    rku(i,rkorder) = h(i)*u(i)
    rkv(i,rkorder) = h(i)*v(i)
    rkh(i,rkorder) = h(i)
  END DO

  DO

    ! Stop criterion: maximum number of iterations is reached.
    iter = iter + 1
    IF (iter > istop) THEN
      iter = iter - 1
      EXIT
    END IF

    DO i = 1,n_elems  ! Set up initial values for RK stepping.
      rku(i,0) = rku(i,rkorder)
      rkv(i,0) = rkv(i,rkorder)
      rkh(i,0) = rkh(i,rkorder)
    END DO

    ! Set-up spill from 1D channels.
    !CALL flood_chan(stime)

    ! Compute time step size.
    IF (.NOT. steady) THEN
      IF (fixed_dt) THEN
        delt = cfl
      ELSE
        delt = deltatFVT(dt_op,cfl)
      END IF
    END IF
    ttime = ttime + delt

    ! SSP Runge-Kutta time stepping.
    DO rkstep = 1,rkorder

      ! Reset some variables.
      cqin = zero   ! Debugging variables containing computed inflow disch.
      cqin2 = zero
      cqin3 = zero
      cqout = zero  ! Debugging variables containing computed outflow disch.
      cqout2 = zero
      cqout3 = zero

      ! Adjustment of the bed elevation of Brufau and Garcia-Navarro (2003).
      ! The following two lines may be commented out to turn the algorithm off.
      ! If this algorithm is used, there is no need for heffective().
      !CALL adjust_bed_elevation(zchanged)
      !IF (zchanged) CALL bed_slope_VTX

      CALL find_dry_cells

      ! The wave celerities are  computed and stored in array 'wavec'.
      wavec = zero
      WHERE (.NOT. drycell) wavec = SQRT(g*h)

      ! Compute the cell gradients of the dependent variables.
      CALL cell_grads
      ! Sets zero grad of stage at wall triangles.
      !CALL hgradfix(rkstep,opt_solver)
      ! Sets zero grad of stage at partially dry triangles.
      CALL hgradfix2(rkstep,opt_solver)

      ! Interpolate solution to the triangle vertices.
      CALL to_vertices

      ! Compute viscous fluxes.
      phi = zero  ! Set cell residuals to zero.
      CALL invisc_fluxes

      ! Compute bed friction explicitely.  Some of the info computed here is
      ! used in the eddy viscosity computations.
      !CALL fric_terms

      ! Parabolic eddy viscosity model based on shear velocity.
      !IF (parabev) THEN
      !  ev = zero
      !  CALL eddy_viscosity
      !END IF

      ! Compute the viscous fluxes.
      !CALL visc_fluxes

      ! Compute the source terms due to bed slope.
      !CALL bed_slope_DFB
      CALL bed_slope_STD
      !CALL bed_slope_ETA(opt_solver)

      ! Add wind forcing terms.  This can be done here, but probably it should
      ! be implemented in the same step as the bed friction, after the Runge-
      ! Kutta cycle, in an implicit manner.
      IF (wind_terms) THEN
        DO i = 1,n_elems
          IF (drycell(i)) CYCLE
          phi(i,2) = phi(i,2) + wcoef1(i)
          phi(i,3) = phi(i,3) + wcoef2(i)
        END DO
      END IF

      ! Add culvert discharges as source/sink terms in the continuity equation.
      IF (culvert) THEN
        CALL cvdisch
        DO i = 1,nculvert
          phi(cvtrigin(i),1) = phi(cvtrigin(i),1) - cvsrc(i)
          phi(cvtrigout(i),1) = phi(cvtrigout(i),1) + cvsrc(i)
        END DO
      END IF

      ! If desired, use the residual smoothing method of Jameson and Mavriplis
      ! (1986).  This is a convergence acceleration technique that should be
      ! used only for steady state runs.
      IF (activateRSAVG) THEN
        phi0 = phi
        ! In this implementation Gauss-Seidel iterations are used instead of
        ! the Jacobi iterations suggested in the original reference.  This
        ! allows to save memory and slightly improve the rate of convergence.
        DO k = 1,nRSAVG
          DO i = 1,n_elems
            phi(i,1) = fac1*phi0(i,1)
            phi(i,2) = fac1*phi0(i,2)
            phi(i,3) = fac1*phi0(i,3)
            DO m = 2,t2t(i,1) + 1
              j = t2t(i,m)  ! Neighbor element used in the Laplacian operator.
              phi(i,1) = phi(i,1) + fac2*phi(j,1)
              phi(i,2) = phi(i,2) + fac2*phi(j,2)
              phi(i,3) = phi(i,3) + fac2*phi(j,3)
            END DO
          END DO
        END DO
      END IF

      ! Update dependent variables.
      DO i = 1,n_elems
        IF (steady) THEN  ! Local time step for each element.
          delt = u_mag(i) + wavec(i)
          IF (delt < small) delt = 1.044_mp
          ! Use CFL criterion from Brufau et al. (2004).
          delt = cfl*grid(i)%area/(delt*grid(i)%perim)
        END IF

        ! Set-up the pseudo-time marching scheme.
        ! aeff = aeffective(i)
        aeff = grid(i)%area
        dt = delt/aeff
        rkh(i,rkstep) = zero
        rku(i,rkstep) = zero
        rkv(i,rkstep) = zero
        DO j = 0,rkstep - 1
          rkh(i,rkstep) = rkh(i,rkstep) + rkalpha(rkstep,j)*rkh(i,j)
          rku(i,rkstep) = rku(i,rkstep) + rkalpha(rkstep,j)*rku(i,j)
          rkv(i,rkstep) = rkv(i,rkstep) + rkalpha(rkstep,j)*rkv(i,j)
        END DO
        rkh(i,rkstep) = rkh(i,rkstep) + rkbeta(rkstep)*dt*phi(i,1)
        heff = MAX(h_dry,rkh(i,rkstep))  ! Prevents division by zero.
        rku(i,rkstep) = (rku(i,rkstep) + rkbeta(rkstep)*dt*phi(i,2))
        rkv(i,rkstep) = (rkv(i,rkstep) + rkbeta(rkstep)*dt*phi(i,3))
      END DO

      ! Perform velocity clipping within the Runge-Kutta cycle.
      IF (vclip) THEN
        DO i = 1,n_elems
          heff = MAX(h_dry,rkh(i,rkstep))
          umag = SQRT(rku(i,rkstep)*rku(i,rkstep) + &
                 rkv(i,rkstep)*rkv(i,rkstep))/heff
          IF (umag > vceiling) THEN
            rku(i,rkstep) = vceiling*rku(i,rkstep)/umag
            rkv(i,rkstep) = vceiling*rkv(i,rkstep)/umag
          END IF
        END DO
      END IF

      ! Set-up the variables for the next Runge-Kutta step and process some
      ! wetting and drying.
      DO i = 1,n_elems
        rkh(i,rkstep) = MAX(rkh(i,rkstep),zero)
        h(i) = rkh(i,rkstep)

        ! 'h_delta' is the threshold for wetting and drying.  Hysteresis
        ! determines that if a cell is dry, it will only get wet if the depth
        ! becomes larger than 'h_wet'; if the cell is wet, it will only become
        ! dry if the depth becomes smaller that 'h_dry'.
        IF (drycell(i)) THEN
          h_delta = h_wet
        ELSE
          h_delta = h_dry
        END IF

        IF (h(i) > h_delta) THEN
          !heff = heffective(i)  ! Effective depth for triangle i.
          heff = h(i)
          u(i) = rku(i,rkstep)/heff
          v(i) = rkv(i,rkstep)/heff
        ELSE
          rku(i,rkstep) = zero
          rkv(i,rkstep) = zero
          u(i) = zero
          v(i) = zero
        END IF
      END DO

      ! Correct the water depth if the bed adjustment algorithm was used.  No
      ! need to comment this line out if the algorithm is not used...
      !IF (zchanged) CALL adjust_water_depth

    END DO  ! End Runge-Kutta step.


    ! Implicit implementation of the friction terms.
    CALL find_dry_cells
    CALL fric_terms_impl(delt)

    ! These variables are used to compute the RMS change in water depth.  This
    ! is one of the criteria used to determine convergence to a steady state
    ! solution (steady state convergence criterion #1).
    sumtop = zero
    sumbot = zero
    DO i = 1,n_elems
      sumtop = sumtop + (h(i) - rkh(i,0))*(h(i) - rkh(i,0))
      sumbot = sumbot + rkh(i,0)*rkh(i,0)
    END DO

    !CALL overdraft
    !DO i = 1,n_elems  ! Make sure that h > 0 everywhere to avoid math problems.
    !   h(i) = MAX(h(i),zero)  ! This is faster...
    !  !IF (h(i) < zero) THEN  !.. but this provides a diagnose.
    !  !  h(n_elems+1) = h(n_elems+1) + h(i)*grid(i)%area/grid(n_elems+1)%area
    !  !  h(i) = zero
    !  !END IF
    !END DO

    IF (hvalve) THEN
      DO j = 1,n_bcedges
        IF (.NOT. hbar(j)) CYCLE  ! Stage boundary edges.
        k = bcedges(j)  ! Edge.
        i = edge_in_element(edges(k))  ! Element.
        m = local_edge(grid(i),k)  ! Element i local index to edge k.
        ksign = SIGN(one,REAL(grid(i)%edge(m),mp))
        dotprod = ksign*(u(i)*edges(k)%normal(1) + v(i)*edges(k)%normal(2))
        IF (dotprod < 0) THEN
          u(i) = zero
          v(i) = zero
        END IF
      END DO
    END IF
    ! Note about the code above: the 'no inflow' condition is better if
    ! implemented directly at the edges, in subroutine 'flux_ac'.

    ! The values of the residual must be reset to zero for dry nodes and in
    ! elements where the kinematic overland flow model is used.
    DO i = 1,n_elems
      IF (drycell(i)) THEN
        phi(i,1) = zero
        phi(i,2) = zero
        phi(i,3) = zero
      END IF
      !IF (h(i) < h_shallow) THEN
      !  phi(i,2) = zero
      !  phi(i,3) = zero
      !END IF
    END DO

    ! Check if intermediary output is requested and accomplish it.
    IF (allow_out) THEN
      IF (select_out(iter)) THEN
        ! Make solution file pretty: interpolate solution to triangle vertices
        ! and enforce boundary conditions at inflow/outflow strings and channel
        ! margins.
        CALL to_vertices

        CALL error_est(2)
        IF (cgns) THEN
          CALL ctr2vtx(nctr2vtx,cd,cdvtx)
          CALL output_data(iter,ttime,rmax,fn,title)
        ELSE
          CALL spit_out(iter,fn,title)
          !CALL spit_outv(iter,fn,title)
        ENDIF
      END IF
    END IF

    IF (steady) THEN  ! Steady state stopping criteria.
      SELECT CASE (stop_cr)

      CASE (0)  ! Use the cell residuals and stop if the maximum cell residual
                ! in the is smaller than a defined threshold.
        rmax = zero
        DO i = 1,n_elems
          rmax = MAX(rmax,ABS(phi(i,1)),ABS(phi(i,2)),ABS(phi(i,3)))
        END DO

      CASE (1)  ! Use the relative error method, eq. (36) of Lai et al. (2005).
        IF (equals(sumbot,zero)) sumbot = one  ! To avoid division by zero.
        rmax = SQRT(sumtop/sumbot)

      END SELECT
      tmass = zero  ! Total mass.
      DO i = 1,n_elems
        tmass = tmass + h(i)*grid(i)%area
      END DO
      !WRITE (*,'(I8,4ES14.5)') iter,rmax,tmass
      WRITE (*,'(I8,4ES14.5)') iter,rmax,tmass,ABS(cqin3(1)),cqout3(1)
      WRITE (*,'(8X,ES14.5)') (ABS(cqin3(i)),i=2,n_inflowbdr)
      WRITE (*,'(8X,ES14.5)') (cqout3(i),i=2,n_outflowbdr)

      IF (rmax < rstop) EXIT

    ELSE  ! March in time.
      stime = stime + delt
      ! Compute and print the maximum cell residual for informative purposes.
      rmax = zero
      DO i = 1,n_elems
        rmax = MAX(rmax,ABS(phi(i,1)),ABS(phi(i,2)),ABS(phi(i,3)))
      END DO
      tmass = zero  ! Total mass.
      DO i = 1,n_elems
        tmass = tmass + h(i)*grid(i)%area
      END DO
      IF (equals(sumbot,zero)) sumbot = one  ! To avoid division by zero.
      ! The number of ES13.4 values in the FORMAT statement in the next line
      ! must be equal to 3 + nmax.
      sumbot = SQRT(sumtop/sumbot)
      WRITE (*,'(I8,11ES13.4)') iter,stime,sumbot,tmass, &  ! Statement for degugging.
                                ABS(qin(1)),cqout3(1)
      !WRITE (*,'(I8,11ES13.4)') iter,stime,sumbot,tmass, &  !!!*THIS IS THE RIGHT STATEMENT*!!!
      !                          ABS(cqin3(1)),cqout3(1)
      !WRITE (*,'(I8,11ES13.4)') iter,stime,SQRT(sumtop/sumbot),tmass, &
      !                          (ABS(cqin3(i)),i=1,nimin), &
      !                          (cqout3(i),i=1,nomin)
      !WRITE (*,'(I8,4ES14.5)') iter,stime,rmax,SQRT(sumtop/sumbot),tmass
      !WRITE (*,'(I8,6ES14.5)') iter,cqin(1),cqin2(1),cqin3(1),cqout(1), &
      !  cqout2(1),cqout3(1)
      !WRITE (*,'(8X,3ES14.5)') (cqin(i),cqin2(i),cqin3(i),i=2,n_inflowbdr)
      !WRITE (*,'(8X,3ES14.5)') (cqout(i),cqout2(i),cqout3(i),i=2,n_outflowbdr)

      ! Get new boundary condition values.
      CALL get_time_bcs(stime)

      ! Write iterative solution parameters to cgns file.
      !IF (cgns) THEN
        !CALL cg_iric_write_sol_time_f(ttime,i)
        !CALL gc_iric_write_sol_baseiterative_real_f('Inflow1',cqin3(1),i)
        !CALL gc_iric_write_sol_baseiterative_real_f('Outflow1',cqout3(1),i)
        !CALL gc_iric_write_sol_baseiterative_real_f('TotalVolume',tmass,i)
        !CALL gc_iric_write_sol_baseiterative_real_f('Residual',sumbot,i)
      !END IF
    END IF

    ! Write probe output.
    IF (steady) stime = REAL(iter,mp)
    IF (MOD(iter,probeskip) == 0) CALL probefileo(stime,n_pts,hvtx,uvtx,vvtx)

  END DO

! Make solution file pretty.  The cell gradients of the dependent variables are
! needed for subroutine 'wetdryfix'.  This needs to be done every time that is
! necessary to output the vertex values (for iRIC, for example).
  CALL to_vertices
  CALL ctr2vtx(nctr2vtx,cd,cdvtx)

  ! Compute error estimates.
  CALL error_est(2)

  tecplot_print = .TRUE.  ! Debug: print to Tecplot file on exit.

END SUBROUTINE solverFVT
