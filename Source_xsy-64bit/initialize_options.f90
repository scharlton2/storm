SUBROUTINE initialize_options
  USE constants
  USE options
  USE io
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This subroutines sets the default values for all the optional parameters   !
!  needed by SToRM.  Essentially, this means that SToRM will run even if no   !
!  options are specified, because all options have a default value and that   !
!  default value is set in this subroutine.                                   !
!                                                                             !
!  F. Simoes, March 2004                                                      !
!  Last updated (mm-dd-yyyy): 04-25-2013 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables:
  REAL(mp), EXTERNAL :: viscosity

! Printing options.
  output_file = .TRUE.  ! Regurgitate the contents of the run file.
  output_geom = .TRUE.  ! Print computational grid.
  output_head = .TRUE.  ! Print header in output file.
  output_icnd = .TRUE.  ! Print initial conditions.
  output_opts = .TRUE.  ! Print options table.

! Solver type.
  opt_solver = 2  ! FVT scheme.

! Friction coefficient.
  opt_friction = 0  ! Manning's n.

! Water temperature (degrees Centigrade).
  temp = REAL(18,mp)
  visc = viscosity(temp)

! Computation of the residuals.
  opt_residual = 0  ! Finite element interpolation.

! Residual distribution scheme.
  opt_rdscheme = 0  ! The linear N scheme (narrow).

! Dry cell threshold.
  h_dry = 0.001_mp  ! 1 mm

! Threshold for wetting cells.
  h_wet = h_dry*one

! Threshold for shallow flow simplification.
  h_shallow = 0.01_mp

! Stagnation point threshold.
  u_stag = 0.001_mp  ! 1 mm/s

! Cut-of for velocity (magnitude).
  vclip = .FALSE.  ! Option is not active.
  vceiling = 35.0_mp

! CFL time stepping restriction.
  cfl = 1.0_mp

! Stopping criteria for the iterative solution procedure.
  stop_cr = 0
  istop = 0
  rstop = REAL(1.0e-4,mp)

! Solid walls' boundary treatment.
  btype = 1  ! Default is free velocity with zero flux across.
  btype = 0

! Outflow boundary condition.
  hvalve = .TRUE.  ! Activate numerical valve.

! Decide whether to apply the algorithm of Brufau an Garcia-Navarro (2003)
! (zadjust=.TRUE.) or not (zadjust = .FALSE.).
  zadjust = .FALSE.

! Parameters used in the code that determines intermediate solution output.
! These values must not be read from an external file.
  n_iout = 0  ! = 0 before memory allocation.
  p_iout = 1

! Relaxation parameters used in enforcing the stage boundary condition.
  alpha = 0.100_mp
  kappa = 0.0100_mp

! Steady vs. unsteady computations.
  steady = .TRUE.
  dt_op = 0

! Order of the Runge-Kutta method.
  rkorder = 2  ! SSPRK 2nd order.

! Reading and writing mesh connectivity tables.
  read_conn = .FALSE.
  write_conn = .FALSE.

! Define the type of treatment of the source/sink terms.
  opt_src = 0

! Limiter used in the cell gradient reconstruction.
  opt_limiter = 2

! Parabolic eddy viscosity constant.
  parabev = .FALSE.
  omega = 0.068_mp  ! Theoretical value.

! Activation of high-definition least-squares gradient reconstruction.
  activateHDLS = .FALSE.

! Activation of residual averaging.
  activateRSAVG = .FALSE.
  epsRSAVG = half  ! Residual smoothing constant.
  nRSAVG = 2  ! Number of smoothing iterations.

! Variables for wind forcing terms.
  wind_terms = .FALSE.

! This is a relative threshold value below which the equivalent depth of Kuiry
! et al. (2008) is used for the water depth at edges.
  delta_hequiv = 0.020_mp  ! Set at 2% difference between depths on each side.

! This is used as a shock detector to compute the fluxes using a better shock
! capturing method.
  delta_hshock = 0.10_mp  ! Set at 10% difference between depths on each side.

! Special output at "probe" points.
  nprobe = 0
  probeskip = 1  ! Number of iterations skiped between outputs.
  probe_file = ''

! This variable sets the center-to-vertex interpolation technique.
  nctr2vtx = 3

! On start, when using CGNS I/O, printing to Tecplot file is turned off.
  tecplot_print = .FALSE.

! Culvert computations.
  culvert = .FALSE.

END SUBROUTINE initialize_options
