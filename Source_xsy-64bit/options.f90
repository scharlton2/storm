!-----------------------------------------------------------------------------!
!                                                                             !
!  This file contains variables used in setting-up the options necessary for  !
!  a SToRM run.                                                               !
!                                                                             !
!  F. Simoes, March 2004                                                      !
!  Last updated (mm-dd-yyyy): 04-25-2013 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

MODULE options
  USE parameters
  IMPLICIT NONE
  SAVE

! Printing options.
  LOGICAL :: output_file,output_geom,output_head,output_icnd,output_opts

! Solver type: = 1 for RDS (default), = 2 for finite volumes based (triangles),
! = 3 for another finite volume solver based on triangles (FVZ).
  INTEGER :: opt_solver

! Friction coefficient.
  INTEGER :: opt_friction
  CHARACTER (LEN=40) :: f_fric_facts

! Water temperature (degrees Centigrade) and viscosity.
  REAL (KIND=mp) :: temp,visc

! Technique to compute the residual (= 0 for finite element interpolation
! (default), = 1 for approximate finite element interpolation, or = 2 using the
! conservative residual distribution approach, i.e., the Jacobian.
  INTEGER :: opt_residual

! Residual distribution scheme.
  INTEGER :: opt_rdscheme

! Threshold to determine if a node or cell is dry (or wet).
  REAL (KIND=mp) :: h_dry,h_wet

! Defines the depth to discard the convective terms in the momentum equations.
  REAL (KIND=mp) :: h_shallow

! Threshold for a stagnation point.
  REAL (KIND=mp) :: u_stag

! Cut-of for velocity (magnitude).
  LOGICAL :: vclip  ! = .TRUE. if option is activated.
  REAL (KIND=mp) ::vceiling

! Courrant-Friedrichs-Levy number.
  REAL (KIND=mp) :: cfl

! Stopping criteria for the iterative solution procedure.  If stop_cr = 0, the
! stopping criterion is the maximum value of the residuals; if stop_cr = 1, the
! stoping criterion is the relative error of the water depth, as defined in eq.
! (36) of Lai et al. (2005).
  INTEGER :: stop_cr
  INTEGER :: istop  ! Maximum number of iterations.
  REAL (KIND=mp) :: rstop  ! Maximum allowed value of the selected quantity.

! Type of boundary condition applied to solid walls: btype = 0 for zero
! velocity, = 1 for zero flux but free velocity.
  INTEGER :: btype

! Decide whether to apply a numerical valve at the boundary where the stage
! boundary condition is enforced.  If this boundary coincides with the
! computational domain's outlet, it should be set to .TRUE., otherwise it may
! have to be set to .FALSE..  A typical exemple is subcritical flow (.TRUE.)
! and supercritical flow (.FALSE.). Default is .TRUE..
  LOGICAL :: hvalve

! Decide whether to apply the algorithm of Brufau an Garcia-Navarro (2003)
! (zadjust=.TRUE.) or not (zadjust = .FALSE.).
  LOGICAL :: zadjust

! Relaxation parameters used in enforcing the stage boundary condition.
  REAL (KIND=mp) :: alpha,kappa

! Steay vs. unsteady computations.  Variable dt_op is used to select the method
! to compute the time step: dt_op = 0 for fixed time step size (set in variable
! cfl); = 1 for automatic time step size.  For now it is only used in unsteady
! computations, but its use can be expanded to steady state computations
! later.
  LOGICAL :: steady  ! = .TRUE. for steady computations.
  INTEGER :: dt_op

! Time stepping method chosen.  Currently this is done with a general
! implementation of the Runge-Kutta method and this parameter is the order of
! the method.
  INTEGER :: rkorder

! Options to read and write the mesh connectivity to external files.  read_conn
! =.TRUE. to read the mesh connectivity from file f_read_conn; otherwise, set
! read_conn=.FALSE., which is the default.  To write the connectivity tables to
! file f_write_conn set write_conn=.TRUE.
  LOGICAL :: read_conn,write_conn
  CHARACTER (LEN=40) ::f_read_conn,f_write_conn

! Type of treatment of the source/sink terms of the governing equations.
! opt_src = 0 for the residual-based approach, = 1 for the finite volume
! formulation.
  INTEGER :: opt_src

! Limiter used in the cell gradient reconstruction process:  = 0 for zero
! gradient; = 1 for no limiter (default); = 2 for Barth and Jespersen (1989).
  INTEGER :: opt_limiter

! Parabolic eddy viscosity model.
  REAL (KIND=mp) :: omega  ! Eddy viscosity constant.
  LOGICAL :: parabev  ! = .TRUE. if model is activated, = .FALSE. otherwise.

! Activation of high-definition least-squares gradient reconstruction.
  LOGICAL :: activateHDLS

! Activation of residual averaging.
  LOGICAL :: activateRSAVG
  INTEGER :: nRSAVG
  REAL (KIND=mp) :: epsRSAVG

! Culvert computations.
  LOGICAL :: culvert  ! .TRUE. if culvert computations are active.
  CHARACTER (LEN=256) :: cvfile  ! Path and name of file with culvert data.

! Variables for wind forcing terms.
  LOGICAL :: wind_terms
  CHARACTER (LEN=40) :: wind_file

! Path for Triangle edge file output
   CHARACTER(LEN=250) :: edgeFilePath, vedgeFilepath !rmcd mod

! Threshold value to apply the equivalent depth at cell edges.  It is a
! relative difference between the depths at each side of the edge and must be
! in the range [0,1], but should typically be < 0.05 (i.e., 5%).  Set to zero
! to deactivate the usage of the equivalent depth.
  REAL (KIND=mp) :: delta_hequiv

! This is used as a shock detector to trigger the use of a shock capturing
! method at cell edges that have a relative depth difference larger than the
! threshold specified by this parameter.  Must be in the range [0,1], but
! should typically > 0.05 (i.e., 5%).  Set to zero to always use the shock
! capturing method, set to 1 to deactivate it.
  REAL (KIND=mp) :: delta_hshock

! Variables for probe output.  nprobe is the number of probe points; probeskip
! is the number of iterations to skip between output; probepts(nprobe) is an
! array containing the node numbers (vertices) where output is collected;
! probe_format is the FORMAT used in the probed data WRITE statements; and
! probe_file is the output file name.
  INTEGER :: nprobe,probeskip
  INTEGER, ALLOCATABLE, DIMENSION(:) :: probepts
  CHARACTER (LEN=40) :: probe_format,probe_file

! This variable sets the center-to-vertex interpolation, which is argumet
! 'type' in subroutines 'c2v_setup' and 'ctr2vtx': = 1 area-weighted
! interpolation,  = 2 for the least-squares, and = 3 for the inverse distance
! weighting.
  INTEGER :: nctr2vtx

! A variable the allows printing the solution in a Tecplot-compatible ASCII
! file even when the CGNS I/O is active.
  LOGICAL :: tecplot_print

END MODULE options
