SUBROUTINE read_options(funit)
  USE parameters
  USE options
  USE constants
  USE io
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This routine implements a parser to read through the options file and      !
!  initialize the appropriate variables.                                      !
!                                                                             !
!  Francisco Simoes, October 2004                                             !
!  Last updated (mm-dd-yyyy): 04-22-2013 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables:
  INTEGER, INTENT(IN) :: funit

! Local variables:
  INTEGER :: c,counter,i,ierror,lineno,n,p
  REAL(KIND=mp) :: coef
  CHARACTER (LEN=40) :: line,opt
  LOGICAL :: flag,fnexist

  REAL(KIND=mp), EXTERNAL :: viscosity
  INTEGER, EXTERNAL :: readp
  CHARACTER, EXTERNAL :: to_upper
  LOGICAL, EXTERNAL :: equals,readline

  lineno = 0
  coef = zero

  DO

    flag = readline(funit,line,lineno)
    IF (.NOT. flag) EXIT  ! Error or end-of-file encountered.

    ! * * * * * * * * * * Water temperature * * * * * * * * * * !
    IF (line(1:12) == "TEMPERATURE") THEN
      n = LEN_TRIM(line)
      c = IACHAR(line(n:n))  ! ASCII decimal value of line(n:n).
      IF ((c >= 48 .AND. c <=57) .OR. c == 46) THEN  ! Default: Centigrade.
        READ(line(13:40),*) temp
      ELSE IF (c == 67) THEN  ! Also Centigrade, but with explicit unit symbol.
        READ(line(13:n-1),*) temp
      ELSE IF (c == 70) THEN  ! Degrees Fahrenheit.
        READ(line(13:n-1),*) temp
        temp = REAL(5,mp)*(temp - REAL(32,mp))/REAL(9,mp)
      ELSE
        PRINT *,''
        PRINT *,'ERROR: illegal value in option temperature (units)'
        PRINT *,'in line',lineno,' of the options file.'
        CALL byebye('Program SToRM stopped.')
      END IF
      visc = viscosity(temp)
      IF (equals(visc,zero)) CALL err_opts('temperature',lineno)

    ! * * * * * * * * * * * * * * * Viscosity * * * * * * * * * * * * * * * !
    ELSE IF (line(1:9) == "VISCOSITY") THEN
      READ(line(10:40),*) visc
      IF (visc < 0) CALL err_opts('viscosity',lineno)

    ! * * * * * * * * * * * * * * * * Solver * * * * * * * * * * * * * * * * !
    ELSE IF (line(1:6) == "SOLVER") THEN
      opt = ADJUSTL(line(7:))
      IF (opt(1:3) == 'RDS') THEN
        opt_solver = 1
      ELSE IF (opt(1:3) == 'FVT') THEN
        opt_solver = 2
      ELSE IF (opt(1:3) == 'FVZ') THEN
        opt_solver = 3
      ELSE
        CALL err_opts('solver',lineno)
      END IF

    ! * * * * * * * * * * Computation of the residual * * * * * * * * * * !
    ELSE IF (line(1:8) == "RESIDUAL") THEN
      opt = ADJUSTL(line(9:))
      IF (opt(1:2) == 'FE') THEN
        opt_residual = 0
      ELSE IF (opt(1:3) == 'AFE') THEN
        opt_residual = 1
      ELSE IF (opt(1:3) == 'RDS') THEN
        opt_residual = 2
      ELSE
        CALL err_opts('residual',lineno)
      END IF

    ! * * * * * * * * * * Residual distribution scheme * * * * * * * * * * !
    ELSE IF (line(1:9) == "RD SCHEME") THEN
      opt = ADJUSTL(line(10:))
      IF (opt(1:1) == 'N') THEN
        opt_rdscheme = 0
      ELSE IF (opt(1:3) == 'LDA') THEN
        opt_rdscheme = 1
      ELSE IF (opt(1:1) == 'B') THEN
        opt_rdscheme = 2
      ELSE IF (opt(1:3) == 'PSI') THEN
        opt_rdscheme = 3
      ELSE
        CALL err_opts('rd scheme',lineno)
      END IF

    ! * * * * * * * * * * Type of friction coefficient * * * * * * * * * * !
    ELSE IF (line(1:9) == "ROUGHNESS") THEN
      READ(line(10:40),*) opt_friction
      IF (opt_friction < 0 .OR. opt_friction > 5) &
        CALL err_opts('roughness',lineno)
      ! Read file name where friction coefficients are stored.  Use line below
      ! for the file name.
      IF (opt_friction == 5) THEN
        flag = readline(funit,line,lineno)
        READ(line(1:40),'(A)') f_fric_facts
      END IF

    ! * * * * * * * Option to regurgitate the main input file * * * * * * * !
    ELSE IF (line(1:10) == "MASTERFILE") THEN
      READ(line(11:40),*) i
      IF (i <= 0) THEN
        output_file = .FALSE.
      ELSE
        output_file = .TRUE.
      END IF

    ! Option to regurgitate the geometry data (nodal coordinates, control
    ! volume connectivity, etc.).
    ELSE IF (line(1:8) == "GEOMETRY") THEN
      READ(line(9:40),*) i
      IF (i <= 0) THEN
        output_geom = .FALSE.
      ELSE
        output_geom = .TRUE.
      END IF

    ! Option to create a header in the output file with the run's time stamp.
    ELSE IF (line(1:6) == "HEADER") THEN
      READ(line(7:40),*) i
      IF (i <= 0) THEN
        output_head = .FALSE.
      ELSE
        output_head = .TRUE.
      END IF

    ! * * * * * * * * Option to echo the initial conditions * * * * * * * * !
    ELSE IF (line(1:13) == "INITIAL CONDS") THEN
      READ(line(14:40),*) i
      IF (i <= 0) THEN
        output_icnd = .FALSE.
      ELSE
        output_icnd = .TRUE.
      END IF

    ! * * * * * * * * * * Option to echo all the options * * * * * * * * * * !
    ELSE IF (line(1:7) == "OPTIONS") THEN
      READ(line(8:40),*) i
      IF (i <= 0) THEN
        output_opts = .FALSE.
      ELSE
        output_opts = .TRUE.
      END IF

    ! * * * * * * * * * * * Threshold for shallow cells * * * * * * * * * * * !
    ELSE IF (line(1:9) == "H SHALLOW") THEN
      !READ(line(10:40),*) h_shallow
      !h_shallow = MIN(one,ABS(h_shallow))
      WRITE (*,'(/"WARNING: option H SHALLOW in options file is an obsolete", &
        " feature.")')
      WRITE (*,'("H SHALLOW value ignored in line",I4," of options file."/)') &
        lineno

    ! * * * * * * * * * * * * Threshold for dry cells * * * * * * * * * * * * !
    ELSE IF (line(1:11) == "H THRESHOLD") THEN
      READ(line(12:40),*) h_dry
      h_dry = MIN(one,ABS(h_dry))

    ! * * * * * * * * * * * * Threshold for wet cells * * * * * * * * * * * * !
    ELSE IF (line(1:11) == "C THRESHOLD") THEN
      READ(line(12:40),*) coef
      coef = ABS(coef)
      coef = MAX(coef,one)

    ! * * * * * * * * * * Threshold for stagnant cells * * * * * * * * * * !
    ELSE IF (line(1:11) == "U THRESHOLD") THEN
      READ(line(12:40),*) u_stag
      u_stag = MIN(one,ABS(u_stag))

    ! * * * * * * * * * * * Velocity magnitude ceiling * * * * * * * * * * * !
    ELSE IF (line(1:12) == "MAX VELOCITY") THEN
      READ(line(13:40),*) vceiling
      vceiling = ABS(vceiling)
      vclip = .TRUE.

    ! * * * * * * * CFL number for time step size calculation * * * * * * * !
    ELSE IF (line(1:3) == "CFL") THEN
      READ(line(4:40),*) cfl
      cfl = ABS(cfl)

    ! * * * * * * * * * * Steady state stopping criterion * * * * * * * * * * !
    ELSE IF (line(1:13) == "STOPPING CRIT") THEN
      READ(line(14:40),*) stop_cr
      IF (stop_cr < 0 .OR. stop_cr > 1) CALL err_opts('stopping crit',lineno)

    ! * * * * * * * * * * Maximum number of iterations * * * * * * * * * * !
    ELSE IF (line(1:13) == "MAX ITERATION") THEN
      READ(line(14:40),*) istop

    ! * * * * * * * * * * Convergence criterion threshold * * * * * * * * * * !
    ELSE IF (line(1:12) == "MAX RESIDUAL") THEN
      READ(line(13:40),*) rstop
      rstop = ABS(rstop)

    ! * * * * * * * * * Treatment of the solid boundaries * * * * * * * * * !
    ELSE IF (line(1:5) == "WALLS") THEN
      opt = ADJUSTL(line(6:))
      IF (opt(1:3) == 'NOS') THEN
        btype = 0
      ELSE IF (opt(1:4) == 'FREE') THEN
        btype = 1
      ELSE
        CALL err_opts('walls',lineno)
      END IF

    ! * * * * * * * * * * * * Reflux valve at outlet * * * * * * * * * * * * !
    ELSE IF (line(1:6) == "HVALVE") THEN
      READ(line(7:40),*) i
      IF (i <= 0) THEN
        hvalve = .FALSE.
      ELSE
        hvalve = .TRUE.
      END IF

    ! * * * * * * Algorithm of Brufau an Garcia-Navarro (2003) * * * * * * !
    ELSE IF (line(1:13) == "BGN ALGORITHM") THEN
      READ(line(14:40),*) i
      IF (i <= 0) THEN
        zadjust = .FALSE.
      ELSE
        zadjust = .TRUE.
      END IF

    ! * * * * * Select where to print-out intermediate solutions * * * * * !
    ELSE IF (line(1:12) == "INT SOLUTION") THEN
      READ(line(13:40),*) n_iout
      IF (n_iout <= 0) CALL err_opts('int solution',lineno)
      ALLOCATE(iout(n_iout),STAT=ierror)
      IF (ierror /= 0) CALL alloc_err(ierror)
      CALL mem_add(4*n_iout)

      ! Read all the values (iteration numbers) where intermediate output is
      ! desired.  Note that there may be multiple nodes per line, but they must
      ! be separated by one or more blank spaces.
      counter = 0
      DO WHILE (counter < n_iout)
        flag = readline(funit,line,lineno)  ! Read each line.
        IF (.NOT.flag) EXIT

        p = readp(line)
        DO WHILE (p > 0)  ! Read all the points in the same line.
          counter = counter + 1
          iout(counter) = p
          p = readp(line)
          IF (counter >= n_iout) EXIT
        END DO
      END DO
      max_iout = MAXVAL(iout)

    ! * * * * * Relaxation parameters for stage boundary condition * * * * * !
    ELSE IF (line(1:9) == "STAGEBC A") THEN
      READ(line(10:40),*) alpha
      alpha = MIN(alpha,one)
      IF (alpha <= 0) CALL err_opts('stagebc a',lineno)

    ELSE IF (line(1:9) == "STAGEBC K") THEN
      READ(line(10:40),*) kappa
      IF (kappa <= 0) CALL err_opts('stagebc k',lineno)

    ! * * * * * * * * * Steady vs. unsteady computations * * * * * * * * * !
    ELSE IF (line(1:13) == "TIME STEPPING") THEN
      READ(line(14:40),*) i
      IF (i <= 0) THEN
        steady = .TRUE.
      ELSE
        steady = .FALSE.
      END IF

    ! * * * * * * * * * * Criterion for time step size * * * * * * * * * * !
    ELSE IF (line(1:9) == "DELT TYPE") THEN
      READ(line(10:40),*) dt_op
      IF (dt_op < 0 .OR. dt_op > 1) CALL err_opts('delt type',lineno)

    ! * * * * * * * * * * * Order of Runge-Kutta scheme * * * * * * * * * * * !
    ELSE IF (line(1:8) == "RK ORDER") THEN
      READ(line(9:40),*) rkorder
      IF (rkorder < 1 .OR. rkorder > 3) CALL err_opts('rk order',lineno)

    ! * Option to read-in the mesh connectivity tables from an external file *
    ELSE IF (line(1:9) == "READ CONN") THEN
      read_conn = .TRUE.
      READ(line(10:40),'(A)') f_read_conn
      f_read_conn = ADJUSTL(f_read_conn)

    ! * Option to write the mesh connectivity tables to an external file * !
    ELSE IF (line(1:10) == "WRITE CONN") THEN
      write_conn = .TRUE.
      READ(line(11:40),'(A)') f_write_conn
      f_write_conn = ADJUSTL(f_write_conn)

    ! * * * * * * * * * Treatment of the source/sink terms * * * * * * * * * !
    ELSE IF (line(1:9) == "SRC TERMS") THEN
      opt = ADJUSTL(line(10:))
      IF (opt(1:3) == 'RDS') THEN
        opt_src = 0
      ELSE IF (opt(1:2) == 'FV') THEN
        opt_src = 1
      ELSE
        CALL err_opts('src terms',lineno)
      END IF

    ! * * * * * * * * Limiter in the gradient reconstruction * * * * * * * * !
    ELSE IF (line(1:7) == "LIMITER") THEN
      opt = ADJUSTL(line(8:))
      IF (opt(1:4) == 'ZERO') THEN
        opt_limiter = 0
      ELSE if (opt(1:4) == 'NONE') THEN
        opt_limiter = 1
      ELSE IF (opt(1:2) == 'BJ') THEN
        opt_limiter = 2
      ELSE
        CALL err_opts('limiter',lineno)
      END IF

    ! * * * * * * * * * * Parabolic eddy viscosity model * * * * * * * * * * !
    ELSE IF (line(1:12) == "PARABOLIC EV") THEN
      parabev = .TRUE.
      IF (LEN_TRIM(line(13:40)) > 0) THEN
        READ(line(13:40),*) omega
        omega = ABS(omega)
      END IF

    ! * * * * * * * * Higher order least squares gradients * * * * * * * * !
    ELSE IF (line(1:13) == "ACTIVATE HDLS") THEN
      activateHDLS = .TRUE.

    ! * * * * * * * * * * * * * Residual smoothing * * * * * * * * * * * * * !
    ELSE IF (line(1:14) == "ACTIVATE RSAVG") THEN
      activateRSAVG = .TRUE.
      IF (LEN_TRIM(line(15:40)) > 0) THEN
        READ(line(13:40),*) epsRSAVG,nRSAVG
        epsRSAVG = MIN(ABS(epsRSAVG),one)
      END IF

    ! * * * * * * * Set-up for wind-driven circulation terms * * * * * * * !
    ELSE IF (line(1:4) == "WIND") THEN
      wind_terms = .TRUE.
      READ(line(5:40),'(A)') wind_file
      wind_file = ADJUSTL(wind_file)
      IF (LEN_TRIM(wind_file) == 0) CALL err_opts('wind',lineno)
      INQUIRE (FILE=wind_file,EXIST=fnexist)
      IF (.NOT.fnexist) THEN
        WRITE (*,'(2A/)') 'ERROR: wind datafile not found, ', TRIM(wind_file)
        CALL err_opts('wind',lineno)
      END IF

    ! * * * * * * * * * Use equivalent depth at cell edges * * * * * * * * * !
    ELSE IF (line(1:10) == "DEL HEQUIV") THEN
      READ(line(11:40),*) delta_hequiv
      IF (delta_hequiv < zero .OR. delta_hequiv > one) &
        CALL err_opts('del hequiv',lineno)

    ! * * * * * * * * Used for shock detection at cell edges * * * * * * * * !
    ELSE IF (line(1:9) == "DEL SHOCK") THEN
      READ(line(10:40),*) delta_hshock
      IF (delta_hshock < zero .OR. delta_hshock > one) &
        CALL err_opts('del shock',lineno)

    ! * * * * * * * * * * Set-up for culvert computations * * * * * * * * * * !
    ELSE IF (line(1:4) == "CULV") THEN
      culvert = .TRUE.
      READ(line(5:40),'(A)') cvfile
      cvfile = ADJUSTL(cvfile)
      IF (LEN_TRIM(cvfile) == 0) CALL err_opts('culv',lineno)
      INQUIRE (FILE=cvfile,EXIST=fnexist)
      IF (.NOT.fnexist) THEN
        WRITE (*,'(2A/)') 'ERROR: culvert datafile not found, ', TRIM(cvfile)
        CALL err_opts('culv',lineno)
      END IF

    ! * * * * * * * * * File name of the probe output file * * * * * * * * * !
    ELSE IF (line(1:8) == "GAGES FN") THEN
      READ(line(9:40),'(A)') probe_file
      probe_file = ADJUSTL(probe_file)

    ! * * * * * * * * * File name of the probe output file * * * * * * * * * !
    ELSE IF (line(1:10) == "GAGES SKIP") THEN
      READ(line(11:40),*) probeskip
      IF (probeskip < 0) CALL err_opts('gages skip',lineno)

    ! * * * * * * * * * * * * * Probed grid nodes * * * * * * * * * * * * * !
    ELSE IF (line(1:5) == "GAGES") THEN
      READ(line(6:40),*) nprobe
      IF (nprobe <= 0) CALL err_opts('gages',lineno)
      ALLOCATE(probepts(nprobe),STAT=ierror)
      IF (ierror /= 0) CALL alloc_err(ierror)
      CALL mem_add(4*nprobe)

      ! Read all the values (node numbers) where probed output is desired.
      ! There may be multiple nodes per line, separated by one or more blank
      ! spaces.
      counter = 0
      DO WHILE (counter < nprobe)
        flag = readline(funit,line,lineno)  ! Read each line.
        IF (.NOT.flag) EXIT

        p = readp(line)
        DO WHILE (p > 0)  ! Read all the points in the same line.
          counter = counter + 1
          probepts(counter) = p
          p = readp(line)
          IF (counter >= nprobe) EXIT
        END DO
      END DO

    ! * * * * * * * * * * * * * * Obsolete feature * * * * * * * * * * * * * !
    ELSE IF (line(1:7) == "VEL BCS") THEN
      WRITE (*,'(/"WARNING: option VEL BCS in options file an obsolete ", &
        "feature.")')
      WRITE (*,'("VEL BCS value ignored in line",I4," of options file."/)') &
        lineno

    ! * * * * * * * * * * * * * * Error condition * * * * * * * * * * * * * * !
    ELSE
      WRITE (*,'("ERROR: unrecognized option in line",I4, &
        & " of options file.")') lineno
      CALL byebye('Program SToRM stopped.')
    END IF

  END DO

  IF (coef >= one) h_wet = coef*h_dry

  IF (nprobe > 0 .AND. LEN_TRIM(probe_file) == 0) THEN
    PRINT *,''
    PRINT *,'ERROR: undefined file name for gaging points output.'
    CALL byebye('Incomplete options file. Program SToRM stopped.')
  END IF

END SUBROUTINE read_options
