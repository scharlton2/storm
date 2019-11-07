SUBROUTINE header(id,filename,title)
  USE options
  USE geometry
  USE dep_vars
  USE vbc_arrays
  USE io
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Main print-out subroutine for the run parameters before the solution run.  !
!                                                                             !
!  Francisco Simoes, March 2004                                               !
!  Last updated (mm-dd-yyyy): 05-13-2014 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments:
  INTEGER, INTENT(IN) :: id  ! Unit number for the external WRITE statements.
  CHARACTER (LEN=*), INTENT(IN) :: filename  ! Main filename, used as case id.
  CHARACTER (LEN=*), INTENT(OUT) :: title  ! Case id for output to Tecplot.

! Local variables:
  INTEGER :: i,ierror,j,k,value(8)
  CHARACTER :: hour*10,today*8,zone*5  ! Used in intrinsic subroutine.
  CHARACTER (LEN=2) :: day
  CHARACTER (LEN=3) :: month(12),res(3),rds(4)
  CHARACTER (LEN=4) :: what(2),year
  CHARACTER (LEN=5) :: limitr(3)
  CHARACTER (LEN=9) :: c(5),req(6),stype(3),tm(2),walltr(2)
  CHARACTER (LEN=80) :: buffer,fn,fntrim

  INTEGER, EXTERNAL :: mem_call,query
  LOGICAL, EXTERNAL :: get_iounit
  INTERFACE
    INTEGER FUNCTION npoints(varno,datav,n_pts,n_elems)
    INTEGER, INTENT(IN) :: varno,n_elems,n_pts
    LOGICAL, DIMENSION(:), INTENT(IN) :: datav
    END FUNCTION npoints
  END INTERFACE

  limitr = (/' ZERO',' NONE','BARTH'/)
  month = (/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct', &
            'Nov','Dec'/)
  req = (/"Manning's",". . Chezy","Drag Coef","Colebrook","Strickler","Vegetated"/)
  res = (/" FE","AFE","RDS"/)
  rds = (/". N","LDA",". B","PSI"/)
  stype = (/"Res Distr","Fin Vol T","Fin Vol Z"/)
  what = (/" OFF","  ON"/)
  walltr = (/" ZERO VEL","ZERO FLUX"/)
  tm = (/" UNSTEADY",".  STEADY"/)

! Program header.
  IF (output_head) THEN
    WRITE (id,'(A)')" ____  _____  ___   ____   __  __"
    WRITE (id,'(A)')"/ ___||_   _|/ _ \ |  _ \ |  \/  |"
    WRITE (id,'(A)')"\___ \  | | | | | || |_) )| |\/| |"
    WRITE (id,'(A)')" ___) ) | | | |_| ||  _ < | |  | |"
    WRITE (id,'(A)')"|____/  |_|  \___/ |_| \_\|_|  |_|"
    WRITE (id,*)
    WRITE (id,'(A)')"System for Transport and River Modeling [Version 0.3]"
    WRITE (id,'(A)')"U.S. Geological Survey                  [Track   xsy]"
    WRITE (id,'(A)')"National Research Program"  !           [Variant xsy]
    WRITE (id,'(A)')"Golden, CO 80403 - U.S.A."  !           [Tag     xsy]
    CALL DATE_AND_TIME(today,hour,zone,value)
    WRITE (id,'(/A,I2,1X,A,1X,I4)')"Run started at "//hour(1:2)//":"// &
        & hour(3:4)//":"//hour(5:6)//" on ",value(3),month(value(2)),value(1)
    WRITE (id,'(A,A)')"Run data file is ",TRIM(filename)
  END IF

! Data for output later in Tecplot file.
  WRITE (day,'(I2)') value(3)
  WRITE (year,'(I4)') value(1)
  title = "SToRM run [" // TRIM(filename) // "] started " // hour(1:2) // ":" &
    // hour(3:4) // ":" // hour(5:6) // " on " // day // "-" // &
    month(value(2)) // "-" // year

! Regurgitate input file.
  IF (output_file .AND. .NOT. cgns) THEN
    WRITE (id,'(/A/21("-"))')"Contents of run file:"
    IF (.NOT. get_iounit(i)) THEN
      WRITE (id,'(" ERROR: unable to open file ",A)') TRIM(filename)
      WRITE (id,'(" SToRM stopped.")')
      CALL byebye('SToRM stopped.')
    END IF
    OPEN (i,FILE=filename,STATUS='OLD',IOSTAT=ierror)
    IF (ierror /= 0) THEN
      WRITE (id,'(" ERROR: unable to open file ",A)') TRIM(filename)
      WRITE (id,'(" SToRM stopped.")')
      CALL byebye('SToRM stopped.')
    END IF
    DO
      READ (i,'(A)',IOSTAT=ierror) buffer
      IF (ierror /= 0) EXIT
      WRITE (id,'(A)') TRIM(buffer)
    END DO
    CLOSE (i)
  END IF

! Options set for the present run.
  IF (output_opts) THEN
    WRITE (id,'(/A/14("-"))')"Print options:"
    WRITE (id,'(A,31(" ."),A)')"Print run file",what(query(output_file))
    WRITE (id,'(A,32(" ."),A)')"Print header",what(query(output_head))
    WRITE (id,'(A,31(" ."),A)')"Print options.",what(query(output_opts))
    WRITE (id,'(A,26(" ."),A)')"Print computational grid", &
      what(query(output_geom))
    WRITE (id,'(A,26(" ."),A)')"Print initial conditions", &
      what(query(output_icnd))

    WRITE (id,'(/A/16("-"))')"Problem options:"
    WRITE (id,'(A,23(" ."),ES12.5)')"Water temperature (oC)",temp
    WRITE (id,'(A,26(" ."),ES12.5)')"Water viscosity.",visc
    WRITE (id,'(A,26(" ."),1X,A9)')"Roughness equation",req(opt_friction + 1)
    WRITE (id,'(A,29(" ."),1X,A9)')"Solver type ", stype(opt_solver)
    SELECT CASE (opt_solver)
    CASE (1)  ! Residual distribution scheme options.
      WRITE (id,'(A,26(" ."),1X,A3)')"Residual computation by.", &
        res(opt_residual + 1)
      WRITE (id,'(A,24(" ."),1X,A3)')"Residual distribution scheme", &
        rds(opt_rdscheme + 1)
      WRITE (id,'(A,29(" ."),A)')"Use BGN algorithm.",what(query(zadjust))

    CASE (2)  ! Finite volume method options.
      WRITE (id,'(A,29(" ."),1X,A5)')"Limiter function",limitr(opt_limiter+1)
      IF (parabev) &
        WRITE (id,'(A,22(" ."),ES12.5)')"Eddy viscosity constant.",omega
      WRITE (id,'(A,18(" ."),A)')"Use higher-order least squares gradients", &
        what(query(activateHDLS))
      WRITE (id,'(A,21(" ."),ES12.5)')"Shock detection threshold.",delta_hshock
      WRITE (id,'(A,21(" ."),ES12.5)')"Equivalent depth threshold",delta_hequiv

    END SELECT
    WRITE (id,'(A,25(" ."),ES12.5)')"Dry cell threshold",h_dry
    WRITE (id,'(A,25(" ."),ES12.5)')"Wet cell threshold",h_wet
    !WRITE (id,'(A,22(" ."),ES12.5)')"Shallow depth threshold.",h_shallow
    WRITE (id,'(A,21(" ."),ES12.5)')"Stagnation cell threshold.",u_stag
    IF (vclip) &
      WRITE (id,'(A,21(" ."),ES12.5)')"Maximum velocity clipping.",vceiling
    WRITE (id,'(A,27(" ."),ES12.5)')"CFL parameter.",cfl
    IF (steady) &
      WRITE (id,'(A,21(" ."),I4)')"Steady state convergence criterion",stop_cr
    WRITE (id,'(A,22(" ."),I8)')"Maximum number of iterations",istop
    WRITE (id,'(A,19(" ."),ES12.5)')"Cut-off value of the residual.",rstop
    WRITE (id,'(A,16(" ."),1X,A)')"Type of solid-walls boundary treatment", &
      walltr(btype + 1)
    WRITE (id,'(A,26(" ."),ES12.5)')"Stage bc: alpha.",alpha
    WRITE (id,'(A,26(" ."),ES12.5)')"Stage bc: kappa.",kappa
    WRITE (id,'(A,28(" ."),A)')"No-inflow at h-nodes",what(query(hvalve))
    WRITE (id,'(A,28(" ."),I4)')"Time stepping scheme",rkorder
    WRITE (id,'(A,19(" ."),1X,A)')"Type of time-stepping treatment.", &
      tm(query(steady))
    IF (.NOT. steady) WRITE (id,'(A,26(" ."),I4)')"Time step size criterion", &
      dt_op
    WRITE (id,'(A,17(" ."),I4)')"Type of treatment of the source/sink terms", &
      opt_src
    WRITE (id,'(A,27(" ."),A)')"Use residual smoothing", &
      what(query(activateRSAVG))
    IF (activateRSAVG) THEN
      WRITE (id,'(A,20(" ."),ES12.5)')"Residual smoothing constant.",epsRSAVG
      WRITE (id,'(A,23(" ."),I4)')"Number of smoothing iterations",nRSAVG
    END IF
    WRITE (id,'(A,19(" ."),I4)')"Type of center-to-vertex interpolation",nctr2vtx

    ! Read/write mesh connectivity using external files.
    WRITE (id,'(A,17(" ."),A)')"Read mesh connectivity from external file.",&
      what(query(read_conn))
    IF (read_conn) WRITE (id,'(A,A,A)')"(Mesh connectivity file name: ", &
      TRIM(f_read_conn),")"
    WRITE (id,'(A,18(" ."),A)')"Write mesh connectivity to external file", &
      what(query(write_conn))
    IF (write_conn) WRITE (id,'(A,A,A)')"(Mesh connectivity file name: ", &
      TRIM(f_write_conn),")"

    ! Write where intermediate solution steps will be printed.
    IF (n_iout > 0) THEN
      WRITE (id,'(/5X,A)')"Intermediate solution print-out"
      WRITE (id,'(8X,A)')"N     ITER   FILE NAME"
      WRITE (id,'(40("-"))')
      fntrim = filename
      CALL strip_ext(fntrim)
      DO i = 1,n_iout
        CALL ifname(iout(i),max_iout,fntrim,fn)
        WRITE (id,'(2I9,3X,A)') i,iout(i),TRIM(fn)
      END DO
    END IF

    ! Write the probed values for output.
    IF (nprobe > 0) THEN
      WRITE (id,'(/A)')"Location of gaging stations:"
      WRITE (id,'(28("-"))')
      DO i = 1,nprobe
        WRITE (id,'(I10)') probepts(i)
      END DO
      WRITE (id,'(" Output filename: ",A)') TRIM(probe_file)
    END IF

  END IF

! Geometry (i.e., computational grid) data base.  Supports up to 99,999,999
! data points (limited by the formatting instructions).
  IF (output_geom) THEN
    WRITE (id,'(/10X,A)')"Node coordinates"
    WRITE (id,'(7X,A)')"N          X                Y"
    WRITE (id,'(43("-"))')
    DO i = 1,n_pts
      WRITE (id,'(I8,2(2X,ES15.7))') i,nodes(i)
    END DO

    WRITE (id,'(/22X,A)')"Control volume connectivity table"
    WRITE (id,'(A)')"        N       N1       N2       N3       E1       &
    &E2       E3       AREA"
    WRITE (id,'(78("-"))')
    DO i = 1,n_elems
      WRITE (id,'(7I9,2X,ES12.5)') i,grid(i)%vertex(1),grid(i)%vertex(2), &
        grid(i)%vertex(3),grid(i)%edge(1),grid(i)%edge(2),grid(i)%edge(3), &
        grid(i)%area
    END DO

    WRITE (id,'(/31X,A)')"Edge data base"
    WRITE (id,'(A)')"       N      P1      P2      E1      E2      NX    " // &
      "      NY         LENGTH"
    WRITE (id,'(77("-"))')
    DO i = 1,n_edges
      WRITE (id,'(5I8,3(2X,ES10.3))') i,edges(i)%p(1),edges(i)%p(2), &
        edges(i)%e(1),edges(i)%e(2),edges(i)%normal(1),edges(i)%normal(2), &
        edges(i)%length
    END DO

    WRITE (id,'(/17X,A)')"Edge arrays"
    WRITE (id,'(A)')" FLOWEDGE FLOWEDGE WALLEDGE WALLEDGE BOUNDARY"
    WRITE (id,'(46("-"))')
    DO i = 1,MAX(n_bpolygon,flow_edges,wall_edges,flow_edges1,wall_edges1)
      c(1) = "        -"
      IF (i <= flow_edges) WRITE (c(1),'(I9)') flowedg(i)
      c(2) = "        -"
      IF (i <= flow_edges1) WRITE (c(2),'(I9)') flowedg1(i)
      c(3) = "        -"
      IF (i <= wall_edges) WRITE (c(3),'(I9)') walledg(i)
      c(4) = "        -"
      IF (i <= wall_edges1) WRITE (c(4),'(I9)') walledg1(i)
      c(5) = "        -"
      IF (i <= n_bpolygon) WRITE (c(5),'(I9)') bpolygon(i)
      WRITE (id,'(5A9)') c(1),c(2),c(3),c(4),c(5)
    END DO

    WRITE (id,'(/35X,A)')"Wall nodes"
    WRITE (id,'(80("-"))')
    i = 1
    DO
      IF (i > n_wall) EXIT
      j = MIN(i+9,n_wall)
      WRITE (id,'(10I8)')(wall_pts(k),k=i,j)
      i = j + 1
    END DO

    k = SIZE(n2t,2)
    WRITE (id,'(/A,I2,A)')"Node-to-triangle connectivity table (",k - 1,")"
    WRITE (id,'(A)')"     N  # of Ts  N2T:"
    WRITE (id,'(77("-"))')
    DO i = 1,n_pts
      WRITE (id,'(11I7)') i,(n2t(i,j),j=1,MIN(n2t(i,1) + 1,10))
    END DO

    k = SIZE(n2t2,2)
    WRITE (id,'(/A,I2,A)')"Node-to-triangle connectivity table (",k - 1,")"
    WRITE (id,'(A)')"     N  # of Ts N2T2:"
    WRITE (id,'(77("-"))')
    DO i = 1,n_pts
      WRITE (id,'(11I7)') i,(n2t2(i,j),j=1,MIN(n2t2(i,1) + 1,10))
    END DO

    k = SIZE(n2n,2)
    WRITE (id,'(/A,I2,A)')"Node-to-node connectivity table (",k - 1,")"
    WRITE (id,'(A)')"     N  # of Ns  N2N:"
    WRITE (id,'(77("-"))')
    DO i = 1,n_pts
      WRITE (id,'(11I7)') i,(n2n(i,j),j=1,MIN(n2n(i,1) + 1,10))
    END DO

    WRITE (id,'(/A)')"  Element-to-element connectivity table"
    WRITE (id,'(A)')"        N    T2T3(1)   T2T3(2)   T2T3(3)"
    WRITE (id,'(41("-"))')
    DO i = 1,n_elems
      WRITE (id,'(4I10)') i,(t2t3(j,i),j=1,3)
    END DO

    k = SIZE(t2t,2)
    WRITE (id,'(/A,I2,A)')"Element-to-element connectivity table (",k - 1,")"
    WRITE (id,'(A,I2,A)')"[For least-squares computational molecule]"
    WRITE (id,'(A)')"     N  # of Ns  T2T:"
    WRITE (id,'(77("-"))')
    DO i = 1,n_elems
      WRITE (id,'(11I7)') i,(t2t(i,j),j=1,MIN(t2t(i,1) + 1,10))
    END DO

    IF (ALLOCATED(t2tHD)) THEN
      k = SIZE(t2tHD,2)
      WRITE (id,'(/A,I2,A)')"Element-to-element connectivity table (",k - 1,")"
      WRITE (id,'(A,I2,A)')"[For higher order least-squares computational &
        &molecule]"
      WRITE (id,'(A)')"       N    # of Ns   T2THD:"
      WRITE (id,'(77("-"))')
      DO i = 1,n_elems
        WRITE (id,'(30I8)') i,(t2tHD(i,j),j=1,t2tHD(i,1) + 1)
        !WRITE (id,'(30I8)') i,(t2tHD(i,j),j=1,MIN(t2tHD(i,1) + 1,10))
      END DO
    END IF

    WRITE (id,'(/9X,A)')"Control volume size"
    WRITE (id,'(A)')"        N      AREA        PERIMETER"
    WRITE (id,'(38("-"))')
    DO i = 1,n_pts
      WRITE (id,'(I9,2(2X,ES12.5))') i,cv_area(i),cv_perim(i)
    END DO

  END IF

! Solution variables (initial conditions).
  IF (output_icnd) THEN

    ! Print inflow boundary data.
    DO j = 1,n_inflowbdr
      WRITE (id,'(/A,I4)')"Inflow boundary number",j
      WRITE (id,'(26("-"))')

      ! Print the line with the varaible names.
      SELECT CASE (vtype(j))
      CASE (0)  ! Normal velocity.
        WRITE (id,'(/A)')"Boundary conditions table"
        WRITE (id,'(A)')"        N       Time        Velocity"

      CASE (1,5)  ! Discharge.
        WRITE (id,'(/A)')"Boundary conditions table"
        WRITE (id,'(A)')"        N       Time       Discharge"

      CASE (3)  ! Full velocity vector.
        WRITE (id,'(/A)')"Boundary conditions table"
        WRITE (id,'(A)')"       N        Time           U             V"

      CASE (4)  ! Stage at both ends.
        WRITE (id,'(/A)')"Boundary conditions table"
        WRITE (id,'(A)')"        N       Time         Stage"
      END SELECT

      IF (vtype(j) == 3) THEN
        DO i = 1,n_timeserq(j)
          WRITE (id,'(I9,3(2X,ES12.5))') i,timelocq(i,j),timeserq(i,j), &
                                         timeserq2(i,j)
        END DO
      ELSE
        DO i = 1,n_timeserq(j)
          WRITE (id,'(I9,3(2X,ES12.5))') i,timelocq(i,j),timeserq(i,j)
        END DO
      END IF

      WRITE (id,'(/A)')"Inflow boundary node numbers:"
      DO i = 1,n_qin(j)
        WRITE (id,'(I10)') qin_nodes(i,j)
      END DO
    END DO

    ! Print outflow boundary data.
    DO j = 1,n_outflowbdr
      WRITE (id,'(/A,I4)')"Outflow boundary number",j
      WRITE (id,'(27("-"))')

      ! Print the line with the varaible names.
      SELECT CASE (htype(j))
      CASE (1)  ! Subcritical outflow.
        WRITE (id,'(/A)')"Boundary conditions table"
        WRITE (id,'(A)')"        N       Time         Stage"
        DO i = 1,n_timeserh(j)
          WRITE (id,'(I9,3(2X,ES12.5))') i,timeloch(i,j),timeserh(i,j)
        END DO

      CASE (0)
        WRITE (id,'(/A)')"Free flow boundary."

      CASE (2)
        WRITE (id,'(/A)')"Critical flow boundary (Fr = 1)."
      END SELECT

      WRITE (id,'(/A)')"Outflow boundary node numbers:"
      DO i = 1,n_hbc(j)
        WRITE (id,'(I10)') hbc_nodes(i,j)
      END DO
    END DO

    ! Print channel data.
    !DO j = 1,n_chanbdr
    !  WRITE (id,'(/A,I4)')"Channel number",j
    !  WRITE (id,'(18("-"))')
    !  WRITE (id,'(/A)')"        #     Node     Edge Triangle     Param t"
    !  DO i = 1,n_channel(j) - 1
    !    WRITE (id,'(4I9,2X,ES12.5)') i,channel(i,j),chanedges(i,j), &
    !      chantrigs(i,j),chant(i,j)
    !  END DO
    !  WRITE (id,'(2I9)') n_channel(j),channel(i,n_channel(j))
    !END DO

    WRITE (id,'(/31X,A)')"Initial conditions"
    WRITE (id,'(8X,"N",7X,"U",13X,"V",13X,"H",13X,"Z",13X,"CD")')
    WRITE (id,'(80("-"))')
    ! NOTE: the following DO-loop assumes that all the dependent variables have
    ! the same cardinality of u().  For now, this solves the problem of
    ! printing vertex or cell-centered dependent variables.
    DO i = 1,npoints(5,datav,n_pts,n_elems)
      WRITE (id,'(I9,5(2X,ES12.5))') i,u_tcp(i),v_tcp(i),h_tcp(i),z_tcp(i),cd_tcp(i)
    END DO

  END IF

END SUBROUTINE header
