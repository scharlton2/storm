SUBROUTINE read_cgns(filename,title)
  USE parameters
  USE geometry
  USE dep_vars
  USE io
  USE vbc_arrays
  USE options
  USE constants
  IMPLICIT NONE
  INCLUDE "Header Files\cgnslib_f.h"
!!DEC$ ATTRIBUTES REFERENCE, C, VARYING :: cg_goto_f
!!DEC$ ATTRIBUTES REFERENCE, C, VARYING :: cg_array_read_f

!-----------------------------------------------------------------------------!
!                                                                             !
!  This is the data loader for CGNS data files passed from iRIC 2.0.  Upon    !
!  exit from this subroutine, the dummy argument variable 'title' contains    !
!  a title to be written to the output file.  The dummy argument variable     !
!  'filename' contains the name of the main input file (with extension        !
!  '.cgn' or '.cgns') that contains the grid, optional parameters, initial    !
!  conditions, and boundary conditions for the run.                           !
!                                                                             !
!  Francisco Simoes, 4 May 2012                                               !
!  Last updated (mm-dd-yyyy): 06-02-2017 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  CHARACTER (LEN=*), INTENT(IN) :: filename
  CHARACTER (LEN=*), INTENT(OUT) :: title

! Local variables.
  INTEGER :: celldim,datatype,edgefile,i,ierror,isize(3,3), &
             j,k,maxnin,maxqin,narrays,nbases,ncoord,ncount,nsections, &
             nuser_data,nzones,physdim
  INTEGER, ALLOCATABLE, DIMENSION(:) :: elements,itempvar,parentdata
  REAL(8), ALLOCATABLE, DIMENSION(:) :: work1,work2,x,y
  DOUBLE PRECISION :: u0,v0,wse0
  CHARACTER (LEN=32) :: basename,coordname,username,zonename,strg
  INTEGER, EXTERNAL :: findtt
  LOGICAL, EXTERNAL :: get_iounit,sortEdges

  TITLE = 'CGNS'

! Open CGNS file for read.
  CALL cg_open_f(filename,CG_MODE_MODIFY,findex,ierror)
  IF (ierror /= 0) THEN
    CALL cg_error_print_f()
    CALL byebye('SToRM ERROR: cannot open CGNS input file.')
  END IF

! Initialize iRIC CGNS file.
  CALL cg_iric_init_f(findex,ierror)
  IF (ierror /= 0) THEN
    CALL cg_error_print_f()
    CALL byebye('SToRM iRIC ERROR: cannot start CGNS file.')
  END IF

! Read computational grid.
  CALL cg_nbases_f(findex,nbases,ierror)
  DO i = 1,nbases
    CALL cg_base_read_f(findex,i,basename,celldim,physdim,ierror)
    IF (ierror /= 0) THEN
      CALL cg_error_print_f()
      CALL byebye('SToRM ERROR: cannot read CGNS base.')
    END IF
    IF (TRIM(basename) == 'iRIC') THEN  ! iRIC base name.
      CALL cg_nzones_f(findex,i,nzones,ierror)
      IF (ierror /= 0) THEN
        CALL cg_error_print_f()
        CALL byebye('SToRM ERROR: cannot find number of CGNS zones.')
      END IF
      DO j = 1,nzones
        CALL cg_zone_read_f(findex,i,j,zonename,isize,ierror)
        IF (ierror /= 0) THEN
          CALL cg_error_print_f()
          CALL byebye('SToRM ERROR: cannot read CGNS zone.')
        END IF
        IF (TRIM(zonename) == 'iRICZone') THEN  ! iRIC zone name.
          ! At this point it would be good to have a check to make sure that
          ! this zone is really for unstructured grids.

          ! Grid size.
          n_pts = isize(1,1)
          n_elems = isize(2,1)
          IF (n_pts < 1 .OR. n_elems < 1) THEN  ! Error checking.
            CALL cg_error_print_f()
            CALL byebye('SToRM CGNS ERROR: no points or nodes in CGNS file.')
          END IF

          ! Allocate space for some of the grid-related global variables.
          ALLOCATE(nodes(n_pts),STAT=ierror)  ! Nodal information.
          IF (ierror /= 0) CALL alloc_err(ierror)
          ALLOCATE(grid(n_elems+1),STAT=ierror)  ! Element connectivity.
          IF (ierror /= 0) CALL alloc_err(ierror)

          ! Set datav array for vertex data comming from the CGNS file.
          datav = .TRUE.
          datav(8)  = .FALSE.  ;  datav(9)  = .FALSE.
          datav(10) = .FALSE.  ;  datav(11) = .FALSE.

          ! Read nodal coordinates.
          ALLOCATE(x(n_pts),y(n_pts),STAT=ierror)  ! Local working arrays.
          IF (ierror /= 0) CALL alloc_err(ierror)
          CALL cg_ncoords_f(findex,i,j,ncoord,ierror)
          IF (ierror /= 0) THEN
            CALL cg_error_print_f()
          END IF
          DO k = 1,ncoord
            CALL cg_coord_info_f(findex,i,j,k,datatype,coordname,ierror)
            IF (ierror /= 0) THEN
              CALL cg_error_print_f()
              CALL byebye('SToRM CGNS ERROR: cannot read coord info.')
            END IF

            SELECT CASE(TRIM(coordname))
            CASE('CoordinateX')
              CALL cg_coord_read_f(findex,i,j,coordname,RealDouble,1, &
                                   n_pts,x,ierror)
              IF (ierror /= 0) THEN
                CALL cg_error_print_f()
                CALL byebye('SToRM CGNS ERROR: cannot read x-coords.')
              END IF
            CASE('CoordinateY')
              CALL cg_coord_read_f(findex,i,j,coordname,RealDouble,1, &
                                   n_pts,y,ierror)
              IF (ierror /= 0) THEN
                CALL cg_error_print_f()
                CALL byebye('SToRM CGNS ERROR: cannot read y-coords.')
              END IF
            END SELECT
          END DO
          DO k = 1,n_pts
            nodes(k)%x = x(k)
            nodes(k)%y = y(k)
          END DO

          ! Read element connectivity.
          ALLOCATE(elements(n_elems*3),parentdata(n_elems),STAT=ierror)
          IF (ierror /= 0) CALL alloc_err(ierror)
          CALL cg_nsections_f(findex,i,j,nsections,ierror)
          IF (ierror /= 0) THEN
            CALL cg_error_print_f()
            CALL byebye('SToRM CGNS ERROR: cannot read # sections.')
          END IF
          DO k = 1,nsections
            CALL cg_elements_read_f(findex,i,j,k,elements,parentdata, &
                                    ierror)
            IF (ierror /= 0) THEN
              CALL cg_error_print_f()
              CALL byebye('SToRM CGNS ERROR: cannot read connectivity table.')
            END IF
          END DO
          DO k = 1,n_elems
            grid(k)%vertex(1) = elements(k*3 - 2)
            grid(k)%vertex(2) = elements(k*3 - 1)
            grid(k)%vertex(3) = elements(k*3)
          END DO
          DEALLOCATE(elements,parentdata)

          ! Read bed elevation and roughness.
          ALLOCATE(z_tcp(n_pts),cd_tcp(n_pts),STAT=ierror) ! Bed elev.
          IF (ierror /= 0) CALL alloc_err(ierror)
          CALL cg_goto_f(findex,i,ierror,'Zone_t',j,'end')
          CALL cg_nuser_data_f(nuser_data,ierror)
          CALL cg_user_data_read_f(nuser_data,username,ierror)
          CALL cg_goto_f(findex,i,ierror,'Zone_t',j,'UserDefinedData_t',1, &
                         'end')
          CALL cg_goto_f(findex,i,ierror,'Zone_t',j,'UserDefinedData_t',1, &
                         'UserDefinedData_t',1,'end')
          CALL cg_narrays_f(narrays,ierror)
          CALL cg_array_read_f(1,z_tcp,ierror)  ! Read bed elevation.
          IF (ierror /= 0) THEN
            CALL cg_error_print_f()
            CALL byebye('SToRM CGNS ERROR: cannot read bed elevation.')
          END IF
          CALL cg_goto_f(findex,i,ierror,'Zone_t',j,'UserDefinedData_t',1, &
                         'UserDefinedData_t',2,'end')
          CALL cg_array_read_f(1,cd_tcp,ierror)  ! Read friction coefficient.
          IF (ierror /= 0) THEN
            CALL cg_error_print_f()
            CALL byebye('SToRM CGNS ERROR: cannot read friction coefficients.')
          END IF

          ! Read initial conditions from coverage polygons.
          ALLOCATE(u_tcp(n_pts),v_tcp(n_pts),h_tcp(n_pts),STAT=ierror)
          IF (ierror /= 0) CALL alloc_err(ierror)
          CALL cg_goto_f(findex,i,ierror,'Zone_t',j,'UserDefinedData_t',1, &
                         'UserDefinedData_t',3,'end')
          CALL cg_array_read_f(1,u_tcp,ierror)  ! Read x-component of velocity.
          IF (ierror /= 0) THEN
            CALL cg_error_print_f()
            CALL byebye('SToRM CGNS ERROR: cannot read initial U-velocity.')
          END IF
          CALL cg_goto_f(findex,i,ierror,'Zone_t',j,'UserDefinedData_t',1, &
                         'UserDefinedData_t',4,'end')
          CALL cg_array_read_f(1,v_tcp,ierror)  ! Read y-component of velocity.
          IF (ierror /= 0) THEN
            CALL cg_error_print_f()
            CALL byebye('SToRM CGNS ERROR: cannot read initial V-velocity.')
          END IF
          CALL cg_goto_f(findex,i,ierror,'Zone_t',j,'UserDefinedData_t',1, &
                         'UserDefinedData_t',5,'end')
          CALL cg_array_read_f(1,h_tcp,ierror)  ! Read water surface elevation.
          IF (ierror /= 0) THEN
            CALL cg_error_print_f()
            CALL byebye('SToRM CGNS ERROR: cannot read initial stage.')
          END IF
          DO k = 1,n_pts
            h_tcp(k) = MAX(zero,h_tcp(k) - z_tcp(k))
          END DO

          ! Read wind forcing terms from coverage polygons.
          ALLOCATE(w_fric(n_pts),w_mag(n_pts),w_dir(n_pts),STAT=ierror)
          IF (ierror /= 0) CALL alloc_err(ierror)
          CALL cg_goto_f(findex,i,ierror,'Zone_t',j,'UserDefinedData_t',1, &
                         'UserDefinedData_t',6,'end')
          CALL cg_array_read_f(1,w_fric,ierror)  ! Read wind friction.
          IF (ierror /= 0) THEN
            CALL cg_error_print_f()
            CALL byebye('SToRM CGNS ERROR: cannot read wind friction.')
          END IF
          CALL cg_goto_f(findex,i,ierror,'Zone_t',j,'UserDefinedData_t',1, &
                         'UserDefinedData_t',7,'end')
          CALL cg_array_read_f(1,w_mag,ierror)  ! Read wind magnitude.
          IF (ierror /= 0) THEN
            CALL cg_error_print_f()
            CALL byebye('SToRM CGNS ERROR: cannot read wind magnitude.')
          END IF
          CALL cg_goto_f(findex,i,ierror,'Zone_t',j,'UserDefinedData_t',1, &
                         'UserDefinedData_t',8,'end')
          CALL cg_array_read_f(1,w_dir,ierror)  ! Read wind direction (azimuth).
          IF (ierror /= 0) THEN
            CALL cg_error_print_f()
            CALL byebye('SToRM CGNS ERROR: cannot read wind direction.')
          END IF

        END IF  ! END iRIC zone name.
      END DO  ! END loop over zones.
    END IF  ! END iRIC base name.
  END DO  ! END loop over bases

! Read edge and vertex information.  Filenames are hardcoded. First read the
! .edge file.
  IF (.NOT. get_iounit(edgefile)) THEN
    WRITE (*,'(" ERROR: unable to open triangle .edge file.")')
    CALL byebye('SToRM stopped.')
  END IF
  OPEN(edgefile,FILE="iRICZone.1.edge",STATUS='OLD',IOSTAT=ierror)
  IF (ierror /= 0) THEN
  WRITE (*,'(" ERROR: unable to open file iRICZone.1.edge.")')
    CALL byebye('SToRM stopped.')
  END IF
  READ (edgefile,*) n_edges
  ! Allocate storage for global array edges.
  ALLOCATE(edges(n_edges),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  DO i = 1,n_edges
    READ (edgefile,*) j,edges(i)%p(1),edges(i)%p(2)
  END DO
  CLOSE (edgefile)

! Now read the .v.edge file.
  IF (.NOT. get_iounit(edgefile)) THEN
    WRITE (*,'(" ERROR: unable to open triangle .v.edge file.")')
    CALL byebye('SToRM stopped.')
  END IF
  OPEN(edgefile,FILE="iRICZone.1.v.edge",STATUS='OLD',IOSTAT=ierror)
  IF (ierror /= 0) THEN
  WRITE (*,'(" ERROR: unable to open file iRICZone.1.v.edge.")')
    CALL byebye('SToRM stopped.')
  END IF
  READ (edgefile,*) i
  IF (i /= n_edges) THEN
    WRITE (*,'(" ERROR: inconsistent data in .v.edge file.")')
    CALL byebye('SToRM stopped.')
  END IF
  DO i = 1,n_edges
    READ (edgefile,*) j,edges(i)%e(1),edges(i)%e(2)
    IF (edges(i)%e(1) < 0) THEN  ! Swap nodes.
      edges(i)%e(1) = edges(i)%e(2)
      edges(i)%e(2) = -1
    END IF
  END DO
  CLOSE (edgefile)

! Initial conditions.
  CALL cg_iric_read_integer_f("uniform_vinit",k,ierror)
  usehsfile = (k == 2)
  IF (k == 1) THEN
    ! Read uniform initial conditions.  If activated, these conditions superseed
    ! those defined by coverage polygons.
    CALL cg_iric_read_real_f("u_init",u0,ierror)
    CALL cg_iric_read_real_f("v_init",v0,ierror)
    CALL cg_iric_read_real_f("w_init",wse0,ierror)
    DO i = 1,n_pts  ! For a start with uniform conditions.
      u_tcp(i) = u0
      v_tcp(i) = v0
      h_tcp(i) = MAX(zero,wse0 - z_tcp(i))
    END DO

  ELSE IF (k == 2) THEN
    ! Read existing solution from file.  If activated, these conditions superseed
    ! those defined by coverage polygons.
    CALL cg_iric_read_string_f("startup_file",hsfile,ierror)
    CALL hot_start_cgns(hsfile,n_pts,h_tcp,u_tcp,v_tcp)
  END IF

! Set-up default on all options.
  CALL initialize_options

! Read solution parameters from CGNS file.
  steady = .FALSE.
  CALL cg_iric_read_integer_f("tstep_method",rkorder,ierror)  ! RK order.
  IF (ierror /= 0) THEN
    CALL cg_error_print_f()
    CALL byebye('SToRM CGNS ERROR: cannot read time-stepping order.')
  END IF
  CALL cg_iric_read_real_f("tstep_size",u0,ierror)  ! Time step size.
  IF (ierror /= 0 .OR. u0 < zero) THEN
    IF (ierror /= 0) CALL cg_error_print_f()
    CALL byebye('SToRM CGNS ERROR: invalid time step size.')
  END IF
  cfl = u0
  !CALL cg_iric_read_integer_f("tstep_numbr",istop,ierror)  ! # of time steps.
  CALL cg_iric_read_string_f("tstep_numbr",strg,ierror)
  IF (ierror /= 0 .OR. istop < 0) THEN
    IF (ierror /= 0) CALL cg_error_print_f()
    CALL byebye('SToRM CGNS ERROR: value of the number of time steps.')
  END IF
  strg = ADJUSTL(strg)
  istop = findtt(strg,cfl)
  IF (istop < 0) THEN
    WRITE (*,'("String in number of time steps: ",A)') TRIM(strg)
    CALL byebye("SToRM CGNS ERROR: format is not valid.")
  END IF
  !CALL cg_iric_read_integer_f("tstep_interval",k,ierror)  ! Output interval.
  CALL cg_iric_read_string_f("tstep_interval",strg,ierror)
  strg = ADJUSTL(strg)
  k = findtt(strg,cfl)
  IF (k < 0) THEN
    WRITE (*,'("String in interval for plotting: ",A)') TRIM(strg)
    CALL byebye("SToRM CGNS ERROR: format is not valid.")
  END IF
  IF (k < 1) k = MAX(istop,1)
  !IF (k > istop) k = istop
  n_iout = (istop + 1)/k

  ! This is to set-up variables for CGNS solution output.
  TSPrintCount = 0
  TSCellPrintCount = 0
  ALLOCATE(SolNames(n_iout*10+1),CellSolNames(n_iout*10+1), &
           TimeIncrements(n_iout*10+1),CellTimeIncrements(n_iout*10+1), &
           iout(n_iout),STAT=ierror)
  ncount = 0
  DO i = 1,istop
    IF (MOD(i,k) == 0) THEN
      ncount = ncount + 1
      iout(ncount) = i
    ENDIF
  ENDDO
  max_iout = MAXVAL(iout)
  k = 0
  CALL cg_iric_read_integer_f("tstep_large",k,ierror)  ! Velocity clipping.
  IF (k == 1) THEN
    vclip = .TRUE.
    CALL cg_iric_read_real_f("tstep_large_val",u0,ierror)
    IF (u0 < zero) &
      CALL byebye('SToRM CGNS ERROR: negative magnitude of maximum velocity.')
    vceiling = u0
  END IF

! Read wetting and drying parameters.
  CALL cg_iric_read_real_f("dry_thresh",u0,ierror)  ! Threshold for dry cells
  IF (ierror /= 0 .OR. u0 <= zero) THEN
    IF (ierror /= 0) CALL cg_error_print_f()
    CALL byebye('SToRM ERROR: value of threshold for dry cells must be > 0.')
  END IF
  h_dry = u0
  CALL cg_iric_read_real_f("wet_thresh",u0,ierror)  ! Threshold for wet cells.
  IF (ierror /= 0 .OR. u0 < one) THEN
    IF (ierror /= 0) CALL cg_error_print_f()
    CALL byebye('SToRM ERROR: value of threshold for wet cells must be > 1.')
  END IF
  u0 = MAX(u0,one)
  h_wet = h_dry*u0

! Read type of resistance coefficient used.
  CALL cg_iric_read_integer_f("friction_option",k,ierror)
  IF (ierror /= 0 .OR. k < 0 .OR. k > 2) THEN
    IF (ierror /= 0) CALL cg_error_print_f()
    CALL byebye('SToRM CGNS ERROR: friction coefficient has invalid type.')
  END IF
  opt_friction = k

! Activate culvert computations.
  CALL cg_iric_read_integer_f("culv_option",k,ierror)
  IF (ierror /= 0) THEN
    culvert = .FALSE.  ! Legacy from the previous version.
  ELSE
    culvert = .FALSE.
    IF (k == 1) THEN
      culvert = .TRUE.
      CALL cg_iric_read_string_f("culv_file",cvfile,ierror)
    END IF
  END IF

! ACtivate wind forcing terms.
  CALL cg_iric_read_integer_f("wind_option",k,ierror)
  IF (ierror /= 0) THEN
    wind_terms = .FALSE.  ! Legacy from the previous version.
  ELSE
    wind_terms = .FALSE.
    IF (k == 1) wind_terms = .TRUE.
  END IF

! Prepare edge information and data structure.
  CALL find_edges

! Compute triangle areas.
  CALL t_areas(grid,n_elems,nodes,n_pts,edges,n_edges)

!-----------------------------------------------------------------------------!
!                                                                             !
!  The following code replaces the code in SUBROUTINE ioflow.  Subroutine     !
!  ioflow is used when reading from text flat files, and the following code   !
!  does the same thing, but reading from a CGNS iRIC 2.0 data file.           !
!                                                                             !
!-----------------------------------------------------------------------------!

! Read boundary conditions (see SUBROUTINE ioflow).
! Inflow boundaries are read first.
  CALL cg_iric_read_bc_count_f('inflow',n_inflowbdr)
  IF (n_inflowbdr < 1) &
    CALL byebye('SToRM CGNS ERROR: no inflow boundaries found.')
  maxqin = 0
  DO i = 1,n_inflowbdr  ! Find inflow table with the most entries.
    CALL cg_iric_read_bc_functionalsize_f('inflow',i,'qin',k,ierror)
    maxqin = MAX(k,maxqin)
  END DO
  ALLOCATE(qin(n_inflowbdr),qin2(n_inflowbdr),vtype(n_inflowbdr), &
    n_timeserq(n_inflowbdr),timelocq(maxqin,n_inflowbdr),n_qin(n_inflowbdr), &
    timeserq(maxqin,n_inflowbdr),timeserq2(maxqin,n_inflowbdr),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
! Find the maximum number of nodes in the inflow node strings.
  maxnin = 0
  DO i = 1,n_inflowbdr
    CALL cg_iric_read_bc_indicessize_f('inflow',i,n_qin(i),ierror)
    IF (maxnin < n_qin(i)) maxnin = n_qin(i)
  END DO
  ALLOCATE(qin_nodes(maxnin,n_inflowbdr),itempvar(maxnin),work1(maxqin), &
    work2(maxqin),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  DO j = 1,n_inflowbdr  ! Initialize array to -1.  This may be inportant later,
    DO i = 1,maxnin     ! when using function in_array2 to check wall edges.
      qin_nodes(i,j) = -1
    END DO
  END DO

! Loop through all inflow boundaries to load inflow node strings and boundary
! condition values.
  DO i = 1,n_inflowbdr
    IF (n_qin(i) < 1) &
      CALL byebye('SToRM CGNS ERROR: inflow boundary has no nodes assigned.')
    ! Load node string:
    CALL cg_iric_read_bc_indices_f('inflow',i,itempvar,ierror)
    ! Load time-discharge table:
    CALL cg_iric_read_bc_functional_f('inflow',i,'qin',work1,work2,ierror)
    qin2(i) = zero
    call cg_iric_read_bc_integer_f ('inflow',i,'qin_type',k,ierror)
    IF (ierror /= 0) THEN
      vtype(i) = 1  ! Lecacy from the previous version.
    ELSE
      vtype(i) = k
    END IF
    IF (k /= 1 .AND. k /= 5) vtype(i) = 1  ! Error catch.
    CALL cg_iric_read_bc_functionalsize_f('inflow',i,'qin',k,ierror)
    n_timeserq(i) = k
    IF (k < 1) &
      CALL byebye('SToRM CGNS ERROR: inflow boundary has no discharge table.')
    DO j = 1,k
      timelocq(j,i) = work1(j)
      timeserq(j,i) = work2(j)
      timeserq2(j,i) = zero
    END DO
    IF (.NOT. sortEdges(n_qin(i),itempvar)) THEN
      PRINT *,itempvar
      PRINT *,"ERROR: inflow boundary edges not contiguous."
      CALL byebye("Program SToRM stopped.")
    END IF
    DO j = 1,n_qin(i)
      qin_nodes(j,i) = itempvar(j)
    END DO
  END DO
  DEALLOCATE(itempvar,work1,work2)

! Now read outflow boundaries.
  CALL cg_iric_read_bc_count_f('outflow',n_outflowbdr)
  IF (n_outflowbdr < 1) &
    CALL byebye('SToRM CGNS ERROR: no outflow boundaries found.')
  maxqin = 0
  DO i = 1,n_outflowbdr  ! Find outflow table with the most entries.
    CALL cg_iric_read_bc_functionalsize_f('outflow',i,'etaout',k,ierror)
    maxqin = MAX(k,maxqin)
  END DO
  ALLOCATE(hbc(n_outflowbdr),htype(n_outflowbdr),n_timeserh(n_outflowbdr), &
    timeloch(maxqin,n_outflowbdr),timeserh(maxqin,n_outflowbdr), &
    n_hbc(n_outflowbdr),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
! Find the maximum number of nodes in the outflow node strings.
  maxnin = 0
  DO i = 1,n_outflowbdr
    CALL cg_iric_read_bc_indicessize_f('outflow',i,n_hbc(i),ierror)
    IF (maxnin < n_hbc(i)) maxnin = n_hbc(i)
  END DO
  ALLOCATE(hbc_nodes(maxnin,n_outflowbdr),itempvar(maxnin),work1(maxnin), &
    work2(maxnin),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  DO j = 1,n_outflowbdr  ! Initialize array to -1.  This may be inportant later
    DO i = 1,maxnin      ! when using function in_array2 to check wall edges.
      hbc_nodes(i,j) = -1
    END DO
  END DO

! Loop through all outflow boundaries to load outflow node strings and boundary
! condition values.
  DO i = 1,n_outflowbdr
    IF (n_hbc(i) < 1) &
      CALL byebye('SToRM CGNS ERROR: outflow boundary has no nodes assigned.')
    ! Load node string:
    CALL cg_iric_read_bc_indices_f('outflow',i,itempvar,ierror)
    ! Load time-stage table:
    CALL cg_iric_read_bc_functional_f('outflow',i,'etaout',work1,work2,ierror)
    htype(i) = 1  ! Subcritical outlet.
    CALL cg_iric_read_bc_functionalsize_f('outflow',i,'etaout',k,ierror)
    n_timeserh(i) = k
    IF (k < 1) &
      CALL byebye('SToRM CGNS ERROR: outflow boundary has no stage table.')
    DO j = 1,k
      timeloch(j,i) = work1(j)
      timeserh(j,i) = work2(j)
    END DO
    IF (.NOT. sortEdges(n_hbc(i),itempvar)) THEN
      PRINT *,itempvar
      PRINT *,"ERROR: outflow boundary edges not contiguous."
      CALL byebye("Program SToRM stopped.")
    END IF
    DO j = 1,n_hbc(i)
      hbc_nodes(j,i) = itempvar(j)
    END DO
  END DO
  DEALLOCATE(itempvar,work1,work2)

! Allocate some of the arrays used to keep track of computed inflows and
! outflows during the time stepping of the solution.
  ALLOCATE(cqin(n_inflowbdr),cqin2(n_inflowbdr),cqin3(n_inflowbdr), &
    cqout(n_outflowbdr),cqout2(n_outflowbdr),cqout3(n_outflowbdr),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)

  DEALLOCATE(x,y)

!-----------------------------------------------------------------------------!
!                                                                             !
!  End of ioflow CGNS code.                                                   !
!                                                                             !
!-----------------------------------------------------------------------------!

! Find boundary edges.
  CALL walls

! Set-up arrays for inflow boundary condition computations.
  CALL ibc_arrays

! Set-up node data structure.
  CALL node_db

! Set-up element connectivity data structure.
  CALL elemnt_db

END SUBROUTINE read_cgns
