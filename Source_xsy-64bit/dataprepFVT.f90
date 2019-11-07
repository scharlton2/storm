SUBROUTINE dataprepFVT
  USE parameters
  USE geometry
  USE dep_vars
  USE options
  USE constants
  USE memory
  USE vbc_arrays
  USE RKparams
  USE io
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Prepares the data for use by the finite volume solver (FVT, as in Finite   !
!  Volume, Triangles).  It allocates needed array space, moves the variables  !
!  that need it from the vertices of the cell to the cell geometric centers,  !
!  and releases the memory used for temporary arrays.                         !
!                                                                             !
!  Francisco Simoes, February 2007                                            !
!  Last updated (mm-dd-yyyy): 05-02-2013 by F. Simões                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables.
  INTEGER :: aux,i,ierror,j,k,l
  REAL (KIND=mp) :: xc,yc,zc

! Allocate space for some of the global variables.
  ! Bed elevation.
  ALLOCATE(z(n_elems),zvtx(n_pts),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(8*(n_elems + n_pts))
  ! Solution variabless.
  ALLOCATE(h(n_elems+1),hvtx(n_pts),u(n_elems),uvtx(n_pts),v(n_elems), &
    vvtx(n_pts),u_mag(n_elems),wavec(n_elems),zeta(n_elems),zetavtx(n_pts), &
    cdvtx(n_pts),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(8*(6*n_elems + 4*n_pts) + 8)
  ALLOCATE(cd(n_elems),rough(n_elems),STAT=ierror)  ! Bed friction coefficient.
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(8*2*n_elems)
  ! Gradients.
  ALLOCATE(delh(n_elems),delu(n_elems),delv(n_elems),delz(n_elems),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(4*v_size*n_elems)
  ALLOCATE(ev(n_pts),STAT=ierror)  ! Turbulent eddy viscosity.
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(8*n_pts)
  ALLOCATE(phi(n_elems,3),STAT=ierror)  ! Cell residual.
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(8*3*n_elems)
  ALLOCATE(e_est(n_elems),STAT=ierror)  ! Error estimate.
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(8*n_elems)
  ! Dry cell flags.
  ALLOCATE(bdrycell(n_elems),drycell(n_elems),partdry(n_elems),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(3*4*n_elems)
  ALLOCATE(z2(n_elems),delz2(n_elems),STAT=ierror)  ! Used in 'bed_slope_ETA'.
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add((8 + v_size)*n_elems)

! If CGNS hot start is used, load cell center values from Tecplot ASCII file.
  IF (cgns .AND. usehsfile) CALL hot_start_ASCII

! Set-up arrays for interpolation from center of cells to vertices of cells
! (nodes).
  CALL c2v_setup(nctr2vtx)

! For the data that needs it, interpolate from cell-centered to cell-vertex
! positions.
  IF (datav(3)) THEN
    CALL vtx2ctr(z_tcp,z)
    DO i = 1,n_pts
      zvtx(i) = z_tcp(i)
    END DO
  ELSE
    DO i = 1,n_elems
      z(i) = z_tcp(i)
    END DO
    CALL ctr2vtx(nctr2vtx,z_tcp,zvtx)
    ! The bed elevations at the centers of the computational cells are
    ! calculated from the vertex values to ensure consistency and accuracy of
    ! the fluxes at the cell interfaces.  This means that output of z should
    ! be done at the vertices, not at the cell centers.
    CALL vtx2ctr(zvtx,z)
  ENDIF

  IF (datav(4)) THEN
    CALL vtx2ctr(h_tcp,h)
    DO i = 1,n_pts
      hvtx(i) = h_tcp(i)
      zetavtx(i) = hvtx(i) + zvtx(i)
    END DO
    CALL vtx2ctr(zetavtx,zeta)
  ELSE
    DO i = 1,n_elems
      h(i) = h_tcp(i)
    END DO
    ! Note: interpolation from cell centers to cell vertices should always be
    ! done using the water surface elevation, not the water depth.
    zeta = z + h
    CALL ctr2vtx(nctr2vtx,zeta,zetavtx)
    hvtx = zetavtx - zvtx
  ENDIF

  IF (datav(5)) THEN
    CALL vtx2ctr(u_tcp,u)
    DO i = 1,n_pts
      uvtx(i) = u_tcp(i)
    END DO
  ELSE
    DO i = 1,n_elems
      u(i) = u_tcp(i)
    END DO
    CALL ctr2vtx(nctr2vtx,u_tcp,uvtx)
  ENDIF

  IF (datav(6)) THEN
    CALL vtx2ctr(v_tcp,v)
    DO i = 1,n_pts
      vvtx(i) = v_tcp(i)
    END DO
  ELSE
    DO i = 1,n_elems
      v(i) = v_tcp(i)
    END DO
    CALL ctr2vtx(nctr2vtx,v_tcp,vvtx)
  ENDIF

  IF (datav(7)) THEN
    CALL vtx2ctr(cd_tcp,cd)
    cdvtx = cd_tcp
  ELSE
    DO i = 1,n_elems
      cd(i) = cd_tcp(i)
    END DO
    CALL ctr2vtx(nctr2vtx,cd_tcp,cdvtx)
  ENDIF

! Free the memory used for the Tecplot temporary arrays.
  DEALLOCATE(z_tcp,h_tcp,u_tcp,v_tcp,cd_tcp)

! Allocate storage variables for the bed elevation adjustment algorithm.
  ALLOCATE(zstore(n_elems),zvtxstore(n_pts),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(8*(n_elems + n_pts)*3)
  zstore = z
  zvtxstore = zvtx

! Reset datav() so that it can be used to write to Tecplot-formatted files.
  datav = .FALSE.
  datav(1) = .TRUE.  ;  datav(2) = .TRUE.
  datav(3) = .TRUE.
  ! To use datav(3) = .TRUE. properly, the output in subroutine
  ! 'tecplot_write' was correspondingly modified to print 'zvtx' instead of
  ! 'z'.

! Set-up vectors from center of cells to mid-points of edges.
  ALLOCATE(rc(n_elems,3),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(v_size*n_elems*3)
  DO i = 1,n_elems
    DO j = 1,3
      k = ABS(grid(i)%edge(j))  ! Edge being processed.
      ! Coordinates of the mid-point of edge k.
      xc = half*(nodes(edges(k)%p(1))%x + nodes(edges(k)%p(2))%x)
      yc = half*(nodes(edges(k)%p(1))%y + nodes(edges(k)%p(2))%y)
      rc(i,j)%x = xc - grid(i)%xc! - xc
      rc(i,j)%y = yc - grid(i)%yc! - yc
    END DO
  END DO

! Normalize normals and tangents in triangle database.
  DO i = 1,n_edges
    edges(i)%normal(1) = edges(i)%normal(1)/edges(i)%length
    edges(i)%normal(2) = edges(i)%normal(2)/edges(i)%length
    edges(i)%tang(1) = edges(i)%tang(1)/edges(i)%length
    edges(i)%tang(2) = edges(i)%tang(2)/edges(i)%length
  END DO

! For each triangle, find the element where the overdraft will be taken from in
! the drying cycles.
  DO i = 1,n_elems
    zc = z(i)
    grid(i)%o_cell = n_elems + 1  ! Catch-all ghost cell.
    DO j = 1,3
      k = ABS(grid(i)%edge(j))
      l = edges(k)%e(1)
      IF (l == i .OR. l < 0) l = edges(k)%e(2)
      IF (l == i .OR. l < 0) CYCLE
      ! Now l is the element connected to element i by edge k.
      IF (z(l) < zc) THEN
        zc = z(l)
        grid(i)%o_cell = l
      END IF
    END DO
  END DO

! The overdraft accounting method requires these two special storage locations
! to be properly initialized:
  h(n_elems+1) = zero
  grid(n_elems+1)%area = one

! ALLOCATE space for the quantities necessary to solve the 2D kinematic wave
! equation for shallow flows.
  IF (h_shallow > h_dry) THEN
    ALLOCATE(magzb(n_elems),ngzb(n_elems),STAT=ierror)
    IF (ierror /= 0) CALL alloc_err(ierror)
    CALL mem_add(n_elems*(v_size + 8))
  END IF

! If activated, ALLOCATE and prepare the arrays necessary to use a higher order
! least squares gradient reconstruction using all the first neighbors as
! computational molecule.  Otherwise, gradients are also calculated using
! least-squares reconstruction, but a smaller computational molecule is used.
  IF (activateHDLS) THEN
    CALL hdlsq_setup(1)
  ELSE
    CALL lsq_setup(opt_solver,2)
  END IF
  !CALL gaussgrad_setup

! If activated, ALLOCATE space for use by the residual smoothing operations.
  IF (activateRSAVG) THEN
    ALLOCATE(phi0(n_elems,3),STAT=ierror)  ! Unsmoothed cell residual.
    IF (ierror /= 0) CALL alloc_err(ierror)
    CALL mem_add(8*3*n_elems)
  END IF

! Read the data for the vegetation flow resistance model.
  IF (opt_friction == 5) CALL rvegdata

! Read the data for the wind forcing terms.  Do this for the text input only;
! CGNS input uses different quantities and procedure.
  IF (wind_terms) THEN
    IF (cgns) THEN
      CALL windforcing2
    ELSE
      CALL windforcing
    END IF
  END IF

! Set-up the Runge-Kutta time stepping procedure.  This follows a technique
! similar to the one in Cockburn (1998), pages 170-171.  Currently, there are
! three methods: the usual explicit Euler stepping (rkorder = 1), an optimal
! 2nd order SSPRK method (rkoder = 2), and an optimal 3rd order SSPRK method
! (rkorder = 3). See page 100 of Gottlieb et al. (2001).
  ALLOCATE(rkalpha(rkorder,0:rkorder-1),rkbeta(rkorder),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  ALLOCATE(rkh(n_elems,0:rkorder),rku(n_elems,0:rkorder), &
    rkv(n_elems,0:rkorder),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(8*n_elems*(rkorder + 1))
  SELECT CASE (rkorder)
  CASE (1)
    rkalpha(1,0) = one;  rkbeta(1) = one

  CASE (2)
    rkalpha(1,0) = one;                       rkbeta(1) = one
    rkalpha(2,0) = half; rkalpha(2,1) = half; rkbeta(2) = half

  CASE(3)
    rkalpha(1,0) = one
    rkalpha(2,0) = 0.750_mp; rkalpha(2,1) = 0.250_mp
    rkalpha(3,0) = one_third; rkalpha(3,1) = zero; rkalpha(3,2) = 2.0_mp/3.0_mp
    rkbeta(1) = one;  rkbeta(2) = 0.250_mp; rkbeta(3) = 2.0_mp/3.0_mp

  CASE DEFAULT
    CALL byebye('Incorrect Runge-Kutta order. SToRM stopped.')

  END SELECT

! Find the highest value of the bed elevation in the string of nodes belonging
! to the inflow boundary.  This information is used in a safety-net procedure
! for when all the inflow triangles get dry.
  bch_max = zero
  DO l = 1,n_bcedges
    IF (.NOT. qbar(l)) CYCLE
    k = bcedges(l)
    i = edges(k)%p(1)
    bch_max = MAX(bch_max,zvtx(i))
    i = edges(k)%p(2)
    bch_max = MAX(bch_max,zvtx(i))
  END DO

! For each triangle, build a table pointing to its vertices such that p(1) has
! the lowest elevation and p(3) has the highest.
  ALLOCATE(csortedz(n_elems),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(4*n_elems*3)
  DO i = 1,n_elems
    csortedz(i)%p(1) = grid(i)%vertex(1)
    csortedz(i)%p(2) = grid(i)%vertex(2)
    csortedz(i)%p(3) = grid(i)%vertex(3)
    DO j = 1,2  ! Two passes (worst case scenario).
      DO k = 1,2  ! Bubble-sort a list with 3 elements.
        IF (zvtx(csortedz(i)%p(k)) > zvtx(csortedz(i)%p(k+1))) THEN
          aux = csortedz(i)%p(k)  ! Swap elements.
          csortedz(i)%p(k) = csortedz(i)%p(k+1)
          csortedz(i)%p(k+1) = aux
        END IF
      END DO
    END DO
  END DO

! Do culvert preliminary computations.
  IF (culvert) CALL cvprep

! Build an array with the partial areas needed for the interpolation procedure
! used to compute the effective area of partially wet cells.
  !ALLOCATE(areap(n_elems,2),STAT=ierror)
  !IF (ierror /= 0) CALL alloc_err(ierror)
  !CALL mem_add(8*n_elems*2)
  !DO i = 1,n_elems

  !END DO

! Set-up arrays for one-dimensional channels in the computational domain.
  CALL chan_arrays
! Print-out for debugging SUBROUTINE chan_arrays.
  !DO j = 1,n_chanbdr
  !  WRITE (*,'(/A,I4)')"Channel number",j
  !  WRITE (*,'(18("-"))')
  !  WRITE (*,'(/A)')"        #     Node     Edge Triangle     Param t"
  !  DO i = 1,n_channel(j) - 1
  !    WRITE (*,'(4I9,2X,ES12.5)') i,channel(i,j),chanedges(i,j), &
  !      chantrigs(i,j),chant(i,j)
  !  END DO
  !  WRITE (*,'(2I9)') n_channel(j),channel(n_channel(j),j)
  !END DO
  ! Print-out coordinates of the channel nodes, useful for plotting in Tecplot.
  !DO j = 1,n_chanbdr
  !  DO i = 1,n_channel(j)
  !    WRITE (*,'(2I9,2(2X,ES12.5))') i,channel(i,j),nodes(channel(i,j))%x, &
  !      nodes(channel(i,j))%y
  !  END DO
  !END DO
  !CALL byebye('Stopped in dataprepFVT for debugging.')

END SUBROUTINE dataprepFVT
