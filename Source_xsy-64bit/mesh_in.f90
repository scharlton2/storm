SUBROUTINE mesh_in(filename)
  USE geometry
  USE dep_vars
  USE vbc_arrays
  USE memory
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  All the mesh connectivity arrays are read in this suroutine, except those  !
!  set-up in subroutine ioflow (for programming convenience reasons).         !
!  Subroutine ioflow must be CALLed before or after a CALL to mesh_in.  The   !
!  data is read from file filename in binary format.  The format and read     !
!  order must be identical to the ones in subroutine mesh_out.                !
!                                                                             !
!  Francisco Simoes, October 2006                                             !
!  Last updated (mm-dd-yyyy): 03-17-2011 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  CHARACTER (LEN=*) :: filename

! Local variables.
  INTEGER :: dim1,dim2,f_unit,i,ierror,j,k

  LOGICAL, EXTERNAL :: get_iounit

  IF (.NOT. get_iounit(f_unit)) THEN
    PRINT *,'ERROR in get_iounit().'
    CALL byebye('Program SToRM stopped.')
  END IF

  OPEN (f_unit,FILE=filename,POSITION='REWIND',FORM='UNFORMATTED')

  READ (f_unit) (grid(i),i=1,n_elems)  ! Allocated in Tecplot.

! From 'find_edges'.
  READ (f_unit) n_edges
  IF (n_edges > 0) THEN
    IF (.NOT. ALLOCATED(edges)) THEN
      ALLOCATE(edges(n_edges),STAT=ierror)
      IF (ierror /= 0) CALL alloc_err(ierror)
    END IF
    CALL mem_add(e_size*n_edges)
    READ (f_unit) (edges(i),i=1,n_edges)
  END IF

! From 'walls'.
  READ (f_unit) wall_edges,flow_edges,wall_edges1,flow_edges1,n_bpolygon,n_wall
  IF (wall_edges > 0) THEN
    ALLOCATE(walledg(wall_edges),STAT=ierror)
    IF (ierror /= 0) CALL alloc_err(ierror)
    CALL mem_add(4*wall_edges)
    READ (f_unit) (walledg(i),i=1,wall_edges)
  END IF
  IF (flow_edges > 0) THEN
    ALLOCATE(flowedg(flow_edges),STAT=ierror)
    IF (ierror /= 0) CALL alloc_err(ierror)
    CALL mem_add(4*flow_edges)
    READ (f_unit) (flowedg(i),i=1,flow_edges)
  END IF
  IF (wall_edges1 > 0) THEN
    ALLOCATE(walledg1(wall_edges1),STAT=ierror)
    IF (ierror /= 0) CALL alloc_err(ierror)
    CALL mem_add(4*wall_edges1)
    READ (f_unit) (walledg1(i),i=1,wall_edges1)
  END IF
  IF (flow_edges1 > 0) THEN
    ALLOCATE(flowedg1(flow_edges1),STAT=ierror)
    IF (ierror /= 0) CALL alloc_err(ierror)
    CALL mem_add(4*flow_edges1)
    READ (f_unit) (flowedg1(i),i=1,flow_edges1)
  END IF
  IF (n_bpolygon > 0) THEN
    ALLOCATE(bpolygon(n_bpolygon),STAT=ierror)
    IF (ierror /= 0) CALL alloc_err(ierror)
    CALL mem_add(4*n_bpolygon)
    READ (f_unit) (bpolygon(i),i=1,n_bpolygon)
  END IF
  IF (n_wall > 0) THEN
    ALLOCATE(wall_pts(n_wall),STAT=ierror)
    IF (ierror /= 0) CALL alloc_err(ierror)
    CALL mem_add(4*n_wall)
    READ (f_unit) (wall_pts(i),i=1,n_wall)
  END IF

! From 'ibc_arrays'.
  READ (f_unit) n_inflowbdr,dim1,n_outflowbdr,dim2,n_wall,wall_edges1,n_bcedges
  IF (n_inflowbdr > 0) THEN
    ALLOCATE(n_qin(n_inflowbdr),STAT=ierror)
    IF (ierror /= 0) CALL alloc_err(ierror)
    CALL mem_add(8*n_inflowbdr)
    READ (f_unit) (n_qin(i),i=1,n_inflowbdr)
    ALLOCATE(pnormals(dim1,n_inflowbdr),STAT=ierror)
    IF (ierror /= 0) CALL alloc_err(ierror)
    CALL mem_add(v_size*dim1*n_inflowbdr)
    ! Arrays used in FUNCTION vbc_by_h are allocated here.
    ALLOCATE(l_inflow(dim1,n_inflowbdr),h_inflow(dim1),a_inflow(dim1), &
      u_inflow(dim1),v_inflow(dim1),STAT=ierror)
    IF (ierror /= 0) CALL alloc_err(ierror)
    CALL mem_add(8*5*n_inflowbdr)
    READ (f_unit) ((l_inflow(i,j),i=1,dim1),j=1,n_inflowbdr)
    READ (f_unit) ((pnormals(i,j),i=1,dim1),j=1,n_inflowbdr)
  END IF
  IF (n_outflowbdr > 0) THEN
    ALLOCATE(n_hbc(n_outflowbdr),STAT=ierror)
    IF (ierror /= 0) CALL alloc_err(ierror)
    CALL mem_add(8*n_outflowbdr)
    READ (f_unit) (n_hbc(i),i=1,n_outflowbdr)
    ALLOCATE(hnormals(dim2,n_outflowbdr),STAT=ierror)
    IF (ierror /= 0) CALL alloc_err(ierror)
    CALL mem_add(v_size*dim2*n_outflowbdr)
    READ (f_unit) ((hnormals(i,j),i=1,dim2),j=1,n_outflowbdr)
  END IF
  IF (n_wall > 0) THEN
    ALLOCATE(wtangs(n_wall),STAT=ierror)
    IF (ierror /= 0) CALL alloc_err(ierror)
    CALL mem_add(v_size*n_wall)
    READ (f_unit) (wtangs(i),i=1,n_wall)
  END IF
  IF (wall_edges1 > 0) THEN
    ALLOCATE(wtangs1(wall_edges1),STAT=ierror)
    IF (ierror /= 0) CALL alloc_err(ierror)
    CALL mem_add(v_size*wall_edges1)
    READ (f_unit) (wtangs1(i),i=1,wall_edges1)
  END IF
  IF (n_bcedges > 0) THEN
    ALLOCATE(bcedges(n_bcedges),hbar(n_bcedges),qbar(n_bcedges),STAT=ierror)
    IF (ierror /= 0) CALL alloc_err(ierror)
    CALL mem_add(4*3*n_bcedges)
    READ (f_unit) (bcedges(i),i=1,n_bcedges)
    READ (f_unit) (hbar(i),i=1,n_bcedges)
    READ (f_unit) (qbar(i),i=1,n_bcedges)
  END IF

! From 'node_db'.
  ALLOCATE(cv_area(n_pts),cv_perim(n_pts),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(8*2*n_pts)
  READ (f_unit) (cv_area(i),cv_perim(i),i=1,n_pts)
  READ (f_unit) k
  ALLOCATE(n2t(n_pts,k),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(4*n_pts*k)
  READ (f_unit) ((n2t(i,j),j=1,k),i=1,n_pts)
  READ (f_unit) k
  ALLOCATE(n2t2(n_pts,k),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(4*n_pts*k)
  READ (f_unit) ((n2t2(i,j),j=1,k),i=1,n_pts)
  READ (f_unit) k
  ALLOCATE(n2n(n_pts,k),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(4*n_pts*k)
  READ (f_unit) ((n2n(i,j),j=1,k),i=1,n_pts)

! From 'elment_db'.
  READ (f_unit) k
  ALLOCATE(t2t(n_elems,k),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(4*n_elems*k)
  READ (f_unit) ((t2t(i,j),j=1,k),i=1,n_elems)
  READ (f_unit) k
  IF (k > 0) THEN
    ALLOCATE(t2tHD(n_elems,k),STAT=ierror)
    IF (ierror /= 0) CALL alloc_err(ierror)
    CALL mem_add(4*n_elems*k)
    READ (f_unit) ((t2tHD(i,j),j=1,k),i=1,n_elems)
  END IF
  ALLOCATE(t2t3(3,n_elems),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(4*n_elems*3)
  READ (f_unit) ((t2t3(j,i),j=1,3),i=1,n_elems)

  CLOSE (f_unit)

END SUBROUTINE mesh_in
