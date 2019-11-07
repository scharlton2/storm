SUBROUTINE mesh_out(filename)
  USE geometry
  USE dep_vars
  USE vbc_arrays
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  All the mesh connectivity arrays are written in this suroutine.  The data  !
!  format and writing sequence must be identical to those in subroutine       !
!  mesh_in.                                                                   !
!                                                                             !
!  Francisco Simoes, November 2006                                            !
!  Last updated (mm-dd-yyyy): 03-17-2011 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  CHARACTER (LEN=*) :: filename

! Local variables.
  INTEGER :: f_unit,i,j,k

  LOGICAL, EXTERNAL :: get_iounit

  IF (.NOT. get_iounit(f_unit)) THEN
    PRINT *,'ERROR in get_iounit().'
    CALL byebye('Program SToRM stopped.')
  END IF

  OPEN (f_unit,FILE=filename,POSITION='REWIND',FORM='UNFORMATTED')

  WRITE (f_unit) (grid(i),i=1,n_elems)

! From 'find_edges'.
  WRITE (f_unit) n_edges
  IF (n_edges > 0) WRITE (f_unit) (edges(i),i=1,n_edges)

! From 'walls'.
  WRITE (f_unit) wall_edges,flow_edges,wall_edges1,flow_edges1,n_bpolygon, &
    n_wall
  IF (wall_edges > 0) WRITE (f_unit) (walledg(i),i=1,wall_edges)
  IF (flow_edges > 0) WRITE (f_unit) (flowedg(i),i=1,flow_edges)
  IF (wall_edges1 > 0) WRITE (f_unit) (walledg1(i),i=1,wall_edges1)
  IF (flow_edges1 > 0) WRITE (f_unit) (flowedg1(i),i=1,flow_edges1)
  IF (n_bpolygon > 0) WRITE (f_unit) (bpolygon(i),i=1,n_bpolygon)
  IF (n_wall > 0) WRITE (f_unit) (wall_pts(i),i=1,n_wall)

! From 'ibc_arrays'.
  WRITE (f_unit) n_inflowbdr,SIZE(pnormals,1),n_outflowbdr,SIZE(hnormals,1), &
    n_wall,wall_edges1,n_bcedges
  IF (n_inflowbdr > 0) THEN
    WRITE (f_unit) (n_qin(i),i=1,n_inflowbdr)
    k = SIZE(l_inflow,1)
    WRITE (f_unit) ((l_inflow(i,j),i=1,k),j=1,n_inflowbdr)
    k = SIZE(pnormals,1)
    WRITE (f_unit) ((pnormals(i,j),i=1,k),j=1,n_inflowbdr)
  END IF
  IF (n_outflowbdr > 0) THEN
    WRITE (f_unit) (n_hbc(i),i=1,n_outflowbdr)
    k = SIZE(hnormals,1)
    WRITE (f_unit) ((hnormals(i,j),i=1,k),j=1,n_outflowbdr)
  END IF
  IF (n_wall > 0) WRITE (f_unit) (wtangs(i),i=1,n_wall)
  IF (wall_edges1 > 0) WRITE (f_unit) (wtangs1(i),i=1,wall_edges1)
  IF (n_bcedges > 0) THEN
    WRITE (f_unit) (bcedges(i),i=1,n_bcedges)
    WRITE (f_unit) (hbar(i),i=1,n_bcedges)
    WRITE (f_unit) (qbar(i),i=1,n_bcedges)
  END IF

! From 'node_db'.
  WRITE (f_unit) (cv_area(i),cv_perim(i),i=1,n_pts)
  k = SIZE(n2t,2)
  WRITE (f_unit) k
  WRITE (f_unit) ((n2t(i,j),j=1,k),i=1,n_pts)
  k = SIZE(n2t2,2)
  WRITE (f_unit) k
  WRITE (f_unit) ((n2t2(i,j),j=1,k),i=1,n_pts)
  k = SIZE(n2n,2)
  WRITE (f_unit) k
  WRITE (f_unit) ((n2n(i,j),j=1,k),i=1,n_pts)

! From 'elemnt_db'.
  k = SIZE(t2t,2)
  WRITE (f_unit) k
  WRITE (f_unit) ((t2t(i,j),j=1,k),i=1,n_elems)
  k = 0
  IF (ALLOCATED(t2tHD)) k = SIZE(t2tHD,2)
  WRITE (f_unit) k
  IF (k > 0) WRITE (f_unit) ((t2tHD(i,j),j=1,k),i=1,n_elems)
  WRITE (f_unit) ((t2t3(j,i),j=1,3),i=1,n_elems)

  CLOSE (f_unit)

END SUBROUTINE mesh_out
