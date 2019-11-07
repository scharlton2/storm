SUBROUTINE tecplot_read(id,varloc)
  USE geometry
  USE dep_vars
  USE memory
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Read in nodal coordinates and element connectivity from an external data   !
!  file in Tecplot format.  Read also all the variables commonly used by      !
!  STORM, regardless of being located at the cell centers or at the cell      !
!  vertices.  This subroutine only reads the data itself; the Tecplot file    !
!  headers must be read before the call to this subroutine is made.  Note     !
!  also that the array datav() is set here so that it can be used also in     !
!  the output routines.                                                       !
!                                                                             !
!  INPUT:                                                                     !
!    id      file unit where the data is read from, as obtained in the OPEN   !
!            Fortran statement;                                               !
!    varloc  the CELLCENTERED string from the Tecplot file header.            !
!                                                                             !
!  Francisco Simoes, February 2007                                            !
!  Last updated (mm-dd-yyyy): 04-05-2008 by F. Simões                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables.
  INTEGER, INTENT(IN) :: id
  CHARACTER (LEN=*), INTENT(IN) :: varloc

! Local variables.
  INTEGER :: i,ierror,j,lineno
  REAL (KIND=mp) :: var
  !CHARACTER (LEN=120) :: buffer
  CHARACTER (LEN=40) :: buffer  ! Size set for readline().
  LOGICAL :: line
  LOGICAL, EXTERNAL :: cell_centered,readline

  INTERFACE
    INTEGER FUNCTION npoints(varno,datav,n_pts,n_elems)
    INTEGER, INTENT(IN) :: varno,n_elems,n_pts
    LOGICAL, DIMENSION(:), INTENT(IN) :: datav
    END FUNCTION npoints
  END INTERFACE

! Read point coordinates, bed elevation (x,y,z), and friction coefficients.  An
! initial solution must be included here, too, but it may be set to zero if
! necessary.
  datav = .TRUE.
  DO i = 1,n_pts  ! Variable #1 in the Tecplot order,
    READ (id,*) nodes(i)%x
  END DO
  DO i = 1,n_pts  ! variable #2 in the Tecplot order,
    READ (id,*) nodes(i)%y
  END DO
! Variable #3 in the Tecplot order, etc.
  IF (cell_centered('03',varloc)) datav(3) = .FALSE.
  DO i = 1,npoints(3,datav,n_pts,n_elems)
    READ (id,*) z_tcp(i)
  END DO
  IF (cell_centered('04',varloc)) datav(4) = .FALSE.
  DO i = 1,npoints(4,datav,n_pts,n_elems)
    READ (id,*) h_tcp(i)
  END DO
  IF (cell_centered('05',varloc)) datav(5) = .FALSE.
  DO i = 1,npoints(5,datav,n_pts,n_elems)
    READ (id,*) u_tcp(i)
  END DO
  IF (cell_centered('06',varloc)) datav(6) = .FALSE.
  DO i = 1,npoints(6,datav,n_pts,n_elems)
    READ (id,*) v_tcp(i)
  END DO
  IF (cell_centered('07',varloc)) datav(7) = .FALSE.
  DO i = 1,npoints(7,datav,n_pts,n_elems)
    READ (id,*) cd_tcp(i)
  END DO

! Read the other parameters that are written for every STORM output.
  IF (cell_centered('08',varloc)) datav(8) = .FALSE.
  DO i = 1,npoints(8,datav,n_pts,n_elems)
    READ (id,*) var
  END DO
  IF (cell_centered('09',varloc)) datav(9) = .FALSE.
  DO i = 1,npoints(9,datav,n_pts,n_elems)
    READ (id,*) var
  END DO
  IF (cell_centered('10',varloc)) datav(10) = .FALSE.
  DO i = 1,npoints(10,datav,n_pts,n_elems)
    READ (id,*) var
  END DO
  IF (cell_centered('11',varloc)) datav(11) = .FALSE.
  DO i = 1,npoints(11,datav,n_pts,n_elems)
    READ (id,*) var
  END DO

! Read element connectivity, counterclockwise.
  DO i = 1,n_elems
    READ (id,*) grid(i)%vertex(1),grid(i)%vertex(2),grid(i)%vertex(3)
  END DO

! For enhanced performance, some of the grid quantities may be read here.  This
! area of the file will be ignored by Tecplot when reading the file.  The
! format used is that taken directly from the .edge and .v.edge files generated
! by the mesh generator 'triangle' and may be added to the SToRM grid input
! file by a simple copy-and-paste operation, one after the other, respectively.
! Blank lines and comment lines starting with # (such as the ones usually put
! at the end of those files by program 'triangle') may be used, but only
! between the places where the data is pasted onto.

! I don't know if the data from 'triangle' work with the RDS solver because its
! code has the requirement that the sequence in which edges are set in
! grid()%edge() must follow this rule: each node i must be followed by edge i
! so that the edge opposite to node i can easily be found from
! MAX(1,MOD(i+1,4)).  This needs to be checked at a later time...

! First read the .edge data.
  lineno = 0
  DO
    line = readline(id,buffer,lineno)
    IF (.NOT. line) RETURN  ! EOF reached.
    READ (buffer,*) n_edges
    ! Allocate storage for global array edges.
    ALLOCATE(edges(n_edges),STAT=ierror)
    IF (ierror /= 0) CALL alloc_err(ierror)
    CALL mem_add(e_size*n_edges)
    DO i = 1,n_edges
      READ (id,*) j,edges(i)%p(1),edges(i)%p(2)
    END DO
    EXIT
  END DO

! Now read the .v.edge data.
  lineno = 0
  DO
    line = readline(id,buffer,lineno)
    IF (.NOT. line) CALL byebye("Error reading grid file: edge data missing.")
    DO i = 1,n_edges
      READ (id,*) j,edges(i)%e(1),edges(i)%e(2)
      IF (edges(i)%e(1) < 0) THEN  ! Swap nodes.
        edges(i)%e(1) = edges(i)%e(2)
        edges(i)%e(2) = -1
      END IF
    END DO
    EXIT
  END DO

END SUBROUTINE tecplot_read
