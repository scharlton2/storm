SUBROUTINE tecplot_write(id)
  USE parameters
  USE constants
  USE geometry
  USE dep_vars
  USE options
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Output of the nodal coordinates and element connectivity to an external    !
!  datafile in Tecplot format.  Also writes all the variables commonly used   !
!  by STORM, regardless of being located at the cell centers or at the cell   !
!  vertices.  This subroutine only writes the data itself; the Tecplot file   !
!  headers must be written before the call to this subroutine is made.        !
!                                                                             !
!  INPUT:                                                                     !
!    id      file unit where the data is written to, as obtained in the OPEN  !
!            Fortran statement.                                               !
!                                                                             !
!  Francisco Simoes, February 2007                                            !
!  Last updated (mm-dd-yyyy): 02-12-2009 by F. Simões                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables:
  INTEGER, INTENT(IN)  :: id

! Local variables:
  INTEGER :: i,j

  INTERFACE
    INTEGER FUNCTION npoints(varno,datav,n_pts,n_elems)
    INTEGER, INTENT(IN) :: varno,n_elems,n_pts
    LOGICAL, DIMENSION(:), INTENT(IN) :: datav
    END FUNCTION npoints
  END INTERFACE

! Write the nodal coordinates and solution information.
  DO i = 1,npoints(1,datav,n_pts,n_elems)  ! X-coordinate.
    WRITE (id,'(ES24.15)') nodes(i)%x
  END DO
  DO i = 1,npoints(2,datav,n_pts,n_elems)  ! Y-coordinate.
    WRITE (id,'(ES24.15)') nodes(i)%y
  END DO
! Bed elevation and other variables, in the proper order.
  DO i = 1,npoints(3,datav,n_pts,n_elems)
     !WRITE (id,'(ES24.15)') z(i)
     ! Write 'zvtx' instead of 'z' to be consistent with datav(3) = .TRUE.
     WRITE (id,'(ES24.15)') zvtx(i)
  END DO
  DO i = 1,npoints(4,datav,n_pts,n_elems)
    WRITE (id,'(ES24.15)') h(i)
  END DO
  DO i = 1,npoints(5,datav,n_pts,n_elems)
    WRITE (id,'(ES24.15)') u(i)
  END DO
  DO i = 1,npoints(6,datav,n_pts,n_elems)
    WRITE (id,'(ES24.15)') v(i)
  END DO
  DO i = 1,npoints(7,datav,n_pts,n_elems)
    WRITE (id,'(ES24.15)') cd(i)
  END DO
  DO j = 1,3
    DO i = 1,npoints(7+j,datav,n_pts,n_elems)
      WRITE (id,'(ES24.15)') phi(i,j)
    END DO
  END DO
  DO i = 1,npoints(11,datav,n_pts,n_elems)
    WRITE (id,'(ES24.15)') e_est(i)
  END DO

! Old code for DATAPACKING=POINT in the Tecplot ZONE field.
!  DO i = 1,n_pts
!    ! First, set velocity to zero on dry nodes (for pretty plots).
!    uvel = u(i)
!    vvel = v(i)
!    IF (h(i) < h_dry) THEN
!      uvel = zero
!      vvel = zero
!    END IF
!    WRITE (id,'(10(ES24.15))') nodes(i)%x,nodes(i)%y,z(i),h(i),uvel, &
!      vvel,cd(i),phi(i,1),phi(i,2),phi(i,3)
!  END DO

! Finally, write the element connectivity table.
  DO i = 1,n_elems
    WRITE (id,*) grid(i)%vertex(1),grid(i)%vertex(2),grid(i)%vertex(3)
  END DO

! For enhanced performance, some of the grid quantities are written here.  This
! is equivalent to similar code in 'tecplot_read' for reading 'triangle' data.

! First write the equivalent .edge data.
  WRITE (id,*)
  WRITE (id,*) n_edges
  DO i = 1,n_edges
    WRITE (id,*) i,edges(i)%p(1),edges(i)%p(2)
  END DO

! Now write the equivalent .v.edge data.
  WRITE (id,*)
  WRITE (id,*) n_edges
  DO i = 1,n_edges
    WRITE (id,*) i,edges(i)%e(1),edges(i)%e(2)
  END DO

END SUBROUTINE tecplot_write
