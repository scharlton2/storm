INTEGER FUNCTION npoints(varno,datav,n_pts,n_elems)
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This function is designed for in the I/O process via Tecplot files, to     !
!  determine if a particular variable is cell-centered or vertex-centered.    !
!  It uses the information in array datav() to determine that.                !
!                                                                             !
!  INPUT:                                                                     !
!    varno      the variable order, as it is read and written in the Tecplot  !
!               file. Currently, it goes x (varno = 1), y (2), z (3), h (4),  !
!               u (5), v(6), cd (7), etc.;                                    !
!    datav      the array that indicates if the variable is located at the    !
!               vertices (datav(varno) = .TRUE.) or at the element centers    !
!               (datav(varno) = .FALSE.);                                     !
!    n_pts      number of points (triangle vertices) in the computational     !
!               grid;                                                         !
!    n_elems    number of elements in the computational grid.                 !
!                                                                             !
!  OUTPUT:                                                                    !
!    npoints  = n_pts if the variable is vertex-centered, = n_elems if the    !
!               variable is cell-centered.                                    !
!                                                                             !
!  NOTE: this function must be used in an INTERFACE construct in the calling  !
!  program:                                                                   !
!                                                                             !
!    INTERFACE                                                                !
!      INTEGER FUNCTION npoints(varno,datav,n_pts,n_elems)                    !
!      INTEGER, INTENT(IN) :: varno,n_elems,n_pts                             !
!      LOGICAL, DIMENSION(:), INTENT(IN) :: datav                             !
!      END FUNCTION npoints                                                   !
!    END INTERFACE                                                            !
!                                                                             !
!  Francisco Simoes, February 2007                                            !
!  Last updated (mm-dd-yyyy): 08-20-2007 by F. Simões                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables.
  INTEGER, INTENT(IN) :: varno,n_elems,n_pts
  LOGICAL, DIMENSION(:), INTENT(IN) :: datav

  IF (varno > SIZE(datav) .OR. varno < 1) THEN
    PRINT *,"ERROR: bad call to function npoints."
    CALL byebye('SToRM stopped.')
  END IF

  npoints = n_elems
  IF (datav(varno)) npoints = n_pts

END FUNCTION npoints
