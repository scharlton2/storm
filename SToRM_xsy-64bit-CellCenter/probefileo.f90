SUBROUTINE probefileo(t,nn,hh,uu,vv)
  USE options
  USE io
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Output of probe data to an external file using a nice format.  The way     !
!  the arguments work facilitates the use of this subroutine with multiple    !
!  solvers: the h, u, and v are always defined at the vertices of the         !
!  triangles, therefore the arrays must be cell-vertex arrays.  By passing    !
!  these these quantities in the argument the code is not tied down to        !
!  particular solvers, COMMON BLOCKs, or variable names.                      !
!                                                                             !
!  INPUT:                                                                     !
!    t      time;                                                             !
!    hh     array containing the values of the water depth;                   !
!    uu     array containing the u velocity;                                  !
!    vv     array containing the v velocity;                                  !
!    nn     dimension of all the arrays.                                      !
!                                                                             !
!  Francisco Simoes, October 2008                                             !
!  Last updated (mm-dd-yyyy): 10-21-2008 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  INTEGER, INTENT(IN) :: nn
  REAL (KIND=mp), DIMENSION(nn), INTENT(IN) :: hh,uu,vv
  REAL (KIND=mp) :: t

! Local variables.
  INTEGER :: i

  WRITE (p_unit,probe_format) t,(hh(probepts(i)),uu(probepts(i)), &
                              vv(probepts(i)),i=1,nprobe)

END SUBROUTINE probefileo
