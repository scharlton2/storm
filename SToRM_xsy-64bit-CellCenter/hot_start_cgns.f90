SUBROUTINE hot_start_cgns(cgnsfile,n,h,u,v)
  USE parameters
  IMPLICIT NONE
  INCLUDE "Header Files\cgnslib_f.h"

!-----------------------------------------------------------------------------!
!                                                                             !
!  Reads an existing solution from an iRIC 2.0 CGNS file.                     !
!                                                                             !
!  INPUT:                                                                     !
!    cgnsfile   filename of CGNS file (iRIC 2.0);                             !
!    n          dimension of arrays (number of triangle vertices in the       !
!               computational grid.                                           !
!                                                                             !
!  OUTPUT:                                                                    !
!    h          water depth;                                                  !
!    u          x-component of the flow velocity;                             !
!    v          y-component of the flow velocity.                             !
!                                                                             !
!  Francisco Simoes, 13 July 2012                                             !
!  Last updated (mm-dd-yyyy): 11-03-2016 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables.
  INTEGER, INTENT(IN) :: n
  REAL (KIND=mp), INTENT(OUT), DIMENSION(n) :: h,u,v
  CHARACTER(*), INTENT(IN) :: cgnsfile

! Local variables.
  INTEGER :: fid,ierror,iloc,index_base,index_zone,isize(3,3),npts,nsols
  CHARACTER (LEN=32) :: solname,zonename

! Open CGNS file for read.
  CALL cg_open_f(cgnsfile,CG_MODE_READ,fid,ierror)
  IF (ierror /= 0) THEN
    CALL cg_error_print_f()
    CALL byebye('SToRM ERROR: cannot open CGNS solution file.')
  END IF

  index_base = 1  ! These two values are hardcoded for iRIC 2.0.
  index_zone = 1
  CALL cg_zone_read_f(fid,index_base,index_zone,zonename,isize,ierror)
  npts = isize(1,1)
  IF (npts /= n) THEN  ! Error checking.
    PRINT *,'Error reading starting solution:'
    PRINT *,'number of nodes does not match current grid size.'
    CALL byebye('Program SToRM stopped.')
  END IF

  CALL cg_nsols_f(fid,index_base,index_zone,nsols,ierror)
  CALL cg_sol_info_f(fid,index_base,index_zone,nsols,solname,iloc,ierror)
  CALL cg_field_read_f(fid,index_base,index_zone,nsols,'VelocityX', &
                       RealDouble,1,n,u,ierror)
  IF (ierror /= 0) THEN
    CALL cg_error_print_f()
    CALL byebye('SToRM CGNS ERROR: cannot read initial U-velocity.')
  END IF
  CALL cg_field_read_f(fid,index_base,index_zone,nsols,'VelocityY', &
                       RealDouble,1,n,v,ierror)
  IF (ierror /= 0) THEN
    CALL cg_error_print_f()
    CALL byebye('SToRM CGNS ERROR: cannot read initial V-velocity.')
  END IF
  CALL cg_field_read_f(fid,index_base,index_zone,nsols,'Depth', &
                       RealDouble,1,n,h,ierror)
  IF (ierror /= 0) THEN
    CALL cg_error_print_f()
    CALL byebye('SToRM CGNS ERROR: cannot read initial water depth.')
  END IF

  CALL cg_close_f(fid,ierror)

END SUBROUTINE hot_start_cgns
