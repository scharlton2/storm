LOGICAL FUNCTION cell_centered(var,string)
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This function is used to check if a particular variable is located at the  !
!  cell centers or vertices of the elements when reading a Tecplot text       !
!  file.                                                                      !
!                                                                             !
!  INPUT:                                                                     !
!    var     a 2-character variable containing the location (order) of the    !
!            variable of interest in the Tecplot file.  For example, the      !
!            velocity component u is the 4th variable in the order used by    !
!            STORM, therefore var = '04';                                     !
!    string  the CELLCENTERED string from the Tecplot file header.            !
!                                                                             !
!  OUTPUT:                                                                    !
!    cell_centered  returns .TRUE. if 'string' contains 'var', .FALSE. if     !
!                   not.                                                      !
!                                                                             !
!  For example, the velocity component u is the 4th variable in the order     !
!  used by STORM, therefore var = '04'.                                       !
!                                                                             !
!  Francisco Simoes, February 2007                                            !
!  Last updated (mm-dd-yyyy): 02-27-2007 by F. Simões                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables.
  CHARACTER (LEN=*), INTENT(IN) :: var,string

  cell_centered = INDEX(string,var) /= 0

END FUNCTION cell_centered
