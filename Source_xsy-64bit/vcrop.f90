SUBROUTINE vcrop(vclip,vceiling)
  USE parameters
  USE dep_vars
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Points which have a flow velocity magnitude greater than vceiling get      !
!  their velocity magnitude clipped to the value specified in vceiling,       !
!  keeping the direction of the original velocity vector.                     !
!                                                                             !
!  INPUT:                                                                     !
!    vclip     = .TRUE. if clipping is activated, .FALSE. otherwise;          !
!    vceiling  maximum magnitude of velocity, above which value the velocity  !
!              is reduced to vceiling.                                        !
!                                                                             !
!  Francisco Simoes, April 2011                                               !
!  Last updated (mm-dd-yyyy): 04-14-2011 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables.
  LOGICAL, INTENT(IN) :: vclip
  REAL(KIND=mp), INTENT(IN) :: vceiling

  IF (.NOT. vclip) RETURN

  WHERE (u_mag > vceiling)
    u = vceiling*u/u_mag
    v = vceiling*v/u_mag
    u_mag = vceiling
  END WHERE

END SUBROUTINE vcrop
