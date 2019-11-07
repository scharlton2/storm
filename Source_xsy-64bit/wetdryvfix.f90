SUBROUTINE wetdryvfix
  USE parameters
  USE constants
  USE geometry
  USE options
  USE dep_vars
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  The purpose of this subroutine is to set to zero the velocity at the       !
!  vertices of triangles for which at least one neighboring triangle is dry.  !
!  In other words, it sets to zero the velocity of vertices located at the    !
!  wet-dry boundaries.  This subroutine should only be used when the          !
!  no-slip wall boundary condition is chosen by the user.                     !
!                                                                             !
!  NOTE: it may set to zero the velocity at nodes at the inflow and outflow   !
!  boundaries...                                                              !
!                                                                             !
!  Francisco Simoes, January 2008                                             !
!  Last updated (mm-dd-yyyy): 01-22-2008 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables.
  INTEGER :: i,j,k
  LOGICAL :: dryt,wett

  DO i = 1,n_pts
    wett = .FALSE.
    dryt = .FALSE.
    ! First determine if the node is at the boundary between wet and dry
    ! triangles.
    DO j = 2,n2t(i,1) + 1
      k = n2t(i,j)
      IF (drycell(k)) THEN
        dryt = .TRUE.  ! Dry triangle.
      ELSE
        wett = .TRUE.  ! Wet triangle.
      END IF
    END DO

    IF (wett .AND. dryt) THEN  ! Yes, vertex is at a wet-dry boundary.
      uvtx(i) = zero
      vvtx(i) = zero
    END IF

  END DO

END SUBROUTINE wetdryvfix
