FUNCTION heffective(cell_no)
  USE parameters
  USE constants
  USE geometry
  USE dep_vars
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This function computes the effective water depth for a cell that is        !
!  partially wet.  It does not work if the cell is dry, therefore it must     !
!  only be used for wetted cells!                                             !
!                                                                             !
!  INPUT:                                                                     !
!    cell_no      element number of the cell where the effective water depth  !
!                 is sought.                                                  !
!  OUTPUT:                                                                    !
!    heffective   effective water depth.                                      !
!                                                                             !
!  F. Simoes, 2 February 2009                                                 !
!  Last updated (mm-dd-yyyy): 02-02-2009 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  REAL (KIND=mp) :: heffective
  INTEGER, INTENT(IN) :: cell_no

! Local variables.
  REAL (KIND=mp) :: depth,eta

  depth = h(cell_no)
  eta = depth + z(cell_no)
  IF (eta > zvtx(csortedz(cell_no)%p(3))) THEN
    ! If all vertices are wet, do nothing.
    heffective = depth
    RETURN
  END IF

  IF (eta > zvtx(csortedz(cell_no)%p(2))) THEN
    ! Two vertices wet.
    heffective = one_third*(2.0_mp*eta - zvtx(csortedz(cell_no)%p(1)) - &
                 zvtx(csortedz(cell_no)%p(2)))
  ELSE
    ! One vertex wet.
    heffective = one_third*(eta - zvtx(csortedz(cell_no)%p(1)))
  END IF

END FUNCTION heffective
