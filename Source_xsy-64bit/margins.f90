SUBROUTINE margins
  USE constants
  USE geometry
  USE dep_vars
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This subroutine sets to zero the velocity vector in cells at river         !
!  margins.  River margins are defined as having at least one neighboring     !
!  cell fully wet and another fully dry.  Uses array drycell to check for     !
!  dry cells.                                                                 !
!                                                                             !
!  Francisco Simoes, July 2012                                                !
!  Last updated (mm-dd-yyyy): 07-13-2012 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables.
  INTEGER :: count1,count2,i,j,n

  DO i = 1,n_elems
    IF (drycell(i)) CYCLE
    IF (zeta(i) < zvtx(csortedz(i)%p(3))) THEN
      DO j = 1,3
        n = t2t3(j,i)
        count1 = 0  ! Counts the number of neighboring dry cells.
        count2 = 0  ! Counts the number of neighboring fully wet cells.
        IF (n < 1) THEN
          CYCLE
        ELSE IF (drycell(n)) THEN  ! Dry cell.
          count1 = count1 + 1
        ELSE IF (zeta(n) > zvtx(csortedz(n)%p(3))) THEN  ! Fully wet cell.
          count2 = count2 + 1
        END IF
      END DO
      IF ((count1 > 0) .AND. (count2 > 0)) THEN
        u(i) = zero
        v(i) = zero
      END IF
    END IF
  END DO
END SUBROUTINE margins

! This subroutine may be triggered unintentionally in wetting fronts.  Its
! impact to the behavior of wetting needs to be looked into with more detail.
