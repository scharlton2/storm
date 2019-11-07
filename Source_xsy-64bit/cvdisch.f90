SUBROUTINE cvdisch
  USE dep_vars
  USE constants
  USE options
  USE geometry
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This SUBROUTINE computes the source/sink terms in the continuity equation  !
!  due to the presence of culverts.                                           !
!                                                                             !
!  F. Simoes, February 2013                                                   !
!  Last updated (mm-dd-yyyy): 05-02-2013 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables:
  INTEGER :: i,j
  REAL (KIND=mp) :: head,t

  cvsrc = zero
  DO i = 1,nculvert
    IF (drycell(cvtrigin(i))) CYCLE  ! Cell is dry.
    ! Compute headwater elevation at culvert inlet = zeta + grad_zeta dot rcv.
    head = zeta(cvtrigin(i)) + (delh(cvtrigin(i))%x + &
      delz(cvtrigin(i))%x)*cvrin(1,i) + (delh(cvtrigin(i))%y + &
      delz(cvtrigin(i))%y)*cvrin(2,i)
    head = head - cvptoin(3,i)
    IF (head < h_dry) CYCLE  ! Inlet of culvert is dry.

    ! Use culvert performance curve to find the discharge through it using
    ! linear interpolation.
    IF (head < cvhead(1,i)) THEN
      cvsrc(i) = cvq(1,i)
      CYCLE
    ELSE IF (head >= cvhead(ncvtable(i),i)) THEN
      cvsrc(i) = cvq(ncvtable(i),i)
      CYCLE
    ELSE
      DO j = 1,ncvtable(i)-1
        IF (head < cvhead(j+1,i)) THEN
          t = (head - cvhead(j,i))/(cvhead(j+1,i) - cvhead(j,i))
          cvsrc(i) = cvq(j+1,i)*t + cvq(j,i)*(one - t)
          EXIT
        END IF
      END DO
    END IF
  END DO

END SUBROUTINE cvdisch
