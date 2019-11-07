SUBROUTINE adjust_bed_elevation(zchange)
  USE parameters
  USE dep_vars
  USE geometry
  USE options
  USE constants
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This subroutine implements the bed elevation adjustment of Brufau and      !
!  Garcia-Navarro (2003), whereby the dry nodes of a partially wet triangle   !
!  are lowered to the level of the water surface elevation.                   !
!                                                                             !
!  OUTPUT:                                                                    !
!    zchange     .TRUE. if bed changes took place, .FALSE. if not.            !
!                                                                             !
!  Francisco Simoes, 3 February 2009                                          !
!  Last updated (mm-dd-yyyy): 03-03-2009 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  LOGICAL, INTENT(OUT) :: zchange
! Local avriables.
  INTEGER :: i,j

  zchange = .FALSE.
  DO i = 1,n_elems
    IF (h(i) < h_dry) CYCLE  ! Dry cell.
    IF (zeta(i) < zvtx(csortedz(i)%p(3))) THEN
      zchange = .TRUE.
      DO j = 1,3
        zvtx(grid(i)%vertex(j)) = MIN(zeta(i),zvtx(grid(i)%vertex(j)))
      END DO
      z(i) = (zvtx(grid(i)%vertex(1)) + zvtx(grid(i)%vertex(2)) + &
             zvtx(grid(i)%vertex(3)))*one_third
      h(i) = zeta(i) - z(i)
      DO j = 1,3
        hvtx(grid(i)%vertex(j)) = zetavtx(grid(i)%vertex(j)) - &
                                  zvtx(grid(i)%vertex(j))
      END DO
    END IF
  END DO

END SUBROUTINE adjust_bed_elevation
