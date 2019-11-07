SUBROUTINE bed_slope_ETA(solver)
  USE parameters
  USE constants
  USE geometry
  USE dep_vars
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This subroutine computes the source term due to the bed slope using the    !
!  stage instead of the water depth.  It follows the approach in eq. (15) of  !
!  Bradford and Sanders (2002).                                               !
!                                                                             !
!  INPUT:                                                                     !
!    solver = 2 for FVT and other cell-centered solvers;                      !
!             3 for FVZ and other cell-centered solvers that use the free     !
!               surface elevation as dependent variable instead of water      !
!               depth.                                                        !
!                                                                             !
!  Francisco Simoes, October 2007                                             !
!  Last updated (mm-dd-yyyy): 11-26-2007 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  INTEGER, INTENT(IN) :: solver

! Local variables.
  REAL (KIND=mp) :: eta
  TYPE(vector) :: gradzb,gradzb2
  INTEGER :: i

  SELECT CASE (solver)

!-----------------------------------------------------------------------------!
!                                 FVT solver                                  !
!-----------------------------------------------------------------------------!

  CASE (2)
    DO i = 1,n_elems
      IF (drycell(i)) CYCLE

      eta = z(i) + h(i)

      ! The following code/correction is not needed anymore because the free
      ! surface gradient is enforced exactly at the boundaries between dry and
      ! wet cells.
      !IF (partdry(i)) THEN
      !  CALL adjustS0(.TRUE.,i,gradzb,gradzb2)
      !ELSE
      !  gradzb%x = delz(i)%x
      !  gradzb%y = delz(i)%y
      !  gradzb2%x = delz2(i)%x
      !  gradzb2%y = delz2(i)%y
      !END IF
      !phi(i,2) = phi(i,2) - g*(eta*gradzb%x - half*gradzb2%x)*grid(i)%area
      !phi(i,3) = phi(i,3) - g*(eta*gradzb%y - half*gradzb2%y)*grid(i)%area

      phi(i,2) = phi(i,2) - g*(eta*delz(i)%x - half*delz2(i)%x)*grid(i)%area
      phi(i,3) = phi(i,3) - g*(eta*delz(i)%y - half*delz2(i)%y)*grid(i)%area
    END DO

!-----------------------------------------------------------------------------!
!                                 FVZ solver                                  !
!-----------------------------------------------------------------------------!

  CASE (3)
    DO i = 1,n_elems
      IF (drycell(i)) CYCLE

      ! The following code/correction is not needed anymore because the free
      ! surface gradient is enforced exactly at the boundaries between dry and
      ! wet cells.
      !IF (partdry(i)) THEN
      !  CALL adjustS0(.TRUE.,i,gradzb,gradzb2)
      !ELSE
      !  gradzb%x = delz(i)%x
      !  gradzb%y = delz(i)%y
      !  gradzb2%x = delz2(i)%x
      !  gradzb2%y = delz2(i)%y
      !END IF
      !phi(i,2) = phi(i,2) - g*(zeta(i)*gradzb%x - half*gradzb2%x)*grid(i)%area
      !phi(i,3) = phi(i,3) - g*(zeta(i)*gradzb%y - half*gradzb2%y)*grid(i)%area

      phi(i,2) = phi(i,2) - g*(zeta(i)*delz(i)%x - half*delz2(i)%x)*grid(i)%area
      phi(i,3) = phi(i,3) - g*(zeta(i)*delz(i)%y - half*delz2(i)%y)*grid(i)%area
    END DO

  CASE DEFAULT
    PRINT *,"ERROR: invalid option in bed_slope_ETA."
    CALL byebye('Program STORM stopped.')

  END SELECT

END SUBROUTINE bed_slope_ETA
