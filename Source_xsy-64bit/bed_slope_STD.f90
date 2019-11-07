SUBROUTINE bed_slope_STD
  USE parameters
  USE constants
  USE geometry
  USE dep_vars
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Compute the source terms due to the bed slope using the standard least     !
!  squares computation of the bed gradient.                                   !
!                                                                             !
!  Francisco Simoes, October 2007                                             !
!  Last updated (mm-dd-yyyy): 02-18-2009 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables.
  INTEGER :: i,j,k
  TYPE(vector) :: dummy,gradzb
  LOGICAL :: set_grad

  DO i = 1,n_elems
    IF (drycell(i)) CYCLE
    !IF (bdrycell(i)) CYCLE

    !--------------------------------------------------------------------------
    ! The following code/correction is not needed anymore because the free
    ! surface gradient is enforced exactly at the boundaries between dry and
    ! wet cells.
    !IF (partdry(i)) THEN
    !  CALL adjustS0(.FALSE.,i,gradzb,dummy)
    !ELSE
    !  gradzb%x = delz(i)%x
    !  gradzb%y = delz(i)%y
    !END IF

    !--------------------------------------------------------------------------
    ! This code sets the bed gradient temporarily to zero if there are
    ! neighboring dry cells.
    !set_grad = .FALSE.
    !DO j = 1,3
    !  k = t2t3(j,i)
    !  IF (k < 1) CYCLE  ! Edge triangle.
    !  IF (drycell(k)) set_grad = .TRUE.
    !END DO
    !IF (set_grad) THEN
    !  gradzb%x = zero
    !  gradzb%y = zero
    !ELSE
    !  gradzb%x = delz(i)%x
    !  gradzb%y = delz(i)%y
    !END IF

    !--------------------------------------------------------------------------
    ! Yet another variation: this method uses the gradient of the water depth
    ! to determine if any of the vertices is dry and adjusts the bed slope
    ! accordingly.
    CALL adjustZb(.FALSE.,i,gradzb,dummy)

    phi(i,2) = phi(i,2) - g*h(i)*gradzb%x*grid(i)%area
    phi(i,3) = phi(i,3) - g*h(i)*gradzb%y*grid(i)%area

    !phi(i,2) = phi(i,2) - g*h(i)*delz(i)%x*grid(i)%area
    !phi(i,3) = phi(i,3) - g*h(i)*delz(i)%y*grid(i)%area

    !WRITE (*,'(I8,4(3X,ES12.5))') i,source(i,1),g*h(i)*delz(i)%x*grid(i)%area, &
    !  source(i,2),g*h(i)*delz(i)%y*grid(i)%area

  END DO
  !CALL byebye('bed_slope_STD')

END SUBROUTINE bed_slope_STD
