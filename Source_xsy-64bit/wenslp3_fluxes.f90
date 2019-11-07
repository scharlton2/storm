SUBROUTINE wenslp3_fluxes(hdry)
  USE parameters
  USE constants
  USE geometry
  USE dep_vars
  USE vbc_arrays
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Computes the contributions to the residual due to the inviscid wall        !
!  convective fluxes for a no-slip boundary.  It enforces the proper          !
!  inviscid fluxes exactly without using Roe's flux, using eq. (29) of Kuiry  !
!  et al. (2008) to compute the water depth at the wall edge.                 !
!                                                                             !
!  INPUT:                                                                     !
!    hdry   water depth threshold below which the cell is considered dry.     !
!                                                                             !
!  Francisco Simoes, April 2008                                               !
!  Last updated (mm-dd-yyyy): 04-15-2008 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  REAL (KIND=mp), INTENT(IN) :: hdry

! Local variables.
  INTEGER :: i,k,l,m
  REAL (KIND=mp) :: eps,h1,h2,hedge
  TYPE(edge) :: ework  ! Working edge.
  INTEGER, EXTERNAL :: edge_in_element,local_edge
  REAL (KIND=mp), EXTERNAL :: invert_edge_dirs

  DO l = 1,wall_edges1
    k = walledg1(l)
    i = edge_in_element(edges(k))
    CALL edge_copy(edges(k),ework)
    eps = one

    IF (drycell(i)) CYCLE

    ! Compute the water depth at the element's edge.
    h1 = hvtx(ework%p(1))
    h2 = hvtx(ework%p(2))
    hedge = SQRT((h1*h1 + h1*h2 + h2*h2)*one_third)

    IF (hedge < zero) THEN  ! Flag cell for bed slope treatment.
      hedge = zero
      partdry(i) = .TRUE.
    END IF
    IF (hedge < hdry) CYCLE

    ! Make sure that the working normal points out from element i, i.e., it
    ! points out of the domain.
    m = local_edge(grid(i),k)
    IF (grid(i)%edge(m) < 0) eps = invert_edge_dirs(ework)

    ! Add the edge flux to the residual's array for element i.  Note that the
    ! flux of mass through the wall must be zero, therefore f(1) = 0.
    phi(i,2) = phi(i,2) - half*g*hedge*hedge*ework%normal(1)*ework%length
    phi(i,3) = phi(i,3) - half*g*hedge*hedge*ework%normal(2)*ework%length

    !IF (k == 4131) THEN
    !  WRITE (*,'(2I6,3ES15.6)') k,i,half*g*hedge*hedge*ework%normal(1), &
    !    half*g*hedge*hedge*ework%normal(2)
    !END IF
  END DO

END SUBROUTINE wenslp3_fluxes
