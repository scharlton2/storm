SUBROUTINE bed_slope
  USE parameters
  USE geometry
  USE dep_vars
  USE constants
  USE options
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Compute the source terms due to the bed slope.                             !
!                                                                             !
!  Francisco Simoes, November 2005                                            !
!  Last updated (mm-dd-yyyy): 03-22-2007 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables:
  INTEGER :: i
  REAL(KIND=mp) :: h_cv,s0x,s0y,x1,x2,x3,y1,y2,y3,znew(3)
  TYPE(triangle) :: cv

  DO i = 1,n_elems

    ! Skip dry or stagnant cells.
    IF (u_avg(i,1) < h_dry) CYCLE
    IF ((u_avg(i,2)*u_avg(i,2) + u_avg(i,3)*u_avg(i,3)) < u_stag*u_stag) CYCLE

    cv = grid(i)  ! Control volume.

    ! Compute bed slope within the element.
    x1 = nodes(cv%vertex(1))%x; y1 = nodes(cv%vertex(1))%y
    x2 = nodes(cv%vertex(2))%x; y2 = nodes(cv%vertex(2))%y
    x3 = nodes(cv%vertex(3))%x; y3 = nodes(cv%vertex(3))%y
    ! Apply the algorithm of Brufau and Garcia-Navarro (2003) to partially wet
    ! cells.
    CALL drynwet(cv,znew)
    ! The bed slope is constant within the control volume, therefore it can be
    ! kept outside of the integral over the control volume.
    s0x = -half*(znew(1)*(y2 - y3) + znew(2)*(y3 - y1) + znew(3)*(y1 - y2))
    s0y = -half*(znew(1)*(x3 - x2) + znew(2)*(x1 - x3) + znew(3)*(x2 - x1))
    ! Integrate h within the control volume.
    h_cv = u_avg(i,1)
    ! Math note: in the operations above the term in cv%area droped because
    ! it appears in the numerator and in the denominator.

    ! Source term due to bed slope.  Note that we don't multiply by the area
    ! of the control volume because it is not necessary.
    source(i,1) = source(i,1) + zero        ! Continuity.
    source(i,2) = source(i,2) + g*h_cv*s0x  ! x-momentum.
    source(i,3) = source(i,3) + g*h_cv*s0y  ! y-momentum.
    !PRINT *,i,s0x,s0y
  END DO

END SUBROUTINE bed_slope
