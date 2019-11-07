SUBROUTINE divide(i,k,grid,nelems,npts)
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Given a triangle and a point on an edge, this subroutine divides it into   !
!  two triangles using that point and keeping the circulation of the          !
!  vertices.                                                                  !
!                                                                             !
!  INPUT:       i  element in array grid to be subdivided;                    !
!               k  face where point npts is located;                          !
!            grid  array with the triangle information;                       !
!          nelems  number of elements in array grid;                          !
!            npts  point to be used as the new vertex in the subdivision.     !
!                                                                             !
!  OUTPUT:   grid  this array is now modified (the original triangle has new  !
!                  connectivity and a new triangle is appended to the         !
!                  array).  Make sure that the array dimensions are large     !
!                  enough to accomodate the new triangle;                     !
!          nelems  number of elements in array grid (= nelems + 1).           !
!                                                                             !
!  Note that this subroutine has a Fortran 90-style assumed shape array       !
!  (containing dimension(:)) that requires the subroutine to be place in a    !
!  MODULE or to have an explicit interface wherever it is used.               !
!                                                                             !
!  Francisco Simoes, Mar. 2004                                                !
!                                                                             !
!-----------------------------------------------------------------------------!
!  Date of last change: 02-27-2007, F. Simoes

! Machine precision.
  INTEGER, PARAMETER :: mp = KIND(1.0D0) ! = KIND(1.0) for single precision.

! Computational grids (as defined in the calling SUBROUTINE corner()).
  TYPE :: triangle
    INTEGER :: vertex(3)        ! Points to coordinates of vertices.
    INTEGER :: edge(3)          ! Points to edge information (can be < 0).
    INTEGER :: s                ! Sum of the frequency of all the edges.
!   REAL(KIND=mp) :: area       ! Area of element.
  END TYPE

! Dummy arguments.
  INTEGER, INTENT(IN) :: i,k,npts
  INTEGER, INTENT(INOUT) :: nelems
  TYPE(triangle), DIMENSION(:), INTENT(INOUT) :: grid

! Local variables.
  INTEGER :: m
  TYPE(triangle) :: oldt

  oldt = grid(i)
  nelems = nelems + 1
  DO m = 1,3
    IF (k == ABS(oldt%edge(m))) EXIT
  END DO

  SELECT CASE (m)
  CASE (1)
    grid(i)%vertex(1) = oldt%vertex(1)
    grid(i)%vertex(2) = npts
    grid(i)%vertex(3) = oldt%vertex(3)
    grid(nelems)%vertex(1) = npts
    grid(nelems)%vertex(2) = oldt%vertex(2)
    grid(nelems)%vertex(3) = oldt%vertex(3)
  CASE (2)
    grid(i)%vertex(1) = oldt%vertex(1)
    grid(i)%vertex(2) = oldt%vertex(2)
    grid(i)%vertex(3) = npts
    grid(nelems)%vertex(1) = npts
    grid(nelems)%vertex(2) = oldt%vertex(3)
    grid(nelems)%vertex(3) = oldt%vertex(1)
  CASE (3)
    grid(i)%vertex(1) = oldt%vertex(1)
    grid(i)%vertex(2) = oldt%vertex(2)
    grid(i)%vertex(3) = npts
    grid(nelems)%vertex(1) = npts
    grid(nelems)%vertex(2) = oldt%vertex(2)
    grid(nelems)%vertex(3) = oldt%vertex(3)
  CASE DEFAULT
    WRITE (*,'(" ERROR: invalid edge in subdivision, triangle",I7)') i
    WRITE (*,'(" Program CORNER stopped.")')
    STOP
  END SELECT

END SUBROUTINE divide
