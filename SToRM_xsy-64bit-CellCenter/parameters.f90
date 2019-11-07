!-----------------------------------------------------------------------------!
!                                                                             !
!  This file contains the definition of basic variable types and parameters   !
!  used in SToRM.                                                             !
!                                                                             !
!  F. Simoes, December 2003                                                   !
!  Last updated (mm-dd-yyyy): 02-02-2009                                      !
!                                                                             !
!-----------------------------------------------------------------------------!

MODULE parameters
  IMPLICIT NONE
  SAVE

! Machine precision:
  INTEGER, PARAMETER :: mp = KIND(1.0D0)  ! = KIND(1.0) for single precision.

! Points:
  TYPE :: point
    REAL (KIND=mp) :: x,y  ! Point coordinates.
  END TYPE point

! Edges (element sides):
  TYPE :: edge
    INTEGER :: p(2)  ! Points to beginning and end coordinate points.
    INTEGER :: e(2)  ! Points to elements on each side of the edge.
    REAL (KIND=mp) :: normal(2)  ! Components of the normal to the edge.
    REAL (KIND=mp) :: tang(2)    !     "      "   "  tangent "  "   "
    REAL (KIND=mp) :: length     ! Length of edge.
  END TYPE

! Computational grids:
  TYPE :: triangle
    INTEGER :: vertex(3)        ! Points to coordinates of vertices.
    INTEGER :: edge(3)          ! Points to edge information (can be < 0).
    INTEGER :: o_cell           ! Points to the cell where overdraft is taken.
    REAL (KIND=mp) :: area      ! Area of element.
    REAL (KIND=mp) :: perim     ! Perimeter of element.
    REAL (KIND = mp) :: xc,yc   ! Center coordinates.
  END TYPE

! Vectors:
  TYPE :: vector
    REAL (KIND=mp) :: x,y  ! Vector components.
  END TYPE vector

! Array with three INTEGER components.
  TYPE :: vector3
    INTEGER :: p(3)
  END TYPE

END MODULE parameters
