!-----------------------------------------------------------------------------!
!                                                                             !
!  This file contains the definition of basic variable types and parameters   !
!  used in SToRM.                                                             !
!                                                                             !
!  F. Simoes, December 2003                                                   !
!  Last updated (mm-dd-yyyy): 11-15-2004                                      !
!                                                                             !
!-----------------------------------------------------------------------------!

MODULE parameters
  IMPLICIT NONE
  SAVE

! Machine precision:
  INTEGER, PARAMETER :: mp = KIND(1.0D0)  ! = KIND(1.0) for single precision.

! Points:
  TYPE :: point
    REAL(KIND=mp) :: x,y  ! Point coordinates.
  END TYPE point

! Edges (element sides):
  TYPE :: edge
    INTEGER :: p(2)  ! Points to beginning and end coordinate points.
    REAL(KIND=mp) :: normal(2)  ! Components of the normal to the edge.
    REAL(KIND=mp) :: length     ! Length of edge.
  END TYPE

! Computational grids:
  TYPE :: triangle
    INTEGER :: vertex(3)        ! Points to coordinates of vertices.
    INTEGER :: edge(3)          ! Points to edge information (can be < 0).
    REAL(KIND=mp) :: area       ! Area of element.
  END TYPE

END MODULE parameters
