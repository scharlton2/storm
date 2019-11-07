MODULE parameters
  IMPLICIT NONE
  SAVE

! Machine precision:
  INTEGER, PARAMETER :: mp = KIND(1.0D0)  ! = KIND(1.0) for single precision.

END MODULE parameters


SUBROUTINE ghost_cell(xc,yc,x1,y1,x2,y2,xg,yg)
  USE parameters
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Computes the center of a ghost cell, i.e., of a cell that is the           !
!  symmetric (reflection) of the current cell.  The line of symmetry is the   !
!  prescribed edge.                                                           !
!                                                                             !
!  INPUT:                                                                     !
!    xc,yc     coordinates of the center of the cell;                         !
!    x1,y1     coordinates of one of the end points of the edge that is       !
!              going to be used as line of symmetry;                          !
!    x2,y2     coordinates of the other end point of the same edge.           !
!                                                                             !
!  OUTPUT:                                                                    !
!    xg,yg     center of the ghost cell, placed symmetrically to (xc,yc).     !
!                                                                             !
!  NOTE: the order of point 1 and point 2 does not matter.                    !
!                                                                             !
!  Francisco Simoes, February 2009                                            !
!  Last updated (mm-dd-yyyy): 02-03-2009 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  REAL(KIND=mp), INTENT(IN) :: xc,yc,x1,y1,x2,y2
  REAL(KIND=mp), INTENT(OUT) :: xg,yg

! Local variables.
  REAL(KIND=mp) :: ax,ay,bx,by,cx,cy,denom

  ax = xc - x1;  ay = yc - y1
  bx = x2 - x1;  by = y2 - y1
  denom = bx*bx + by*by
  cx = (bx*(ax*bx + ay*by) + by*(ay*bx - ax*by))/denom
  cy = (bx*(ax*by - ay*bx) + by*(ax*bx + ay*by))/denom

  xg = cx + x1
  yg = cy + y1

END SUBROUTINE ghost_cell


PROGRAM boo
  USE parameters
  IMPLICIT NONE

  REAL(KIND=mp) :: xc,yc,xg,yg,x1,y1,x2,y2,x3,y3

  x1 = 3.0;  y1 = 3.0
  x2 = 5.0;  y2 = 5.0
  xc = 8.0_mp/3.0_mp
  yc = 4.0
  CALL ghost_cell(xc,yc,x1,y1,x2,y2,xg,yg)
  PRINT *,xg,yg

  x2 = 3.0;  y2 = 3.0
  x1 = 5.0;  y1 = 5.0
  CALL ghost_cell(xc,yc,x1,y1,x2,y2,xg,yg)
  PRINT *,xg,yg

  x1 = 5.0;  y1 = 2.0
  x2 = 5.0;  y2 = 5.0
  x3 = 1.0;  y3 = 1.0
  xc = (x1 + x2 + x3)/3.0_mp
  yc = (y1 + y2 + y3)/3.0_mp
  PRINT *,xc,yc
  x3 = 9.0; y3 = 1.0
  xg = (x1 + x2 + x3)/3.0_mp
  yg = (y1 + y2 + y3)/3.0_mp
  PRINT *,xg,yg
  CALL ghost_cell(xc,yc,x1,y1,x2,y2,xg,yg)
  PRINT *,xg,yg
  CALL ghost_cell(xc,yc,x2,y2,x1,y1,xg,yg)
  PRINT *,xg,yg

END PROGRAM boo
