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
  REAL (KIND=mp), INTENT(IN) :: xc,yc,x1,y1,x2,y2
  REAL (KIND=mp), INTENT(OUT) :: xg,yg

! Local variables.
  REAL (KIND=mp) :: ax,ay,bx,by,cx,cy,denom

  ax = xc - x1;  ay = yc - y1
  bx = x2 - x1;  by = y2 - y1
  denom = bx*bx + by*by
  cx = (bx*(ax*bx + ay*by) + by*(ay*bx - ax*by))/denom
  cy = (bx*(ax*by - ay*bx) + by*(ax*bx + ay*by))/denom

  xg = cx + x1
  yg = cy + y1

END SUBROUTINE ghost_cell
