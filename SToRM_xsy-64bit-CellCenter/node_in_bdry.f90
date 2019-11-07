LOGICAL FUNCTION node_in_bdry(pto,e1,e2)
  USE geometry
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Function node_in_bdry() checks if a given node (pto) is located on the     !
!  bounding polygon of the computational grid.                                !
!                                                                             !
!    INPUT:                                                                   !
!      pto           index of the point to be checked.                        !
!    OUTPUT:                                                                  !
!      e1            index of one of the boundary edges containig point pto;  !
!      e2            index of the other boundary edge containing point pto;   !
!      node_in_bdry  .TRUE. if point pto is located on the bounding polygon   !
!                    of the computational grid; .FALSE. otherwise.            !
!                                                                             !
!  Francisco Simoes, November 2005                                            !
!  Last updated (mm-dd-yyyy): 11-11-2005 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dymmy arguments.
  INTEGER, INTENT(IN) :: pto
  INTEGER, INTENT(OUT) :: e1,e2

! Local variables.
  INTEGER :: counter,i
  LOGICAL, EXTERNAL :: pto_in_edge

  node_in_bdry = .FALSE.
  counter = 0

! Cycle over all the edges in the bounding polygon.
  DO i = 1,n_bpolygon
    IF (pto_in_edge(pto,edges(bpolygon(i)))) THEN
      e2 = bpolygon(i)
      counter = counter + 1
      IF (counter == 1) e1 = bpolygon(i)
    END IF
  END DO

  IF (counter == 0) RETURN  ! Node not in boundary.

  IF (counter == 2) THEN
    node_in_bdry = .TRUE.  ! Node in boundary.
    RETURN
  END IF

  ! This happens if counter = 1 or > 2, which are illegal values.
  PRINT *,"FATAL ERROR in FUNCTION node_in_bdry(); counter =",counter
  PRINT *,"Please report error to SToRM's development team."
  CALL byebye('Program SToRM stopped.')

END FUNCTION node_in_bdry
