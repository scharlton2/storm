LOGICAL FUNCTION sortEdges(n,iarray)
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Get the right order for the boundary edges.  This FUNCTION is needed       !
!  because iRIC gives the inflow and outflow boundary nodes in a sorted       !
!  array, while SToRM needs it in an array with the order of the nodes that   !
!  forms the chain of edges in the right sequence.                            !
!                                                                             !
!  INPUT:                                                                     !
!    n           dimension of the array;                                      !
!    iarray      array with the indeces to the boundary points randomly       !
!                arranged (from iRIC).                                        !
!                                                                             !
!  OUTPUT:                                                                    !
!    iarray      array with the indeces of the boundary points arranged in    !
!                the right sequence that forms the contiguous grid edges;     !
!    sortEdges   error checking: .TRUE. if the sorting process is             !
!                successful, .FALSE. if it isn't.                             !
!                                                                             !
!  F. Simoes, 1 May 2013                                                      !
!  Last updated (mm-dd-yyyy): 05-07-2013 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  INTEGER, INTENT(IN) :: n
  INTEGER, DIMENSION(n), intent(INOUT) :: iarray

! Local variables.
  INTEGER :: i,j,k,m,mm
  INTEGER, DIMENSION(n) :: bndpts
  INTEGER, EXTERNAL :: is_edge

  bndpts = 0  ! Output array.
  mm = 0  ! Dimension of output array.
  m = n
  sortEdges = .TRUE.

! First sweep: add edges to the right of the first point in the input array.
  i = 1
  CALL add_array_ptoH(iarray(i),mm,bndpts)  ! Add point 1 to sorted array.
  CALL rmv_array_pto(i,m,iarray)  ! Remove point 1 from random array.
  DO
    k = 0
    DO j = 1,m  ! Find point in "iarray" that defines an edge with bndpts(i).
      k = is_edge(bndpts(i),iarray(j))
      IF (k > 0) EXIT
    END DO
    IF (k == 0) THEN  ! Fail: no more edges to be found.
      EXIT
    ELSE IF (k > 0) THEN
      CALL add_array_ptoT(iarray(j),mm,bndpts)  ! Add point to sorted array.
      CALL rmv_array_pto(j,m,iarray)  ! Remove point from random array.
      i = i + 1
    ELSE  ! Error condition: RETURN to calling program with error flag set.
      sortEdges = .FALSE.
      RETURN
    END IF
  END DO

! Second sweep: add edges to the left of the first point in the input array.
  DO
    k = 0
    DO j = 1,m  ! Find point in "iarray" that defines an edge with bndpts(1).
      k = is_edge(bndpts(1),iarray(j))
      IF (k > 0) EXIT
    END DO
    IF (k == 0) THEN  ! Fail: no more edges to be found.
      EXIT
    ELSE IF (k > 0) THEN
      CALL add_array_ptoH(iarray(j),mm,bndpts)  ! Add point to sorted array.
      CALL rmv_array_pto(j,m,iarray)  ! Remove point from random array.
    ELSE  ! Error condition: RETURN to calling program with error flag set.
      sortEdges = .FALSE.
      RETURN
    END IF
  END DO

! Last error check.
  IF (mm /= n) sortEdges = .FALSE.
  IF (m /= 0) sortEdges = .FALSE.

  iarray = bndpts

END FUNCTION sortEdges
