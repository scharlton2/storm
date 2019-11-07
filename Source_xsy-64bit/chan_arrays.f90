SUBROUTINE chan_arrays
  USE parameters
  USE constants
  USE geometry
  USE dep_vars
  USE memory
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Calculation of preliminary arrays needed later for implementing the        !
!  one-dimensional channel computations.                                      !
!                                                                             !
!  F. Simoes, September 2011                                                  !
!  Last updated (mm-dd-yyyy): 10-02-2011 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables.
  INTEGER :: d1,d2,e,i,i0,ierror,j,t1,t2
  INTEGER, EXTERNAL :: is_edge
  REAL (KIND=mp) :: length,rlength

  d1 = SIZE(channel,1)  ! Max number of points per boundary.
  d2 = SIZE(channel,2)  ! Number of outflow boundaries.
  IF (d2 /= n_chanbdr) THEN
    PRINT *,"ERROR in 1D channel arrays:,"
    PRINT *,"Dimension does not match count."
    CALL byebye('Program SToRM stopped.')
  END IF

  ALLOCATE(chanedges(0:d1,d2),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(4*(d1 + 1)*d2)

! Find the strings of outflow boundary edges.
  length = 0  ! Used to parametrize the channels.
  DO i0 = 1,n_chanbdr
    chanedges(0,i0) = n_channel(i0) - 1  ! Number of edges for current channel.
    DO i = 1,n_channel(i0)-1

      ! Find the edge defined by the two nodes.
      e = is_edge(channel(i,i0),channel(i+1,i0))
      IF (e == 0) THEN
        PRINT *,"ERROR in 1D channel node assignement,"
        PRINT *,"nodes",channel(i,i0),channel(i+1,i0)
        PRINT *,"These nodes are not connected by an edge."
        PRINT *,"Please revise your channel link data."
        CALL byebye('Program SToRM stopped.')
      ELSE IF (e < 0) THEN
        PRINT *,"ERROR in 1D channel node assignement,"
        PRINT *,"nodes",channel(i,i0)," or ",channel(i+1,i0),"(or both)."
        PRINT *,"Nodes out of range.  Please revise your channel link data."
        CALL byebye('Program SToRM stopped.')
      END IF
      chanedges(i,i0) = e
      length = length + edges(e)%length
    END DO
  END DO

! Set the array with the source and sink triangles.  Note that this step
! requires knowing the bed elevation at the center of the triangles, which
! means that this step must be computed in "dataprep" for fixed bed runs,
! but from within the solver for movable bed runs.
  ALLOCATE(chantrigs(0:d1,d2),chant(d1,d2),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(4*(d1 + 1)*d2 + 8*d1*d2)
  DO j = 1,n_chanbdr
    chantrigs(0,j) = chanedges(0,j)
    rlength = zero  ! Running length.
    DO i = 1,chanedges(0,j)
      t1 = edges(chanedges(i,j))%e(1)
      t2 = edges(chanedges(i,j))%e(2)
      IF (t1 < 0) THEN  ! In case the edge is a boundary edge.
        chantrigs(i,j) = t2
      ELSE IF (t2 < 0) THEN
        chantrigs(i,j) = t1
      ELSE  ! Choose the triangle with lowest elevation.
        IF (z(t1) < z(t2)) THEN
          chantrigs(i,j) = t1
        ELSE
          chantrigs(i,j) = t2
        END IF
      END IF

      ! Parametrize the channels' edges.
      rlength = rlength + half*(edges(chanedges(i,j))%length)
      IF (i > 1) rlength = rlength + half*(edges(chanedges(i-1,j))%length)
      chant(i,j) = rlength/length
    END DO
  END DO

END SUBROUTINE chan_arrays
