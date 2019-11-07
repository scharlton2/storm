SUBROUTINE elemnt_db
  USE geometry
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This subroutine prepares element-oriented database arrays.  For example,   !
!  it initializes the table that connects each element to its neighbors for   !
!  the least-squares gradient computations.                                   !
!                                                                             !
!  IMPORTANT NOTE: at present, table t2t is set-up for least-squares          !
!  computation of the gradients.  Consequently, it must be viewed as a        !
!  computational molecule and not as a true connectivity table.  The          !
!  computational molecule set up in this subroutine has 3 or more nodes for   !
!  each element.                                                              !
!                                                                             !
!  Table t2t3 does contain a true triangle-to-triangle connectivity table     !
!  that gives, for each triangle i, all the triangles connected to it via     !
!  the edge j (j=1,2,3).  t2t3(j,i) = -1 when edge j is a bounday edge.       !
!                                                                             !
!  F. Simoes, March 2007                                                      !
!  Last updated (mm-dd-yyyy): 08-15-2016 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables.
  INTEGER :: e1,e2,i,ierror,j,k,l,m,maxdim,n,nn
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: t2tw  ! Temporary working array.

! The code below requires that the computational domain has AT LEAST 4
! elements.
  IF (n_elems < 4) THEN
    WRITE (*,'("ERROR: too few triangles in computational domain.",/, &
      & "Number of elements (triangles) in computational domain must be >= &
      &4.")')
    CALL byebye('Program SToRM stopped.')
  END IF

  WRITE (*,'("  Element-to-element connectivity table...",$)')

! ALLOCATE space for working array.
  maxdim = 20
  ALLOCATE(t2tw(n_elems,maxdim),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  t2tw = 0

! First sweep.  This sets-up the table for all the triangles inside the domain.
! Each triangle is linked to the 3 neighbors that have a common edge.

  DO i = 1,n_edges
    e1 = edges(i)%e(1)
    e2 = edges(i)%e(2)
    IF (e1 < 1 .OR. e2 < 1) CYCLE  ! Boundary edge.
    k = t2tw(e1,1) + 1
    t2tw(e1,k+1) = e2
    t2tw(e1,1) = k
    k = t2tw(e2,1) + 1
    t2tw(e2,k+1) = e1
    t2tw(e2,1) = k
  END DO

! Boundary triangles have less than three neighbors.  The connectivity table
! needs to be augmented in order to ensure appropriate rank of the least
! squares matrices.

  DO i = 1,n_elems
    n = t2tw(i,1)
    IF (n == 3) CYCLE  ! Connectivity is ok.
    IF (n < 1 .OR. n > 3) THEN  ! Error condition.
      WRITE (*,'("PROGRAM ERROR: triangle connectivity table incorrect.",/, &
           & "Please check your mesh connectivity data.")')
      DEALLOCATE(t2tw)
      CALL byebye('Program stopped.')
    END IF

! In the following DO-loop, i is the triangle of interest, j are the closest
! neighbors of i (i.e., the triangles that share one edge with triangle i), and
! k are the closest neighbors of j.
    DO l = 1,n
      j = t2tw(i,l+1)  ! Neighboring triangle.
      DO m = 1,MIN(3,t2tw(j,1))
        k = t2tw(j,m+1)
        IF (k == i) CYCLE
        nn = t2tw(i,1) + 1  ! Augment connectivity table with triangle k.
        t2tw(i,nn+1) = k
        t2tw(i,1) = nn
      END DO
    END DO

    ! Clean-up duplicate entries.
    CALL cde_t2t(t2tw,i,n_elems,maxdim)
  END DO

! In very particular cases, the above augmentation may not be enough.  This may
! happen in long channels with one triangle in width, for example.

  DO i = 1,n_elems
    n = t2tw(i,1)
    IF (n >= 3) CYCLE  ! Connectivity is ok.
    ! Augment computational cell using the last entry in the t2t table.
    j = t2tw(i,n+1)
    DO m = 1,t2tw(j,1)
      k = t2tw(j,m+1)
      IF (k == i) CYCLE
      nn = t2tw(i,1) + 1  ! Augment connectivity table with triangle k.
      t2tw(i,nn+1) = k
      t2tw(i,1) = nn
    END DO

    ! Clean-up duplicate entries (again).
    CALL cde_t2t(t2tw,i,n_elems,maxdim)
  END DO

! Set-up the final t2t table and clean-up t2tw.
  maxdim = 0
  DO i = 1,n_elems  ! Find maximum dimension of array.
    maxdim = MAX(maxdim,t2tw(i,1))
  END DO
  maxdim = maxdim + 1
  ALLOCATE(t2t(n_elems,maxdim),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(4*n_elems*maxdim)
  DO i = 1,n_elems
    DO j = 1,maxdim
      t2t(i,j) = t2tw(i,j)
    END DO
  END DO
  DEALLOCATE(t2tw)

! Prepare the t2tHD connectivity table.  The t2tHD table contains all the first
! neighbors of each triangle.  Note that t2tHD is built using the n2t table
! (the n2t table contains strictly the node-to-triangle connectivity table).
! Pathologic cases may appear at corners and long channels with one triangle in
! width, for example.  In those cases, t2tHD is augmented just like t2t is.

! ALLOCATE space for working array.
  maxdim = 0
  DO i = 1,n_pts
    maxdim = MAX(maxdim,n2t(i,1))
  END DO
  maxdim = 4*maxdim + 3
  ALLOCATE(t2tw(n_elems,maxdim),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  t2tw = 0

  DO i = 1,n_elems
    DO j = 1,3
      k = grid(i)%vertex(j)  ! Vertex being processed.

      DO l = 2,n2t(k,1) + 1
        m = n2t(k,l)  ! Triangle connected to vertex k.
        IF (m == i) CYCLE
        t2tw(i,t2tw(i,1)+2) = m
        t2tw(i,1) = t2tw(i,1) + 1
      END DO
    END DO

    ! Clean-up duplicate entries.
    CALL cde_t2t(t2tw,i,n_elems,maxdim)
  END DO

! Boundary triangles may have less than three neighbors.  The connectivity
! table needs to be augmented in order to ensure appropriate rank of the least
! squares matrices.

  DO i = 1,n_elems
    n = t2tw(i,1)
    IF (n > 2) CYCLE  ! Connectivity is ok.
    IF (n < 1) THEN  ! Error condition.
      WRITE (*,'("PROGRAM ERROR: triangle connectivity table incorrect.",/, &
           & "Please check your mesh connectivity data.")')
      DEALLOCATE(t2tw)
      CALL byebye('Program stopped.')
    END IF

    nn = MIN(3,t2tw(i,1) + 1)
    DO j = 2,nn  ! Elements used to expand the computational molecule.
      DO k = 1,3
        l = grid(j)%vertex(k)  ! Vertex being processed.

        DO m = 2,n2t(l,1) + 1
          n = n2t(l,m)  ! Triangle connected to vertex l.
          IF (n == i .OR. m == j) CYCLE
          IF (t2tw(i,1)+2 > maxdim) THEN
            WRITE (*,'("PROGRAM ERROR: allocation bounds in T2TW.",/, &
              & "Please report the error to the SToRM programming team.")')
            DEALLOCATE(t2tw)
            CALL byebye('Program stopped.')
          END IF
          t2tw(i,t2tw(i,1)+2) = n
          t2tw(i,1) = t2tw(i,1) + 1
        END DO
      END DO
    END DO

    ! Clean-up duplicate entries (again).
    CALL cde_t2t(t2tw,i,n_elems,maxdim)
  END DO

! Set-up the final t2tHD table and clean-up t2tw.
  maxdim = 0
  DO i = 1,n_elems  ! Find maximum dimension of array.
    maxdim = MAX(maxdim,t2tw(i,1))
  END DO
  maxdim = maxdim + 1
  ALLOCATE(t2tHD(n_elems,maxdim),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(4*n_elems*maxdim)
  DO i = 1,n_elems
    DO j = 1,maxdim
      t2tHD(i,j) = t2tw(i,j)
    END DO
  END DO
  DEALLOCATE(t2tw)

  ALLOCATE(t2t3(3,n_elems),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(4*n_elems*3)

! Array t2t3 gives the three next-of-kin triangles, i.e., the triangles that
! share an edge with the triangle of interest.  In other words, t2t3(j,i)
! points to the triangle that shares local edge j (j=1,2,3) with the triangle
! of interest (triangle i).  It takes the value -1 (the "no element" state)
! when the edge is a boundary edge.

  DO i = 1,n_elems
    DO j = 1,3
      k = ABS(grid(i)%edge(j))
      e1 = edges(k)%e(1)
      IF (e1 == i) e1 = edges(k)%e(2)
      t2t3(j,i) = e1  ! Will be < 1 for boundary edges.
    END DO
  END DO

  PRINT *,'Done.'

END SUBROUTINE elemnt_db
