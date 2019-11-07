SUBROUTINE correct_vtx(depth,wselev,wetv,n)
  USE parameters
  USE constants
  USE geometry
  USE dep_vars
  USE options
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Corrects water surface elevation of vertices located at the river          !
!  margins.                                                                   !
!                                                                             !
!  Francisco Simoes, 20 Sep 2012                                              !
!  Last updated (mm-dd-yyyy): 10-23-2012 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  INTEGER, INTENT(IN) :: n  ! Number of vertices (should be = n_pts).
  REAL(8), INTENT(OUT) :: depth(n),wselev(n)  ! Depth and wse at vertices.
  LOGICAL, INTENT(OUT) :: wetv(n)  ! .TRUE. if vertex is wet, .FALSE. if not.

! Local variables.
  INTEGER :: counter,i,j,k
  REAL (KIND=mp) :: a,id,vmag,zetaref
  !LOGICAL :: pdry(n_elems)  ! Part-dry cells.

! Initialize variables.
  DO k = 1,n_pts
    depth(k) = hvtx(k)
    wselev(k) = zetavtx(k)
  END DO

! Mask dry vertices
  wetv = .FALSE.
  DO k = 1,n_pts
    vmag = sqrt((uvtx(k)*uvtx(k))+(vvtx(k)*vvtx(k)))
    IF (hvtx(k) > h_wet) THEN
      wetv(k) = .TRUE.
    ELSE IF (hvtx(k) > h_dry .AND. vmag > 1.0e-5) THEN
      wetv(k) = .TRUE.
    END IF
  END DO

! Jump out of this SUBROUTINE here if the old approach is to be used.  The new
! approach doesn't use partially dry cells in the interpolation.
  !RETURN

! Find partially dry cells.
  !pdry = .FALSE.
  !DO k = 1,n_elems
  !  IF ((.NOT. wetv(grid(k)%vertex(1))) .OR. (.NOT. wetv(grid(k)%vertex(2))) &
  !    .OR. (.NOT. wetv(grid(k)%vertex(3)))) pdry(k) = .TRUE.
  !END DO

! Interpolation of water surface elevation using inverse distance wighting
! coefficients, which is right only if nctr2vtx = 3 in SUBROUTINE
! initialize_options, because the weights are computed at the data
! pre-processing stage to save CPU time.
  DO i = 1,n_pts
    id = zero
    a = zero
    counter = zero
    DO j = 2,n2t(i,1)+1
      k = n2t(i,j)
      IF (drycell(k)) CYCLE  ! Do not use dry-cell values...
      !IF (pdry(k)) CYCLE  ! ... nor partially-dry-cell values.
      IF (zeta(k) < zvtx(csortedz(k)%p(3))) CYCLE
      counter = counter + 1
      id = id + c2v_weights_id(i,j-1)
      a = a + c2v_weights_id(i,j-1)*zeta(k)
    END DO
    IF (counter > 0) THEN
      wselev(i) = a/id
    ELSE
      wselev(i) = zvtx(i)
    END IF
  END DO
  depth = wselev - zvtx

! Limit vertex water surface elevation to the highest center value of any
! neighboring triangles, just as it is done in SUBROUTINE to_vertices.
  DO i = 1,n_pts
    zetaref = zvtx(i)
    DO j = 2,n2t(i,1) + 1
      k = n2t(i,j)
      IF (drycell(k)) CYCLE
      zetaref = MAX(zeta(k),zetaref)
    END DO
    wselev(i) = MIN(wselev(i),zetaref)
    depth(i) = wselev(i) - zvtx(i)  ! Corrected water depth.
    IF (depth(i) > h_wet) THEN
      wetv(i) = .TRUE.
    ELSE
      wetv(i) = .FALSE.
    END IF
  END DO

END SUBROUTINE correct_vtx
