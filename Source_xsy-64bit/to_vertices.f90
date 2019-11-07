SUBROUTINE to_vertices
  USE constants
  USE geometry
  USE dep_vars
  USE options
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This subroutine replaces the code documented in file                       !
!  solverFVT.document.f90, which used to do the interpolation of the          !
!  solution variables from the cell centers to the cell vertices.             !
!                                                                             !
!  Francisco Simoes, July 2012                                                !
!  Last updated (mm-dd-yyyy): 06-02-2017 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables.
  INTEGER :: i,j,k
  REAL (KIND=mp) :: zetaref
  LOGICAL :: flag

! First standard interpolation.
  !WHERE (drycell) h = zero
  CALL ctr2vtxWSE(nctr2vtx,zeta,zetavtx)
  hvtx = zetavtx - zvtx

! Limit vertex water surface elevation to the highest center value of any
! neighboring triangles.
  DO i = 1,n_pts
    zetaref = zvtx(i)
    DO j = 2,n2t(i,1) + 1
      k = n2t(i,j)
      IF (drycell(k)) CYCLE
      zetaref = MAX(zeta(k),zetaref)
    END DO
    zetavtx(i) = MIN(zetavtx(i),zetaref)
    hvtx(i) = zetavtx(i) - zvtx(i)  ! Corrected water depth.
  END DO
  CALL positivity_vtx  ! Make sure water depth is >= zero.

! Velocity interpolation.  Set velocity to zero in all partially-wet cells.
  DO i = 1,n_elems
    IF (drycell(i)) CYCLE
    flag = .TRUE.
    DO j = 1,3
      k = grid(i)%vertex(j)
      flag = flag .AND. (hvtx(k) > h_dry)
    END DO
    IF (.NOT. flag) THEN
      u(i) = zero
      v(i) = zero
    END IF
  END DO
  CALL ctr2vtx(nctr2vtx,u,uvtx)
  CALL ctr2vtx(nctr2vtx,v,vvtx)
  IF (btype == 0) CALL wetdryvfix  ! Enforce zero velocity at walls.
  WHERE (hvtx < h_dry)  ! Enforce zero velocity at dry vertices.
    uvtx = zero
    vvtx = zero
  END WHERE

! Enforce boundary conditions at inflow/outflow strings.
  CALL enforce_bcs(hvtx,uvtx,vvtx,zvtx,n_pts)
  !CALL zeta_bcs

END SUBROUTINE to_vertices
