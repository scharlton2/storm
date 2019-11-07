SUBROUTINE enforce_bcs(h0,u0,v0,z0,ndim)
  USE constants
  USE geometry
  USE dep_vars
  USE options
  USE vbc_arrays
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Enforce stage and discharge boundary conditions, plus the solid wall       !
!  boundary conditions for the non-slip (btype = 0) option.                   !
!                                                                             !
!  INPUT:                                                                     !
!    h0    water depth   at vertices of triangles;                            !
!    u0    U-velocity    "     "     "      "    ;                            !
!    v0    V-velocity    "     "     "      "    ;                            !
!    z0    bed elevation "     "     "      "    ;                            !
!    ndim  dimension of arrays (used because the bounds checking option of    !
!          many compilers do not work for arrays of unknown length).          !
!                                                                             !
!  Francisco Simoes, May 2005                                                 !
!  Last updated (mm-dd-yyyy): 04-18-2013 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  INTEGER, INTENT(IN) :: ndim
  REAL (KIND=mp), INTENT(INOUT), DIMENSION(ndim) :: h0,u0,v0,z0

! Local variables.
  INTEGER :: i,i0,j
  REAL (KIND=mp) :: dotprod,hnot,incr,q
  REAL (KIND=mp), EXTERNAL :: vbc_by_h

! Enforce free surface elevation at the selected nodes.
  DO i0 = 1,n_outflowbdr
    SELECT CASE (htype(i0))

    CASE (1)  ! Subcritical outflow boundary, stage specified.
      DO i = 1,n_hbc(i0)
      ! In this implementation, the prescribed stage at the boundary nodes is
      ! enforced in a gradual manner.  The relaxation parameter is 'alpha' and
      ! 'kappa' is the minimum increment allowed.  Set kappa = [0.001,0.01] for
      ! laboratory channels, and kappa = [0.01,0.1] for larger channels and
      ! natural rivers. alpha = [0,1].
      j = hbc_nodes(i,i0)
      hnot = MAX(zero,hbc(i0) - z0(j))
      incr = alpha*(hnot - h0(j))
      IF (ABS(incr) < kappa) incr = hnot - h0(j)  ! Avoid infinitesimal steps.
      h0(j) = h0(j) - incr

      ! In the following code, the prescribed h is enforced instantaneously.
      j = hbc_nodes(i,i0)
      h0(j) = MAX(zero,hbc(i0) - z0(j))
    END DO

      ! Prevent flow reversal into the domain.  This is a numerical valve that
      ! sets the velocity to zero whenever the velocity vector points into the
      ! domain.  It is identical to what is done at the free boundary nodes.
      IF (hvalve) THEN
        DO i = 1,n_hbc(i0)
          dotprod = u0(hbc_nodes(i,i0))*hnormals(i,i0)%x + &
                    v0(hbc_nodes(i,i0))*hnormals(i,i0)%y
          IF (dotprod < zero) THEN
            u0(hbc_nodes(i,i0)) = zero
            v0(hbc_nodes(i,i0)) = zero
          END IF
        END DO
      END IF

    CASE (0)  ! Free flow boundary.

      ! Enforce outflow at freeboundary nodes.  This is simply a numerical
      ! valve that prevents flux into the domain by setting the node velocity
      ! to zero when inflow is detected.
      DO i = 1,n_hbc(i0)
        dotprod = u0(hbc_nodes(i,i0))*hnormals(i,i0)%x + &
                  v0(hbc_nodes(i,i0))*hnormals(i,i0)%y
        IF (dotprod < zero) THEN
          u0(hbc_nodes(i,i0)) = zero
          v0(hbc_nodes(i,i0)) = zero
        END IF
      END DO

    CASE (2)  ! Critical flow, Froude number equal to 1.

    ! Do nothing for now.

    CASE DEFAULT
      PRINT *,''
      PRINT *,'ERROR: invalid value in htype.'
      CALL byebye('Program SToRM stopped.')

    END SELECT
  END DO

! Enforce non-slip boundary conditions at solid walls.
  IF (btype == 0) THEN
    DO i = 1,n_wall
      u0(wall_pts(i)) = zero
      v0(wall_pts(i)) = zero
    END DO

  ELSE
    ! Make sure velocity vectors are tangent to solid walls.
    DO i = 1,n_wall
      dotprod = u0(wall_pts(i))*wtangs(i)%x + v0(wall_pts(i))*wtangs(i)%y
      u0(wall_pts(i)) = dotprod*wtangs(i)%x
      v0(wall_pts(i)) = dotprod*wtangs(i)%y
    END DO

  END IF

! Enforce zero velocity over dry nodes.  Dry nodes are set with zero velocity,
! but not with zero water depth.  The water depth h0() is kept as it is to
! preserve conservation of mass when wetting-drying cycles take place.
  DO i = 1,n_pts
    IF (h0(i) < h_dry) THEN
      u0(i) = zero
      v0(i) = zero
    END IF
  END DO

  DO i0 = 1,n_inflowbdr
    ! Inflow boundary conditions.
    SELECT CASE (vtype(i0))

    CASE (0)
      ! Enforce specified velocity at the selected nodes, if they are wet; set
      ! to zero the dry nodes (which was already done by part of the code
      ! above).
      DO i = 1,n_qin(i0)
        IF (h0(qin_nodes(i,i0)) > h_dry) THEN
          u0(qin_nodes(i,i0)) = qin(i0)*pnormals(i,i0)%x
          v0(qin_nodes(i,i0)) = qin(i0)*pnormals(i,i0)%y
        END IF
      END DO

    CASE (1,5)
      ! Enforce discharge along the specified chain of grid edges according to
      ! flow depth.
      q = vbc_by_h(i0,h0,ndim)
      DO i = 1,n_qin(i0)
        IF (h0(qin_nodes(i,i0)) > h_dry) THEN
          u0(qin_nodes(i,i0)) = v_inflow(i)*pnormals(i,i0)%x
          v0(qin_nodes(i,i0)) = v_inflow(i)*pnormals(i,i0)%y
        END IF
      END DO

    CASE (2)
      ! Enforce discharge along the specified chain of grid edges according to
      ! conveyance.
      PRINT *,"VTYPE = 2: option not implemented yet."
      CALL byebye('Program SToRM stopped.')

    CASE (3)
      ! Enforce velocity vector along the specified chain.
      DO i = 1,n_qin(i0)
      IF (h0(qin_nodes(i,i0)) > h_dry) THEN
          u0(qin_nodes(i,i0)) = qin(i0)
          v0(qin_nodes(i,i0)) = qin2(i0)
        END IF
      END DO

    CASE (4)
      ! This is for when stage is specified at the subcritical inflow. In this
      ! case do nothing.

    CASE DEFAULT
      PRINT *,''
      PRINT *,'ERROR: invalid value in vtype.'
      CALL byebye('Program SToRM stopped.')

    END SELECT
  END DO

END SUBROUTINE enforce_bcs
