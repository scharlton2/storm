!-----------------------------------------------------------------------------!
!                                                                             !
!  This module contains the arrays used to implement several types of         !
!  velocity boundary conditions in SToRM.                                     !
!                                                                             !
!  F. Simoes, November 2005                                                   !
!  Last updated (mm-dd-yyyy): 03-12-2013 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

MODULE vbc_arrays
  USE parameters
  IMPLICIT NONE
  SAVE

! vtype    Type of inflow boundary
!--------------------------------------
!   0      Velocity magnitude, velocity normal to edge.
!   1      Discharge distributed according to flow depth, with all of inflow
!          edges considered wet.
!   2      Discharge distributed according to conveyance.
!   3      Full velocity vector.
!   4      Stage is specified, subcritical inflow.
!   5      Discharge distributed according to flow depth, but only using the
!          edges belonging to wet triangles.
!   6      Culverts, weirs, and other structures defined by a rating curve.
  INTEGER, ALLOCATABLE, DIMENSION (:) :: vtype

! Array pnormals contains the unit normals at each node of the inflow array,
! pointing into the computational domain.
  TYPE (vector), ALLOCATABLE, DIMENSION(:,:) :: pnormals

! The following arrays are not necessay and could do away with, but they make
! implementing the velocity boundary conditions much, much easier...
  ! Length of the edges containing the inflow.
  REAL (KIND=mp), ALLOCATABLE, DIMENSION(:,:) :: l_inflow
  REAL (KIND=mp), ALLOCATABLE, DIMENSION(:) :: &
    h_inflow, &  ! Depth of the nodes in the inflow boundary;
    a_inflow, &  ! flow area of the inflow element faces;
    u_inflow, &  ! flow velocity at the center of the inflow element faces; and
    v_inflow     ! velocity at the inflow nodes.

! Array with the different types of outflow boundary condition: htype = 0 for
! free boundaries; = 1 for subcritical outflow; = 2 for critical (Fr = 1)
! outflow.
  INTEGER, ALLOCATABLE, DIMENSION (:) :: htype

! Array that contains the unit normals at each node where the stage boundary
! condition is enforced.  The normals point out of the computational domain.
  TYPE (vector), ALLOCATABLE, DIMENSION(:,:) :: hnormals

! Array that contains the unit tangent vectors at each node of the solid walls.
  TYPE (vector), ALLOCATABLE, DIMENSION(:) :: wtangs,wtangs1

! Logical arrays to help enforcing the inflow/outflow boundary conditions.
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: hbar,qbar

! Array that contains the unit normal at each edge that is an inflow, outflow,
! or frr flow boundary.
!  TYPE(vector), ALLOCATABLE, DIMENSION(:) :: abcnormals

! This variable is used in a fail-safe pocedure when all the inflow triangles
! get dry.  It is set in SUBROUTINE 'dataprepFVT'.
  REAL (KIND=mp) :: bch_max

END MODULE vbc_arrays
