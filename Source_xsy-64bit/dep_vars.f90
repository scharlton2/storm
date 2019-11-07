!-----------------------------------------------------------------------------!
!                                                                             !
!  This file contains the dependent flow variables used in SToRM.             !
!                                                                             !
!  F. Simoes, March 2003                                                      !
!  Last updated (mm-dd-yyyy): 05-02-2013 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

MODULE dep_vars
  USE parameters
  IMPLICIT NONE
  SAVE

! Main flow variables and respective gradients.
  REAL (KIND=mp), ALLOCATABLE, DIMENSION(:) :: h,hvtx,u,uvtx,v,vvtx,z,z2, &
                                               zeta,zetavtx
  TYPE(vector), ALLOCATABLE, DIMENSION(:) :: delh,delu,delv,delz,delz2

! Used for storing the bed elevations.  Needed when using the bed adjustment
! technique of Brufau and Garcia-Navarro (2003).
  REAL (KIND=mp), ALLOCATABLE, DIMENSION(:) :: zstore,zvtxstore

! Bed elevation arrays.
  REAL (KIND=mp), ALLOCATABLE, DIMENSION(:) :: zbx,zby,zvtx

! Friction coefficients, including extra parameters for vegetated flow
! resistance models.
  REAL (KIND=mp), ALLOCATABLE, DIMENSION(:) :: cd,cdvtx,rough
  REAL (KIND=mp), ALLOCATABLE, DIMENSION(:) :: rcoef1,rcoef2

! Wind forcing terms.
  REAL (KIND=mp), ALLOCATABLE, DIMENSION(:) :: wcoef1,wcoef2
  ! Temporary variables for reading data from CGNS file.
  REAL (KIND=mp), ALLOCATABLE, DIMENSION(:) :: w_dir,w_fric,w_mag

! Turbulent eddy viscosity.
  REAL (KIND=mp), ALLOCATABLE, DIMENSION(:) :: ev

! Cell velocities and celerities.
  REAL (KIND=mp), ALLOCATABLE, DIMENSION(:,:) :: u_avg
  REAL (KIND=mp), ALLOCATABLE, DIMENSION(:) :: u_mag,wavec

! Convective flux through each computational edge.
  REAL (KIND=mp), ALLOCATABLE, DIMENSION(:,:) :: flux

! Source terms due to bed and friction slope.
  REAL (KIND=mp), ALLOCATABLE, DIMENSION(:,:) :: source

! Cell residual.
  REAL (KIND=mp), ALLOCATABLE, DIMENSION(:,:) :: residual

! Nodal residual (resulting from applying the residual distribution scheme to
! the cell residual in the RDS solver) or element residual (resulting from
! adding all the edge and source contributions to each triangle in the FVT
! solver).
  REAL (KIND=mp), ALLOCATABLE, DIMENSION(:,:) :: phi,phi0

! Boundary conditions.
  INTEGER :: n_inflowbdr,n_outflowbdr,n_wall
  INTEGER, ALLOCATABLE, DIMENSION (:) :: n_hbc,n_timeserh,n_timeserq,n_qin, &
                                         wall_pts
  INTEGER, ALLOCATABLE, DIMENSION (:,:) :: hbc_nodes,qin_nodes
  REAL (KIND=mp), ALLOCATABLE, DIMENSION (:) :: hbc,qin,qin2
  REAL (KIND=mp), ALLOCATABLE, DIMENSION (:,:) :: timeloch,timelocq,timeserh, &
                                                  timeserq,timeserq2

! Culvert computations.
  INTEGER :: nculvert  ! Number of culverts.
  ! Location (triangle number) of inflow and outflow.
  INTEGER, ALLOCATABLE, DIMENSION (:) :: cvtrigin,cvtrigout
  ! Location (spatial coordinates) of inflow and outflow.
  REAL (KIND=mp), ALLOCATABLE, DIMENSION (:,:) :: cvptoin,cvptoout
  ! Performance curve for each culvert: number of rows in table.
  INTEGER, ALLOCATABLE, DIMENSION (:) :: ncvtable
  ! Performance curve for each culvert: headwater elevation and discharge.
  REAL (KIND=mp), ALLOCATABLE, DIMENSION (:,:) :: cvhead,cvq
  ! Vector pointing from center of triangle to location of culvert inlet.
  REAL (KIND=mp), ALLOCATABLE, DIMENSION (:,:) :: cvrin
  ! Source term for continuity equation.
  REAL (KIND=mp), ALLOCATABLE, DIMENSION (:) :: cvsrc

! TEMPORARY boundary condition type for flood inundation in South Fork river &
! Buckeye Lake, Columbus, OH.
  INTEGER :: n_chanbdr
  INTEGER, ALLOCATABLE, DIMENSION (:) :: n_channel
  INTEGER, ALLOCATABLE, DIMENSION (:,:) :: channel

! rmcd mod Boundary Conditions
  INTEGER, ALLOCATABLE, DIMENSION (:) :: inflow_type, outflow_type, &
                                         inflow_td_type, outflow_td_type, &
                                         inflow_td_index, outflow_td_index
  REAL (KIND=mp), ALLOCATABLE, DIMENSION (:) :: inflow_constval, &
                                                outflow_constval

! rmcd mod Time-dependent Vals
  INTEGER, ALLOCATABLE, DIMENSION (:) :: RCNumVals, TSNumVals, RCType, TSType, TSIndexUsed
  REAL, ALLOCATABLE, DIMENSION (:) :: RCVals, TSVals

! These variables store the results from various forms of computing the inflow
! and outflow discharge, for debugging.
  REAL (KIND=mp), ALLOCATABLE, DIMENSION (:) :: cqin,cqin2,cqin3,cqout, &
                                                cqout2,cqout3

! Quantities related to the bed gradient necessary to solve the 2D kinematic
! wave equation. magzb is the magnitide of the bed gradient and ngzv is the
! unit vector in the direction of grad zb.
  REAL (KIND=mp), ALLOCATABLE, DIMENSION(:) :: magzb
  TYPE(vector), ALLOCATABLE, DIMENSION(:) :: ngzb

! `A posteriori error estimate.
  REAL (KIND=mp), ALLOCATABLE, DIMENSION(:) :: e_est

! Temporary arrays to facilitate reading data from Tecplot-formated ASCII data
! files.  NOTE: datav is a logical array that determines if the data in the
! variables read from the Tecplot file is located at cell vertices (.TRUE.) or
! at cell centers (.FALSE.).  Its dimension is the number of variables in the
! Tecplot file.
  LOGICAL :: datav(11)  ! 11 variables used by STORM.
  REAL (KIND=mp), ALLOCATABLE, DIMENSION(:) :: cd_tcp,h_tcp,u_tcp,v_tcp,z_tcp

END MODULE dep_vars
