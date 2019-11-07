!-----------------------------------------------------------------------------!
!                                                                             !
!  This file contains the global variables used in the description of the     !
!  grid system used in SToRM.                                                 !
!                                                                             !
!  F. Simoes, March 2003                                                      !
!  Last updated (mm-dd-yyyy): 08-22-2012 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

MODULE geometry
  USE parameters
  IMPLICIT NONE
  SAVE

! Computational nodal coordinates.
  INTEGER :: n_pts  ! Number of points.
  TYPE(point), ALLOCATABLE, DIMENSION(:) :: nodes

! Computational mesh.
  INTEGER :: n_elems  ! Number of triangles.
  TYPE(triangle), ALLOCATABLE, DIMENSION(:) :: grid

! Triangle edge structure.
  INTEGER :: n_edges  ! Overall number of edges, i.e., of triangle sides.
  TYPE(edge), ALLOCATABLE, DIMENSION(:) :: edges

! Vectors pointing from center of triangle to midpoint of each edge.
  TYPE(vector), ALLOCATABLE, DIMENSION(:,:) :: rc

! Node and element data base.
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: n2t,n2t2  ! Node-to-triangle tables.
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: n2n  ! Node-to-node table.

  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: t2t, & ! Element-to-element tables.
                                          t2tHD,t2t3
  REAL (KIND=mp), ALLOCATABLE, DIMENSION(:) :: cv_area  ! Control volume areas.
  REAL (KIND=mp), ALLOCATABLE, DIMENSION(:) :: cv_perim  !  "       "    perim.

! Edge classification arrays (wall edges, flow edges, 1D channel edges, etc.).
  INTEGER :: n_bcedges,n_bpolygon,flow_edges,flow_edges1,wall_edges,wall_edges1
  INTEGER, ALLOCATABLE, DIMENSION(:) :: bcedges,bpolygon,flowedg,flowedg1, &
                                        walledg,walledg1
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: bcedgesinflow,bcedgesoutflow, &
                                          chanedges,chantrigs
  REAL (KIND=mp), ALLOCATABLE, DIMENSION(:,:) :: chant

! Mesh-based arrays needed for least-squares gradient construction.
  REAL (KIND=mp), ALLOCATABLE, DIMENSION(:) :: ai,bi,ci,det
  REAL (KIND=mp), ALLOCATABLE, DIMENSION(:,:) :: dx,dy,lsqweight

! Mesh-based arrays needed for cell-to-vertex interpolation.
  REAL (KIND=mp), ALLOCATABLE, DIMENSION(:,:) :: c2v_weights,c2v_weights_id

! Masking arrays that indicate if each cell is dry (.TRUE.) or wet (.FALSE.).
! 'drycell' indicates if the entire cell is wet or dry and 'partdy' concerns
! only partially wet cells.  'partdry' is .TRUE. for partially wet cells and is
! .FALSE. for totally wet or totally dry cells.  'bdrycell' indicates wet cells
! that are located at the wet/dry boundary.
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: bdrycell,drycell,partdry

! Array containing links to the vertices of the triangles.  The vertex
! information is sorted by elevation within each triangle, p(1) being the
! lowest vertex and p(3) being the highest.  Array 'areap' contains the area
! of each triangle divided in two parts and is used to compute the effective
! area in partially wet triangles.
  TYPE(vector3), ALLOCATABLE, DIMENSION(:) :: csortedz
  REAL (KIND=mp), DIMENSION(:,:), ALLOCATABLE ::areap

! Mesh-based arrays for gradient computation based on the Gaussian method.
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: gconn
  REAL(KIND=mp), ALLOCATABLE, DIMENSION(:,:) :: gcoefx,gcoefy,gweight
  TYPE(vector), ALLOCATABLE, DIMENSION(:) :: gwgrad  ! Working array...

END MODULE geometry
