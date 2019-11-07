!-----------------------------------------------------------------------------!
!                                                                             !
!  Parameters for integration of one-dimensional model in SToRM.              !
!  All units are metric.                                                      !
!                                                                             !
!                                                                             !
!  F. Simoes, December 2011                                                   !
!  Last updated (mm-dd-yyyy): 09-05-2012                                      !
!                                                                             !
!-----------------------------------------------------------------------------!

MODULE oned_grid
  IMPLICIT NONE
  SAVE

! Machine precision:
  INTEGER, PARAMETER :: mp = KIND(1.0D0)  ! = KIND(1.0) for single precision.

! Points:
  TYPE :: point
    REAL (KIND=mp) :: x,y  ! Point coordinates.
  END TYPE point

! Basic geometry data.
  INTEGER :: d1_nxs  ! Number of cross sections.
  ! Number of points defining each cross section:
  INTEGER, ALLOCATABLE, DIMENSION(:) :: d1_xspts
  ! Location of cross sections (georeferenced to same datum as the SToRM data):
  REAL (KIND=point), ALLOCATABLE, DIMENSION(:) :: d1_xsecloc
  ! Cross section points, transverse and vertical coordinates.  d1_xscoordy
  ! must be georeferenced to same datum as the SToRM data.
  REAL (KIND=mp), ALLOCATABLE, DIMENSION(:,:) d1_xscoordx,d1_xscoordy

! Numerical parameters.
  REAL (KIND=mp) :: d1_delt, &   ! Time step size.
                    d1_theta, &  ! Time weighting factor.
                    d1_maxit     ! Maximum number of iterations per time step.

! Dependent variables.
  REAL (KIND=mp), ALLOCATABLE, DIMENSION(:) :: &
                  d1_wse, &    ! Stage for each cross section.
                  d1_disch, &  ! Discharge for each cross section.
                  d1_rough     ! Roughness for each cross section.

END MODULE oned_grid

! Usage:
!
! d1_xspts(i), i = 1,...,d1_nxs
! d1_xsecloc(i), i = 1,...,d1_nxs
! d1_xscoordx(i,j), i = 1,...,d1_xspts(j), j = 1,...,d1_nxs
! d1_xscoordy(i,j), i = 1,...,d1_xspts(j), j = 1,...,d1_nxs
! d1_wse(i), i = 1,...,d1_nxs
! d1_disch(i), i = 1,...,d1_nxs
! d1_rough(i), i = 1,...,d1_nxs
!
! Cross sections run from upstream to downstream, i.e., cross section #1 is the
! upstream-most cross section (where d1_disch(1) is the boundary condition) and
! cross section #d1_nxs is the downstream-most cross section (where
! d1_wse(d1_nxs) is the boundary condition).
