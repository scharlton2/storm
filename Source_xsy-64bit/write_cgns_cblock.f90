!-----------------------------------------------------------------------------!
!                                                                             !
!  This file contains the variables used in SUBROUTINEs "write_cgns" and      !
!  "write_cgns_cell.  Variables with 'cell' in their names belong to the      !
!  latter SUBROUTINE. These variables are local to these subroutines, but     !
!  this implementation was done to try to fix a stack overflow problem        !
!  suffered by this subroutine in cases with very large datasets.             !
!                                                                             !
!  F. Simoes, 22 August 2014                                                  !
!  Last updated (mm-dd-yyyy): 06-02-2017 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

MODULE write_cgns_cblock
  USE parameters
  IMPLICIT NONE
  SAVE

  INTEGER, ALLOCATABLE, DIMENSION(:) :: ibc,ibcell
  REAL(8), ALLOCATABLE, DIMENSION(:) :: depth,dShields,dtau,dtauv,tau,taubx, &
                                        tauby,taucell,tauv,wselev
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: wetv
  TYPE(vector), ALLOCATABLE, DIMENSION(:) :: deltaubx,deltauby

END MODULE write_cgns_cblock
