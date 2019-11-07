!-----------------------------------------------------------------------------!
!                                                                             !
!  File unit id's for READ and WRITE statements within SToRM.  Also contains  !
!  the necessary information to generate intermediary output solution files.  !
!                                                                             !
!  F. Simoes, April 2004                                                      !
!  Last updated (mm-dd-yyyy): 06-02-2017 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

MODULE io
  IMPLICIT NONE
  SAVE

! Unit id of the main output file.
  INTEGER :: o_unit

! Unit id of the file where the solution is stored, in Tecplot format.
  INTEGER :: s_unit

! Unit id of the file where probe output is stored.
  INTEGER :: p_unit

! Array containing the iteration information for which intermediate solution
! output is required.
  INTEGER, DIMENSION(:), ALLOCATABLE :: iout
  INTEGER :: n_iout
! Counter pointing to next element in array iout, used for fast searching.
  INTEGER :: p_iout
! Variable used for calling subroutine ifname.
  INTEGER :: max_iout

! Basic parameters for CGNS I/O:

! For CGNS output; .TRUE. if CGNS I/O is used.
  LOGICAL :: cgns

! CGNS file quantities.
  INTEGER :: findex

! For intermediate printing.
  INTEGER :: TSPrintCount      ! For cell vertex solutions.
  INTEGER :: TSCellPrintCount  ! For cell-center solutions.

! Unsteady solution output names in CGNS file.
  CHARACTER*32, ALLOCATABLE, DIMENSION(:) :: SolNames
  CHARACTER*32, ALLOCATABLE, DIMENSION(:) :: CellSolNames

  REAL(8), ALLOCATABLE, DIMENSION(:) :: TimeIncrements
  REAL(8), ALLOCATABLE, DIMENSION(:) :: CellTimeIncrements

! Path and name of file with CGNS hot-start data.
  CHARACTER (LEN=256) :: hsfile
  LOGICAL :: usehsfile  ! .TRUE. if CGNS hot start is used.

END MODULE io
