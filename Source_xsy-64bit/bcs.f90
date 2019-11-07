SUBROUTINE bcs(funit)
  IMPLICIT NONE

! Dummy variables:
  INTEGER, INTENT(IN) :: funit

! Read inflow/outflow and other boundary data.
  CALL ioflow(funit)

! Find boundary edges.
  CALL walls

! Set-up arrays for inflow boundary condition computations.
  CALL ibc_arrays

END SUBROUTINE bcs
