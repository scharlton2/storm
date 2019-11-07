SUBROUTINE zeta_bcs
  USE parameters
  USE dep_vars
  USE vbc_arrays
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Enforces stage at the appropriate boundary nodes.  It assumes that a call  !
!  to subroutine 'enforce_bcs' has already occured so that the values of      !
!  array 'hvtx()' are correct.                                                !
!                                                                             !
!  Francisco Simoes, April 2008                                               !
!  Last updated (mm-dd-yyyy): 02-27-2011 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables.
  INTEGER :: i,i0,j

  DO i0 = 1,n_outflowbdr
    IF (htype(i0) /= 1) CYCLE
    DO i = 1,n_hbc(i0)
      j = hbc_nodes(i,i0)
      zetavtx(j) = hvtx(j) + zvtx(j)
    END DO
  END DO

END SUBROUTINE zeta_bcs
