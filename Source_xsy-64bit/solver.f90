SUBROUTINE solver(iter,ttime,rmax,fn,title) !rmcd mod
  USE parameters
  USE options
  USE io
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This subroutines manages the solution method used.  It takes care of       !
!  deciding what solver is used and prepares the data for it.  On output it   !
!  passes the number of iterations performed (iter) and the value of the      !
!  maximum residual (rmax).                                                   !
!                                                                             !
!  INPUT:                                                                     !
!    fn     filename for the output files;                                    !
!    title  title used in the Tecplot files.                                  !
!                                                                             !
!  OUTPUT:                                                                    !
!    iter   maximum iteration before stopping criteria was reached;           !
!    rmax   maximum residual or other value used in the stopping criteria     !
!           (i.e., quantity that is checked for to determine convergence).    !
!                                                                             !
!  Francisco Simoes, February 2007                                            !
!  Last updated (mm-dd-yyyy): 04-16-2008 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables.
  INTEGER, INTENT(OUT) :: iter
  REAL(KIND=mp), INTENT(OUT) :: rmax, ttime !rmcd mod
  CHARACTER (LEN=*), INTENT(IN) :: fn,title

! Local variables.
  INTEGER, EXTERNAL :: mem_call

  SELECT CASE (opt_solver)
  CASE (1)  ! RDS solver.
    WRITE (*,'("  Preparing solver-dependent data...",$)')
    CALL dataprepRDS
    PRINT *,'Done.'
    WRITE (*,'(A,ES13.5)')"  Estimated memory allocated (MB):", &
      mem_call()/1048576.0_mp
    IF (output_head) WRITE (o_unit,'(//A,ES13.5)') &
      "Estimated memory allocated (MB):",mem_call()/1048576.0_mp
    PRINT *,' Solver started...'
    CALL solverRDS(iter,rmax,fn,title)

  CASE (2)  ! Finite volume solver using triangles for control volumes (FVT).
    WRITE (*,'("  Preparing solver-dependent data...",$)')
    CALL dataprepFVT
    PRINT *,'Done.'
    WRITE (*,'(A,ES13.5)')"  Estimated memory allocated (MB):", &
      mem_call()/1048576.0_mp
    IF (output_head) WRITE (o_unit,'(//A,ES13.5)') &
      "Estimated memory allocated (MB):",mem_call()/1048576.0_mp
    PRINT *,' Solver started...'
    CALL solverFVT(iter,ttime,rmax,fn,title) !rmcd mod ('ttime')

  CASE DEFAULT
    PRINT *,' STORM error: solver type is not implemented in the'
    PRINT *,' current version.'
    CALL byebye('Program stopped.')

  END SELECT

END SUBROUTINE solver
