SUBROUTINE write_cgns(ttime)
  USE geometry
  USE dep_vars
  USE options
  USE io
  USE constants
  USE write_cgns_cblock
  IMPLICIT NONE
  INCLUDE "Header Files\cgnslib_f.h"
!!DEC$ ATTRIBUTES REFERENCE, C, VARYING :: cg_goto_f
!!DEC$ ATTRIBUTES REFERENCE, C, VARYING :: cg_array_write_f

!-----------------------------------------------------------------------------!
!                                                                             !
!  This SUBROUTINE writes solution data to the CGNS datafile used by iRIC     !
!  2.0.  The dummy argument variable is the time step number.                 !
!                                                                             !
!  Francisco Simoes, 20 May 2012                                              !
!  Last updated (mm-dd-yyyy): 06-02-2017 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy argument.
  REAL (KIND=mp), INTENT(IN) :: ttime

! Local variables.
  INTEGER :: celldim,fsol_id,fsol_iter,i,ierror,isize(3,3),j,k,nbases, &
             nzones,physdim
  INTEGER :: idata(2)
  REAL(8) :: mb,vmag
  CHARACTER(LEN=5) :: tsindex
  CHARACTER(LEN=32) :: basename,tmp_flowsol2D,zonename

  TSCellPrintCount = TSCellPrintCount + 1
  WRITE (tsindex,'(I5)') TSCellPrintCount
  tsindex = ADJUSTL(tsindex)
  tmp_flowsol2D = 'FlowCellSolution'//TRIM(tsindex)
  tmp_flowsol2D = TRIM(tmp_flowsol2D)
  CellSolNames(TSCellPrintCount) = tmp_flowsol2D
  CellTimeIncrements(TSCellPrintCount) = ttime

! ALLOCATE variable space the first time the subroutine is called.
  IF (.NOT. ALLOCATED(ibc)) THEN
    ALLOCATE(ibc(n_elems),wetv(n_elems),STAT=ierror)
    IF (ierror /= 0) CALL alloc_err(ierror)
    CALL mem_add(4*n_pts)
    ALLOCATE(depth(n_pts),dShields(n_pts),dtauv(n_pts),tauv(n_pts), &
      wselev(n_pts),STAT=ierror)
    IF (ierror /= 0) CALL alloc_err(ierror)
    CALL mem_add(8*n_pts*5)
    ALLOCATE(dtau(n_elems),tau(n_elems),taubx(n_elems),tauby(n_elems), &
      STAT=ierror)
    IF (ierror /= 0) CALL alloc_err(ierror)
    CALL mem_add(8*n_elems*4)
    ALLOCATE(deltaubx(n_elems),deltauby(n_elems),STAT=ierror)
    IF (ierror /= 0) CALL alloc_err(ierror)
    CALL mem_add(8*n_elems*4)
  END IF

  CALL correct_vtx(depth,wselev,wetv,n_pts)

  CALL cg_nbases_f(findex,nbases,ierror)
  DO i = 1,nbases
    CALL cg_base_read_f(findex,i,basename,celldim,physdim,ierror)
    IF (ierror /= 0) THEN
      CALL cg_error_print_f()
      CALL byebye('SToRM ERROR: cannot read CGNS base.')
    END IF
    IF (TRIM(basename) == 'iRIC') THEN  ! iRIC base name.
      CALL cg_nzones_f(findex,i,nzones,ierror)
      IF (ierror /= 0) THEN
        CALL cg_error_print_f()
        CALL byebye('SToRM ERROR: cannot find number of CGNS zones.')
      END IF
      CALL cg_biter_write_f(findex,i,'TimeIterValues',TSCellPrintCount,ierror)
      call cg_goto_f(findex,i,ierror,'BaseIterativeData_t',1,'end')
      call cg_array_write_f('TimeValues',RealDouble,1,TSCellPrintCount,CellTimeIncrements,ierror)
      DO j = 1,nzones
        CALL cg_zone_read_f(findex,i,j,zonename,isize,ierror)
        IF (ierror /= 0) THEN
          CALL cg_error_print_f()
          CALL byebye('SToRM ERROR: cannot read CGNS zone.')
        END IF
        IF (TRIM(zonename) == 'iRICZone') THEN  ! iRIC zone name.
          CALL cg_sol_write_f(findex,i,j,tmp_flowsol2D,CellCenter,fsol_iter,ierror)
          CALL cg_goto_f(findex,i,ierror,'Zone_t',j,'ZoneIterativeData_t',1,'end')
          idata(1)=32
          idata(2)=TSCellPrintCount
          CALL cg_array_write_f('FlowCellSolutionPointers',Character,2,idata,CellSolNames,ierror)

          ! Write velocity field.
          call cg_field_write_f(findex,i,j,fsol_iter,RealDouble,'VelocityX',u,fsol_id,ierror)
          IF (ierror /= 0) THEN
            CALL cg_error_print_f()
            CALL byebye('SToRM ERROR: error writing VelocityX to CGNS file.')
          END IF
          call cg_field_write_f(findex,i,j,fsol_iter,RealDouble,'VelocityY',v,fsol_id,ierror)
          IF (ierror /= 0) THEN
            CALL cg_error_print_f()
            CALL byebye('SToRM ERROR: error writing VelocityY to CGNS file.')
          END IF

          ! Mask dry nodes.
          ibc = 0
          DO k = 1,n_elems
            IF (h(k) > h_dry) ibc(k) = 1
          END DO
          call cg_field_write_f(findex,i,j,fsol_iter,Integer,'IBC',ibc,fsol_id,ierror)
          IF (ierror /= 0) THEN
            CALL cg_error_print_f()
            CALL byebye('SToRM ERROR: error writing IBC to CGNS file.')
          END IF
          ! Write water surface elevation and depth.
          !depth = 0.0d0
          !wselev = 0.0d0
          call cg_field_write_f(findex,i,j,fsol_iter,RealDouble,'Depth',h,fsol_id,ierror)
          IF (ierror /= 0) THEN
            CALL cg_error_print_f()
            CALL byebye('SToRM ERROR: error writing Depth to CGNS file.')
          END IF
          ! Write friction coefficients.
          call cg_field_write_f(findex,i,j,fsol_iter,RealDouble,'Cd',cd,fsol_id,ierror)
          IF (ierror /= 0) THEN
            CALL cg_error_print_f()
            CALL byebye('SToRM ERROR: error writing Cd to CGNS file.')
          END IF

        END IF
      END DO
    END IF
  ENDDO

END SUBROUTINE write_cgns
