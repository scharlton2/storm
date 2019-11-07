SUBROUTINE get_time_bcs(time)
  USE parameters
  USE dep_vars
  USE vbc_arrays
  USE options
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This subroutine uses the time-dependent boundary values tables and         !
!  interpolates q and h for the current time.                                 !
!                                                                             !
!  INPUT:                                                                     !
!    time    physical time, in seconds.  The simulation always starts at      !
!            time zero.                                                       !
!                                                                             !
!  F. Simoes, February 2008                                                   !
!  Last updated (mm-dd-yyyy): 03-25-2011 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables:
  REAL (KIND=mp), INTENT(IN) :: time

! Local variables:
  INTEGER :: i,j,k
  REAL (KIND=mp) :: delta,delta2
  LOGICAL, EXTERNAL :: equals

! Interpolate for outflow boundary (stage).  Do this only for htype() = 1.
  DO k = 1,n_outflowbdr
    IF (htype(k) /= 1) CYCLE
    IF (n_timeserh(k) > 0) THEN
      IF (time < timeloch(1,k)) THEN
        hbc(k) = timeserh(1,k)
      ELSE IF (time >= timeloch(n_timeserh(k),k)) THEN
        hbc(k) = timeserh(n_timeserh(k),k)
      ELSE
        DO i = 1,n_timeserh(k) - 1
          ! Use linear interpolation for regular table time intervals.
          IF (time >= timeloch(i,k) .AND. time < timeloch(i+1,k)) THEN
            delta = time - timeloch(i,k)
            delta2 = timeloch(i+1,k) - timeloch(i,k)
            hbc(k) = timeserh(i,k) + delta*(timeserh(i+1,k) - &
                          timeserh(i,k))/delta2
            EXIT
          END IF
        END DO
      END IF
    END IF
  END DO

! Interpolate for inflow boundary (velocity or discharge).
  DO k = 1,n_inflowbdr
    IF (n_timeserq(k) > 0) THEN
      IF (time < timelocq(1,k)) THEN
        qin(k) = timeserq(1,k)
        IF (vtype(k) == 3) qin2(k) = timeserq2(1,k)
      ELSE IF (time >= timelocq(n_timeserq(k),k)) THEN
        qin(k) = timeserq(n_timeserq(k),k)
        IF (vtype(k) == 3) qin2(k) = timeserq2(n_timeserq(k),k)
      ELSE
        DO j = 1,n_timeserq(k) - 1
          ! Use linear interpolation for regular table time intervals.
          IF (time >= timelocq(j,k) .AND. time < timelocq(j+1,k)) THEN
            delta = time - timelocq(j,k)
            delta2 = timelocq(j+1,k) - timelocq(j,k)
            qin(k) = timeserq(j,k) + delta*(timeserq(j+1,k) - &
                     timeserq(j,k))/delta2
            IF (vtype(k) == 3) qin2(k) = timeserq2(j,k) + delta* &
                                    (timeserq2(j+1,k) - timeserq2(j,k))/delta2
            EXIT
          END IF
        END DO
      END IF
    END IF
  END DO

END SUBROUTINE get_time_bcs
