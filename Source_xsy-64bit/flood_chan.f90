SUBROUTINE flood_chan(t)
  USE parameters
  USE constants
  USE geometry
  USE dep_vars
  USE memory
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This subroutine is where the spill from the one-dimensional channels to    !
!  the two-dimensional computational domain of SToRM.                         !
!                                                                             !
!  INPUT:                                                                     !
!    t     time, in seconds, corresponting to the current time step.          !
!                                                                             !
!  F. Simoes, September 2011                                                  !
!  Last updated (mm-dd-yyyy): 10-04-2011 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy subroutine arguments:
  REAL (KIND=mp) :: t  ! Time of the simulation, in seconds.

! Local variables.
  INTEGER :: i,j
  REAL (KIND=mp) :: clength,hsink,hwater,s,scut,s_i(2),t0,vwave

! Define basic parameters.
! Initial time, needed for SToRM restarts.  If this is a continuation run, t0
! should be the ending time of the previous run; otherwise set to zero.
  t0 = zero
! Velocity of the propagating wave.
  vwave = 1.0_mp  ! m/s
! Define inflow, up to where there's spill.  This is done in terms of the
! channel length parameter s.
  scut = 0.60_mp  ! 60%
  clength = 28841.599_mp  ! Channel length.
  s = vwave*(t + t0)/clength
  s = MIN(s,scut)
  IF (t > 86400) s = zero  ! Cut inflow off at 24 hours.
! Define the region of sink back into the 1D channel.
  s_i(1) = 0.90_mp; s_i(2) = one  ! 90% and 100%
  IF (t > 93600_mp) s_i(1) = zero  ! Start pumping out more water at hour 26.
! Amount of water to add to cell for each time step (depth, m).
  hwater = 0.10_mp
! Amount of water to remove from cell at each time step (% from depth).
  hsink = 0.20_mp

  s = zero  ! For 2 days after.
  s_i(1) = zero
  hsink = one

  DO j = 1,n_chanbdr
    DO i = 1,chantrigs(0,j)
      ! Add water to cell:
      IF (chant(i,j) < s) h(chantrigs(i,j)) = h(chantrigs(i,j)) + hwater
      ! Remove water from cell:
      IF (chant(i,j) > s_i(1) .AND. chant(i,j) < s_i(2)) &
        h(chantrigs(i,j)) = h(chantrigs(i,j))*(one - hsink)
    END DO
  END DO

END SUBROUTINE flood_chan
