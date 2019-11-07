SUBROUTINE ifname(v,vmax,fn,fnout)
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This subroutine creates a file name from a root file name and a number.    !
!  The name generation is as follows: the root of the name is passed in       !
!  variable 'fn', to which the numerical value in variable 'v' is added,      !
!  followed by extension '.dat'.  Variable 'vmax' is used to add zeroes       !
!  between 'fn' and 'v' if 'v' has less digits than 'vmax'.  For example, if  !
!  fn = 'Agent', v = 300, and vmax = 500, then fnout = 'Agent300.dat'. If,    !
!  however, v = 7, then fnout = 'Agent007.dat'.                               !
!                                                                             !
!  INPUT:                                                                     !
!    v        value (string) to be added to the filename;                     !
!    vmax     value to be used for spacing purposes;                          !
!    fn       root filename.                                                  !
!                                                                             !
!  OUTPUT:                                                                    !
!    fnout    filename created by concatenating 'fn' with 'v', using 'vmax'   !
!             as a space template.  Returns '' if 'v' has more digits than    !
!             'vmax', which is an error condition.                            !
!                                                                             !
!  Francisco Simoes, November 2005                                            !
!  Last updated (mm-dd-yyyy): 11-28-2005 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables.
  INTEGER, INTENT(IN) :: v,vmax
  CHARACTER (LEN=*), INTENT(IN) :: fn
  CHARACTER (LEN=*), INTENT(OUT) :: fnout

! Local variables.
  INTEGER :: i,j,k
  CHARACTER (LEN=80) :: buffer

! Determine the size of vmax (i.e., how many characters in contains).
  WRITE (buffer,*) vmax
  buffer = ADJUSTL(buffer)
  i = LEN_TRIM(buffer)

! Determine the size of V.
  WRITE (buffer,*) v
  buffer = ADJUSTL(buffer)
  j = LEN_TRIM(buffer)

  IF (j > i) THEN  ! ERROR: v cannot be larger than vmax.
    fnout = ''
    RETURN
  END IF

  fnout = TRIM(ADJUSTL(fn))
  DO k = j,i - 1
    fnout = TRIM(fnout) // "0"
  END DO
  fnout = TRIM(fnout) // TRIM(buffer) // ".dat"

END SUBROUTINE ifname
