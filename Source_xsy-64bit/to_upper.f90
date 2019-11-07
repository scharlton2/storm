CHARACTER FUNCTION to_upper(ch)
  IMPLICIT NONE

! Converts a single character from lower case to upper case.  Leaves character
! unchanged if the character is not a lower case letter.
! F. Simoes, Dec. 2001

  ! Dummy argument.
  CHARACTER, INTENT(IN) :: ch

  ! Local variables.
  INTEGER, PARAMETER :: strand = IACHAR("A") - IACHAR("a")

  IF ("a" <= ch .AND. ch <= "z") THEN
    ! Lower case: change to upper case.
    to_upper = ACHAR(IACHAR(ch) + strand)
  ELSE
    ! Do nothing.
    to_upper = ch
  END IF

END FUNCTION to_upper
