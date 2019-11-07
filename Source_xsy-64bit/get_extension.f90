CHARACTER (LEN=3) FUNCTION get_extension(fname)
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  A function that returns the extension of a given file name.  The           !
!  extension is defined as being the string after the last dot.  The          !
!  returned extension is in upper case. The first three characters of the     !
!  extension are returned.  Example: get_extension("testfile.data") returns   !
!  "DAT"; get_extension("bull.1010.cgns") returns "CGN";                      !
!  get_extension("filena.me") returns "ME".                                   !
!                                                                             !
!  INPUT:                                                                     !
!    fname    filename with extension, in upper case.                         !
!                                                                             !
!  OUTPUT:                                                                    !
!                                                                             !
!  Francisco Simoes, April 2006                                               !
!  Last updated (mm-dd-yyyy): 06-21-2006 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables.
  CHARACTER (LEN=*), INTENT(IN) :: fname

! Local variables.
  INTEGER :: i,j
  CHARACTER, EXTERNAL :: to_upper

  get_extension = ''
  i = INDEX(fname,'.',.TRUE.)  ! Find last dot.
  IF (i == 0) RETURN  ! Extension not found.
  j = LEN_TRIM(fname)
  IF (i == j) RETURN  ! Extension not found after the dot.
  j = MIN(j,i+3)  ! The extension has 3 characters or less.
  get_extension = fname(i+1:j)

  ! Convert to upper case.
  DO i = 1,3
    get_extension(i:i) = to_upper(get_extension(i:i))
  END DO

END FUNCTION get_extension
