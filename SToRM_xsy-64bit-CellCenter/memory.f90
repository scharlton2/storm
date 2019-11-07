!-----------------------------------------------------------------------------!
!                                                                             !
!  This file contains the physical memory size used by the main variables in  !
!  SToRM (i.e., those that are ALLOCATEd).                                    !
!                                                                             !
!  F. Simoes, July 2007                                                       !
!  Last updated (mm-dd-yyyy): 08-20-2007 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

MODULE memory
  IMPLICIT NONE
  SAVE

  INTEGER :: mem_used = 0  ! Number of bytes ALLOCATEd by SToRM.
  INTEGER, PARAMETER :: &
    e_size = 56, &  ! Size of TYPE(edge)    variables, in bytes.
    p_size = 16, &  !  "   "  TYPE(point)       "      "   "   .
    t_size = 60, &  !  "   "  TYPE(triangle)    "      "   "   .
    v_size = 16     !  "   "  TYPE(vector)      "      "   "   .

END MODULE memory

! Note about Fortran: a local variable that is initialized when it is declared
! has an implicit SAVE attribute.  This means that the variable 'mem_used' is
! initialized only the first time the MODULE 'memory' is used.  Subsequent uses
! of this MODULE do not reinitialize 'mem_used', and the value of its contents
! (i.e., of the contents of variable 'mem_used') is retained.  This is also
! true when initializing variables in their declaration statements in functions
! and subroutines: the variable is initialized on the first call to the
! function/subroutine, but on subsequent calls the old value  of the variable
! is retained.
