SUBROUTINE cde_t2t(t2tw,i,dim1,dim2)
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This subroutine is used in subroutine 'node_db' to clean-up duplicate      !
!  entries in the construction of arrays t2t and t2tHD.                       !
!                                                                             !
!  INPUT:                                                                     !
!    t2tw       working array containing the t2t connectivity table;          !
!    i          entry to be processed in table t2t;                           !
!    dim1,dim2  dimension of array t2tw.                                      !
!                                                                             !
!  OUTPUT:                                                                    !
!    t2tw       array without duplicate entries in row i.                     !
!                                                                             !
!  F. Simoes, September 2007                                                  !
!  Last updated (mm-dd-yyyy): 09-25-2007 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables.
  INTEGER, INTENT(IN) :: i,dim1,dim2
  INTEGER, INTENT(INOUT), DIMENSION(dim1,dim2) :: t2tw

! Local variables.
  INTEGER :: j,k,l,m,n,nn

  ! Clean-up duplicate entries.
  n = t2tw(i,1)
  l = 1
  DO
    j = t2tw(i,l+1)
    DO m = l+1,n
      k = t2tw(i,m+1)
      IF (j == k) THEN ! Duplicate entry found.
        n = n - 1
        DO nn = m,n
          t2tw(i,nn+1) = t2tw(i,nn+2)
        END DO
        t2tw(i,1) = n
      END IF
    END DO

    l = l + 1
    IF (l > n) EXIT
  END DO

END SUBROUTINE cde_t2t
