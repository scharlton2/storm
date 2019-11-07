FUNCTION edge_in_array(edg,edges_array,n)
  USE parameters
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This subroutine returns .TRUE. if variable edg is in array edges_array,    !
!  and returns .FALSE. if it is not.  The order of the end points in edg      !
!  does not need to be the same as the same edge in array edges_array.  n is  !
!  the number of edges in edges_array (n >= 0).                               !
!                                                                             !
!  Francisco Simoes, March 2004                                               !
!  Last updated (mm-dd-yyyy): 03-11-2004                                      !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments:
  LOGICAL :: edge_in_array
  INTEGER, INTENT(IN) :: n
  TYPE(edge), INTENT(IN) :: edg,edges_array(n+1)

! Local variables:
  INTEGER :: i

  edge_in_array = .FALSE.
  IF (n > 0) THEN
    DO i = 1,n
      IF ((edg%p(1) == edges_array(i)%p(1)  .AND. &
           edg%p(2) == edges_array(i)%p(2)) .OR.  &
          (edg%p(2) == edges_array(i)%p(1)  .AND. &
           edg%p(1) == edges_array(i)%p(2))) THEN
        edge_in_array = .TRUE.
        RETURN
      END IF
    END DO
  END IF

END FUNCTION edge_in_array
