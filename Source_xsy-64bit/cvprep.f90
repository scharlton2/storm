SUBROUTINE cvprep
  USE geometry
  USE dep_vars
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This SUBROUTINE prepares some variables needed to perform the culvert      !
!  computations.                                                              !
!                                                                             !
!  F. Simoes, January 2013                                                    !
!  Last updated (mm-dd-yyyy): 04-19-2013 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables:
  INTEGER :: i,j
  LOGICAL, EXTERNAL :: pto_in_triangle

! Find the triangles where the inlet and outlet of each culvert are located.
  DO i = 1,nculvert
    ! Find the inlet triangle.
    cvtrigin(i) = 0
    DO j = 1,n_elems
      IF (pto_in_triangle(grid(j),cvptoin(1,i),cvptoin(2,i))) THEN
        cvtrigin(i) = j
        EXIT
      END IF
    END DO
    IF (cvtrigin(i) == 0) THEN  ! Error checking.
      WRITE(*,'(" ERROR: inlet of culvert ",I5, &
        " not in computational domain.")') i
      CALL byebye(' Program SToRM stopped.')
    END IF

    ! Compute the position vector of the culvert inlet in the triangle.
    cvrin(1,i) = cvptoin(1,i) - grid(j)%xc
    cvrin(2,i) = cvptoin(2,i) - grid(j)%yc

    ! Find the outlet triangle.
    cvtrigout(i) = 0
    DO j = 1,n_elems
      IF (pto_in_triangle(grid(j),cvptoout(1,i),cvptoout(2,i))) THEN
        cvtrigout(i) = j
        EXIT
      END IF
    END DO
    IF (cvtrigout(i) == 0) THEN  ! Error checking.
      WRITE(*,'(" ERROR: outlet of culvert ",I5, &
        " not in computational domain.")') i
      CALL byebye(' Program SToRM stopped.')
    END IF

  END DO

END SUBROUTINE cvprep
