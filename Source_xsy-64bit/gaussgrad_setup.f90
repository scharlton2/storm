SUBROUTINE gaussgrad_setup
  USE parameters
  USE constants
  USE geometry
  USE memory
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This subroutine computes all the coefficients necessary to compute the     !
!  dependent variables' cell gradients using the method described in          !
!  Anastasiou and Chan (1997) for the inviscid terms, eqs. (11-13).  This     !
!  subroutine must be CALLed before subroutine 'gaussgrad()' is ever used.    !
!                                                                             !
!  Francisco Simoes, February 2009                                            !
!  Last updated (mm-dd-yyyy): 02-28-2011 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables.
  INTEGER :: i,ierror,j,k,l
  REAL (KIND=mp) :: area,denom,x1,x2,xpto(3),y1,y2,ypto(3)
  TYPE(vector) :: n(3)

  ALLOCATE(gconn(n_elems,3),gcoefx(n_elems,3),gcoefy(n_elems,3), &
           gweight(n_elems,3),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(4*n_elems*3 + 8*n_elems*9)
  ALLOCATE(gwgrad(n_elems),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(v_size*n_elems)

  DO i = 1,n_elems
    DO j = 1,3
      k = t2t3(i,j)  ! Triangle adjacent to local edge j.
      IF (k > 0) THEN
        ! Normal triangle.
        xpto(j) = grid(i)%xc
        ypto(j) = grid(i)%yc
        gconn(i,j) = k
      ELSE
        ! Wall: use a ghost cell to mirror the contents of the currente cell.
        l = grid(i)%edge(j)
        x1 = nodes(edges(l)%p(1))%x
        y1 = nodes(edges(l)%p(1))%y
        x2 = nodes(edges(l)%p(2))%x
        y2 = nodes(edges(l)%p(2))%y
        CALL ghost_cell(grid(i)%xc,grid(i)%yc,x1,y1,x2,y2,xpto(j),ypto(j))
        gconn(i,j) = i
      END IF
    END DO

    ! Area is positive if the vertices are oriented in a clockwise manner,
    ! negative otherwise.
    area = half*((xpto(2) - xpto(1))*(ypto(3) - ypto(1)) - &
                 (xpto(3) - xpto(1))*(ypto(2) - ypto(1)))
    IF (area < zero) &
      CALL byebye("Error is Gauss gradient set-up (gaussgrad_setup).")

    n(1)%x = ypto(1) - ypto(2)  ;  n(1)%y = xpto(2) - xpto(1)
    n(2)%x = ypto(2) - ypto(3)  ;  n(2)%y = xpto(3) - xpto(2)
    n(3)%x = ypto(3) - ypto(1)  ;  n(3)%y = xpto(1) - xpto(3)

    ! Coefficients to compute the primary cell gradient.
    gcoefx(i,1) = half*(n(1)%x + n(3)%x)/area
    gcoefx(i,2) = half*(n(1)%x + n(2)%x)/area
    gcoefx(i,3) = half*(n(2)%x + n(3)%x)/area
    gcoefy(i,1) = half*(n(1)%y + n(3)%y)/area
    gcoefy(i,2) = half*(n(1)%y + n(2)%y)/area
    gcoefy(i,3) = half*(n(2)%y + n(3)%y)/area

    ! Weights for the final gradient computation.
    denom = (xpto(2) - xpto(1))*(ypto(3) - ypto(2)) - &
            (xpto(3) - xpto(2))*(ypto(2) - ypto(1))
    gweight(i,1) = ((xpto(2) - grid(i)%xc)*(ypto(3) - ypto(2)) - &
                    (xpto(3) - xpto(2))*(ypto(2) - grid(i)%yc))/denom
    gweight(i,2) = ((xpto(3) - grid(i)%xc)*(ypto(1) - ypto(3)) - &
                    (xpto(1) - xpto(3))*(ypto(3) - grid(i)%yc))/denom
    gweight(i,3) = ((xpto(1) - grid(i)%xc)*(ypto(2) - ypto(1)) - &
                    (xpto(2) - xpto(1))*(ypto(1) - grid(i)%yc))/denom
  END DO

END SUBROUTINE gaussgrad_setup
