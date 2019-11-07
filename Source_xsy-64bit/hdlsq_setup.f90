SUBROUTINE hdlsq_setup(wtype)
  USE parameters
  USE geometry
  USE constants
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Set-up the grid geometry quatities needed to use the least-squares         !
!  gradient construction technique of Mavriplis (2003).  This only impacts    !
!  the FVT solver, because the RDS solver already uses the largest            !
!  computational molecule possible.                                           !
!                                                                             !
!  INPUT:                                                                     !
!    wtype  = 1 for unweighted least-squares construction;                    !
!           = 2 for inverse distance weighting.                               !
!                                                                             !
!  F. Simoes, September 2007                                                  !
!  Last updated (mm-dd-yyyy): 09-25-2007 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments:
  INTEGER, INTENT(IN) :: wtype

! Local variables:
  INTEGER :: i,ierror,j,k,l

! Allocate array memory.
  k = 0
  DO i = 1,n_elems
    k = MAX(k,t2tHD(i,1))
  END DO
  ALLOCATE(ai(n_elems),bi(n_elems),ci(n_elems),det(n_elems),dx(n_elems,k), &
    dy(n_elems,k),lsqweight(n_elems,k),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(8*(4 + 3*k)*n_elems)

! Compute mesh connectivity-related quantities.  Note that the dx and dy arrays
! are related to the t2tHD table.
  DO i = 1,n_elems
    j = t2tHD(i,1)
    DO k = 2,j+1
      l = t2tHD(i,k)
      dx(i,k-1) = grid(l)%xc - grid(i)%xc
      dy(i,k-1) = grid(l)%yc - grid(i)%yc
    END DO
  END DO

! Compute the weights.
  SELECT CASE (wtype)
  CASE (1)  ! Unweighted least squares.
    lsqweight = one

  CASE (2)            ! Weighted least squares.  Imparts better conditioning
    DO i = 1,n_elems  ! to the system of equations.
      j = t2tHD(i,1)
      DO k = 1,j
        ! Weight squared.
        lsqweight(i,k) = one/(dx(i,k)*dx(i,k) + dy(i,k)*dy(i,k))
      END DO
    END DO

  CASE DEFAULT
    PRINT *,"ERROR: least-squares procedure, invalid option."
    CALL byebye('Program STORM stopped.')

  END SELECT

! Calculate ai, bi, ci, and det used in Cramer's rule to compute cell
! gradients.
  ai = zero; bi = zero; ci = zero
  DO i = 1,n_elems
    j = t2tHD(i,1)
    DO k = 1,j
      ai(i) = ai(i) + lsqweight(i,k)*dx(i,k)*dx(i,k)
      bi(i) = bi(i) + lsqweight(i,k)*dx(i,k)*dy(i,k)
      ci(i) = ci(i) + lsqweight(i,k)*dy(i,k)*dy(i,k)
    END DO
    det(i) = ai(i)*ci(i) - bi(i)*bi(i)
  END DO

END SUBROUTINE hdlsq_setup
