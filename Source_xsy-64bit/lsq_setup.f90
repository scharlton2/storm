SUBROUTINE lsq_setup(solver,type)
  USE parameters
  USE geometry
  USE constants
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Set-up the grid geometry quatities needed to use the least-squares         !
!  gradient construction technique of Mavriplis (2003).                       !
!                                                                             !
!  INPUT:                                                                     !
!    solver = 1 for RDS and other vertex-oriented solvers;                    !
!             2 for FVT and other cellcentered solvers;                       !
!    type   = 1 for unweighted least-squares construction;                    !
!           = 2 for inverse distance weighting.                               !
!                                                                             !
!  F. Simoes, January 2007                                                    !
!  Last updated (mm-dd-yyyy): 11-02-2007 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments:
  INTEGER, INTENT(IN) :: type,solver

! Local variables:
  INTEGER :: i,ierror,j,k,l

  SELECT CASE (solver)

!-----------------------------------------------------------------------------!
!                                 RDS solver                                  !
!-----------------------------------------------------------------------------!
  CASE (1)
    ! Allocate array memory.
    k = 0
    DO i = 1,n_pts
      k = MAX(k,n2n(i,1))
    END DO
    ALLOCATE(ai(n_pts),bi(n_pts),ci(n_pts),det(n_pts),dx(n_pts,k), &
      dy(n_pts,k),lsqweight(n_pts,k),STAT=ierror)
    IF (ierror /= 0) CALL alloc_err(ierror)
    CALL mem_add(8*(4 + 3*k)*n_pts)

    ! Compute mesh connectivity-related quantities.  Note that the dx and dy
    ! arrays are related to the n2n table.
    DO i = 1,n_pts
      j = n2n(i,1)
      DO k = 2,j+1
        l = n2n(i,k)
        dx(i,k-1) = nodes(l)%x - nodes(i)%x
        dy(i,k-1) = nodes(l)%y - nodes(i)%y
      END DO
    END DO

    ! Compute the weights.
    SELECT CASE (type)
    CASE (1)  ! Unweighted least squares.
      lsqweight = one

    CASE (2,3)        ! Weighted least squares.  Imparts better conditioning to
      DO i = 1,n_pts  ! the system of equations.
        j = n2n(i,1)
        DO k = 1,j
          ! Weight squared.
          lsqweight(i,k) = one/(dx(i,k)*dx(i,k) + dy(i,k)*dy(i,k))
        END DO
      END DO

    CASE DEFAULT
      PRINT *,"ERROR: least-squares procedure, invalid option."
      CALL byebye('Program STORM stopped.')

    END SELECT

    ! Calculate ai, bi, ci, and det used in Cramer's rule to compute node
    ! gradients.
    ai = zero; bi = zero; ci = zero
    DO i = 1,n_pts
      j = n2n(i,1)
      DO k = 1,j
        ai(i) = ai(i) + lsqweight(i,k)*dx(i,k)*dx(i,k)
        bi(i) = bi(i) + lsqweight(i,k)*dx(i,k)*dy(i,k)
        ci(i) = ci(i) + lsqweight(i,k)*dy(i,k)*dy(i,k)
      END DO
      det(i) = ai(i)*ci(i) - bi(i)*bi(i)
    END DO

!-----------------------------------------------------------------------------!
!                                 FVT solver                                  !
!-----------------------------------------------------------------------------!
  CASE (2)
    ! Allocate array memory.
    k = 0
    DO i = 1,n_elems
      k = MAX(k,t2t(i,1))
    END DO
    ALLOCATE(ai(n_elems),bi(n_elems),ci(n_elems),det(n_elems),dx(n_elems,k), &
      dy(n_elems,k),lsqweight(n_elems,k),STAT=ierror)
    IF (ierror /= 0) CALL alloc_err(ierror)
    CALL mem_add(8*(4 + 3*k)*n_elems)

    ! Compute mesh connectivity-related quantities.  Note that the dx and dy
    ! arrays are related to the t2t table.
    DO i = 1,n_elems
      j = t2t(i,1)
      DO k = 2,j+1
        l = t2t(i,k)
        dx(i,k-1) = grid(l)%xc - grid(i)%xc
        dy(i,k-1) = grid(l)%yc - grid(i)%yc
      END DO
    END DO

    ! Compute the weights.
    SELECT CASE (type)
    CASE (1)  ! Unweighted least squares.
      lsqweight = one

    CASE (2)            ! Weighted least squares.  Imparts better conditioning
      DO i = 1,n_elems  ! to the system of equations.
        j = t2t(i,1)
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
      j = t2t(i,1)
      DO k = 1,j
        ai(i) = ai(i) + lsqweight(i,k)*dx(i,k)*dx(i,k)
        bi(i) = bi(i) + lsqweight(i,k)*dx(i,k)*dy(i,k)
        ci(i) = ci(i) + lsqweight(i,k)*dy(i,k)*dy(i,k)
      END DO
      det(i) = ai(i)*ci(i) - bi(i)*bi(i)
    END DO

  CASE DEFAULT
    PRINT *,"ERROR: least-squares procedure, invalid option."
    CALL byebye('Program STORM stopped.')

  END SELECT

END SUBROUTINE lsq_setup
