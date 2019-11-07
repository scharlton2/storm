SUBROUTINE c2v_setup(type)
  USE parameters
  USE geometry
  USE constants
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Set-up interpolating weights in the cell-to-vertex least-squares           !
!  interpolation technique following section 3.3 of Coudiere et al (1999).    !
!  Note that this construction uses the augmented table n2t2.                 !
!                                                                             !
!  INPUT:                                                                     !
!    type      numerical method used in the interpolation: = 1 for Maisano    !
!              et al. (2006); = 2 for the least-squares technique of          !
!              Coudiere et al. (1999); = 3 for the inverse distance.  This    !
!              value must be the same as the equivalent parameter in the      !
!              call to subroutine 'ctr2vtx'.                                  !
!                                                                             !
!  F. Simoes, September 2007                                                  !
!  Last updated (mm-dd-yyyy): 08-22-2012 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  INTEGER, INTENT(IN) :: type

! Local variables:
  INTEGER :: i,ierror,j,k,l
  REAL (KIND=mp) :: d,ixx,ixy,iyy,lx,ly,rx,ry,tdist,xk,yk
  LOGICAL :: clipping

  k = SIZE(n2t2,2)
  ALLOCATE(c2v_weights(n_pts,k-1),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)
  CALL mem_add(8*n_pts*(k - 1))

  SELECT CASE (type)
  CASE (1)
    ! No precomputed coefficients are used: do nothing.

  CASE (2)
    ! Set-up interpolating weights in the cell-to-vertex least-squares
    ! interpolation technique following section 3.3 of Coudiere et al (1999).
    ! Note that this construction uses the augmented table n2t2.

    DO i = 1,n_pts
      rx = zero
      ry = zero
      ixx = zero
      ixy = zero
      iyy = zero
      k = n2t2(i,1)
      DO j = 2,k+1
        l = n2t2(i,j)
        xk = grid(l)%xc - nodes(i)%x
        yk = grid(l)%yc - nodes(i)%y
        rx = rx + xk
        ry = ry + yk
        ixx = ixx + xk*xk
        ixy = ixy + xk*yk
        iyy = iyy + yk*yk
      END DO
      d = ixx*iyy - ixy*ixy
      lx = (ixy*ry - iyy*rx)/d
      ly = (ixy*rx - ixx*ry)/d
      DO j = 2,k+1
        l = n2t2(i,j)
        xk = grid(l)%xc - nodes(i)%x
        yk = grid(l)%yc - nodes(i)%y
        c2v_weights(i,j-1) = (one + lx*xk + ly*yk)/(REAL(k,mp) + lx*rx + ly*ry)
      END DO

      ! If weight clipping is desired, uncomment the next lines.
      !clipping = .FALSE.
      !DO j = 2,k+1
      !  IF (c2v_weights(i,j-1) < zero) THEN  ! Clip below.
      !    c2v_weights(i,j-1) = zero
      !    clipping = .TRUE.
      !  END IF
      !  IF (c2v_weights(i,j-1) > 2.0_mp) THEN  ! Clip above.
      !    c2v_weights(i,j-1) = 2.0_mp
      !    clipping = .TRUE.
      !  END IF
      !END DO
      !IF (clipping) THEN  ! Normalize weights (not sure if this is needed).
      !  d = zero
      !  DO j = 2,k+1
      !    d = d + c2v_weights(i,j-1)
      !  END DO
      !  DO j = 2,k+1
      !    c2v_weights(i,j-1) = c2v_weights(i,j-1)/d
      !  END DO
      !END IF
      ! End of clipping procedure.
    END DO

  CASE (3)
    ! Set-up the inverse distance weighting coefficients -- see eqs. (2-3) of
    ! Maisano et al. (2006), for example.

    ALLOCATE(c2v_weights_id(n_pts,k-1),STAT=ierror)
    IF (ierror /= 0) CALL alloc_err(ierror)
    CALL mem_add(8*n_pts*(k - 1))

    DO i = 1,n_pts
      tdist = zero
      DO j = 2,n2t(i,1)+1
        k = n2t(i,j)
        xk = grid(k)%xc - nodes(i)%x
        yk = grid(k)%yc - nodes(i)%y
        c2v_weights(i,j-1) = one/SQRT(xk*xk + yk*yk)
        c2v_weights_id(i,j-1) = c2v_weights(i,j-1)
        tdist = tdist + c2v_weights(i,j-1)
      END DO
      DO j = 2,n2t(i,1)+1
        c2v_weights(i,j-1) = c2v_weights(i,j-1)/tdist
      END DO
    END DO

  CASE DEFAULT
    PRINT *,"ERROR: invalid option in c2v_setup."
    CALL byebye('Program STORM stopped.')

  END SELECT

END SUBROUTINE c2v_setup
