FUNCTION vbc_by_h(ibdr,wdepth,ndim)
  USE parameters
  USE constants
  USE dep_vars
  USE options
  USE vbc_arrays
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This subroutine computes the velocity distribution along the inflow        !
!  boundary, distributing the discharge as a function of the local water      !
!  depth.  The discharge computed using the thus calculated velocity distri-  !
!  bution is given by q.  All tests using this subroutine concluded that      !
!  q = qin always.                                                            !
!                                                                             !
!  INPUT:                                                                     !
!    ibdr      inflow boundary number belonging to [1,n_inflowbdr];           !
!    wdepth    main array containing the water depth at the vertices of the   !
!              triangles;                                                     !
!    ndim      size of array wdepth.  It's used here just because the bounds  !
!              checking option of many compilers do not work for arrays of    !
!              unknown length.                                                !
!                                                                             !
!  OUTPUT:                                                                    !
!    vbc_by_h  discharge at inflow cross section calculated from the          !
!              velocity vectors computed in this subroutine.  This value is   !
!              passed out mainly for debugging purposes.  It should come out  !
!              equal to variable 'qin(ibdr)' in MODULE 'dep_vars'.            !
!                                                                             !
!  Francisco Simoes, November 2005                                            !
!  Last updated (mm-dd-yyyy): 02-25-2011 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables.
  INTEGER, INTENT(IN) :: ibdr,ndim
  REAL (KIND=mp), INTENT(IN), DIMENSION(ndim) :: wdepth
  REAL (KIND=mp) :: vbc_by_h

! Local variables.
  INTEGER :: i
  REAL (KIND=mp) :: at,ht,q

! Initialize needed arrays.
  DO i = 1,n_qin(ibdr)
    h_inflow(i) = wdepth(qin_nodes(i,ibdr))
  END DO
  ht = zero
  at = zero
  DO i = 1,n_qin(ibdr)-1
    a_inflow(i) = l_inflow(i,ibdr)*(h_inflow(i) + h_inflow(i+1))*half
    ht = ht + (h_inflow(i) + h_inflow(i+1))*half
    at = at + a_inflow(i)
  END DO

! Compute velocity at mid-points of faces.
  DO i = 1,n_qin(ibdr)-1
    u_inflow(i) = half*(h_inflow(i) + h_inflow(i+1))*qin(ibdr)/ht/at
  END DO

! First pass to compute v_inflow(i).
  v_inflow = zero
  DO i = 2,n_qin(ibdr)-1
    IF (h_inflow(i) > h_dry) v_inflow(i) = half*(u_inflow(i-1) + u_inflow(i))
  END DO

! Now compute the new u_inflow(i).
  DO i = 1,n_qin(ibdr)-1
    u_inflow(i) = half*(v_inflow(i) + v_inflow(i+1))
  END DO

! Compute the value of the new (inaccurate) discharge.
  q = zero
  DO i = 1,n_qin(ibdr)-1
    q = q + u_inflow(i)*a_inflow(i)
  END DO

! Second pass: adjust v_inflow(i) according to the difference between q and
! qin().
  !v_inflow(1) = zero  ! No need: this is already enforced in the code above.
  !v_inflow(n_qin(ibdr)) = zero
  DO i = 2,n_qin(ibdr)-1
    IF (h_inflow(i) > h_dry) v_inflow(i) = v_inflow(i)*qin(ibdr)/q
  END DO

! Compute the adjusted discharge (for debugging purposes only).
  DO i = 1,n_qin(ibdr)-1
    u_inflow(i) = half*(v_inflow(i) + v_inflow(i+1))
  END DO
  q = zero
  DO i = 1,n_qin(ibdr)-1
    q = q + u_inflow(i)*a_inflow(i)
  END DO

  vbc_by_h = q

END FUNCTION vbc_by_h
