INTEGER FUNCTION distribution()
  USE parameters
  USE geometry
  USE dep_vars
  USE constants
  USE options
  !USE io  ! Needed only for debugging, delete after removing all WRITE statms.
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Applies the selected residual distribution scheme to the cell residual to  !
!  come up with the nodal residual, PHI, used in the equations to update the  !
!  solution variables.  Upon return, it contains the number of failed         !
!  residual distribution schemes.                                             !
!                                                                             !
!  Francisco Simoes, October 2004                                             !
!  Last updated (mm-dd-yyyy): 11-16-2005 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Local variables:
  INTEGER :: counter,i,j,k,m,n
  REAL(KIND=mp) :: hcell,nxi,nyi,ucell,vcell
  REAL(KIND=mp) :: kpls(3,3),kmns(3,3),kplus(3,3,3),kminus(3,3,3),r(3,3)
  LOGICAL :: flag
  LOGICAL, EXTERNAL :: lda,narrow,nonlinearb,psi

! Reset nodal quantities PHI to zero.
  phi = zero  ! Note that this is an implicit DO-loop.

! Counter counts the number of failed residual distribution schemes.
  counter = 0

! Main DO-loop over all the elements.
  DO j = 1,n_elems

    ! Cell quantities, linearized as required by the residual distribution
    ! theory.
    hcell = u_avg(j,1)
    ucell = u_avg(j,2)
    vcell = u_avg(j,3)

    IF (hcell < h_dry) CYCLE  ! Dry cell.
    IF ((ucell*ucell + vcell*vcell) < u_stag*u_stag) THEN  ! Stagnant cell.
      CALL centered(j,r)

    ELSE
      ! Compute the kplus and kminus matrices for all of the element's nodes.
      DO i = 1,3
        k = MAX(1,MOD(i+1,4))  ! Element index of edge opposite to node i.
        nxi = edges(ABS(grid(j)%edge(k)))%normal(1)
        nyi = edges(ABS(grid(j)%edge(k)))%normal(2)
        ! We want the normal pointing inwards, therefore a sign change may be
        ! needed:
        IF (grid(j)%edge(k) > 0) THEN
          nxi = -nxi
          nyi = -nyi
        END IF

        CALL kmatrices(kpls,kmns,nxi,nyi,hcell,ucell,vcell)

        DO m = 1,3
          DO n = 1,3
            kplus(i,m,n) = kpls(m,n)
            kminus(i,m,n) = kmns(m,n)
          END DO
        END DO

      END DO

      ! Compute the nodal residuals according to the chosen residual
      ! distribution scheme.
      SELECT CASE (opt_rdscheme)

      CASE (0)  ! The N scheme (narrow).  This is the default option.
        flag = narrow(kplus,kminus,j,r)
        IF (.NOT. flag) THEN  ! N scheme failure.
          flag = lda(kplus,j,r)
          IF (.NOT. flag) counter = counter + 1  ! Both schemes failed
          ! Note: the computations can proceed normally if the resudual schemes
          ! fail, because in that case the residual contriburions to the nodes,
          ! r(.,.), are simply zero and can be added to phi(.,.) without ill
          ! effect.
        END IF

      CASE (1)  ! The LDA scheme (low diffusion A).
        flag = lda(kplus,j,r)
        IF (.NOT. flag) THEN  ! LDA scheme failure.
          flag = narrow(kplus,kminus,j,r)
          IF (.NOT. flag) counter = counter + 1  ! Both schemes failed.
        END IF

      CASE (2)  ! The B scheme.
        flag = nonlinearb(kplus,kminus,j,r)
        IF (.NOT. flag) counter = counter + 1

      CASE (3)  ! The PSI scheme.
        PRINT *,''
        PRINT *,'ERROR: the PSI scheme not yet implemented.'
        CALL byebye('Program SToRM stopped.')

      CASE DEFAULT
        PRINT *,''
        PRINT *,'ERROR: invalid value in opt_rdscheme.'
        CALL byebye('Program SToRM stopped.')
      END SELECT

    END IF

    ! r(m,n) is the element's residual contribution to node m (m is the local
    ! node indicator, m = 1,2,3). Phi(i,n) computed in this way is done for
    ! each node and has the contribution of all the adjacent triangles (i is
    ! the global node number).
    !WRITE (*,*) j
    !DO m = 1,3
    !  WRITE (*,'(3(2X,ES12.5))') (r(m,n),n=1,3)
    !END DO
    !WRITE (*,*)
    DO m = 1,3
      DO n = 1,3  ! h for n = 1, u for n = 2, and v for n = 3.
        phi(grid(j)%vertex(m),n) = phi(grid(j)%vertex(m),n) + r(m,n)
      END DO
    END DO

    ! Note about the debugging lines that follow: if the residual distribution
    ! scheme is properly implemented we should have
    !
    !                 r(1,m) + r(2,m) + r(3,m) = residual(j,m)
    !
    ! because of the conservative nature of the schemes.
    !WRITE (o_unit,'(/A)')'Nodal residual contribution:'
    !DO m = 1,3
    !  WRITE (o_unit,'(5(1X,ES12.5))') r(1,m),r(2,m),r(3,m), &
    !    r(1,m)+r(2,m)+r(3,m),residual(j,m)
    !END DO
    !WRITE (o_unit,*)

  END DO  ! End of sweep over all the elements.

  distribution = counter

END FUNCTION distribution
