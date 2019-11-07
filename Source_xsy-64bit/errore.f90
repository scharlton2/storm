SUBROUTINE errore
  USE parameters
  USE geometry
  USE dep_vars
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Computes `a posteriori error estimates based on the ideas of Zienkiewicz   !
!  and Zhu (1987).  Call before the main solution output subroutines.         !
!                                                                             !
!  Francisco Simoes, July 2006                                                !
!  Last updated (mm-dd-yyyy): 09-04-2007 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

  ! Local variables:
  INTEGER :: i,ierror
  REAL (KIND=mp) :: l2h,l2u,l2v
  REAL (KIND=mp), ALLOCATABLE, DIMENSION(:,:) :: sigmah,sigmahs,sigmau, &
                                                 sigmaus,sigmav,sigmavs
  ! Allocate local variables. Not very efficient, but ok because this
  ! subroutine is only called a few times (i.e., once before each solution
  ! output request).
  ALLOCATE(sigmah(n_elems,2),sigmahs(n_elems,2),sigmau(n_elems,2), &
    sigmaus(n_elems,2),sigmav(n_elems,2),sigmavs(n_elems,2),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)

  ! Compute the finite element triangle gradients.
  CALL grad(u,n_pts,sigmau,n_elems)
  CALL grad(v,n_pts,sigmav,n_elems)
  CALL grad(h,n_pts,sigmah,n_elems)

  ! Compute the higher-order gradients.
  CALL gradhordr(sigmau,sigmaus)
  CALL gradhordr(sigmav,sigmavs)
  CALL gradhordr(sigmah,sigmahs)

!  DO i = 1,n_elems
!    WRITE (*,'(I6,6(ES19.10))') i,sigmau(i,1),sigmau(i,2),sigmav(i,1), &
!      sigmav(i,2),sigmah(i,1),sigmah(i,2)
!    WRITE (*,'(I6,6(ES19.10))') i,sigmaus(i,1),sigmaus(i,2),sigmavs(i,1), &
!      sigmavs(i,2),sigmahs(i,1),sigmahs(i,2)
!  END DO

  ! Use the L2 error norm of difference between the two gradients.
  DO i = 1,n_elems
    l2u = SQRT((sigmaus(i,1) - sigmau(i,1))**2 + (sigmaus(i,2) - sigmau(i,2))**2)
    l2v = SQRT((sigmavs(i,1) - sigmav(i,1))**2 + (sigmavs(i,2) - sigmav(i,2))**2)
    l2h = SQRT((sigmahs(i,1) - sigmah(i,1))**2 + (sigmahs(i,2) - sigmah(i,2))**2)
    e_est(i) = MAX(l2u,l2v,l2h) ! Uses euclidean norm (Maisano et al (2006)).
  END DO

  DEALLOCATE(sigmah,sigmahs,sigmau,sigmaus,sigmav,sigmavs)

END SUBROUTINE errore
