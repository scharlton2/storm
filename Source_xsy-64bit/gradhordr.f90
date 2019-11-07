SUBROUTINE gradhordr(fegrad,gradphi)
  USE parameters
  USE geometry
  USE constants
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Computes the gradient of function phi using higher order interpolants.     !
!                                                                             !
!  INPUT:                                                                     !
!    fegrad    finite element gradient of phi, defined at the element         !
!              centroids (computed, for example, with subroutine grad()).     !
!                                                                             !
!  OUTPUT:                                                                    !
!    gradphi   gradient of phi, defined at the centroid of the elements.      !
!                                                                             !
!  Francisco Simoes, July 2006                                                !
!  Last updated (mm-dd-yyyy): 03-22-2007 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables:
  REAL (KIND=mp), INTENT(IN) :: fegrad(n_elems,2)
  REAL (KIND=mp), INTENT(OUT) :: gradphi(n_elems,2)

! Local variables:
  INTEGER :: i,ierror,j,k,n1,n2,n3
  REAL (KIND=mp) :: tarea,weight
  REAL (KIND=mp), ALLOCATABLE, DIMENSION(:) :: phix,phiy

  ! Allocate local variables. Not very efficient, but ok because this
  ! subroutine is only called a few times (i.e., three times before each
  ! solution output request).
  ALLOCATE(phix(n_pts),phiy(n_pts),STAT=ierror)
  IF (ierror /= 0) CALL alloc_err(ierror)

  ! Compute the gradient at each node using eq. (2) of Maisano et al (2006) and
  ! the areas of the neighboring triangles as weights.
  DO i = 1,n_pts
    tarea = zero
    phix(i) = zero
    phiy(i) = zero
    DO j = 2,n2t(i,1) + 1
      k = n2t(i,j)
      weight = one/grid(k)%area
      tarea = tarea + weight
      phix(i) = phix(i) + weight*fegrad(k,1)
      phiy(i) = phiy(i) + weight*fegrad(k,2)
    END DO
    phix(i) = phix(i)/tarea
    phiy(i) = phiy(i)/tarea
  END DO

  ! Now take the average to find the gradient at the centroid of the triangle.
  DO i = 1,n_elems
    n1 = grid(i)%vertex(1)
    n2 = grid(i)%vertex(2)
    n3 = grid(i)%vertex(3)
    gradphi(i,1) = (phix(n1) + phix(n2) + phix(n3))/3.0_mp
    gradphi(i,2) = (phiy(n1) + phiy(n2) + phiy(n3))/3.0_mp
  END DO

  DEALLOCATE(phix,phiy)

END SUBROUTINE gradhordr

! While the finite element gradient is only first-order, I have oberved this
! technique to be second-order using the analytic solutions of Simoes (2005).
! See program gaccuracy and data in aanalysis.txt.
