SUBROUTINE grad(phi,nvdim,gradphi,nedim)
  USE parameters
  USE geometry
  USE constants
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Use linear shape functions to compute the gradient of variable 'phi' in    !
!  the triangles of the computational domain.  'Phi' is any variable defined  !
!  at the nodes (u, v, h, etc.).  'Gradphi' contains the components of grad   !
!  'phi' and is defined at the centroids of the triangular elements.          !
!                                                                             !
!  INPUT:                                                                     !
!    nvdim    number of vertices in the mesh (must be the size of array       !
!             phi);                                                           !
!    phi      array with the variable to compute the gradient of.  This       !
!             variable must be defined at the vertices of triangles;          !
!    nedim    number of elements (triangles) in the grid.                     !
!                                                                             !
!  OUTPUT:                                                                    !
!    gradphi  an array of dimension 'nedim' containing the two components     !
!             (in the x and y directions) of the gradient.                    !
!                                                                             !
!  Francisco Simoes, July 2006                                                !
!  Last updated (mm-dd-yyyy): 09-04-2007 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables:
  INTEGER, INTENT(IN) :: nedim,nvdim
  REAL (KIND=mp), INTENT(IN) :: phi(nvdim)
  REAL (KIND=mp), INTENT(OUT) :: gradphi(nedim,2)

! Local variables:
  INTEGER :: i,n1,n2,n3

  DO i = 1,nedim
    n1 = grid(i)%vertex(1)
    n2 = grid(i)%vertex(2)
    n3 = grid(i)%vertex(3)
    gradphi(i,1) = (phi(n1)*(nodes(n2)%y - nodes(n3)%y) + &  ! X-component.
                    phi(n2)*(nodes(n3)%y - nodes(n1)%y) + &
                    phi(n3)*(nodes(n1)%y - nodes(n2)%y))*half/grid(i)%area
    gradphi(i,2) = (phi(n1)*(nodes(n3)%x - nodes(n2)%x) + &  ! Y-component.
                    phi(n2)*(nodes(n1)%x - nodes(n3)%x) + &
                    phi(n3)*(nodes(n2)%x - nodes(n1)%x))*half/grid(i)%area
  END DO

END SUBROUTINE grad
