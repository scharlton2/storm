SUBROUTINE drynwet(e,zout)
  USE parameters
  USE dep_vars
  USE options
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  Applies the Brufau and Garcia-Navarro (2003) algorithm--sometimes called   !
!  the BGN algorithm in this code--to those triangles that fall within the    !
!  criterion.  It takes element 'e' as input and outputs the modified bed     !
!  elevations in array zout().                                                !
!                                                                             !
!  Francisco Simoes, October 2005                                             !
!  Last updated (mm-dd-yyyy): 10-31-2005 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy variables.
  TYPE(triangle), INTENT(IN) :: e
  REAL(KIND=mp), INTENT(OUT) :: zout(3)

! Local variables.
  INTEGER :: n,i,j,k
  REAL(KIND=mp) :: dj,dk,dmax
  LOGICAL, EXTERNAL :: part_dry

  zout(1) = z(e%vertex(1))
  zout(2) = z(e%vertex(2))
  zout(3) = z(e%vertex(3))

  IF (.NOT.zadjust) RETURN  ! Do not use the BGN algorithm.

  IF (.NOT.part_dry(e,n,i,j)) RETURN  ! None of the nodes is dry.

! Apply the BGN algorithm to the nodes that lye abobe the maximum stage in the
! element.
  IF (n == 1) THEN  ! Only one node is dry.
    j = MAX(1,MOD(i+1,4))
    k = MAX(1,MOD(i+2,4))
    dj = z(e%vertex(j)) + h(e%vertex(j))
    dk = z(e%vertex(k)) + h(e%vertex(k))
    dmax = MAX(dj,dk)
    zout(i) = MIN(dmax,z(e%vertex(i)))  ! Adjust the dry node.

  ELSE IF (n == 2) THEN  ! Two dry nodes.
    k = 6 - i - j
    dk = z(e%vertex(k)) + h(e%vertex(k))
    zout(i) = MIN(dk,z(e%vertex(i)))  ! Adjust the dry nodes.
    zout(j) = MIN(dk,z(e%vertex(j)))

  END IF

END SUBROUTINE drynwet
