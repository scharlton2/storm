SUBROUTINE wetdryfix(itype,hv,hc,gradhc)
  USE parameters
  USE constants
  USE geometry
  USE options
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This subroutine applies a fix to the water depth computed at the vertices  !
!  of triangles for the vertices that are located at wet-dry boundaries.      !
!  This subroutine should be called after a call to subroutine 'ctr2vtx'.     !
!                                                                             !
!  INPUT:                                                                     !
!    itype     numerical method used in the interpolation from cell centers   !
!              to the vertices, as used in subroutine 'ctr2vtx';              !
!    hv        water depth at the vertices of the triangular grid;            !
!    hc        water depth at the centers of the triangules;                  !
!    gradhc    gradient of the water depth for each triangle.                 !
!                                                                             !
!  OUTPUT:                                                                    !
!    hv        corrected values of the water depth.                           !
!                                                                             !
!  NOTE: after further research, it was verified that a call to this          !
!  subroutine with itype = 2 can cause problems when 'n2t2' has wet and dry   !
!  nodes, but 'n2t' only has dry nodes (can happen in corners, for example).  !
!  In this case, counter is equal to 0 and a division by zero occurs.         !
!  Therefore, the recommendation is to use itype = 1 for all cases.           !
!                                                                             !
!  Francisco Simoes, January 2008                                             !
!  Last updated (mm-dd-yyyy): 01-25-2008 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  REAL (KIND=mp), INTENT(IN) :: hc(n_elems)
  REAL (KIND=mp), INTENT(INOUT) :: hv(n_pts)
  TYPE(vector), INTENT(IN) :: gradhc(n_elems)
  INTEGER, INTENT(IN) :: itype

! Local variables.
  INTEGER :: counter,i,j,k
  REAL (KIND=mp) :: newh
  LOGICAL :: dryt,wett

  SELECT CASE (itype)
  CASE (1)

    DO i = 1,n_pts
      wett = .FALSE.
      dryt = .FALSE.
      ! First determine if the node is at the boundary between wet and dry
      ! triangles.
      DO j = 2,n2t(i,1) + 1
        k = n2t(i,j)
        IF (drycell(k)) THEN
          dryt = .TRUE.  ! Dry triangle.
        ELSE
          wett = .TRUE.  ! Wet triangle.
        END IF
      END DO

      IF (wett .AND. dryt) THEN  ! Yes, vertex is at a wet-dry boundary.
        counter = 0  ! Counts the number of wet nodes.
        newh = zero
        DO j = 2,n2t(i,1) + 1
          k = n2t(i,j)
          IF (drycell(k)) CYCLE  ! Dry triangle.
          counter = counter + 1
          newh = newh + hc(k) + (nodes(i)%x - grid(k)%xc)*gradhc(k)%x + &
                                (nodes(i)%y - grid(k)%yc)*gradhc(k)%y
        END DO
        hv(i) = newh/REAL(counter,mp)
      END IF

    END DO

  CASE (2)

    DO i = 1,n_pts
      wett = .FALSE.
      dryt = .FALSE.
      ! First determine if the node is at the boundary between wet and dry
      ! triangles.
      DO j = 2,n2t2(i,1) + 1
        k = n2t2(i,j)
        IF (drycell(k)) THEN
          dryt = .TRUE.  ! Dry triangle.
        ELSE
          wett = .TRUE.  ! Wet triangle.
        END IF
      END DO

      IF (wett .AND. dryt) THEN  ! Yes, vertex is at a wet-dry boundary.
        counter = 0  ! Counts the number of wet nodes.
        newh = zero
        DO j = 2,n2t(i,1) + 1
          k = n2t(i,j)
          IF (drycell(k)) CYCLE  ! Dry triangle.
          counter = counter + 1
          newh = newh + hc(k) + (nodes(i)%x - grid(k)%xc)*gradhc(k)%x + &
                                (nodes(i)%y - grid(k)%yc)*gradhc(k)%y
        END DO
        hv(i) = newh/MAX(one,REAL(counter,mp))  ! Avoids division by zero.
      END IF

    END DO

  CASE DEFAULT
    PRINT *,"ERROR: invalid option in wetdryfix."
    CALL byebye('Program STORM stopped.')

  END SELECT

END SUBROUTINE wetdryfix
