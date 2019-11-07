SUBROUTINE ctr2vtxWSE(type,carray,varray)
  USE parameters
  USE constants
  USE geometry
  USE dep_vars
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This subroutine interpolates the values from the cell centers to the       !
!  vertices of a triangular grid.  Yhis particular version of the ctr2vtx     !
!  SUBROUTINE was designed specifically to interpolate the water surface      !
!  elevation: it does not use values in dry cells.  Logical variable dycell   !
!  is used to check for dry cells.                                            !
!                                                                             !
!  INPUT:                                                                     !
!    type      numerical method used in the interpolation: = 1 for Maisano    !
!              et al. (2006); = 2 for the least-squares technique of          !
!              Coudiere et al. (1999); = 3 for the inverse distance.  This    !
!              value must be the same used in the call to subroutine          !
!              'c2v_setup';                                                   !
!    carray    array with the quantity of interest located at the geometric   !
!              centers of the triangles (dimension = n_elems).                !
!                                                                             !
!  OUTPUT:                                                                    !
!    varray    array with the same quantity interpolated to the vertices of   !
!              the triangles (dimension = n_pts).                             !
!                                                                             !
!  NOTE: when the least squares interpolation technique is used, option       !
!  type = 2, the desired benefits of this SUBROUTINE are not realized.  In    !
!  other words, dry cells are used for the interpolation and spureous wet     !
!  vertices my appear in holes and other concave regions.                     !
!                                                                             !
!  Francisco Simoes, February 2007                                            !
!  Last updated (mm-dd-yyyy): 08-22-2012 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy arguments.
  REAL (KIND=mp), INTENT(IN) :: carray(n_elems)
  REAL (KIND=mp), INTENT(OUT) :: varray(n_pts)
  INTEGER, INTENT(IN) :: type

! Local variables.
  INTEGER :: counter,i,j,k
  REAL (KIND=mp) :: a,id,tarea,weight

  SELECT CASE (type)
  CASE (1)
    ! The interpolation follows very closely the technique used to interpolate
    ! the gradients -- eq. (2) of Maisano et al. (2006) -- and uses the areas
    ! of neighboring triangles as weights.
    DO i = 1,n_pts
      tarea = zero
      a = zero
      counter = 0
      DO j = 2,n2t(i,1)+1
        k = n2t(i,j)
        IF (drycell(k)) CYCLE
        counter = counter + 1
        weight = grid(k)%area
        tarea = tarea + weight
        a = a + weight*carray(k)
      END DO
      IF (counter > 0) THEN
        varray(i) = a/tarea
      ELSE
        varray(i) = zvtx(i)  ! Hardcoded for water surface elevation.
      END IF
    END DO

  CASE (2)
    ! The second-order accurate least-squares technique in eq. (21) of Coudiere
    ! et al. (1999).  The weights are computed at the data pre-processing stage
    ! to save CPU time.
    DO i = 1,n_pts
      a = zero
      DO j = 2,n2t2(i,1)+1
        k = n2t2(i,j)
        a = a + c2v_weights(i,j-1)*carray(k)
      END DO
      varray(i) = a
    END DO

  CASE (3)
    ! This option uses inverse distance wighting coefficients -- see eqs. (2-3)
    ! of Maisano et al. (2006), for example.  The weights are computed at the
    ! data pre-processing stage to save CPU time.
    DO i = 1,n_pts
      id = zero
      a = zero
      counter = zero
      DO j = 2,n2t(i,1)+1
        k = n2t(i,j)
        IF (drycell(k)) CYCLE
        counter = counter + 1
        id = id + c2v_weights_id(i,j-1)
        a = a + c2v_weights_id(i,j-1)*carray(k)
      END DO
      IF (counter > 0) THEN
        varray(i) = a /id
      ELSE
        varray(i) = zvtx(i)  ! Hardcoded for water surface elevation.
      END IF
    END DO

  CASE DEFAULT
    PRINT *,"ERROR: invalid option in ctr2vtx."
    CALL byebye('Program STORM stopped.')

  END SELECT

END SUBROUTINE ctr2vtxWSE
