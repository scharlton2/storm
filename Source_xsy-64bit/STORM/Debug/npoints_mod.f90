        !COMPILER-GENERATED INTERFACE MODULE: Fri Jan 23 13:46:48 2009
        MODULE NPOINTS_mod
          INTERFACE 
            FUNCTION NPOINTS(VARNO,DATAV,N_PTS,N_ELEMS) RESULT(NPOINTS_0&
     &)
              INTEGER(KIND=4), INTENT(IN) :: VARNO
              LOGICAL(KIND=4), INTENT(IN) :: DATAV(:)
              INTEGER(KIND=4), INTENT(IN) :: N_PTS
              INTEGER(KIND=4), INTENT(IN) :: N_ELEMS
              INTEGER(KIND=4) :: NPOINTS_0
            END FUNCTION NPOINTS
          END INTERFACE 
        END MODULE NPOINTS_mod
