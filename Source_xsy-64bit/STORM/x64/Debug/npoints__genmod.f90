        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 19 14:12:04 2018
        MODULE NPOINTS__genmod
          INTERFACE 
            FUNCTION NPOINTS(VARNO,DATAV,N_PTS,N_ELEMS)
              INTEGER(KIND=4), INTENT(IN) :: VARNO
              LOGICAL(KIND=4), INTENT(IN) :: DATAV(:)
              INTEGER(KIND=4), INTENT(IN) :: N_PTS
              INTEGER(KIND=4), INTENT(IN) :: N_ELEMS
              INTEGER(KIND=4) :: NPOINTS
            END FUNCTION NPOINTS
          END INTERFACE 
        END MODULE NPOINTS__genmod
