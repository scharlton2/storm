        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 19 14:13:01 2018
        MODULE SOLVERRDS__genmod
          INTERFACE 
            SUBROUTINE SOLVERRDS(ITER,RMAX,FN,TITLE)
              INTEGER(KIND=4), INTENT(OUT) :: ITER
              REAL(KIND=8), INTENT(OUT) :: RMAX
              CHARACTER(*), INTENT(IN) :: FN
              CHARACTER(*), INTENT(IN) :: TITLE
            END SUBROUTINE SOLVERRDS
          END INTERFACE 
        END MODULE SOLVERRDS__genmod
