        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 07 13:43:35 2014
        MODULE SOLVER__genmod
          INTERFACE 
            SUBROUTINE SOLVER(ITER,TTIME,RMAX,FN,TITLE)
              INTEGER(KIND=4), INTENT(OUT) :: ITER
              REAL(KIND=8), INTENT(OUT) :: TTIME
              REAL(KIND=8), INTENT(OUT) :: RMAX
              CHARACTER(*), INTENT(IN) :: FN
              CHARACTER(*), INTENT(IN) :: TITLE
            END SUBROUTINE SOLVER
          END INTERFACE 
        END MODULE SOLVER__genmod
