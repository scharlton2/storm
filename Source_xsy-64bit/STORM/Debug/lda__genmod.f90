        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug 22 11:39:34 2014
        MODULE LDA__genmod
          INTERFACE 
            FUNCTION LDA(KPLUS,N,R)
              REAL(KIND=8), INTENT(IN) :: KPLUS(3,3,3)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(OUT) :: R(3,3)
              LOGICAL(KIND=4) :: LDA
            END FUNCTION LDA
          END INTERFACE 
        END MODULE LDA__genmod
