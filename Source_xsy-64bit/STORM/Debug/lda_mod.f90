        !COMPILER-GENERATED INTERFACE MODULE: Tue Feb 03 18:02:55 2009
        MODULE LDA_mod
          INTERFACE 
            FUNCTION LDA(KPLUS,N,R) RESULT(LDA_0)
              REAL(KIND=8), INTENT(IN) :: KPLUS(3,3,3)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(OUT) :: R(3,3)
              LOGICAL(KIND=4) :: LDA_0
            END FUNCTION LDA
          END INTERFACE 
        END MODULE LDA_mod
