        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 07 13:43:20 2014
        MODULE PROBEFILEO__genmod
          INTERFACE 
            SUBROUTINE PROBEFILEO(T,NN,HH,UU,VV)
              INTEGER(KIND=4), INTENT(IN) :: NN
              REAL(KIND=8) :: T
              REAL(KIND=8), INTENT(IN) :: HH(NN)
              REAL(KIND=8), INTENT(IN) :: UU(NN)
              REAL(KIND=8), INTENT(IN) :: VV(NN)
            END SUBROUTINE PROBEFILEO
          END INTERFACE 
        END MODULE PROBEFILEO__genmod
