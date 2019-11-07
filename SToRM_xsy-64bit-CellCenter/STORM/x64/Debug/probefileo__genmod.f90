        !COMPILER-GENERATED INTERFACE MODULE: Fri Jun 02 13:48:29 2017
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
