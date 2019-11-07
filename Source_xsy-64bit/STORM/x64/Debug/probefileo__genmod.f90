        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 19 14:12:18 2018
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
