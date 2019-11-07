        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug 22 11:39:50 2014
        MODULE JACOBIAN__genmod
          INTERFACE 
            SUBROUTINE JACOBIAN(KI,NXI,NYI,H,U,V)
              REAL(KIND=8), INTENT(OUT) :: KI(3,3)
              REAL(KIND=8), INTENT(IN) :: NXI
              REAL(KIND=8), INTENT(IN) :: NYI
              REAL(KIND=8), INTENT(IN) :: H
              REAL(KIND=8), INTENT(IN) :: U
              REAL(KIND=8), INTENT(IN) :: V
            END SUBROUTINE JACOBIAN
          END INTERFACE 
        END MODULE JACOBIAN__genmod
