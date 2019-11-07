        !COMPILER-GENERATED INTERFACE MODULE: Mon Feb 02 12:45:40 2009
        MODULE JACOBIAN_mod
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
        END MODULE JACOBIAN_mod
