        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 19 14:12:16 2018
        MODULE KMATRICES__genmod
          INTERFACE 
            SUBROUTINE KMATRICES(KPLUS,KMINUS,NXI,NYI,H,U,V)
              REAL(KIND=8), INTENT(OUT) :: KPLUS(3,3)
              REAL(KIND=8), INTENT(OUT) :: KMINUS(3,3)
              REAL(KIND=8), INTENT(IN) :: NXI
              REAL(KIND=8), INTENT(IN) :: NYI
              REAL(KIND=8), INTENT(IN) :: H
              REAL(KIND=8), INTENT(IN) :: U
              REAL(KIND=8), INTENT(IN) :: V
            END SUBROUTINE KMATRICES
          END INTERFACE 
        END MODULE KMATRICES__genmod
