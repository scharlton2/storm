        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug 22 11:39:36 2014
        MODULE NONLINEARB__genmod
          INTERFACE 
            FUNCTION NONLINEARB(KPLUS,KMINUS,N,R)
              REAL(KIND=8), INTENT(IN) :: KPLUS(3,3,3)
              REAL(KIND=8), INTENT(IN) :: KMINUS(3,3,3)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(OUT) :: R(3,3)
              LOGICAL(KIND=4) :: NONLINEARB
            END FUNCTION NONLINEARB
          END INTERFACE 
        END MODULE NONLINEARB__genmod
