        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug 22 11:39:49 2014
        MODULE NARROW__genmod
          INTERFACE 
            FUNCTION NARROW(KPLUS,KMINUS,N,R)
              REAL(KIND=8), INTENT(IN) :: KPLUS(3,3,3)
              REAL(KIND=8), INTENT(IN) :: KMINUS(3,3,3)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(OUT) :: R(3,3)
              LOGICAL(KIND=4) :: NARROW
            END FUNCTION NARROW
          END INTERFACE 
        END MODULE NARROW__genmod
