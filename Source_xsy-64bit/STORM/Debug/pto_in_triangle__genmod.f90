        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug 22 11:39:29 2014
        MODULE PTO_IN_TRIANGLE__genmod
          INTERFACE 
            FUNCTION PTO_IN_TRIANGLE(E,X,Y)
              USE CONSTANTS
              TYPE (TRIANGLE) :: E
              REAL(KIND=8), INTENT(IN) :: X
              REAL(KIND=8), INTENT(IN) :: Y
              LOGICAL(KIND=4) :: PTO_IN_TRIANGLE
            END FUNCTION PTO_IN_TRIANGLE
          END INTERFACE 
        END MODULE PTO_IN_TRIANGLE__genmod
