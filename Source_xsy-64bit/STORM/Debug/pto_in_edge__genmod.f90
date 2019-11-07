        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 07 13:43:45 2014
        MODULE PTO_IN_EDGE__genmod
          INTERFACE 
            FUNCTION PTO_IN_EDGE(I,E)
              USE PARAMETERS
              INTEGER(KIND=4), INTENT(IN) :: I
              TYPE (EDGE), INTENT(IN) :: E
              LOGICAL(KIND=4) :: PTO_IN_EDGE
            END FUNCTION PTO_IN_EDGE
          END INTERFACE 
        END MODULE PTO_IN_EDGE__genmod
