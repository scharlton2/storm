        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 07 13:43:32 2014
        MODULE LOCAL_EDGE__genmod
          INTERFACE 
            FUNCTION LOCAL_EDGE(ELEMENT,EDG)
              USE PARAMETERS
              TYPE (TRIANGLE), INTENT(IN) :: ELEMENT
              INTEGER(KIND=4), INTENT(IN) :: EDG
              INTEGER(KIND=4) :: LOCAL_EDGE
            END FUNCTION LOCAL_EDGE
          END INTERFACE 
        END MODULE LOCAL_EDGE__genmod
