        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 07 13:43:19 2014
        MODULE ADD_EDGE_TO_ARRAY__genmod
          INTERFACE 
            SUBROUTINE ADD_EDGE_TO_ARRAY(EDG,EDGES_ARRAY,N)
              USE PARAMETERS
              INTEGER(KIND=4), INTENT(INOUT) :: N
              TYPE (EDGE), INTENT(IN) :: EDG
              TYPE (EDGE), INTENT(INOUT) :: EDGES_ARRAY(N+1)
            END SUBROUTINE ADD_EDGE_TO_ARRAY
          END INTERFACE 
        END MODULE ADD_EDGE_TO_ARRAY__genmod
