        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug 22 11:39:21 2014
        MODULE T_AREAS__genmod
          INTERFACE 
            SUBROUTINE T_AREAS(GRID,NGRID,PTO,NPTO,EDG,NEDG)
              USE PARAMETERS
              INTEGER(KIND=4), INTENT(IN) :: NEDG
              INTEGER(KIND=4), INTENT(IN) :: NPTO
              INTEGER(KIND=4), INTENT(IN) :: NGRID
              TYPE (TRIANGLE) :: GRID(NGRID)
              TYPE (POINT) :: PTO(NPTO)
              TYPE (EDGE) :: EDG(NEDG)
            END SUBROUTINE T_AREAS
          END INTERFACE 
        END MODULE T_AREAS__genmod
