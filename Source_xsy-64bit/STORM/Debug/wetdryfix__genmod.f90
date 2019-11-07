        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug 22 11:39:18 2014
        MODULE WETDRYFIX__genmod
          INTERFACE 
            SUBROUTINE WETDRYFIX(ITYPE,HV,HC,GRADHC)
              USE GEOMETRY
              USE PARAMETERS
              INTEGER(KIND=4), INTENT(IN) :: ITYPE
              REAL(KIND=8), INTENT(INOUT) :: HV(N_PTS)
              REAL(KIND=8), INTENT(IN) :: HC(N_ELEMS)
              TYPE (VECTOR), INTENT(IN) :: GRADHC(N_ELEMS)
            END SUBROUTINE WETDRYFIX
          END INTERFACE 
        END MODULE WETDRYFIX__genmod
