        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 19 14:12:25 2018
        MODULE ADJUSTZB__genmod
          INTERFACE 
            SUBROUTINE ADJUSTZB(FLAG,E,GRADZB,GRADZB2)
              USE PARAMETERS
              LOGICAL(KIND=4), INTENT(IN) :: FLAG
              INTEGER(KIND=4), INTENT(IN) :: E
              TYPE (VECTOR), INTENT(OUT) :: GRADZB
              TYPE (VECTOR), INTENT(OUT) :: GRADZB2
            END SUBROUTINE ADJUSTZB
          END INTERFACE 
        END MODULE ADJUSTZB__genmod
