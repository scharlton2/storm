        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 19 14:12:16 2018
        MODULE PART_DRY__genmod
          INTERFACE 
            FUNCTION PART_DRY(E,N,I,J)
              USE OPTIONS
              TYPE (TRIANGLE), INTENT(IN) :: E
              INTEGER(KIND=4), INTENT(OUT) :: N
              INTEGER(KIND=4), INTENT(OUT) :: I
              INTEGER(KIND=4), INTENT(OUT) :: J
              LOGICAL(KIND=4) :: PART_DRY
            END FUNCTION PART_DRY
          END INTERFACE 
        END MODULE PART_DRY__genmod
