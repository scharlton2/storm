        !COMPILER-GENERATED INTERFACE MODULE: Thu Jun 01 12:41:34 2017
        MODULE CTR2VTXWSE__genmod
          INTERFACE 
            SUBROUTINE CTR2VTXWSE(TYPE,CARRAY,VARRAY)
              USE GEOMETRY
              INTEGER(KIND=4), INTENT(IN) :: TYPE
              REAL(KIND=8), INTENT(IN) :: CARRAY(N_ELEMS)
              REAL(KIND=8), INTENT(OUT) :: VARRAY(N_PTS)
            END SUBROUTINE CTR2VTXWSE
          END INTERFACE 
        END MODULE CTR2VTXWSE__genmod
