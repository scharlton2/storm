        !COMPILER-GENERATED INTERFACE MODULE: Wed Feb 18 13:31:09 2009
        MODULE GHOST_CELL_mod
          INTERFACE 
            SUBROUTINE GHOST_CELL(XC,YC,X1,Y1,X2,Y2,XG,YG)
              REAL(KIND=8), INTENT(IN) :: XC
              REAL(KIND=8), INTENT(IN) :: YC
              REAL(KIND=8), INTENT(IN) :: X1
              REAL(KIND=8), INTENT(IN) :: Y1
              REAL(KIND=8), INTENT(IN) :: X2
              REAL(KIND=8), INTENT(IN) :: Y2
              REAL(KIND=8), INTENT(OUT) :: XG
              REAL(KIND=8), INTENT(OUT) :: YG
            END SUBROUTINE GHOST_CELL
          END INTERFACE 
        END MODULE GHOST_CELL_mod
