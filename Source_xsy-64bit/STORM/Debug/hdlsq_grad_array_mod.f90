        !COMPILER-GENERATED INTERFACE MODULE: Tue Feb 03 18:03:17 2009
        MODULE HDLSQ_GRAD_ARRAY_mod
          INTERFACE 
            SUBROUTINE HDLSQ_GRAD_ARRAY(PHI,GRADPHI,NDIM)
            INTEGER(KIND=4), INTENT(IN) :: NDIM
            REAL(KIND=8), INTENT(IN) :: PHI(NDIM)
            TYPE (VECTOR), INTENT(OUT) :: GRADPHI(NDIM)
          END SUBROUTINE HDLSQ_GRAD_ARRAY
        END INTERFACE 
      END MODULE HDLSQ_GRAD_ARRAY_mod
