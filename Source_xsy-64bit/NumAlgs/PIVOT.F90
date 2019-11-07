![  {The Matrix Inverse via Exchange Steps}
![  {The Matrix Inverse via Exchange Steps}*)
      SUBROUTINE PIVOT (A, LDA, N, B, S1, S2, IERR, MX, MY, VAL)
!
!*****************************************************************
!                                                                *
!  This subroutine calculates the inverse of a real square NxN   *
!  matrix  A  by applying exchange steps, also called the method *
!  of pivotization.                                              *
!                                                                *
!                                                                *
!  INPUT PARAMETERS:                                             *
!  =================                                             *
!  A        : 2-dimensional array A(1:LDA,1:N) containing the    *
!             matrix A that is to be inverted.                   *
!  LDA      : leading dimension of A as defined in the calling   *
!             program.                                           *
!  N        : order of the matrix A.                             *
!                                                                *
!                                                                *
!  OUTPUT PARAMETERS:                                            *
!  ==================                                            *
!  B        : 2-dimensional array B(1:LDA,1:N) containing the    *
!             inverse matrix of A.                               *
!  S1,S2    : control parameters;  S1 is the sum of the absolute *
!             values of the diagonal entries of the matrix A*B-I,*
!             where I is the nxn identity matrix.                *
!             S2 is the sum of the absolute values of the off-   *
!             diagonal entries in A*B-I. Theoretically A*B-I = 0.*
!  IERR     : = 1, inverse of A has been found.                  *
!             = 2, the matrix A is numerically singular, no      *
!                  inverse exists.                               *
!  VAL      : last pivot element, if A is numerically singular.  *
!                                                                *
!                                                                *
!  AUXILIARY VECTORS:                                            *
!  ==================                                            *
!  MX,MY    : N-vectors of INTEGER type.                         *
!                                                                *
!----------------------------------------------------------------*
!                                                                *
!  subroutines required: MACHPD                                  *
!                                                                *
!*****************************************************************
!                                                                *
!  author     : Gisela Engeln-Muellges                           *
!  date       : 05.18.1987                                       *
!  source     : FORTRAN 77                                       *
!                                                                *
!*****************************************************************
!
!  declarations.
!
      IMPLICIT DOUBLEPRECISION (A - H, O - Z)
      DOUBLEPRECISION A (LDA, N), B (LDA, N)
      INTEGER MX (N), MY (N)
!
!  calculating the machine constant FMACHP.
!
      FMACHP = 1.0D0
    5 FMACHP = 0.5D0 * FMACHP
      IF (MACHPD (1.0D0 + FMACHP) .EQ.1) GOTO 5
      FMACHP = FMACHP * 2.0D0
!
!  store matrix A in array B.
!
      DO 10 I = 1, N
         DO 20 J = 1, N
            B (I, J) = A (I, J)
   20    END DO
   10 END DO
!
!  initialize the pivot vectors MX and MY to zero.
!
      DO 30 I = 1, N
         MX (I) = 0
         MY (I) = 0
   30 END DO
!
!  determine the pivot element.
!
      DO 40 I = 1, N
         PIVO = 0.0D0
         DO 50 IX = 1, N
            IF (MX (IX) .EQ.0) THEN
               DO 60 IY = 1, N
                  IF (MY (IY) .EQ.0) THEN
                     IF (DABS (B (IX, IY) ) .GT.DABS (PIVO) ) THEN
                        PIVO = B (IX, IY)
                        NX = IX
                        NY = IY
                     ENDIF
                  ENDIF
   60          END DO
            ENDIF
   50    END DO
!
!  if the pivot element is nearly zero, the matrix is numerically singul
!
         IF (DABS (PIVO) .LT.4.0D0 * FMACHP) THEN
            VAL = PIVO
            IERR = 2
            RETURN
         ENDIF
!
!  saving the indices of the pivot element.
!
         MX (NX) = NY
         MY (NY) = NX
!
!  calculation of the matrix elements according to the
!  rules for an exchange step.
!
         DUMMY = 1.0D0 / PIVO
         DO 70 J = 1, N
            IF (J.NE.NX) THEN
               FACTOR = B (J, NY) * DUMMY
               DO 80 K = 1, N
                  B (J, K) = B (J, K) - B (NX, K) * FACTOR
   80          END DO
               B (J, NY) = - FACTOR
            ENDIF
   70    END DO
         DO 90 K = 1, N
            B (NX, K) = B (NX, K) * DUMMY
   90    END DO
         B (NX, NY) = DUMMY
   40 END DO
!
!  reverse row and column permutations.
!
      DO 100 I = 1, N - 1
         DO 110 M = I, N
            IF (MX (M) .EQ.I) GOTO 120
  110    END DO
  120    J = M
         IF (J.NE.I) THEN
            DO 130 K = 1, N
               H = B (I, K)
               B (I, K) = B (J, K)
               B (J, K) = H
  130       END DO
            MX (J) = MX (I)
            MX (I) = I
         ENDIF
         DO 140 M = I, N
            IF (MY (M) .EQ.I) GOTO 150
  140    END DO
  150    J = M
         IF (J.NE.I) THEN
            DO 160 K = 1, N
               H = B (K, I)
               B (K, I) = B (K, J)
               B (K, J) = H
  160       END DO
            MY (J) = MY (I)
            MY (I) = I
         ENDIF
  100 END DO
!
!  Forming the difference S= A*B-I, where I is the identity matrix.
!  Forming the sum S1 of the absolute values of the diagonal elements of
!  and the sum S2 of remaining elements. Theoretically
!  S1 and S2 should both equal zero.
!
      S1 = 0.0D0
      S2 = 0.0D0
      DO 170 I = 1, N
         DO 180 J = 1, N
            H = 0.0D0
            DO 190 K = 1, N
               H = H + A (I, K) * B (K, J)
  190       END DO
            IF (I.EQ.J) THEN
               S1 = S1 + DABS (H - 1.0D0)
            ELSE
               S2 = S2 + DABS (H)
            ENDIF
  180    END DO
  170 END DO
      IERR = 1
      RETURN
      END SUBROUTINE PIVOT
