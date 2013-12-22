        !COMPILER-GENERATED INTERFACE MODULE: Fri Nov 16 13:49:44 2012
        MODULE CG__genmod
          INTERFACE 
            SUBROUTINE CG(NP,IA,JA,A_SPA,AD,B,X,TOL,ITMAX,ITER,ERR)
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: IA(*)
              INTEGER(KIND=4) :: JA(*)
              REAL(KIND=8) :: A_SPA(*)
              REAL(KIND=8) :: AD(*)
              REAL(KIND=8) :: B(*)
              REAL(KIND=8) :: X(*)
              REAL(KIND=8) :: TOL
              INTEGER(KIND=4) :: ITMAX
              INTEGER(KIND=4) :: ITER
              REAL(KIND=8) :: ERR
            END SUBROUTINE CG
          END INTERFACE 
        END MODULE CG__genmod
