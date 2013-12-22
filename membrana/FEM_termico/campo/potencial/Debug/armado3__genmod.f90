        !COMPILER-GENERATED INTERFACE MODULE: Fri Nov 16 13:49:43 2012
        MODULE ARMADO3__genmod
          INTERFACE 
            SUBROUTINE ARMADO3(NLE,X,Y,NS,NOPE,ESM,EF,SIGMA_EL,QE,SOL,  &
     &EX_EL,EY_EL)
              INTEGER(KIND=4) :: NOPE
              INTEGER(KIND=4) :: NLE
              REAL(KIND=8) :: X(NOPE)
              REAL(KIND=8) :: Y(NOPE)
              INTEGER(KIND=4) :: NS(NOPE)
              REAL(KIND=8) :: ESM(NOPE,NOPE)
              REAL(KIND=8) :: EF(NOPE)
              REAL(KIND=8) :: SIGMA_EL
              REAL(KIND=8) :: QE
              REAL(KIND=8) :: SOL(NOPE)
              REAL(KIND=8) :: EX_EL
              REAL(KIND=8) :: EY_EL
            END SUBROUTINE ARMADO3
          END INTERFACE 
        END MODULE ARMADO3__genmod
