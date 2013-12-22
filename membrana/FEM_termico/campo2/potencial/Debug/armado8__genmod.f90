        !COMPILER-GENERATED INTERFACE MODULE: Fri Nov 16 13:49:43 2012
        MODULE ARMADO8__genmod
          INTERFACE 
            SUBROUTINE ARMADO8(NCASE,NLE,X,Y,Z,NS,NOPE,ESM,EF,SIGMA_EL, &
     &QE,SOL,EX_EL,EY_EL,EZ_EL)
              INTEGER(KIND=4) :: NOPE
              INTEGER(KIND=4) :: NCASE
              INTEGER(KIND=4) :: NLE
              REAL(KIND=8) :: X(NOPE)
              REAL(KIND=8) :: Y(NOPE)
              REAL(KIND=8) :: Z(NOPE)
              INTEGER(KIND=4) :: NS(NOPE)
              REAL(KIND=8) :: ESM(NOPE,NOPE)
              REAL(KIND=8) :: EF(NOPE)
              REAL(KIND=8) :: SIGMA_EL
              REAL(KIND=8) :: QE
              REAL(KIND=8) :: SOL(NOPE)
              REAL(KIND=8) :: EX_EL
              REAL(KIND=8) :: EY_EL
              REAL(KIND=8) :: EZ_EL
            END SUBROUTINE ARMADO8
          END INTERFACE 
        END MODULE ARMADO8__genmod
