        !COMPILER-GENERATED INTERFACE MODULE: Fri Nov 16 13:49:43 2012
        MODULE ARMADO__genmod
          INTERFACE 
            SUBROUTINE ARMADO(NCASE,NLE,X,Y,Z,NS,NOPE,ESM,EF,SIGMA_EL,QE&
     &)
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
            END SUBROUTINE ARMADO
          END INTERFACE 
        END MODULE ARMADO__genmod
