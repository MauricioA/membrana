        !COMPILER-GENERATED INTERFACE MODULE: Fri Nov 16 13:49:43 2012
        MODULE ARMADO32D__genmod
          INTERFACE 
            SUBROUTINE ARMADO32D(NLE,X,Y,NS,NOPE,ESM,EF,EM,PR,NTENSION, &
     &ET)
              INTEGER(KIND=4) :: NTENSION
              INTEGER(KIND=4) :: NOPE
              INTEGER(KIND=4) :: NLE
              REAL(KIND=8) :: X(NOPE)
              REAL(KIND=8) :: Y(NOPE)
              INTEGER(KIND=4) :: NS(2*NOPE)
              REAL(KIND=8) :: ESM(2*NOPE,2*NOPE)
              REAL(KIND=8) :: EF(2*NOPE)
              REAL(KIND=8) :: EM
              REAL(KIND=8) :: PR
              REAL(KIND=8) :: ET(NTENSION)
            END SUBROUTINE ARMADO32D
          END INTERFACE 
        END MODULE ARMADO32D__genmod
