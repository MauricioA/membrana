        !COMPILER-GENERATED INTERFACE MODULE: Fri Nov 16 13:49:44 2012
        MODULE LIN_BCG__genmod
          INTERFACE 
            SUBROUTINE LIN_BCG(NP,IA,JA,A_SPA,AD,B,X,UNIT_CONT)
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: IA(NP+1)
              INTEGER(KIND=4) :: JA(*)
              REAL(KIND=8) :: A_SPA(*)
              REAL(KIND=8) :: AD(NP)
              REAL(KIND=8) :: B(NP)
              REAL(KIND=8) :: X(NP)
              INTEGER(KIND=4) :: UNIT_CONT
            END SUBROUTINE LIN_BCG
          END INTERFACE 
        END MODULE LIN_BCG__genmod
