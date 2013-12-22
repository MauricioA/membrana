        !COMPILER-GENERATED INTERFACE MODULE: Fri Nov 16 13:49:44 2012
        MODULE MATXVECSIM__genmod
          INTERFACE 
            SUBROUTINE MATXVECSIM(IA,JA,AN,AD,B,C,NP)
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: IA(NP+1)
              INTEGER(KIND=4) :: JA(*)
              REAL(KIND=8) :: AN(*)
              REAL(KIND=8) :: AD(NP)
              REAL(KIND=8) :: B(NP)
              REAL(KIND=8) :: C(NP)
            END SUBROUTINE MATXVECSIM
          END INTERFACE 
        END MODULE MATXVECSIM__genmod
