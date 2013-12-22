        !COMPILER-GENERATED INTERFACE MODULE: Fri Nov 16 13:49:43 2012
        MODULE FUNCORDENALYA__genmod
          INTERFACE 
            SUBROUTINE FUNCORDENALYA(S,T,Z,PHI,DPHIX,DPHIY,DPHIZ,AJACO, &
     &AJACOI,DXHI,DTHE,DPSI)
              REAL(KIND=8) :: S
              REAL(KIND=8) :: T
              REAL(KIND=8) :: Z
              REAL(KIND=8) :: PHI(27)
              REAL(KIND=8) :: DPHIX(27)
              REAL(KIND=8) :: DPHIY(27)
              REAL(KIND=8) :: DPHIZ(27)
              REAL(KIND=8) :: AJACO(3,3)
              REAL(KIND=8) :: AJACOI(3,3)
              REAL(KIND=8) :: DXHI(27)
              REAL(KIND=8) :: DTHE(27)
              REAL(KIND=8) :: DPSI(27)
            END SUBROUTINE FUNCORDENALYA
          END INTERFACE 
        END MODULE FUNCORDENALYA__genmod
