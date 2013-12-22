        !COMPILER-GENERATED INTERFACE MODULE: Fri Nov 16 13:49:43 2012
        MODULE CAMPO__genmod
          INTERFACE 
            SUBROUTINE CAMPO(NNODES,NELEMENTS,NODPEL,SOLUTION,MATERIAL, &
     &CONECT,COOR_X,COOR_Y,COOR_Z,SIGMA1,SIGMA2,SIGMA3,SIGMA4,GRAD_X,   &
     &GRAD_Y,GRAD_Z,GRADXEL_X,GRADXEL_Y,GRADXEL_Z)
              INTEGER(KIND=4) :: NODPEL
              INTEGER(KIND=4) :: NELEMENTS
              INTEGER(KIND=4) :: NNODES
              REAL(KIND=8) :: SOLUTION(NNODES)
              INTEGER(KIND=4) :: MATERIAL(NELEMENTS)
              INTEGER(KIND=4) :: CONECT(NELEMENTS,NODPEL)
              REAL(KIND=8) :: COOR_X(NNODES)
              REAL(KIND=8) :: COOR_Y(NNODES)
              REAL(KIND=8) :: COOR_Z(NNODES)
              REAL(KIND=8) :: SIGMA1
              REAL(KIND=8) :: SIGMA2
              REAL(KIND=8) :: SIGMA3
              REAL(KIND=8) :: SIGMA4
              REAL(KIND=8) :: GRAD_X(NNODES)
              REAL(KIND=8) :: GRAD_Y(NNODES)
              REAL(KIND=8) :: GRAD_Z(NNODES)
              REAL(KIND=8) :: GRADXEL_X(NELEMENTS)
              REAL(KIND=8) :: GRADXEL_Y(NELEMENTS)
              REAL(KIND=8) :: GRADXEL_Z(NELEMENTS)
            END SUBROUTINE CAMPO
          END INTERFACE 
        END MODULE CAMPO__genmod
