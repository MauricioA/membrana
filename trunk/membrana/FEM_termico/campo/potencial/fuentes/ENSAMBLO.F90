!   CALCULA LAS MATRICES DE CADA ELEMENTO
SUBROUTINE ARMADO3(NLE,X,Y,ns,NOPE,ESM,EF,sigma_el,qe,sol,Ex_el,Ey_el)

!CALL ARMADO3(JEL,X,Y,ns,nodpel,ESM,EF,sigma_el,qe,sol,gradxel_x(jel),gradxel_y(jel))

implicit none
INTEGER :: nope,NS(NOPE),NLE
DOUBLE PRECISION :: X(NOPE),Y(NOPE),EF(NOPE),ESM(NOPE,NOPE),sigma_el,qe,sol(nope),Ex_el,Ey_el

! NOPE = 3

DOUBLE PREcIsION:: B(3),C(3),DETER,A,RMED
DOUBLE PREcIsION:: PI=3.14159
integer i,j
	 
B(1)=Y(2)-Y(3)
B(2)=Y(3)-Y(1)
B(3)=Y(1)-Y(2)
C(1)=X(3)-X(2)
C(2)=X(1)-X(3)
C(3)=X(2)-X(1)

DETER=X(2)*Y(3)+X(3)*Y(1)+X(1)*Y(2)-X(2)*Y(1)-X(3)*Y(2)-X(1)*Y(3)
RMED=(X(1)+X(2)+X(3))/3.

IF(ABS(DETER)<0.0000001) then
    write(6,*) 'Area del elemento : ',  NLE,'   es cero ', DETER 
   stop ' '
endif

!borrar
if (NLE.eq.1) then
	write(6,*) "deter ", DETER
	write(6,*) "rmed ", RMED
	write(6,*) "b ", B
	write(6,*) "c ", C
	write(6,*) "sigma ", sigma_el
	write(6,*) 200e-9
endif
 
do I=1,3
  EF(I)=0.0
  do J=1,3
      A=1.0
      IF(I.EQ.J) A=2.0
      EF(I)=EF(I)+(QE*DETER*PI/12.)*(A*X(J))
      ESM(I,J)=(sigma_el*B(I)*B(J)+sigma_el*C(I)*C(J))*PI*RMED/DETER
  enddo     
enddo      

end subroutine armado3


subroutine ARMADO32d(NLE,X,Y,ns,nope,ESM,EF,EM,PR,ntension,et)
implicit none
INTEGER :: nope,NS(2*NOPE),NLE,ntension
DOUBLE PRECISION :: X(NOPE),Y(NOPE),EF(2*NOPE),ESM(2*NOPE,2*NOPE),EM,PR,et(ntension)

DOUBLE PREcIsION:: B(ntension,3*nope),C(2*nope,ntension),D(ntension,ntension),aux(nope),DETER,A,RMED,zmed
DOUBLE PREcIsION:: PI=3.14159, R,sum
integer i,j,k
	 
 B=0.0
 c=0.0

 B(1,1)=Y(2)-Y(3)
 B(1,3)=Y(3)-Y(1)
 B(1,5)=Y(1)-Y(2)
 B(3,2)=X(3)-X(2)
 B(3,4)=X(1)-X(3)
 B(3,6)=X(2)-X(1)
 B(4,1)=B(3,2)
 B(4,2)=B(1,1)
 B(4,3)=B(3,4)
 B(4,4)=B(1,3)
 B(4,5)=B(3,6)
 B(4,6)=B(1,5)
 
 deter=X(2)*Y(3)+X(3)*Y(1)+X(1)*Y(2)-X(2)*Y(1)-X(3)*Y(2)-X(1)*Y(3)     
 AUX(1)=X(2)*Y(3)-X(3)*Y(2)
 AUX(2)=X(3)*Y(1)-X(1)*Y(3)
 AUX(3)=X(1)*Y(2)-X(2)*Y(1)
 RMED=(X(1)+X(2)+X(3))/3.
 ZMED=(Y(1)+Y(2)+Y(3))/3.
      
 B(2,1)=(AUX(1)+B(1,1)*RMED+B(3,2)*ZMED)/RMED
 B(2,3)=(AUX(2)+B(1,3)*RMED+B(3,4)*ZMED)/RMED
 B(2,5)=(AUX(3)+B(1,5)*RMED+B(3,6)*ZMED)/RMED
 



      R=EM/(1.+PR)
      D(1,1)=R*(1.-PR)/(1-2*PR)
      D(2,2)=D(1,1)
      D(3,3)=D(1,1)
      D(4,4)=R/2.
      D(1,2)=PR*R/(1.-2*PR)
      D(2,1)=D(1,2)
      D(1,3)=D(1,2)
      D(3,1)=D(1,2)
      D(2,3)=D(1,2)
      D(3,2)=D(1,2)
      D(1,4)=0.0
      D(2,4)=0.0
      D(3,4)=0.0
      D(4,1)=0.0
      D(4,2)=0.0
      D(4,3)=0.0
  
!    %MATRIZ [C]=[BT]*[D]
      DO I= 1,2*nope
        DO J= 1,ntension
          C(I,J) = 0.0
          DO K= 1,ntension
              C(I,J) = C(I,J)+B(K,I)*D(K,J)
          enddo
        enddo
      enddo

!    %MATRIZ [ESM]=[BT]*[D]*[B]=[C]*[B]
     DO  I=1,2*nope
        DO  J=1,2*nope
          SUM=0.0
          DO K=1,ntension
            SUM=SUM + C(I,K)*B(K,J)
          enddo
         
          ESM(I,J)=SUM*RMED*PI/deter
       enddo
    enddo
   
   
 !  C    %MATRIZ FUERZA TERMICA [FT]=[C]*[ET]
      DO I=1,2*nope
         EF(I)=(C(I,1)*ET(1) + C(I,2)*ET(2) + C(I,3)*ET(3) + C(I,4)*ET(4))*RMED*PI
      enddo
   
    
end subroutine armado32d

SUBROUTINE ARMADO8(NCASE,NLE,X,Y,Z,ns,NOPE,ESM,EF,sigma_el,qe,sol,Ex_el,Ey_el,Ez_el)
implicit none
INTEGER :: nope,NS(NOPE),NCASE
DOUBLE PRECISION :: X(NOPE),Y(NOPE),Z(NOPE),EF(NOPE),ESM(NOPE,NOPE),sigma_el,qe,sol(nope),Ex_el,Ey_el,Ez_el

DOUBLE PREcIsION:: PHI(8),DPHIX(8),DPHIY(8),DPHIZ(8),AJACO(3,3),AJACOI(3,3),DXHI(8), &
     DTHE(8),DPSI(8),GAUSSPT(2),GAUSSWT(2),XHI,THE,PSI,DETER,S11,S22,S33,CNST

integer kk,jj,i,j,ii,K,NLE,pgaus,I2
	 
DATA GAUSSPT/ -0.57735027, 0.57735027/
DATA GAUSSWT/ 1.0, 1.0 /


ESM=0.0
EF=0.0
  pgaus=0
  DO KK=1,2
    DO JJ=1,2
      DO II=1,2
        pgaus=pgaus+1

		XHI = GAUSSPT(KK)
        THE = GAUSSPT(JJ)
        PSI = GAUSSPT(II)
      
!  CALCULO LA FUNCION PHI, SU DERIVADA EN FUNCCION DE XHI Y DE THE

        call funciones8(XHI,THE,PSI,PHI,DPHIX,DPHIY,DPHIZ,AJACO,AJACOI,DXHI,DTHE,DPSI ) !!ojo orden distinto!!
       ! call funcORDENALYA(XHI,THE,PSI,PHI,DPHIX,DPHIY,DPHIZ,AJACO,AJACOI,DXHI,DTHE,DPSI )
  
!   JACOBIANO Y SU INVERSA
        AJACO=0.0

          DO K=1,NOPE
              AJACO(1,1)=AJACO(1,1)+DXHI(K)*X(K)
              AJACO(1,2)=AJACO(1,2)+DXHI(K)*Y(K)
              AJACO(1,3)=AJACO(1,3)+DXHI(K)*Z(K)
              AJACO(2,1)=AJACO(2,1)+DTHE(K)*X(K)
              AJACO(2,2)=AJACO(2,2)+DTHE(K)*Y(K)
              AJACO(2,3)=AJACO(2,3)+DTHE(K)*Z(K)
              AJACO(3,1)=AJACO(3,1)+DPSI(K)*X(K)
              AJACO(3,2)=AJACO(3,2)+DPSI(K)*Y(K)
              AJACO(3,3)=AJACO(3,3)+DPSI(K)*Z(K)
          ENDDO

          CALL DETERM(AJACO, AJACOI, DETER)

          IF(DETER.EQ.0.0) THEN
              WRITE(6,*) ' SUBRUTINE: ARMADO '
              WRITE(6,*) ' ELEMENTO: ',NLE,' PUNTO ',KK,JJ,II,' TIENE DETERMINANTE CERO',DETER
              STOP ' '
          ENDIF

 
!   DERIVADAS DE LAS FUNCIONES DE INTERPOLACION EN FUNCION DE X-Y-Z

          DO I=1,NOPE
            DPHIX(I)=AJACOI(1,1)*DXHI(I) + AJACOI(1,2)*DTHE(I) + AJACOI(1,3)*DPSI(I)
            DPHIY(I)=AJACOI(2,1)*DXHI(I) + AJACOI(2,2)*DTHE(I) + AJACOI(2,3)*DPSI(I)
            DPHIZ(I)=AJACOI(3,1)*DXHI(I) + AJACOI(3,2)*DTHE(I) + AJACOI(3,3)*DPSI(I)
          ENDDO

!       DIFERENCIAL DE INTEGRACION
          CNST = DETER * GAUSSWT(JJ)*GAUSSWT(KK)*GAUSSWT(II)

!     MATRIZ DE RIGIDEZ LOCAL Y VETOR INDEPENDIENTE

          DO I=1,nope
            DO J=1,nope
              S11= DPHIX(I)*DPHIX(J)
              S22= DPHIY(I)*DPHIY(J)
              S33= DPHIZ(I)*DPHIZ(J)

              ESM(I,J)=ESM(I,J) + (S11+S22+S33)*CNST*sigma_el
            ENDDO
            Ex_el=Ex_el- DPHIX(i)*sol(i)/real(nope)
            Ey_el=Ey_el- DPHIY(i)*sol(i)/real(nope)
            Ez_el=Ez_el- DPHIZ(i)*sol(i)/real(nope)

            !EF(I)=EF(I) + qe*phi(i)*cnst

          ENDDO
        ENDDO
      ENDDO
    ENDDO


    

end subroutine armado8



SUBROUTINE ARMADO(NCASE,NLE,X,Y,Z,ns,NOPE,ESM,EF,sigma_el,qe)
implicit none
INTEGER :: nope,NS(NOPE),NCASE
DOUBLE PRECISION :: X(NOPE),Y(NOPE),Z(NOPE),EF(NOPE),ESM(NOPE,NOPE),sigma_el,qe

DOUBLE PREcIsION:: PHI(27),DPHIX(27),DPHIY(27),DPHIZ(27),AJACO(3,3),AJACOI(3,3),DXHI(27), &
     DTHE(27),DPSI(27),GAUSSPT(3),GAUSSWT(3),XHI,THE,PSI,DETER,S11,S22,S33,CNST

integer kk,jj,i,j,ii,K,NLE,pgaus,I2
	 
DATA GAUSSPT/ -0.774596669241483377035853079956, 0.0, 0.774596669241483377035853079956/
DATA GAUSSWT/ 0.5555555556, 0.8888888889 ,0.5555555556 /

ESM=0.0
EF=0.0
  pgaus=0
  DO KK=1,3
    DO JJ=1,3
      DO II=1,3
        pgaus=pgaus+1

		XHI = GAUSSPT(KK)
        THE = GAUSSPT(JJ)
        PSI = GAUSSPT(II)
      
!  CALCULO LA FUNCION PHI, SU DERIVADA EN FUNCCION DE XHI Y DE THE

        call funciones(XHI,THE,PSI,PHI,DPHIX,DPHIY,DPHIZ,AJACO,AJACOI,DXHI,DTHE,DPSI ) !!ojo orden distinto!!
       ! call funcORDENALYA(XHI,THE,PSI,PHI,DPHIX,DPHIY,DPHIZ,AJACO,AJACOI,DXHI,DTHE,DPSI )
  
!   JACOBIANO Y SU INVERSA
        AJACO=0.0

          DO K=1,NOPE
              AJACO(1,1)=AJACO(1,1)+DXHI(K)*X(K)
              AJACO(1,2)=AJACO(1,2)+DXHI(K)*Y(K)
              AJACO(1,3)=AJACO(1,3)+DXHI(K)*Z(K)
              AJACO(2,1)=AJACO(2,1)+DTHE(K)*X(K)
              AJACO(2,2)=AJACO(2,2)+DTHE(K)*Y(K)
              AJACO(2,3)=AJACO(2,3)+DTHE(K)*Z(K)
              AJACO(3,1)=AJACO(3,1)+DPSI(K)*X(K)
              AJACO(3,2)=AJACO(3,2)+DPSI(K)*Y(K)
              AJACO(3,3)=AJACO(3,3)+DPSI(K)*Z(K)
          ENDDO

          CALL DETERM(AJACO, AJACOI, DETER)

          IF(DETER.EQ.0.0) THEN
              WRITE(6,*) ' SUBRUTINE: ARMADO '
              WRITE(6,*) ' ELEMENTO: ',NLE,' PUNTO ',KK,JJ,II,' TIENE DETERMINANTE CERO',DETER
              STOP ' '
          ENDIF

 
!   DERIVADAS DE LAS FUNCIONES DE INTERPOLACION EN FUNCION DE X-Y-Z

          DO I=1,27
            DPHIX(I)=AJACOI(1,1)*DXHI(I) + AJACOI(1,2)*DTHE(I) + AJACOI(1,3)*DPSI(I)
            DPHIY(I)=AJACOI(2,1)*DXHI(I) + AJACOI(2,2)*DTHE(I) + AJACOI(2,3)*DPSI(I)
            DPHIZ(I)=AJACOI(3,1)*DXHI(I) + AJACOI(3,2)*DTHE(I) + AJACOI(3,3)*DPSI(I)
          ENDDO

!       DIFERENCIAL DE INTEGRACION
          CNST = DETER * GAUSSWT(JJ)*GAUSSWT(KK)*GAUSSWT(II)

!     MATRIZ DE RIGIDEZ LOCAL Y VETOR INDEPENDIENTE

          DO I=1,nope
            DO J=1,nope
              S11= DPHIX(I)*DPHIX(J)
              S22= DPHIY(I)*DPHIY(J)
              S33= DPHIZ(I)*DPHIZ(J)

              ESM(I,J)=ESM(I,J) + (S11+S22+S33)*CNST*sigma_el

            ENDDO

            EF(I)=EF(I) + qe*phi(i)*cnst

          ENDDO
        ENDDO
      ENDDO
    ENDDO


end subroutine armado

subroutine funciones8(s,t,z,PHI,DPHIX,DPHIY,DPHIZ,AJACO,AJACOI,DXHI,DTHE,DPSI )   
IMPLICIT none
double precision PHI(8),DPHIX(8),DPHIY(8),DPHIZ(8),AJACO(3,3),AJACOI(3,3),DXHI(8), &
          DTHE(8),DPSI(8),s,t,z
double precision sm,tm,zm,sq,tp,zp


             
             
                 sm = 0.5*(1.0-s)
                 tm = 0.5*(1.0-t)
                 zm = 0.5*(1.0-z)
                 sq = 0.5*(1.0+s)
                 tp = 0.5*(1.0+t)
                 zp = 0.5*(1.0+z)
                 

                  phi(1)    = sm*tm*zm
                  DXHI(1)  =-0.5*tm*zm
                  DTHE(1)  =-0.5*sm*zm
                  DPSI(1)=-0.5*sm*tm


                 phi(2)     = sq*tm*zm
                 dxhi( 2) = 0.5*tm*zm
                 dthe( 2) =-0.5*sq*zm
                 dpsi( 2) =-0.5*sq*tm
                 
                 phi(3)     = sq*tp*zm
                 DXHI( 3) = 0.5*tp*zm
                 dthe( 3) = 0.5*sq*zm
                 dpsi( 3) =-0.5*sq*tp
                 
                 phi(4)     = sm*tp*zm
                 DXHI( 4) =-0.5*tp*zm
                 dthe( 4) = 0.5*sm*zm
                 dpsi( 4) =-0.5*sm*tp
                 
                 phi( 5)    = sm*tm*zp
                 DXHI( 5) =-0.5*tm*zp
                 dthe( 5) =-0.5*sm*zp
                 dpsi( 5) = 0.5*sm*tm
                 
                 phi(6)    = sq*tm*zp 
                 DXHI(6) = 0.5*tm*zp
                 dthe( 6) =-0.5*sq*zp
                 dpsi( 6) = 0.5*sq*tm
                 
                 phi(  7)   = sq*tp*zp
                 DXHI(7) = 0.5*tp*zp
                 dthe( 7) = 0.5*sq*zp
                 dpsi( 7) = 0.5*sq*tp

                 phi(8)     = sm*tp*zp
                 DXHI( 8) =-0.5*tp*zp
                 dthe( 8) = 0.5*sm*zp
                 dpsi( 8) = 0.5*sm*tp

	 
	 
return
end


subroutine funciones(s,t,z,PHI,DPHIX,DPHIY,DPHIZ,AJACO,AJACOI,DXHI,DTHE,DPSI )   
IMPLICIT none
double precision PHI(27),DPHIX(27),DPHIY(27),DPHIZ(27),AJACO(3,3),AJACOI(3,3),DXHI(27), &
          DTHE(27),DPSI(27),s,t,z
double precision s1,z1,t1,sl,tl,zl,sq,tp,zp,s2,t2,z2,s3,t3,z3,s4,t4,z4


	 sl=s*(s-1.0)
     tl=t*(t-1.0)
     zl=z*(z-1.0)
     sq=s*(s+1.0)
     tp=t*(t+1.0)
     zp=z*(z+1.0)
     s1= 2.0*s-1.0
     t1= 2.0*t-1.0
     z1= 2.0*z-1.0
     s2= 1.0-s*s
     t2= 1.0-t*t
     z2= 1.0-z*z
     s3= 1.0+2.0*s
     t3= 1.0+2.0*t
     z3= 1.0+2.0*z
     s4=-2.0*s
     t4=-2.0*t
     z4=-2.0*z
     
	 phi (1) = 0.125*sl*tl*zl
     DXHI(1) = 0.125*s1*tl*zl
     DTHE(1) = 0.125*sl*t1*zl
     DPSI(1) = 0.125*sl*tl*z1
     
	 phi (2) = 0.25*s2*tl*zl
     DXHI(2) = 0.25*s4*tl*zl
     DTHE(2) = 0.25*s2*t1*zl
     DPSI(2) = 0.25*s2*tl*z1
     
     phi (3) = 0.125*sq*tl*zl
     DXHI(3) = 0.125*s3*tl*zl
     DTHE(3) = 0.125*sq*t1*zl
     DPSI(3) = 0.125*sq*tl*z1
     
 	 phi (4) = 0.25*sq*t2*zl
     DXHI(4) = 0.25*s3*t2*zl
     DTHE(4) = 0.25*sq*t4*zl
     DPSI(4) = 0.25*sq*t2*z1
     

	 phi (5) = 0.125*sq*tp*zl
     DXHI(5) = 0.125*s3*tp*zl
     DTHE(5) = 0.125*sq*t3*zl
     DPSI(5) = 0.125*sq*tp*z1
     
     phi (6) = 0.25*s2*tp*zl
     DXHI(6) = 0.25*s4*tp*zl
     DTHE(6) = 0.25*s2*t3*zl
     DPSI(6) = 0.25*s2*tp*z1

	 phi (7) = 0.125*sl*tp*zl
     DXHI(7) = 0.125*s1*tp*zl
     DTHE(7) = 0.125*sl*t3*zl
     DPSI(7) = 0.125*sl*tp*z1
       
     phi (8) = 0.25*sl*t2*zl
     DXHI(8) = 0.25*s1*t2*zl
     DTHE(8) = 0.25*sl*t4*zl
     DPSI(8) = 0.25*sl*t2*z1
     
     phi (9) = 0.5*s2*t2*zl
     DXHI(9) = 0.5*s4*t2*zl
     DTHE(9) = 0.5*s2*t4*zl
     DPSI(9) = 0.5*s2*t2*z1

		 
     phi (10) = 0.25*sl*tl*z2
     DXHI(10) = 0.25*s1*tl*z2
     DTHE(10) = 0.25*sl*t1*z2
     DPSI(10) = 0.25*sl*tl*z4

     phi (11) = 0.5*s2*tl*z2
     DXHI(11) = 0.5*s4*tl*z2
     DTHE(11) = 0.5*s2*t1*z2
     DPSI(11) = 0.5*s2*tl*z4
     
	 	      
	 phi (12) = 0.25*sq*tl*z2
     DXHI(12) = 0.25*s3*tl*z2
     DTHE(12) = 0.25*sq*t1*z2
     DPSI(12) = 0.25*sq*tl*z4
     
	 phi (13) = 0.5*sq*t2*z2
     DXHI(13) = 0.5*s3*t2*z2
     DTHE(13) = 0.5*sq*t4*z2
     DPSI(13) = 0.5*sq*t2*z4
     
	 
	 phi (14) = 0.25*sq*tp*z2
     DXHI(14) = 0.25*s3*tp*z2
     DTHE(14) = 0.25*sq*t3*z2
     DPSI(14) = 0.25*sq*tp*z4
     
	 phi (15) = 0.5*s2*tp*z2
     DXHI(15) = 0.5*s4*tp*z2
     DTHE(15) = 0.5*s2*t3*z2
     DPSI(15) = 0.5*s2*tp*z4
     
	 phi (16) = 0.25*sl*tp*z2
     DXHI(16) = 0.25*s1*tp*z2
     DTHE(16) = 0.25*sl*t3*z2
     DPSI(16) = 0.25*sl*tp*z4
     
	 phi (17) = 0.5*sl*t2*z2
     DXHI(17) = 0.5*s1*t2*z2
     DTHE(17) = 0.5*sl*t4*z2
     DPSI(17) = 0.5*sl*t2*z4
     
	 phi (18) = s2*t2*z2
     DXHI(18) = s4*t2*z2
     DTHE(18) = s2*t4*z2
     DPSI(18) = s2*t2*z4
	 
	 phi (19) = 0.125*sl*tl*zp
     DXHI(19) = 0.125*s1*tl*zp
     DTHE(19) = 0.125*sl*t1*zp
     DPSI(19) = 0.125*sl*tl*z3

     phi (20) = 0.25*s2*tl*zp
     DXHI(20) = 0.25*s4*tl*zp
     DTHE(20) = 0.25*s2*t1*zp
     DPSI(20) = 0.25*s2*tl*z3

	 phi (21) = 0.125*sq*tl*zp
     DXHI(21) = 0.125*s3*tl*zp
     DTHE(21) = 0.125*sq*t1*zp
     DPSI(21) = 0.125*sq*tl*z3
     
     phi (22) = 0.25*sq*t2*zp
     DXHI(22) = 0.25*s3*t2*zp
     DTHE(22) = 0.25*sq*t4*zp
     DPSI(22) = 0.25*sq*t2*z3
     
	 phi (23) = 0.125*sq*tp*zp
     DXHI(23) = 0.125*s3*tp*zp
     DTHE(23) = 0.125*sq*t3*zp
     DPSI(23) = 0.125*sq*tp*z3
     
	 phi (24) = 0.25*s2*tp*zp
     DXHI(24) = 0.25*s4*tp*zp
     DTHE(24) = 0.25*s2*t3*zp
     DPSI(24) = 0.25*s2*tp*z3
     
	 phi (25) = 0.125*sl*tp*zp
     DXHI(25) = 0.125*s1*tp*zp
     DTHE(25) = 0.125*sl*t3*zp
     DPSI(25) = 0.125*sl*tp*z3
     
	 phi (26) = 0.25*sl*t2*zp
     DXHI(26) = 0.25*s1*t2*zp
     DTHE(26) = 0.25*sl*t4*zp
     DPSI(26) = 0.25*sl*t2*z3
     
	 phi (27) = 0.5*s2*t2*zp
     DXHI(27) = 0.5*s4*t2*zp
     DTHE(27) = 0.5*s2*t4*zp
     DPSI(27) = 0.5*s2*t2*z3
	
	
return
end

subroutine funcORDENALYA(s,t,z,PHI,DPHIX,DPHIY,DPHIZ,AJACO,AJACOI,DXHI,DTHE,DPSI )   
IMPLICIT none
double precision PHI(27),DPHIX(27),DPHIY(27),DPHIZ(27),AJACO(3,3),AJACOI(3,3),DXHI(27), &
          DTHE(27),DPSI(27),s,t,z
double precision s1,z1,t1,sl,tl,zl,sq,tp,zp,s2,t2,z2,s3,t3,z3,s4,t4,z4

	 sl=s*(s-1.0)
     tl=t*(t-1.0)
     zl=z*(z-1.0)
     sq=s*(s+1.0)
     tp=t*(t+1.0)
     zp=z*(z+1.0)
     s1= 2.0*s-1.0
     t1= 2.0*t-1.0
     z1= 2.0*z-1.0
     s2= 1.0-s*s
     t2= 1.0-t*t
     z2= 1.0-z*z
     s3= 1.0+2.0*s
     t3= 1.0+2.0*t
     z3= 1.0+2.0*z
     s4=-2.0*s
     t4=-2.0*t
     z4=-2.0*z
     
	 phi (1) = 0.125*sl*tl*zl
     DXHI(1) = 0.125*s1*tl*zl
     DTHE(1) = 0.125*sl*t1*zl
     DPSI(1) = 0.125*sl*tl*z1
     
     phi (2) = 0.125*sq*tl*zl
     DXHI(2) = 0.125*s3*tl*zl
     DTHE(2) = 0.125*sq*t1*zl
     DPSI(2) = 0.125*sq*tl*z1
     
	 phi (3) = 0.125*sq*tp*zl
     DXHI(3) = 0.125*s3*tp*zl
     DTHE( 3) = 0.125*sq*t3*zl
     DPSI( 3) = 0.125*sq*tp*z1
     
	 phi ( 4) = 0.125*sl*tp*zl
     DXHI( 4) = 0.125*s1*tp*zl
     DTHE( 4) = 0.125*sl*t3*zl
     DPSI( 4) = 0.125*sl*tp*z1
     
	 phi ( 5) = 0.125*sl*tl*zp
     DXHI( 5) = 0.125*s1*tl*zp
     DTHE( 5) = 0.125*sl*t1*zp
     DPSI( 5) = 0.125*sl*tl*z3
     
	 phi ( 6) = 0.125*sq*tl*zp
     DXHI( 6) = 0.125*s3*tl*zp
     DTHE( 6) = 0.125*sq*t1*zp
     DPSI( 6) = 0.125*sq*tl*z3
     
	 phi ( 7) = 0.125*sq*tp*zp
     DXHI( 7) = 0.125*s3*tp*zp
     DTHE( 7) = 0.125*sq*t3*zp
     DPSI( 7) = 0.125*sq*tp*z3
     
	 phi ( 8) = 0.125*sl*tp*zp
     DXHI( 8) = 0.125*s1*tp*zp
     DTHE( 8) = 0.125*sl*t3*zp
     DPSI( 8) = 0.125*sl*tp*z3
     
	 phi ( 9) = 0.25*s2*tl*zl
     DXHI( 9) = 0.25*s4*tl*zl
     DTHE( 9) = 0.25*s2*t1*zl
     DPSI( 9) = 0.25*s2*tl*z1
     
	 phi (10) = 0.25*sq*t2*zl
     DXHI(10) = 0.25*s3*t2*zl
     DTHE(10) = 0.25*sq*t4*zl
     DPSI(10) = 0.25*sq*t2*z1
     
	 phi (11) = 0.25*s2*tp*zl
     DXHI(11) = 0.25*s4*tp*zl
     DTHE(11) = 0.25*s2*t3*zl
     DPSI(11) = 0.25*s2*tp*z1
     
	 phi (12) = 0.25*sl*t2*zl
     DXHI(12) = 0.25*s1*t2*zl
     DTHE(12) = 0.25*sl*t4*zl
     DPSI(12) = 0.25*sl*t2*z1
     
	 phi (13) = 0.25*sl*tl*z2
     DXHI(13) = 0.25*s1*tl*z2
     DTHE(13) = 0.25*sl*t1*z2
     DPSI(13) = 0.25*sl*tl*z4
     
	 phi (14) = 0.25*sq*tl*z2
     DXHI(14) = 0.25*s3*tl*z2
     DTHE(14) = 0.25*sq*t1*z2
     DPSI(14) = 0.25*sq*tl*z4
     
	 phi ( 15) = 0.25*sq*tp*z2
     DXHI(15) = 0.25*s3*tp*z2
     DTHE(15) = 0.25*sq*t3*z2
     DPSI(15) = 0.25*sq*tp*z4
     
	 phi ( 16) = 0.25*sl*tp*z2
     DXHI(16) = 0.25*s1*tp*z2
     DTHE(16) = 0.25*sl*t3*z2
     DPSI(16) = 0.25*sl*tp*z4
     
	 phi ( 17) = 0.25*s2*tl*zp
     DXHI(17) = 0.25*s4*tl*zp
     DTHE(17) = 0.25*s2*t1*zp
     DPSI(17) = 0.25*s2*tl*z3
     
	 phi ( 18) = 0.25*sq*t2*zp
     DXHI(18) = 0.25*s3*t2*zp
     DTHE(18) = 0.25*sq*t4*zp
     DPSI(18) = 0.25*sq*t2*z3
     
	 phi ( 19) = 0.25*s2*tp*zp
     DXHI(19) = 0.25*s4*tp*zp
     DTHE(19) = 0.25*s2*t3*zp
     DPSI(19) = 0.25*s2*tp*z3
     
	 phi ( 20) = 0.25*sl*t2*zp
     DXHI(20) = 0.25*s1*t2*zp
     DTHE(20) = 0.25*sl*t4*zp
     DPSI(20) = 0.25*sl*t2*z3
     
	 phi ( 21) = 0.5*s2*t2*zl
     DXHI(21) = 0.5*s4*t2*zl
     DTHE(21) = 0.5*s2*t4*zl
     DPSI(21) = 0.5*s2*t2*z1
     
	 phi ( 22) = 0.5*s2*tl*z2
     DXHI(22) = 0.5*s4*tl*z2
     DTHE(22) = 0.5*s2*t1*z2
     DPSI(22) = 0.5*s2*tl*z4
     
	 phi (23) = 0.5*sq*t2*z2
     DXHI(23) = 0.5*s3*t2*z2
     DTHE(23) = 0.5*sq*t4*z2
     DPSI(23) = 0.5*sq*t2*z4
     
	 phi (  24) = 0.5*s2*tp*z2
     DXHI(24) = 0.5*s4*tp*z2
     DTHE(24) = 0.5*s2*t3*z2
     DPSI(24) = 0.5*s2*tp*z4
     
	 phi ( 25) = 0.5*sl*t2*z2
     DXHI(25) = 0.5*s1*t2*z2
     DTHE(25) = 0.5*sl*t4*z2
     DPSI(25) = 0.5*sl*t2*z4
     
	 phi (  26) = 0.5*s2*t2*zp
     DXHI(26) = 0.5*s4*t2*zp
     DTHE(26) = 0.5*s2*t4*zp
     DPSI(26) = 0.5*s2*t2*z3
     
	 phi (27) = s2*t2*z2
     DXHI(27) = s4*t2*z2
     DTHE(27) = s2*t4*z2
     DPSI(27) = s2*t2*z4
   
return
end


subroutine DETERM(AJACO, AJACOI, DETER)   
IMPLICIT none
double precision :: deter, AJACO(3,3),AJACOI(3,3)
!local
double precision :: t1,t2,t3,denom
 
     t1  = AJACO(2,2)*AJACO(3,3) - AJACO(3,2)*AJACO(2,3)
     t2  =-AJACO(2,1)*AJACO(3,3) + AJACO(3,1)*AJACO(2,3)
     t3  = AJACO(2,1)*AJACO(3,2) - AJACO(3,1)*AJACO(2,2)
     deter = AJACO(1,1)*t1 + AJACO(1,2)*t2 + AJACO(1,3)*t3
     if(deter .eq. 0.0) return
     
	 denom = 1.0/deter
     AJACOI(1,1) = t1*denom
     AJACOI(2,1) = t2*denom
     AJACOI(3,1) = t3*denom
     AJACOI(2,2) = ( AJACO(1,1)*AJACO(3,3) - AJACO(3,1)*AJACO(1,3))*denom
     AJACOI(3,2) = (-AJACO(1,1)*AJACO(3,2) + AJACO(1,2)*AJACO(3,1))*denom
     AJACOI(3,3) = ( AJACO(1,1)*AJACO(2,2) - AJACO(2,1)*AJACO(1,2))*denom
     AJACOI(1,2) = (-AJACO(1,2)*AJACO(3,3) + AJACO(3,2)*AJACO(1,3))*denom
     AJACOI(1,3) = ( AJACO(1,2)*AJACO(2,3) - AJACO(2,2)*AJACO(1,3))*denom
     AJACOI(2,3) = (-AJACO(1,1)*AJACO(2,3) + AJACO(2,1)*AJACO(1,3))*denom


end subroutine determ

