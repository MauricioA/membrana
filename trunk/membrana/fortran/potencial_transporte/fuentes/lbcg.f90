
SUBROUTINE CG(NP,IA,JA,A_SPA,AD,B,X,TOL,ITMAX,ITER,ERR)
      implicit none 
	DOUBLE PRECISION A_SPA(*),AD(*),X(*),B(*)
      INTEGER IA(*),JA(*)
      DOUBLE PRECISION bnrm, Znrm, bknum, bkden, XANT, akden,bk,ak,zm1nrm,dxnrm,xnrm,alfa,b_norma,r_norma,aknum
	DOUBLE PRECISION TOL,EPS,ERR,FUNNORM
	INTEGER ITMAX,ITER,NP,j,kk

      DOUBLE PRECISION, ALLOCATABLE ::Z(:),P(:),PP(:),R(:),ZZ(:),RR(:)

	ALLOCATE (P(NP),PP(NP),ZZ(NP),Z(NP),R(NP),RR(NP)) 

      
	CALL MATXVECSIM(IA,JA,A_SPA,AD,X,R,NP)

      b_norma = 0.0
	DO KK=1,NP
        R(KK)  = B(KK) - R(KK)
        b_norma =b_norma  + b(kk)*b(kk) 
	ENDDO
      b_norma=sqrt(b_norma) 
      

	DO KK=1,NP
	  Z(KK)=R(KK)/AD(KK)
        rr(kk)=z(kk)
	ENDDO

      ZNRM=FUNNORM(NP,Z)
      
      do while(err.gt.tol .and. iter.lt.itmax)

        iter=iter+1

	  bknum=0.d0
  	  do j=1,NP
          bknum=bknum + Z(j)*R(j)
        enddo
	  CALL MATXVECSIM(IA,JA,A_SPA,AD,rr,Z,NP)
        
	  bkden=0.0
	  do j=1,NP
          bkden=bkden + rr(j)*Z(j)
        enddo
	   
        alfa= bknum/bkden    	  

        do kk=1,np
          x(kk)=x(kk)+alfa*rr(kk)          
	    r(kk)= r(kk) - alfa * z(kk)
	  enddo

     	  DO KK=1,NP
	    Z(KK)=R(KK)/AD(KK)
        ENDDO

    	  akden=bknum
	  aknum=0.0
        do  j=1,NP
          aknum = aknum+ z(j)*r(j)
        enddo

        ak=aknum/akden
          
        DO J=1,NP
          rr(j) = z(j)  + ak*rr(j)
        ENDDO

        r_norma = 0.0
	  DO KK=1,NP
          r_norma =r_norma  + r(kk)*r(kk)    
	  ENDDO
        r_norma=sqrt(r_norma) 

        err=r_norma/b_norma
        !write(2,*) iter,err,r_norma,ak,alfa

      ENDDO


	DEALLOCATE (P,PP,ZZ,Z,R,RR) 
	
	end subroutine cg
      
   DOUBLE PRECISION FUNCTION FUNNORM(N,Y)
   implicit none
   DOUBLE PRECISION:: Y(*)
   integer:: N
   !local
   integer :: isamax,i
   DOUBLE PRECISION:: sum2
   
        !ISAMAX = 1
        !DO  I=1,N
        !  IF(ABS(Y(I)).GT.ABS(Y(ISAMAX))) ISAMAX=I
        !ENDDO


	  !FUNNORM=ABS(Y(ISAMAX))

	  sum2=0.0

!$OMP PARALLEL DO 
!$omp& shared(Y), private(i), reduction(+: sum2) 
        DO  I=1,N
          sum2=sum2 + Y(i)*Y(i)
        END DO
!$OMP END PARALLEL DO

	  FUNNORM=sqrt(sum2)


	RETURN
	END


subroutine BCG(NP,IA,JA,A_SPA,AD,B,X,TOL,ITMAX,ITER,ERR)
implicit none 
DOUBLE PRECISION A_SPA(*),AD(*),X(*),B(*),err,tol
INTEGER IA(*),JA(*)
integer :: np,iter,itmax
!local
DOUBLE PRECISION ::bnrm, Znrm, bknum, bkden, XANT, EPS,akden,bk,ak,zm1nrm,dxnrm,xnrm,FUNNORM
DOUBLE PRECISION, ALLOCATABLE ::Z(:),P(:),PP(:),R(:),ZZ(:),RR(:)
INTEGER j,kk



aLLOCATE (P(NP),PP(NP),ZZ(NP),Z(NP),R(NP),RR(NP)) 

 EPS= 1.d-14

call MATXVEC(IA,JA,A_SPA,AD,X,R,NP)

	DO KK=1,NP
      R(KK)  = B(KK) - R(KK)
	  RR(KK) = R(KK)   
	ENDDO

      CALL MATXVEC(IA,JA,A_SPA,AD,R,RR,NP)


 	DO KK=1,NP
	  Z(KK)=B(KK)/AD(KK)
    ENDDO

      bnrm=FUNNORM(NP,Z)

	DO KK=1,NP
	  Z(KK)=R(KK)/AD(KK)
    ENDDO

      ZNRM=FUNNORM(NP,Z)

      do while(err.gt.tol .and. iter.lt.itmax)

        iter=iter+1
      
 	  DO KK=1,NP
	    ZZ(KK)=RR(KK)/AD(KK)
        ENDDO
   	  
	  bknum=0.d0
  	  do j=1,NP
          bknum=bknum + Z(j)*RR(j)
        enddo

        if(iter.eq.1) then
          do  j=1,NP
             P (j)=Z(j)
             PP(j)=ZZ(j)
          enddo
        else
          bk=bknum/bkden
          do  j=1,NP
            P(j) = bk*P(j)+Z(j)
            PP(j)= bk*PP(j)+ZZ(j)
          enddo
        endif

        bkden=bknum

	  CALL MATXVEC(IA,JA,A_SPA,AD,P,Z,NP)
      
    	  akden=0.d0
        do  j=1,NP
          akden = akden+ z(j)*pp(j)
        enddo

        ak=bknum/akden
          
        CALL MATXVEC(IA,JA,A_SPA,AD,PP,ZZ,NP)

        DO J=1,NP
          x(j) = x(j)  + ak*p(j)
          r(j) = r(j)  - ak*z(j)
          RR(j)= RR(j) - ak*zz(j)
        ENDDO

     	  DO KK=1,NP
	    Z(KK)=R(KK)/AD(KK)
        ENDDO

        zm1nrm=znrm
        znrm=FUNNORM(NP,Z)
        
	  if(abs(zm1nrm-znrm).gt.EPS*znrm) then
          
		 dxnrm=abs(ak)*FUNNORM(NP,P)
         err=znrm/abs(zm1nrm-znrm)*dxnrm
                  
      	 xnrm=FUNNORM(NP,X)
          
	     if(err.le.0.5d0*xnrm) then
             err=err/xnrm
         else
             err=znrm/bnrm
         endif

	  else
        
	    err=znrm/bnrm
       
	  endif

    ENDDO


dEALLOCATE (P,PP,ZZ,Z,R,RR) 
	

end subroutine BCG

      

SUBROUTINE LIN_BCG(NP,IA,JA,A_SPA,AD,B,X,unit_cont)
use def_constantes
integer :: np,unit_cont
double precision :: a_spa(*),AD(np),X(np),B(np)
integer :: IA(np+1),JA(*)
!local
DOUBLE PRECISION ::bnrm, Znrm, bknum, bkden, XANT, EPS,akden,bk,ak,err,zm1nrm,dxnrm,xnrm,FUNNORM
DOUBLE PRECISION, ALLOCATABLE ::Z(:),P(:),PP(:),R(:),ZZ(:),RR(:)



aLLOCATE (P(NP),PP(NP),ZZ(NP),Z(NP),R(NP),RR(NP)) 

iter=0
EPS= 1.d-14
err = 1.0
      
call MATXVECSIM(IA,JA,A_SPA,AD,X,R,NP)

	DO KK=1,NP
      R(KK)  = B(KK) - R(KK)
	  RR(KK) = R(KK)   
	ENDDO

      CALL MATXVECSIM(IA,JA,A_SPA,AD,R,RR,NP)


 	DO KK=1,NP
	  Z(KK)=B(KK)/AD(KK)
      if(abs(Z(kk))>1.0e-9) then
         write(6,*) kk,zz(kk),b(kk)/ad(kk)
      endif
    ENDDO

      bnrm=FUNNORM(NP,Z)

	DO KK=1,NP
	  Z(KK)=R(KK)/AD(KK)
      ENDDO

      ZNRM=FUNNORM(NP,Z)

      do while(err.gt.toler .and. iter.lt.itmax)

        iter=iter+1
      
 	  DO KK=1,NP
	    ZZ(KK)=RR(KK)/AD(KK)
        ENDDO
   	  
	  bknum=0.d0
  	  do j=1,NP
          bknum=bknum + Z(j)*RR(j)
        enddo

        if(iter.eq.1) then
          do  j=1,NP
             P (j)=Z(j)
             PP(j)=ZZ(j)
          enddo
        else
          bk=bknum/bkden
          do  j=1,NP
            P(j) = bk*P(j)+Z(j)
            PP(j)= bk*PP(j)+ZZ(j)
          enddo
        endif

        bkden=bknum

	  CALL MATXVECSIM(IA,JA,A_SPA,AD,P,Z,NP)
      
    	  akden=0.d0
        do  j=1,NP
          akden = akden+ z(j)*pp(j)
        enddo

        ak=bknum/akden
          
        CALL MATXVECSIM(IA,JA,A_SPA,AD,PP,ZZ,NP)

        DO J=1,NP
          x(j) = x(j)  + ak*p(j)
          r(j) = r(j)  - ak*z(j)
          RR(j)= RR(j) - ak*zz(j)
        ENDDO

     	  DO KK=1,NP
	    Z(KK)=R(KK)/AD(KK)
        ENDDO

        zm1nrm=znrm
        znrm=FUNNORM(NP,Z)
        
	  if(abs(zm1nrm-znrm).gt.EPS*znrm) then
          
		 dxnrm=abs(ak)*FUNNORM(NP,P)
         err=znrm/abs(zm1nrm-znrm)*dxnrm
                  
      	 xnrm=FUNNORM(NP,X)
          
	     if(err.le.0.5d0*xnrm) then
             err=err/xnrm
         else
             err=znrm/bnrm
         endif

	  else
        
	    err=znrm/bnrm
       
	  endif

    ENDDO

    WRITE(unit_cont,*) 'iteraciones internas ', ITER,err  

dEALLOCATE (P,PP,ZZ,Z,R,RR) 
	
end subroutine LIN_BCG

!DOUBLE PRECISION FUNCTION FUNNORM(N,Y) 
!integer :: n
!DOUBLE PRECISION :: Y(N)
! local
!integer :: isamax,i

!  isamax=1
!  do  i=1,N
!      if(abs(Y(i)).gt.abs(Y(isamax))) isamax=i
!  enddo
        
!  FUNNORM=abs(Y(isamax))

!end function FUNNORM
	
SUBROUTINE MATXVECSIM(IA,JA,AN,AD,B,C,Np)
implicit none
integer :: np
INTEGER :: IA(np+1),JA(*)
DOUBLE PRECISION :: B(np),AN(*),AD(np),C(np)
! local
integer :: k,i,iaf,iai,j

DO K=1,np
     C(K)= AD(K)*B(K)
enddo
      
DO I=1,np
     IAI = IA(I)
     IAF = IA(I+1)- 1


     IF(IAF.GE.IAI) THEN
          DO J=IAI,IAF
            C(I) = C(I) + AN(J)*B(JA(J))
            C(JA(J))=   C(JA(J)) + AN(J)*B(I)
          ENDDO
      ENDIF
enddo
     
end subroutine MATXVECSIM
      

SUBROUTINE MATXVEC(IA,JA,AN,AD,B,C,Np)  ! no simetrica o completa
implicit none
integer :: np
INTEGER :: IA(np+1),JA(*)
DOUBLE PRECISION :: B(np),AN(*),AD(np),C(np)
! local
integer :: k,i,iaf,iai,j

DO K=1,np
     C(K)= AD(K)*B(K)
enddo
      
DO I=1,np
     IAI = IA(I)
     IAF = IA(I+1)- 1

     DO J=IAI,IAF
         C(I) = C(I) + AN(J)*B(JA(J))
     ENDDO

enddo
     
end subroutine MATXVEC
         