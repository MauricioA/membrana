subroutine grad_concentra(nnodes,nelements,nodpel,solution,material,conect,coor_x,coor_y,grad_x,grad_y)


implicit none
integer :: nnodes,nelements,nodpel,conect(nelements,nodpel),material(nelements)
double precision :: grad_x(nnodes),grad_y(nnodes)

double precision :: solution(nnodes),coor_x(nnodes),coor_y(nnodes)
! local
INTEGER :: nope,NS(nodpel),jel,mat
DOUBLE PRECISION :: X(nodpel),Y(nodpel),sol(nodpel),sigma_el,Ex_el,ey_el,ez_el,ex(nodpel),ey(nodpel)
DOUBLE PREcIsION:: B(3),C(3),DETER,A,RMED
DOUBLE PREcIsION:: PI=3.14159
DOUBLE PRECISION,allocatable :: gausspt(:),gausswt(:)
DOUBLE PRECISION,allocatable :: phi(:,:),dphi(:,:,:),gxcod(:,:),PHIdX(:,:,:),cteI(:) 

integer kk,jj,i,j,ii,K,NLE,kgaus,I2,ngaus,ndimension,ndime,ndi
dOUBLE PREcIsION:: t,s,sm,tm,sq,tp,AJACO(2,2),AJACOI(2,2),DPHIX(nodpel),DPHIY(nodpel)
dOUBLE PREcIsION:: denomJ(nnodes), gradxel_x(nelements) , gradxel_y(nelements) 


     ngaus=2
     ndimension=2

      allocate(gausspt(Ngaus),gausswt(Ngaus),phi(2*Ngaus,nodpel),dphi(ndimension,2*Ngaus,nodpel),gxcod(ndimension,2*Ngaus),PHIdX(ndimension,2*Ngaus,nodpel),cteI(Ngaus*2)) 


grad_x=0.0
grad_y=0.0


DO JEL=1,nelements

      
   
   DO I=1,nodpel
        ns(I)=conect(JEL,I)
	    j=NS(I)
        X(i)=coor_x(j)
        y(i)=coor_y(j)
        sol(i)=solution(j)
        Ex(i)=0.0
        Ey(i)=0.0
   ENDDO
        
   Ex_el=0
   Ey_el=0


  if(nodpel==3) then    
    !  EVALUACION DEL GRADIENTE EN LOS ELEMENTOS
       B(1)=Y(2)-Y(3)
       B(2)=Y(3)-Y(1)
       B(3)=Y(1)-Y(2)
       C(1)=X(3)-X(2)
       C(2)=X(1)-X(3)
       C(3)=X(2)-X(1)
   
       DETER=X(2)*Y(3)+X(3)*Y(1)+X(1)*Y(2)-X(2)*Y(1)-X(3)*Y(2)-X(1)*Y(3)
   
       
    !   gradxel_x(jel)=sigma_el* gradxel_x(jel)*(-1.)
    !   gradxel_y(jel)=sigma_el* gradxel_y(jel)*(-1.)   

   
       do i=1,nodpel
            j=ns(i)
       !     grad_x(j) =  grad_x(j) + gradxel_x(i)/real(nodpel)
       !     grad_y(j) =  grad_y(j) + gradxel_y(i)/real(nodpel)
       enddo
    
 else


      GAUSSPT(1)=-0.57735027
      GAUSSPT(2)= 0.57735027
      GAUSSWT(1)= 1.0
      GAUSSWT(2)= 1.0 
 
       kgaus=0
   DO KK=1,ngaus
     DO JJ=1,ngaus
         kgaus=kgaus+1

		 t = GAUSSPT(KK)
         s = GAUSSPT(JJ)
      

             
              sm = 0.5*(1.0-s)
              tm = 0.5*(1.0-t)
              sq = 0.5*(1.0+s)
              tp = 0.5*(1.0+t)
              
              phi(kgaus,1)    = sm*tm
              dphi(1,kgaus,1) =-0.5*tm
              dphi(2,kgaus,1) =-0.5*sm
              phi(kgaus,2)     = sq*tm
              dphi(1,kgaus, 2) = 0.5*tm
              dphi(2,kgaus, 2) =-0.5*sq
              phi(kgaus,3)     = sq*tp
              dphi(1,kgaus, 3) = 0.5*tp
              dphi(2,kgaus, 3) = 0.5*sq
              phi(kgaus,4)     = sm*tp
              dphi(1,kgaus, 4) =-0.5*tp
              dphi(2,kgaus, 4) = 0.5*sm
              

!   JACOBIANO Y SU INVERSA
        AJACO=0.0

          DO K=1,nodpel
              AJACO(1,1)=AJACO(1,1)+dphi(1,kgaus,k)*X(K)
              AJACO(1,2)=AJACO(1,2)+dphi(1,kgaus,k)*Y(K)
              AJACO(2,1)=AJACO(2,1)+dphi(2,kgaus,k)*X(K)
              AJACO(2,2)=AJACO(2,2)+dphi(2,kgaus,k)*Y(K)
          ENDDO

          DETER= AJACO(1,1)*AJACO(2,2)-AJACO(1,2)*AJACO(2,1)
          AJACOI(1,1)=AJACO(2,2)/DETER
          AJACOI(2,2)=AJACO(1,1)/DETER
          AJACOI(1,2)=-AJACO(1,2)/DETER
          AJACOI(2,1)=-AJACO(2,1)/DETER

          DO I=1,nodpel
            DPHIX(I)=AJACOI(1,1)*dphi(1,kgaus,i) + AJACOI(1,2)*dphi(2,kgaus,i)
            DPHIY(I)=AJACOI(2,1)*dphi(1,kgaus,i) + AJACOI(2,2)*dphi(2,kgaus,i) 
          ENDDO

         DO I=1,nodpel
        
             Ex(i)=Ex(i)+ DPHIX(i)*sol(i)
             Ey(i)=Ey(i)+ DPHIY(i)*sol(i)
             
             Ex_el=Ex_el+ DPHIX(i)*sol(i)
             Ey_el=Ey_el+ DPHIY(i)*sol(i)
         ENDDO


      ENDDO

    ENDDO
      
    do i=1,nodpel
  !      j=ns(i)

   !     grad_x(j) =  grad_x(j) - Ex(i)*0.25
   !     grad_y(j) =  grad_y(j) - Ey(i)*0.25
             
     !   Ex_el=Ex_el+ DPHIX(i)*sol(i)
     !   Ey_el=Ey_el+ DPHIY(i)*sol(i)
    ENDDO



 endif  

    gradxel_x(jel) = Ex_el/real(nodpel)
    gradxel_y(jel) = Ey_el/real(nodpel)
    
enddo



grad_x=0.0
grad_y=0.0
denomJ=0
     
  do kk=1,nelements
     DO I=1,nodpel
        ns(I)=conect(kk,I)
	    j=NS(I)
       
        denomJ(j)=denomJ(j)+1
      
    enddo

enddo

do kk=1,nelements
     DO I=1,nodpel
        ns(I)=conect(kk,I)
	    j=NS(I)

        
        grad_x(j)=grad_x(j) + gradxel_x(kk)
        grad_y(j)=grad_y(j) + gradxel_y(kk)
    enddo

enddo





end subroutine grad_concentra

