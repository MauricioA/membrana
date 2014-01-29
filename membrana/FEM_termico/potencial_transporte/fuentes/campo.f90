subroutine campo(nnodes,nelements,nodpel,solution,material,conect,coor_x,coor_y,coor_z,sigma1,sigma2,sigma3,sigma4,  &
   &        grad_x,grad_y,grad_z,gradxel_x,gradxel_y,gradxel_z)
implicit none
integer :: nnodes,nelements,nodpel,conect(nelements,nodpel),material(nelements)
double precision :: grad_x(nnodes),grad_y(nnodes),grad_z(nnodes) 
double precision :: gradxel_x(nelements),gradxel_y(nelements),gradxel_z(nelements)

double precision :: solution(nnodes),coor_x(nnodes),coor_y(nnodes),coor_z(nnodes),sigma1,sigma2 ,sigma3,sigma4
! local
INTEGER :: nope,NS(nodpel),jel,mat
DOUBLE PRECISION :: X(nodpel),Y(nodpel),Z(nodpel),sol(nodpel),sigma_el,Ex_el,ey_el,ez_el,ex(nodpel),ey(nodpel),ez(nodpel)


DOUBLE PREcIsION:: PHI(nodpel),DPHIX(nodpel),DPHIY(nodpel),DPHIZ(nodpel),AJACO(3,3),AJACOI(3,3),DXHI(nodpel), &
     DTHE(nodpel),DPSI(nodpel),XHI,THE,PSI,DETER,S11,S22,S33,CNST,denom
DOUBLE PREcIsION,allocatable :: GAUSSPT(:),GAUSSWT(:)

integer kk,jj,i,j,ii,K,NLE,pgaus,I2,ngaus
	 


if(nodpel==8) then
  allocate(gausspt(2),gausswt(2))
  ngaus=2
  GAUSSPT(1)=-0.57735027
  GAUSSPT(2)= 0.57735027
  GAUSSWT(1)= 1.0
  GAUSSWT(2)= 1.0 

elseif(nodpel==27) then
  allocate(gausspt(3),gausswt(3))
  ngaus=3
  GAUSSPT(1)=-0.774596669241483377035853079956
  GAUSSPT(2)= 0.0
  GAUSSPT(3)= 0.774596669241483377035853079956
  GAUSSWT(1)= 0.5555555556
  GAUSSWT(2)= 0.8888888889 
  GAUSSWT(3)=0.5555555556 

endif

grad_x=0.0
grad_y=0.0
grad_z=0.0
DO JEL=1,nelements

   mat=material(jel)        
   if(mat==4) then
       sigma_el = sigma1
   elseif(mat==2 .or. mat==3) then
       sigma_el = sigma2
   elseif(mat==1 ) then
       sigma_el = sigma3
   elseif(mat==0.or.mat==5 ) then
       sigma_el = sigma4
   elseif(mat==10 .or. mat==11 ) then
       sigma_el = 1.0
   endif
   
   
   DO I=1,nodpel
        ns(I)=conect(JEL,I)
	    j=NS(I)
        X(i)=coor_x(j)
        y(i)=coor_y(j)
        z(i)=coor_z(j)
        sol(i)=solution(j)
        Ex(i)=0.0
        Ey(i)=0.0
        Ez(i)=0.0
   ENDDO
        
   Ex_el=0
   Ey_el=0
   Ez_el=0

   pgaus=0
   DO KK=1,ngaus
     DO JJ=1,ngaus
       DO II=1,ngaus
         pgaus=pgaus+1

		 XHI = GAUSSPT(KK)
         THE = GAUSSPT(JJ)
         PSI = GAUSSPT(II)
      
      if(ngaus==2) then
        call funciones8(XHI,THE,PSI,PHI,DPHIX,DPHIY,DPHIZ,AJACO,AJACOI,DXHI,DTHE,DPSI ) 
     else         

        call funciones(XHI,THE,PSI,PHI,DPHIX,DPHIY,DPHIZ,AJACO,AJACOI,DXHI,DTHE,DPSI ) 
    endif
!   JACOBIANO Y SU INVERSA
        AJACO=0.0

          DO K=1,nodpel
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

          DO I=1,nodpel
            DPHIX(I)=AJACOI(1,1)*DXHI(I) + AJACOI(1,2)*DTHE(I) + AJACOI(1,3)*DPSI(I)
            DPHIY(I)=AJACOI(2,1)*DXHI(I) + AJACOI(2,2)*DTHE(I) + AJACOI(2,3)*DPSI(I)
            DPHIZ(I)=AJACOI(3,1)*DXHI(I) + AJACOI(3,2)*DTHE(I) + AJACOI(3,3)*DPSI(I)
          ENDDO

!          CNST = DETER * GAUSSWT(JJ)*GAUSSWT(KK)*GAUSSWT(II)
         ! gradientes
         DO I=1,nodpel
        
             Ex(i)=Ex(i)+ DPHIX(i)*sol(i)
             Ey(i)=Ey(i)+ DPHIY(i)*sol(i)
             Ez(i)=Ez(i)+ DPHIZ(i)*sol(i)
             Ex_el=Ex_el+ DPHIX(i)*sol(i)/real(nodpel)
             Ey_el=Ey_el+ DPHIY(i)*sol(i)/real(nodpel)
             Ez_el=Ez_el+ DPHIZ(i)*sol(i)/real(nodpel)

         ENDDO
      ENDDO
      ENDDO
    ENDDO
      
    do i=1,nodpel
        j=ns(i)
        if(i<9) denom=8.0
        if(i>=9 .and. i<=20) denom=4.0
        if(i>=21 .and. i<=26) denom=2.0
        if(i==27) denom=1.0


        grad_x(j) =  grad_x(j) - Ex(i)/denom
        grad_y(j) =  grad_y(j) - Ey(i)/denom
        grad_z(j) =  grad_z(j) - Ez(i)/denom
    enddo
   gradxel_x(jel) = -Ex_el
   gradxel_y(jel) = -Ey_el
   gradxel_z(jel) = -Ez_el
    
enddo


grad_x=0
grad_y=0
grad_z=0
 do jel=1,nelements
     do i=1,nodpel
        
        grad_x(conect(jel,i)) =  grad_x(conect(jel,i)) + gradxel_x(jel)/real(nodpel) 
        grad_y(conect(jel,i)) =  grad_y(conect(jel,i)) + gradxel_y(jel)/real(nodpel) 
        grad_z(conect(jel,i)) =  grad_z(conect(jel,i)) + gradxel_z(jel)/real(nodpel) 
    enddo
enddo

deallocate(gausspt,gausswt)

end subroutine campo




subroutine campo2d(nnodes,nelements,nodpel,solution,material,conect,coor_x,coor_y,sigmaext,sigmaint,sigmamem,grad_x,grad_y,  &
      &    gradxel_x,gradxel_y)


implicit none
integer :: nnodes,nelements,nodpel,conect(nelements,nodpel),material(nelements)
double precision :: grad_x(nnodes),grad_y(nnodes)
double precision :: gradxel_x(nelements),gradxel_y(nelements)

double precision :: solution(nnodes),coor_x(nnodes),coor_y(nnodes),sigmaext,sigmaint,sigmamem
! local
INTEGER :: nope,NS(nodpel),jel,mat
DOUBLE PRECISION :: X(nodpel),Y(nodpel),sol(nodpel),sigma_el,Ex_el,ey_el,ez_el,ex(nodpel),ey(nodpel)
DOUBLE PREcIsION:: B(3),C(3),DETER,A,RMED
DOUBLE PREcIsION:: PI=3.14159
DOUBLE PRECISION,allocatable :: gausspt(:),gausswt(:)
DOUBLE PRECISION,allocatable :: phi(:,:),dphi(:,:,:),gxcod(:,:),PHIdX(:,:,:),cteI(:) 

integer kk,jj,i,j,ii,K,NLE,kgaus,I2,ngaus,ndimension,ndime,ndi
dOUBLE PREcIsION:: t,s,sm,tm,sq,tp,AJACO(2,2),AJACOI(2,2),DPHIX(nodpel),DPHIY(nodpel)


     ngaus=2
     ndimension=2

      allocate(gausspt(Ngaus),gausswt(Ngaus),phi(2*Ngaus,nodpel),dphi(ndimension,2*Ngaus,nodpel),gxcod(ndimension,2*Ngaus),PHIdX(ndimension,2*Ngaus,nodpel),cteI(Ngaus*2)) 


grad_x=0.0
grad_y=0.0
gradxel_x=0.0
gradxel_y=0.0

DO JEL=1,nelements

   mat=material(jel)        
   if(mat==3) then
       sigma_el = sigmaint
   elseif(mat==2) then
       sigma_el = sigmamem
   elseif(mat==1 ) then
       sigma_el = sigmaext
   endif
   
   
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
   
       gradxel_x(jel)=-(B(1)*sol(1)+B(2)*sol(2)+B(3)*sol(3))/DETER
       gradxel_y(jel)=-(C(1)*sol(1)+C(2)*sol(2)+C(3)*sol(3))/DETER
   
    !   gradxel_x(jel)=sigma_el* gradxel_x(jel)*(-1.)
    !   gradxel_y(jel)=sigma_el* gradxel_y(jel)*(-1.)   

   
       do i=1,nodpel
            j=ns(i)
            grad_x(j) =  grad_x(j) + gradxel_x(i)/real(nodpel)
            grad_y(j) =  grad_y(j) + gradxel_y(i)/real(nodpel)
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
      
   ! do i=1,nodpel
   !     j=ns(i)

   !     grad_x(j) =  grad_x(j) - Ex(i)*0.25
   !     grad_y(j) =  grad_y(j) - Ey(i)*0.25
   ! enddo
    gradxel_x(jel) = -Ex_el*0.25
    gradxel_y(jel) = -Ey_el*0.25
    



 endif  





enddo



end subroutine campo2d




subroutine campo1d(nnodes,nelements,nodpel,solution,material,conect,coor_x,sigmaext,sigmaint,sigmamem,grad_x,gradxel_x)
implicit none
integer :: nnodes,nelements,nodpel,conect(nelements,nodpel),material(nelements)
double precision :: grad_x(nnodes)
double precision :: gradxel_x(nelements)

double precision :: solution(nnodes),coor_x(nnodes),sigmaext,sigmaint,sigmamem
! local
INTEGER :: nope,NS(nodpel),jel,mat
DOUBLE PRECISION :: X(nodpel),sol(nodpel),sigma_el,Ex_el,ex(nodpel)
DOUBLE PREcIsION:: PI=3.14159
DOUBLE PRECISION,allocatable :: gauspt(:),gauswt(:)
DOUBLE PRECISION,allocatable :: phi(:,:),dphi(:,:,:),gxcod(:,:),PHIdX(:,:,:),cteI(:) 

integer kk,jj,i,j,ii,K,NLE,kgaus,I2,ngaus,ndimension,ndime,ndi
dOUBLE PREcIsION:: t,s,sm,tm,sq,tp,AJACO,AJACOI,DPHIX(nodpel),DPHIY(nodpel),xhi


     ngaus=nodpel-1
     ndimension=1

 allocate(gauspt(Ngaus),gauswt(Ngaus),phi(Ngaus,nope),dphi(ndimension,Ngaus,nope),gxcod(ndimension,Ngaus),PHIdX(ndimension,Ngaus,nope),cteI(Ngaus)) 

	      gauspt(1) = -.906179845938664
          gauspt(2) = -.538469310105683
          gauspt(3) = 0.0
          gauspt(4) = .538469310105683
          gauspt(5) = .906179845938664
		      
          gauswt(1) = .236926885056189
          gauswt(2) = .478628670499366
          gauswt(3) = .568888888888889
          gauswt(4) = .478628670499366
          gauswt(5) = .236926885056189



grad_x=0.0
gradxel_x=0.0

DO JEL=1,nelements

   mat=material(jel)        
   if(mat==3) then
       sigma_el = sigmaint
   elseif(mat==2) then
       sigma_el = sigmamem
   elseif(mat==1 ) then
       sigma_el = sigmaext
   endif
   
   
   DO I=1,nodpel
        ns(I)=conect(JEL,I)
	    j=NS(I)
        X(i)=coor_x(j)
        sol(i)=solution(j)
        Ex(i)=0.0
   ENDDO
        
   Ex_el=0
         
	do kgaus=1,Ngaus
             XHI = GAUSPT(kgaus)
	         PHI(kgaus,1) = (-625./768.)*( XHI**5-XHI**4-(10.0/25.)*XHI**3+       &
                    &     (10.0/25.)*XHI**2+(9./625.0)*XHI-(9./625.))
             PHI(kgaus,2) = (3125./768.)*( XHI**5-3./5.*XHI**4-(26.0/25.)*XHI**3+  &
                    &     (78.0/125.)*XHI**2+(1./25.)*XHI-3./125.)
             PHI(kgaus,3) = (-3125./384.)*( XHI**5-1./5*XHI**4-(34.0/25.)*XHI**3+  &
                    &     (34.0/125.)*XHI**2+(9./25)*XHI-9./125)
             PHI(kgaus,4)= (3125./384.)*( XHI**5+1./5*XHI**4-(34.0/25.)*XHI**3-    &
                    &     (34.0/125.)*XHI**2+(9./25.)*XHI+9./125.)                      
             PHI(kgaus,5) =(-3125./768.)*( XHI**5+(3./5.)*XHI**4-                  &
                    &    (26.0/25.)*XHI**3 -(78.0/125.)*XHI**2+(1./25.)*XHI+3./125.)
             PHI(kgaus,6) = (625./768.)*( XHI**5+XHI**4-(10.0/25.)*XHI**3-         &
                    &      (10.0/25.)*XHI**2+(9./625.0)*XHI+(9./625.))


              DpHI(1,kgaus,1)  = (-625./768.)*( 5*XHI**4-4*XHI**3-(30.0/25.)*XHI**2+  &
                      &               (20.0/25.)*XHI+(9./625.0))
              DpHI(1,kgaus,2)  = (3125./768.)*( 5*XHI**4-12./5.*XHI**3-               &
                      &                (78.0/25.)*XHI**2+(156.0/125)*XHI+(1./25.))
              DpHI(1,kgaus,3)  = (-3125./384.)*( 5*XHI**4- 4./5. *XHI**3-             &
                      &                (102.0/25.)*XHI**2+ (68.0/125.)*XHI+(9./25.))
              DpHI(1,kgaus,4)  = (3125./384.)*( 5*XHI**4+4./5*XHI**3-                  &  
                      &                (102.0/25.)*XHI**2-(68.0/125.)*XHI+(9./25.))
              DpHI(1,kgaus,5)  =(-3125./768.)*( 5*XHI**4+(12./5.)*XHI**3-             &
                      &              (78.0/25.)*XHI**2-(156.0/125)*XHI+(1./25))
              DpHI(1,kgaus,6)  = (625./768.)*( 5*XHI**4+4*XHI**3-(30.0/25)*XHI**2-     &
                      &              (20.0/25.)*XHI+(9./625.0))
	 
             do ndime=1,ndimension
                 gxcod(ndime,kgaus)=0.0
                 do ii=1,nodpel
                     gxcod(ndime,kgaus) = gxcod(ndime,kgaus) + x(ii)*phi(kgaus,ii)
                 enddo
             enddo


             AJACO=0.0
             DO ii=1,nodpel
                 AJACO = AJACO + DPHI(1,kgaus,ii)*x(ii)
             ENDDO
          

             IF(AJACO.EQ.0.0) THEN
                  WRITE(6,*) ' ELEMENTO: con determinante cero en el PUNTO (shapes)',Kgaus
                  print*,'Dionisio va a cerrarse!'
                  STOP ' '
             ENDIF
    
             AJACOI = 1.0/AJACO


             DO ii=1,nodpel
                PHIdx(1,kgaus,ii)=AJACOI*DPHI(1,kgaus,ii)
             ENDDO

             
	         cteI(kgaus)= AJACO * GAUSWT(kgaus)*2*PI*GXCOD(1,kgaus) 
          

             DO I=1,nodpel
        
               Ex(i)=Ex(i)+ PHIdx(1,kgaus,i)*sol(i)
             
               Ex_el=Ex_el+ PHIdx(1,kgaus,i)*sol(i)
             ENDDO


	 enddo


      
    do i=1,nodpel
        j=ns(i)

        grad_x(j) =  grad_x(j) - Ex(i)*0.25
    enddo
    gradxel_x(jel) = -Ex_el*0.25
    

enddo



end subroutine campo1d
