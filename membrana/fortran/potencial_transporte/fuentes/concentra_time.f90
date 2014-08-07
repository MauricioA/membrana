subroutine concentraH_time(time,deltat,ch_ant,coh_ant,cna_ant,ccl_ant)
use def_transpor
use def_solver
use def_variables
use def_constantes
implicit none
double precision :: time,deltat,ch_ant(nnodes),coh_ant(nnodes),cna_ant(nnodes),ccl_ant(nnodes)
!local
integer :: npasos,nlimite,inode,jnode,kk,ncase,i,j,ipoin,ii,jj,jj2,iaux,keje,ke
double precision :: epsil,error,x(nodpel),y(nodpel),ef(nodpel),esm(nodpel,nodpel),chdes(nodpel),cohdes(nodpel),ef_ant(nodpel),mas(nodpel)
double precision :: ymed,tmed,dx,qe,fundx,funqe,adiag,denom,burn,potvol,  &
   &                hconduc,funconduc,Pint,funpint

double precision :: th2,Acoef2,Acoef1,funrho,funcp,cp,rho,landa,rinvaina,flujomed,Tvai,f0,b0,dista,mu

double precision, allocatable :: ch_ini(:),coh_ini(:)
integer :: mat,NS(nodpel),iter
double precision :: chmed,cohmed,pot_med,Dh_el,err,hmed,current,current2

allocate(ch_ini(nnodes),coh_ini(nnodes))


th2=0.5

Acoef2 = (1.0-th2)*deltat
Acoef1 = th2*deltat

NCASE=1
epsil=0.01
npasos=0
error=1.0
nlimite=100

ncase=1

!do while(error>epsil .and. npasos<nlimite)

  an=0
  ad=0
  rhs=0
  
  do kk=1,nelements
       
     mat=material(kk)        
     ymed=0.0
	 Chmed=0.0
	 Cohmed=0.0
     pot_med=0.0
     qe=0.0

     do i=1,nodpel
         ns(I)=conect(kk,I)
	     j=NS(I)
         X(i)=coor_x(j)
         if(ndimension==2) then
            y(i)=coor_y(j)
            ymed=ymed+y(i)/real(nodpel)
         else

            ymed=ymed+x(i)/real(nodpel)

         endif
         mas(i)=masa(j)

         chdes(i)=ch_ant(j)
         chmed=chmed+ch_ant(j)/real(nodpel)
         cohmed=cohmed+coh_ant(j)/real(nodpel)
         pot_med=pot_med + solucion(j)/real(nodpel)
	 enddo

     if(ndimension==2) then
       hmed = y(nodpel)-y(1)
     else
       hmed = (x(nodpel)-x(1))/(real(nodpel))
     endif 

     if(nmode==2) then
         if(mat==1 ) then
               Dh_el = D_h*exp(-zh*Number_clave*pot_med)
         elseif(mat==2 ) then
               Dh_el =  D_h*exp(-zh*Number_clave*pot_med)*1e-3
         elseif(mat==3) then
               Dh_el = D_h*exp(-zh*Number_clave*pot_med)
         endif 
     
     
         landa = exp(-zh*Number_clave*pot_med)
     else
         if(mat==1 ) then
               Dh_el = D_h
         elseif(mat==2 ) then
               Dh_el =  D_h*1e-3
         elseif(mat==3) then
               Dh_el = D_h
         endif 
     
     
         landa =1.0

         if(ndimension==2) then
             mu = -Dh_el*Number_clave*zh*gradxel_y(kk)
         elseif(ndimension==1) then
             mu = -Dh_el*Number_clave*zh*gradxel_x(kk)
         endif

     endif
     
    ! if(ymed>=2.0 .and. ymed <= 98.0 ) then
    !    QE =  kwb*Ch2o - kwf* chmed* cohmed 
    ! endif
     
     
     
     call ARMADO_t(kk,X,Y,ns,nodpel,ESM,EF,Dh_el,qe,landa,mu,chdes,Acoef1,Acoef2,mas)

     do inode=1,nodpel
            ipoin=NS(inode)
            if(  vec_tierra(ipoin)/=-1 ) then

              ! calculo ch_catodo mediante ecuacion con I
               current = 0.0

              ! ch_catodo = (D_h*chmed + current/Faraday)/(mu+D_h/(hmed*0.5))
            

              adiag=ESM(inode,inode)
              do jnode=1,nodpel
                 ESM(inode,jnode)=0.0
                 EF (jnode)=EF(jnode)-ESM(jnode,inode)*ch_catodo
                 ESM(jnode,inode)=0.0
              end do
              ESM(inode,inode)= adiag
              EF(inode)       = adiag* ch_catodo
            end if
            if(  vec_poten(ipoin)/=-1 ) then


              ! calculo ch_anodo mediante ecuacion con I
               current = I_eq* ( exp(-Faraday*(solucion(ipoin)+E_eq)/(2*R_cte*T_cte)) - ch_ant(ipoin)/ch_inicial*exp(Faraday*(solucion(ipoin)+E_eq)/(2*R_cte*T_cte))) 
               current2 = I2_eq* ( ccl_ant(ipoin)/ccl_inicial * exp(-Faraday*(solucion(ipoin)+E2_eq)/(2*R_cte*T_cte)) - exp(Faraday*(solucion(ipoin)+E2_eq)/(2*R_cte*T_cte))) 

            !  ch_anodo = (D_h*chmed + (current+current2)/Faraday)/(D_h/(hmed)- mu)
              
            !  if(ipoin==1301) then
            !     write(1111,*)'anodo H: ', time,ipoin,ch_anodo,current+current2
            !  endif

              adiag=ESM(inode,inode)
              do jnode=1,nodpel
                 ESM(inode,jnode)=0.0
                 EF (jnode)=EF(jnode)-ESM(jnode,inode)*ch_anodo
                 ESM(jnode,inode)=0.0
              end do
              ESM(inode,inode)= adiag
              EF(inode)       = adiag* ch_anodo
            end if
        end do

! ENSAMBLO
       if(nmode==2) then
         DO II=1,nodpel
            RHS(NS(II))=RHS(NS(II)) + EF(II)
	        AD(NS(II)) = AD(NS(II)) + ESM(II,II)
            	     
	        DO JJ=1,nodpel


                JJ2=NS(JJ)+1-NS(II)
	            IF(JJ2.GT.0) THEN
		 
	                DO IAUX = 1,CX(NS(II))
         
		                KEJE = IA(NS(II))+IAUX-1  
	          	  	        
       		            IF( JA(KEJE) .EQ. NS(II) + JJ2 -1) THEN           
 		       
		                    AN( KEJE ) = AN(KEJE ) +  ESM(II,JJ)
                
		                ENDIF
          
		           ENDDO

	            ENDIF
	
	        ENDDO
         ENDDO	  
      elseif(nmode==3) then
         DO II=1,nodpel
            RHS(NS(II))=RHS(NS(II)) + EF(II)
	        AD(NS(II)) = AD(NS(II)) + ESM(II,II)
            	     
	        DO JJ=1,nodpel

!                JJ2=NS(JJ)+1-NS(II)
!	            IF(JJ2.NE.0) THEN
	            IF(JJ.NE.II) THEN
		 
	                DO IAUX = 1,CX(NS(II))
         
		                KEJE = IA(NS(II))+IAUX-1  
	          	  	        
       		            IF( JA(KEJE) .EQ. NS(jj)) THEN           
 		       
		                    AN( KEJE ) = AN(KEJE ) +  ESM(II,JJ)
                
		                ENDIF
          
		           ENDDO

	            ENDIF
	
	        ENDDO
         ENDDO	  


      endif

    ENDDO

     
    ch=RHS
    iter=0
    err=1.0
    if(nmode==2) then
       call CG(nnodes,IA,JA,AN,AD,RHS,ch,toler,itermax,ITER,ERR)
    else
       call BCG(nnodes,IA,JA,AN,AD,RHS,ch,toler_dos,itermax,ITER,ERR)
    endif

!  error=0.0 
!  denom=0.0
!  do ke=1,nnodes
!      error=error+ (ch(ke)-ch_ini(ke))**2
!	  denom = denom + ch(ke)**2
!      ch_ini(ke)=ch(ke)
!  enddo   

!  error = sqrt(error/denom)
!  npasos=npasos+1
  
  
!enddo ! fin while

!deallocate(ch_ini)


end subroutine concentraH_time


subroutine concentraOH_time(time,deltat,ch_ant,coh_ant,cna_ant,ccl_ant)
use def_transpor
use def_solver
use def_variables
use def_constantes

implicit none
double precision :: time,deltat,ch_ant(nnodes),coh_ant(nnodes),cna_ant(nnodes),ccl_ant(nnodes)
!local
integer :: npasos,nlimite,inode,jnode,kk,ncase,i,j,ipoin,ii,jj,jj2,iaux,keje,ke
double precision :: epsil,error,x(nodpel),y(nodpel),ef(nodpel),esm(nodpel,nodpel),chdes(nodpel),cohdes(nodpel),ef_ant(nodpel),mas(nodpel)
double precision :: ymed,tmed,dx,qe,fundx,funqe,adiag,denom,burn,potvol,  &
   &                hconduc,funconduc,Pint,funpint

double precision :: th2,Acoef2,Acoef1,funrho,funcp,cp,rho,landa,rinvaina,flujomed,Tvai,f0,b0,dista,mu
double precision, allocatable :: ch_ini(:),coh_ini(:)
integer :: mat,ns(nodpel),iter
double precision :: chmed,cohmed,Doh_el,err,pot_med,hmed,current

!allocate(coh_ini(nnodes))


th2=0.5

Acoef2 = (1.0-th2)*deltat
Acoef1 = th2*deltat


NCASE=1
epsil=0.01
npasos=0
error=1.0
nlimite=100

ncase=1

!do while(error>epsil .and. npasos<nlimite)

  an=0
  ad=0
  rhs=0


  do kk=1,nelements
       
     mat=material(kk)        
     ymed=0.0
	 Chmed=0.0
	 Cohmed=0.0
     pot_med=0.0
     qe=0.0

     do i=1,nodpel
         ns(I)=conect(kk,I)
	     j=NS(I)
         X(i)=coor_x(j)
         if(ndimension==2) then
            y(i)=coor_y(j)
            ymed=ymed+y(i)/real(nodpel)
         else

            ymed=ymed+x(i)/real(nodpel)

         endif
        
         mas(i)=masa(j)

         cohdes(i)=coh_ant(j)
         chmed=chmed+ch_ant(j)/real(nodpel)
         cohmed=cohmed+coh_ant(j)/real(nodpel)
         pot_med=pot_med + solucion(j)/real(nodpel)


	 enddo

          if(ndimension==2) then
       hmed = y(nodpel)-y(1)
     else
       hmed = x(nodpel)-x(1)
     endif 



     if(nmode==2) then

         if(mat==1 ) then
               Doh_el = D_oh*exp(-zoh*Number_clave*pot_med)
         elseif(mat==2 ) then
               Doh_el = D_oh*exp(-zoh*Number_clave*pot_med)*1e-3
         elseif(mat==3) then
               Doh_el = D_oh*exp(-zoh*Number_clave*pot_med)
         endif 
     
         landa=exp(-zoh*Number_clave*pot_med)
     else

         if(mat==1 ) then
               Doh_el = D_oh
         elseif(mat==2 ) then
               Doh_el = D_oh*1e-3
         elseif(mat==3) then
               Doh_el = D_oh
         endif 
     
         landa=1.0
          if(ndimension==2) then
             mu = -Doh_el*Number_clave*zoh*gradxel_y(kk)
         elseif(ndimension==1) then
             mu = -Doh_el*Number_clave*zoh*gradxel_x(kk)
         endif


     endif

    ! if(ymed>=2.0 .and. ymed <= 98.0 ) then
        QE =  kwb*Ch2o - kwf* chmed* cohmed 
     !   if(qe<0) qe=0
     !endif
     

     call ARMADO_t(kk,X,Y,ns,nodpel,ESM,EF,Doh_el,qe,landa,mu,cohdes,Acoef1,Acoef2,mas)

     do inode=1,nodpel
            ipoin=NS(inode)
            if(  vec_tierra(ipoin)/=-1 ) then
             
                ! calculo coh_catodo mediante ecuacion con I
               current = I_eq* ( exp(-Faraday*(solucion(ipoin)+E_eq)/(2*R_cte*T_cte)) - ch_ant(ipoin)/ch_inicial*exp(Faraday*(solucion(ipoin)+E_eq)/(2*R_cte*T_cte))) 

              !coh_catodo = (-D_oh*cohmed + current/Faraday)/(-mu - D_oh/(hmed*0.5))
               



              
             
              adiag=ESM(inode,inode)
              do jnode=1,nodpel
                 ESM(inode,jnode)=0.0
                 EF (jnode)=EF(jnode)-ESM(jnode,inode)*coh_catodo
                 ESM(jnode,inode)=0.0
              end do
              ESM(inode,inode)= adiag
              EF(inode)       = adiag* coh_catodo
            end if
            if(  vec_poten(ipoin)/=-1 ) then

              ! calculo coh_anodo mediante ecuacion con I
               current = 0.0

             !  coh_anodo = (D_oh*cohmed + current/Faraday)/(D_oh/(hmed*0.5)- mu)

              adiag=ESM(inode,inode)
              do jnode=1,nodpel
                 ESM(inode,jnode)=0.0
                 EF (jnode)=EF(jnode)-ESM(jnode,inode)*coh_anodo
                 ESM(jnode,inode)=0.0
              end do
              ESM(inode,inode)= adiag
              EF(inode)       = adiag* coh_anodo
            end if
        end do

! ENSAMBLO
       if(nmode==2) then
         DO II=1,nodpel
            RHS(NS(II))=RHS(NS(II)) + EF(II)
	        AD(NS(II)) = AD(NS(II)) + ESM(II,II)
            	     
	        DO JJ=1,nodpel


                JJ2=NS(JJ)+1-NS(II)
	            IF(JJ2.GT.0) THEN
		 
	                DO IAUX = 1,CX(NS(II))
         
		                KEJE = IA(NS(II))+IAUX-1  
	          	  	        
       		            IF( JA(KEJE) .EQ. NS(II) + JJ2 -1) THEN           
 		       
		                    AN( KEJE ) = AN(KEJE ) +  ESM(II,JJ)
                
		                ENDIF
          
		           ENDDO

	            ENDIF
	
	        ENDDO
         ENDDO	  
      elseif(nmode==3) then
         DO II=1,nodpel
            RHS(NS(II))=RHS(NS(II)) + EF(II)
	        AD(NS(II)) = AD(NS(II)) + ESM(II,II)
            	     
	        DO JJ=1,nodpel

!                JJ2=NS(JJ)+1-NS(II)
!	            IF(JJ2.NE.0) THEN
	            IF(JJ.NE.II) THEN
		 
	                DO IAUX = 1,CX(NS(II))
         
		                KEJE = IA(NS(II))+IAUX-1  
	          	  	        
       		            IF( JA(KEJE) .EQ. NS(jj)) THEN           
 		       
		                    AN( KEJE ) = AN(KEJE ) +  ESM(II,JJ)
                
		                ENDIF
          
		           ENDDO

	            ENDIF
	
	        ENDDO
         ENDDO	  


      endif
    

    ENDDO

     
    coh=RHS
    iter=0
    err=1.0
    if(nmode==2) then
        call CG(nnodes,IA,JA,AN,AD,RHS,coh,toler,itermax,ITER,ERR)
    else
        
        call BCG(nnodes,IA,JA,AN,AD,RHS,coh,toler_dos,itermax,ITER,ERR)

    endif

!  error=0.0 
!  denom=0.0
!  do ke=1,nnodes
!      error=error+ (coh(ke)-coh_ini(ke))**2
!	  denom = denom + coh(ke)**2
!      coh_ini(ke)=coh(ke)
!  enddo   

!  error = sqrt(error/denom)
!  npasos=npasos+1
  
  
!enddo ! fin while

!deallocate(coh_ini)

end subroutine concentraOH_time



subroutine concentraNa_time(time,deltat,ch_ant,coh_ant,cna_ant,ccl_ant)
use def_transpor
use def_solver
use def_variables
use def_constantes
implicit none
double precision :: time,deltat,ch_ant(nnodes),coh_ant(nnodes),cna_ant(nnodes),ccl_ant(nnodes)
!local
integer :: npasos,nlimite,inode,jnode,kk,ncase,i,j,ipoin,ii,jj,jj2,iaux,keje,ke
double precision :: epsil,error,x(nodpel),y(nodpel),ef(nodpel),esm(nodpel,nodpel),chdes(nodpel),cohdes(nodpel),ef_ant(nodpel),mas(nodpel)
double precision :: ymed,tmed,dx,qe,fundx,funqe,adiag,denom,burn,potvol,  &
   &                hconduc,funconduc,Pint,funpint

double precision :: th2,Acoef2,Acoef1,funrho,funcp,cp,rho,landa,rinvaina,flujomed,Tvai,f0,b0,dista,mu

double precision, allocatable :: ch_ini(:),coh_ini(:),cna_ini(:),ccl_ini(:)
integer :: mat,NS(nodpel),iter
double precision :: chmed,cohmed,cnamed,cclmed,pot_med,DNa_el,err,hmed,current

allocate(ch_ini(nnodes),coh_ini(nnodes),cna_ini(nnodes),ccl_ini(nnodes))


th2=0.5

Acoef2 = (1.0-th2)*deltat
Acoef1 = th2*deltat

NCASE=1
epsil=0.01
npasos=0
error=1.0
nlimite=100

ncase=1

!do while(error>epsil .and. npasos<nlimite)

  an=0
  ad=0
  rhs=0
  
  do kk=1,nelements
       
     mat=material(kk)        
     ymed=0.0
	 Chmed=0.0
	 Cohmed=0.0
	 Cnamed=0.0
	 Cclmed=0.0
     pot_med=0.0
     qe=0.0

     do i=1,nodpel
         ns(I)=conect(kk,I)
	     j=NS(I)
         X(i)=coor_x(j)
         if(ndimension==2) then
            y(i)=coor_y(j)
            ymed=ymed+y(i)/real(nodpel)
         else

            ymed=ymed+x(i)/real(nodpel)

         endif
         mas(i)=masa(j)

         chdes(i)=cna_ant(j)
         
         chmed=chmed+ch_ant(j)/real(nodpel)
         cohmed=cohmed+coh_ant(j)/real(nodpel)
         cnamed=cnamed+cna_ant(j)/real(nodpel)
         cclmed=cclmed+ccl_ant(j)/real(nodpel)
         pot_med=pot_med + solucion(j)/real(nodpel)
	 enddo

         if(ndimension==2) then
       hmed = y(nodpel)-y(1)
     else
       hmed = x(nodpel)-x(1)
     endif 


     if(nmode==2) then
         if(mat==1 ) then
               Dna_el = D_na*exp(-zh*Number_clave*pot_med)
         elseif(mat==2 ) then
               Dna_el =  D_na*exp(-zh*Number_clave*pot_med)*1e-3
         elseif(mat==3) then
               Dna_el = D_na*exp(-zh*Number_clave*pot_med)
         endif 
     
     
         landa = exp(-zh*Number_clave*pot_med)
     else
         if(mat==1 ) then
               Dna_el = D_na
         elseif(mat==2 ) then
               Dna_el =  D_na*1e-3
         elseif(mat==3) then
               Dna_el = D_na
         endif 
     
     
         landa =1.0

          if(ndimension==2) then
             mu = -Dna_el*Number_clave*zna*gradxel_y(kk)
         elseif(ndimension==1) then
             mu = -Dna_el*Number_clave*zna*gradxel_x(kk)
         endif

     endif
     
        QE = 0.0 
     
     
     
     call ARMADO_t(kk,X,Y,ns,nodpel,ESM,EF,Dna_el,qe,landa,mu,chdes,Acoef1,Acoef2,mas)

     do inode=1,nodpel
            ipoin=NS(inode)
            if(  vec_tierra(ipoin)/=-1 ) then

              ! calculo cna_catodo mediante ecuacion con I
               current = 0.0

              ! cna_catodo = -(D_na*cnamed + current/Faraday)/(mu-D_na/(hmed*0.5))
            

              adiag=ESM(inode,inode)
              do jnode=1,nodpel
                 ESM(inode,jnode)=0.0
                 EF (jnode)=EF(jnode)-ESM(jnode,inode)*cna_catodo
                 ESM(jnode,inode)=0.0
              end do
              ESM(inode,inode)= adiag
              EF(inode)       = adiag* cna_catodo
            end if
            if(  vec_poten(ipoin)/=-1 ) then


              ! calculo ch_anodo mediante ecuacion con I
               current = 0.0

            !   cna_anodo = (D_na*cnamed + current/Faraday)/(D_na/(hmed*0.5)+ mu)
              
            

              adiag=ESM(inode,inode)
              do jnode=1,nodpel
                 ESM(inode,jnode)=0.0
                 EF (jnode)=EF(jnode)-ESM(jnode,inode)*cna_anodo
                 ESM(jnode,inode)=0.0
              end do
              ESM(inode,inode)= adiag
              EF(inode)       = adiag* cna_anodo
            end if
        end do

! ENSAMBLO
       if(nmode==2) then
         DO II=1,nodpel
            RHS(NS(II))=RHS(NS(II)) + EF(II)
	        AD(NS(II)) = AD(NS(II)) + ESM(II,II)
            	     
	        DO JJ=1,nodpel


                JJ2=NS(JJ)+1-NS(II)
	            IF(JJ2.GT.0) THEN
		 
	                DO IAUX = 1,CX(NS(II))
         
		                KEJE = IA(NS(II))+IAUX-1  
	          	  	        
       		            IF( JA(KEJE) .EQ. NS(II) + JJ2 -1) THEN           
 		       
		                    AN( KEJE ) = AN(KEJE ) +  ESM(II,JJ)
                
		                ENDIF
          
		           ENDDO

	            ENDIF
	
	        ENDDO
         ENDDO	  
      elseif(nmode==3) then
         DO II=1,nodpel
            RHS(NS(II))=RHS(NS(II)) + EF(II)
	        AD(NS(II)) = AD(NS(II)) + ESM(II,II)
            	     
	        DO JJ=1,nodpel

!                JJ2=NS(JJ)+1-NS(II)
!	            IF(JJ2.NE.0) THEN
	            IF(JJ.NE.II) THEN
		 
	                DO IAUX = 1,CX(NS(II))
         
		                KEJE = IA(NS(II))+IAUX-1  
	          	  	        
       		            IF( JA(KEJE) .EQ. NS(jj)) THEN           
 		       
		                    AN( KEJE ) = AN(KEJE ) +  ESM(II,JJ)
                
		                ENDIF
          
		           ENDDO

	            ENDIF
	
	        ENDDO
         ENDDO	  


      endif

    ENDDO

     
    cna=RHS
    iter=0
    err=1.0
    if(nmode==2) then
       call CG(nnodes,IA,JA,AN,AD,RHS,cna,toler,itermax,ITER,ERR)
    else
       call BCG(nnodes,IA,JA,AN,AD,RHS,cna,toler_dos,itermax,ITER,ERR)
    endif

!  error=0.0 
!  denom=0.0
!  do ke=1,nnodes
!      error=error+ (ch(ke)-ch_ini(ke))**2
!	  denom = denom + ch(ke)**2
!      ch_ini(ke)=ch(ke)
!  enddo   

!  error = sqrt(error/denom)
!  npasos=npasos+1
  
  
!enddo ! fin while

!deallocate(ch_ini)


end subroutine concentraNa_time



subroutine concentracl_time(time,deltat,ch_ant,coh_ant,cna_ant,ccl_ant)
use def_transpor
use def_solver
use def_variables
use def_constantes
implicit none
double precision :: time,deltat,ch_ant(nnodes),coh_ant(nnodes),cna_ant(nnodes),ccl_ant(nnodes)
!local
integer :: npasos,nlimite,inode,jnode,kk,ncase,i,j,ipoin,ii,jj,jj2,iaux,keje,ke
double precision :: epsil,error,x(nodpel),y(nodpel),ef(nodpel),esm(nodpel,nodpel),chdes(nodpel),cohdes(nodpel),ef_ant(nodpel),mas(nodpel)
double precision :: ymed,tmed,dx,qe,fundx,funqe,adiag,denom,burn,potvol,  &
   &                hconduc,funconduc,Pint,funpint

double precision :: th2,Acoef2,Acoef1,funrho,funcp,cp,rho,landa,rinvaina,flujomed,Tvai,f0,b0,dista,mu

double precision, allocatable :: ch_ini(:),coh_ini(:),cna_ini(:),ccl_ini(:)
integer :: mat,NS(nodpel),iter
double precision :: chmed,cohmed,cnamed,cclmed,pot_med,Dcl_el,err,hmed,current,current2

allocate(ch_ini(nnodes),coh_ini(nnodes),cna_ini(nnodes),ccl_ini(nnodes))


th2=0.5

Acoef2 = (1.0-th2)*deltat
Acoef1 = th2*deltat

NCASE=1
epsil=0.01
npasos=0
error=1.0
nlimite=100

ncase=1

!do while(error>epsil .and. npasos<nlimite)

  an=0
  ad=0
  rhs=0
  
  do kk=1,nelements
       
     mat=material(kk)        
     ymed=0.0
	 Chmed=0.0
	 Cohmed=0.0
	 Cnamed=0.0
	 Cclmed=0.0
     pot_med=0.0
     qe=0.0

     do i=1,nodpel
         ns(I)=conect(kk,I)
	     j=NS(I)
         X(i)=coor_x(j)
          if(ndimension==2) then
            y(i)=coor_y(j)
            ymed=ymed+y(i)/real(nodpel)
         else

            ymed=ymed+x(i)/real(nodpel)

         endif
         mas(i)=masa(j)

         chdes(i)=ccl_ant(j)
         
         chmed=chmed+ch_ant(j)/real(nodpel)
         cohmed=cohmed+coh_ant(j)/real(nodpel)
         cnamed=cnamed+cna_ant(j)/real(nodpel)
         cclmed=cclmed+ccl_ant(j)/real(nodpel)
         pot_med=pot_med + solucion(j)/real(nodpel)
	 enddo

          if(ndimension==2) then
       hmed = y(nodpel)-y(1)
     else
       hmed = x(nodpel)-x(1)
     endif 


     if(nmode==2) then
         if(mat==1 ) then
               Dcl_el = D_cl*exp(-zcl*Number_clave*pot_med)
         elseif(mat==2 ) then
               Dcl_el =  D_cl*exp(-zcl*Number_clave*pot_med)*1e-3
         elseif(mat==3) then
               Dcl_el = D_cl*exp(-zcl*Number_clave*pot_med)
         endif 
     
     
         landa = exp(-zcl*Number_clave*pot_med)
     else
         if(mat==1 ) then
               Dcl_el = D_cl
         elseif(mat==2 ) then
               Dcl_el =  D_cl*1e-3
         elseif(mat==3) then
               Dcl_el = D_cl
         endif 
     
     
         landa =1.0
         if(ndimension==2) then
             mu = -Dcl_el*Number_clave*zcl*gradxel_y(kk)
         elseif(ndimension==1) then
             mu = -Dcl_el*Number_clave*zcl*gradxel_x(kk)
         endif


     endif
     
     QE = 0.0 
     
          
     call ARMADO_t(kk,X,Y,ns,nodpel,ESM,EF,Dcl_el,qe,landa,mu,chdes,Acoef1,Acoef2,mas)

     do inode=1,nodpel
            ipoin=NS(inode)
            if(  vec_tierra(ipoin)/=-1 ) then

              ! calculo ccl_catodo mediante ecuacion con I
               current = 0.0
               
            !   ccl_catodo = -(D_cl*cclmed + current/Faraday)/(mu-D_cl/(hmed*0.5))
            

              adiag=ESM(inode,inode)
              do jnode=1,nodpel
                 ESM(inode,jnode)=0.0
                 EF (jnode)=EF(jnode)-ESM(jnode,inode)*ccl_catodo
                 ESM(jnode,inode)=0.0
              end do
              ESM(inode,inode)= adiag
              EF(inode)       = adiag* ccl_catodo
            end if
            if(  vec_poten(ipoin)/=-1 ) then


              ! calculo ccl_anodo mediante ecuacion con I
               current = I_eq* ( exp(-Faraday*(solucion(ipoin)+E_eq)/(2*R_cte*T_cte)) - ch_ant(ipoin)/ch_inicial*exp(Faraday*(solucion(ipoin)+E_eq)/(2*R_cte*T_cte))) 
               current2 = I2_eq* ( ccl_ant(ipoin)/ccl_inicial * exp(-Faraday*(solucion(ipoin)+E2_eq)/(2*R_cte*T_cte)) - exp(Faraday*(solucion(ipoin)+E2_eq)/(2*R_cte*T_cte))) 

           !   ccl_anodo = (D_cl*cclmed + (current+current2)/Faraday)/(D_cl/(hmed*0.5)- mu)

              
            
              adiag=ESM(inode,inode)
              do jnode=1,nodpel
                 ESM(inode,jnode)=0.0
                 EF (jnode)=EF(jnode)-ESM(jnode,inode)*ccl_anodo
                 ESM(jnode,inode)=0.0
              end do
              ESM(inode,inode)= adiag
              EF(inode)       = adiag* ccl_anodo
            end if
        end do

! ENSAMBLO
       if(nmode==2) then
         DO II=1,nodpel
            RHS(NS(II))=RHS(NS(II)) + EF(II)
	        AD(NS(II)) = AD(NS(II)) + ESM(II,II)
            	     
	        DO JJ=1,nodpel


                JJ2=NS(JJ)+1-NS(II)
	            IF(JJ2.GT.0) THEN
		 
	                DO IAUX = 1,CX(NS(II))
         
		                KEJE = IA(NS(II))+IAUX-1  
	          	  	        
       		            IF( JA(KEJE) .EQ. NS(II) + JJ2 -1) THEN           
 		       
		                    AN( KEJE ) = AN(KEJE ) +  ESM(II,JJ)
                
		                ENDIF
          
		           ENDDO

	            ENDIF
	
	        ENDDO
         ENDDO	  
      elseif(nmode==3) then
         DO II=1,nodpel
            RHS(NS(II))=RHS(NS(II)) + EF(II)
	        AD(NS(II)) = AD(NS(II)) + ESM(II,II)
            	     
	        DO JJ=1,nodpel

!                JJ2=NS(JJ)+1-NS(II)
!	            IF(JJ2.NE.0) THEN
	            IF(JJ.NE.II) THEN
		 
	                DO IAUX = 1,CX(NS(II))
         
		                KEJE = IA(NS(II))+IAUX-1  
	          	  	        
       		            IF( JA(KEJE) .EQ. NS(jj)) THEN           
 		       
		                    AN( KEJE ) = AN(KEJE ) +  ESM(II,JJ)
                
		                ENDIF
          
		           ENDDO

	            ENDIF
	
	        ENDDO
         ENDDO	  


      endif

    ENDDO

     
    ccl=RHS
    iter=0
    err=1.0
    if(nmode==2) then
       call CG(nnodes,IA,JA,AN,AD,RHS,ccl,toler,itermax,ITER,ERR)
    else
       call BCG(nnodes,IA,JA,AN,AD,RHS,ccl,toler_dos,itermax,ITER,ERR)
    endif

!  error=0.0 
!  denom=0.0
!  do ke=1,nnodes
!      error=error+ (ch(ke)-ch_ini(ke))**2
!	  denom = denom + ch(ke)**2
!      ch_ini(ke)=ch(ke)
!  enddo   

!  error = sqrt(error/denom)
!  npasos=npasos+1
  
  
!enddo ! fin while

!deallocate(ch_ini)


end subroutine concentracl_time
