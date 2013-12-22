subroutine poisson()
use def_variables
use def_constantes
use def_solver
implicit none
double precision :: x(nodpel),y(nodpel),z(nodpel),esm(nodpel,nodpel),ef(nodpel),adiag,sigma_el,err,qe
double precision, allocatable  :: esm_tot(:,:)

integer :: ns(nodpel)
integer  NPASO,KK,JEL,I,II,ncase,inode,ipoin,jnode,jj,iaux,keje,jj2,mat,j,iter
double precision, allocatable :: solucion_ant(:)
double precision :: error,epsil,denom,numer,sol(nodpel),funsigma1,campoxl,tierra
integer :: nconta,NCOTA



if(nopcion==1) then     
  write(6,*) '      voy a control!'
  call control()
else
   write(6,*) '      leo el sistema!'
   call lee_sistema(archi_sistema,nnodes)
endif

allocate(solucion_ant(nnodes),esm_tot(nnodes,nnodes))

esm_tot=0.0

tierra=0.0
error=1.0
nconta=0
epsil = 1e-3
NCOTA=10
solucion_ant=0.0

do while(error>epsil .and. nconta< NCOTA )
    
    nconta=nconta+1
    an=0
    ad=0
    rhs=0
        
    DO JEL=1,nelements

        mat=material(jel)        
        if(mat==1 ) then
            sigma_el = sigmaext
        elseif(mat==2 ) then
            sigma_el = sigmamem
        elseif(mat==3) then
            sigma_el = sigmaint
        endif
   
        qe=0.0
        DO I=1,nodpel
            ns(I)=conect(JEL,I)
	        j=NS(I)
            X(i)=coor_x(j)
            y(i)=coor_y(j)
!            z(i)=coor_z(j)
            sol(i)=solucion_ant(j)
        ENDDO
        
        gradxel_x(jel)=0.0
        gradxel_y(jel)=0.0
!        gradxel_z(jel)=0.0

        if(nodpel==27) then
            CALL ARMADO(NCASE,JEL,X,Y,Z,ns,nodpel,ESM,EF,sigma_el,qe)
        elseif(nodpel==8) then
            CALL ARMADO8(NCASE,JEL,X,Y,Z,ns,nodpel,ESM,EF,sigma_el,qe,sol,gradxel_x(jel),gradxel_y(jel),gradxel_z(jel))
        elseif(nodpel==3) then
            CALL ARMADO3(JEL,X,Y,ns,nodpel,ESM,EF,sigma_el,qe,sol,gradxel_x(jel),gradxel_y(jel))
        endif
! INTRODUZCO EN LAS MATRICES LAS CONDICIONES DE CONTORNO

        do inode=1,nodpel
            ipoin=NS(inode)
            if(  vec_tierra(ipoin)/=-1 ) then
              adiag=ESM(inode,inode)
              do jnode=1,nodpel
                 ESM(inode,jnode)=0.0
                 EF (jnode)=EF(jnode)-ESM(jnode,inode)*tierra
                 ESM(jnode,inode)=0.0
              end do
              ESM(inode,inode)= adiag
              EF(inode)       = adiag* tierra
            end if
            if(  vec_poten(ipoin)/=-1 ) then
              adiag=ESM(inode,inode)
              do jnode=1,nodpel
                 ESM(inode,jnode)=0.0
                 EF (jnode)=EF(jnode)-ESM(jnode,inode)* potencial
                 ESM(jnode,inode)=0.0
              end do
              ESM(inode,inode)= adiag
              EF(inode)       = adiag* potencial
            end if
        end do


! ENSAMBLO
        DO II=1,nodpel
            RHS(NS(II))=RHS(NS(II)) + EF(II)
	        AD(NS(II)) = AD(NS(II)) + ESM(II,II)
            	     
	        DO JJ=1,nodpel

                esm_tot(NS(II),NS(jj))=esm_tot(NS(II),NS(jj)) + ESM(II,jj)

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
    

    ENDDO

    do ii=1,nnodes
    do jj=1,nnodes
         write(111,*) 'MR ',II,jj,esm_tot(ii,jj)
         enddo
         enddo
     
    solucion=RHS
    iter=0
    err=1.0
    call CG(nnodes,IA,JA,AN,AD,RHS,solucion,toler,itermax,ITER,ERR)
       do jj=1,nnodes
         write(111,*) 'Sol ',jj,solucion(jj)
        
         
         enddo

    if(problema=='CEREBRO') then
       error=0.0
       denom=0
       numer=0.0
       do kk=1,nnodes
          numer=numer + (solucion(kk)-solucion_ant(kk))*(solucion(kk)-solucion_ant(kk))
          denom = denom +  solucion(kk)*solucion(kk)
          solucion_ant(kk)=solucion(kk)
       enddo
       error=dsqrt(numer/denom)
    else
       error=epsil*0.5
    endif

    WRITE(unit_cont,*) 'iteraciones internas ', nconta, error, ITER,err  
    WRITE(6,*) 'iteraciones internas del loop poisson ', nconta, error, ITER,err  


enddo ! end while

!call campo(nnodes,nelements,nodpel,solucion,material,conect,coor_x,coor_y,coor_z,sigma1,sigma2,sigma3,sigma4,grad_x,grad_y,grad_z,  &
!      &    gradxel_x,gradxel_y,gradxel_z)
call campo2d(nnodes,nelements,nodpel,solucion,material,conect,coor_x,coor_y,sigmaext,sigmaint,sigmamem,grad_x,grad_y,  &
      &    gradxel_x,gradxel_y)

call salida_sol(solucion)


deallocate(solucion_ant)

end subroutine poisson	
      


double precision function funsigma1(mat,E)
implicit none
integer :: mat
double precision :: E

if(mat==4 .or. mat==2 .or. mat==3) then
   if(E>= 46.0 .and. E<=70.0) then
      funsigma1 = 1 + 2.5  
   else
      funsigma1 = 1
   endif
else
   funsigma1 = 1
endif

end function funsigma1

