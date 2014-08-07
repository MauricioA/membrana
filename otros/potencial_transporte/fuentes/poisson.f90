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




allocate(solucion_ant(nnodes))

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
        
        !if(nmode==3) then
        !    sigma_el=1.0
        !endif


! jel = elemento actual
! nodpel = 3
! x = arreglo con 3 coordenadas x, 1 por cada nodo del elem
! y = arreglo con 3 coordenadas y, 1 por cada nodo del elem
! sol = arreglo con 3 soluciones anteriores, una por cada nodo
        qe=0.0
        DO I=1,nodpel
            ns(I)=conect(JEL,I)
	        j=NS(I)
            X(i)=coor_x(j)
            if(ndimension==2) then
               y(i)=coor_y(j)
            endif
            sol(i)=solucion_ant(j)
            
            !qe=qe+carga(j)/real(nodpel)

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
        elseif(nodpel==4) then
            CALL ARMADO4(JEL,X,Y,ns,nodpel,ESM,EF,sigma_el,qe,sol,gradxel_x(jel),gradxel_y(jel))
        elseif(nodpel==6) then
            CALL ARMADO1d(JEL,X,Y,ns,nodpel,ESM,EF,sigma_el,qe,sol,gradxel_x(jel),gradxel_y(jel))
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

     
    solucion=RHS
    iter=0
    err=1.0
    call CG(nnodes,IA,JA,AN,AD,RHS,solucion,toler,itermax,ITER,ERR)

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
    if(nmode==1) then
       WRITE(6,*) 'iteraciones internas del loop poisson ', nconta, error, ITER,err  
    endif


enddo ! end while

!call campo(nnodes,nelements,nodpel,solucion,material,conect,coor_x,coor_y,coor_z,sigma1,sigma2,sigma3,sigma4,grad_x,grad_y,grad_z,  &
!      &    gradxel_x,gradxel_y,gradxel_z)
if(ndimension==2) then
    call campo2d(nnodes,nelements,nodpel,solucion,material,conect,coor_x,coor_y,sigmaext,sigmaint,sigmamem,grad_x,grad_y,  &
      &    gradxel_x,gradxel_y)
elseif(ndimension==1) then
    call campo1d(nnodes,nelements,nodpel,solucion,material,conect,coor_x,sigmaext,sigmaint,sigmamem,grad_x,gradxel_x)

endif


if(nmode==1) then
   call salida_sol(solucion)
endif

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

