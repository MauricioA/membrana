subroutine transporte()
use def_variables
use def_constantes
use def_solver
use def_transpor
implicit none
! local
double precision,allocatable :: ch_ant(:),phHaux(:)
double precision,allocatable :: coh_ant(:),phOHaux(:)
double precision,allocatable :: ccl_ant(:),cna_ant(:)
integer :: npaso_kk,nsale,i,ki,jj,kk,j,ii,cadena(10)
double precision :: tcero,deltat,tt,numer, denom,error,rsa,xmed,ymed


allocate(ch(nnodes),coh(nnodes),masa(nnodes))
allocate(cna(nnodes),ccl(nnodes),phOHaux(nnodes),phHaux(nnodes))
allocate(ch_ant(nnodes),coh_ant(nnodes))
allocate(cna_ant(nnodes),ccl_ant(nnodes))

allocate(grad_ccl_x(nnodes),grad_ccl_y(nnodes))
allocate(grad_cna_x(nnodes),grad_cna_y(nnodes))
allocate(grad_coh_x(nnodes),grad_coh_y(nnodes))
allocate(grad_ch_x(nnodes),grad_ch_y(nnodes))

grad_ccl_x=0.0
grad_ccl_y=0.0
grad_cna_x=0.0
grad_cna_y=0.0

grad_coh_x=0.0
grad_coh_y=0.0

grad_ch_x=0.0
grad_ch_y=0.0


write(unit_histo,*) nnodes
write(unit_ph,*) nnodes
write(unit_camp,*) nelements
       
call masadiag2d()

! inicializo mejor
ch = ch_inicial
coh = coh_inicial


!! ***CAMBIO***   LOs escenarios

do kk=1,nelements
   
   do jj=1,nodpel
      j=conect(kk,jj)
      cna(j) = cna_inicial
      ccl(j) = ccl_inicial

    !  if(material(kk)>2) then   ! escenario de cero concetracion dentro
    !     ch(j) = 0.0   
    !     coh(j) =0.0
    !     cna(j) = 0.0
    !     ccl(j) = 0.0
    !  endif

    !  if(material(kk)>2) then   ! escenario de concetracion menor
      !   ch(j) =ch_inicial*0.5  
      !   coh(j) =coh_inicial*0.5
      !   cna(j) = 0.0
      !   ccl(j) = 0.0
    !  endif

      if(material(kk)>2) then   ! escenario de concetracion igual  se espera acumulacion en los bordes
         cna(j) = cna_inicial*0.5
         ccl(j) = ccl_inicial*0.5
      endif


   enddo
enddo


ch_ant=ch
coh_ant=coh
cna_ant=cna
ccl_ant=ccl

cadena=0

do jj=1,nnodes
   phHaux(jj)=-log10((ch(jj)+1e-18)*1e+15/6.02E23)
   phOHaux(jj)=-log10((coh(jj)+1e-18)*1e+15/6.02E23)
enddo


if(nmode==2) then
    Number_clave= Faraday/(R_cte*T_cte)*0.6
    Ch_anodo=Ch_anodo*exp(zh*Number_clave*potencial)
    Coh_anodo=Coh_anodo*exp(zoh*Number_clave*potencial)
else
    Number_clave= Faraday/(R_cte*T_cte)
endif


npaso_kk=0
nsale=100

tcero=0.020
deltat=1e-7 ! 0.2e-5    !! ***CAMBIO*** bajo el dt
rsa = 0.50

tt=0.0

do while( tt<tcero)

   tt=tt+deltat
   npaso_kk=npaso_kk+1

   call calculo_carga()
   call poisson()


   if(nmode==2) then
       do ki=1,nnodes
          ch(ki) = ch(ki)*exp(zh*Number_clave*solucion(ki))
          coh(ki)= coh(ki)*exp(zoh*Number_clave*solucion(ki))
          ch_ant(ki) = ch_ant(ki)*exp(zh*Number_clave*solucion(ki))
          coh_ant(ki)= coh_ant(ki)*exp(zoh*Number_clave*solucion(ki))
       enddo
   endif
  
   
   call concentraH_time(tt,deltat,ch_ant,coh_ant,cna_ant,ccl_ant)
   call concentraOH_time(tt,deltat,ch_ant,coh_ant,cna_ant,ccl_ant)
   call concentraNa_time(tt,deltat,ch_ant,coh_ant,cna_ant,ccl_ant)
   call concentraCl_time(tt,deltat,ch_ant,coh_ant,cna_ant,ccl_ant)
   

   if(nmode==2) then
      do ki=1,nnodes
         ch(ki) = ch(ki)*exp(-zh*Number_clave*solucion(ki))
         coh(ki)= coh(ki)*exp(-zoh*Number_clave*solucion(ki))
         ch_ant(ki) = ch_ant(ki)*exp(-zh*Number_clave*solucion(ki))
         coh_ant(ki)= coh_ant(ki)*exp(-zoh*Number_clave*solucion(ki))
         phHaux(ki)=-log10((ch(ki)+1e-18)*1e+15/6.02E23)
         phOHaux(ki)=-log10((coh(ki)+1e-18)*1e+15/6.02E23)
      enddo
   else
      
      do ki=1,nnodes
        
         if(abs(ch(ki))<1.0e-8) then
             ch(ki)  =0.0
             ch_ant(ki)=ch(ki)
         endif
         if(abs(coh(ki))<1.0e-8) then
             coh(ki)=0.0
             coh_ant(ki)=0.0
         endif
         if(cna(ki)<1.0e-8)  then
             cna(ki)=0.0
             cna_ant(ki)=0.0
         endif

         if(ccl(ki)<1.0e-8) then 
             ccl(ki)=0.0 
             ccl_ant(ki)=0.0 
         endif
      
         phHaux(ki)=-log10((ch(ki)+ 1e-18)*1e+15/6.02E23)
         if(phHaux(ki)<=0.1)  phHaux(ki)=0.1
         phOHaux(ki)=-log10((coh(ki)+1e-18)*1e+15/6.02E23)
         if(phOHaux(ki)<=0.1) phOHaux(ki)=0.1
      enddo

   endif


   numer=0.0
   denom=0.0
   
   rresist=0.0

   do jj=1,nnodes
        numer=numer + (ch(jj)-ch_ant(jj))*(ch(jj)-ch_ant(jj))+ (coh(jj)-coh_ant(jj))*(coh(jj)-coh_ant(jj))  + (cna(jj)-cna_ant(jj))*(cna(jj)-cna_ant(jj))+ (ccl(jj)-ccl_ant(jj))*(ccl(jj)-ccl_ant(jj))
        denom = denom +  ch(jj)*ch(jj) + coh(jj)*coh(jj) +  cna(jj)*cna(jj) + ccl(jj)*ccl(jj)
        ch(jj) = rsa * ch(jj) + (1-rsa)*ch_ant(jj)
        coh(jj) = rsa * coh(jj) + (1-rsa)*coh_ant(jj)
        cna(jj) = rsa * cna(jj) + (1-rsa)*cna_ant(jj)
        ccl(jj) = rsa * ccl(jj) + (1-rsa)*ccl_ant(jj)
        
        ch_ant(jj)=ch(jj)
        coh_ant(jj)=coh(jj)
        cna_ant(jj)=cna(jj)
        ccl_ant(jj)=ccl(jj)

        rresist = rresist + ch(jj) + coh(jj) +  cna(jj) +  ccl(jj)

   enddo

   rresist  = rresist/nnodes

   rresist = rresist*124.2/1000.0 *1e12/6.02e+23 * (1.57e-4)
   ccurrent = potencial*rresist /(3.14159*50*50)  ! divido por el area en micrones

   write(unit_celda,*) tt,potencial,rresist,ccurrent,currentH,currentCl

   error=dsqrt(numer/denom)

   if(error>=1e+3) then
       write(6,*) 'paso: ',tt,error
       stop ' '
   endif

   if(mod(npaso_kk,nsale)==0) then
      
       write(6,*) 'paso: ',tt,error
   
       write(unit_histo,*) 'paso: ',tt
       write(unit_ph,*) 'paso: ',tt
       write(unit_camp2,*) 'paso: ',tt
       do jj=1,nnodes 
          
          write(unit_histo,'(i5,7e15.5)') jj,coor_x(jj),coor_y(jj),solucion(jj),ch(jj),coh(jj),cna(jj),ccl(jj)
          write(unit_ph,'(i5,5e15.5)') jj,coor_x(jj),coor_y(jj),solucion(jj),phHaux(jj),phOHaux(jj)
       enddo

       do jj=1,nelements
          xmed=0
          ymed=0
          do ii=1,nodpel
             xmed=xmed + coor_x(conect(jj,ii))/nodpel
             ymed=ymed + coor_y(conect(jj,ii))/nodpel
          enddo

          write(unit_camp2,'(4e15.5)') xmed,ymed,gradxel_x(jj),gradxel_y(jj)
       enddo



       ! reparo paso de tiempo
       !if(npaso_kk > 1000 .and. cadena(1)==0) then
       !     deltat=deltat*2
       !     cadena(1)=1
       !elseif(npaso_kk > 2000 .and. cadena(2)==0) then
       !     deltat=deltat*2
       !     cadena(2)=1
       !elseif(npaso_kk > 4000 .and. cadena(3)==0) then
       !     deltat=deltat*2
       !     cadena(3)=1
       !elseif(npaso_kk > 10000 .and. cadena(4)==0) then
       !     deltat=deltat*2
       !     cadena(4)=1
       !endif
   endif


enddo


call salida_concentra(ch,coh,cna,ccl,phHaux,phOHaux)

call salida_sol(solucion)


endsubroutine transporte



subroutine masadiag2d()
use def_solver
use def_variables
use def_constantes
implicit none
double precision x(nodpel),y(nodpel),gpdet,gpvol

integer ns(nodpel),ipoin,inode,jnode,ielem,kk,j,k,pdime,pnode,mdime,igaus,ilocs,jlocs
  !    calculamos la matriz de masa diagonal usando una regla cerrada de integracion

double precision ::  weigc(nodpel),posgl(4),weigl(4),B(nodpel),C(nodpel),deter,rmed,AJACO(2,2),deriv(2,nodpel,nodpel)
DOUBLE PREcIsION:: PI=3.14159



        

  masa=0.0
  
  do ielem = 1,nelements

      rmed=0.0   
      do inode = 1,nodpel
         ns(inode)= conect(ielem,inode)
         x(inode) = coor_x(ns(inode))
         y(inode) = coor_y(ns(inode))
         rmed=rmed+x(inode)/real(nodpel)
      end do
     
     
     if(nodpel==3) then
          B(1)=Y(2)-Y(3)
          B(2)=Y(3)-Y(1)
          B(3)=Y(1)-Y(2)
          C(1)=X(3)-X(2)
          C(2)=X(1)-X(3)
          C(3)=X(2)-X(1)

          DETER=X(2)*Y(3)+X(3)*Y(1)+X(1)*Y(2)-X(2)*Y(1)-X(3)*Y(2)-X(1)*Y(3)
          RMED=(X(1)+X(2)+X(3))/3.
      
      
          do inode=1,nodpel
          
             gpvol=  DETER*PI*RMED  
             masa(ns(inode))=masa(ns(inode))+gpvol

          enddo
    
    
    
     elseif(nodpel==4) then
      
      
        call armotodo(nodpel,x,y,deriv,weigc)

        do inode=1,nodpel
             
             AJACO(1,1) = 0.0
             AJACO(1,2) = 0.0
             AJACO(2,1) = 0.0
             AJACO(2,2) = 0.0
             do jnode = 1,nodpel
                AJACO(1,1) = AJACO(1,1) + x(jnode) * deriv(1,jnode,inode)
                AJACO(1,2) = AJACO(1,2) + x(jnode) * deriv(2,jnode,inode)
                AJACO(2,1) = AJACO(2,1) + y(jnode) * deriv(1,jnode,inode)
                AJACO(2,2) = AJACO(2,2) + y(jnode) * deriv(2,jnode,inode)
             end do
 
             gpdet = AJACO(1,1) * AJACO(2,2) - AJACO(2,1) * AJACO(1,2)
 
             gpvol=  weigc(inode)*gpdet
             
             masa(ns(inode))=masa(ns(inode))+gpvol*2*PI*rmed

      enddo
        
     elseif(nodpel==6) then
        
         do inode=1,nodpel
         
            masa(ns(inode))= 0.25

         enddo
              
     endif

   enddo

  
     !
     ! Loop over nodes to control zero-volume points
     !
     do ipoin=1,nnodes
       if(masa(ipoin)<1e-12) then
          write(6,*) 'nodo  ',ipoin,' posee matrix de masa cero  ',masa(ipoin)    
		  stop ' '
        end if
    ! debugdebugdebugdebugdebugdebugdebugdebugdebugdebugdebug
    !    write(12,*) ' *masa(',ipoin,')= ', masa(ipoin)
    ! debugdebugdebugdebugdebugdebugdebugdebugdebugdebugdebug

     enddo 


    


end subroutine masadiag2d
 


subroutine armotodo(nope,x,y,deriv,weigc)
implicit none
integer :: nope
double precision x(nope),y(nope),weigc(nope),deriv(2,nope,nope),shapf(nope)
double precision posgc(2,nope),posgl(nope),weigl(nope),s,t,st
integer inoga(nope),pnode,pdime,mdime,nlocs,igaus,ilocs,jlocs,klocs,inode


 !
  ! Element shape function and derivatives SHAPC,DERIC,HESLC,WEIGC 
  ! for a close rule
  !- For each element type, using a closed integration rule:
  !      WEIGC(nnode)
  !      SHAPC(nnode,nnode)
  !      DERIC(ndime,nnode,nnode)
  !      HESLC(ntens,nnode,nnode)
 
 if(nope==4) then
     inoga(1)= 1
     inoga(2)= 4
     inoga(3)= 2
     inoga(4)= 3
     
     pdime=2
     mdime=2
	 
     nlocs=2 
     posgl(1)=-1.0
     posgl(2)= 1.0
     weigl(1)= 1.0
     weigl(2)= 1.0
 
  
     igaus=0
     do ilocs=1,nlocs
        do jlocs=1,nlocs
              igaus=igaus+1
              weigc(  inoga(igaus))=weigl(ilocs)*weigl(jlocs)
              posgc(1,inoga(igaus))=posgl(ilocs)
              posgc(2,inoga(igaus))=posgl(jlocs)
        end do
     end do
     
      

     do inode=1,nope
        
         s=posgc(1,inode)
         t=posgc(2,inode)
         st=s*t                                           
      !   shapf(1)=(1.0-t-s+st)*0.25                     !  4         3
      !   shapf(2)=(1.0-t+s-st)*0.25                     !
      !   shapf(3)=(1.0+t+s+st)*0.25                     !      
      !   shapf(4)=(1.0+t-s-st)*0.25                     !
         
         deriv(1,1,inode)=(-1.0+t)*0.25                       !  1         2
         deriv(1,2,inode)=(+1.0-t)*0.25
         deriv(1,3,inode)=(+1.0+t)*0.25
         deriv(1,4,inode)=(-1.0-t)*0.25
         deriv(2,1,inode)=(-1.0+s)*0.25
         deriv(2,2,inode)=(-1.0-s)*0.25
         deriv(2,3,inode)=(+1.0+s)*0.25
         deriv(2,4,inode)=(+1.0-s)*0.25
    	 
	 end do

 elseif(nope==6) then
    
     
	 
     posgl(1)=-1.0
     posgl(2)= -2.0/3.0
     posgl(3)= -1.0/3.0
     posgl(4)= 1/3.
     posgl(5)= 2/3.
     posgl(6)= 1.0


     weigl(1)= 1.0
     weigl(2)= 1.0
     weigl(3)= 1.0
     weigl(4)= 1.0
     weigl(5)= 1.0
     weigl(6)= 1.0
 
  
     igaus=0
     do ilocs=1,nlocs
        do jlocs=1,nlocs
              igaus=igaus+1
              weigc(  inoga(igaus))=weigl(ilocs)*weigl(jlocs)
              posgc(1,inoga(igaus))=posgl(ilocs)
              posgc(2,inoga(igaus))=posgl(jlocs)
        end do
     end do
     
      

     do inode=1,nope
        
         s=posgc(1,inode)
         
         deriv(1,1,inode)=(-1.0+t)*0.25                       
         deriv(1,2,inode)=(+1.0-t)*0.25
         deriv(1,3,inode)=(+1.0+t)*0.25
         deriv(1,4,inode)=(-1.0-t)*0.25
         deriv(2,1,inode)=(-1.0+s)*0.25
         deriv(2,2,inode)=(-1.0-s)*0.25
         deriv(2,3,inode)=(+1.0+s)*0.25
         deriv(2,4,inode)=(+1.0-s)*0.25
    	 
	 end do


 endif

end subroutine armotodo
       
      
     