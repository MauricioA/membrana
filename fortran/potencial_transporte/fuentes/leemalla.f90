subroutine leemalla()
use def_variables
implicit none
! local
character(50) :: tilte
character(5) :: tri
integer :: kk,kk1,kk2,kk3,jj,ii,j,n,ns(3),nelet,npiso,jj1,jj2

double precision :: rad,dista,rmed
double precision ::PI=3.14159265358979323846


open(unit=1111,file=archi_malla)

read(1111,*) tilte
read(1111,*) nnodes

allocate(coor_x(nnodes),coor_y(nnodes))
do kk=1,nnodes
    read(1111,*) n,coor_x(kk),coor_y(kk)
enddo

read(1111,*) tilte
read(1111,*) num_mat
allocate(matelem(num_mat))
nelements=0
do kk=1,num_mat
    read(1111,*) n,matelem(kk),tri
    nelements=nelements+matelem(kk)
enddo
  
read(1111,*) tilte


allocate(conect(nelements,nodpel),material(nelements))
allocate(tension(nelements,ntension),deforma(nelements,ntension))

npiso=0
do jj=1,num_mat
   do kk=1, matelem(jj)
       read(1111,*) n,(conect(npiso+kk,j),j=1,nodpel)
       material(npiso+kk)=jj
   enddo
   npiso=npiso+matelem(jj)
enddo          

close(1111)



mat_externo=0
mat_interno=0
mat_membrana=0

! determino el material de cada elemento
do kk=1,nelements
   
   
   if(material(kk)==1) then !! exterior
      mat_externo=mat_externo+1
   elseif(material(kk)==3) then ! interno
      mat_interno=mat_interno+1
   elseif(material(kk)==2) then ! membrana
      mat_membrana=mat_membrana+1
   endif

enddo

allocate(vec_tierra(nnodes),vec_poten(nnodes),vec_quieto(2*nnodes))

vec_tierra=-1
vec_poten=-1
vec_quieto=-1
ndirichV=0
ndirichT=0


do kk=1,nnodes
   if(coor_y(kk)>=largo*0.999) then
       ndirichV=ndirichV+1
   endif

   if(coor_y(kk)<=0.0001) then
       ndirichT=ndirichT+1
   endif

enddo

!deallocate(nod_dirichV,nod_dirichT)
allocate(nod_dirichV(ndirichV),nod_dirichT(ndirichT))
ndirichV=0
ndirichT=0

do kk=1,nnodes
   if(coor_y(kk)>=largo*0.999) then
       ndirichV=ndirichV+1
       nod_dirichV(ndirichV)=kk
   endif

   if(coor_y(kk)<=0.0001) then
       ndirichT=ndirichT+1
       nod_dirichT(ndirichT)=kk
   endif

enddo



 
do kk=1,ndirichV
    vec_poten(nod_dirichV(kk))=potencial
    vec_quieto(2*nod_dirichV(kk))=0
    vec_quieto(2*nod_dirichV(kk)-1)=0
enddo
 
do kk=1,ndirichT
    vec_tierra(nod_dirichT(kk))=0
    vec_quieto(2*nod_dirichT(kk))=0
    vec_quieto(2*nod_dirichT(kk)-1)=0
enddo



nod_tierra=ndirichT
nod_poten=ndirichV
nquieto=2*ndirichT+2*ndirichV

! continuo con las condiciones de borde para el probelma elastico
do kk=1,nelements
   if(material(kk)==2) then
       
       do jj=1,nodpel
           if(coor_x(conect(kk,jj)) == 0.0) then
              nquieto=nquieto +1
              vec_quieto(2*conect(kk,jj)-1)=0
           endif
           if(conect(kk,jj) ==6 .or. conect(kk,jj) ==9 ) then
              nquieto=nquieto +1
              vec_quieto(2*conect(kk,jj))=0
           endif
       enddo
   endif
enddo


allocate(nod_quieto(nquieto))
nquieto=0
do kk=1,2*nnodes
   if(vec_quieto(kk)==0) then
       nquieto=nquieto+1
       nod_quieto(nquieto)=kk
   endif
enddo

write(unit_cc,* )' materiales ', nelements
write(unit_cc,* )'      externo:   ', mat_externo
write(unit_cc,* )'      interno:   ', mat_interno
write(unit_cc,* )'      memebrana: ', mat_membrana

write(unit_cc,* )'  ************************************************  '
write(unit_cc,* )'  ************************************************  '


write(unit_cc,* )'nodos tierra ', nod_tierra
      
do kk=1,nod_tierra
      write(unit_cc,'(2i6,2e15.5)') kk,nod_dirichT(kk),coor_x(nod_dirichT(kk)), coor_y(nod_dirichT(kk))
enddo

write(unit_cc,* )'nodos poten ', nod_poten
do kk=1,nod_poten
      write(unit_cc,'(2i6,2e15.5)') kk,nod_dirichV(kk), coor_x(nod_dirichV(kk)),coor_y(nod_dirichV(kk))
enddo

write(unit_cc,* )'nodos quietos ', nquieto
do kk=1,nquieto
      jj=kk/2
      if(kk - 2*jj == 1) jj=jj+1
      
      write(unit_cc,'(2i6,2e15.5)') kk,nod_quieto(kk), coor_x(jj),coor_y(jj)
enddo
 
allocate(grad_x(nnodes),grad_y(nnodes),grad_z(nnodes),gradxel_x(nelements),gradxel_y(nelements),gradxel_z(nelements))




! determino nodos de la membranas

nod_mem_ext=64
nod_mem_int=64


!do kk=1,nnodes
!   rad = sqrt(coor_x(kk)*coor_x(kk) + (coor_y(kk)-50.0)*(coor_y(kk)-50.0))
!   if(rad > radio_int*0.999 .and. rad < radio_int*1.001) then
!      nod_mem_int=nod_mem_int+1
!   endif
!   if(rad > radio_ext*0.999 .and. rad < radio_ext*1.001 ) then
!      nod_mem_ext=nod_mem_ext+1
!   endif
!enddo


       
   allocate(sigma_new(nod_mem_ext))
allocate(nodos_mem(2,nod_mem_ext),densi_poros(nod_mem_ext),poten_tm(nod_mem_ext),radio_poro(nod_mem_ext,1),curr_poro(nod_mem_ext,1),poros_totales(nod_mem_ext))
allocate(area_zona(nod_mem_ext),conductance(nod_mem_ext,1))
allocate(porosxradio(nod_mem_ext,1),number_poros(nod_mem_ext),jporo(nod_mem_ext))

allocate( currxzona(nod_mem_ext),conducxzona(nod_mem_ext))

sigma_new=0.0
currxzona=0.0
conducxzona=0.0

densi_poros=0.0
radio_poro=0.0
curr_poro=0.0
poros_totales=0.0
conductance=0.0
jporo=0
porosxradio=0
number_poros=0


do kk=1,nod_mem_ext
    nodos_mem(2,kk) = kk+20
    nodos_mem(1,kk) = kk+1253
enddo


write(unit_cc,* )'nodos membrana interna', nod_mem_int
do kk=1,nod_mem_int
   jj=nodos_mem(1,kk)
   rad = sqrt(coor_x(jj)*coor_x(jj) + coor_y(jj)*coor_y(jj))
   write(unit_cc,*) kk,nodos_mem(1,kk),coor_x(jj),coor_y(jj)
enddo

! calculo area de la zona asociada a cada nodo externo
do kk=1,nod_mem_ext-1
   jj1=nodos_mem(2,kk)
   jj2=nodos_mem(2,kk+1)
   
   dista = sqrt((coor_x(jj1)-coor_x(jj2))**2 + (coor_y(jj1)-coor_y(jj2))**2 )
   
   rmed = (coor_x(jj1)+coor_x(jj2))*0.5


   if(kk==1) then
      area_zona(kk) = dista*PI*rmed
   else

      area_zona(kk) = dista*2*PI*rmed

      if(kk==nod_mem_ext-1) then
          area_zona(kk+1) = dista*PI*rmed
      endif

   endif

enddo


write(unit_cc,* )'nodos membrana externa', nod_mem_ext
do kk=1,nod_mem_ext
   jj=nodos_mem(2,kk)
   rad = sqrt(coor_x(jj)*coor_x(jj) + coor_y(jj)*coor_y(jj))
   
   
   write(unit_cc,*) kk,nodos_mem(2,kk),coor_x(jj),coor_y(jj)




enddo




end subroutine leemalla

