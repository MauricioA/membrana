subroutine salida_sol(solucion)
use def_variables
use def_constantes
implicit none
double precision :: solucion(nnodes)
! local 
integer :: kk,jj,ii,ngrid,nzgrid2,casex,ncamp,j
double precision :: camp,xmed,ymed,zmed,col,z=0,solmed
character*(1)::coma=','

write(unit_sal,*) 'X,    Y,   V  ' ! nnodes

do kk=1,nnodes
   write(unit_sal,'(e15.5,a,e15.5,a,e15.5)') coor_x(kk),coma,coor_y(kk),coma,solucion(kk)
enddo


write(unit_camp,*)  ' x,     y,   camp'

do kk=1,nelements
   
   camp = dsqrt(gradxel_x(kk)*gradxel_x(kk) + gradxel_y(kk)*gradxel_y(kk))
   xmed=0
   ymed=0
   zmed=0
   solmed=0
   do jj=1,nodpel
        j=conect(kk,jj)
        xmed = xmed + coor_x(j)/real(nodpel)
        ymed = ymed + coor_y(j)/real(nodpel)
        solmed=solmed + solucion(j)/real(nodpel)
   enddo

   write(unit_camp,'(e15.5,2(a,e15.5))')   xmed,coma,ymed,coma,camp

   if(material(kk)==2) then
      write(unit_2d,'(4e15.5)') xmed,ymed,solmed, camp
   endif



enddo




end subroutine salida_sol



subroutine salida_tension(solucion)
use def_variables
use def_constantes
implicit none
double precision :: solucion(2*nnodes)
! local 
integer :: kk,jj,ii,ngrid,nzgrid2,casex,ncamp,j
double precision :: camp,xmed,ymed,zmed,col,z=0,solmed,tenmed
character*(1)::coma=','

!write(unit_sal,*) 'X,    Y,  Z,   V  ' ! nnodes

do kk=1,nnodes
   write(unit_sal2d,'(e15.5,a,e15.5,a,e15.5,a,e15.5,a,e15.5)') coor_x(kk),coma,coor_y(kk),coma,solucion(2*kk-1),coma,solucion(2*kk)
enddo



do kk=1,nelements
    xmed=0.0
    ymed=0.0   
    do jj=1,nodpel
       xmed = xmed + coor_x(conect(kk,jj))/real(nodpel)
       ymed = ymed + coor_y(conect(kk,jj))/real(nodpel)
   enddo
   tenmed = sqrt(tension(kk,1)*tension(kk,1)+tension(kk,2)*tension(kk,2))
   write(unit_el2d,'(e15.5,a,e15.5,a,e15.5,a,e15.5)') xmed,coma,ymed,coma,z,coma,tenmed
enddo

end subroutine salida_tension
