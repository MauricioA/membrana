subroutine salida_sol(solucion)
use def_variables
use def_constantes
implicit none
double precision :: solucion(nnodes)
! local 
integer :: kk,jj,ii,ngrid,nzgrid2,casex,ncamp,j
double precision :: camp,xmed,ymed,zmed,col,z=0,solmed
character*(1)::coma=','



if(ndimension==2) then

    write(unit_sal,*) 'X,    Y,   V  ' ! nnodes

    do kk=1,nnodes
       write(unit_sal,'(e15.5,a,e15.5,a,e15.5)') coor_x(kk),coma,coor_y(kk),coma,solucion(kk)
    enddo

    write(unit_gra,*)  ' x,     y,    campx, campy'

    write(unit_camp,*)  ' x,     y,    camp'
    write(unit_2d,*)  ' x,     y,    camp'

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
          write(unit_2d,'(e15.5,2(a,e15.5))') xmed,coma,ymed,coma,camp
       endif

       write(unit_gra,'(e15.5,3(a,e15.5))')   xmed,coma,ymed,coma,gradxel_x(kk),coma,gradxel_y(kk)


    enddo

elseif(ndimension==1) then

    write(unit_sal,*) 'X,     V  ' ! nnodes

    do kk=1,nnodes
       write(unit_sal,'(e15.5,a,e15.5)') coor_x(kk),coma,solucion(kk)
    enddo

    write(unit_gra,*)  ' x,   campx'

    write(unit_camp,*)  ' x,     camp'
    write(unit_2d,*)  ' x,      camp'

    do kk=1,nelements
   
       camp = dsqrt(gradxel_x(kk)*gradxel_x(kk) )
       xmed=0
       solmed=0
       do jj=1,nodpel
            j=conect(kk,jj)
            xmed = xmed + coor_x(j)/real(nodpel)
            solmed=solmed + solucion(j)/real(nodpel)
       enddo

       write(unit_camp,'(e15.5,1(a,e15.5))')   xmed,coma,camp

       if(material(kk)==2) then
          write(unit_2d,'(e15.5,(a,e15.5))') xmed,coma,camp
       endif

       write(unit_gra,'(e15.5,1(a,e15.5))')   xmed,coma,gradxel_x(kk)


    enddo


endif



end subroutine salida_sol



subroutine salida_concentra(ch,coh)
use def_variables
use def_constantes
implicit none
double precision :: ch(nnodes),coh(nnodes)
! local 
integer :: kk
character*(1)::coma=','



    write(unit_ch,*) 'X,    Y,   CH  ' ! nnodes

    do kk=1,nnodes
       write(unit_ch,'(e15.5,a,e15.5,a,e15.5)') coor_x(kk),coma,coor_y(kk),coma,ch(kk)
    enddo

    write(unit_coh,*) 'X,    Y,   COH  ' ! nnodes

    do kk=1,nnodes
       write(unit_coh,'(e15.5,a,e15.5,a,e15.5)') coor_x(kk),coma,coor_y(kk),coma,coh(kk)
    enddo



end subroutine salida_concentra



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