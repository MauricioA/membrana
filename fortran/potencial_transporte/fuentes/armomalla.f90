subroutine armo_malla()

use def_variables

! local
integer, allocatable :: grilla(:,:)


if(ndimension==2) then


    nodpel=4

    ladoZ=100.0
    ladoR=50.0


    size=2.0


    Nel_z=LadoZ/size
    Nel_r=LadoR/size

    nnodes = (Nel_z+1)*(Nel_r+1)

    nelements = Nel_z*Nel_r

    allocate(coor_x(nnodes),coor_y(nnodes),conect(nelements,nodpel),material(nelements),vec_tierra(nnodes),vec_poten(nnodes))
    allocate(grilla(nnodes,nnodes))
    allocate(grad_x(nnodes),grad_y(nnodes),grad_z(nnodes),gradxel_x(nelements),gradxel_y(nelements),gradxel_z(nelements))

    grilla=0

    nno=0

    do kk=1,Nel_z+1
    
        do jj=1,Nel_r+1
       
           nno=nno+1 
       
           coor_y(nno) = (kk-1)*size
           coor_x(nno) = (jj-1)*size
       
           grilla(kk,jj)=nno

        enddo


    enddo


    ne=0
    do kk=1,Nel_z
    
        do jj=1,Nel_r
       
           ne=ne+1 
       
           conect(ne,1)=grilla(kk,jj)
           conect(ne,2)=grilla(kk,jj+1)
           conect(ne,3)=grilla(kk+1,jj+1)
           conect(ne,4)=grilla(kk+1,jj)
           material(ne)=1
        enddo

    enddo


    deallocate(grilla)


    ! ahora CC

    nod_tierra=0
    nod_poten=0

    vec_tierra=-1
    vec_poten=-1

    do kk=1,nnodes
       if(coor_y(kk)==0.0) then
          nod_tierra=nod_tierra+1
          vec_tierra(kk)=1
       endif
       if(coor_y(kk)==LadoZ) then
          nod_poten=nod_poten+1
          vec_poten(kk)=1
       endif
    enddo


elseif(ndimension==1) then


    nodpel=6
    nope1=nodpel-1

    ladoZ=50.0

    size=0.5


    Nel_z=LadoZ/size

    nnodes = Nel_z*(nodpel-1)+1

    nelements = Nel_z

    allocate(coor_x(nnodes),coor_y(nnodes),conect(nelements,nodpel),material(nelements),vec_tierra(nnodes),vec_poten(nnodes))
    allocate(grad_x(nnodes),grad_y(nnodes),gradxel_x(nelements),gradxel_y(nelements))
    
    coor_y=0.0
    gradxel_y=0.0
    grad_y=0.0

    nno=0

    base=0.0 

    do kk=1,nelements
    
        do jj=1,nope1
       
           nno=nno+1 
       
           coor_x(nno) = base + (jj-1)*size/real(nope1)
       

        enddo

        base = kk*size

    enddo

    nno=nno+1
    coor_x(nno) = ladoZ


    nnodes=nno

    ne=0
    
    do kk=1,nelements
       material(kk)=1
  
       do jj=1,nodpel
          conect(kk,jj) = nope1*(kk-1)+ jj
       enddo
    enddo


    ! ahora CC

    nod_tierra=0
    nod_poten=0

    vec_tierra=-1
    vec_poten=-1

    nod_tierra=nod_tierra+1
    vec_tierra(1)=1
    nod_poten=nod_poten+1
    vec_poten(nnodes)=1

endif

end subroutine armo_malla