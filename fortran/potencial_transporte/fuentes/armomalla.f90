subroutine armo_malla()

use def_variables
implicit none
! local
integer, allocatable :: grilla(:,:)
double precision, allocatable ::radio(:)
double precision ::Ladoz,Lador,size,angulo,pi,base
integer :: nsuperf,mater,Nel_z,Nel_r,nno,kk,jj,ne,nnod,nope1,nrads


PI=3.14159
if(ndimension==2) then
    mater=1
    if(mater==1) then

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


         write(111,*) nnodes
        
        do kk=1,nnodes
            write(111,*) kk,coor_x(kk),coor_y(kk)
        enddo

        write(111,*) nelements
        do kk=1,nelements
            write(111,*) kk,conect(kk,1),conect(kk,2),conect(kk,3),conect(kk,4)
        enddo



        deallocate(grilla)
    
    elseif(mater==3) then
        nodpel=4
        
        nsuperf=9
        nrads=22
        ladoZ=100.0
        ladoR=50.0
        
        nnodes=10000
        allocate(coor_x(nnodes),coor_y(nnodes),radio(nrads),conect(10000,nodpel))
        allocate(grilla(nnodes,nnodes))
        
        do kk=1,10
           radio(kk) = kk*1.0
        enddo
        do kk=11,12
           radio(kk) = 10.0 + (kk-10)*0.5
        enddo
        do kk=13,nrads
           radio(kk) = 11.0 + (kk-12)*(ladoR-11.0)/(nrads-12)
        enddo



       
        coor_x=0.0
        coor_z=0.0
        nnod=0

        do kk=1,nsuperf+1
            angulo=(kk-1)*pi*0.5/real(nsuperf)
            do jj=1,nrads
                nnod=nnod+1
                
                grilla(jj,kk)=nnod

                coor_x(nnod) = radio(jj)*cos(angulo)
                coor_y(nnod) = radio(jj)*sin(angulo)
                if(jj==nrads) then
                   if((ladoR-coor_x(nnod)) > (ladoR-coor_y(nnod)) ) then
            
                        coor_y(nnod) = ladoR
                   else
                        coor_x(nnod) = ladoR
                   endif
               endif
               if(jj==1) then
                  if((coor_x(nnod)) < 0.75 ) then
            
                        coor_y(nnod) = radio(jj)
                   else
                        coor_x(nnod) = radio(jj)
                   endif
               endif
            enddo



        enddo


        nelements=0
        do kk=1,nsuperf
            
            do jj=1,nrads-1
               nelements=nelements+1
               conect(nelements,1)= grilla(jj,kk)
               conect(nelements,2)= grilla(jj,kk+1)
               conect(nelements,3)= grilla(jj+1,kk+1)
               conect(nelements,4)= grilla(jj+1,kk)
            enddo


        enddo

        ! creo elementos adicionales


        write(111,*) nnod
        
        do kk=1,nnod
            write(111,*) kk,coor_x(kk),coor_y(kk)
        enddo

        write(111,*) nelements
        do kk=1,nelements
            write(111,*) kk,conect(kk,1),conect(kk,2),conect(kk,3),conect(kk,4)
        enddo


        Nel_z=LadoZ/size
        Nel_r=LadoR/size

        nnodes = (Nel_z+1)*(Nel_r+1)

        nelements = Nel_z*Nel_r

        allocate(coor_x(nnodes),coor_y(nnodes),conect(nelements,nodpel),material(nelements),vec_tierra(nnodes),vec_poten(nnodes))
        
        allocate(grad_x(nnodes),grad_y(nnodes),grad_z(nnodes),gradxel_x(nelements),gradxel_y(nelements),gradxel_z(nelements))

        grilla=0

    endif

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