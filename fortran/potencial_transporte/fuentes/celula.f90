program celula
use def_variables
use def_transpor
implicit none


    write(6,*) ' inicio del programa celula'
    call inicio()
    write(6,*) ' leo datos'
    call lectura()
    
    if(nopcion==1) then     
       write(6,*) '      voy a control!'
       call control()
    else
       write(6,*) '      leo el sistema!'
       call lee_sistema(archi_sistema,nnodes)
       allocate(carga(nnodes))

       carga=0.0

    endif



    if(nmode==1) then
       write(6,*) ' resuelvo problema electromagnetico'
       call poisson()
    elseif(nmode==2 .or. nmode==3) then
       call transporte()
    endif
    
    write(6,*) ' cierro todo y me voy'
    call finalice()

    write(6,*) ' chau!'

    
end program celula