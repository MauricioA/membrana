program celula
use def_variables
implicit none


    write(6,*) ' inicio del programa celula'
    call inicio()
    write(6,*) ' leo datos'
    call lectura()
    
    write(6,*) ' resuelvo problema electromagnetico'
    call poisson()
    
    
    write(6,*) ' cierro todo y me voy'
    call finalice()

    write(6,*) ' chau!'

    
end program celula