subroutine calculo_carga()
use def_transpor
use def_variables
implicit none
! local
integer :: kk
double precision :: cte,cte_dilucion

cte = 1e+6/6.03e+23  !paso a Mol/micron2
cte_dilucion = ch_inicial / cna_inicial

do kk=1,nnodes
    
   carga(kk) = Faraday/(epsilon*epsilon_0)* (zh*ch(kk)*cte + zoh*coh(kk)*cte +zna*cna(kk)*cte*cte_dilucion + zcl*ccl(kk)*cte*cte_dilucion)

enddo

end subroutine calculo_carga

