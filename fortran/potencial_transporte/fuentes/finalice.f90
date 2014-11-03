subroutine finalice()
use def_variables
implicit none

    close(unit_data)
    close(unit_malla)
    close(unit_cc)
    close(unit_cont)
    close(unit_sal)
    close(unit_2d)
    close(unit_gra)
    close(unit_grid)
    close(unit_camp)
    close(unit_sal2d)
    close(unit_el2d)
    close(unit_ch)
    close(unit_coh)
    close(unit_histo)
    close(unit_ph)
    close(unit_den)
    close(unit_por)
    close(unit_esta)
    close(unit_celda)
    close(unit_ph2)
    close(unit_camp2)


    
end subroutine finalice