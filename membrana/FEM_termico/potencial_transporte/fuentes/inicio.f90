subroutine inicio()
use def_variables
implicit none

    open(unit=unit_data,file=filedata)
    open(unit=unit_malla,file=filemalla)
    open(unit=unit_cc,file=filecc)
    open(unit=unit_cont,file=file_aux)
    open(unit=unit_sal,file=file_sal)
    open(unit=unit_2d,file=file_2D)
    open(unit=unit_gra,file=file_gra)
    open(unit=unit_grid,file=file_grid)
    open(unit=unit_camp,file=file_camp)
    open(unit=unit_sal2d,file=file_sal2d)
    open(unit=unit_el2d,file=file_el2d)
    open(unit=unit_ch,file=file_ch2d)
    open(unit=unit_coh,file=file_coh2d)
    open(unit=unit_histo,file=file_histo)
    open(unit=unit_ph,file=file_ph)
    
end subroutine inicio