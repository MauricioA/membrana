subroutine lectura()
use def_variables
implicit none
character*(1)::Z 
character*(120):: textinput,opcion
integer :: leng,last,suma,kk


write(Z,'(A1)') Z'09'

do while(textinput /= 'END_DATA')
     read(unit_data,'(A120)') textinput
	 call upcase(textinput)
	 leng=len_trim(textinput) 
     
	 last=0
	 suma=0
	 opcion='  '
	 do while(last<leng)
        last=last+1
		if(textinput(last:last)/=' '.and.textinput(last:last)/=':'.and.textinput(last:last)/='='.and.textinput(last:last)/=Z) then
		  if(textinput(last:last)=='#') then
		    last=leng
		  else
		    suma=suma+1
		    opcion(suma:suma)= textinput(last:last)
		  endif
		endif

	 enddo
     
	 leng=len_trim(opcion)
	 last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='MODO') then
				read(opcion(last+1:leng),'(i3)') nmode
				write(unit_cont,*) ' MODO: ', nmode
				last=leng+1
			endif
        enddo
       last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='DIMEN') then
				read(opcion(last+1:leng),'(i3)') ndimension
				write(unit_cont,*) ' DIMEN: ', ndimension
				last=leng+1
			endif
        enddo


	 last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='ARCHIVO1') then
				read(opcion(last+1:leng),'(a20)') archi_malla
				write(unit_cont,*) ' ARCHIVO1: ', archi_malla
				last=leng+1
			endif
        enddo

 
      last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='NODPEL') then
				read(opcion(last+1:leng),'(i3)') nodpel
				write(unit_cont,*) ' NODPEL: ', nodpel
				last=leng+1
			endif
        enddo

       last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='OPCION') then
				read(opcion(last+1:leng),'(i3)') nopcion
				write(unit_cont,*) ' OPCION: ', nopcion
				last=leng+1
			endif
        enddo
        
         last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='ARCHIVO3') then
				read(opcion(last+1:leng),'(a20)') archi_sistema
				write(unit_cont,*) ' ARCHIVO3: ', archi_sistema
				last=leng+1
			endif
        enddo
	    last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='SIGINT') then
				read(opcion(last+1:leng),'(e15.5)') sigmaint
				write(unit_cont,*) ' SIGINT: ', sigmaint
				last=leng+1
			endif
        enddo
        last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='SIGEXT') then
				read(opcion(last+1:leng),'(e15.5)') sigmaext
				write(unit_cont,*) ' SIGEXT: ', sigmaext
				last=leng+1
			endif
        enddo
        last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='SIGMEM') then
				read(opcion(last+1:leng),'(e15.5)') sigmamem
				write(unit_cont,*) ' SIGMEM: ', sigmamem
				last=leng+1
			endif
        enddo
        last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='PERMIT') then
				read(opcion(last+1:leng),'(e15.5)') permit
				write(unit_cont,*) ' PERMIT: ', permit
				last=leng+1
			endif
        enddo
        last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='POTENCIAL') then
				read(opcion(last+1:leng),'(e15.5)') potencial
				write(unit_cont,*) ' POTENCIAL: ', potencial
				last=leng+1
			endif
        enddo
        last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='FRECUENCIA') then
				read(opcion(last+1:leng),'(e15.5)') frecuencia
				write(unit_cont,*) ' FRECUENCIA: ', frecuencia
				last=leng+1
			endif
        enddo
        last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='RADIOEXT') then
				read(opcion(last+1:leng),'(e15.5)') radio_ext
				write(unit_cont,*) ' RADIOEXT: ', radio_ext
				last=leng+1
			endif
        enddo
        last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='RADIOINT') then
				read(opcion(last+1:leng),'(e15.5)') radio_int
				write(unit_cont,*) ' RADIOINT: ', radio_int
				last=leng+1
			endif
        enddo
        last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='LARGO') then
				read(opcion(last+1:leng),'(e15.5)') largo
				write(unit_cont,*) ' LARGO: ', largo
				last=leng+1
			endif
        enddo
          last=0	 
	   do while(last<leng)
			last=last+1
			if(opcion(1:last)=='ANCHO') then
				read(opcion(last+1:leng),'(e15.5)') ancho
				write(unit_cont,*) ' ANCHO: ', ancho
				last=leng+1
			endif
        enddo
   enddo
   
   

   if(archi_malla/='  ') then
       write(6,*) ' armo mallado desde los archivos! '
       call leemalla()
   else
      
       write(6,*) ' armo malla interna de elementos de cuatro nodos! '
       call armo_malla()


   endif
       
   close(unit_data)
end subroutine lectura

subroutine lee_sistema(archi_sistema,nnodes)
use def_solver
implicit none
character(20)::archi_sistema
integer :: nnodes
!local
integer::kk,n,j
integer :: unit_sist

open(unit=unit_sist,file=archi_sistema)



read(unit_sist,*) n
allocate(ia(nnodes+1),cx(nnodes+1))
if(n/=nnodes+1) then
   write(6,*) 'sistema incorrecto!!'
   stop ' '
endif

do kk=1,n

   read(unit_sist,*) j,ia(kk),cx(kk)

enddo

read(unit_sist,*) nonull
 allocate(ad(nnodes),rhs(nnodes),solucion(nnodes))
 solucion=0.0

allocate(ja(nonull),an(nonull))
do kk=1,nonull

   read(unit_sist,*) j,ja(kk)

enddo

read(unit_sist,*) n
allocate(ia2d(2*nnodes+1),cx2d(2*nnodes+1))
if(n/=2*nnodes+1) then
   write(6,*) 'sistema incorrecto!!'
   stop ' '
endif

do kk=1,n

   read(unit_sist,*) j,ia2d(kk),cx2d(kk)

enddo

read(unit_sist,*) nonull2d
 allocate(ad2d(2*nnodes),rhs2d(2*nnodes),solucion2d(2*nnodes))

allocate(ja2d(nonull2d),an2d(nonull2d))
do kk=1,nonull2d

   read(unit_sist,*) j,ja2d(kk)

enddo


close(unit_sist)



end subroutine lee_sistema

subroutine upcase(word)
    !-----------------------------------------------------------------------
    !
    !     This routine converts wopos to upper case 
    !
    !-----------------------------------------------------------------------
    implicit none
    character(*), intent(inout) :: word
    integer                     :: iposi,ioctv

    do iposi=1,60                                   ! process all positions
       ioctv=ichar(word(iposi:iposi))               ! octal value
       if(o'141'<=ioctv.and.o'172'>=ioctv) then ! it is a lower case
          ioctv=ioctv-o'40'                          ! equivalent upper case
          word(iposi:iposi)=char(ioctv)              ! convert it to upcase
       end if
    end do ! iposi=1,5

  end subroutine upcase


