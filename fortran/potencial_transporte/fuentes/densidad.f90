subroutine densidad_poros()
use def_variables
implicit none

double precision:: DENSIDAD_INICIAL	= 0;
double precision:: RADIO_INICIAL 		= 510e-6;		!// r* 0.51 nm
double precision:: RADIO_MIN_ENERGIA	= 800e-6;		!// rm 0.80 nm
double precision:: ALPHA				= 1e-3;			!// Coeficiente de creación 1e9 m**-2 s**-1
double precision:: V_EP				= 0.258;		!// Voltaje característico [V]
double precision:: DENSIDAD_EQ		= 1.5e-3;		!// N0 Densidad de poros en equilibrio 1.5e9 m**-2
double precision:: CONST_Q			
double precision:: DIFF_POROS			= 50e-3;		!// D Coeficiente de diffusión para poros 5e-14 m**2 s**-1
double precision:: DELTA_T_POROS		= 1e-7;   ! seg
double precision:: F_MAX				= 0.7e-9;		!// Max fuerza electrica [N V**-2]
double precision:: R_H				= 970e-6;		!// 0.97e-9 m
double precision:: R_T				= 310e-6;		!// 0.31e-9 m
double precision:: BETA 				= 1.4e-19;		!// Repulsión estérica [J]
double precision:: GAMA 				= 1.8e-11;		!// pasado a N         1.8e-11 J m**-1
double precision:: SIGMA_P			= 2e-8;		!// pasada a N/micron           2e-2 J m**-2
double precision:: SIGMA_0			= 1e-12;		!// pasada a N/micron       1e-6 J m**-2
double precision :: TEMPERATURA 		= 310;			!// 37ºC
double precision :: K_bolztman        = 1.38e-23    ! J/K
double precision :: Pot_rest          = -0.08 ! V
double precision :: sporo             = 2e-6        ! S/micron

double precision :: Cm    =1e-14               ! F/micron**2 
double precision :: gl    = 2e-12                ! S/micron**2

double precision :: PI = 3.14159265358979323846

double precision ::  h_membrana
double precision :: time,t_final,rad,alfa,n_eq,term4,term3,term2,term1,sigmaeff,denom1, denom2,denom,pot_ant,normalx,normaly,  &
    &               ra_aux,radio_average
integer :: nfinal,kk,kt,jj,nsale,poros_ant,indice,ii,Tot_por,jmem

integer :: jporomax=1

double precision,allocatable :: aux1(:,:),aux2(:,:),aux3(:,:)
integer, allocatable :: naux(:,:)


CONST_Q = (RADIO_MIN_ENERGIA / RADIO_INICIAL)**2

radio_poro=RADIO_INICIAL

! loop temporal

t_final=2e-3  ! sef

nfinal = t_final/DELTA_T_POROS

nsale=100

! sigmaint,sigmaext,sigmamem

h_membrana = radio_ext-radio_int
area_poros = 0.0
corriente_poros=0.0


do kt=1,nfinal

   time  = DELTA_T_POROS*kt
   if(mod(kt,nsale)==0)  then
       write(unit_den,*) time
   endif




   do kk=1,nod_mem_ext
       
       jj=nodos_mem(2,kk)
       rad = sqrt(coor_x(jj)*coor_x(jj)+(coor_y(jj)-50.0)*(coor_y(jj)-50.0))
       alfa = atan(coor_x(jj)/(coor_y(jj)-50.0))
       
       if(alfa<0.0) then
          alfa = alfa + 3.14159
       endif

       normalx =  rad*dsin(alfa)/rad
       normaly =  rad*dcos(alfa)/rad
      
       denom = -(grad_x(jj)*normalx+grad_y(jj)*normaly)
       pot_ant = poten_tm(kk)
      ! sigma_new(kk) = (Cm*(poten_tm(kk)-pot_ant)/DELTA_T_POROS + gl*(poten_tm(kk)-pot_rest) + currxzona(kk)/area_zona(kk))/denom
       sigma_new(kk) = currxzona(kk)/area_zona(kk)/denom

       if(sigma_new(kk) == 0.0) sigma_new(kk) = sigmamem



       radio_average=0.0
       currxzona(kk) =0.0
       conducxzona(kk)=0.0

       n_eq = DENSIDAD_EQ * exp(CONST_Q * (poten_tm(kk) / V_EP)**2 )

       densi_poros(kk) = densi_poros(kk) + DELTA_T_POROS*ALPHA * exp((poten_tm(kk) / V_EP)**2) * (1 - densi_poros(kk)  / n_eq) 
       
       poros_ant = poros_totales(kk)
       poros_totales(kk) = densi_poros(kk)*area_zona(kk)

       indice = (poros_totales(kk) -poros_ant)

       
       if( indice >= 1 ) then
          
          
          ! zona de reallocacion de lo anterior
          !if(jporo(kk)>jporomax) then
             
             allocate(aux1(nod_mem_ext,jporomax),aux2(nod_mem_ext,jporomax),aux3(nod_mem_ext,jporomax),naux(nod_mem_ext,jporomax))


          do ii=1,nod_mem_ext
             do jj=1,jporomax
                 aux1(ii,jj)=radio_poro(ii,jj)   
                 aux2(ii,jj)=curr_poro(ii,jj) 
                 aux3(ii,jj)=conductance(ii,jj) 
                 naux(ii,jj)=porosxradio(ii,jj) 
             enddo
          enddo


          
             deallocate(radio_poro,curr_poro,conductance,porosxradio)
                    
          !endif
 
          ! sumo nuevos poros
          jporo(kk)=jporo(kk)+1      
          number_poros(kk)=number_poros(kk) + indice
          jmem=jporomax
          if(jporo(kk)>jporomax) then
             jmem = jporo(kk)
          endif
                 
         ! if(jporo(kk)>jporomax) then
              allocate(radio_poro(nod_mem_ext,jmem),curr_poro(nod_mem_ext,jmem),conductance(nod_mem_ext,jmem),porosxradio(nod_mem_ext,jmem))
         ! endif
          radio_poro=0.0
          curr_poro=0.0
          conductance=0.0
          porosxradio=0

          do ii=1,nod_mem_ext
             do jj=1,jporomax
                radio_poro(ii,jj) = aux1(ii,jj)   
                curr_poro(ii,jj) = aux2(ii,jj)
                conductance(ii,jj) = aux3(ii,jj)
                porosxradio(ii,jj) = naux(ii,jj)
             enddo
          enddo
          radio_poro(kk,jporo(kk)) =RADIO_INICIAL
          curr_poro(kk,jporo(kk)) = 0.0
          conductance(kk,jporo(kk)) = 0.0
          porosxradio(kk,jporo(kk)) = indice
          
          deallocate(aux1,aux2,aux3,naux)
          
          if(jporo(kk)>jporomax) then
              jporomax = jporo(kk)
          endif
          

       endif
          
       do jj=1,jporo(kk)
             ra_aux = radio_poro(kk,jj)

             term1 = poten_tm(kk)**2 * F_MAX/ (1 +R_H/(ra_aux +  R_T) )
             term2 = 4*BETA*1.e+6*(RADIO_INICIAL/ra_aux)**4 * (1.0/ra_aux)
             term3 = -2*PI*GAMA
              
             sigmaeff = 2*SIGMA_P - (2*SIGMA_P-SIGMA_0)/(1-area_poros/(4*PI*radio_ext**2))**2

             term4 = 2*PI*sigmaeff* ra_aux
       
             radio_poro(kk,jj) = ra_aux + DELTA_T_POROS * DIFF_POROS*1e-6/(K_bolztman*TEMPERATURA) * ( term1 + term2 + term3 + term4)

            denom1 = h_membrana/(PI*sporo* radio_poro(kk,jj) * radio_poro(kk,jj) )
            denom2 = 1/(2*sporo* radio_poro(kk,jj))
           
            curr_poro(kk,jj) = poten_tm(kk)/(denom1 + denom2)   ! unidades V/S = A
            
            conductance(kk,jj) = 1.0/(denom1 + denom2)  ! S
            
            currxzona(kk) = currxzona(kk) + porosxradio(kk,jj) * curr_poro(kk,jj)
            conducxzona(kk) = conducxzona(kk) + porosxradio(kk,jj)* conductance(kk,jj)

            radio_average = radio_average +  radio_poro(kk,jj)/jporo(kk)

       enddo

       
         if(mod(kt,nsale)==0)  then
            write(unit_den,'(9E15.5)') alfa,area_zona(kk),densi_poros(kk),poros_totales(kk),radio_average,currxzona(kk),conducxzona(kk),poten_tm(kk),sigma_new(kk)  
         endif

   enddo ! loop sobre nodos de membrana

   ! estadistica de poros: poros de tamaño menor y mayor al limite. 
   if(mod(kt,nsale)==0) then
       write(unit_esta,*) time,Tot_por
       write(6,*) time,Tot_por
   endif

   ! calculo corriente total y area total !
   area_poros = 0.0
   corriente_poros=0.0
   Tot_por=0
   conductance_total=0.0
   do kk=1,nod_mem_ext

       jj=nodos_mem(2,kk)
       rad = sqrt(coor_x(jj)*coor_x(jj)+(coor_y(jj)-50.0)*(coor_y(jj)-50.0))
       alfa = atan(coor_x(jj)/(coor_y(jj)-50.0))
       
       if(alfa<0.0) alfa = alfa + 3.14159

       if(mod(kt,nsale)==0)  then

           write(unit_esta,*) kk,alfa,number_poros(kk),jporo(kk)
       endif
       do jj=1,jporo(kk)
          
          area_poros = area_poros + PI*radio_poro(kk,jj)**2 *porosxradio(kk,jj)
          Tot_por=Tot_por+porosxradio(kk,jj)

          if(mod(kt,nsale)==0) then
              write(unit_esta,*) jj,porosxradio(kk,jj),radio_poro(kk,jj)
          endif
      enddo
      corriente_poros=corriente_poros + currxzona(kk)/area_zona(kk)
      conductance_total = conductance_total + conducxzona(kk)

   enddo
   
   write(unit_por,'(4e15.5,i6)') time,area_poros,corriente_poros,conductance_total,Tot_por
   

   


enddo




end subroutine densidad_poros