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
double precision :: gl    = 2e-6                ! S/micron

double precision :: PI = 3.14159265358979323846

double precision ::  h_membrana,Tot_por
double precision :: dt,rad,alfa,n_eq,term4,term3,term2,term1,sigmaeff,denom1, denom2,denom,pot_ant,normalx,normaly,sigma_new
integer :: t_final,kk,kt,jj,nsale


CONST_Q = (RADIO_MIN_ENERGIA / RADIO_INICIAL)**2

radio_poro=RADIO_INICIAL

! loop temporal

t_final=1000
nsale=100

! sigmaint,sigmaext,sigmamem

h_membrana = radio_ext-radio_int
area_poros = 0.0
corriente_poros=0.0


do kt=1,t_final

   dt  = DELTA_T_POROS*kt
   if(mod(kt,nsale)==0)  then
       write(unit_den,*) dt
   endif

   do kk=1,nod_mem_ext
       
       jj=nodos_mem(2,kk)
       rad = sqrt(coor_x(jj)*coor_x(jj)+(coor_y(jj)-50.0)*(coor_y(jj)-50.0))
       alfa = atan(coor_x(jj)/(coor_y(jj)-50.0))
       
       if(alfa<0.0) alfa = alfa + 3.14159

       n_eq = DENSIDAD_EQ * exp(CONST_Q * (poten_tm(kk) / V_EP)**2 )

       densi_poros(kk) = densi_poros(kk) + DELTA_T_POROS*ALPHA * exp((poten_tm(kk) / V_EP)**2) * (1 - densi_poros(kk)  / n_eq) 
       
       poros_totales(kk) = densi_poros(kk)*area_zona(kk)

       term1 = poten_tm(kk)**2 * F_MAX/ (1 +R_H/(radio_poro(kk) +  R_T) )
       term2 = 4*BETA*1.e+6*(RADIO_INICIAL/radio_poro(kk))**4 * (1/radio_poro(kk))
       term3 = -2*PI*GAMA
              
       sigmaeff = 2*SIGMA_P - (2*SIGMA_P-SIGMA_0)/(1-area_poros/(4*PI*radio_ext**2))**2

       term4 = 2*PI*sigmaeff* radio_poro(kk)
       
       radio_poro(kk) = radio_poro(kk) + DELTA_T_POROS * DIFF_POROS*1e-6/(K_bolztman*TEMPERATURA) * ( term1 + term2 + term3 + term4)

       denom1 = h_membrana/(PI*sporo* radio_poro(kk) * radio_poro(kk) )
       denom2 = 1/(2*sporo* radio_poro(kk))
       curr_poro(kk) = poten_tm(kk)/(denom1 + denom2)   ! unidades V/S = A


       conductance(kk) = 1.0/(denom1 + denom2)/area_zona(kk) ! S/micron**2
       


       ! saco la membrana
       normalx =  rad*dsin(alfa)
       normaly =  rad*dcos(alfa)
       if(alfa<0.0) alfa = alfa + 3.14159
       
       denom = grad_x(jj)*normalx+grad_y(jj)*normaly
       
       pot_ant = poten_tm(kk)
       sigma_new = (Cm*(poten_tm(kk)-pot_ant)/DELTA_T_POROS + gl*(poten_tm(kk)-pot_rest) + curr_poro(kk)/area_zona(kk))/denom  ! ecuacion (3)


       if(mod(kt,nsale)==0)  then
           write(unit_den,*) alfa, densi_poros(kk),poros_totales(kk),radio_poro(kk),curr_poro(kk),conductance(kk),poten_tm(kk),sigma_new 
       endif

   enddo


   ! calculo corriente total y area total !
   area_poros = 0.0
   corriente_poros=0.0
   Tot_por=0.0
   do kk=1,nod_mem_ext
       area_poros = area_poros + PI*radio_poro(kk)**2 * poros_totales(kk)
       corriente_poros=corriente_poros + poros_totales(kk)* curr_poro(kk)/area_zona(kk)
       Tot_por=Tot_por+poros_totales(kk)
   enddo


   write(unit_por,*) dt,area_poros,corriente_poros,Tot_por
   

enddo




end subroutine densidad_poros