function borde(pot,CH)

E_eq=1.23; % V
C0_h=0.16*6.02e23/1000/(1e+4*1e+4*1e+4); %at/micro**3
I_eq=1e-18; % A/micro**2

F=96485.34;  %C/mol
R=8.314;     %J/K/mol
T=310.0;      %K


I = I_eq* ( exp(-F*(pot+E_eq)/(2*R*T)) - CH/C0_h*exp(F*(pot+E_eq)/(2*R*T))) 


