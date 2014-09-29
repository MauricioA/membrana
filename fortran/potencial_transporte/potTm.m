fid=fopen('grid2D.csv','r+');

s=fscanf(fid,'%s',1);

N=fscanf(fid,'%i',1);

for kk=1:N
    xe(kk)=fscanf(fid,'%f',1);
    ye(kk)=fscanf(fid,'%f',1);
    
    rade(kk)=sqrt(xe(kk)*xe(kk)+(ye(kk)-50.0)*(ye(kk)-50.0) );
    alfae(kk) = atan(xe(kk)/(ye(kk)-50.0));
    rare(kk)=fscanf(fid,'%f',1);
    alfe(kk)=fscanf(fid,'%f',1);
    Vex(kk)=fscanf(fid,'%f',1);
    Eex(kk)=fscanf(fid,'%f\n',1);
end

s=fscanf(fid,'%s',1);

N=fscanf(fid,'%i',1);

for kk=1:N
    x(kk)=fscanf(fid,'%f',1);
    y(kk)=fscanf(fid,'%f',1);
    rar(kk)=fscanf(fid,'%f',1);
    alf(kk)=fscanf(fid,'%f',1);

    Vin(kk)=fscanf(fid,'%f',1);
    Ein(kk)=fscanf(fid,'%f\n',1);
end

N=fscanf(fid,'%i',1);

for kk=1:N
    Vm(kk)=fscanf(fid,'%f\n',1);
end

fclose(fid);

subplot(2,3,1);plot(alf,Vex,'.')
subplot(2,3,2);plot(alf,Vin,'.')
subplot(2,3,3);plot(alf,Eex,'.')
subplot(2,3,4);plot(alf,Ein,'.')
subplot(2,3,5);plot(alf,Vm,'.')


pause


% ahora evoluciono poros

RADIO_INICIAL 		= 510e-6;		%!// r* 0.51 nm
RADIO_MIN_ENERGIA	= 800e-6;		%!// rm 0.80 nm
ALPHA				= 1e-3;			%!// Coeficiente de creación 1e9 m**-2 s**-1
V_EP				= 0.258;		%!// Voltaje característico [V]
DENSIDAD_EQ		= 1.5e-3;		%!// N0 Densidad de poros en equilibrio 1.5e9 m**-2

DELTA_T_POROS		= 10000e-9;


CONST_Q = (RADIO_MIN_ENERGIA / RADIO_INICIAL)^2;


%! loop temporal

t_final=1000;

d_poros(1:N)=0.0;

for kt=1:t_final

   dt  = DELTA_T_POROS*kt
   

   for kk=1:N
       
       n_eq = DENSIDAD_EQ * exp(CONST_Q * (Vm(kk)/ V_EP)^2 )

       
       d_poros(kk) = d_poros(kk) + DELTA_T_POROS*ALPHA * exp((Vm(kk) / V_EP)^2) * (1 - d_poros(kk)  / n_eq); 

       

   end

   plot(alfe,d_poros,'.');
   pause;%(0.1)
   
end




