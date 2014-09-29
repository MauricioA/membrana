function veo_densi(N)

fi=fopen('densidad.dat','r+');

kt=0;

while(feof(fi) == 0)

     t = fscanf(fi,'%f',1)
     kt=kt+1;
     
     tt(kt)=t;
     
     for kk=1:N
        alf(kk)= fscanf(fi,'%f',1);
        area(kk)= fscanf(fi,'%f',1);
        dens(kk)= fscanf(fi,'%f',1);
        Ptot(kk)= fscanf(fi,'%f',1);
        rad(kk)= fscanf(fi,'%f',1);
        curr(kk)= fscanf(fi,'%f',1);
        cond(kk)= fscanf(fi,'%f',1);
        ptm(kk)= fscanf(fi,'%f',1);
        sig(kk)= fscanf(fi,'%f\n',1);
     end
     
     subplot(3,3,1);plot(alf,dens,'.');title('K/mcr2')
     subplot(3,3,2);plot(alf,Ptot,'.');title('K')
     subplot(3,3,3);plot(alf,rad,'.');title('mcr')
     subplot(3,3,4);plot(alf,curr,'.');title('A')
     subplot(3,3,5);plot(alf,cond,'.');title('S')
     subplot(3,3,6);plot(alf,sig,'.');title('S/mcr')
     subplot(3,3,7);plot(alf,area,'.');title('mcr2')
     subplot(3,3,8);plot(alf,ptm,'.');title('V')

     pause(0.1);

end

fclose(fi);


end