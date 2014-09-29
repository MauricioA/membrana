function veo_porototi()

fi=fopen('poros.dat','r+');

kt=0;

while(feof(fi) == 0)
     
     kt=kt+1;
     
     t = fscanf(fi,'%f',1)
     
     tt(kt)=t;
     
     area(kt)= fscanf(fi,'%f',1);
     curr(kt)= fscanf(fi,'%f',1);
     conduc(kt)= fscanf(fi,'%f',1);
     Ptot(kt)= fscanf(fi,'%f\n',1);
end


subplot(2,2,1);plot(tt,area,'.');title('mcr2')
subplot(2,2,2);plot(tt,Ptot,'.');title('K')
subplot(2,2,4);plot(tt,conduc,'.');title('S')
subplot(2,2,3);plot(tt,curr,'.');title('A')

fclose(fi);




end