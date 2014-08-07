function veo_trans()

fi=fopen('historia.dat','r+');

N=fscanf(fi,'%i',1)

set(0, 'defaultfigurevisible', 'off');

ancho=50.0;
largo=100.0;

med=0.0;

k=0;
while(feof(fi)==0)

     k=k+1;
     s=fscanf(fi,'%s',1);
     t = fscanf(fi,'%f',1)
     tt(k)=t;
     jy=0;
     for j=1:N
        n=fscanf(fi,'%i',1);
        x(j)=fscanf(fi,'%f',1);
        y(j)=fscanf(fi,'%f',1);
        sol(j)=fscanf(fi,'%f',1);
        ch(j)=fscanf(fi,'%f',1);
        coh(j)=fscanf(fi,'%f',1);
        cna(j)=fscanf(fi,'%f',1);
        ccl(j)=fscanf(fi,'%f\n',1);
        
        if(x(j)==med)
            jy=jy+1;
            yy(jy)=y(j);
            ch_t(k,jy)=ch(j);
            coh_t(k,jy)=coh(j);
            cna_t(k,jy)=cna(j);
            ccl_t(k,jy)=ccl(j);
            sol_t(k,jy)=sol(j);
        end
        
     end
        
     subplot(2,3,1);plot3(x,y,sol,'k.');
     subplot(2,3,2);plot3(x,y,ch,'b.');title('H')
     subplot(2,3,3);plot3(x,y,coh,'r.');title('OH')
     subplot(2,3,4);plot3(x,y,cna,'k.');title('Na')
     subplot(2,3,5);plot3(x,y,ccl,'k.');title('Cl')
     
     filename=sprintf('output/%05d.png',k);
     print(filename);
     
     %~ pause
     %pause(0.1);
end
fclose(fi);

%subplot(2,2,1);plot(yy,sol_t,'.');
subplot(2,2,2);plot(yy,ch_t,'.');title('H')
subplot(2,2,3);plot(yy,coh_t,'.');title('OH')
subplot(2,2,4);plot(yy,cna_t,'.');title('Na')
subplot(2,2,1);plot(yy,ccl_t,'.');title('Cl')
