function veo_trans(caso)

if(caso==1)

    fi=fopen('historia.dat','r+');

    N=fscanf(fi,'%i',1)

    ancho=250.0;
    largo=500.0;
    radio_ext=50.005;
    radio_int=50.0;

    posmas(1) = largo*0.5+radio_ext;
    posmas(2) = largo*0.5+radio_ext;
    posmenos(1) = largo*0.5-radio_ext;
    posmenos(2) = largo*0.5-radio_ext;

    con(1)=0.0;
    con(2)=100.0;
    con2(1)=0.0;
    con2(2)=20E+7;

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

         %subplot(2,3,1);plot3(x,y,sol,'k.');
        % subplot(2,2,1);plot3(x,y,ch,'b.');title('H')
        % subplot(2,2,2);plot3(x,y,coh,'r.');title('OH')
        % subplot(2,3,4);plot3(x,y,cna,'k.');title('Na')
        % subplot(2,3,5);plot3(x,y,ccl,'k.');title('Cl')
         subplot(2,2,1);plot(yy,cna_t(k,:),'.' ,posmas,con2,posmenos,con2);title('Na+')
         subplot(2,2,2);plot(yy,ccl_t(k,:),'.' ,posmas,con2,posmenos,con2);title('Cl-')
         subplot(2,2,3);plot(yy,ch_t(k,:),'.'  ,posmas,con,posmenos,con);title('H+')
         subplot(2,2,4);plot(yy,coh_t(k,:),'.' ,posmas,con,posmenos,con);title('OH-')




        % pause
        pause(0.1);
    end
    fclose(fi);

    figure
    %plot(yy,sol_t,'.');
    plot3(x,y,sol,'k.');
    %subplot(2,2,2);plot(yy,ch_t,'.');title('H')
    %subplot(2,2,3);plot(yy,coh_t,'.');title('OH')
    %subplot(2,2,4);plot(yy,cna_t,'.');title('Na')
    %subplot(2,2,1);plot(yy,ccl_t,'.');title('Cl')

else
    
   fi2=fopen('campo.dat','r+');
   NE=fscanf(fi2,'%i',1)
   
   k=0;
   while(feof(fi2)==0)

         k=k+1;
         s=fscanf(fi2,'%s',1);
         t = fscanf(fi2,'%f',1)
         tt(k)=t;
         jy=0;
         for j=1:NE
            x(j)=fscanf(fi2,'%f',1);
            y(j)=fscanf(fi2,'%f',1);
            solx(j)=fscanf(fi2,'%f',1);
            soly(j)=fscanf(fi2,'%f\n',1);
         end
         
         subplot(2,2,1);plot3(x,y,solx,'b.');title('X')
         subplot(2,2,2);plot3(x,y,soly,'r.');title('Y')

         pause
    end
    fclose(fi2);
   
end