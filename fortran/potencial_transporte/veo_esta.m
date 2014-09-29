function veo_esta(N)

fi=fopen('estadis.dat','r+');

kt=0;

rmin = 0.005; % 2 nanometros


while(feof(fi) == 0)

    kt=kt+1;

     t = fscanf(fi,'%f',1)
     tt(kt)=t;
     Ptot(kt)=fscanf(fi,'%f\n',1);
    
     for kk=1:N
         
         posi(kk)=kk;
         nk=fscanf(fi,'%i',1);
         alfa(kk)=fscanf(fi,'%f',1);

         Np(kk,kt)=fscanf(fi,'%i',1);
         jp(kk,kt)=fscanf(fi,'%i\n',1);
         
         rad_ave(kk)=0.0;
         rad_max(kk)=0.0;
         rad_min(kk)=1.0;
         Nrmin(kk)=0;
         Nrmax(kk)=0;
         
         for jj=1:jp(kk,kt)
             n = fscanf(fi,'%i',1);
             pxr(jj) = fscanf(fi,'%f',1);
             rad(jj) = fscanf(fi,'%f\n',1);
             
             rad_ave(kk)=rad_ave(kk) + rad(jj)* pxr(jj)/Np(kk,kt);
             
             if(rad(jj)> rad_max(kk))
                 rad_max(kk)=rad(jj);
             end
             if(rad(jj)< rad_min(kk))
                 rad_min(kk)=rad(jj);
             end
             
                 
             
             if(rad(jj)< rmin)
                 Nrmin(kk)=Nrmin(kk)+1;
             else
                 Nrmax(kk)= Nrmax(kk)+1;
             end
                  
         end
         
            
         
     end
     
     subplot(2,2,1);plot(alfa,rad_ave,'.');title('<r>')
     subplot(2,2,2);plot(alfa,rad_max,'.');title('rmax')
     %subplot(2,2,2);plot(posi,rad_min);title('rmin')
     
     subplot(2,2,3);plot(alfa,Nrmin,'-');title('Nmin')
     subplot(2,2,4);plot(alfa,Nrmax,'-');title('Nmax')
     
     pause(0.1);
     
end


fclose(fi);