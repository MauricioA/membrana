function Vel_2d(ndimy,ndimx,dir,ntop);
%function Vel_2d(ii,jj)  ejey ejex

archivo1=strcat(dir,'VeloU.dat');
archivo2=strcat(dir,'VeloV.dat');
archivo3=strcat(dir,'vortic.dat');
archivo4=strcat(dir,'Corrien.dat');
archivo5=strcat(dir,'dbm.dat')

archivo6=strcat(dir,'SalU.dat')
archivo7=strcat(dir,'SalV.dat')
archivo8=strcat(dir,'SalW.dat')

fid1 = fopen(archivo1,'rb+');
fid2 = fopen(archivo2,'rb+');
fid3 = fopen(archivo3,'rb+');
fid4 = fopen(archivo4,'rb+');
dbm = fopen(archivo5,'rb+');

saleu = fopen(archivo6,'wb+');
salev = fopen(archivo7,'wb+');
salew = fopen(archivo8,'wb+');

tt = fscanf(dbm,'%f\n',1)

db = fscanf(dbm,'%i\n',[ndimx,ndimy]);

hold
for jj=1:ndimy
  for kk=1:ndimx
      xx(kk,jj)=(kk-1);
      yy(kk,jj)=(jj-1);
      if(jj==ntop)
         xx2(kk)=xx(kk,ntop)/ndimx;
      end
            
      if(db(kk,jj)==-1)
         plot(xx(kk,jj),yy(kk,jj),'k.');axis([0 ndimx 0 ndimy])
       end
   end
end
pause

is=0;
while(feof(fid1)==0)
    is=is+1
    
%    VeU = fscanf(fid1,'%f\n',[ndimy,ndimx]);
%    VeV = fscanf(fid2,'%f\n',[ndimy,ndimx]);
%    Wor = fscanf(fid3,'%f\n',[ndimy,ndimx]);
%    PHI = fscanf(fid4,'%f\n',[ndimy,ndimx]);
    VeU = fscanf(fid1,'%f\n',[ndimx,ndimy]);
    VeV = fscanf(fid2,'%f\n',[ndimx,ndimy]);
    Wor = fscanf(fid3,'%f\n',[ndimx,ndimy]);
    PHI = fscanf(fid4,'%f\n',[ndimx,ndimy]);
VeU = VeU'; 
VeV = VeV'; 
Wor = Wor'; 
PHI = PHI'; 

for k=1:ndimx
  Uef(k) = VeU(ntop,k);
  Vef(k) = VeV(ntop,k);
  Wef(k) = Wor(ntop,k);
end

for k=1:ndimx
  Uef1(k) = VeU(ntop,k);
  Vef1(k) = VeV(ntop,k);
  Wef1(k) = Wor(ntop+2,k);
end

for k=1:ndimx
  Uef2(k) = VeU(ntop+5,k);
  Vef2(k) = VeV(ntop+5,k);
  Wef2(k) = Wor(ntop+5,k);
end

for k=1:ndimx
  Uef3(k) = VeU(ntop+10,k);
  Vef3(k) = VeV(ntop+10,k);
  Wef3(k) = Wor(ntop+10,k);
end
for k=1:ndimx
  Uef4(k) = VeU(ntop+15,k);
  Vef4(k) = VeV(ntop+15,k);
  Wef4(k) = Wor(ntop+15,k);
end
for k=1:ndimx
  Uef5(k) = VeU(ntop+20,k);
  Vef5(k) = VeV(ntop+20,k);
  Wef5(k) = Wor(ntop+20,k);

end
for k=1:ndimx
  Uef6(k) = VeU(ntop+25,k);
  Vef6(k) = VeV(ntop+25,k);
  Wef6(k) = Wor(ntop+25,k);

end
for k=1:ndimx
  Uef7(k) = VeU(ntop+30,k);
  Vef7(k) = VeV(ntop+30,k);
  Wef7(k) = Wor(ntop+30,k);

end
for k=1:ndimx
  Uef8(k) = VeU(ntop+35,k);
  Vef8(k) = VeV(ntop+35,k);
  Wef8(k) = Wor(ntop+35,k);

end

subplot(2,2,1);contour(VeU,30);xlabel('U')  
subplot(2,2,2);contour(VeV,30);xlabel('V')  
%subplot(2,2,1);quiver(VeU,VeV)
subplot(2,2,3);contour(Wor,60);%xlabel('W')  
subplot(2,2,4);contour(PHI,30);xlabel('Phi')  
%subplot(2,2,3);;plot(xx2,Wef)  
%subplot(2,2,4);plot(xx2,Vef,xx2,Uef)
%pause(0.2);
%plot(xx2,Uef,xx2,Uef1,xx2,Uef2,xx2,Uef3,xx2,Uef4,xx2,Uef5,xx2,Uef6,xx2,Uef7,xx2,Uef8)
pause();
cla;
%pause(0.3);
%plot(xx2,Vef,xx2,Vef1,xx2,Vef2,xx2,Vef3,xx2,Vef4,xx2,Vef5,xx2,Vef6,xx2,Vef7,xx2,Vef8)
%pause(0.2);
%cla;


if(is==1047)
    
for jj=1:ndimx
    fprintf(saleu,'%f %f %f %f %f %f %f %f %f %f\n',xx2(jj),Uef(jj),Uef1(jj),Uef2(jj),Uef3(jj),Uef4(jj),Uef5(jj),Uef6(jj),Uef7(jj),Uef8(jj));
    fprintf(salev,'%f %f %f %f %f %f %f %f %f %f\n',xx2(jj),Vef(jj),Vef1(jj),Vef2(jj),Vef3(jj),Vef4(jj),Vef5(jj),Vef6(jj),Vef7(jj),Vef8(jj));
    fprintf(salew,'%f %f %f %f %f %f %f %f %f %f\n',xx2(jj),Wef(jj),Wef1(jj),Wef2(jj),Wef3(jj),Wef4(jj),Wef5(jj),Wef6(jj),Wef7(jj),Wef8(jj));
end
fclose(saleu);
fclose(salev);
fclose(salew);


end



end

fclose(fid1);
fclose(fid2);
fclose(fid3);
fclose(fid4);





 