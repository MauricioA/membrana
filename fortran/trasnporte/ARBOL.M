function arbol(caso)
% caso = 1 LAM ecd 
% caso = 2 celular

if(caso==1)
   fid = fopen('dla.dat','rb+');
   Nx = fscanf(fid,'%i',1)
   Ny = fscanf(fid,'%i',1)
elseif(caso==2)
   fid = fopen('fracmc.txt','rb+');
end

for jj=1:Nx+1
  for kk=1:Ny+1
      xx(kk,jj)=(kk-1)/(Nx+1);
      yy(kk,jj)=(jj-1)/(Ny+1);
  end
end

%VeU = fscanf(fid1,'%f\n',[Ny+1,Nx+1]);
%VeV = fscanf(fid2,'%f\n',[Ny+1,Nx+1]);

while(feof(fid)==0)
    [re,cont]=fscanf(fid,'%i %i\n');
end
dim=cont/2

for k=1:dim,
    x(k)=re(2*k-1);
    y(k)=re(2*k);
end    

fclose(fid);
whitebg('white') 

plot(y,x,'k.');axis([0 Ny 0 Nx])



 