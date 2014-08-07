fi=fopen('results.3D','r+');

%~ N=fscanf(fi,'%i',1);
%~ 
%~ for kk=1:N   
   %~ x(kk)=fscanf(fi,'%f',1);
   %~ y(kk)=fscanf(fi,'%f',1);
   %~ sol(kk)=fscanf(fi,'%f',1);
%~ end
%~ fclose(fi);

k=0;
while(feof(fi)==0)
   k=k+1;
   x(k)=fscanf(fi,'%f',1);
   y(k)=fscanf(fi,'%f',1);
   sol(k)=fscanf(fi,'%f',1);
end
fclose(fi);

plot3(x,y,sol,'.');
