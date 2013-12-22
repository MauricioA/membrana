fi=fopen('results.3D','r+');

N=fscanf(fi,'%i',1);

for kk=1:N
   x(kk)=fscanf(fi,'%f',1);
   y(kk)=fscanf(fi,'%f',1);
   sol(kk)=fscanf(fi,'%f',1);
end
fclose(fi);


plot3(x,y,sol,'.');