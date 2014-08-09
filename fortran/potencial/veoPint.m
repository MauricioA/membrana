f=fopen('saleplano.dat','r+');

k=0;
while(feof(f)==0)
   k=k+1;
    x(k)=fscanf(f,'%f',1);
    y(k)=fscanf(f,'%f',1);
    V(k)=fscanf(f,'%f',1);
    C(k)=fscanf(f,'%f\n',1);
end

fclose(f);

plot3(x,y,V,'.')
pause
plot3(x,y,C,'.');xlabel('x');ylabel('y')
