

fid=fopen('celula4.fem','r+');

s=fscanf(fid,'%s',1);

N=fscanf(fid,'%i',1);

for kk=1:N
    
   nn=fscanf(fid,'%i',1);
   x(kk)=fscanf(fid,'%f',1);
   y(kk)=fscanf(fid,'%f\n',1);
     
   Gx(kk) = 2*x(kk)*y(kk)+y(kk)*y(kk);
   Gy(kk) = 2*x(kk)*y(kk)+x(kk)*x(kk);
   
   
end

s=fscanf(fid,'%s',1);
ng=fscanf(fid,'%i',1);


ne=0;
for jj=1:ng
    
    nn=fscanf(fid,'%i',1);
    group(jj)=fscanf(fid,'%i',1);
    s=fscanf(fid,'%s\n',1);
    
    for ii=ne+1:ne+group(jj)
        materia(ii)=nn;
    end
    
    ne=ne+group(jj);
    
    
end

    s=fscanf(fid,'%s\n',1);


    
for kk=1:ne
    n=fscanf(fid,'%i',1);
    nel(kk,1)=fscanf(fid,'%i',1);
    nel(kk,2)=fscanf(fid,'%i',1);
    nel(kk,3)=fscanf(fid,'%i',1);
    nel(kk,4)=fscanf(fid,'%i\n',1);

    xm=0.0;
    ym=0.0;
    for ii=1:4
       xm=xm+x( nel(kk,ii))*0.25;
       ym=ym+y( nel(kk,ii))*0.25;
    end
    xme(kk)=xm;
    yme(kk)=ym;
    Gx_e(kk) = 2*xm*ym+ym*ym;
    Gy_e(kk) = 2*xm*ym+xm*xm;
    
end


fclose(fid);

fid=fopen('fort.333','r+');

while(feof(fid)==0)
   
    n=fscanf(fid,'%i',1);
    x_1(n)=fscanf(fid,'%f',1);
    y_1(n)=fscanf(fid,'%f',1);
    
    grx_ne(n)=fscanf(fid,'%f',1);
    gry_ne(n)=fscanf(fid,'%f\n',1);
    
    
end

fclose(fid)



fid=fopen('fort.444','r+');

while(feof(fid)==0)
   
    n=fscanf(fid,'%i',1);
    x_2(n)=fscanf(fid,'%f',1);
    y_2(n)=fscanf(fid,'%f',1);
    
    grx_nn(n)=fscanf(fid,'%f',1);
    gry_nn(n)=fscanf(fid,'%f\n',1);
    
    
end

fclose(fid)


plot3(x,y,Gx,x_2,y_2,grx_nn)

pause

plot3(x,y,Gy,x_2,y_2,gry_nn)

pause

plot3(xme,yme,Gx_e,x_1,y_1,grx_ne)

pause

plot3(xme,yme,Gy_e,x_1,y_1,gry_ne)





