function veo_malla(archivo,mate)

fid=fopen(archivo,'r+');
N=fscanf(fid,'%i',1);


nt=0;
np=0;
for kk=1:N
    n=fscanf(fid,'%i',1);
    x(kk)=fscanf(fid,'%f',1);
    y(kk)=fscanf(fid,'%f\n',1);
    
    if(y(kk)==0.0)
        nt=nt+1;
        tierra(nt)=kk;
    elseif(y(kk)==100.0)
        np=np+1;
        poten(np)=kk;
    end
end
nt
tierra
np
poten

ne=fscanf(fid,'%i',1)

for kk=1:ne
    n=fscanf(fid,'%i',1);
    nel(kk,1)=fscanf(fid,'%i',1);
    nel(kk,2)=fscanf(fid,'%i',1);
    nel(kk,3)=fscanf(fid,'%i',1);
    nel(kk,4)=fscanf(fid,'%i\n',1);
    mat(kk)=1;
    if(kk>=4761 & kk < 5013)
        mat(kk)=2;
    elseif(kk>=5013 )
        mat(kk)=3;
    end
end



fclose(fid);

%plot(x,y,'*');axis([0 100 0 100])

%pause
%cla;

for kk=1:ne
    i=nel(kk,1);
    j=nel(kk,2);
    k=nel(kk,3);
    ii=nel(kk,4);
    x1=[x(i),x(j)];
    y1=[y(i),y(j)];
    x2=[x(j),x(k)];
    y2=[y(j),y(k)];
    x3=[x(k),x(ii)];
    y3=[y(k),y(ii)];
    x4=[x(ii),x(i)];
    y4=[y(ii),y(i)];

    if(mat(kk)==mate) 
       line(x1,y1);
       line(x2,y2);
       line(x3,y3);
       line(x4,y4);
    end

end




