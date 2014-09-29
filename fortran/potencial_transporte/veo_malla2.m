function veo_malla2(archivo,nnod)

fid=fopen(archivo,'r+');

s=fscanf(fid,'%s',1);

N=fscanf(fid,'%i',1);

for kk=1:N
    
   nn=fscanf(fid,'%i',1);
   x(kk)=fscanf(fid,'%f',1);
   y(kk)=fscanf(fid,'%f\n',1);
    
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
    
    if(materia(kk)==2) 
       line(x1,y1);
       line(x2,y2);
       line(x3,y3);
       line(x4,y4);
   
    if(nnod==1)
      
        xno = x(i);
        yno = y(i);
        if(xno< 11 & yno >39 & yno < 61)
        s = sprintf('%i',i);
        text(xno,yno,s);
        end
        xno = x(j);
        yno = y(j);
        if(xno< 11 & yno >39 & yno < 61)
        s = sprintf('%i',j);
        text(xno,yno,s);
        end
        xno = x(k);
        yno = y(k);
        if(xno< 11 & yno >39 & yno < 61)
        s = sprintf('%i',k);
        text(xno,yno,s);
        end
        xno = x(ii);
        yno = y(ii);
        if(xno< 11 & yno >39 & yno < 61)
        s = sprintf('%i',ii);
        text(xno,yno,s);
        end 
    end
    
    end
    
    
end




