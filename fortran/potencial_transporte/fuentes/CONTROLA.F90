SUBROUTINE CONTROL()
use def_solver
use def_constantes
use def_variables
implicit none
double precision :: x(nodpel),y(nodpel),z(nodpel),esm(nodpel,nodpel),ef(nodpel),adiag,sigma_el,qe
integer :: ns(nodpel),ns2d(2*nodpel)
INTEGER  ::INBWT,nbwt,i,ij,j,nb,kk,NCASE,inode,ipoin,jnode,ii,jj,j1,nicio,k,nno,npas,ne,INBWE,nbwE
double precision, allocatable :: ud(:),un(:)
integer, allocatable :: cont(:),ip(:),iu(:),iup(:),ju(:),iut(:),consim(:)
integer, allocatable :: cont2d(:),ip2d(:),iu2d(:),iup2d(:),ju2d(:),iut2d(:),consim2d(:)
integer:: unit_sist


! ARMO LAS ESTRUCTURAS LOGICAS SPARCE
ALLOCATE(AD(nnodes),IA(nnodes+1),CONT(nnodes+1),CX(nnodes+1),solucion(nnodes),rhs(nnodes),consim(nnodes))
ALLOCATE(AD2d(2*nnodes),IA2d(2*nnodes+1),CONT2d(2*nnodes+1),CX2d(2*nnodes+1),solucion2d(2*nnodes),rhs2d(2*nnodes),consim2d(2*nnodes))
allocate(carga(nnodes))



RHS=0.0
ad=0.0
solucion=0.0
carga=0.0
ia=0
cx=0
cont=0
consim=0
NE=nelements




if(nmode<=2) then

      INBWT=0
      NBWT =0
      DO KK=1,NE
        DO I=1,nodpel
          NS(I)=conect(KK,I)
        ENDDO
        DO  I=1,nodpel-1
		  IJ=I+1
          DO J=IJ,nodpel
            NB=IABS(NS(I)-NS(J))
            IF(NB.EQ.0) THEN 
			   WRITE(unit_cont,*) 'ELEMENTO  ',KK,' TIENE DOS NODOS IDENTICOS'
			   WRITE(6,*) 'ELEMENTO  ',KK,' TIENE DOS NODOS IDENTICOS'
			   STOP
			ENDIF    
            IF(NB.GT.NBWT) THEN
               INBWT=KK
               NBWT =NB
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      NBWT=NBWT+1

      WRITE(unit_cont,*) ' BANDWIDTH: ',NBWT,'  EN ELEMENTO  ',INBWT
   
   
! DETERMINO LA FORMA SIMBOLICA DEL PROBLEMA SPARCE. PARA ESO ENSAMBLO UN SISTEMA DE UNOS Y CEROS



    DO KK=1,nnodes+1
      CX(KK)=0.0
      CONT(KK)=0.0
      IA(KK) = 0
    ENDDO


    NCASE=0
    sigma_el=1.0
    qe=1.0
    IA(1)=1

    do nno=1,nnodes
      consim=0
      DO KK=1,nelements
         npas=0
         DO I=1,nodpel
           NS(I)=conect(KK,I)
           if(nno==ns(i)) npas=1
         ENDDO
         if(npas==1) then   
           ESM=1.0
 

           DO I=1,nodpel
              II=NS(I)
                DO J=1,nodpel
                  JJ=NS(J)+1-II
                  IF(JJ.GT.0) THEN
	                 if(nno==ns(i).and.nno/=ns(j)) consim(ns(j)) = consim(ns(j))+1          
                  ENDIF  
	            ENDDO
            ENDDO

          endif
       ENDDO

       do i=1,nnodes
          if(consim(i)/=0) then
             cont(nno)=cont(nno)+1
          endif
       enddo
        IA(nno+1) = IA(nno) + CONT(nno)
    enddo




    NONULL = IA(nnodes+1)-1 !- nnodes


    ALLOCATE (an(NONULL),ja(NONULL))
    an=0.0
    ja=0

    write(unit_cont,*) 'nonulos del sistema  ',nnodes,nonull,nonull*(1.0/real(nnodes))*(1/real(nnodes))

    do nno=1,nnodes
      consim=0
      DO KK=1,nelements
         npas=0
         DO I=1,nodpel
           NS(I)=conect(KK,I)
           if(nno==ns(i)) npas=1
         ENDDO
         if(npas==1) then   
           ESM=1.0
   
         DO I=1,nodpel
           II=NS(I)
           DO J=1,nodpel
              JJ=NS(J)+1-II
              IF(JJ.GT.0) THEN
                if(nno==ns(i).and.nno/=ns(j)) consim(ns(j)) = consim(ns(j))+1          
              ENDIF  
	       ENDDO
         ENDDO

       endif
      ENDDO
      do i=1,nnodes
      
          if(consim(i)/=0) then
              JA(IA(nno)+CX(nno)) =   i
              CX(nno)=CX(nno)+1
          endif
       enddo
    enddo

elseif(nmode==3) then



    DO KK=1,nnodes+1
      CX(KK)=0.0
      CONT(KK)=0.0
      IA(KK) = 0
    ENDDO


    NCASE=0
    sigma_el=1.0
    qe=1.0
    IA(1)=1

    do nno=1,nnodes
      consim=0
      DO KK=1,nelements
         npas=0
         DO I=1,nodpel
           NS(I)=conect(KK,I)
           if(nno==ns(i)) npas=1
         ENDDO
         if(npas==1) then   
           ESM=1.0
 

           DO I=1,nodpel
              II=NS(I)
                DO J=1,nodpel
!                  JJ=NS(J)+1-II
!                  IF(JJ.NE.0) THEN
                  IF(J.NE.I) THEN
	                 if(nno==ns(i).and.nno/=ns(j)) consim(ns(j)) = consim(ns(j))+1          
                  ENDIF  
	            ENDDO
            ENDDO

          endif
       ENDDO

       do i=1,nnodes
          if(consim(i)/=0) then
             cont(nno)=cont(nno)+1
          endif
       enddo
        IA(nno+1) = IA(nno) + CONT(nno)
    enddo




    NONULL = IA(nnodes+1)-1 !- nnodes


    ALLOCATE (an(NONULL),ja(NONULL))
    an=0.0
    ja=0

    write(unit_cont,*) 'nonulos del sistema  ',nnodes,nonull,nonull*(1.0/real(nnodes))*(1/real(nnodes))

    do nno=1,nnodes
      consim=0
      DO KK=1,nelements
         npas=0
         DO I=1,nodpel
           NS(I)=conect(KK,I)
           if(nno==ns(i)) npas=1
         ENDDO
         if(npas==1) then   
           ESM=1.0
   
         DO I=1,nodpel
           II=NS(I)
           DO J=1,nodpel
!              JJ=NS(J)+1-II
!              IF(JJ.NE.0) THEN
              IF(J.NE.i) THEN
                if(nno==ns(i).and.nno/=ns(j)) consim(ns(j)) = consim(ns(j))+1          
              ENDIF  
	       ENDDO
         ENDDO

       endif
      ENDDO
      do i=1,nnodes
      
          if(consim(i)/=0) then
              JA(IA(nno)+CX(nno)) =   i
              CX(nno)=CX(nno)+1
          endif
       enddo
    enddo

endif
! caso mecanico

RHS2d=0.0
ad2d=0.0
solucion2d=0.0
ia2d=0
cx2d=0
cont2d=0
consim2d=0
NE=nelements

      INBWE=0
      NBWE =0
      DO KK=1,NE
        DO I=1,nodpel
          NS2d(2*I-1)= 2*conect(KK,I)-1
          NS2d(2*I)  = 2*conect(KK,I)
        ENDDO
        DO  I=1,2*nodpel-1
		  IJ=I+1
          DO J=IJ,2*nodpel
            NB=IABS(NS2d(I)-NS2d(J))
            IF(NB.EQ.0) THEN 
			   WRITE(unit_cont,*) 'ELEMENTO  ',KK,' TIENE DOS NODOS IDENTICOS'
			   WRITE(6,*) 'ELEMENTO  ',KK,' TIENE DOS NODOS IDENTICOS'
			   STOP
			ENDIF    
            IF(NB.GT.NBWE) THEN
               INBWe=KK
               NBWE =NB
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      NBWE=2*NBWE+1

      WRITE(unit_cont,*) ' BANDWIDTH: ',NBWE,'  EN ELEMENTO  ',INBWE
   
   
! DETERMINO LA FORMA SIMBOLICA DEL PROBLEMA SPARCE. PARA ESO ENSAMBLO UN SISTEMA DE UNOS Y CEROS



DO KK=1,2*nnodes+1
  CX2d(KK)=0.0
  CONT2d(KK)=0.0
  IA2d(KK) = 0
ENDDO


NCASE=0
sigma_el=1.0
qe=1.0
IA2d(1)=1

do nno=1,2*nnodes
  consim2d=0
  DO KK=1,nelements
     npas=0
     DO I=1,nodpel
       NS2d(2*I-1)=2*conect(KK,I)-1
       NS2d(2*I)=2*conect(KK,I)
       if(nno==ns2d(2*i-1) .or. nno==ns2d(2*i) ) npas=1
     ENDDO
     if(npas==1) then   
       DO I=1,2*nodpel
            II=NS2d(I)
            DO J=1,2*nodpel
              JJ=NS2d(J)+1-II
              IF(JJ.GT.0) THEN
	             if(nno==ns2d(i) .and. nno/=ns2d(j)) consim2d(ns2d(j)) = consim2d(ns2d(j))+1          
              ENDIF  
	        ENDDO
        ENDDO

      endif
   ENDDO

   do i=1,2*nnodes
      if(consim2d(i)/=0) then
         cont2d(nno)=cont2d(nno)+1
      endif
   enddo
    IA2d(nno+1) = IA2d(nno) + CONT2d(nno)
enddo




NONULL2d = IA2d(2*nnodes+1) !- nnodes


ALLOCATE (an2d(NONULL2d),ja2d(NONULL2d))
an2d=0.0
ja2d=0

write(unit_cont,*) 'nonulos del sistema  ',2*nnodes,nonull2d,nonull2d*(1.0/real(2*nnodes))*(1/real(2*nnodes))

do nno=1,2*nnodes
  consim2d=0
  DO KK=1,nelements
     npas=0
     DO I=1,nodpel
       NS2d(2*I-1)=2*conect(KK,I)-1
       NS2d(2*I)=2*conect(KK,I)
       if(nno==ns2d(2*i-1) .or. nno==ns2d(2*i) ) npas=1
     ENDDO
     if(npas==1) then   
   
     DO I=1,2*nodpel
       II=NS2d(I)
       DO J=1,2*nodpel
          JJ=NS2d(J)+1-II
          IF(JJ.GT.0) THEN
            if(nno==ns2d(i).and.nno/=ns2d(j)) consim2d(ns2d(j)) = consim2d(ns2d(j))+1          
          ENDIF  
	   ENDDO
     ENDDO

   endif
  ENDDO
  
  do i=1,2*nnodes
     if(consim2d(i)/=0) then
          JA2d(IA2d(nno)+CX2d(nno)) = i
          CX2d(nno)=CX2d(nno)+1
     endif
  enddo

enddo


open(unit=unit_sist,file=archi_sistema)

write(unit_sist,*) nnodes+1
do kk=1,nnodes+1

   write(unit_sist,*) kk,ia(kk),cx(kk)

enddo

write(unit_sist,*) nonull
do kk=1,nonull

   write(unit_sist,*) kk,ja(kk)

enddo

write(unit_sist,*) 2*nnodes+1
do kk=1,2*nnodes+1

   write(unit_sist,*) kk,ia2d(kk),cx2d(kk)

enddo

write(unit_sist,*) nonull2d
do kk=1,nonull2d

   write(unit_sist,*) kk,ja2d(kk)

enddo



close(unit_sist)
RETURN
END

