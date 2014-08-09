#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <conio.h>

#define  ii  100   //eje y
#define  jj  100   //exe x
#define  TTOTAL 25

#define GRABAR 100
#define TOL   1e-6   /* min. valor res. convergencia */
#define error 1e-5   /* min. valor res. convergencia */

#define tita0 0.0   /* Generacion de protones en el Anodo (PAra cationes es 1-tita0)*/
#define titaP 1.0   /* Inicial de Cationes*/

#define  ndimen2 5000
#define  eta1    5
#define  eta2    5

//#define  VB    0.1
//#define  EY    0.2


int **masc();
long double  **matrix();
void free_matrix(long double **m);
void free_masc(int **m);

void inicializo();
void calculo();
void calcular_V(int i,int j);
float salida();

void DBM(int);
void DLA(int);

void gradiente();
int probabilidad();
void correlacion();
void imprimir(float);

//archivos de salida
FILE *anio,*cati,*proto,*pote,*veloU,*veloV,*vorti,*depos,*corri,*datos,*frentes,*dlaf,ca;

FILE *dbmf,*densf,*corref,*vxb1,*vxb2,*campox,*campoy;


/* definiciones de variables de agregacion*/
long double *grad,sum_gra,**corre;
int num_gra,**ngrad,nlimite,NMAX,ncontac;
float **etha;

/* definicion de variables de calculo */

const long double   k = 0.001;              /* Paso de tiempo (s) */
const long double   wt = 0.8;	 /* Parametro de relajacion SOR */
const long double   L=1;                 /* Largo de la celda (cm)*/
const long double   c_0 = 1E9;           /* Concentracion inicial (1 /cm3)*/
const long double   pot_0 = 1.0;         /* Potencial inicial */
const long double   pot_an = 1.0;         /* Potencial real */
const long double   B_0 = -100000.0;         /* campo magnetico */
const long double   epsilon = 80.0;	 /* Constante dielectrica del agua */
const long double   epsilon_0 =0.07953983; /* Constante de permitividad (uc2 / dyna cm2) */
const long double   za = 2.0;           /* Carga de los aniones (-e) */
const long double   zc = 2.0;	          /* Carga de los cationes (e) */
const long double   zp = 1.0;	          /* Carga de los protones (e) */

const long double   kb = 1.3805E-16;/* Constante de Boltzmann (dyna cm / K)*/
const long double   t = 293.0;          /* Temperatura (Kelvin) */
const long double   e = 4.8E-10;	/* Carga del proton (uc)*/
const long double   Da = 1E-5;          /* Coef de dif aniones (cm2 / s)*/
const long double   Dc = 7E-6;          /* Coef de  "  cationes */
const long double   Dp = 4.6E-5;          /* Coef de  "  protones */
const long double   mua = 8.29E-4;	/* Mov ionica sulfato (cm2 / V s)*/
const long double   muc = 5.37E-4;	/* Mov ionica cobre (cm2 / V s) */
const long double	  mup= 36.23E-4;	/* Mov ionica proton (cm2 / V s) */
const long double   cte= 299.792458; /* constante de proporcionalidad cgs*/
const long double   rho = 1;         /* dens agua */
const long double   nu = 0.01;           /* visc. cinematica fluido cm^2/seg*/ //(nu=0.01== Re=10, nu=0.001==Re=100)
const long double   g=980.665;           /* aceleracion gravedad cm/seg^2*/
const long double   ALPHA=79.8;           /* coef. tipo Bousinesq aniones cm^3*/
const long double   BETA=79.8;           /* coef. tipo Bousinesq cationes cm^3*/
const long double   BETAP=79.8;           /* coef. tipo Bousinesq protones cm^3*/

int **mask;
long double **a0, **a1, **c0, **c1, **pr0, **pr1, **pot, **w0, **w1, **psi, **uu, **vv, **campoEx, **campoEy,**permi;

unsigned long indx=0;
long double F, h, Poa, Poc, Pop, rsa, rsc, rsp, rma, rmc, rmp, rescon, Fra, Frc, Frp;
float tiempo=0;

int i,j,ind,count,m,npart;

long double x_0,E_0,u_0,mu;    /* ahora las escalas se calculan */

/*-- definicion de numeros adimensionales ------------------*/
long double Re,Pea,Pec,Pep,Ma,Mc,Mp,Gra,Grc,Grp,GraB,GrcB,GrpB;
long double cbA,cbC,cbP;

long double Denn=0, Rtotal=0, RtotalP=0, RtotalA=0, RtotalC=0, resistencia=0,I_i=0;

main()
{
//long double t_incuba,t_agrega;
float Frente;
int det=0;

/* tamanio de la celda en direccion x (largo)*/
x_0=L;

/* definicion de paso espacial*/
h = x_0 / jj;

/* Otras constantes utiles */
F = e/epsilon_0;

u_0=0.1;

Poa = epsilon * pot_0 / x_0 / x_0 / F / za / c_0 /cte;
Poc = epsilon * pot_0 / x_0 / x_0 / F / zc / c_0 /cte;
Pop = epsilon * pot_0 / x_0 / x_0 / F / zp / c_0 /cte;

Re=u_0*x_0/nu;

Pea=x_0*u_0/Da;
Pec=x_0*u_0/Dc;
Pep=x_0*u_0/Dp;
Ma=u_0*x_0 /mua /pot_0;
Mc=u_0*x_0 /muc /pot_0;
Mp=u_0*x_0 /mup /pot_0;


//Grasshof magnetico

GraB=x_0*c_0*e*za*B_0/(u_0  * rho * cte) ;
GrcB=x_0*c_0*e*zc*B_0/(u_0  * rho * cte) ;
GrpB=x_0*c_0*e*zp*B_0/(u_0  * rho * cte );


rsa = k / h / h / Pea;
rma = k / h / h / Ma;
rsc = k / h / h / Pec;
rmc = k / h / h / Mc;
rsp = k / h / h / Pep;
rmp = k / h / h / Mp;

Fra = u_0 * u_0 * rho * cte / pot_0 / c_0 / e / za ;
Frc = u_0 * u_0 * rho * cte / pot_0 / c_0 / e / zc ;
Frp = u_0 * u_0 * rho * cte / pot_0 / c_0 / e / zp ;

// numero de grassohf gravitatorio
/*
Gra=x_0*c_0*g*ALPHA /u_0 /u_0 /cte;
Grc=x_0*c_0*g*BETA  /u_0 /u_0 /cte;
Grp=x_0*c_0*g*BETAP /u_0 /u_0 /cte;
*/


Gra=1/Fra;
Grc=1/Frc;
Grp=1/Frp;


/* aloco las matrices de potencial */
mask=masc();
pot=matrix();
c0=matrix();
c1=matrix();
a0=matrix();
a1=matrix();
pr0=matrix();
pr1=matrix();
w0=matrix();
w1=matrix();
psi=matrix();
uu=matrix();
vv=matrix();
campoEx=matrix();
campoEy=matrix();
permi=matrix();

/* aloco vectores del modelo de agregacion */
   corre=(long double**) malloc((unsigned) ((ii+1)*sizeof(long double*)));
	for(i=0;i<=ii;i++)
    	corre[i]=(long double*) malloc((unsigned) ((jj+1)*sizeof(long double)));

   etha=(float**) malloc((unsigned) ((ii+1)*sizeof(float*)));
	for(i=0;i<=ii;i++)
    	etha[i]=(float*) malloc((unsigned) ((jj+1)*sizeof(float)));

   ngrad=(int**) malloc((unsigned) ((ndimen2+1)*sizeof(int*)));
	for(i=0;i<=ndimen2;i++)
    	ngrad[i]=(int*) malloc((unsigned) (3*sizeof(int)) );

   grad=(long double*) malloc((unsigned) ((ndimen2+1)*sizeof(long double)));

//crea archivos de salida
   if( (anio  = fopen("anion.dat","w+"))==NULL ||
       (cati  = fopen("catio.dat","w+"))==NULL ||
       (pote  = fopen("poten.dat","w+"))==NULL ||
       (proto  = fopen("proto.dat","w+"))==NULL ||
       (depos = fopen("deposito.dat","w+"))==NULL ||
       (datos = fopen("datos.dat","w+"))==NULL ||
       (veloU = fopen("veloU.dat","w+"))==NULL ||
       (veloV = fopen("veloV.dat","w+"))==NULL ||
       (vorti= fopen("vortic.dat","w+"))==NULL||
       (corri= fopen("corrien.dat","w+"))==NULL  ||
       (frentes= fopen("frentes.dat","w+"))==NULL  ||
       (dbmf = fopen("dbm.dat","w+"))       ==NULL ||
       (densf  = fopen("densidad.dat","w+"))==NULL ||
       (dlaf  = fopen("dla.dat","w+"))==NULL ||
       (vxb1  = fopen("vxb1.dat","w+"))==NULL ||
       (vxb2  = fopen("vxb2.dat","w+"))==NULL ||
       (corref = fopen("correla.dat","w+")) ==NULL  ||
       (campox = fopen("campoEx.dat","w+")) ==NULL   ||
       (campoy = fopen("campoEy.dat","w+")) ==NULL)
     abort();

// para el agregado
//    randomize();

/*    for(j=0;j<=jj;j++)
    {
        for(i=ii*3/5;i<=ii;i++)
        {
                etha[i][j]=eta1;
        }
        for(i=0;i<ii*3/5;i++)
        {
                etha[i][j]=eta2;
        }
    }
*/

/* condiciones iniciales */
   inicializo();

   fprintf(datos,"Froude: %Lf %Lf %Lf\n",Fra,Frc,Frp);
   fprintf(datos,"Grasssofh: %Lf %Lf %Lf\n",Gra,Grc,Grp);
   fprintf(datos,"Grasssofh magneticos: %Lf %Lf %Lf\n",GraB,GrcB,GrpB);
   fprintf(datos,"Poisson: %Lf %Lf %Lf\n",Poa,Poc,Pop);
   fprintf(datos,"Migra: %Lf %Lf %Lf\n",Ma,Mc,Mp);
   fprintf(datos,"Peclet: %Lf %Lf %Lf\n",Pea,Pec,Pep);
   fprintf(datos,"Reynolls: %Lf\n",Re);

   fclose(datos);
/* grabo condiciones iniciales*/
   salida();

   printf("Inicio: \n");

   rescon=1;
   nlimite=50;
   NMAX = 20;
   ncontac=0;
   Frente=1;
   while(tiempo<=TTOTAL)// && rescon>TOL)
	{

       tiempo=tiempo+k;
	    indx++;


       calculo();

// aca tengo que calcular resistencia, corriente y posicion del frente anodico de protones

       Rtotal=RtotalC+RtotalP+RtotalA;
       Denn = (zc*muc*RtotalC+za*mua*RtotalA+zp*mup*RtotalP);
       resistencia = x_0/Denn;

//       I_i=pot_0/resistencia;
       I_i=pot_an/resistencia;

// Cuando el frente alcanza el catodo, comienza DBM

/*       if(Frente < 1)
       {
          det++;
//          if(det==(GRABAR*0.5))
          if(det==GRABAR)
          {
//             DBM(nlimite);
             DLA(nlimite);
             ncontac=ncontac+nlimite;
             imprimir(tiempo);
             printf("fractal Tiempo:%f  indice: %i ncontact: %i Frente: %f\n",tiempo, det, ncontac, Frente);
             det=0;
          }
       }
*/
       if(fmod(indx,GRABAR)==0)
       {
	       Frente = salida();
          printf("tiempo: %f  indice: %i  error: %Lf Frente: %f\n",tiempo,indx,rescon,Frente);
       }

       //if(ncontac>5000)
       //    tiempo=TTOTAL+1;

	}

   //correlacion();

   printf("tiempo: %Lf  indice: %i   error: %Lf\n",tiempo,indx,rescon);
   salida();
   printf("fin!\n");

/* desaloco todos los vectores y matrices
free_masc(mask);
free_matrix(pot);
free_matrix(c0);
free_matrix(c1);
free_matrix(a0);
free_matrix(a1);
free_matrix(pr0);
free_matrix(pr1);
free_matrix(psi);
free_matrix(w1);
free_matrix(w0);
free_matrix(uu);
free_matrix(vv);
*/

//fcloseall();
return 0;
}

// calculo!!!
void calculo ()

{
static	long double  aux_a,aux_c,aux_p,aux_v,aux_w,aux_psi;
static	long double r, p1, p2, p3, p4, delta, pot1, pot2, pot3, pot4,
		   res,err_tot,resa,resc,resp,respot,k1,k2,qt;
static  long double resw,respsi,q,alfa,beta,betap,gamma,gammap,eta,gammaG,
                    gammapG,etaG;
static	int i, j, ind, imas,imenos,jmas,jmenos,iimax;

RtotalA=0;
RtotalC=0;
RtotalP=0;

qt=k/h;
err_tot=2*error;
ind=1;



/*	Iteraciones sucesivas    */
while (err_tot > error)
{

/*	Calculo SOR de omega y Gauss-Seidel de phi */
   err_tot = 0.0;

//   iimax=NMAX;
//   if(tiempo<0.2)
      iimax=0;

	for (i=iimax; i<=ii; i++)
	{
   	for (j = 0; j<=jj; j++)
		{
        imas=i+1;
        imenos=i-1;
        jmas=j+1;
        jmenos=j-1;

        respot=0.0;
        resa=0.0;
        resc=0.0;
        resp=0.0;
        resw=0;
        respsi=0;
        aux_a  = a1[i][j];
        aux_v  = pot[i][j];
     	  aux_c  = c1[i][j];
   	  aux_p  = pr1[i][j];
	     aux_w  = w1[i][j];
	     aux_psi= psi[i][j];

//        if(mask[i][j]==-1) pot[i][j]=0;

        if( mask[i][j]!=-1 )
        {

/*	Calculo del potencial electrico  */
//	       if(i==ii) pot[i][j]=pot_0;//1-kb*t*cte/e*log(c1[i][j]); /*anodo*/
	       if(i==ii) pot[i][j]=pot_an;//1-kb*t*cte/e*log(c1[i][j]); /*anodo*/
		    else if (j==0) pot[i][j]=(4*pot[i][j+1]-pot[i][j+2])/3;/*BORDE dereho*/
          else if (j==jj) pot[i][j]=(4*pot[i][j-1]-pot[i][j-2])/3;/*BORDE izquierdo*/
	       else if(i==0) pot[i][j]=0;//-kb*t*cte/e*log(c1[i][j]); /*catodo*/
		    else
		    {
             if(mask[imas][j]==-1) pot[imas][j]= 0.0;//(4*pot[i-1][j]-pot[i-2][j])/3;
             if(mask[imenos][j]==-1) pot[imenos][j]= 0.0;//(4*pot[i+1][j]-pot[i+2][j])/3;
             if(mask[i][jmas]==-1) pot[i][jmas]= 0.0;//(4*pot[i][j-1]-pot[i][j-2])/3;
             if(mask[i][jmenos]==-1) pot[i][jmenos]=0.0;// (4*pot[i][j+1]-pot[i][j+2])/3;

             pot[i][j]=0.25*(pot[i+1][j]+pot[i][j+1]+ pot[i-1][j]+ pot[i][j-1]) - h*h/4.0 *(a1[i][j]/Poa - c1[i][j]/Poc - pr1[i][j]/Pop);
             pot[i][j]=(1.0 - wt) * aux_v + wt * pot[i][j];
          }

			 res=fabs(pot[i][j]-aux_v);
		    respot=(res > respot) ? res: respot;

/*	Calculo de aniones */

          if (i==ii)/*anodo*/
		    {
		       k1=(3*pot[i][j]-4*pot[i-1][j]+pot[i-2][j])/(Re*Ma);
         	 k2=1/(Pea*Re);
		       a1[i][j]=k2/(3*k2-k1)*(4*a1[i-1][j]-a1[i-2][j]);
		       a1[i][j]=wt*a1[i][j]+(1-wt)*aux_a;
	       }
		    else if(j==0)   a1[i][j]=(4*a1[i][j+1]-a1[i][j+2])/3;/*borde derecho*/
		    else if (j==jj) a1[i][j]=(4*a1[i][j-1]-a1[i][j-2])/3;/*borde izquierdo*/
          else if (i==0)/*catodo*/
		    {
		   	 k1=(3*pot[i][j]-4*pot[i+1][j]+pot[i+2][j])/(Re*Ma);
				 k2=1/(Pea*Re);
				 a1[i][j]=k2/(3*k2-k1)*(4*a1[i+1][j]-a1[i+2][j]);
				 a1[i][j]=wt*a1[i][j]+(1-wt)*aux_a;
  				 a1[i][j]=a1[i][j];
	       }
          /*primeros vecinos al agregado (catodo)*/
		    else
			 {
 			   if(mask[imas][j]==-1)
			   {
			     k2=1/(Pea);
			     k1=(3*pot[i+1][j]-4*pot[i][j]+pot[i-1][j])/(Ma);

	           if (imas+1==ii) uu[i+1][j]=0;
              else uu[i+1][j]=(psi[i+1][j+1]-psi[i+1][j-1])/2/h;

			     a1[i+1][j]=(4*a1[i][j]-a1[i-1][j])*k2/(3*k2-k1-uu[i+1][j]);
			     a1[i+1][j]=wt*a1[i+1][j]+(1-wt)*aux_a;
			   }
			   if (mask[imenos][j]==-1)
			   {
			     k2=1/(Pea);
			     k1=(3*pot[i-1][j]-4*pot[i][j]+pot[i+1][j])/(Ma);

	           if (imenos-1==0) uu[i-1][j]=0;
              else uu[i-1][j]=(psi[i-1][j+1]-psi[i-1][j-1])/2/h;

			     a1[i-1][j]=(4*a1[i][j]-a1[i+1][j])*k2/(3*k2-k1+uu[i-1][j]);
			     a1[i-1][j]=wt*a1[i-1][j]+(1-wt)*aux_a;
			   }
			   else if (mask[i][jmenos]==-1 && i > 0 )
			   {
			     k2=1/(Pea);
			     k1=(3*pot[i][j-1]-4*pot[i][j]+pot[i][j+1])/(Ma);

			     if (jmenos==0) vv[i][j-1]=0;
              else  vv[i][j-1]=-(psi[i+1][j-1]-psi[i-1][j-1])/2/h;

			     a1[i][j-1]=(4*a1[i][j]-a1[i][j+1])*k2/(3*k2-k1+vv[i][j-1]);
			     a1[i][j-1]=.1*a1[i][j-1]+(1-.1)*aux_a;
			   }
			   else if (mask[i][jmas]==-1 && i > 0)
			   {
			     k2=1/(Pea);
			     k1=(3*pot[i][j+1]-4*pot[i][j]+pot[i][j-1])/(Ma);

	           if (jmas==jj) vv[i][j+1]=0;
              else vv[i][j+1]=-(psi[i+1][j+1]-psi[i-1][j+1])/2/h;

			     a1[i][j+1]=(4*a1[i][j]-a1[i][j-1])*k2/(3*k2-k1-vv[i][j+1]);
			     a1[i][j+1]=.1*a1[i][j+1]+(1-.1)*aux_a;
			   }
            pot1=pot[i+1][j]-pot[i][j];
	         pot2=pot[i][j]-pot[i-1][j];
			   pot3=pot[i][j+1]-pot[i][j];
			   pot4=pot[i][j]-pot[i][j-1];
			   delta=(1+4*rsa+rma*pot1/2-rma*pot2/2+rma*pot3/2-rma*pot4/2);
			   calcular_V(i,j);
			   p1=(rsa-rma*pot1*0.5-qt*0.5*uu[i+1][j])/delta;
			   p2=(rsa-rma*pot3*0.5-qt*0.5*vv[i][j+1])/delta;
			   p3=(rsa+rma*pot2*0.5+qt*0.5*uu[i-1][j])/delta;
			   p4=(rsa+rma*pot4*0.5+qt*0.5*vv[i][j-1])/delta;
			   a1[i][j]=p1*a1[i+1][j]+p2*a1[i][j+1]+p3*a1[i-1][j]+p4*a1[i][j-1]+a0[i][j]/delta;
   		   a1[i][j]=(1.0 - wt) * aux_a + wt * a1[i][j];
          }
			 res=fabs(a1[i][j]-aux_a);
			 resa=(res > resa) ? res: resa;

/*	Calculo de cationes */

		    if(i==ii) c1[i][j]=(1-tita0)*a1[i][j];   /* anodo */
		    else if (j==0)  c1[i][j]=(4*c1[i][j+1]-c1[i][j+2])/3;  /*borde derecho*/
		    else if (j==jj) c1[i][j]=(4*c1[i][j-1]-c1[i][j-2])/3;/*borde izquierdo*/
          else if (i==0) c1[i][j]=(4*c1[i+1][j]-c1[i+2][j])/3;
      /*primeros vecinos al agregado (catodo)*/
          else
			 {
//			   if(mask[imas][j]==-1)   c1[i+1][j]=c1[i][j];
//			   else if(mask[imenos][j]==-1) c1[i-1][j]=c1[i][j];
//			   else if(mask[i][jmenos]==-1) c1[i][j-1]=c1[i][j];
//			   else if(mask[i][jmas]==-1)   c1[i][j+1]=c1[i][j];
			   if(mask[imas][j]==-1) c1[imas][j]=c1[i-1][j];
//			   if(mask[imas][j]==-1) c1[i][j]=(4*c1[i-1][j]-c1[i-2][j])/3;
			   else if(mask[imenos][j]==-1) c1[imenos][j]=(4*c1[i+1][j]-c1[i+2][j])/3;
			   else if(mask[i][jmenos]==-1) c1[i][jmenos]=(4*c1[i][j+1]-c1[i][j+2])/3;
			   else if(mask[i][jmas]==-1)   c1[i][jmas]=(4*c1[i][j-1]-c1[i][j-2])/3;
			   pot1=pot[i+1][j]-pot[i][j];
			   pot2=pot[i][j]-pot[i-1][j];
			   pot3=pot[i][j+1]-pot[i][j];
			   pot4=pot[i][j]-pot[i][j-1];
			   delta=(1+4*rsc-rmc*pot1*0.5+rmc*pot2*0.5-rmc*pot3*0.5+rmc*pot4*0.5);
			   p1=(rsc+rmc*pot1*0.5-qt*0.5*uu[i+1][j])/delta;
			   p2=(rsc+rmc*pot3*0.5-qt*0.5*vv[i][j+1])/delta;
			   p3=(rsc-rmc*pot2*0.5+qt*0.5*uu[i-1][j])/delta;
			   p4=(rsc-rmc*pot4*0.5+qt*0.5*vv[i][j-1])/delta;

			   c1[i][j]=p1*c1[i+1][j]+p2*c1[i][j+1]+p3*c1[i-1][j]+p4*c1[i][j-1]+c0[i][j]/delta;
   		   c1[i][j] = (1.0 - wt) * aux_c + wt * c1[i][j];
			 }

			 res=fabs(c1[i][j]-aux_c);
			 resc=(res > resc) ? res: resc;

/*	Calculo de protones */
		    if(j==0) pr1[i][j]=(4*pr1[i][j+1]-pr1[i][j+2])/3;  /*borde derecho*/
		    else if (i==ii) pr1[i][j]=tita0*a1[i][j];   /* anodo */
		    else if (j==jj) pr1[i][j]=(4*pr1[i][j-1]-pr1[i][j-2])/3;/*borde izquierdo*/
		    else if(i==0) pr1[i][j]=(4*pr1[i+1][j]-pr1[i+2][j])/3;  /*catodo*/
		    else
			 {
//			   if(mask[imas][j]==-1) pr1[i][j]=(4*pr1[i-1][j]-pr1[i-2][j])/3;
			   if(mask[imas][j]==-1) pr1[imas][j]=pr1[i-1][j];
			   else if(mask[imenos][j]==-1) pr1[imenos][j]=(4*pr1[i+1][j]-pr1[i+2][j])/3;
			   else if(mask[i][jmenos]==-1) pr1[i][jmenos]=(4*pr1[i][j+1]-pr1[i][j+2])/3;
			   else if(mask[i][jmas]==-1)   pr1[i][jmas]=(4*pr1[i][j-1]-pr1[i][j-2])/3;
			     pot1=pot[i+1][j]-pot[i][j];
			     pot2=pot[i][j]-pot[i-1][j];
			     pot3=pot[i][j+1]-pot[i][j];
			     pot4=pot[i][j]-pot[i][j-1];
			     delta=(1+4*rsp-rmp*pot1/2+rmp*pot2/2-rmp*pot3/2+rmp*pot4/2);
			     p1=(rsp+rmp*pot1/2-qt/2*uu[i+1][j])/delta;
			     p2=(rsp+rmp*pot3/2-qt/2*vv[i][j+1])/delta;
			     p3=(rsp-rmp*pot2/2+qt/2*uu[i-1][j])/delta;
			     p4=(rsp-rmp*pot4/2+qt/2*vv[i][j-1])/delta;

			     pr1[i][j]=p1*pr1[i+1][j]+p2*pr1[i][j+1]+p3*pr1[i-1][j]+p4*pr1[i][j-1]+pr0[i][j]/delta;
			     pr1[i][j] = (1.0 - wt) * aux_p + wt * pr1[i][j];
			 }

			 res=fabs(pr1[i][j]-aux_p);
			 resp=(res > resp) ? res: resp;

/*	Calculo de la funcion de corriente */

		    if(i!=0 && i!=ii && j!=0 && j!=jj)
		    {
		       psi[i][j] = (w1[i][j]*h*h+psi[i+1][j]+psi[i-1][j]+psi[i][j+1]+psi[i][j-1])/4.0;
		       psi[i][j] = (1.0 - wt) * aux_psi + wt * psi[i][j];
		    }

		    res = fabs (psi[i][j] - aux_psi);
		    respsi = (res > respsi) ? res : respsi;

/*	Calculo de la vorticidad */

		    if(j==0) w1[i][j]=-2*psi[i][j+1] /h /h;
		    else if (j==jj) w1[i][j]=-2*psi[i][j-1]/h /h;
		    else if (i==ii) w1[i][j]=-2*psi[i-1][j] /h /h;
		    else if (i==0) w1[i][j]=-2*psi[i+1][j] /h /h;
		    else
			 {
			   if (mask[i][jmas] == -1)
				    w1[i][jmas] = - (2 * psi[i][j-1]) / h / h;
			   if (mask[imenos][j] == -1)
				    w1[imenos][j] = -(2 * psi[imas][j]) / h / h;
			   if (mask[i][jmenos] == -1)
				    w1[i][jmenos] = -(2 * psi[i][j]) / h / h;
			   if (mask[imas][j] == -1)
				    w1[imas][j] = -(2 * psi[i][j])/h/h;

	         r = k / h / h / Re;
			   alfa=((pot[i][j+1]-pot[i][j-1])*(a1[i+1][j]-a1[i-1][j])-
				     (pot[i+1][j] - pot[i-1][j]) * (a1[i][j+1] - a1[i][j-1]))/4.0/h/h;
			   beta=((pot[i+1][j]-pot[i-1][j])*(c1[i][j+1]-c1[i][j-1])-
				       (pot[i][j+1] - pot[i][j-1]) * (c1[i+1][j] - c1[i-1][j]))/4.0/h/h;
            betap=((pot[i+1][j]-pot[i-1][j])*(pr1[i][j+1]-pr1[i][j-1])-
				       (pot[i][j+1] - pot[i][j-1]) * (pr1[i+1][j] - pr1[i-1][j]))/4.0/h/h;

               q= k / h / 2;
               delta = 1.0 + 4.0 * r;
			      p1 = (r-q*uu[i+1][j])/delta;
			      p2 = (r-q*vv[i][j+1])/delta;
			      p3 = (r+q*uu[i-1][j])/delta;
			      p4 = (r+q*vv[i][j-1])/delta;


              	// side view
              // gamma=(c1[i+1][j]-c1[i][j])/h;
              // gammap=(pr1[i+1][j]-pr1[i][j])/h;
              // eta=(a1[i+1][j]-a1[i][j])/h;

	// LAs derivadas van en sentido j (celda horizontal top view)
 				   gamma=(c1[i][j+1]-c1[i][j])/h;
				  	gammap=(pr1[i][j+1]-pr1[i][j])/h;
				  	eta=(a1[i][j+1]-a1[i][j])/h;

	// LAs derivadas van en sentido i (celda vertical top view)
//               gammaG=(c1[i+1][j]-c1[i][j])/h;
//               gammapG=(pr1[i+1][j]-pr1[i][j])/h;
//               etaG=(a1[i+1][j]-a1[i][j])/h;


// terminos magneticos
                 cbA = (1./h)*(a1[i][j]*(uu[i+1][j]+vv[i][j+1]-uu[i][j]-vv[i][j]) + uu[i][j]*(a1[i+1][j]-a1[i][j])  + vv[i][j]*(a1[i][j+1]-a1[i][j]));
                 cbC = (1./h)*(c1[i][j]*(uu[i+1][j]+vv[i][j+1]-uu[i][j]-vv[i][j]) + uu[i][j]*(c1[i+1][j]-c1[i][j])  + vv[i][j]*(c1[i][j+1]-c1[i][j]));
                 cbP = (1./h)*(pr1[i][j]*(uu[i+1][j]+vv[i][j+1]-uu[i][j]-vv[i][j])+ uu[i][j]*(pr1[i+1][j]-pr1[i][j])+ vv[i][j]*(pr1[i][j+1]-pr1[i][j]));

           /* Termino electroconv.  */
               w1[i][j] = p1*w1[i+1][j]+p2*w1[i][j+1]+p3*w1[i-1][j]+p4*w1[i][j-1]+w0[i][j]/delta - k*alfa/delta/Fra + k*beta/delta/Frc+ k*betap/delta/Frp;
//  termino gravitovconvectivo
               w1[i][j] = w1[i][j] +k*Grc*gammaG/delta+k*Grp*gammapG/delta+k*Gra*etaG/delta;
//  termino Magnetoconvectivo
               w1[i][j] = w1[i][j] + k*GraB*cbA/delta - k*GrcB*cbC/delta - k*GrpB*cbP/delta;
			      w1[i][j] = wt*w1[i][j] + (1.0 - wt)*aux_w;

			      w1[i][j]=wt*w1[i][j]+(1.0 - wt)*aux_w;
          }

			 res = fabs (w1[i][j] - aux_w);
			 resw = (res > resw) ? res : resw;

		/* Actualizacion de residuos */
		    if (respot>err_tot) err_tot=respot;
		    if (resc>err_tot) err_tot=resc;
		    if (resp>err_tot) err_tot=resp;
		    if (resa>err_tot) err_tot=resa;
		    if (resw>err_tot) err_tot=resw;
		    if (respsi>err_tot) err_tot=respsi;

        } /* end if de mask*/
     } /* end for */
  } /* end for */

  ind=ind+1;

  if (ind>50000)
  {
    printf("no converge...salgo \n");
	 printf("error_tot= %f",err_tot);
    salida();
    getch();
	 exit(0);
  }
}  /* end while iteracion interna*/

/* Calculo residuo de convergencia*/
rescon=0.0;
for (i=0;i<=ii;i++)
  for(j=0;j<=jj;j++)
  {
   	if (fabs(c1[i][j]-c0[i][j])>rescon) rescon=fabs(c1[i][j]-c0[i][j]);
   	if (fabs(pr1[i][j]-pr0[i][j])>rescon) rescon=fabs(pr1[i][j]-pr0[i][j]);
		if (fabs(a1[i][j]-a0[i][j])>rescon) rescon=fabs(a1[i][j]-a0[i][j]);
		if (fabs(w1[i][j]-w0[i][j])>rescon) rescon=fabs(w1[i][j]-w0[i][j]);
  }

/* transferencia de valores */
  for (i = 0; i <= ii; i++)
  { for (j = 0; j <= jj; j++)
	 {
	 	w0[i][j] = w1[i][j];
	 	a0[i][j] = a1[i][j];
	 	c0[i][j] = c1[i][j];
	 	pr0[i][j]= pr1[i][j];
      RtotalA += a1[i][j]/(ii*jj);
      RtotalC += c1[i][j]/(ii*jj);
      RtotalP += pr1[i][j]/(ii*jj);
	 }
  }

}  /* end rutina calculo */


void calcular_V(int i,int j)
{
	if (i+1==ii) uu[i+1][j]=0;
	else uu[i+1][j]=(psi[i+1][j+1]-psi[i+1][j-1])/2/h;
	if (j+1==jj) vv[i][j+1]=0;
	else vv[i][j+1]=-(psi[i+1][j+1]-psi[i-1][j+1])/2/h;
	if (i-1==0) uu[i-1][j]=0;
	else uu[i-1][j]=(psi[i-1][j+1]-psi[i-1][j-1])/2/h;
	if (j-1==0) vv[i][j-1]=0;
	else  vv[i][j-1]=-(psi[i+1][j-1]-psi[i-1][j-1])/2/h;
}


void inicializo()
{
int i,j,j_medio,s,s2;

RtotalA=0;
RtotalC=0;
RtotalP=0;


/*	Condiciones iniciales */
for (i = 0; i <= ii; i++)
   for (j = 0; j <= jj; j++)
   {
   /* Condiciones iniciales uniformes */
   a0[i][j] = 1.0;
   a1[i][j] = 1.0;
   c0[i][j] = titaP;
   c1[i][j] = titaP;
   pr0[i][j] = 2*(1.0-titaP);
   pr1[i][j] = 2*(1.0-titaP);
//   pot[i][j] = i * pot_0 / ii;
   pot[i][j] = i * pot_an / ii;
   uu[i][j]=0;
   vv[i][j]=0;
   psi[i][j]=0;
   w0[i][j]=0;
   w1[i][j]=0;
   RtotalA=RtotalA+ a1[i][j];
   RtotalC=RtotalC+ c1[i][j];
   RtotalP=RtotalP+ pr1[i][j];
   }

   fprintf(dlaf,"%i %i\n",ii+1,jj+1);

   for(i=0;i<=ii;i++)
	  for (j=0;j<=jj;j++)
		 mask[i][j]=0;

  i=1;
  for (j=0;j<=jj;j++)
  {
       mask[0][j]=-1;
       fprintf(dlaf,"%i %i\n",i,j+1);
  }



   j_medio = jj/2;

///* esto es un spike recto ancho
for(i=1;i<20;i++)
{
   mask[i][j_medio]=-1;
   mask[i][j_medio+1]=-1;
   mask[i][j_medio-1]=-1;
   mask[i][j_medio+2]=-1;
   mask[i][j_medio-2]=-1;
   fprintf(dlaf,"%i %i\n",i+1,j_medio+1);
   fprintf(dlaf,"%i %i\n",i+1,j_medio+2);
   fprintf(dlaf,"%i %i\n",i+1,j_medio+3);
   fprintf(dlaf,"%i %i\n",i+1,j_medio);
   fprintf(dlaf,"%i %i\n",i+1,j_medio-1);
}
//*/


   /* esto es un spike recto
   for(i=1;i<20;i++)
   {
      mask[i][j_medio]=-1;
      fprintf(dlaf,"%i %i\n",i+1,j_medio+1);
   }
*/


/* esto son dos spikes uno junto a otro, uno de ellos mas grande , o iguales!
   for(i=1;i<18;i++)
   {
      mask[i][j_medio-2]=-1;
      mask[i][j_medio+2]=-1;
      fprintf(dlaf,"%i %i\n",i+1,j_medio+1);
   }
//   mask[18][j_medio+2]=-1;
//   mask[19][j_medio+2]=-1;
   mask[18][j_medio-2]=-1;
   mask[19][j_medio-2]=-1;

*/


/* esto es un spike recto con una pequeña V en la punta
   s=0;
   s2=0;
   for(i=1;i<15;i++)
   {
      if(i<10)
      {
        mask[i][j_medio-s]=-1;
      }
      else if(i>=10)
      {
        s2++;
        s++;
        mask[i][j_medio-s]=-1;
        mask[i][j_medio+s2]=-1;
      }
      fprintf(dlaf,"%i %i\n",i+1,j_medio+1);
   }
*/

/* esto son 7 iniciadores de spike
   for(i=1;i<10;i++)
   {
      mask[i][j_medio]=-1;
      mask[i][j_medio-5]=-1;
      mask[i][j_medio-10]=-1;
      mask[i][j_medio-15]=-1;
      mask[i][j_medio+5]=-1;
      mask[i][j_medio+10]=-1;
      mask[i][j_medio+15]=-1;
      fprintf(dlaf,"%i %i\n",i+1,j_medio+1);
      fprintf(dlaf,"%i %i\n",i+1,j_medio+5+1);
      fprintf(dlaf,"%i %i\n",i+1,j_medio+10+1);
      fprintf(dlaf,"%i %i\n",i+1,j_medio+15+1);
      fprintf(dlaf,"%i %i\n",i+1,j_medio-5+1);
      fprintf(dlaf,"%i %i\n",i+1,j_medio-10+1);
      fprintf(dlaf,"%i %i\n",i+1,j_medio-15+1);
   }
*/


/* esto es un spike inclinado
   s=0;
   for(i=1;i<20;i++)
   {
      if(fmod(i,2)==0) s++;
      mask[i][j_medio+s]=-1;   //hacia la derecha (B<0)
//      mask[i][j_medio-s]=-1; //hacia la izqueirda (B>0)
      fprintf(dlaf,"%i %i\n",i+1,j_medio+1);
   }
*/


/* esto son 16 spikes inclinados
   s=0;
   for(i=1;i<20;i++)
   {
      if(fmod(i,2)==0) s++;

      if(B_0>0.0)
      {
         mask[i][j_medio+s]=-1;   //hacia la derecha (B>0)
         mask[i][j_medio-5+s]=-1;
         mask[i][j_medio-10+s]=-1;
         mask[i][j_medio-15+s]=-1;
         mask[i][j_medio-20+s]=-1;
         mask[i][j_medio-25+s]=-1;
         mask[i][j_medio-30+s]=-1;
         mask[i][j_medio-35+s]=-1;
         mask[i][j_medio-40+s]=-1;
         mask[i][j_medio+5+s]=-1;
         mask[i][j_medio+10+s]=-1;
         mask[i][j_medio+15+s]=-1;
         mask[i][j_medio+20+s]=-1;
         mask[i][j_medio+25+s]=-1;
         mask[i][j_medio+30+s]=-1;
         mask[i][j_medio+35+s]=-1;
         fprintf(dlaf,"%i %i\n",i+1,j_medio-5+s+1);
         fprintf(dlaf,"%i %i\n",i+1,j_medio-10+s+1);
         fprintf(dlaf,"%i %i\n",i+1,j_medio-15+s+1);
         fprintf(dlaf,"%i %i\n",i+1,j_medio-20+s+1);
         fprintf(dlaf,"%i %i\n",i+1,j_medio-25+s+1);
         fprintf(dlaf,"%i %i\n",i+1,j_medio-30+s+1);
         fprintf(dlaf,"%i %i\n",i+1,j_medio-35+s+1);
         fprintf(dlaf,"%i %i\n",i+1,j_medio-40+s+1);
         fprintf(dlaf,"%i %i\n",i+1,j_medio+5+s+1);
         fprintf(dlaf,"%i %i\n",i+1,j_medio+10+s+1);
         fprintf(dlaf,"%i %i\n",i+1,j_medio+15+s+1);
         fprintf(dlaf,"%i %i\n",i+1,j_medio+20+s+1);
         fprintf(dlaf,"%i %i\n",i+1,j_medio+25+s+1);
         fprintf(dlaf,"%i %i\n",i+1,j_medio+30+s+1);
         fprintf(dlaf,"%i %i\n",i+1,j_medio+35+s+1);
      }
      else
      {mask[i][j_medio-s]=-1; //hacia la izqueirda (B<0)
      mask[i][j_medio-5-s]=-1;
      mask[i][j_medio-10-s]=-1;
      mask[i][j_medio-15-s]=-1;
      mask[i][j_medio-20-s]=-1;
      mask[i][j_medio-25-s]=-1;
      mask[i][j_medio-30-s]=-1;
      mask[i][j_medio-35-s]=-1;
      mask[i][j_medio+5-s]=-1;
      mask[i][j_medio+10-s]=-1;
      mask[i][j_medio+15-s]=-1;
      mask[i][j_medio+20-s]=-1;
      mask[i][j_medio+25-s]=-1;
      mask[i][j_medio+30-s]=-1;
      mask[i][j_medio+35-s]=-1;
      mask[i][j_medio+40-s]=-1;
      fprintf(dlaf,"%i %i\n",i+1,j_medio-5-s+1);
      fprintf(dlaf,"%i %i\n",i+1,j_medio-10-s+1);
      fprintf(dlaf,"%i %i\n",i+1,j_medio-15-s+1);
      fprintf(dlaf,"%i %i\n",i+1,j_medio-20-s+1);
      fprintf(dlaf,"%i %i\n",i+1,j_medio-25-s+1);
      fprintf(dlaf,"%i %i\n",i+1,j_medio-30-s+1);
      fprintf(dlaf,"%i %i\n",i+1,j_medio-35-s+1);
      fprintf(dlaf,"%i %i\n",i+1,j_medio+5-s+1);
      fprintf(dlaf,"%i %i\n",i+1,j_medio+10-s+1);
      fprintf(dlaf,"%i %i\n",i+1,j_medio+15-s+1);
      fprintf(dlaf,"%i %i\n",i+1,j_medio+20-s+1);
      fprintf(dlaf,"%i %i\n",i+1,j_medio+25-s+1);
      fprintf(dlaf,"%i %i\n",i+1,j_medio+30-s+1);
      fprintf(dlaf,"%i %i\n",i+1,j_medio+35-s+1);
      fprintf(dlaf,"%i %i\n",i+1,j_medio+40-s+1);
      }
   }

*/


/* esto es un spike inclinado con una pequeña V en la punta
   s=0;
   s2=0;
   for(i=1;i<15;i++)
   {
      if(fmod(i,2)==0) s++;
      mask[i][j_medio-s]=-1;

      if(i>=10)
      {
        s2++;
        mask[i][j_medio-s+s2]=-1;
      }
      fprintf(dlaf,"%i %i\n",i+1,j_medio+1);
   }
*/


}


float salida()
{
int i,j,kk,jmedio;
int pasoP1,pasoP2,pasoP3;
float frenteP1,frenteP2,frenteP3;

int mascara[ii+1][jj+1],imas,jmas,imenos, jmenos;


   fprintf(dbmf,"%f\n",tiempo);
   for(i=0;i<=ii;i++)
      for(j=0;j<=jj;j++)
          mascara[i][j] = mask[i][j];



   for(i=0;i<ii;i++)
   {
      imas=i+1;
      imenos=i-1;
      if(i==0) imenos = 0;
      for(j=0;j<=jj;j++)
      {
          jmas=j+1;
          if(j==jj) jmas=0;
          jmenos=j-1;
          if(j==0) jmenos=jj;

          if(mask[i][j]!= -1 && mask[imas][j]  ==-1) mascara[i][j]=1;
          if(mask[i][j]!= -1 && mask[imenos][j]==-1) mascara[i][j]=1;
          if(mask[i][j]!= -1 && mask[i][jmas]  ==-1) mascara[i][j]=1;
          if(mask[i][j]!= -1 && mask[i][jmenos]==-1) mascara[i][j]=1;
      }
   }

   for(i=0;i<=ii;i++)
   {
      for(j=0;j<=jj;j++)
          fprintf(dbmf,"%i\t",mascara[i][j]);
      fprintf(dbmf,"\n");
   }

for (i = 0; i <= ii; i++)
{
   for (j = 0; j <= jj; j++)
   {
       fprintf(pote,"%Lf\t",pot[i][j]);
       fprintf(anio,"%Lf\t",a1[i][j]);
       fprintf(cati,"%Lf\t",c1[i][j]);
       fprintf(proto,"%Lf\t",pr1[i][j]);
       fprintf(vorti,"%Lf\t",w1[i][j]);
       fprintf(corri,"%Lf\t",psi[i][j]);
       fprintf(veloU,"%Lf\t",uu[i][j]);
       fprintf(veloV,"%Lf\t",vv[i][j]);
       fprintf(vxb1,"%Lf\t",vv[i][j]*B_0);
       fprintf(vxb2,"%Lf\t",-uu[i][j]*B_0);
   }
   fprintf(pote,"\n");
   fprintf(anio,"\n");
   fprintf(cati,"\n");
   fprintf(proto,"\n");
   fprintf(vorti,"\n");
   fprintf(corri,"\n");
   fprintf(veloU,"\n");
   fprintf(veloV,"\n");
   fprintf(vxb1,"\n");
   fprintf(vxb2,"\n");
//   fprintf(depos,"\n");
}
//calculo y saco el campo electrico 
for (i = 1; i < ii; i++)
{
   for (j = 1; j < jj; j++)
   {
       campoEx[i][j]= - (pot[i][j]-pot[i-1][j])/h;
       campoEy[i][j]= - (pot[i][j]-pot[i][j-1])/h;
       
       fprintf(campox,"%Lf\t",campoEx[i][j]);
       fprintf(campoy,"%Lf\t",campoEy[i][j]);
      
   }
   fprintf(campox,"\n");
   fprintf(campoy,"\n");
   
}


  pasoP1=0;
  pasoP2=0;
  pasoP3=0;
  frenteP1=0;
  frenteP2=0;
  frenteP3=0;
  jmedio = (int) jj*0.5;

  for(kk=ii;kk>0;kk--)
  {
      if(a1[kk][jmedio] < 1.5 && pasoP1==0)
      {
         frenteP1=(float)kk/ii;
         pasoP1=1;
      }
      if(a1[kk][jmedio]< 1.05 && pasoP2==0)
      {
         frenteP2=(float)kk/ii;
         pasoP2=1;
      }
      if(a1[kk][jmedio] < 1.005 && pasoP3==0)
      {
         frenteP3=(float)kk/ii;
         pasoP3=1;
      }

  }

/*  for(kk=ii;kk>0;kk--)
  {
      if(pr1[kk][jmedio] > (2.0*(1.0 - titaP)*0.5) && pasoP1==0)
      {
         frenteP1=(float)kk/ii;
         pasoP1=1;
      }
      if(pr1[kk][jmedio]> (2.0*(1.0 - titaP)*0.05) && pasoP2==0)
      {
         frenteP2=(float)kk/ii;
         pasoP2=1;
      }
      if(pr1[kk][jmedio] > (2.0*(1.0 - titaP)*0.005) && pasoP3==0)
      {
         frenteP3=(float)kk/ii;
         pasoP3=1;
      }

  }
*/


fprintf(frentes,"%f %f %f %f %Lf %Lf\n",tiempo,frenteP1,frenteP2,frenteP3,resistencia,I_i);
//printf("%f %f %f %f %Lf %Lf\n",tiempo,frenteP1,frenteP2,frenteP3,resistencia,I_i);

  return (frenteP2);
}


long double **matrix()
{
	int i;
	long double **m;
	m=(long double**) malloc((unsigned) ((ii+1)*sizeof(long double*)));
	if (!m) printf("allocation failure 1 in matrix()");
	for(i=0;i<=ii;i++) {
	m[i]=(long double*) malloc((unsigned) ((jj+1)*sizeof(long double)));
	if(!m[i]) printf("allocation failure 2 in matrix()");
				 }
	return m;
}


int **masc()
{
	int i;
	int **m;
	m=( int **) malloc((int) ((ii+1)*sizeof( int*)));
	if (!m) printf("allocation failure 1 in matrix()");
	for(i=0;i<=ii;i++) {
	m[i]=( int*) malloc((int) ((jj+1)*sizeof( int)));
	if(!m[i]) printf("allocation failure 2 in matrix()");
				 }
	return m;

}


void free_masc( int **m)
{
int i;
for(i=ii;i>=0;i--) free((char*) (m[i]));
free((char*) (m));
}

void free_matrix(long double **m)
{
int i;
for(i=ii;i>=0;i--) free((char*) (m[i]));
}


/*Algo asi: en DBM se pegan nlimite particulas por vez. La subrutina va a calcular
la distibucion de flujo o de potencial a eleccion del usuario. luego el gradiente y luego
la probabilidad. Luego pega nlimite veces. despues regresa a recalcular todo.*/

void DBM(int nlimite)
{
int ntotal=1,proba=0;

    while(ntotal<=nlimite)
    {

        gradiente();
//        fprintf(datos,"..........Viene de Gradiente:%Lf  %i\n",sum_gra,num_gra);
//        printf("..........Viene de Gradiente:%Lf  %i\n",sum_gra,num_gra);
        proba=probabilidad();
//        fprintf(datos,"..........Viene de Probabilidad: %i\n",proba);
//        printf("..........Viene de Probabilidad: %i\n",proba);
        ntotal++;
    }

}


void gradiente()
{
int kx,ky,kxmas,kxmenos,kymas,kymenos;
long double A1,A2,A3,A4,flujx,flujy;


   num_gra=0;
   sum_gra=0;

   for(kx=0;kx<ndimen2;kx++)
     grad[kx]=0.0;
/*Gradiente calculado con el Potencial*/

   for(ky=0;ky<=jj;ky++)
   {
      for(kx=1;kx<ii;kx++)
      {

          kymas=ky+1;
          kymenos=ky-1;
          if(ky==0)
             kymenos=jj;
          if(ky==jj)
             kymas=0;

          kxmas  =kx+1;
          kxmenos=kx-1;

          if(mask[kx][ky]!=-1 && (mask[kx][kymas]==-1||mask[kx][kymenos]==-1||mask[kxmas][ky]==-1||mask[kxmenos][ky]==-1))
          {
              A1 = pot[kxmas][ky];
              A2 = pot[kx][kymas];
              A3 = pot[kx][kymenos];
              A4 = pot[kxmenos][ky];

              if(mask[kxmas][ky]==-1)   A1=0;
              if(mask[kx][kymas]==-1)   A2=0;
              if(mask[kx][kymenos]==-1) A3=0;
              if(mask[kxmenos][ky]==-1) A4=0;

              grad[num_gra] = pow(sqrt( (A1-A4)*(A1-A4)+(A2-A3)*(A2-A3) ),etha[kx][ky]);

              ngrad[num_gra][0]=kx;
              ngrad[num_gra][1]=ky;

//              printf("%i %Lf %i %i\n",num_gra,grad[num_gra],ngrad[num_gra][0],ngrad[num_gra][1]);

              sum_gra = sum_gra + grad[num_gra];
              num_gra=num_gra+1;

          }

      }
   }
/*fin del calculo utilizanto potencial */

/*Gradiente calculado con el flujo

   for(ky=0;ky<=jj;ky++)
   {
      for(kx=1;kx<ii;kx++)
      {
          kymas=ky+1;
          kymenos=ky-1;
          if(ky==0)
             kymenos=jj;
          if(ky==jj)
             kymas=0;

          kxmas=kx+1;
          kxmenos=kx-1;
          if(mask[kx][kymas]==-1||mask[kx][kymenos]==-1||mask[kxmas][ky]==-1||mask[kxmenos][ky]==-1)
          {

              A1 = mask[kxmas][ky];
              A2 = mask[kx][kymas];
              A3 = mask[kx][kymenos];
              A4 = mask[kxmenos][ky];

              if(A1==-1) A1=0;
              if(A2==-1) A2=0;
              if(A3==-1) A3=0;
              if(A4==-1) A4=0;

              flujx=-(A1-A4)/Mc/2/h-(c1[kxmas][ky]-c1[kxmenos][ky])/Pec/2/h+c1[kx][ky]*uu[kx][ky];
              flujy=-(A2-A3)/Mc/2/h-(c1[kx][kymas]-c1[kx][kymenos])/Pec/2/h+c1[kx][ky]*vv[kx][ky];

              if(A1==0 && flujx<0) flujx=0;
              if(A4==0 && flujx>0) flujx=0;
              if(A2==0 && flujy<0) flujy=0;
              if(A3==0 && flujy>0) flujy=0;

              grad[num_gra] = pow(sqrt(pow(flujx,2)+pow(flujy,2)),etha[kx][ky]);

              ngrad[num_gra][0]=kx;
              ngrad[num_gra][1]=ky;

              sum_gra = sum_gra + grad[num_gra];
              num_gra=num_gra+1;
          }
      }
   }

fin del calculo con flujo*/

}

int probabilidad()
{
int kx=0,paso=0,probabilidad;
double arand,probable=0;

    arand = rand()/(double)RAND_MAX;

    while(paso==0 && kx<= num_gra)
    {

        if(sum_gra==0)
           probable = probable + 1;
        else
           probable = probable + grad[kx]/sum_gra;

        if(probable>arand)
        {
          paso=1;
          mask[ngrad[kx][0]][ngrad[kx][1]]=-1;
          pot[ngrad[kx][0]][ngrad[kx][1]]=0;
          if(ngrad[kx][0]==ii)
              probabilidad = 0;
          else
              probabilidad = 1;
        }

        kx++;

    }

    fprintf(datos,"probabilidad: %f %f\n",probable,arand);
//    printf("probabilidad: %f %f\n",probable,arand);

return (probabilidad);
}


void correlacion()
{
int kx,ky,suma=0,j;
long double rho;
int AFUN(long double);

//   DENSIDAD DE DISTRIBUCON E PARTICULAS. INOHOMGENEA EN LA DIRECCION PEREPEDICULAR
//   AL SUSTRATO. ESTO SE ESTUDIA CALCULANDO EL NUMERO NORMALIZADO DE PARTICULAS A UNA
//   ALTURA H

   for(kx=0;kx<=ii;kx++)
   {
      for(ky=0;ky<=jj-1;ky++)
      {
          if(mask[kx][ky]==-1)
           suma=suma+1;
      }
      rho=(long double) suma/((ii+1)*(jj+1));
      fprintf(densf,"%i %Lf\n",kx,rho);
   }
//  DENSIDAD DE CORRELACION PARALELA AL DEPOSITO

   for(kx=0;kx<=ii+1;kx=kx+5)
   {
      for(ky=0;ky<=jj;ky++)
      {

         corre[kx][ky]=0;

         for(j=0;j<=jj;j++)
              if(ky+j <= jj) corre[kx][ky]=corre[kx][ky]+ AFUN(mask[kx][j])*AFUN(mask[kx][ky+j]);

         corre[kx][ky] = corre[kx][ky] / ((long double) jj+1);

         fprintf(corref,"%Lf\t",corre[kx][ky]);

      }
      fprintf(corref,"\n");
   }

}

int AFUN(long double xx)
{

   if(xx==-1)
     return 1;

return 0;
}

void imprimir(float tiempo)
{
int i,j;
int mascara[ii+1][jj+1],imas,jmas,imenos, jmenos;

   fprintf(dbmf,"%f\n",tiempo);
   for(i=0;i<=ii;i++)
      for(j=0;j<=jj;j++)
          mascara[i][j] = mask[i][j];

   for(i=1;i<ii;i++)
   {
      imas=i+1;
      imenos=i-1;
      for(j=1;j<jj;j++)
      {
          jmas=j+1;
          jmenos=j-1;
          if(mask[i][j]!= -1 && mask[imas][j]==-1  ) mascara[i][j]=1;
          if(mask[i][j]!= -1 && mask[imenos][j]==-1) mascara[i][j]=1;
          if(mask[i][j]!= -1 && mask[i][jmas]==-1  ) mascara[i][j]=1;
          if(mask[i][j]!= -1 && mask[i][jmenos]==-1) mascara[i][j]=1;
      }
   }


   for(i=0;i<=ii;i++)
   {
      for(j=0;j<=jj;j++)
          fprintf(dbmf,"%i\t",mascara[i][j]);
      fprintf(dbmf,"\n");
   }

/*   for(i=0;i<=ii;i++)
   {
      for(j=0;j<=jj;j++)
          fprintf(dbmf,"%i\t",mask[i][j]);
      fprintf(dbmf,"\n");
   }*/

}


void DLA(int nlimite)
{
float PI=3.14159,xran,zmov,VUP,VRIGHT,VLEFT,VDOWN,umax,VB,EY,nlado,Vtotal;
int k,NI,NJ,NTESTI,NTESTJ,NLIM,pegado,naz2;
double angulo;

int    funB(float,float,float,float,float);
double angul(double,double,float);

nlado=10;

for(k=1;k<=nlimite;k++)
{
   xran = rand()/(double)RAND_MAX;

   NI=ii;
   NJ=nlado + (jj-2*nlado)*xran;

   pegado=0;

   while(pegado==0)
   {
     zmov = rand()/(double)RAND_MAX;


     angulo = angul(uu[NI][NJ],vv[NI][NJ],PI);
     umax=fabs(uu[NI][NJ]);
     if(umax < fabs(vv[NI][NJ]))
     		  umax = fabs(vv[NI][NJ]);

     EY = (pot[NI][NJ]-pot[NI-1][NJ])/h;
     VB = mua*EY*B_0/fabs(mua*B_0); //0.05;
     //printf("vel  %i %i   %f %f\n",NI,NJ,EY,VB);
     VUP   = 0;
     VDOWN = 0.75;
     VLEFT = 0.25/2.;
     VRIGHT= 0.25/2.;

/*  //E
     VUP   = -EY*cos(angulo)+0.25;
     VDOWN = +EY*cos(angulo);
     VLEFT = -EY*sin(angulo);
     VRIGHT= +EY*sin(angulo);
*/
     // campo MAGNETICO
 //B>0

/*     VUP   = +VB*cos(angulo - PI/2.);//-EY*cos(angulo);
     VDOWN = -VB*cos(angulo - PI/2.);//+EY*cos(angulo);
     VLEFT = +VB*sin(angulo - PI/2.);//-EY*sin(angulo);
     VRIGHT= -VB*sin(angulo - PI/2.);//+EY*sin(angulo);
*/
//     VUP   = +VB*cos(angulo - PI/2.);
//     VDOWN = -VB*cos(angulo - PI/2.)+ EY;
//     VLEFT = +VB*sin(angulo - PI/2.);
//     VRIGHT= -VB*sin(angulo - PI/2.);

// B<0

//     VUP   = +VB*cos(angulo + PI/2.);//-EY*cos(angulo);
//     VDOWN = -VB*cos(angulo + PI/2.);//+EY*cos(angulo);
//     VLEFT = +VB*sin(angulo + PI/2.);//-EY*sin(angulo);
//     VRIGHT= -VB*sin(angulo + PI/2.);//+EY*sin(angulo);

//     VUP   =+VB*cos(angulo + PI/2.);
//     VDOWN =-VB*cos(angulo + PI/2.)+ EY;
//     VLEFT =+VB*sin(angulo + PI/2.);
//     VRIGHT=-VB*sin(angulo + PI/2.);

       //fprintf(datos,"K:%i %f %f %f %f %f\n",k,VUP,VDOWN,VLEFT,VRIGHT,VB);


      Vtotal= fabs(VUP)+fabs(VDOWN)+fabs(VLEFT)+fabs(VRIGHT);
     VUP   =      VUP/Vtotal;
     VDOWN =      VDOWN/Vtotal;
     VLEFT =      VLEFT/Vtotal;
     VRIGHT=      VRIGHT/Vtotal;

       naz2=funB(zmov,VUP,VDOWN,VLEFT,VRIGHT);

       NTESTI=NI;
       NTESTJ=NJ;

       switch(naz2)
       {
          case 0:
            NI=NI-1;
            break;
          case 1:
            NI=NI+1;
            if(NI>=ii)
                NI=ii-nlado;
            break;
          case 2:
            NJ=NJ+1;
            if(NJ >jj-nlado)
               NJ=nlado;
            break;
          case 3:
            NJ=NJ-1;
            if(NJ<nlado)
               NJ=jj-nlado;
            break;
       }

       if(mask[NI][NJ]==-1)
       {
             mask[NTESTI][NTESTJ]=-1;
             pegado=1;
             fprintf(dlaf,"%i %i\n",NTESTI+1,NTESTJ+1);

// BUSCO nLIMITE
            if(NTESTI>NMAX)
            {
              NMAX=NTESTI;
              printf("%i %i\n",NMAX,NTESTJ+1);
            }
       }
   }

}

}


int funB(float XR,float VU,float VD,float VL,float VR)
{
float RD,RU,RR,RL;
int NFUNB4;


RD=VD;
RU=VU;
RL=VL;
RR=VR;

/*
if(RD>0.5)
{
   RD=0.5;
   RU=0.0;
}

if(RU>0.5)
{
   RU=0.5;
   RD=0.0;
}
if(RL>0.5)
{
   RL=0.5;
   RR=0.0;
}
if(RR>0.5)
{
    RR=0.5;
    RL=0.0;
}
*/

if(XR>=0.0 && XR<RD)
    NFUNB4      =0;
if (XR>=RD &&  XR<(RD+RU))
    NFUNB4  =1;
if(XR>=(RD+RU)&& XR<(RD+RU+RR))
    NFUNB4 =2;
if(XR>=(RD+RU+RR) && XR<1.0)
    NFUNB4    =3;

return NFUNB4;

}

double angul(double VEX,double VEY,float PI)
{
double FANGUL;

if(VEX==0.0)
   FANGUL=PI/2.;

if(VEX>0.0 && VEY>0.0)
   FANGUL = atan(VEY/VEX);

if((VEX<0.0 && VEY> 0.0)|| (VEX>0.0 && VEY<=0.0))
   FANGUL = -atan(fabs(VEY/VEX));

if((VEX< 0.0 && VEY<= 0.0))
   FANGUL = atan(VEY/VEX);

return FANGUL;

}


