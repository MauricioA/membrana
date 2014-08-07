#ifndef ARMADO_H_
#define ARMADO_H_

#include "declaraciones.h"

using namespace declaraciones;

class Armado {
public:
	static void armadoPoisson(Double2D pos[], double sigma, int nodpel, double esm[][MAXNPEL]);

	static void armadoTransporte(int nodpel, Double2D pos[], double sigma, double qe, double landa,
			double mu, double mas[], double sol[], double esm[][MAXNPEL], double ef[], double deltaT);

	static double determinante3(Double2D pos[], double b[], double c[]);

	static double iteracion4(int i, int j, int kGauss, Double2D pos[4],
			double phi[2*NGAUSS][4], double dphi[NDIM][2*NGAUSS][4], double phidX[NDIM][2*NGAUSS][4]);

private:
	static void armado3(Double2D pos[], double sigma, double esm[][MAXNPEL]);

	static void armado4(Double2D pos[], double sigma, double qe, bool transp, double landa, double mu,
			double esm[][MAXNPEL], double ef[MAXNPEL], double est[][4], double mas[]);

	static void armadoTransporte3(Double2D pos[], double sigma, double qe, double landa,
			double mu, double mas[], double sol[], double esm[][MAXNPEL], double ef[], double est[][MAXNPEL]);

};

#endif /* ARMADO_H_ */
