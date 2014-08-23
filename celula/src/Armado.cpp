#define _USE_MATH_DEFINES

#include "Armado.h"
#include <cassert>
#include <math.h>

using namespace std;

const double GAUSSPT[] = {
	-1 / sqrt(3.0),
	1 / sqrt(3.0),
};

const double GAUSSWT[] = { 1, 1 };

void Armado::armadoPoisson(Double2D pos[], double sigma, int nodpel, double esm[][MAXNPEL]) {
	switch (nodpel) {
	case 3:
		armado3(pos, sigma, esm);
		break;
	case 4:
		double ef[MAXNPEL];
		armado4(pos, sigma, 0, false, 0, 0, esm, ef, NULL, NULL);
		break;
	default: 
		assert(false);
	}
}

void Armado::armado3(Double2D pos[], double sigma, double esm[][MAXNPEL]) {
	const int NODPEL = 3;
	double b[3], c[3];
	double det = determinante3(pos, b, c);
	double rMed = (pos[0].x + pos[1].x + pos[2].x) / 3;

	assert(abs(det) > TOLER_AREA);

	for (int i = 0; i < NODPEL; i++) for (int j = 0; j < 3; j++) {
		esm[i][j] = sigma * (b[i] * b[j] + c[i] * c[j]) * M_PI * rMed / det;
	}
}

double Armado::determinante3(Double2D pos[], double b[], double c[]) {
	int i = 0;
	b[i++] = pos[1].y - pos[2].y;
	b[i++] = pos[2].y - pos[0].y;
	b[i++] = pos[0].y - pos[1].y;

	i = 0;
	c[i++] = pos[2].x - pos[1].x;
	c[i++] = pos[0].x - pos[2].x;
	c[i++] = pos[1].x - pos[0].x;

	return
		+ pos[1].x * pos[2].y + pos[2].x * pos[0].y + pos[0].x * pos[1].y
		- pos[1].x * pos[0].y - pos[2].x * pos[1].y - pos[0].x * pos[2].y;
}

void Armado::armado4(Double2D pos[], double sigma, double qe, bool transp, double landa, double mu,
		double esm[][MAXNPEL], double ef[MAXNPEL], double est[][4], double mas[]) {

	const int NODPEL = 4;
	double phi[2*NGAUSS][NODPEL];
	double dphi[NDIM][2*NGAUSS][NODPEL];
	double gxCod[NDIM][2*NGAUSS];
	double phidX[NDIM][2*NGAUSS][4];
	double cteI[2*NGAUSS];
	int kGauss = 0;

	for (int i = 0; i < NGAUSS; i++) for (int j = 0; j < NGAUSS; j++) {
		double det = iteracion4(i, j, kGauss, pos, phi, dphi, phidX);

		for (int dim = 0; dim < NDIM; dim++) {
			gxCod[dim][kGauss] = 0;

			for (int k = 0; k < NODPEL; k++) {
				gxCod[dim][kGauss] += pos[k].x * phi[kGauss][k];
			}
		}

		cteI[kGauss] = det * GAUSSWT[i] * GAUSSWT[j] * 2 * M_PI * gxCod[0][kGauss];
		kGauss++;
	}

	for (int i = 0; i < NODPEL; i++) for (int j = 0; j < NODPEL; j++) esm[i][j] = 0;
	for (int i = 0; i < NODPEL; i++) ef[i] = 0;

	for (int kGauss = 0; kGauss < (NGAUSS*NGAUSS); kGauss++) for (int i = 0; i < NODPEL; i++) {
		for (int j = 0; j < NODPEL; j++) {
			for (int d = 0; d < NDIM; d++) {
				esm[i][j] += sigma * phidX[d][kGauss][i] * phidX[d][kGauss][j] * cteI[kGauss];
			}

			if (transp) {	//upwing
				esm[i][j] += mu * phidX[1][kGauss][i] * phi[kGauss][j] * cteI[kGauss];
			}
		}

		ef[i] += cteI[kGauss] * phi[kGauss][i] * qe;
	}

	for (int i = 0; i < NODPEL; i++) {
		if (transp) {
			est[i][i] += mas[i] * landa;
		}
	}
}

double Armado::iteracion4(int i, int j, int kGauss, Double2D pos[4],
		double phi[2*NGAUSS][4], double dphi[NDIM][2*NGAUSS][4], double phidX[NDIM][2*NGAUSS][4]) {

	const int NODPEL = 4;

	double t = GAUSSPT[i];
	double s = GAUSSPT[j];

	double sm = 0.5 * (1.0 - s);
	double tm = 0.5 * (1.0 - t);
	double sq = 0.5 * (1.0 + s);
	double tp = 0.5 * (1.0 + t);

	int k = 0;
	phi[kGauss][k++] = sm * tm;
	phi[kGauss][k++] = sq * tm;
	phi[kGauss][k++] = sq * tp;
	phi[kGauss][k++] = sm * tp;

	k = 0;
	dphi[0][kGauss][k++] = -0.5 * tm;
	dphi[0][kGauss][k++] =  0.5 * tm;
	dphi[0][kGauss][k++] =  0.5 * tp;
	dphi[0][kGauss][k++] = -0.5 * tp;

	k = 0;
	dphi[1][kGauss][k++] = -0.5 * sm;
	dphi[1][kGauss][k++] = -0.5 * sq;
	dphi[1][kGauss][k++] =  0.5 * sq;
	dphi[1][kGauss][k++] =  0.5 * sm;

	double aJacob[2][2];
	for (int k = 0; k < 2; k++) for (int l = 0; l < 2; l++) aJacob[k][l] = 0.0;

	for (k = 0; k < NODPEL; k++) {
		aJacob[0][0] += dphi[0][kGauss][k] * pos[k].x;
		aJacob[0][1] += dphi[0][kGauss][k] * pos[k].y;
		aJacob[1][0] += dphi[1][kGauss][k] * pos[k].x;
		aJacob[1][1] += dphi[1][kGauss][k] * pos[k].y;
	}

	double det =
		aJacob[0][0] * aJacob[1][1] -
		aJacob[0][1] * aJacob[1][0];

	double aJacobInv[2][2] = {
		{  aJacob[1][1] / det, -aJacob[0][1] / det, },
		{ -aJacob[1][0] / det,  aJacob[0][0] / det,	},
	};

	for (int k = 0; k < NODPEL; k++) {
		phidX[0][kGauss][k] =
			aJacobInv[0][0] * dphi[0][kGauss][k] +
			aJacobInv[0][1] * dphi[1][kGauss][k];

		phidX[1][kGauss][k] =
			aJacobInv[1][0] * dphi[0][kGauss][k] +
			aJacobInv[1][1] * dphi[1][kGauss][k];
	}

	return det;
}

void Armado::armadoTransporte(int nodpel, Double2D pos[], double sigma, double qe, double landa,
		double mu, double mas[], double sol[], double esm[][MAXNPEL], double ef[], double deltaT) {

	const double th2 = 0.5;
	const double aCoef1 = deltaT * th2;
	const double aCoef2 = deltaT * (1 - th2);

	double est[MAXNPEL][MAXNPEL];
	for (int i = 0; i < nodpel; i++) for (int j = 0; j < nodpel; j++) est[i][j] = 0;

	switch (nodpel) {
	case 3:
		armadoTransporte3(pos, sigma, qe, landa, mu, mas, sol, esm, ef, est);
		break;
	case 4:
		armado4(pos, sigma, qe, true, landa, mu, esm, ef, est, mas);
		break;
	default: 
		assert(false);
	}

	for (int k = 0; k < nodpel; k++) {
		double sum = 0;
		for (int j = 0; j < nodpel; j++) {
			sum += (est[k][j] - aCoef2 * esm[k][j]) * sol[j];
			esm[k][j] = est[k][j] + aCoef1 * esm[k][j];
		}
		ef[k] = (aCoef1 + aCoef2) * ef[k] + sum;
	}
}

void Armado::armadoTransporte3(Double2D pos[], double sigma, double qe, double landa,
		double mu, double mas[], double sol[], double esm[][MAXNPEL], double ef[], double est[][MAXNPEL]) {

	const int NODPEL = 3;
	double b[NODPEL], c[NODPEL];
	double rMed = 0;

	double det = determinante3(pos, b, c);
	assert(abs(det) > TOLER_AREA);

	for (int i = 0; i < NODPEL; i++) rMed += pos[i].x;
	rMed /= NODPEL;

	for (int i = 0; i < NODPEL; i++) {
		ef[i] = 0;

		for (int j = 0; j < NODPEL; j++) {
			double a = (i==j) ? 2 : 1;
			ef[i] += qe * det * (M_PI / 12) * a * pos[j].x;
			esm[i][j] = (b[i] * b[j] + c[i] * c[j]) * sigma * M_PI  * rMed / det;
		}

		est[i][i] = mas[i] * landa * M_PI * rMed / det;
	}
}
