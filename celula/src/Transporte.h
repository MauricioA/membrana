#ifndef TRANSPORTE_H_
#define TRANSPORTE_H_

#include "Celula.h"

const double CONCENTRACION_INICIAL_INTRA[] = {	// at um**-3
	CONCENT * 0.3978e-7,	// H  .3978e-7 M
	CONCENT * 0.3978e-7,	// OH .3978e-7 M
	CONCENT * 142e-3,		// NA 142 mM
	CONCENT * 108e-3,		// CL 108 mM
};

const double CONCENTRACION_INICIAL_EXTRA[] = {	// at um**-3
	CONCENT *  1e-7,		// H  1e-7 M
	CONCENT *  1e-7,		// OH 1e-7 M
	CONCENT * 14e-3,		// NA 14 mM
	CONCENT *  4e-3,		// CL  4 mM
};

const double CONCENTRACION_ANODO[] = {		// at um**-3
	CONCENT *  1e-7,		// H  1e-7 M
	CONCENT *  1e-7,		// OH 1e-7 M
	CONCENT * 14e-3,		// NA 14 mM
	CONCENT *  4e-3,		// CL  4 mM
};

const double CONCENTRACION_CATODO[] = {		// at um**-3
	CONCENT *  1e-7,		// H  1e-7 M
	CONCENT *  1e-7,		// OH 1e-7 M
	CONCENT * 14e-3,		// NA 14 mM
	CONCENT *  4e-3,		// CL  4 mM
};

const double DIFUSION[] = {	// um**2 / s
	12500,		// H
	7050,		// OH
	1780,		// NA
	3830,		// CL
};

const double CARGA[] = {
	+1,			// H
	-1,			// OH
	+1,			// NA
	-1,			// CL
};

//valores de paper chinos
//const double CONCENTRACION_INICIAL_INTRA[] = {	// at um**-3
//	CONCENT * 220e-9,
//	CONCENT * 220e-9,
//	CONCENT * 220e-9,
//	CONCENT * 220e-9,
//};
//
//const double CONCENTRACION_INICIAL_EXTRA[] = {	// at um**-3
//	CONCENT *  1e-3,
//	CONCENT *  1e-3,
//	CONCENT *  1e-3,
//	CONCENT *  1e-3,
//};
//
//const double CONCENTRACION_ANODO[] = {		// at um**-3
//	CONCENT *  1e-3,
//	CONCENT *  1e-3,
//	CONCENT *  1e-3,
//	CONCENT *  1e-3,
//};
//
//const double CONCENTRACION_CATODO[] = {		// at um**-3
//	CONCENT *  1e-3,
//	CONCENT *  1e-3,
//	CONCENT *  1e-3,
//	CONCENT *  1e-3,
//};
//
//const double CARGA[] = {	//cualquier cosa!!
//	+2,			// H
//	-5,			// OH
//	+2,			// NA
//	-5,			// CL
//};

class Transporte {
public:
	Transporte(Celula& celula);

	void iteracion(double deltaT);

protected:
	Celula* _celula;

	vector<double> masas;

	SparseMatrix<double> matrizTrans[NESPS];

	Transporte() {}

	inline Celula& getCelula();

	void masaDiag2D();

	void concentracion(int esp, double deltaT);

	virtual double difusionMembrana(int iElem, int especie) = 0;

};

#endif /* TRANSPORTE_H_ */
