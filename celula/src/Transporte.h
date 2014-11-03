#ifndef TRANSPORTE_H_
#define TRANSPORTE_H_

#include "Celula.h"
#include "Poros.h"

///* VALORES FORTRAN! */
const double CONC_INICIAL_INTRA[] = {	// at um**-3
	60.20,
	60.20,
	96320000.0 / 2,	//la mitad para intra!
	96320000.0 / 2,	//la mitad para intra!
};

/* VALORES FORTRAN! */
const double CONC_INICIAL_EXTRA[] = {	// at um**-3
	60.20,
	60.20,
	96320000.0,
	96320000.0,
};

/* VALORES FORTRAN! */
const double CONC_ANODO[] = {		// at um**-3
	60.20,
	60.20,
	96320000.0,
	96320000.0,
};

/* VALORES FORTRAN! */
const double CONC_CATODO[] = {		// at um**-3
	60.20,
	60.20,
	96320000.0,
	96320000.0,
};

///* Concentraciones diferentes en citoplasma */	//estaba esta
//const double CONC_INICIAL_INTRA[] = {	// at um**-3
//	CONCENT * 0.3978e-7,	// H  .3978e-7 M
//	CONCENT * 0.3978e-7,	// OH .3978e-7 M
//	CONCENT * 142e-3,		// NA 142 mM
//	CONCENT * 108e-3,		// CL 108 mM
//};
//
///* Concentraciones diferentes en citoplasma */	//estaba esta
//const double CONC_INICIAL_EXTRA[] = {	// at um**-3
//	CONCENT *  1e-7,		// H  1e-7 M
//	CONCENT *  1e-7,		// OH 1e-7 M
//	CONCENT * 14e-3,		// NA 14 mM
//	CONCENT *  4e-3,		// CL  4 mM
//};

/* Concentraciones viejas */
//const double CONCENTRACION_ANODO[] = {			// at um**-3
//	1.5e7,		// H
//	0,			// OH
//	1e12,		// NA
//	0,			// CL
//};
//
//const double CONCENTRACION_CATODO[] = {			// at um**-3
//	0,			// H
//	1.806e7,	// OH
//	0,			// NA
//	0,			// CL
//};

///* Mismos valores en los electrodos que en extra */
//const double CONC_ANODO[] = {		// at um**-3
//	CONCENT *  1e-7,		// H  1e-7 M
//	CONCENT *  1e-7,		// OH 1e-7 M
//	CONCENT * 14e-3,		// NA 14 mM
//	CONCENT *  4e-3,		// CL  4 mM
//};
//
//const double CONC_CATODO[] = {		// at um**-3
//	CONCENT *  1e-7,		// H  1e-7 M
//	CONCENT *  1e-7,		// OH 1e-7 M
//	CONCENT * 14e-3,		// NA 14 mM
//	CONCENT *  4e-3,		// CL  4 mM
//};

const double DIFUSION[] = {	// um**2 / s
	12500,		// H
	7050,		// OH
	1780,		// NA
	2720,		// CL
};

const double CARGA[] = {
	+1,			// H
	-1,			// OH
	+1,			// NA
	-1,			// CL
};

class Transporte {
public:
	Transporte(Celula& celula, bool calcularPoros);

	void iteracion(double deltaT);

	bool _calcularPoros;
	
	Poros* _poros;

	double error = 0;

protected:
	Celula* _celula;

	vector<double> masas;

	SparseMatrix<double> matrizTrans[NESPS];

	Transporte() {}

	inline Celula& getCelula();

	void masaDiag2D();

	void concentracion(int esp, double deltaT);

	inline double difusionMembrana(int iElem, int especie);

};

/* valores de paper chinos
	const double CONCENTRACION_INICIAL_INTRA[] = {	// at um**-3
	CONCENT * 220e-9,
	CONCENT * 220e-9,
	CONCENT * 220e-9,
	CONCENT * 220e-9,
	};

	const double CONCENTRACION_INICIAL_EXTRA[] = {	// at um**-3
	CONCENT *  1e-3,
	CONCENT *  1e-3,
	CONCENT *  1e-3,
	CONCENT *  1e-3,
	};

	const double CONCENTRACION_ANODO[] = {		// at um**-3
	CONCENT *  1e-3,
	CONCENT *  1e-3,
	CONCENT *  1e-3,
	CONCENT *  1e-3,
	};

	const double CONCENTRACION_CATODO[] = {		// at um**-3
	CONCENT *  1e-3,
	CONCENT *  1e-3,
	CONCENT *  1e-3,
	CONCENT *  1e-3,
	};

	const double CARGA[] = {	//cualquier cosa!!
	+2,			// H
	-5,			// OH
	+2,			// NA
	-5,			// CL
	};
*/


#endif /* TRANSPORTE_H_ */
