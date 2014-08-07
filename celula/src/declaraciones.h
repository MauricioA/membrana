#ifndef DECLARACIONES_H_
#define DECLARACIONES_H_

#define _USE_MATH_DEFINES

#include <math.h>
#include <string>

using namespace std;

typedef unsigned int uint;

namespace declaraciones {

	//todo en um!
	const int 	 N_COTA 			= 10;
	const int 	 MAXNPEL 			= 4;
	const int 	 NESPS				= 4;
	const int 	 NGAUSS 			= 2;
	const int 	 NDIM 				= 2;
	const double EPSILON_DIST		= 1e-9;
	const double TOLER_AREA 		= 1e-6;
	const double TOLER_MASA 		= 1e-12;
	const double TOLER_DIST 		= 1e-3;
	const double TIERRA				= 0;
	const double FARADAY 			= 96485.34;		// C/mol cte. de Faraday
	const double R_CTE				= 8.314; 		// J/K/mol cte. de gases
	const double T_CTE				= 310;			// K temperatura
	const double EPSILON_TRANSPORTE	= 78.5;			// cte dieléctrica del agua
	const double EPSILON_0			= 8.85e-12;		// cte de permisividad C**2 / (N m**2)
	const double CONCENT_H2O 		= 3.34e10;		// at um**-3
	const double RSA 				= 0.5;
	const double CONCENT_MINIMO 	= 1e-8;
	const double CONCENT			= 6.02214129e23 * 1e-15;

	const double GAUSSPT[] = {
		-1 / sqrt(3.0),
		 1 / sqrt(3.0),
	};

	const double GAUSSWT[] = { 1, 1 };

	enum Especie {
		H_,	OH,	NA,	CL,
	};

	enum Estado { ON, OFF };

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

	enum Material {
		EXTERNO,
		MEMBRANA,
		INTERNO,
	};

	struct Nodo {
		double x;
		double y;
		bool esTierra;
		bool esPotencia;
	};

	class Elemento {
	private:
		int nodos[MAXNPEL];

	public:
		Material material;
		double sigma;

		Elemento(int nodos_[], int nodpel, Material material_, double _sigma) {
			for (int i = 0; i < nodpel; i++) nodos[i] = nodos_[i];
			material = material_;
			sigma = _sigma;
		}

		inline int operator[](int i) {
			return nodos[i];
		}
	};

	struct Double2D {
		double x;
		double y;
	};
}

#endif /* DECLARACIONES_H_ */
