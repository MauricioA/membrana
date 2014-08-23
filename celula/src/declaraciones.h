#ifndef DECLARACIONES_H_
#define DECLARACIONES_H_

#define _USE_MATH_DEFINES

#include <math.h>
#include <string>

using namespace std;

typedef unsigned int uint;

namespace declaraciones {

	//todo en um! (podrían ir en poisson o transporte)
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

	enum Especie {
		H_,	OH,	NA,	CL,
	};

	enum Material {
		EXTERNO,
		MEMBRANA,
		INTERNO,
	};

	enum Estado { ON, OFF };

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
