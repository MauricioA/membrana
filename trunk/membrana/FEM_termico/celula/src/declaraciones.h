#ifndef DECLARACIONES_H_
#define DECLARACIONES_H_

#define _USE_MATH_DEFINES

#include <math.h>
#include <string>

using namespace std;

#define BREAKPOINT cout << "";

#define uint unsigned int

namespace declaraciones {
	const int 	 N_COTA 			= 10;
	const int 	 MAXNPEL 			= 4;
	const int 	 NESPS				= 4;
	const int 	 NGAUSS 			= 2;
	const int 	 NDIM 				= 2;
	const int	 PASO_DISCO			= 1000;
	const int 	 PASO_CONSOLA 		= 100;
	const double TOLER_AREA 		= 1e-6;
	const double TOLER_MASA 		= 1e-12;
	const double TOLER_DIST 		= 1e-3;
	const double EPSILON_POISSON	= 1e-3;
	const double TIERRA				= 0;
	const double FARADAY 			= 96485.34;		// C/mol
	const double R_CTE				= 8.314; 		// J/K/mol
	const double T_CTE				= 310;			// K
	const double EPSILON_TRANSPORTE	= 78.5;			// cte dieléctrica del agua
	const double EPSILON_0			= 8.85e-12;		// cte de permitividad C**2 / (N m**2)
	const double DELTA_T_TRANSPORTE = 10e-6;
	const double KWF				= 249.16;		// (um**3)/at/s = 1.5e11 1/mol/s
	const double KWB 				= 2.7e-5;		// 1/s
	const double CONCENT_H2O 		= 3.34e10;		// at/(um**3)
	const double RSA 				= 0.5;
	const double CONCENT_MINIMO 	= 1e-8;
	const double CLAVE 				= FARADAY / (R_CTE * T_CTE);
	const string RUTA_HISTORIAL		= "salida/historia.dat";
	const string RUTA_PH			= "salida/ph.dat";

	const double GAUSSPT[] = {
		-1 / sqrt(3.),
		 1 / sqrt(3.),
	};

	const double GAUSSWT[] = { 1, 1 };

	enum Especie {
		H_,	OH,	NA,	CL,
	};

	const double CONCENTRACION_INICIAL[] = {	// En at/(um**3)
		60.2,		// H_ 1e-7 M
		60.2,		// OH 1e-7 M
		96.32e6,	// NA 0.16 M
		96.32e6,	// CL 0.16 M
	};

	const double CONCENTRACION_ANODO[] = {		// En at/(um**3)
		1.5e7,		// H
		0,			// OH
		1e12,		// NA
		0,			// CL
	};

	const double CONCENTRACION_CATODO[] = {		// En at/(um**3)
		0,			// H
		1.806e7,	// OH
		0,			// NA
		0,			// CL
	};

	const double CARGA[] = {
		+1,			// H
		-1,			// OH
		+1,			// NA
		-1,			// CL
	};

	const double DIFUSION[] = {		// (creo)
		12500,		// H
		 7050,		// OH
		 1780,		// NA
		 3830,		// CL
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

		Elemento(int nodos_[], int nodpel, Material material_) {
			for (int i = 0; i < nodpel; i++) nodos[i] = nodos_[i];
			material = material_;
		}

		int operator[](int i) {
			return nodos[i];
		}
	};

	struct Double2D {
		double x;
		double y;
	};
}

#endif /* DECLARACIONES_H_ */
