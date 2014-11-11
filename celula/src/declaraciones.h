#ifndef DECLARACIONES_H_
#define DECLARACIONES_H_

#define _USE_MATH_DEFINES

#include <math.h>
#include <string>
#include <cinttypes>

using namespace std;

typedef unsigned int uint;

namespace declaraciones {
	const int	 MAXNPEL 	= 4;
	const int 	 NESPS		= 4;
	const int 	 NGAUSS 	= 2;
	const int 	 NDIM 		= 2;
	const double TOLER_AREA = 1e-6;
	const double CONCENT	= 6.02214129e23 * 1e-15;

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
		uint16_t nodos[MAXNPEL];

	public:
		Material material;
		double sigma;

		Elemento(int nodos_[], int nodpel, Material material_, double _sigma) {
			for (int i = 0; i < nodpel; i++) nodos[i] = nodos_[i];
			material = material_;
			sigma = _sigma;
		}

		inline uint16_t operator[](int i) {
			return nodos[i];
		}
	};

	struct Double2D {
		double x;
		double y;
	};
}

#endif /* DECLARACIONES_H_ */
