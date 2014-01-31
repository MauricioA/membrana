typedef unsigned int uint;

const int 	 N_COTA 	= 10;
const int 	 MAXNPEL 	= 4;
const int 	 NESPS		= 4;
const double TOLER_AREA = 1e-6;
const double TOLER_MASA = 1e-12;
const double EPSILON 	= 1e-3;
const double TIERRA		= 0.0;

enum Especie {
	H_,	OH,	NA,	CL,
};

const double CONCENTRACION_INICIAL[] = {	// En at/(um**3)
	60.20,		// H_ 1e-7 M
	60.20,		// OH 1e-7 M
	96.32e6,	// NA 0.16 M
	96.32e6,	// CL 0.16 M
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
