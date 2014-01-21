typedef unsigned int uint;

const int 	 NODPEL 	= 3;
const int 	 N_COTA 	= 10;
const double TOLER_AREA = 1e-6;
const double EPSILON 	= 1e-3;
const double TIERRA		= 0.;

enum Material {
	EXTERNO,
	MEMBRANA,
	INTERNO,
};

struct Nodo {
	double x;
	double y;
	double potencial;	//TODO esto se usa?
	bool esTierra;
	bool esPotencia;
};

class Elemento {
private:
	int nodos[NODPEL];

public:
	Elemento(int nodos_[NODPEL], Material materialElem) {
		for (int i = 0; i < NODPEL; i++) nodos[i] = nodos_[i];
		material = materialElem;
	}

	Material material;

	int operator[](int i) {
		return nodos[i];
	}
};
