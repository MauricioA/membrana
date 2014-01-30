typedef unsigned int uint;

const int 	 N_COTA 	= 10;
const double TOLER_AREA = 1e-6;
const double EPSILON 	= 1e-3;
const double TIERRA		= 0.0;

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
	int nodos[4];		//Declaro de sobra porque sí

public:
	Elemento(int nodos_[], int nodpel, Material materialElem) {
		for (int i = 0; i < nodpel; i++) nodos[i] = nodos_[i];
		material = materialElem;
	}

	Material material;

	int operator[](int i) {
		return nodos[i];
	}
};
