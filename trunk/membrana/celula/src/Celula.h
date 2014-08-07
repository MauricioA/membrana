#ifndef CELULA_H_
#define CELULA_H_

#include <Eigen/Sparse>
#include <vector>
#include <memory>
#include "declaraciones.h"

using namespace std;
using namespace Eigen;
using namespace declaraciones;

class EntradaSalida;

//TODO poner en nodos y elementos todo
class Celula {
public:
	int nodpel;
	int nElems;
	int nNodes;
	int threads;
	int pulsos;

	string salida;

	double potencial;
	double radio;
	double alto;
	double ancho;
	double area;
	double time;
	double on_time;
	double times[2];
	double sigmas[3];

	vector<Nodo>     nodos;
	vector<Elemento> elementos;
	vector<Double2D> gradElem;
	vector<Double2D> corrElem;

	VectorXd concentraciones[NESPS];
	VectorXd conc_anteriores[NESPS];
	VectorXd potenciales;

	Celula();

	~Celula();

	void poisson();

	void transportePoros();

	void chequearSimetria();

	inline EntradaSalida& getEntradaSalida() {
		return *_entradaSalida;
	}

	inline Double2D getCenter() {
		return Double2D{ 0, alto / 2 };
	}

private:
	unique_ptr<EntradaSalida> _entradaSalida;

};

#endif
