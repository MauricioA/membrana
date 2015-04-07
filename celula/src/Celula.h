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

class Celula {
public:
	int nodpel;
	int nElems;
	int nNodes;
	int threads;
	int pulsos;

	string salida;

	bool soloPoisson;
	bool calcularPoros;
	bool calcularTransporte;

	double potencial;
	double radio;
	double alto;
	double ancho;
	double area;
	double time;
	double on_time;
	double delta_t;
	double times[2];
	double sigmas[3];

	vector<Nodo>     nodos;
	vector<Elemento> elementos;
	vector<Double2D> gradElem;
	vector<Double2D> corrElem;

	VectorXd concs[NESPS];	//concentraciones actuales
	VectorXd c_ant[NESPS];	//concentraciones iteración anterior
	VectorXd potenciales;

	Celula();

	~Celula();

	void loop();

	void acoplado();

	void poisson();

	void chequearSimetria();

	int valoresExtremos();

	inline EntradaSalida& getEntradaSalida() {
		return *_entradaSalida;
	}

	inline Double2D getCenter() {
		return Double2D { 0, alto / 2 };
	}

private:
	unique_ptr<EntradaSalida> _entradaSalida;

};

#endif
