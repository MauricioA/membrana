#ifndef PROBLEMA_H_
#define PROBLEMA_H_

#include <Eigen/Sparse>
#include <vector>
#include "declaraciones.h"

using namespace std;
using namespace Eigen;
using namespace declaraciones;

class Celula {
public:
	Celula();

	void poisson();

	void transporte();

	Elemento getElement(int i);

	int nodpel;
	int nElems;
	int nNodes;

	double potencial;
	double sigmas[3];

	vector<Nodo>& 	  getNodos();
	vector<Elemento>& getElementos();
	vector<Double2D>& getGradElem();
	vector<Double2D>& getCorrElem();

	SparseMatrix<double>& getMatriz();

	VectorXd& getRhs();
	VectorXd& getSolucion();

private:
	vector<Nodo>     nodos;
	vector<Elemento> elementos;
	vector<Double2D> gradElem;
	vector<Double2D> corrElem;

	SparseMatrix<double> matriz;

	VectorXd rhs;
	VectorXd solucion;

	vector<double>	 masas;
	vector<double> 	 cargas;
	vector<double> 	 phAux[NESPS];

	VectorXd concentraciones[NESPS];
	VectorXd anteriores[NESPS];

	void campo();

	void campo3();

	void campo4();

	void corriente();

	void grabar();

	void chequearSimetria();

	void masaDiag2D();

	void carga();

	void concentracion(int especie);

};

#endif
