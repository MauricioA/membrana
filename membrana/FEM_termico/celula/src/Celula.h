#ifndef CELULA_H_
#define CELULA_H_

#include <Eigen/Sparse>
#include <vector>
#include "declaraciones.h"

using namespace std;
using namespace Eigen;
using namespace declaraciones;

//TODO poner en nodos y elementos todo
class Celula {
public:
	Celula();

	void poisson();

	void poros();

	void transportePoros();

	void chequearSimetria();

	int nodpel;
	int nElems;
	int nNodes;

	double potencial;
	double radio;
	double alto;
	double ancho;
	double area;
	double sigmas[3];

	string salida;

	VectorXd concentraciones[NESPS];
	VectorXd anteriores[NESPS];

	vector<double> phAux[NESPS];

	SparseMatrix<double> matrizTrans[NESPS];

	inline vector<Nodo>& getNodos() {
		return nodos;
	}

	inline vector<Elemento>& getElementos() {
		return elementos;
	}

	inline vector<Double2D>& getGradElem() {
		return gradElem;
	}

	inline vector<Double2D>& getCorrElem() {
		return corrElem;
	}

	inline SparseMatrix<double>& getMatriz() {
		return matriz;
	}

	inline VectorXd& getRhs() {
		return rhs;
	}

	inline VectorXd& getSolucion() {
		return solucion;
	}

	inline vector<double>& getMasas() {
		return masas;
	}

	/* Acá se hace una copia entera! */
	inline void setSolucion(VectorXd sol) {
		solucion = sol;
	}

	inline Double2D getCenter() {
		Double2D center;
		center.x = 0;
		center.y = alto / 2;
		return center;
	}

private:
	vector<Nodo>     nodos;
	vector<Elemento> elementos;
	vector<Double2D> gradElem;
	vector<Double2D> corrElem;

	SparseMatrix<double> matriz;

	VectorXd rhs;
	VectorXd solucion;

	vector<double> masas;

};

#endif
