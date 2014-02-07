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

	void chequearSimetria();

	int nodpel;
	int nElems;
	int nNodes;

	double potencial;
	double sigmas[3];

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

//	inline vector<double>& getCargas() {
//		return cargas;
//	}

	/* Acá se hace una copia entera! */
	inline void setSolucion(VectorXd sol) {
		solucion = sol;
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
//	vector<double> cargas;

};

#endif
