#include <cassert>
#include "Celula.h"
#include "EntradaSalida.h"
#include "Armado.h"
#include "Poisson.h"
#include "Transporte.h"

Celula::Celula() {
	potencial = 0;
	nNodes = nElems = nodpel = 0;

	EntradaSalida::leerInput(*this);
}

void Celula::poisson() {
	Poisson::poisson(*this);
	EntradaSalida::grabarPoisson(*this);
}

void Celula::transporte() {
	Transporte::transporte(*this);
}

void Celula::chequearSimetria() {
	for (int i = 0; i < nNodes; i++) for (int j = i; j < nNodes; j++) {
		assert(abs(matriz.coeff(i, j) - matriz.coeff(j, i)) < 1e-9);
	}
}

//inline vector<Nodo>& Celula::getNodos() {
//	return nodos;
//}
//
//vector<Elemento>& Celula::getElementos() {
//	return elementos;
//}
//
//vector<Double2D>& Celula::getGradElem() {
//	return gradElem;
//}
//
//vector<Double2D>& Celula::getCorrElem() {
//	return corrElem;
//}
//
//SparseMatrix<double>& Celula::getMatriz() {
//	return matriz;
//}
//
//VectorXd& Celula::getRhs() {
//	return rhs;
//}
//
//VectorXd& Celula::getSolucion() {
//	return solucion;
//}
//
//void Celula::setSolucion(VectorXd sol) {
//	solucion = sol;
//}
//
//vector<double>& Celula::getMasas() {
//	return masas;
//}
//
//vector<double>& Celula::getCargas() {
//	return cargas;
//}
