#ifndef PROBLEMA_H_
#define PROBLEMA_H_

#include <Eigen/Sparse>
#include <vector>
#include "declaraciones.cpp"

using namespace std;
using namespace Eigen;

class Problema {
public:
	Problema();

	void poisson();

	Elemento getElement(int i);

private:
	int nElems;
	int nNodes;

	double sigmas[3];
	double potencial;

	vector<Nodo>     nodos;
	vector<Elemento> elementos;
	vector<double>	 corrElemX;
	vector<double>   corrElemY;

	SparseMatrix<double> matriz;

	VectorXd rhs;
	VectorXd solucion;

	void leerMalla(string malla);

	void control();

	void armado(double x[], double y[], double ef[], double qe, double esm[3][3], double sigma);

	void corriente2D();

	void grabar();

	void chequearSimetria();
};

#endif
