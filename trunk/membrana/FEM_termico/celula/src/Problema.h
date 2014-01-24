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
	vector<double>	 campoElemX;
	vector<double>   campoElemY;
	vector<double>	 corrElemX;
	vector<double>   corrElemY;

	SparseMatrix<double> matriz;

	VectorXd rhs;
	VectorXd solucion;

	void leerMalla(string malla);

	void armado3(double x[], double y[], double esm[3][3], double sigma);

	void armado4(double x[], double y[], double esm[4][4], double sigma);

	void campo();

	void corriente();

	void grabar();

	void chequearSimetria();

	double determinante(double x[], double y[], double b[], double c[]);
};

#endif
