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

	void transporte();

	Elemento getElement(int i);

private:
	int nElems;
	int nNodes;
	int nodpel;

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

	void dameLinea(ifstream& archivo, istringstream& iss);

	void armado (double x[], double y[], double esm[][MAXNPEL], double sigma);

	void armado3(double x[], double y[], double esm[][MAXNPEL], double sigma);

	void armado4(double x[], double y[], double esm[][MAXNPEL], double sigma);

	double determinante3(double x[], double y[], double b[], double c[]);

	void campo();

	void corriente();

	void grabar();

	void chequearSimetria();

	void masaDiag2D();

};

#endif
