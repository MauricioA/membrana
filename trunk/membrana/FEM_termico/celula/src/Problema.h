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
	vector<Double2D> gradElem;
	vector<Double2D> corrElem;

	vector<double>	 masas;
	vector<double> 	 cargas;
	vector<double>   cons[NESPS];
	vector<double>   ants[NESPS];
	vector<double> 	phAux[NESPS];

	SparseMatrix<double> matriz;

	VectorXd rhs;
	VectorXd solucion;

	void leerMalla(string malla);

	void dameLinea(ifstream& archivo, istringstream& iss);

	void armado (Double2D pos[], double esm[][MAXNPEL], double sigma);

	void armado3(Double2D pos[], double esm[][MAXNPEL], double sigma);

	void armado4(Double2D pos[], double esm[][MAXNPEL], double sigma,
			bool transp, double landa, double mu, double est[][4], double mas[]);

	double determinante3(Double2D pos[], double b[], double c[]);

	double iteracion4(double phi[2*NGAUSS][4], double dphi[NDIM][2*NGAUSS][4],
			double phidX[NDIM][2*NGAUSS][4], int i, int j, int kGauss, Double2D pos[4]);

	void campo();

	void campo3();

	void campo4();

	void corriente();

	void grabar();

	void chequearSimetria();

	void masaDiag2D();

	void carga();

	void concentracion(int especie);

	void armadoTransporte();
};

#endif
