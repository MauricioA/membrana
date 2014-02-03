#ifndef PROBLEMA_H_
#define PROBLEMA_H_

#include <Eigen/Sparse>
#include <vector>
#include "declaraciones.cpp"

using namespace std;
using namespace Eigen;

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

	vector<Nodo> &getNodos() {
		return nodos;
	}

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

	void leerMalla(string malla);

	void dameLinea(ifstream& archivo, istringstream& iss);

	void armado (Double2D pos[], double esm[][MAXNPEL], double sigma);

	void armado3(Double2D pos[], double esm[][MAXNPEL], double sigma);

	void armado4(Double2D pos[], double esm[][MAXNPEL], double sigma, double ef[MAXNPEL],
			double qe, bool transp, double landa, double mu, double est[][4], double mas[]);

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

	void armadoTransporte(Double2D pos[], double esm[][MAXNPEL], double sigma, double qe, double landa,
			double mu, double sol[], double ef[]);
};

#endif
