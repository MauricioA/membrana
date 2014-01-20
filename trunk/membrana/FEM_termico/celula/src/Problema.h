#ifndef PROBLEMA_H_
#define PROBLEMA_H_

#include <Eigen/Sparse>
#include <vector>
#include "declaraciones.cpp"

using namespace std;
//using Eigen::Triplet;

class Problema {
public:
	Problema();

	void leerMalla(const char* archivo);

	void poisson();

	Elemento getElement(int i);

private:
	int nElems;
	int nNodes;

	double sigmas[3];
	double potencial;

	std::vector<Nodo>     nodos;
	std::vector<Elemento> elementos;
	std::vector<double> 	 solucionAnterior;
	std::vector<double>   gradElemX;
	std::vector<double>   gradElemY;


	Eigen::VectorXd rhs;
	Eigen::SparseMatrix<double> matriz;

	std::vector< Eigen::Triplet<double> > triplets;


	void control();

	void armado(double x[], double y[], double ef[], double qe, double esm[3][3], double sigma);

};

#endif
