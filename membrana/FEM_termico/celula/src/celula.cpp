
#include <Eigen/Sparse>

#include <stdlib.h>

//TODO borrar los que sobren
#include <iostream>
#include <vector>
#include <fstream>
#include "Problema.h"

using namespace std;

//TODO borrar
using namespace Eigen;

int main() {
/*
	SparseMatrix<double> mat(10, 10);
	vector< Triplet<double> > vec;

	vec.push_back(Triplet<double>(1, 5, 10));
	vec.push_back(Triplet<double>(1, 3, 10));
	vec.push_back(Triplet<double>(2, 5, 10));
	vec.push_back(Triplet<double>(4, 1, 10));
	vec.push_back(Triplet<double>(1, 3, 3));

	mat.setFromTriplets(vec.begin(), vec.end());

	cout << mat.coeff(1, 3) << endl;

	vec.clear();

	mat.setFromTriplets(vec.begin(), vec.end()),

	cout << mat.coeff(1, 3) << endl;
*/



	Problema problema;

	problema.poisson();

	return EXIT_SUCCESS;
}



