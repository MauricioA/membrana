#ifndef POISSON_H_
#define POISSON_H_

#include "Celula.h"

class Poisson {
public:
	static void poisson(Celula& celula);

private:
	static VectorXd global_rhs;

	static SparseMatrix<double> matriz;

	static void campo(Celula& celula);

	static void campo3(Celula& celula);

	static void campo4(Celula& celula);

	static void corriente(Celula& celula);
};

#endif /* POISSON_H_ */
