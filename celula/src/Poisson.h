#ifndef POISSON_H_
#define POISSON_H_

#include "Celula.h"

//TODO no debería ser static!

class Poisson {
public:
	static int estado;

	static void iteracion(Celula& celula);

private:
	static VectorXd global_rhs;

	static SparseMatrix<double> matriz;

	static void campo(Celula& celula);

	static void campo3(Celula& celula);

	static void campo4(Celula& celula);

	static void corriente(Celula& celula);

};

#endif /* POISSON_H_ */
