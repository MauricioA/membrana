#include "Celula.h"
#include <iostream>

int main() {
	#ifdef EIGEN_VECTORIZE
	cout << "Vectorization ON\n";
	#else
	cout << "Vectorization OFF\n";
	#endif

	Celula celula;

	//celula.transportePoros();
	celula.poisson();

	return EXIT_SUCCESS;
}
