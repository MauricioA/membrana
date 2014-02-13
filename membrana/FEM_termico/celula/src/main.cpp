#include "Celula.h"
#include <iostream>

//TODO cambiar clases estáticas a singleton

#include "Poros.h"

int main() {

	Poros& poros = Poros::instance();
	Poros& poros3 = Poros::instance();

	cout << poros.s << endl;

	poros.s = "222";

	cout << poros.s << endl;

	Poros& poros2 = Poros::instance();

	cout << poros2.s << endl;
	cout << poros3.s << endl;

//	BREAKPOINT
//
//	Celula problema;
//
//	problema.transporte();
//
//	return EXIT_SUCCESS;
}
