#ifndef ENTRADASALIDA_H_
#define ENTRADASALIDA_H_

#include "Celula.h"

class EntradaSalida {
public:
	static void leerInput(Celula& celula);

	static void grabar(Celula& celula);

	static void printStart(string message, bool verbose=true);

	static void printEnd(int tabs=2, bool verbose=true);

private:
	static clock_t start;

	static void leerMalla(Celula& celula, string malla);

	static void dameLinea(ifstream& archivo, istringstream& iss);

};

#endif /* ENTRADASALIDA_H_ */
