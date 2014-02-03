#ifndef ENTRADASALIDA_H_
#define ENTRADASALIDA_H_

#include "Celula.h"

class EntradaSalida {
public:
	static void leerInput(Celula& celula);

	static void grabar(Celula& celula);

private:
	static void leerMalla(Celula& celula, string malla);

	static void dameLinea(ifstream& archivo, istringstream& iss);

};

#endif /* ENTRADASALIDA_H_ */
