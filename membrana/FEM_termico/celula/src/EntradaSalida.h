#ifndef ENTRADASALIDA_H_
#define ENTRADASALIDA_H_

#include <ctime>
#include "Celula.h"
#include "Poros.h"

class EntradaSalida {
public:
	static void leerInput(Celula& celula);

	static void grabarPoisson(Celula& celula, bool verbose=false);

	static void printStart(string message, bool verbose=true);

	static void printEnd(int tabs=2, bool verbose=true);

	static void grabarTransporte(Celula& cel, double time, bool verbose=true);

	static void grabarRadio(Celula& celula, Poros& radios, double time, bool verbose=true);

	static void grabarITV(Celula& celula, Poros& poros, double time);

private:
	static clock_t start;

	static bool firstWriteTransporte;

	static bool firstWritePoros;

	static bool firstWriteITV;

	static void leerMalla(Celula& celula, string malla);

	static void dameLinea(ifstream& archivo, istringstream& iss);

};

#endif /* ENTRADASALIDA_H_ */
