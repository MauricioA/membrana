#ifndef ENTRADASALIDA_H_
#define ENTRADASALIDA_H_

#include <chrono>
#include "Poros.h"

class Celula;

class EntradaSalida {
public:
	EntradaSalida(Celula& celula);

	~EntradaSalida();

	void leerInput();

	void grabarPoisson(bool verbose=true);

	void printStart(string message, bool verbose=true);

	void printEnd(int tabs=2, bool verbose=true);

	void grabarTransporte(bool verbose=true);

	void grabarPoros(Poros& radios, bool verbose=true);

	void grabarITV(Poros& poros);

private:
	Celula* _celula;
	FILE* fitv;

	int nPoisson = 0;
	int nTransporte = 0;
	int nPoros = 0;

	inline Celula& getCelula();

	void leerMalla(string malla);

	void dameLinea(ifstream& archivo, istringstream& iss);

	chrono::time_point<chrono::high_resolution_clock> start;

};

#endif /* ENTRADASALIDA_H_ */
