#ifndef ENTRADASALIDA_H_
#define ENTRADASALIDA_H_

#include <ctime>
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

	void grabarRadios(Poros& radios, bool verbose=false);

	void grabarITV(Poros& poros);

private:
	Celula* _celula;

	//FILE* ftension;
	//FILE* fcampo;
	FILE* fitv;
	FILE* fporos;
	//FILE* ftransporte;
	//FILE* fph;

	clock_t start;

	int nPoisson = 0;
	int nTransporte = 0;

	inline Celula& getCelula();

	void leerMalla(string malla);

	void dameLinea(ifstream& archivo, istringstream& iss);

};

#endif /* ENTRADASALIDA_H_ */
