#ifndef TRANSPORTE_H_
#define TRANSPORTE_H_

#include "Celula.h"

class Transporte {
public:
	Transporte(Celula& celula);

	void iteracion(double deltaT);

protected:
	Celula* _celula;

	Transporte() {}

	inline Celula& getCelula();

	void masaDiag2D();

	void concentracion(int esp, double deltaT);

	virtual double difusionMembrana(int iElem, int especie) = 0;

};

#endif /* TRANSPORTE_H_ */
