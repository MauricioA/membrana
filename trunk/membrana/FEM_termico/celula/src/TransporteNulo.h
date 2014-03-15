#pragma once

#include "Transporte.h"

class TransporteNulo : public Transporte {
public:
	TransporteNulo(Celula& celula) : Transporte(celula) {};

protected:
	inline double difusionMembrana(int iElem, int especie) {
		return DIFUSION[especie] * 1e-3;
	}

};
