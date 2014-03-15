#pragma once

#include "Transporte.h"
#include "Poros.h"

/* La difusi�n en la membrana depende de las areas de los poros */
class TransporteAreas :	public Transporte {
public:
	TransporteAreas(Celula& celula, Poros& poros);

protected:
	Poros* _poros;

	Poros& getPoros();

	double difusionMembrana(int iElem, int especie);

};
