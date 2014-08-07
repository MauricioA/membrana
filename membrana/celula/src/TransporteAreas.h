#pragma once

#include "Transporte.h"
#include "Poros.h"

/* La difusión en la membrana depende de las areas de los poros */
class TransporteAreas :	public Transporte {
public:
	TransporteAreas(Celula& celula, Poros& poros);

protected:
	Poros* _poros;

	Poros& getPoros();

	double difusionMembrana(int iElem, int especie);

};
