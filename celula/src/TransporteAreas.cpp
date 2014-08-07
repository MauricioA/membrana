#include "TransporteAreas.h"

TransporteAreas::TransporteAreas(Celula& celula, Poros& poros) 
: Transporte(celula) {
	_poros = &poros;
}

Poros& TransporteAreas::getPoros() {
	return *_poros;
}

double TransporteAreas::difusionMembrana(int iElem, int especie) {
	return DIFUSION[especie] * getPoros().getProporsionArea(iElem);
}
