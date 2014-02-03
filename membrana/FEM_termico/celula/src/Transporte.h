#ifndef TRANSPORTE_H_
#define TRANSPORTE_H_

#include "Celula.h"

class Transporte {
public:
	static void transporte(Celula& cel);

private:
	static void masaDiag2D(Celula& cel);

	static void carga(Celula& cel);

	static void concentracion(Celula& cel, int esp);

};

#endif /* TRANSPORTE_H_ */
