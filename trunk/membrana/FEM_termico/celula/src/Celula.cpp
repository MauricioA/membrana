#include <cassert>
#include <iostream>
#include "Celula.h"
#include "EntradaSalida.h"
#include "Armado.h"
#include "Poisson.h"
#include "TransporteAreas.h"
#include "Poros.h"

#include "TransporteNulo.h"

//TODO poner en nodos y elementos todo

Celula::Celula() {
	potencial = 0;
	nNodes = nElems = nodpel  = 0;
	alto = radio = ancho = 0;

	EntradaSalida::leerInput(*this);

	area = 4 * M_PI * pow(radio, 2);
}

void Celula::poisson() {
	Poisson::poisson(*this);
	EntradaSalida::grabarPoisson(*this);
}

void Celula::poros() {
	Poisson::poisson(*this);

	Poros poros(*this);
	poros.loop();
}

void Celula::chequearSimetria() {
	for (int i = 0; i < nNodes; i++) for (int j = i; j < nNodes; j++) {
		assert(abs(matriz.coeff(i, j) - matriz.coeff(j, i)) < 1e-9);
	}
}

void Celula::transporteYPoros() {
	const double TIEMPO_FINAL = 1;
	const double DELTA_T = 1e-6;
	
	Poros poros = Poros(*this);
	TransporteNulo transporte = TransporteNulo(*this);
	
	Poisson::poisson(*this);

	int iter = 0;
	clock_t reloj = 0;

	for (double time = 0; time < TIEMPO_FINAL; time += DELTA_T) {
		transporte.iteracion(DELTA_T);
		//poros.iteracion();

		if (iter % PASO_CONSOLA == 0 && iter != 0) {
			int interv = (clock() - reloj) / (CLOCKS_PER_SEC / 1000);
			cout << time*1e6 << "us\t"
				<< iter << " iters\t"
				<< interv / PASO_CONSOLA << "ms/it" << endl;

			if (iter % PASO_DISCO == 0) {
				EntradaSalida::grabarTransporte(*this, time);
			}

			reloj = clock();
		}
		iter++;
	}
}
