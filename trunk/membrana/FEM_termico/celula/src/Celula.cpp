#include <cassert>
#include <iostream>
#include "Celula.h"
#include "EntradaSalida.h"
#include "Armado.h"
#include "Poisson.h"
#include "TransporteAreas.h"
#include "Poros.h"

#include "TransporteNulo.h"

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

void Celula::transportePoros() {
	const double TIEMPO_FINAL	 = 1;
	const double DELTA_T_POROS	 = 100e-9;
	const int	 PASO_TRANSPORTE = 1000;
	const int	 PASO_DISCO		 = 100000;
	const int 	 PASO_CONSOLA	 = 10000;
	
	Poros poros = Poros(*this);
	TransporteAreas transporte = TransporteAreas(*this, poros);
	
	Poisson::poisson(*this);

	int iter = 0;
	clock_t reloj = 0;

	for (double time = 0; time < TIEMPO_FINAL; time += DELTA_T_POROS) {
		poros.iteracion(DELTA_T_POROS, time);

		if (iter % PASO_TRANSPORTE == 0 && iter != 0) {
			transporte.iteracion(DELTA_T_POROS * PASO_TRANSPORTE);
		}

		if (iter % PASO_CONSOLA == 0 && iter != 0) {
			int interv = (clock() - reloj) / (CLOCKS_PER_SEC / 1000);
			double deltaReal = (double)interv / PASO_CONSOLA * 1000;
			printf("%.0fms %.4f %.0fus/it\n", time*1e3, poros.getRadioMaximo(), deltaReal);

			if (iter % PASO_DISCO == 0) {
				EntradaSalida::grabarTransporte(*this, time);
			}

			reloj = clock();
		}
		iter++;
	}
}

void Celula::actualizarSigmas() {

}
