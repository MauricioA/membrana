#define _CRT_SECURE_NO_WARNINGS

#include <cassert>
#include <iostream>
#include "Celula.h"
#include "EntradaSalida.h"
#include "Armado.h"
#include "Poisson.h"
#include "TransporteAreas.h"
#include "Poros.h"

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

void actualizarSigmas(Celula& celula, Poros& poros) {
	const double SIGMA_PORO = 150e-9;	// sigma de la sol que llena el poro

	for (int i = 0; i < celula.nElems; i++) {
		Elemento& elem = celula.getElementos()[i];
		if (elem.material == MEMBRANA) {
			double areas = poros.getProporsionArea(i);
			elem.sigma = celula.sigmas[MEMBRANA] * (1 - areas) + SIGMA_PORO * areas;
		}
	}
}

//TODO: poisson 
void Celula::transportePoros() {
	const double TIEMPO_FINAL	 = 1e-3;
	const int	 PASO_DISCO_TRANS= 300;
	const int	 PASO_DISCO_PORO = 500;
	const int 	 PASO_CONSOLA	 = 100;
	const int	 PASO_TRANSPORTE = 100;

	Poros poros = Poros(*this);
	TransporteAreas transporte = TransporteAreas(*this, poros);
	
	int iter = 0;
	clock_t reloj = 0;
	double deltaT = 1e-9;

	for (double time = 0; time < TIEMPO_FINAL; time += deltaT) {
		Poisson::poisson(*this, false);

		poros.iteracion(deltaT, time);

		actualizarSigmas(*this, poros);

		if (iter % PASO_TRANSPORTE == 0 && iter != 0) {
			transporte.iteracion(deltaT * PASO_TRANSPORTE);
		}

		if (iter % PASO_CONSOLA == 0 && iter != 0) {
			int interv = (clock() - reloj) / (CLOCKS_PER_SEC / 1000);
			double deltaReal = (double)interv / PASO_CONSOLA;
			printf("%.1fus %.4e %d %d %.0fms/it\n", 
				time*1e6, poros.getRadioMaximo(), poros.getNPoros(), poros.getNPorosChicos(), deltaReal);

			if (iter % PASO_DISCO_TRANS == 0) {
				EntradaSalida::grabarTransporte(*this, time, false);
			}
			if (iter % PASO_DISCO_PORO == 0) {
				EntradaSalida::grabarRadio(*this, poros, time, false);
			}

			reloj = clock();
		}
		iter++;
	}
}
