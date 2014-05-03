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
		assert(abs(matriz.coeff(i, j) - matriz.coeff(j, i)) < 1e-12);
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

// TODO ver valores iniciales y sigma de la sol
void Celula::transportePoros() {
	const double TIEMPO_FINAL	= 20e-3;
	const double DELTA_T		= 1e-9;

	const int PASO_POISSON_1		= 1;
	const int PASO_POISSON_2		= 10;
	const int PASO_POISSON_3		= 50;
	const int PASO_DISCO_POISSON_1	= 10;
	const int PASO_DISCO_POISSON_2	= 1000;
	const int PASO_DISCO_ITV_1		= 10;
	const int PASO_DISCO_ITV_2		= 100;
	const int PASO_DISCO_PORO_1		= 10;
	const int PASO_DISCO_PORO_2		= 1000;
	const int PASO_TRANSPORTE		= 100;
	const int PASO_CONSOLA			= 1000;
	const int PASO_DISCO_TRANSPORTE = 10000;

	int paso_disco_poisson = PASO_DISCO_POISSON_1;
	int paso_disco_itv = PASO_DISCO_ITV_1;
	int paso_disco_poro = PASO_DISCO_PORO_1;
	int paso_poisson = PASO_POISSON_1;

	Poros poros = Poros(*this);
	TransporteAreas transporte = TransporteAreas(*this, poros);
	
	int it_consola = 0, it_disco_trans = 0, it_trans = 0, it_disco_itv = paso_disco_itv; 
	int it_poiss = paso_poisson, it_disco_poro = paso_disco_poro, it_disco_poisson = 0;
	clock_t reloj = 0;
	int fase = 0;

	for (double time = 0; time < TIEMPO_FINAL; time += DELTA_T) {
		
		/* Imprimo por consola */
		if (it_consola == PASO_CONSOLA) {
			int interv = (clock() - reloj) / (CLOCKS_PER_SEC / 1000);
			double deltaReal = (double)interv / PASO_CONSOLA;
			reloj = clock();
			printf("%.1fus %.2e %d %d  %.2f ms/it\n", time*1e6, poros.getRadioMaximo(), 
				poros.getNPoros(), poros.getNPorosChicos(), deltaReal);
			it_consola = 0;
		}
	
		/* Grabo disco poisson */
		if (it_disco_poisson == paso_disco_poisson) {
			EntradaSalida::grabarPoisson(*this, time);
			it_disco_poisson = 0;
		}

		/* Grabo disco itv */
		if (it_disco_poisson == paso_disco_itv) {
			EntradaSalida::grabarITV(*this, poros, time);
			it_disco_poisson = 0;
		}

		/* Grabo disco poros */
		if (it_disco_poro == paso_disco_poro) {
			EntradaSalida::grabarRadio(*this, poros, time, false);
			it_disco_poro = 0;
		}

		/* Grabo disco transporte */
		if (it_disco_trans == PASO_DISCO_TRANSPORTE) {
			EntradaSalida::grabarTransporte(*this, time, false);
			it_disco_trans = 0;
		}
		
		if (it_poiss == paso_poisson) {
			Poisson::poisson(*this);
			it_poiss = 0;
		}

		poros.iteracion(DELTA_T, time);

		actualizarSigmas(*this, poros);

		if (it_trans == PASO_TRANSPORTE) {
			transporte.iteracion(DELTA_T * PASO_TRANSPORTE);
			it_trans = 0;
		}

		/* Actualizo pasos */
		if (fase == 0 && time > 10e-6) {
			fase = 1;
			paso_poisson = PASO_POISSON_2;
			paso_disco_poro = PASO_DISCO_PORO_2;
			paso_disco_poisson = PASO_DISCO_POISSON_2;
			paso_disco_itv = PASO_DISCO_ITV_2;
		} else if (fase == 1 && time > 100e-6) {
			fase = 3;
			paso_poisson = PASO_POISSON_3;
		}

		it_consola++; it_disco_trans++; it_disco_poro++; 
		it_trans++, it_poiss++; it_disco_poisson++; it_disco_itv++;
	}
}
