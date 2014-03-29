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

void Celula::transportePoros() {
	const double TIEMPO_FINAL	 = 20e-3;
	const int	 PASO_DISCO_TRANS= 2000;
	const int	 PASO_DISCO_PORO = 500;
	const int 	 PASO_CONSOLA	 = 500;
	const int	 PASO_TRANSPORTE = 100;
	const int	 PASO_POISSON	 = 10;

	Poros poros = Poros(*this);
	TransporteAreas transporte = TransporteAreas(*this, poros);
	
	int itCons = 0, itTransD = 0, itPoroD = PASO_DISCO_PORO, itTrans = 0, itPoiss = PASO_POISSON;
	clock_t reloj = 0;
	double deltaT = 1e-9;

	for (double time = 0; time < TIEMPO_FINAL; time += deltaT) {
		
		/* Imprimo por cosola */
		if (itCons == PASO_CONSOLA) {
			int interv = (clock() - reloj) / (CLOCKS_PER_SEC / 1000);
			double deltaReal = (double)interv / PASO_CONSOLA;
			reloj = clock();
			printf("%.1fus %.4e %d %d %.3e  %.1f ms/it\n",	time*1e6, poros.getRadioMaximo(), 
				poros.getNPoros(), poros.getNPorosChicos(), deltaT, deltaReal);
			itCons = 0;
		}
		
		/* Grabo disco transporte */
		if (itTransD == PASO_DISCO_TRANS) {
			EntradaSalida::grabarTransporte(*this, time, false);
			itTransD = 0;
		}

		/* Grabo disco poros */
		if (itPoroD == PASO_DISCO_PORO) {
			EntradaSalida::grabarRadio(*this, poros, time, false);
			itPoroD = 0;
		}
		
		if (itPoiss == PASO_POISSON) {
			Poisson::poisson(*this, false);
			itPoiss = 0;
		}

		poros.iteracion(deltaT, time);

		actualizarSigmas(*this, poros);

		if (itTrans == PASO_TRANSPORTE) {
			transporte.iteracion(deltaT * PASO_TRANSPORTE);
			itTrans = 0;
		}

		/* Agrando deltaT */
		//if (time > 30e-6 && deltaT < 0.1e-6) deltaT += 25e-15;

		itCons++; itTransD++; itPoroD++; itTrans++, itPoiss++;
	} //tita 2.86
}
