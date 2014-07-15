#include <cassert>
#include <chrono>
#include <iostream>
#include "Celula.h"
#include "EntradaSalida.h"
#include "Armado.h"
#include "Poisson.h"
#include "TransporteAreas.h"
#include "Poros.h"

using namespace std;

Celula::Celula() {
	potencial = 0;
	nNodes = nElems = nodpel  = 0;
	alto = radio = ancho = 0;
	_entradaSalida = make_unique<EntradaSalida>(*this);
	area = 4 * M_PI * pow(radio, 2);
}

Celula::~Celula() {}

void Celula::poisson() {
	Poisson::poisson(*this);
	getEntradaSalida().grabarPoisson();
}

/* Loop principal */
void Celula::transportePoros() {
	const double TIEMPO_FINAL	= 20e-3;
	const double DELTA_T		= 1e-9;

	const int PASO_POISSON_1		= 1;
	const int PASO_POISSON_2		= 10;
	const int PASO_POISSON_3		= 50;
	const int PASO_DISCO_ITV_1		= 10;
	const int PASO_DISCO_ITV_2		= 100;
	const int PASO_DISCO_ITV_3		= 1000;
	const int PASO_TRANSPORTE		= 1000;
	const int PASO_CONSOLA			= 1000;

	vector<double> paso_disco = { 0, 5e-3, 10e-3, 15e-3, 19.999e-3, 1 };
	int pos_disco = 0;

	int paso_disco_itv = PASO_DISCO_ITV_1;
	int paso_poisson = PASO_POISSON_1;

	Poros poros = Poros(*this);

	TransporteAreas transporte = TransporteAreas(*this, poros);
	
	int it_consola = 0, it_trans = 0, it_disco_itv = paso_disco_itv-1; 
	int it_poiss = paso_poisson;
	auto start = chrono::high_resolution_clock::now();
	int fase = 1;

	for (time = 0; time < TIEMPO_FINAL; time += DELTA_T) {
		
		/* Iteración poisson */
		if (it_poiss == paso_poisson) {
			Poisson::poisson(*this);
			it_poiss = 0;
		}

		/* Iteración poros */
		poros.iteracion(DELTA_T);

		/* Iteración transporte */
		if (it_trans == PASO_TRANSPORTE) {
			transporte.iteracion(DELTA_T * PASO_TRANSPORTE);
			it_trans = 0;
		}

		/* Imprimo por consola */
		if (it_consola == PASO_CONSOLA) {
			auto end = chrono::high_resolution_clock::now();
			long long delta_um = chrono::duration_cast<chrono::microseconds>(end - start).count();
			double deltaReal = delta_um / 1000. / PASO_CONSOLA;
			start = end;

			printf("%.1fus %.2e %d %d  %.2f ms/it\n", time*1e6, poros.getRadioMaximo(),
				poros.getNPoros(), poros.getNPorosChicos(), deltaReal);
			it_consola = 0;
		}

		/* Grabo disco poisson */
		if (time >= paso_disco[pos_disco]) {
			getEntradaSalida().grabarPoisson();
			getEntradaSalida().grabarTransporte();
			getEntradaSalida().grabarPoros(poros);
			pos_disco++;
		}

		/* Grabo disco itv */
		if (it_disco_itv == paso_disco_itv) {
			getEntradaSalida().grabarITV(poros);
			it_disco_itv = 0;
		}

		/* Actualizo pasos */
		if (fase == 1 && time > 10e-6) {
			fase = 2;
			paso_poisson = PASO_POISSON_2;
			paso_disco_itv = PASO_DISCO_ITV_2;
		} else if (fase == 2 && time > 100e-6) {
			fase = 3;
			paso_poisson = PASO_POISSON_3;
		} else if (fase == 3 && time > 5e-3) {
			fase = 4;
			paso_disco_itv = PASO_DISCO_ITV_3;
		}

		it_consola++; it_trans++, it_poiss++; it_disco_itv++;
	}
}
