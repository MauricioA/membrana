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
	Poisson::iteracion(*this);
	getEntradaSalida().grabarPoisson();
}

/* Loop principal */
void Celula::transportePoros() {
	//const double DELTA_T  = 1e-9;
	const double DELTA_T = 0.003e-9;
	
	const int PASO_POISSON_1	= 1;
	const int PASO_POISSON_2	= 10;
	const int PASO_POISSON_3	= 50;
	const int PASO_DISCO_ITV_1	= 10;
	const int PASO_DISCO_ITV_2	= 100;
	const int PASO_DISCO_ITV_3	= 1000;
	const int PASO_TRANSPORTE	= 1000;
	const int PASO_CONSOLA		= 1000;

	auto global_start = chrono::high_resolution_clock::now();

	/* En qué momentos graba a disco */
	/*vector<double> paso_disco = { 
		0, 5e-3, 10e-3, 15e-3, 19.999e-3, 60e-3,
		100e-3, 105e-3, 110e-3, 115e-3, 119.999e-3, 160e-3,	
		200e-3, 205e-3, 210e-3, 215e-3, 219.999e-3, 260e-3,
		300e-3, 305e-3, 310e-3, 315e-3, 319.999e-3, 360e-3,
		DBL_MAX,
	};*/
	
	//vector<double> paso_disco = {
	//	0, 2e-3, 4e-3, 5.999e-3, 8e-3, 10e-3, 12e-3, 15.999e-3, DBL_MAX,
	//};

	/*vector<double> paso_disco = {
		0, 30e-6, 100e-6, .9999e-3, 1,
	};*/

	vector<double> paso_disco = {
		0, 10e-9, 20e-9, 30e-9, 40e-9, 50e-9, 60e-9, 70e-9, 80e-9, 94.99e-9, 1,
	};

	int pos_disco = 0;
	time = 0;

	Poros poros = Poros(*this);
	TransporteAreas transporte = TransporteAreas(*this, poros);
	
	for (int pulso = 0; pulso < pulsos; pulso++) {

		for (int estado = ON; estado < 2; estado++) {

			/* Si comienza nuevo pulso le aviso a poros */
			if (estado == ON) poros.nuevoPulso();

			/* Si terminó el pulso le aviso a poisson */
			if (estado == OFF) Poisson::apagar(*this);

			/* Inicializo variables */
			auto start = chrono::high_resolution_clock::now();
			int paso_disco_itv = PASO_DISCO_ITV_1;
			int paso_poisson = PASO_POISSON_1;
			int it_consola = 0, it_trans = 0, it_disco_itv = paso_disco_itv-1;
			int it_poiss = paso_poisson;
			int fase = 1;
			double start_pulso = time;

			/* Itero todo el pulso */
			for (; time - start_pulso < times[estado]; time += DELTA_T) {
				
				/* Iteración poisson, solo si está ON */
				if ((estado == ON) && (it_poiss == paso_poisson)) {
					Poisson::iteracion(*this);
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
					auto interv = chrono::duration_cast<chrono::microseconds>(end - start);
					double delta_ms = interv.count() / 1000. / PASO_CONSOLA;
					start = end;

					printf("%.1fus %.2e %d %d  %.2f ms/it\n", time*1e6, poros.getRadioMaximo(),
						poros.getNPoros(), poros.getNPorosChicos(), delta_ms);
					it_consola = 0;
				}

				/* Grabo disco poisson */
				if (time >= paso_disco[pos_disco]) {
					getEntradaSalida().grabarPoisson();
					getEntradaSalida().grabarTransporte();
					getEntradaSalida().grabarPoros(poros);
					pos_disco++;
				}

				/* Grabo disco itv si está ON */
				if ((estado == ON) && it_disco_itv == paso_disco_itv) {
					getEntradaSalida().grabarITV(poros);
					getEntradaSalida().grabarPermeabilizacion(poros);
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

	}

	auto global_end = chrono::high_resolution_clock::now();
	auto interv = chrono::duration_cast<chrono::seconds>(global_end - global_start);

	cout << "Tiempo total: " << interv.count() / 60 / 60 << "h " 
		<< (interv.count()/60) % 60 << "m " << interv.count() % 60 << "s\n";
}
