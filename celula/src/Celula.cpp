#include <cassert>
#include <chrono>
#include <iostream>
#include "Celula.h"
#include "Transporte.h"
#include "EntradaSalida.h"
#include "Armado.h"
#include "Poisson.h"
#include "Poros.h"

#include <memory>

using namespace std;

Celula::Celula() {
	potencial = 0;
	nNodes = nElems = nodpel  = 0;
	alto = radio = ancho = 0;
	_entradaSalida = unique_ptr<EntradaSalida>(new EntradaSalida(*this));
	area = 4 * M_PI * pow(radio, 2);
}

Celula::~Celula() {}

void Celula::poisson() {
	Poisson::iteracion(*this);
	getEntradaSalida().grabarPoisson();
}

/* Loop principal */
//TODO mejorar como son los delta_t
void Celula::transportePoros() {
	double DELTA_T = 1e-9;
	
	const int PASO_POROS_1 = 1;
	const int PASO_POROS_2 = 8;
	const int PASO_POROS_3 = 128;
	const int PASO_POROS_4 = 2048;

	const int PASO_POISSON_1 = 1;
	const int PASO_POISSON_2 = 10;
	const int PASO_POISSON_3 = 100;
	const int PASO_POISSON_4 = 200;

	const int PASO_TRANSPORTE = 10;

	const int PASO_DISCO_ITV_1 = 10;
	const int PASO_DISCO_ITV_2 = 100;
	const int PASO_DISCO_ITV_3 = 1000;
	
	const int PASO_CONSOLA = 1000;

	auto global_start = chrono::high_resolution_clock::now();

	vector<double> paso_disco = { 0, };
	for (double t = 100e-6; t < 100e-3; t += 100e-6) paso_disco.push_back(t - 2*DELTA_T);

	int pos_disco = 0;
	time = 0;

	Transporte transporte(*this, calcularPoros);
	unique_ptr<Poros> poros;
	
	if (calcularPoros) {
		poros = unique_ptr<Poros>(new Poros(*this));
		transporte._poros = poros.get();
	}
	
	for (int pulso = 0; pulso < pulsos; pulso++) {

		for (int estado = ON; estado < 2; estado++) {

			/* Si comienza nuevo pulso le aviso a poros */
			if (estado == ON && calcularPoros) poros->nuevoPulso();

			/* Si terminó el pulso le aviso a poisson */
			if (estado == OFF) Poisson::apagar(*this);

			/* Inicializo variables */
			auto start = chrono::high_resolution_clock::now();
			int paso_disco_itv = PASO_DISCO_ITV_1;
			int paso_poisson = PASO_POISSON_1;
			int paso_poros = PASO_POROS_1;
			int it_trans = PASO_TRANSPORTE;
			int it_consola = 0;
			int it_disco_itv = paso_disco_itv-1;
			int it_poiss = paso_poisson;
			int it_poros = paso_poros;
			int fase = 1;
			double start_pulso = time;

			/* Itero todo el pulso */
			for (; time - start_pulso < times[estado]; time += DELTA_T) {
				
				/* Iteración poisson, solo si está ON */
				if ((estado == ON) && (it_poiss == paso_poisson)) {
					Poisson::iteracion(*this);
					if (calcularPoros) it_poiss = 0;		//si no hay poros solo correr poisson 1 vez
				}

				/* Iteración poros */
				if (it_poros == paso_poros && calcularPoros) {
					poros->iteracion(DELTA_T * paso_poros);	//el intervalo está mal cuando cambia de etapa!
					it_poros = 0;
				}

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

					if (calcularPoros) {
						printf("%.1fus %d error: %e  %.2f ms/it\n", time*1e6, 
							poros->getNPoros() + poros->getNPorosChicos(), transporte.error, delta_ms);
					} else {
						printf("%.1fus error: %e  %.2f ms/it\n", time*1e6, transporte.error, delta_ms);
					}

					it_consola = 0;
				}

				/* Grabo disco poisson poros y transporte */
				if (time >= paso_disco[pos_disco]) {
					getEntradaSalida().grabarPoisson();
					getEntradaSalida().grabarTransporte();
					if (calcularPoros) getEntradaSalida().grabarPoros(*poros);
					pos_disco++;
				}

				/* Grabo disco itv si está ON (y calcula poros) */
				if (calcularPoros && (estado == ON) && it_disco_itv == paso_disco_itv) {
					getEntradaSalida().grabarITV(*poros);
					getEntradaSalida().grabarPermeabilizacion(*poros);
					it_disco_itv = 0;
				}

				/* Actualizo pasos */
				if (fase == 1 && time > 10e-6) {
					fase = 2;
					paso_poisson = PASO_POISSON_2;
					paso_disco_itv = PASO_DISCO_ITV_2;
					paso_poros = PASO_POROS_2;
				} else if (fase == 2 && time > 100e-6) {
					fase = 3;
					paso_poisson = PASO_POISSON_3;
					paso_poros = PASO_POROS_3;
					paso_disco_itv = PASO_DISCO_ITV_3;
				} else if (fase == 3 && time > 5e-3) {
					fase = 4;
					paso_poros = PASO_POROS_4;
					paso_poisson = PASO_POISSON_4;
				}

				it_consola++; it_trans++, it_poiss++; it_disco_itv++; it_poros++;
			}

		}

	}

	auto global_end = chrono::high_resolution_clock::now();
	auto interv = chrono::duration_cast<chrono::seconds>(global_end - global_start);

	cout << "Tiempo total: " << interv.count() / 60 / 60 << "h " 
		<< (interv.count()/60) % 60 << "m " << interv.count() % 60 << "s\n";
}
