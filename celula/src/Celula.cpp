#include <cassert>
#include <chrono>
#include <iostream>
#include "Celula.h"
#include "Transporte.h"
#include "EntradaSalida.h"
#include "Armado.h"
#include "Poisson.h"
#include "Poros.h"

using namespace std;

Celula::Celula() {
	potencial = 0;
	nNodes = nElems = nodpel  = 0;
	alto = radio = ancho = 0;
	time = 0;
	_entradaSalida = unique_ptr<EntradaSalida>(new EntradaSalida(*this));
	area = 4 * M_PI * pow(radio, 2);
}

Celula::~Celula() {}

void Celula::poisson() {
	Poisson::iteracion(*this);
	getEntradaSalida().grabarPoisson();
}

double inline getDeltaT(double time_pulso, double multiplier, double delta_t) {
	/* 90% del máximo en 1ms */
	const double K = -1e-3 / log(0.1);
	return multiplier * delta_t * (1 - exp(-time_pulso / K));
}

/* Loop principal */
void Celula::transportePoros() {
	auto global_start = chrono::high_resolution_clock::now();
	const double THRESH_PULSE = 2e-3;
	const double MULTIPLIER_POISSON = 500;
	const double MULTIPLIER_TRANSPORTE = 1000;
	const double MULTIPLIER_CONSOLA = 5000;
	const double MULTIPLIER_ITV = 1000;
	const double PASO_DISCO = 100e-6;

	Transporte transporte(*this, calcularPoros);
	unique_ptr<Poros> poros;

	if (calcularPoros) {
		poros = unique_ptr<Poros>(new Poros(*this));
		transporte._poros = poros.get();
	}

	auto console_start = chrono::high_resolution_clock::now();
	double last_transporte = 0, last_consola = 0;
	double next_poisson, next_transporte, next_disco, next_consola, next_itv;

	for (int pulso = 0; pulso < pulsos; pulso++) {

		for (int estado = ON; estado < 2; estado++) {

			next_poisson = next_disco = next_itv = 0;	//Corren desde la 1º iteración
			next_transporte = next_consola = delta_t;	//No corren en la 1º iteración
			
			/* Si comienza nuevo pulso le aviso a poros */
			if (estado == ON && calcularPoros) poros->nuevoPulso();

			/* Si terminó el pulso le aviso a poisson */
			if (estado == OFF) Poisson::apagar(*this);

			/* Itero todo el pulso */
			for (double pulse_time = 0; pulse_time < times[estado]; pulse_time += delta_t, time += delta_t) {

				/* Iteración poisson, solo si está ON */
				if ((estado == ON) && time >= next_poisson) {
					Poisson::iteracion(*this);
					double step = getDeltaT(pulse_time, MULTIPLIER_POISSON, delta_t);
					next_poisson = time + step;
					
					/* Si no hay poros correr poisson una sola vez */
					if (!calcularPoros) next_poisson = DBL_MAX;
				}

				/* Iteración poros */
				poros->iteracion(delta_t);

				/* Iteración transporte */
				if (time >= next_transporte) {
					transporte.iteracion(time - last_transporte);
					double step = getDeltaT(pulse_time, MULTIPLIER_TRANSPORTE, delta_t);
					last_transporte = time;
					next_transporte = time + step;
				}

				/* Escribo por consola */
				if (time >= next_consola) {
					auto console_end = chrono::high_resolution_clock::now();
					auto interv = chrono::duration_cast<chrono::microseconds>(console_end - console_start);
					double delta_ms = interv.count() / 1000. / ((time - last_consola) / delta_t);
					console_start = console_end;

					if (calcularPoros) {
						printf("%.1fus %d error: %e  %.2f ms/it\n", time*1e6,
							poros->getNPoros() + poros->getNPorosChicos(), transporte.error, delta_ms);
					} else {
						printf("%.1fus error: %e  %.2f ms/it\n", time*1e6, transporte.error, delta_ms);
					}

					if (chequearValoresExtremos()) printf("VALORES EXTREMOS DE pH!!\n");
					double step = getDeltaT(pulse_time, MULTIPLIER_CONSOLA, delta_t);
					last_consola = time;
					next_consola = time + step;
				}

				/* Grabo disco poisson poros y transporte */
				if (time >= next_disco) {
					getEntradaSalida().grabarPoisson();
					getEntradaSalida().grabarTransporte();
					if (calcularPoros) getEntradaSalida().grabarPoros(*poros);
					next_disco = time + PASO_DISCO;
				}

				/* Grabo disco itv */
				if (calcularPoros && time >= next_itv) {
					getEntradaSalida().grabarITV(*poros);
					getEntradaSalida().grabarPermeabilizacion(*poros);
					double step = getDeltaT(pulse_time, MULTIPLIER_ITV, delta_t);
					next_itv = time + step;
				}
			}

		}

	}

	auto global_end = chrono::high_resolution_clock::now();
	auto interv = chrono::duration_cast<chrono::seconds>(global_end - global_start);

	cout << "Tiempo total: " << interv.count() / 60 / 60 << "h ";
	cout << (interv.count() / 60) % 60 << "m " << interv.count() % 60 << "s\n";
}

/* Retorna true si hay nodos con valores extremos de pH */
bool Celula::chequearValoresExtremos() {
	for (int node = 0; node < nNodes; node++) {
		double pH  = -log10((concs[H_][node] + 1e-18) / CONCENT);
		double pOH = -log10((concs[OH][node] + 1e-18) / CONCENT);
		if (pH < 0 || pH > 14 || pOH < 0 || pOH > 14) return true;
	}
	return false;
}
