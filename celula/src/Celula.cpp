#include <cassert>
#include <chrono>
#include <iostream>
#include <cfloat>
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
	Poros poros(*this);				//necesario para tener los ángulos, etc
	getEntradaSalida().grabarPoisson();
	getEntradaSalida().grabarITV(poros);
}

double inline getDeltaT(double time_pulso, double multiplier, double delta_t) {
	/* 90% del máximo en 500e-6 */
	const double K = 500e-6 / log(0.1);
	return multiplier * delta_t * (1 - exp(time_pulso / K));
}

/* Loop principal */
//TODO imprimir por pantalla cada una cantidad de tiempo real (chequear cada x iters)
void Celula::acoplado() {
	auto global_start = chrono::high_resolution_clock::now();
	const double MULTIPLIER_POISSON = 2000;	
	const double MULTIPLIER_TRANSPORTE = 200;	//ponerlo en 20 por lo menos con poros
	const double MULTIPLIER_ITV = 1000;
	const double PASO_DISCO = 100e-6;

	const int CONSOLE_UPDATE_CHECK = 10000;
	const int CONSOLE_INTERVAL_US = 5 * 1000 * 1000;

	Transporte transporte(*this, calcularPoros);
	unique_ptr<Poros> poros;

	if (calcularPoros) {
		poros = unique_ptr<Poros>(new Poros(*this));
		transporte._poros = poros.get();
	}

	auto console_start = chrono::high_resolution_clock::now();
	double last_transporte = 0, last_consola = 0;
	double next_poisson, next_transporte, next_disco, next_itv;
	double multiplier_transporte = MULTIPLIER_TRANSPORTE;

	int console_counter = 0;

	for (int pulso = 0; pulso < pulsos; pulso++) {

		for (int estado = ON; estado < 2; estado++) {

			next_poisson = next_disco = next_itv = 0;	//Corren desde la 1º iteración
			next_transporte = delta_t;					//No corre en la 1º iteración
			Poisson::estado = estado;					//Actualizo valor 
			
			/* Si comienza nuevo pulso le aviso a poros (para que actualize capacitancia) */
			if (estado == ON && calcularPoros) poros->nuevoPulso();

			/* Itero todo el pulso */
			for (double pulse_time = 0; pulse_time < times[estado]; pulse_time += delta_t, time += delta_t) {

				/* Iteración poisson, solo si está ON */
				if (time >= next_poisson) {
					Poisson::iteracion(*this);
					double step = getDeltaT(pulse_time, MULTIPLIER_POISSON, delta_t);
					next_poisson = time + step;
					
					/* Si no hay poros correr poisson una sola vez */
					if (!calcularPoros) next_poisson = DBL_MAX;
				}

				/* Iteración poros */
				if (calcularPoros) poros->iteracion(delta_t);

				/* Iteración transporte */
				if (time >= next_transporte) {
					transporte.iteracion(time - last_transporte);
					last_transporte = time;
					next_transporte = time + multiplier_transporte * delta_t;
				}

				/* Escribo por consola */
				if (console_counter == CONSOLE_UPDATE_CHECK) {
					console_counter = 0;
					auto time_now = chrono::high_resolution_clock::now();
					auto interval_us = chrono::duration_cast<chrono::microseconds>(time_now - console_start);

					if (interval_us.count() > CONSOLE_INTERVAL_US) {
						double intervalo_virtual = time - last_consola;
						double iteraciones = intervalo_virtual / delta_t;
						double tiempo_por_iter = interval_us.count() / iteraciones;

						if (calcularPoros) {
							printf("%.1fus %d error: %e  %.6f us/it", time*1e6,
								poros->getNPoros() + poros->getNPorosChicos(), transporte.error, tiempo_por_iter);
						} else {
							printf("%.1fus error: %e  %.6f us/it", time*1e6, transporte.error, tiempo_por_iter);
						}

						int nodosExtremos = valoresExtremos();
						if (nodosExtremos > 10) {
							assert(false && "poros extremos muchos");
							printf("  pH ext: %d", nodosExtremos);
						}
						printf("\n");

						console_start = time_now;
					}	
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

				console_counter++;
			}

		}

	}

	auto global_end = chrono::high_resolution_clock::now();
	auto interv = chrono::duration_cast<chrono::seconds>(global_end - global_start);

	cout << "Tiempo total: " << interv.count() / 60 / 60 << "h ";
	cout << (interv.count() / 60) % 60 << "m " << interv.count() % 60 << "s\n";
}

int Celula::valoresExtremos() {
	int count = 0;
	for (int node = 0; node < nNodes; node++) {
		double pH  = -log10((concs[H_][node] + 1e-18) / CONCENT);
		double pOH = -log10((concs[OH][node] + 1e-18) / CONCENT);
		if (pH < 0 || pH > 14 || pOH < 0 || pOH > 14) {
			count++;
		}
	}
	return count;
}
