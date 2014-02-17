#include <cmath>
#include <cassert>
#include <algorithm>
#include <iostream> //
#include <cstdio>	//
#include "Poros.h"
#include "EntradaSalida.h"

/* Constantes */
const double DENSIDAD_INICIAL	= 0;
const double RADIO_INICIAL 		= 510e-6;		// r* 0.51 nm
const double RADIO_MIN_ENERGIA	= 800e-6;		// rm 0.80 nm
const double ALPHA				= 1e-3;			// Coeficiente de creación 1e9 m**-2 s**-1
const double V_EP				= 0.258;		// Voltaje característico [V]
const double DENSIDAD_EQ		= 1.5e-3;		// N0 Densidad de poros en equilibrio 1.5e9 m**-2
const double CONST_Q			= pow((RADIO_MIN_ENERGIA / RADIO_INICIAL), 2);
const double DIFF_POROS			= 50e-3;		// D Coeficiente de diffusión para poros 5e-14 m**2 s**-1
const double DELTA_T_POROS		= 100e-9;
const double F_MAX				= 0.7e-9;		// Max fuerza electrica [N V**-2]
const double R_H				= 970e-6;		// 0.97e-9 m
const double R_T				= 310e-6;		// 0.31e-9 m
const double BETA 				= 1.4e-19;		// Repulsión estérica [J]
const double GAMA 				= 1.8e-17;		// 1.8e-11 J m**-1
const double SIGMA_P			= 2e-14;		// 2e-2 J m**-2
const double SIGMA_0			= 1e-18;		// 1e-6 J m**-2
const double TEMPERATURA 		= 310;			// 37ºC

inline bool operator<(const Poros::InfoAngulo& lhs, const Poros::InfoAngulo& rhs){
	return lhs.tita < rhs.tita;
};

Poros::Poros(Celula& celula) {
	assert(celula.nodpel == 4);
	const int NODPEL = 4;

	_celula = &celula;

	for (int iElem = 0; iElem < celula.nElems; iElem++) {
		Elemento elem = celula.getElementos()[iElem];

		/* Chequeo que sea de la membrana y del borde externo */
		if (elem.material == MEMBRANA) {
			bool externo = false;
			for (int i = 0; !externo && i < NODPEL; i++) {
				int iNod = elem[i];
				externo = esNodoExterno(celula.getNodos()[iNod]);
			}

			if (externo) {
				double titas[2];
				int indice = 0;
				InfoAngulo info;

				for (int i = 0; i < NODPEL; i++) {
					int iNod = elem[i];
					Nodo nodo = celula.getNodos()[iNod];
					if (esNodoExterno(nodo)) {
						titas[indice] = getTita(nodo);
						info.nodosExternos[indice] = iNod;
						indice++;
					}
				}
				assert(indice == 2);

				info.densidad = DENSIDAD_INICIAL;
				info.radios = RADIO_INICIAL;
				info.tita = (titas[0] + titas[1]) / 2;
				info.deltaTita = abs(titas[0] - titas[1]);
				info.sin = sin(info.tita);
				valores.push_back(info);
			}
		}
	}
	assert(valores.size() > 0);

	/* Busco los nodos internos de los valores */
	uint ints = 0;
	for (int iElem = 0; iElem < celula.nElems; iElem++) {
		Elemento elem = celula.getElementos()[iElem];

		/* Chequeo que sea de la membrana y del borde interno */
		if (elem.material == MEMBRANA) {
			bool interno = false;
			for (int i = 0; !interno && i < NODPEL; i++) {
				int iNod = elem[i];
				interno = esNodoInterno(celula.getNodos()[iNod]);
			}

			if (interno) {
				double tita = 0;
				int indice = 0;
				int internos[2];

				for (int i = 0; i < NODPEL; i++) {
					int iNod = elem[i];
					Nodo nodo = celula.getNodos()[iNod];
					if (esNodoInterno(nodo)) {
						tita += getTita(nodo);
						internos[indice] = iNod;
						indice++;
					}
				}
				assert(indice == 2);
				tita /= 2;

				/* Busco el info del angulo */
				bool encontrado = false;
				for (uint i = 0; !encontrado && i < valores.size(); i++) {
					if (abs(valores[i].tita - tita) < TOLER_DIST) {
						valores[i].nodosInternos[0] = internos[0];
						valores[i].nodosInternos[1] = internos[1];
						encontrado = true;
					}
				}
				assert(encontrado);
				ints++;
			}
		}
	}
	assert(ints == valores.size());

	sort(valores.begin(), valores.end());
}

inline Celula& Poros::getCelula() {
	return *_celula;
}

bool Poros::esNodoExterno(Nodo nodo) {
	Double2D center = getCelula().getCenter();
	double radio = sqrt(pow(nodo.x - center.x, 2) + pow(nodo.y - center.y, 2));
	return (abs(radio - (getCelula().radio + getCelula().ancho)) < TOLER_DIST);
}

bool Poros::esNodoInterno(Nodo nodo) {
	Double2D center = getCelula().getCenter();
	double radio = sqrt(pow(nodo.x - center.x, 2) + pow(nodo.y - center.y, 2));
	return (abs(radio - getCelula().radio) < TOLER_DIST);
}

double Poros::getTita(Nodo nodo) {
	Double2D center = getCelula().getCenter();
	double radio = sqrt(pow(nodo.x - center.x, 2) + pow(nodo.y - center.y, 2));
	double tita = acos((nodo.y - center.y) / radio);
	assert(tita == tita);
	return tita;
}

void Poros::iteracion() {
	double kPoros = 0;
	double areaPoros = 0;

	for (vector<InfoAngulo>::iterator it = valores.begin(); it != valores.end(); ++it) {
		InfoAngulo info = *it;
		double poros = info.densidad * info.sin * info.deltaTita * 2 * M_PI * pow(getCelula().radio, 2);
		kPoros += poros;
		areaPoros += M_PI * pow(info.radios, 2) * poros;
	}

	for (vector<InfoAngulo>::iterator it = valores.begin(); it != valores.end(); ++it) {
		InfoAngulo& info = *it;
		double itv1 = getCelula().getSolucion()[info.nodosExternos[0]] -
				getCelula().getSolucion()[info.nodosInternos[0]];

		double itv2 = getCelula().getSolucion()[info.nodosExternos[1]] -
				getCelula().getSolucion()[info.nodosInternos[1]];

		double itv = (itv1 + itv2) / 2;

		double n_eq = DENSIDAD_EQ * exp(CONST_Q * pow(itv / V_EP, 2));
		info.densidad = DELTA_T_POROS * ALPHA * exp(pow((itv / V_EP), 2)) * (1 - info.densidad / n_eq) + info.densidad;

		if (kPoros >= 1) {
			double termElectrico = (pow(itv, 2) * F_MAX) / (1 + R_H / (info.radios + R_T));
			double termRepulsion = 4 * BETA * pow(RADIO_INICIAL / info.radios, 4) * (1 / info.radios);
			double termTensionLinea = - 2 * M_PI * GAMA;
			double tensionEfectiva = 2 * SIGMA_P - (2 * SIGMA_P - SIGMA_0) / pow(1 - areaPoros / getCelula().area, 2);
			double termTensionSuperficial = 2 * M_PI * tensionEfectiva * info.radios;

			info.radios += DELTA_T_POROS * DIFF_POROS / (kPoros * TEMPERATURA) *
					(termElectrico + termRepulsion + termTensionLinea + termTensionSuperficial);
		}
	}
}

double Poros::getNPoros() {
	double integral = 0;

	for (vector<InfoAngulo>::iterator it = valores.begin(); it != valores.end(); ++it) {
		InfoAngulo info = *it;
		integral += info.densidad * info.sin * info.deltaTita;
	}

	return integral * 2 * M_PI * pow(getCelula().radio, 2);
}

double Poros::getDensidadPromedio() {
	double densidad = 0;

	for (vector<InfoAngulo>::iterator it = valores.begin(); it != valores.end(); ++it) {
		InfoAngulo info = *it;
		densidad += info.densidad;
	}

	return densidad / valores.size();
}

void Poros::loop() {
	const double T_FINAL = 3;
	const int PASO_CONSOLA = 500000;

	EntradaSalida::printStart("Poros...", false);
	clock_t reloj = 0;
	int it = 0;

	for (double time = 0; time <= T_FINAL; time += DELTA_T_POROS) {
		if (it % PASO_CONSOLA == 0) {
			int interv = (clock() - reloj) / (CLOCKS_PER_SEC / 1000);
			printf("%f, %e, %f, %.0f us/it\n", time, getDensidadPromedio(), getNPoros(), (double)interv / PASO_CONSOLA * 1000);
			fflush(stdout);
			reloj = clock();
		}

		iteracion();

		it++;
	}

	for (vector<InfoAngulo>::iterator it = valores.begin();
			it != valores.end(); ++it) {
		InfoAngulo info = *it;
		printf("%f, %e\n", info.tita, info.radios);
	}

	EntradaSalida::printEnd(3, false);
}
