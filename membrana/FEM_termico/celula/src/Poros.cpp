#define _USE_MATH_DEFINES
#include <math.h>
#include <cassert>
#include <algorithm>
#include <cstdio> 	//
#include "Poros.h"
#include "EntradaSalida.h"

/* Constantes */
const double DENSIDAD_INICIAL	= 0;
const double RADIO_INICIAL 		= 510e-6;		// r* 0.51 nm
const double RADIO_MIN_ENERGIA	= 800e-6;		// rm 0.80 nm
const double ALPHA				= 1e-3;			// Coeficiente de creaci�n 1e9 m**-2 s**-1
const double V_EP				= 0.258;		// Voltaje caracter�stico [V]
const double DENSIDAD_EQ		= 1.5e-3;		// N0 Densidad de poros en equilibrio 1.5e9 m**-2
const double CONST_Q			= pow((RADIO_MIN_ENERGIA / RADIO_INICIAL), 2);
const double DIFF_POROS			= 50e-3;		// D Coeficiente de diffusi�n para poros 5e-14 m**2 s**-1
const double DELTA_T_POROS		= 100e-9;
const double F_MAX				= 0.7e-9;		// Max fuerza electrica [N V**-2]
const double R_H				= 970e-6;		// 0.97e-9 m
const double R_T				= 310e-6;		// 0.31e-9 m
const double BETA 				= 1.4e-19;		// Repulsi�n est�rica [J]
const double GAMA 				= 1.8e-17;		// 1.8e-11 J m**-1
const double SIGMA_P			= 2e-14;		// 2e-2 J m**-2
const double SIGMA_0			= 1e-18;		// 1e-6 J m**-2
const double TEMPERATURA 		= 310;			// 37�C
const double TERM_TENSION_LINEA = - 2 * M_PI * GAMA;

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
				info.tita = (titas[0] + titas[1]) / 2;
				double tita1 = min(titas[0], titas[1]);
				double tita2 = max(titas[0], titas[1]);
				info.constIntegral = 2 * M_PI * pow(getCelula().radio, 2) * (cos(tita1) - cos(tita2));
				assert(info.constIntegral == info.constIntegral);
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
	int kPoros = 0;
	double areaPoros = 0;

	for (uint i = 0; i < valores.size(); i++) {
		InfoAngulo& info = valores[i];
		kPoros += info.poros.size();

		for (uint j = 0; j < info.poros.size(); j++) {
			double radio = info.poros[j];
			areaPoros += M_PI * pow(radio, 2);
		}
	}

	double tensionEfectiva = 2 * SIGMA_P - (2 * SIGMA_P - SIGMA_0) / pow(1 - areaPoros / getCelula().area, 2);

	for (uint i = 0; i < valores.size(); i++) {
		InfoAngulo& info = valores[i];
		double itv1 = getCelula().getSolucion()[info.nodosExternos[0]] -
				getCelula().getSolucion()[info.nodosInternos[0]];

		double itv2 = getCelula().getSolucion()[info.nodosExternos[1]] -
				getCelula().getSolucion()[info.nodosInternos[1]];

		double itv = (itv1 + itv2) / 2;

		double n_eq = DENSIDAD_EQ * exp(CONST_Q * pow(itv / V_EP, 2));
		info.densidad = DELTA_T_POROS * ALPHA * exp(pow(itv / V_EP, 2)) * (1 - info.densidad / n_eq) + info.densidad;

		for (uint j = 0; j < info.poros.size(); j++) {
			double& radio = info.poros[j];
			double termElectrico = (pow(itv, 2) * F_MAX) / (1 + R_H / (radio + R_T));
			double termRepulsion = 4 * BETA * pow(RADIO_INICIAL / radio, 4) * (1 / radio);
			double termTensionSuperficial = 2 * M_PI * tensionEfectiva * radio;

			radio += DELTA_T_POROS * DIFF_POROS / (kPoros * TEMPERATURA) *
					(termElectrico + termRepulsion + TERM_TENSION_LINEA + termTensionSuperficial);

			assert((radio == radio) && radio < 1);
		}

		int porosNuevos = getPorosEnTita(info) - info.poros.size();
		assert(porosNuevos >= 0);

		for (int i = 0; i < porosNuevos; i++) {
			info.poros.push_back(RADIO_INICIAL);
		}
	}
}

int Poros::getPorosEnTita(InfoAngulo& info) {
	return (int) (info.constIntegral * info.densidad);
}

int Poros::getNPoros() {
	int kPoros = 0;
	for (uint i = 0; i < valores.size(); i++) {
		kPoros += valores[i].poros.size();
	}
	return kPoros;
}

double Poros::getDensidadPromedio() {
	double densidad = 0;

	for (uint i = 0; i < valores.size(); i++) {
		densidad += valores[i].densidad;
	}

	return densidad / valores.size();
}

void Poros::loop() {
	const double T_FINAL = 3;
	const int PASO_CONSOLA = 100000;

	EntradaSalida::printStart("Poros...", false);

	printf("tita, itv\n");
	for (uint i = 0; i < valores.size(); i++) {
		InfoAngulo info = valores[i];
		double itv1 = getCelula().getSolucion()[info.nodosExternos[0]] -
				getCelula().getSolucion()[info.nodosInternos[0]];

		double itv2 = getCelula().getSolucion()[info.nodosExternos[1]] -
				getCelula().getSolucion()[info.nodosInternos[1]];

		double itv = (itv1 + itv2) / 2;
		printf("%f, %f\n", info.tita, itv);
	}

	fflush(stdout);

	long reloj = 0;
	int it = 0;

	for (double time = 0; time <= T_FINAL; time += DELTA_T_POROS) {
		if (it % PASO_CONSOLA == 0) {
			//int interv = (clock() - reloj) / (CLOCKS_PER_SEC / 1000);
			//printf("%.2e, %e, %d, %.0fus/it\n", time, getDensidadPromedio(), getNPoros(), (double)interv / PASO_CONSOLA * 1000);
			printf("%.2e, %e, %d\n", time, getDensidadPromedio(), getNPoros());
//			printf("%.2e, %.0f us/it\n", time, (double)interv / PASO_CONSOLA * 1000);
			fflush(stdout);
			//reloj = clock();
		}

		iteracion();

		it++;
	}

	printf("tita, densidad, nPoros\n");
	for (vector<InfoAngulo>::iterator it = valores.begin();	it != valores.end(); ++it) {
		InfoAngulo info = *it;
		printf("%f, %e, %d\n", info.tita, info.densidad, info.poros.size());
	}

	EntradaSalida::printEnd(3, false);
}
