#include <cmath>
#include <cassert>
#include <algorithm>
#include <iostream> //
#include "Poros.h"
#include "EntradaSalida.h"

inline bool operator<(const Poros::InfoAngulo& lhs, const Poros::InfoAngulo& rhs){
	return lhs.tita < rhs.tita;
};

Poros::Poros(Celula& celula) {
	assert(celula.nodpel == 4);
	const int NODPEL = 4;

	_celula = &celula;
	densidadPromedio = 0;

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
				double tita = 0;
				int indice = 0;
				InfoAngulo info;

				for (int i = 0; i < NODPEL; i++) {
					int iNod = elem[i];
					Nodo nodo = celula.getNodos()[iNod];
					if (esNodoExterno(nodo)) {
						tita += getTita(nodo);
						info.nodosExternos[indice] = iNod;
						indice++;
					}
				}
				assert(indice == 2);
				tita /= 2;

				info.densidad = DENSIDAD_INICIAL;
				info.radios = RADIO_INICIAL;
				info.tita = tita;
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
	Celula& celula = *_celula;
	Double2D center;
	center.x = 0;
	center.y = celula.alto / 2;

	double radio = sqrt(pow(nodo.x - center.x, 2) + pow(nodo.y - center.y, 2));
	return (abs(radio - (celula.radio + celula.ancho)) < TOLER_DIST);
}

bool Poros::esNodoInterno(Nodo nodo) {
	Celula& celula = *_celula;
	Double2D center;
	center.x = 0;
	center.y = celula.alto / 2;

	double radio = sqrt(pow(nodo.x - center.x, 2) + pow(nodo.y - center.y, 2));
	return (abs(radio - celula.radio) < TOLER_DIST);
}

double Poros::getTita(Nodo nodo) {
	Celula& celula = *_celula;
	Double2D center;
	center.x = 0;
	center.y = celula.alto / 2;

	double radio = sqrt(pow(nodo.x - center.x, 2) + pow(nodo.y - center.y, 2));
	double tita = acos((nodo.y - center.y) / radio);
	assert(tita == tita);
	return tita;
}

void Poros::iteracion() {
	densidadPromedio = 0;

	for (vector<InfoAngulo>::iterator it = valores.begin(); it != valores.end(); ++it) {
		InfoAngulo& info = *it;
		double itv1 = getCelula().getSolucion()[info.nodosExternos[0]] -
				getCelula().getSolucion()[info.nodosInternos[0]];

		double itv2 = getCelula().getSolucion()[info.nodosExternos[1]] -
				getCelula().getSolucion()[info.nodosInternos[1]];

		double itv = (itv1 + itv2) / 2;
		double n_eq = DENSIDAD_EQ * exp(CONST_Q * pow(itv / V_EP, 2));
		info.densidad = DELTA_T_TRANSPORTE * ALPHA * exp(pow((itv / V_EP), 2)) * (1 - info.densidad / n_eq) + info.densidad;
		densidadPromedio += info.densidad;

//		TODO calular radio de los poros
	}

	densidadPromedio /= valores.size();
}

void Poros::loop() {
	const double T_FINAL = 200e-3;
	EntradaSalida::printStart("Poros...", false);

	int it = 0;
	for (double time = 0; time <= T_FINAL; time += DELTA_T_POROS) {
		if (it % 10000 == 0) {
			cout << time << "\t" << densidadPromedio << endl;
		}

		iteracion();

		it++;
	}

//	for (vector<InfoAngulo>::iterator it = valores.begin();
//			it != valores.end(); ++it) {
//		InfoAngulo info = *it;
//		cout << info.tita << ",\t" << info.densidad << "\n";
//	}

	EntradaSalida::printEnd(3, false);
}
