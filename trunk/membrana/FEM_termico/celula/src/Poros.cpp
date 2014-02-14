#include "Poros.h"
#include <cmath>
#include <cassert>

Poros::Poros(Celula& celula) {
	assert(celula.nodpel == 4);
	const int NODPEL = 4;

	_celula = &celula;

	vector<int> nodosInternos, nodosExternos;
	vector<double> angulosInternos, angulosExternos;
	Double2D center;
	center.x = 0;
	center.y = celula.alto / 2;

	for (int iNodo = 0; iNodo < celula.nNodes; iNodo++) {
		Nodo nodo = celula.getNodos()[iNodo];
		double radio = sqrt(pow(nodo.x - center.x, 2) + pow(nodo.y - center.y, 2));
		double tita = acos(nodo.y - center.y) / celula.radio;

		if (abs(radio - celula.radio) < TOLER_DIST) {
			nodosInternos.push_back(iNodo);
			angulosInternos.push_back(tita);
		} else if (abs(radio - (celula.radio + celula.ancho)) < TOLER_DIST) {
			nodosExternos.push_back(iNodo);
			angulosExternos.push_back(tita);
		}

		assert(nodosExternos.size() == nodosInternos.size());
		assert(nodosExternos.size() > 0);
	}

	for (int iElem = 0; iElem < NODPEL; iElem++) {
		Elemento elem = celula.getElementos()[iElem];
		if (elem.material == MEMBRANA) {
			double tita1, tita2;
			ElementoMembrana infoMemb;

			bool ambos = false;
			int iNod = elem[0];
			tita1 = acos(celula.getNodos()[iNod].y - center.y) / celula.radio;
			for (int i = 1; (!ambos) && i < NODPEL; i++) {
				iNod = elem[i];
				double tita = acos(celula.getNodos()[iNod].y - center.y) / celula.radio;
				if (abs(tita - tita1) > TOLER_DIST) {
					tita2 = tita;
					ambos = true;
				}
			}
			assert(ambos);

			bool tita1Done = false, tita2Done = false;
			for (uint i = 0; !(tita1Done && tita2Done) && i < angulosInternos.size(); i++) {
				if (abs(angulosInternos[i] - tita1) < TOLER_DIST) {
					infoMemb.NodosInternos[0] = nodosInternos[i];
					tita1Done = true;
				} else if (abs(angulosInternos[i] - tita2) < TOLER_DIST) {
					infoMemb.NodosInternos[1] = nodosInternos[i];
					tita2Done = true;
				}
			}
			assert(tita1Done && tita2Done);
			
			tita1Done = false, tita2Done = false;
			for (uint i = 0; !(tita1Done && tita2Done) && i < angulosExternos.size(); i++) {
				if (abs(angulosExternos[i] - tita1) < TOLER_DIST) {
					infoMemb.NodosExternos[0] = nodosExternos[i];
					tita1Done = true;
				} else if (abs(angulosExternos[i] - tita2) < TOLER_DIST) {
					infoMemb.NodosExternos[1] = nodosExternos[i];
					tita2Done = true;
				}
			}
			assert(tita1Done && tita2Done);

			infoMemb.densidad = DENSIDAD_INICIAL;
			infoMemb.radios = RADIO_INICIAL;

			mapaMembrana[iElem] = infoMemb;
		}
	}
}

void Poros::iteracion() {
	for (map<int, ElementoMembrana>::iterator it = mapaMembrana.begin();
			it != mapaMembrana.end(); ++it) {
		ElementoMembrana info = it -> second;
		Celula& cel = *_celula;

		ElementoMembrana nuevo;
		nuevo.NodosInternos[0] = info.NodosInternos[0];
		nuevo.NodosInternos[1] = info.NodosInternos[1];
		nuevo.NodosExternos[0] = info.NodosExternos[0];
		nuevo.NodosExternos[1] = info.NodosExternos[1];

		double itv1 = cel.getSolucion()[info.NodosExternos[0]] - cel.getSolucion()[info.NodosInternos[0]];
		double itv2 = cel.getSolucion()[info.NodosExternos[1]] - cel.getSolucion()[info.NodosInternos[1]];
		double itv = (itv1 + itv2) / 2;

		double n_eq = DENSIDAD_EQ * exp(CONST_Q * pow(itv / V_EP, 2));
		nuevo.densidad = DELTA_T * ALPHA * exp(pow((itv / V_EP), 2)) * (1 - info.densidad / n_eq) + info.densidad;

		//TODO calcular el radio!!
		nuevo.radios = info.radios;

		it -> second = nuevo;
	}
}

void Poros::loop() {
	const double T_FINAL = 1;

//	TODO achicar delta_t!
	for (double time = 0; time <= T_FINAL; time += DELTA_T) {
		iteracion();
	}
}
