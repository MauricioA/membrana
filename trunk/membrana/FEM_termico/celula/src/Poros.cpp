#include "Poros.h"
#include <cmath>
#include <cassert>


Poros::Poros(Celula& celula) {
	assert(celula.nodpel == 4);

	const double TOLER_DIST = 1e-3;
	const int NODPEL = 4;

	vector<int> elementosInternos, elementosExternos;
	vector<double> angulosInternos, angulosExternos;
	Double2D center;
	center.x = 0;
	center.y = celula.alto / 2;

	for (int iNodo = 0; iNodo < celula.nNodes; iNodo++) {
		Nodo nodo = celula.getNodos()[iNodo];
		double radio = sqrt(pow(nodo.x - center.x, 2) + pow(nodo.y - center.y, 2));
		double tita = acos(nodo.y - center.y) / celula.radio;

		if (abs(radio - celula.radio) < TOLER_DIST) {
			elementosInternos.push_back(iNodo);
			angulosInternos.push_back(tita);
		} else if (abs(radio - (celula.radio + celula.ancho)) < TOLER_DIST) {
			elementosExternos.push_back(iNodo);
			angulosExternos.push_back(tita);
		}

		assert(elementosExternos.size() == elementosInternos.size());
		assert(elementosExternos.size() > 0);
	}

	for (int iElem = 0; iElem < NODPEL; iElem++) {
		Elemento elem = celula.getElementos()[iElem];
		if (elem.material == MEMBRANA) {
			double tita1, tita2;

			bool ambos = false;
			tita1 = acos(elem[0].y - center.y) / celula.radio;
			for (int i = 1; (!ambos) && i < NODPEL; i++) {
				double tita = acos(elem[i].y - center.y) / celula.radio;
				if (abs(tita - tita1) > TOLER_DIST) {
					tita2 = tita;
					ambos = true;
				}
			}
			assert(ambos);

			int pos;
			bool anguloEncontrado = false;
			for (int i = 0; !(anguloEncontrado) && i < angulosInternos.size(); i++) {
				if (abs(angulosInternos[i] - tita1) < TOLER_DIST) {
					pos = i;
					anguloEncontrado = true;
				}
			}
			assert(anguloEncontrado);

			//hacer lo mismo para externos y para tita2


//			tita = math.degrees(math.acos((y - centerY) / radio))
//
//			if abs(radio - radioCel) < tolerancia:
//				internos.append((tita, v))
//			elif abs(radio - (radioCel + dMemb)) < tolerancia:
//				externos.append((tita, v))


			ElementoMembrana info;

		}
	}


}

void Poros::iteracion() {

}
