#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS

#include <math.h>
#include <cassert>
#include <algorithm>
#include <cstdio> 
#include "Poros.h"
#include "EntradaSalida.h"

//TODO poros chicos

const double DELTA_T_POROS	= 1.5e-9;	// [s]
const double T_FINAL		= 1e-3;
const int	 PASO_CONSOLA_P	= 100000;
const bool   RADIOS			= true;     // calcular radios

// todo en metros!
const double DENSIDAD_INICIAL	= 0;
const double RADIO_INICIAL 		= 0.51e-9;		// r* 0.51 nm
const double RADIO_MIN_ENERGIA	= 0.80e-9;		// rm 0.80 nm
const double ALPHA				= 1e9;			// Coeficiente de creación 1e9 m**-2 s**-1
const double V_EP				= 0.258;		// Voltaje característico [V]
const double DENSIDAD_EQ		= 1.5e9;		// N0 Densidad de poros en equilibrio 1.5e9 m**-2
const double CONST_Q			= pow((RADIO_MIN_ENERGIA / RADIO_INICIAL), 2);
const double DIFF_POROS			= 50e-14;		// D Coeficiente de diffusión para poros 5e-14 m**2 s**-1
const double F_MAX				= 0.7e-9;		// Max fuerza electrica [N V**-2]
const double R_H				= 0.97e-9;		// 0.97e-9 m
const double R_T				= 0.31e-9;		// 0.31e-9 m
const double BETA 				= 1.4e-19;		// Repulsión estérica [J]
const double GAMA 				= 1.8e-11;		// 1.8e-11 J m**-1
const double SIGMA_P			= 2e-2;			// 2e-2 J m**-2
const double SIGMA_0			= 1e-6;			// 1e-6 J m**-2
const double TEMPERATURA 		= 310;			// 37ºC
const double TERM_TENSION_LINEA = - 2 * M_PI * GAMA;
const double BOLTZMANN			= 1.3806488e-23;// cte de Boltzmann [J K**-1]
const double TOLER_DIST_POROS	= 1e-9;
const double TOLER_ANGULO		= 1e-3;

inline bool operator<(const Poros::InfoAngulo& lhs, const Poros::InfoAngulo& rhs) {
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
				info.area = 2 * M_PI * pow(getCelula().radio, 2) * (cos(tita1) - cos(tita2));
				assert(info.area == info.area);
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
			
			double tita = getTita(elem);

			/* Busco el info del angulo */
			bool encontrado = false;
			for (uint i = 0; !encontrado && i < valores.size(); i++) {
				if (abs(valores[i].tita - tita) < TOLER_ANGULO) {
					encontrado = true;

					/* Lo agrego al mapa */
					mapa[iElem] = &valores[i];

					/* Si es interno agrego los nodos al info */
					if (interno) {
						int indice = 0;
						int internos[2];

						for (int j = 0; j < NODPEL; j++) {
							int iNod = elem[j];
							Nodo nodo = celula.getNodos()[iNod];
							if (esNodoInterno(nodo)) {
								internos[indice] = iNod;
								indice++;
							}
						}
						assert(indice == 2);

						valores[i].nodosInternos[0] = internos[0];
						valores[i].nodosInternos[1] = internos[1];
						ints++;
					}
				}
			}
			assert(encontrado);
		}
	}
	assert(ints == valores.size());

	sort(valores.begin(), valores.end());
}

inline Celula& Poros::getCelula() {
	return *_celula;
}

double Poros::getTita(Nodo nodo) {
	Double2D center = getCelula().getCenter();
	double radio = sqrt(pow(nodo.x - center.x, 2) + pow(nodo.y - center.y, 2));
	double tita = acos((nodo.y - center.y) / radio);
	assert(tita == tita);
	return tita;
}

double Poros::getTita(Elemento elem) {
	const int NODPEL = 4;
	double tita = 0;
	for (int i = 0; i < NODPEL; i++) {
		tita += getTita(getCelula().getNodos()[elem[i]]);
	}
	return tita / NODPEL;
}

bool Poros::esNodoExterno(Nodo nodo) {
	Double2D center = getCelula().getCenter();
	double radio = sqrt(pow(nodo.x - center.x, 2) + pow(nodo.y - center.y, 2));
	return (abs(radio - (getCelula().radio + getCelula().ancho)) < TOLER_DIST_POROS);
}

bool Poros::esNodoInterno(Nodo nodo) {
	Double2D center = getCelula().getCenter();
	double radio = sqrt(pow(nodo.x - center.x, 2) + pow(nodo.y - center.y, 2));
	return (abs(radio - getCelula().radio) < TOLER_DIST_POROS);
}

/* sacar time */
void Poros::iteracion() {
	double areaPoros = 0;

	for (auto& info : valores) for (auto& radio : info.poros) {
		areaPoros += M_PI * pow(radio, 2);
	}

	double tensionEfectiva = 2 * SIGMA_P - (2 * SIGMA_P - SIGMA_0) / pow(1 - areaPoros / getCelula().area, 2);

	for (auto& info : valores) {
		double itv;

		double itv1 = getCelula().getSolucion()[info.nodosExternos[0]] -
			getCelula().getSolucion()[info.nodosInternos[0]];

		double itv2 = getCelula().getSolucion()[info.nodosExternos[1]] -
			getCelula().getSolucion()[info.nodosInternos[1]];

		itv = (itv1 + itv2) / 2;

		double n_eq = DENSIDAD_EQ * exp(CONST_Q * pow(itv / V_EP, 2));
		info.densidad = DELTA_T_POROS * ALPHA * exp(pow(itv / V_EP, 2)) * (1 - info.densidad / n_eq) + info.densidad;

		if (RADIOS) for (auto& radio : info.poros) {
			radio = radio + DELTA_T_POROS * DIFF_POROS / (BOLTZMANN * TEMPERATURA) * (
				(pow(itv, 2) * F_MAX) / (1 + R_H / (radio + R_T)) +
				4 * BETA * pow(RADIO_INICIAL / radio, 4) * (1 / radio) + 
				TERM_TENSION_LINEA + 
				2 * M_PI * tensionEfectiva * radio
			);
		}

		int porosNuevos = getPorosEnTita(info) - info.poros.size();
		assert(porosNuevos >= 0);

		for (int i = 0; i < porosNuevos; i++) {
			info.poros.push_back(RADIO_INICIAL);
		}
	}
}

int Poros::getPorosEnTita(InfoAngulo& info) {
	return (int) (info.area * info.densidad);
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

double Poros::getRadioMaximo() {
	double rMax = 0;

	for (auto& info : valores) for (auto& poro : info.poros) {
		rMax = max(poro, rMax);
	}

	return rMax;
}

/* Radio del primer poro del depolarizado */
double Poros::getUnRadio() {
	InfoAngulo& info = valores[valores.size() - 1];
	if (info.poros.size() > 0) {
		return info.poros[0];
	} else {
		return 0;
	}
}

void Poros::loop() {
	EntradaSalida::printStart("Poros...", false);

	double area = 0;
	for (auto& info : valores) {
		area += info.area;
	}
	printf("area total: %e\n", area);

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

	clock_t reloj = 0;
	int it = 0;

	//FILE* filePoros = fopen("salida/poros.dat", "w");
	
	for (double time = 0; time <= T_FINAL; time += DELTA_T_POROS) {

		if (it % PASO_CONSOLA_P == 0) {
			int interv = (clock() - reloj) / (CLOCKS_PER_SEC / 1000);
			double deltaReal = (double)interv / PASO_CONSOLA_P * 1000;
			printf("%.6f, %d, %.4e, %.0fus/it\n", time, getNPoros(), getRadioMaximo(), deltaReal);
			fflush(stdout);
			
			//fprintf(filePoros, "%e %d\n", time, getNPoros());
			//for (auto& info : valores) for (auto& poro : info.poros) {
			//	fprintf(filePoros, "%e %e\n", info.tita, poro);
			//}

			reloj = clock();
		}
		
		iteracion();

		it++;
	}

	//fclose(filePoros);

	/*FILE* file;
	file = fopen("temp.csv", "w");

	for (auto  &info : valores) {
		for (uint p = 0; p < info.poros.size(); p++) {
			fprintf(file, "%e\n", info.poros[p]);
		}
	}

	fclose(file);*/
	
	EntradaSalida::printEnd(3, false);
}

/* Fracción del area correspondiente a un elemento de la membrana ocupada por poros */
double Poros::getProporsionArea(int iElem) {
	InfoAngulo& info = InfoAngulo(*mapa[iElem]);
	double areaP = 0;
	for (auto& poro : info.poros) {
		areaP += poro;
	}
	return areaP / info.area;
}
