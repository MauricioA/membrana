#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS

#include <math.h>
#include <cassert>
#include <algorithm>
#include <cstdio> 
#include "Poros.h"
#include "EntradaSalida.h"

const bool   CALCULAR_RADIOS		= true;
const bool	 CALCULAR_CAPACITANCIA	= true;
const bool	 POROS_CHICOS			= true;

// um!! J -> 1e12, N -> 1e6 - esta testeado
const double DENSIDAD_INICIAL	= 0;
const double RADIO_INICIAL 		= 0.51e-3;		// r* 0.51 nm
const double RADIO_MIN_ENERGIA	= 0.80e-3;		// rm 0.80 nm
const double ALPHA				= 1e-3;			// Coeficiente de creación 1e9 m**-2 s**-1
const double V_EP				= 0.258;		// Voltaje característico [V]
const double DENSIDAD_EQ		= 1.5e-3;		// N0 Densidad de poros en equilibrio 1.5e9 m**-2
const double DIFF_POROS			= 0.05;			// D Coeficiente de diffusión para poros 5e-14 m**2 s**-1
const double F_MAX				= 0.7e-3;		// Max fuerza electrica [N V**-2]
const double R_H				= 0.97e-3;		// 0.97e-9 m
const double R_T				= 0.31e-3;		// 0.31e-9 m
const double BETA 				= 1.4e-19;		// Repulsión estérica [J]
const double GAMA 				= 1.8e-5;		// 1.8e-11 J m**-1
const double SIGMA_P			= 2e-2;			// 2e-2 J m**-2
const double SIGMA_0			= 1e-6;			// 1e-6 J m**-2
const double SIGMA_PORO			= 2e-6;			// sigma de la sol que llena el poro (2 S/m)
const double TEMPERATURA 		= 310;			// 37ºC
const double CAPACITANCIA		= 1e-14;		// 1e-2 F m**-2
const double BOLTZMANN			= 1.3806488e-11;// cte de Boltzmann 1.3806488e-23 [J K**-1]
const double CONST_Q			= pow((RADIO_MIN_ENERGIA / RADIO_INICIAL), 2);
const double TERM_TENSION_LINEA = - 2 * M_PI * GAMA;
const double TOLER_DIST_POROS	= 1e-3;
const double TOLER_ANGULO		= 1e-3;

inline bool operator<(const Poros::InfoAngulo& lhs, const Poros::InfoAngulo& rhs) {
	return lhs.tita < rhs.tita;
};

Poros::Poros(Celula& celula) {
	assert(celula.nodpel == 4);
	const int NODPEL = 4;

	_celula = &celula;
	tau = celula.radio * CAPACITANCIA * (1 / celula.sigmas[INTERNO] + 1 / (2 * celula.sigmas[EXTERNO]));
	factorPulso = 1;

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
				assert(info.area == info.area && info.area > 0);
				info.porosChicos = 0;
				info.radioChico = RADIO_INICIAL;
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

void Poros::iteracion(double deltaT, double tiempo) {
	double areaPoros = 0;

	for (auto& info : valores) {
		for (auto& poro : info.porosGrandes) {
			areaPoros += M_PI * pow(poro.first, 2);
		}

		areaPoros += info.porosChicos * M_PI * pow(info.radioChico, 2);
	}

	if (CALCULAR_CAPACITANCIA) {
		factorPulso = 1 - exp(-tiempo / tau);
	}

	double tensionEfectiva = 2 * SIGMA_P - (2 * SIGMA_P - SIGMA_0) / pow(1 - areaPoros / getCelula().area, 2);

	for (auto& info : valores) {
		double itv = getITV(info);
		bool neg = false;

		if (CALCULAR_RADIOS) {
			/* Actualizo los radios de los poros grandes */
			for (auto& poro : info.porosGrandes) {
				poro.first = actualizarRadio(poro.first, deltaT, tensionEfectiva, itv);
				if (poro.first < 0) neg = true;
			}

			/* Acá hay gato encerrado */
			if (neg) for (uint i = 0; i < info.porosGrandes.size(); i++) {
				if (info.porosGrandes[i].first < 0) {
					info.porosGrandes.erase(info.porosGrandes.begin() + i);
					i--;
					printf("PORO NEGATIVO BORRADO!\n");
				}
			}

			if (POROS_CHICOS) {
				/* Muevo a chicos los poros grandes con poco radio y bastante antigüedad */
				for (uint i = 0; i < info.porosGrandes.size(); i++) {
					auto poro = info.porosGrandes[i];
					if ((poro.first < 1e-3) && (tiempo - poro.second > 1e-6)) {
						info.porosGrandes.erase(info.porosGrandes.begin() + i);
						info.porosChicos++;
						i--;
					}
				}

				/* Actualizo el radio de los poros chicos */
				if (info.porosChicos > 0) {
					info.radioChico = actualizarRadio(info.radioChico, deltaT, tensionEfectiva, itv);
				}
			}
		}

		/* Calculo densidad */
		double n_eq = DENSIDAD_EQ * exp(CONST_Q * pow(itv / V_EP, 2));
		info.densidad = deltaT * ALPHA * exp(pow(itv / V_EP, 2)) * (1 - info.densidad / n_eq) + info.densidad;

		int porosNuevos = getPorosEnTita(info) - info.porosGrandes.size() - info.porosChicos;

		/* Agrego poros nuevos */
		for (int i = 0; i < porosNuevos; i++) {
			info.porosGrandes.push_back({RADIO_INICIAL, tiempo});
		}

		/* Si se sellaron algunos poros, borro de los poros chicos */
		if (porosNuevos < 0) {
			if (info.porosChicos >= -porosNuevos) {
				info.porosChicos -= -porosNuevos;
			} else {
				assert(false);
			}
		}
	}

	/* Actualizo las conductividades de la membrana */
	actualizarSigmas();
}

double Poros::getITV(InfoAngulo& info) {
	double itv1 = getCelula().getSolucion()[info.nodosExternos[0]] -
		getCelula().getSolucion()[info.nodosInternos[0]];

	double itv2 = getCelula().getSolucion()[info.nodosExternos[1]] -
		getCelula().getSolucion()[info.nodosInternos[1]];

	return factorPulso * (itv1 + itv2) / 2;
}

double inline Poros::actualizarRadio(double radio, double deltaT, double tensionEfectiva, double itv) {
	return radio + deltaT * DIFF_POROS / (BOLTZMANN * TEMPERATURA) * (
		(pow(itv, 2) * F_MAX) / (1 + R_H / (radio + R_T)) +
		4 * BETA * pow(RADIO_INICIAL / radio, 4) * (1 / radio) +
		TERM_TENSION_LINEA +
		2 * M_PI * tensionEfectiva * radio
	);
}

int Poros::getPorosEnTita(InfoAngulo& info) {
	return (int) (info.area * info.densidad);
}

int Poros::getNPorosChicos() {
	int kPoros = 0;
	for (auto& info : valores) {
		kPoros += info.porosChicos;
	}
	return kPoros;
}

int Poros::getNPorosGrandes() {
	int kPoros = 0;
	for (auto& info : valores) {
		kPoros += info.porosGrandes.size();
	}
	return kPoros;
}

int Poros::getNPoros() {
	return getNPorosChicos() + getNPorosGrandes();
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

	for (auto& info : valores) {
		for (auto& poro : info.porosGrandes) {
			rMax = max(poro.first, rMax);
		}
		rMax = max(info.radioChico, rMax);
	}

	return rMax;
}

/* Radio del primer poro del depolarizado */
double Poros::getUnRadio() {
	InfoAngulo& info = valores[valores.size() - 1];
	if (info.porosGrandes.size() > 0) {
		return info.porosGrandes[0].first;
	} else {
		return 0;
	}
}

/* Fracción del area correspondiente a un elemento de la membrana ocupada por poros */
double Poros::getProporsionArea(int iElem) {
	InfoAngulo& info = InfoAngulo(*mapa[iElem]);
	double areaP = 0;
	for (auto& poro : info.porosGrandes) {
		areaP += M_PI * pow(poro.first, 2);
	}
	areaP += info.porosChicos * M_PI * pow(info.radioChico, 2);

	/* Si hubo algún error */
	if (areaP > info.area) {
		for (auto poro : info.porosGrandes) {
			printf("%e %e\n", poro.first, poro.second);
		}
		double itv = getITV(info);
		printf("%f\n", itv);
		printf("%f %f %d %d %f\n", info.area, areaP, info.porosGrandes.size(), info.porosChicos, info.tita);
	}

	assert(areaP < info.area);
	return areaP / info.area;
}

vector<pair<double, double>> Poros::getITVs(double tiempo) {
	vector<pair<double, double>> result;
	for (auto& info : valores) {
		result.push_back(pair<double, double> { info.tita, getITV(info) });
	}
	return result;
}

void Poros::actualizarSigmas() {
	for (int i = 0; i < getCelula().nElems; i++) {
		Elemento& elem = getCelula().getElementos()[i];
		if (elem.material == MEMBRANA) {
			double areas = getProporsionArea(i);
			elem.sigma = getCelula().sigmas[MEMBRANA] * (1 - areas) + SIGMA_PORO * areas;
		}
	}
}
