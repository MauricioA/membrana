#include <cassert>
#include <vector>
#include <iostream>
#include <ctime>
#include <Eigen/Sparse>
#include "Transporte.h"
#include "declaraciones.h"
#include "Poisson.h"
#include "Armado.h"
#include "EntradaSalida.h"

using namespace declaraciones;
using namespace Eigen;
using namespace std;

const double TOLER_MASA 	= 1e-12;
const double FARADAY 		= 96485.34;		// C/mol cte. de Faraday
const double R_CTE			= 8.314; 		// J/K/mol cte. de gases
const double T_CTE			= 310;			// K temperatura
const double RSA 			= 0.5;
const double CONCENT_MINIMO = 1e-12;

Transporte::Transporte(Celula& celula, bool calcularPoros) {
	_celula = &celula;
	_calcularPoros = calcularPoros;

	for (int esp = 0; esp < NESPS; esp++) {
		celula.concs[esp].resize(celula.nNodes);
		celula.c_ant[esp].resize(celula.nNodes);
	}

	for (Elemento elem : celula.elementos) {
		for (int k = 0; k < celula.nodpel; k++) {
			for (int esp = 0; esp < NESPS; esp++) {
				if (elem.material == MEMBRANA || elem.material == INTERNO) {
					celula.concs[esp][elem[k]] = CONC_INICIAL_INTRA[esp];
					celula.c_ant[esp][elem[k]] = CONC_INICIAL_INTRA[esp];
				} else {
					celula.concs[esp][elem[k]] = CONC_INICIAL_EXTRA[esp];
					celula.c_ant[esp][elem[k]] = CONC_INICIAL_EXTRA[esp];
				}
			}
		}
	}

	masaDiag2D();
}

inline Celula& Transporte::getCelula() {
	return *_celula;
}

void Transporte::iteracion(double deltaT) {
	double num = 0, den = 0;
	Celula& celula = getCelula();

	#pragma omp parallel for num_threads(celula.threads)
	for (int esp = 0; esp < NESPS; esp++) {
		concentracion(esp, deltaT);
	}

	for (int kNodo = 0; kNodo < celula.nNodes; kNodo++) {
		for (int esp = 0; esp < NESPS; esp++) {
			num += pow(celula.concs[esp][kNodo] - celula.c_ant[esp][kNodo], 2);
			den += pow(celula.concs[esp][kNodo], 2);

			celula.concs[esp][kNodo] =
				RSA * celula.concs[esp][kNodo] +
				(1 - RSA) * celula.c_ant[esp][kNodo];

			if (celula.concs[esp][kNodo] < CONCENT_MINIMO) {
				celula.concs[esp][kNodo] = 0;
			}

			celula.c_ant[esp][kNodo] = celula.concs[esp][kNodo];
		}
	}

	double error = sqrt(num / den);
	assert(error < 1e3 && error == error);
}

void Transporte::masaDiag2D() {
	Celula& celula = getCelula();
	masas.resize(celula.nNodes);
	for (int i = 0; i < celula.nNodes; i++) masas[i] = 0;

	const int NLOCS = 2;
	const int INOGA[] = { 0, 3, 1, 2 };
	const double POSGL[] = { -1.0, 1.0 };
	const double WEIGL[] = { 1.0, 1.0 };

	for (int iElem = 0; iElem < celula.nElems; iElem++) {
		Elemento elem = celula.elementos[iElem];
		Nodo nodosElem[MAXNPEL];
		double rMed = 0;

		for (int i = 0; i < celula.nodpel; i++) {
			nodosElem[i] = celula.nodos[elem[i]];
			rMed += nodosElem[i].x;
		}
		rMed /= celula.nodpel;

		switch (celula.nodpel) {
			case 3: {
				double b[3], c[3];
				Double2D pos[3];
				for (int i = 0; i < celula.nodpel; i++) {
					pos[i].x = nodosElem[i].x;
					pos[i].y = nodosElem[i].y;
				}

				double det = Armado::determinante3(pos, b, c);

				double gpvol = det * M_PI * rMed;
				for (int i = 0; i < celula.nodpel; i++) {
					masas[elem[i]] += gpvol;
				}

				break;
			} case 4: {
				int iGauss = 0;
				double weigc[MAXNPEL];
				double posgc[2][MAXNPEL];
				double deriv[2][MAXNPEL][MAXNPEL];

				for (int ilocs = 0; ilocs < NLOCS; ilocs++) {
					for (int jlocs = 0; jlocs < NLOCS; jlocs++) {
						weigc[INOGA[iGauss]] = WEIGL[ilocs] * WEIGL[jlocs];
						posgc[0][INOGA[iGauss]] = POSGL[ilocs];
						posgc[1][INOGA[iGauss]] = POSGL[jlocs];
						iGauss++;
					}
				}

				for (int i = 0; i < celula.nodpel; i++) {
					double s = posgc[0][i];
					double t = posgc[1][i];

					int k = 0;
					deriv[0][k++][i] = (-1.0 + t) * 0.25;
					deriv[0][k++][i] = ( 1.0 - t) * 0.25;
					deriv[0][k++][i] = ( 1.0 + t) * 0.25;
					deriv[0][k++][i] = (-1.0 - t) * 0.25;

					k = 0;
					deriv[1][k++][i] = (-1.0 + s) * 0.25;
					deriv[1][k++][i] = (-1.0 - s) * 0.25;
					deriv[1][k++][i] = ( 1.0 + s) * 0.25;
					deriv[1][k++][i] = ( 1.0 - s) * 0.25;
				}

				for (int i = 0; i < celula.nodpel; i++) {
					double aJacob[2][2];
					for (int k = 0; k < 2; k++) for (int l = 0; l < 2; l++) aJacob[k][l] = 0;

					for (int j = 0; j < celula.nodpel; j++) {
						aJacob[0][0] += nodosElem[j].x * deriv[0][j][i];
						aJacob[0][1] += nodosElem[j].x * deriv[1][j][i];
						aJacob[1][0] += nodosElem[j].y * deriv[0][j][i];
						aJacob[1][1] += nodosElem[j].y * deriv[1][j][i];
					}

					double gpDet = aJacob[0][0] * aJacob[1][1] - aJacob[0][1] * aJacob[1][0];
					double gpVol = weigc[i] * gpDet;
					
					masas[elem[i]] += gpVol * 2 * M_PI * rMed;
				}

				break;
			} default: {
				assert(false);
			}
		}
	}

	for (auto& masa : masas) {
		assert(masa > TOLER_MASA);
	}
}

void Transporte::concentracion(int esp, double deltaT) {
	Celula& celula = getCelula();
	double esm[MAXNPEL][MAXNPEL];
	vector<Triplet<double>> triplets;
	triplets.reserve(celula.elementos.size() * celula.nodpel * celula.nodpel);

	VectorXd rhs;
	rhs.resize(celula.nNodes);
	rhs.fill(0);

	for (uint kElem = 0; kElem < celula.elementos.size(); kElem++) {
		Elemento elem = celula.elementos[kElem];
		Double2D pos[MAXNPEL];
		double mas[MAXNPEL];
		double sol[MAXNPEL];
		double ef[MAXNPEL];
		double difusion;

		for (int i = 0; i < celula.nodpel; i++) {
			int iNodo = elem[i];
			Nodo nodo = celula.nodos[iNodo];
			pos[i].x = nodo.x;
			pos[i].y = nodo.y;
			sol[i] = celula.c_ant[esp][iNodo];
			mas[i] = masas[iNodo];
		}

		if (elem.material == MEMBRANA) {
			difusion = difusionMembrana(kElem, esp);
		} else {
			difusion = DIFUSION[esp];
		}

		double mu = -difusion * (FARADAY / (R_CTE * T_CTE)) * CARGA[esp] * celula.gradElem[kElem].y;

		double landa = 1;
		double qe = 0;

		Armado::armadoTransporte(celula.nodpel, pos, difusion, qe, landa, mu, mas, sol, esm, ef, deltaT);

		for (int i = 0; i < celula.nodpel; i++) {
			int iNodo = elem[i];
			Nodo nodo = celula.nodos[iNodo];

			if (nodo.esTierra) {
				double adiag = esm[i][i];
				for (int j = 0; j < celula.nodpel; j++) {
					esm[i][j] = 0;
					ef[j] -= esm[j][i] * CONC_CATODO[esp];
					esm[j][i] = 0;
				}
				esm[i][i] = adiag;
				ef[i] = adiag * CONC_CATODO[esp];
			}

			if (nodo.esPotencia) {
				double adiag = esm[i][i];
				for (int j = 0; j < celula.nodpel; j++) {
					esm[i][j] = 0;
					ef[j] -= esm[j][i] * CONC_ANODO[esp];
					esm[j][i] = 0;
				}
				esm[i][i] = adiag;
				ef[i] = adiag * CONC_ANODO[esp];
			}
		}

		for (int i = 0; i < celula.nodpel; i++) {
			int iNodo = elem[i];
			rhs[iNodo] += ef[i];

			for (int j = 0; j < celula.nodpel; j++) {
				int jNodo = elem[j];
				triplets.push_back(Triplet<double>(iNodo, jNodo, esm[i][j]));
			}
		}
	}

	SparseMatrix<double> matriz(celula.nNodes, celula.nNodes);
	matriz.setFromTriplets(triplets.begin(), triplets.end());

	BiCGSTAB<SparseMatrix<double>> solver(matriz);
	celula.concs[esp] = solver.solveWithGuess(rhs, celula.c_ant[esp]);

	assert(solver.info() == Success);
}

inline double Transporte::difusionMembrana(int iElem, int especie) {
	if (_calcularPoros) {
		return DIFUSION[especie] * (_poros->getProporsionArea(iElem));
	} else {
		return DIFUSION[especie] * 1e-3;
	}
}
