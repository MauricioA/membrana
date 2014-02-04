#include <cassert>
#include <vector>
#include <Eigen/Sparse>
#include "Transporte.h"
#include "declaraciones.h"
#include "Poisson.h"
#include "Armado.h"
#include "EntradaSalida.h"

using namespace declaraciones;
using namespace Eigen;
using namespace std;

//TODO transporte
//TODO calcular error
//TODO FORT: armado_t algo raro en esm con transporte

void Transporte::transporte(Celula& cel) {
	const double T_CERO = 1;

	for (int esp = 0; esp < NESPS; esp++) {
		cel.concentraciones[esp].resize(cel.nNodes);
		cel.anteriores[esp].resize(cel.nNodes);
	}

	cel.phAux[H_].resize(cel.nNodes);
	cel.phAux[OH].resize(cel.nNodes);

	for (int iNode = 0; iNode < cel.nNodes; iNode++) {
		for (int esp = 0; esp < NESPS; esp++) {
			cel.concentraciones[esp][iNode] = CONCENTRACION_INICIAL[esp];
			cel.anteriores[esp][iNode] = CONCENTRACION_INICIAL[esp];
		}

//		phAux[H_][iNode] = -log10(concentraciones[H_][iNode] * 1e15 / 6.02e23);
//		phAux[OH][iNode] = -log10(concentraciones[OH][iNode] * 1e15 / 6.02e23);
	}

	masaDiag2D(cel);
	cel.getCargas().resize(cel.nNodes);

	for (double tt = 0.; tt < T_CERO; tt += DELTA_T) {
		EntradaSalida::printStart("Iteración transporte...");

		Poisson::poisson(cel, false);

		for (int esp = 0; esp < NESPS; esp++) {
			concentracion(cel, esp);
		}

		for (int kNodo = 0; kNodo < cel.nNodes; kNodo++) {
			cel.phAux[H_][kNodo] = -log10(cel.concentraciones[H_][kNodo] * 1e15 / 6.02e23);
			cel.phAux[OH][kNodo] = -log10(cel.concentraciones[OH][kNodo] * 1e15 / 6.02e23);

			for (int esp = 0; esp < NESPS; esp++) {
				cel.concentraciones[esp][kNodo] = RSA * cel.concentraciones[esp][kNodo] + (1-RSA) * cel.anteriores[esp][kNodo];
				if (cel.concentraciones[esp][kNodo] < CONCENT_MINIMO) {
					cel.concentraciones[esp][kNodo] = 0.;
				}
				cel.anteriores[esp][kNodo] = cel.concentraciones[esp][kNodo];
			}
		}

//		TODO escribir por salida
		EntradaSalida::printEnd(1);
	}
}

void Transporte::masaDiag2D(Celula& cel) {
	assert(cel.nodpel == 4);

	cel.getMasas().resize(cel.nNodes);
	for (int i = 0; i < cel.nNodes; i++) cel.getMasas()[i] = 0;

	const int NLOCS = 2;
	const int INOGA[] = {0, 3, 1, 2};
	const double POSGL[] = {-1.0, 1.0};
	const double WEIGL[] = { 1.0, 1.0};

	for (int iElem = 0; iElem < cel.nElems; iElem++) {
		Elemento elem = cel.getElementos()[iElem];
		Nodo nodosElem[cel.nodpel];
		for (int i = 0; i < cel.nodpel; i++) nodosElem[i] = cel.getNodos()[elem[i]];

		/*	subroutine armotodo(nope,x,y,deriv,weigc) */
		int iGauss = 0;
		double weigc[cel.nodpel];
		double posgc[2][cel.nodpel];
		double deriv[2][cel.nodpel][cel.nodpel];

		for (int ilocs = 0; ilocs < NLOCS; ilocs++) {
			for (int jlocs = 0; jlocs < NLOCS; jlocs++) {
				weigc[INOGA[iGauss]] = WEIGL[ilocs] * WEIGL[jlocs];
				posgc[0][INOGA[iGauss]] = POSGL[ilocs];
				posgc[1][INOGA[iGauss]] = POSGL[jlocs];
				iGauss++;
			}
		}

		for (int i = 0; i < cel.nodpel; i++) {
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
		/* end subroutine armotodo */

		for (int iNode = 0; iNode < cel.nodpel; iNode++) {
			double aJacob[2][2];
			for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++) aJacob[i][j] = 0.0;

			for (int jNode = 0; jNode < cel.nodpel; jNode++) {
				aJacob[0][0] += nodosElem[jNode].x * deriv[0][jNode][iNode];
				aJacob[0][1] += nodosElem[jNode].x * deriv[1][jNode][iNode];
				aJacob[1][0] += nodosElem[jNode].y * deriv[0][jNode][iNode];
				aJacob[1][1] += nodosElem[jNode].y * deriv[1][jNode][iNode];
			}

			double gpDet = aJacob[0][0] * aJacob[1][1] - aJacob[0][1] * aJacob[1][0];
			double gpVol = weigc[iNode] * gpDet;
			cel.getMasas()[elem[iNode]] += gpVol;
		}

	}

	for (int iNode = 0; iNode < cel.nNodes; iNode++)	{
		assert(cel.getMasas()[iNode] > TOLER_MASA);
	}
}

void Transporte::carga(Celula& cel) {
	const double CTE = 1e6 / 6.03e23;
	const double CTE_DILUCION = CONCENTRACION_INICIAL[H_] / CONCENTRACION_INICIAL[NA];

	for (int iNodo = 0; iNodo < cel.nNodes; iNodo++) {
		cel.getCargas()[iNodo] = FARADAY / (EPSILON_TRANSPORTE * EPSILON_0) * CTE * (
			CARGA[H_] * cel.concentraciones[H_][iNodo] +
			CARGA[OH] * cel.concentraciones[OH][iNodo] +
			CARGA[NA] * cel.concentraciones[NA][iNodo] * CTE_DILUCION +
			CARGA[CL] * cel.concentraciones[CL][iNodo] * CTE_DILUCION
		);
	}
}

void Transporte::concentracion(Celula& cel, int esp) {
	double esm[cel.nodpel][MAXNPEL];
	vector< Triplet<double> > triplets;

	for (uint kElem = 0; kElem < cel.getElementos().size(); kElem++) {
		Elemento elem = cel.getElementos()[kElem];
		Double2D pos[cel.nodpel];
		double sol[cel.nodpel];
		double ef[cel.nodpel];
		double ch_Med = 0.;
		double cohMed = 0.;

		for (int i = 0; i < cel.nodpel; i++) {
			int iNodo = elem[i];
			Nodo nodo = cel.getNodos()[iNodo];
			pos[i].x = nodo.x;
			pos[i].y = nodo.y;
			sol[i] = cel.getSolucion()[iNodo];

			if (esp == OH) {
				ch_Med += cel.anteriores[H_][iNodo];
				cohMed += cel.anteriores[OH][iNodo];
			}
		}

		if (esp == OH) {
			ch_Med /= cel.nodpel;
			cohMed /= cel.nodpel;
		}

		double sigma = cel.sigmas[elem.material];

		double difElem = DIFUSION[esp];
		if (elem.material == MEMBRANA) difElem *= 1e-3;

		double mu = -difElem * CLAVE * CARGA[esp] * cel.getGradElem()[kElem].y;
		double landa = 1.;
		double qe = 0.;

		if (esp == OH) {
			qe = KWB * CONCENT_H2O - KWF * ch_Med * cohMed;
		}

		Armado::armadoTransporte(pos, sigma, qe, landa, mu, sol, esm, ef);

		for (int i = 0; i < cel.nodpel; i++) {
			int iNodo = elem[i];
			Nodo nodo = cel.getNodos()[iNodo];

			if (nodo.esTierra) {
				double adiag = esm[i][i];
				for (int j = 0; j < cel.nodpel; j++) {
					esm[i][j] = 0.;
					ef[j] -= esm[j][i] * CONCENTRACION_CATODO[esp];
					esm[j][i] = 0.;
				}
				esm[i][i] = adiag;
				ef[i] = adiag * CONCENTRACION_CATODO[esp];
			}

			if (nodo.esPotencia) {
				double adiag = esm[i][i];
				for (int j = 0; j < cel.nodpel; j++) {
					esm[i][j] = 0.;
					ef[j] -= esm[j][i] * CONCENTRACION_ANODO[esp];
					esm[j][i] = 0.;
				}
				esm[i][i] = adiag;
				ef[iNodo] = adiag * CONCENTRACION_ANODO[esp];
			}
		}

		/* Ensamblado */
		for (int i = 0; i < cel.nNodes; i++) cel.getRhs()[i] = 0.;
		triplets.clear();

		for (int i = 0; i < cel.nodpel; i++) {
			int iNodo = elem[i];
			cel.getRhs()[iNodo] += ef[i];

			for (int j = 0; j < cel.nodpel; j++) {
				int jNodo = elem[j];
				triplets.push_back(Triplet<double>(iNodo, jNodo, esm[i][j]));
			}
		}
	}

	cel.getMatriz().resize(cel.nNodes, cel.nNodes);
	cel.getMatriz().setFromTriplets(triplets.begin(), triplets.end());

	/* Resolución */
//	TODO usar LU
	BiCGSTAB< SparseMatrix<double> > solver;

//	TODO usar precondicionador
	solver.compute(cel.getMatriz());

//	concentraciones[esp] = solver.solve(rhs);
	VectorXd guess = cel.concentraciones[esp];
	cel.concentraciones[esp] = solver.solveWithGuess(cel.getRhs(), guess);
}
