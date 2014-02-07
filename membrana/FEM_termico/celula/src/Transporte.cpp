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

//TODO transporte
//TODO calcular error
//TODO FORT: armado_t algo raro en esm con transporte

// -DNDEBUG al compilador

// TODO es siempre la misma matriz para una especie!!

void Transporte::transporte(Celula& cel) {
	const double T_FINAL = 1;

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
	}

	masaDiag2D(cel);
//	cel.getCargas().resize(cel.nNodes);
	long iter = 0;

	clock_t reloj = 0;
	for (double time = 0; time < T_FINAL; time += DELTA_T) {
		double num = 0, den = 0;

		Poisson::poisson(cel, false);

		for (int esp = 0; esp < NESPS; esp++) {
			concentracion(cel, esp);
		}

		for (int kNodo = 0; kNodo < cel.nNodes; kNodo++) {
			cel.phAux[H_][kNodo] = -log10((cel.concentraciones[H_][kNodo] + 1e-18) * 1e15 / 6.02e23);
			cel.phAux[OH][kNodo] = -log10((cel.concentraciones[OH][kNodo] + 1e-18) * 1e15 / 6.02e23);


			for (int esp = 0; esp < NESPS; esp++) {
				num += pow(cel.concentraciones[esp][kNodo] - cel.anteriores[esp][kNodo], 2);
				den += pow(cel.concentraciones[esp][kNodo], 2);

				cel.concentraciones[esp][kNodo] =
						RSA * cel.concentraciones[esp][kNodo] +
						(1-RSA) * cel.anteriores[esp][kNodo];

				if (cel.concentraciones[esp][kNodo] < CONCENT_MINIMO) {
					cel.concentraciones[esp][kNodo] = 0;
				}
				cel.anteriores[esp][kNodo] = cel.concentraciones[esp][kNodo];
			}
		}

		double error = sqrt(num/den);
		if (error > 1e3 || error != error) {
			cerr << "iter " << iter << ", error: " << error << ", den: " << den << endl;
			exit(EXIT_FAILURE);
		}

		if (iter % PASO_CONSOLA == 0) {
			if (iter != 0) {
				int interv = (clock() - reloj) / (CLOCKS_PER_SEC / 1000);
				cout << time*1e6 << "us\t"
					 << iter << " iters\t"
					 << interv/PASO_CONSOLA << " ms/it" <<  endl;
			}

			if (iter % PASO_DISCO == 0) {
				EntradaSalida::grabarTransporte(cel, time);
			}

			reloj = clock();
		}

		iter++;
	}
}

void Transporte::masaDiag2D(Celula& cel) {
	cel.getMasas().resize(cel.nNodes);
	for (int i = 0; i < cel.nNodes; i++) cel.getMasas()[i] = 0;

	const int NLOCS = 2;
	const int INOGA[] = {0, 3, 1, 2};
	const double POSGL[] = {-1.0, 1.0};
	const double WEIGL[] = { 1.0, 1.0};

	for (int iElem = 0; iElem < cel.nElems; iElem++) {
		Elemento elem = cel.getElementos()[iElem];
		Nodo nodosElem[cel.nodpel];
		double rMed = 0;

		for (int i = 0; i < cel.nodpel; i++) {
			nodosElem[i] = cel.getNodos()[elem[i]];
			rMed += nodosElem[i].x;
		}
		rMed /= cel.nodpel;

		switch (cel.nodpel) {
			case 3: {
				double b[3], c[3];
				Double2D pos[3];
				for (int i = 0; i < cel.nodpel; i++) {
					pos[i].x = nodosElem[i].x;
					pos[i].y = nodosElem[i].y;
				}

				double det = Armado::determinante3(pos, b, c);

				double gpvol = det * M_PI * rMed;
				for (int i = 0; i < cel.nodpel; i++) {
					cel.getMasas()[elem[i]] += gpvol;
				}

				break;
			} case 4: {
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

				for (int i = 0; i < cel.nodpel; i++) {
					double aJacob[2][2];
					for (int k = 0; k < 2; k++) for (int l = 0; l < 2; l++) aJacob[k][l] = 0;

					for (int j = 0; j < cel.nodpel; j++) {
						aJacob[0][0] += nodosElem[j].x * deriv[0][j][i];
						aJacob[0][1] += nodosElem[j].x * deriv[1][j][i];
						aJacob[1][0] += nodosElem[j].y * deriv[0][j][i];
						aJacob[1][1] += nodosElem[j].y * deriv[1][j][i];
					}

					double gpDet = aJacob[0][0] * aJacob[1][1] - aJacob[0][1] * aJacob[1][0];
					double gpVol = weigc[i] * gpDet;
					cel.getMasas()[elem[i]] += gpVol * 2 * M_PI * rMed;
				}

				break;
			}
			default: {
				assert(false);
			}
		}
	}

	for (int iNode = 0; iNode < cel.nNodes; iNode++) {
		if (cel.getMasas()[iNode] < TOLER_MASA) {
			cout << iNode << endl;
		}
		assert(cel.getMasas()[iNode] > TOLER_MASA);
	}
}

//void Transporte::carga(Celula& cel) {
//	const double CTE = 1e6 / 6.03e23;
//	const double CTE_DILUCION = CONCENTRACION_INICIAL[H_] / CONCENTRACION_INICIAL[NA];
//
//	for (int iNodo = 0; iNodo < cel.nNodes; iNodo++) {
//		cel.getCargas()[iNodo] = FARADAY / (EPSILON_TRANSPORTE * EPSILON_0) * CTE * (
//			CARGA[H_] * cel.concentraciones[H_][iNodo] +
//			CARGA[OH] * cel.concentraciones[OH][iNodo] +
//			CARGA[NA] * cel.concentraciones[NA][iNodo] * CTE_DILUCION +
//			CARGA[CL] * cel.concentraciones[CL][iNodo] * CTE_DILUCION
//		);
//	}
//}

void Transporte::concentracion(Celula& cel, int esp) {
	double esm[cel.nodpel][MAXNPEL];
	vector< Triplet<double> > triplets;
	for (int i = 0; i < cel.nNodes; i++) cel.getRhs()[i] = 0;

	for (uint kElem = 0; kElem < cel.getElementos().size(); kElem++) {
		Elemento elem = cel.getElementos()[kElem];
		Double2D pos[cel.nodpel];
		double mas[cel.nodpel];
		double sol[cel.nodpel];
		double ef[cel.nodpel];
		double ch_Med = 0.;
		double cohMed = 0.;

		for (int i = 0; i < cel.nodpel; i++) {
			int iNodo = elem[i];
			Nodo nodo = cel.getNodos()[iNodo];
			pos[i].x = nodo.x;
			pos[i].y = nodo.y;
			sol[i] = cel.anteriores[esp][iNodo];
			mas[i] = cel.getMasas()[iNodo];

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

		double landa = 1;
		double qe = 0;

		if (esp == OH) {
			qe = KWB * CONCENT_H2O - KWF * ch_Med * cohMed;
		}

		Armado::armadoTransporte(cel.nodpel, pos, sigma, qe, landa, mu, mas, sol, esm, ef);

		for (int i = 0; i < cel.nodpel; i++) {
			int iNodo = elem[i];
			Nodo nodo = cel.getNodos()[iNodo];

			if (nodo.esTierra) {
				double adiag = esm[i][i];
				for (int j = 0; j < cel.nodpel; j++) {
					esm[i][j] = 0;
					ef[j] -= esm[j][i] * CONCENTRACION_CATODO[esp];
					esm[j][i] = 0;
				}
				esm[i][i] = adiag;
				ef[i] = adiag * CONCENTRACION_CATODO[esp];
			}

			if (nodo.esPotencia) {
				double adiag = esm[i][i];
				for (int j = 0; j < cel.nodpel; j++) {
					esm[i][j] = 0;
					ef[j] -= esm[j][i] * CONCENTRACION_ANODO[esp];
					esm[j][i] = 0;
				}
				esm[i][i] = adiag;
				ef[iNodo] = adiag * CONCENTRACION_ANODO[esp];
			}
		}

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
	cel.getMatriz().makeCompressed();

//	/* SCALING */
//
//	Requiere #include "src/IterativeSolvers/Scaling.h" en eigen/unsupported/Eigen/IterativeSolvers
//	(o acá?) y #include <unsupported/Eigen/IterativeSolvers>
//
//	IterScaling<SparseMatrix<double> > scal;
//
//	scal.computeRef(cel.getMatriz());
//
//	cel.getRhs() = scal.LeftScaling().cwiseProduct(cel.getRhs());
//
//	SparseLU< SparseMatrix<double> > solver;
//
//	solver.analyzePattern(cel.getMatriz());
//	solver.factorize(cel.getMatriz());
//
//	cel.concentraciones[esp] = solver.solve(cel.getRhs());
//
//	cel.concentraciones[esp] = scal.RightScaling().cwiseProduct(cel.concentraciones[esp]);
//
//	/* END SCALING */

	/* Resolución */
//	/* LU */
//	SparseLU< SparseMatrix<double> > solver;
//
//	solver.analyzePattern(cel.getMatriz());
//	solver.factorize(cel.getMatriz());
//
//	cel.concentraciones[esp] = solver.solve(cel.getRhs());
//
//	assert(solver.info() == Success);

	BiCGSTAB< SparseMatrix<double> > solver(cel.getMatriz());
	cel.concentraciones[esp] = solver.solve(cel.getRhs());

//	if (esp == OH) {
//		cout << cel.getRhs()[17] << endl;
//	}

//	if (solver.info() != Success) {
//		for (int i = 0; i < cel.nNodes; i++) {
//			cout << cel.getRhs()[i] << endl;
//		}
//		cout << esp << "\t" << solver.error() << "\t" << solver.iterations() << endl;
//	}

	assert(solver.info() == Success);
}
