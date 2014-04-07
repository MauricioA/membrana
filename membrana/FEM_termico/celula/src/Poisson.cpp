#include <iostream>
#include "Poisson.h"
#include "Armado.h"
#include "EntradaSalida.h"

using namespace std;

//void Poisson::poisson(Celula &celula) {
//	celula.getRhs().resize(celula.nNodes);
//	celula.getRhs().fill(0);
//	celula.getSolucion().resize(celula.nNodes);
//	vector<Triplet<double>> global_triplets(celula.nElems * celula.nodpel * celula.nodpel);
//
//	#pragma omp parallel num_threads(2)
//	{
//		VectorXd local_rhs;
//		local_rhs.resize(celula.nNodes);
//		local_rhs.fill(0);
//		
//		#pragma omp for
//		for (int elemIdx = 0; elemIdx < celula.nElems; elemIdx++) {
//			Elemento elemento = celula.getElementos()[elemIdx];
//			double sigma = elemento.sigma;
//			double ef[MAXNPEL];
//
//			Double2D pos[MAXNPEL];
//			double esm[MAXNPEL][MAXNPEL];
//
//			for (int i = 0; i < celula.nodpel; i++) {
//				int j = elemento[i];
//				pos[i].x = celula.getNodos()[j].x;
//				pos[i].y = celula.getNodos()[j].y;
//				ef[i] = 0.0;
//			}
//
//			Armado::armadoPoisson(pos, sigma, celula.nodpel, esm);
//
//			/* Condiciones de contorno */
//			for (int i = 0; i < celula.nodpel; i++) {
//				Nodo nodo = celula.getNodos()[elemento[i]];
//				double adiag = esm[i][i];
//
//				if (nodo.esTierra) {
//					for (int j = 0; j < celula.nodpel; j++) {
//						esm[i][j] = 0.0;
//						ef[j] -= esm[j][i] * TIERRA;
//						esm[j][i] = 0.0;
//					}
//					esm[i][i] = adiag;
//					ef[i] = adiag * TIERRA;
//				}
//
//				if (nodo.esPotencia) {
//					for (int j = 0; j < celula.nodpel; j++) {
//						esm[i][j] = 0.0;
//						ef[j] -= esm[j][i] * celula.potencial;
//						esm[j][i] = 0.0;
//					}
//					esm[i][i] = adiag;
//					ef[i] = adiag * celula.potencial;
//				}
//			}
//
//			/* Ensamblado */
//			for (int i = 0; i < celula.nodpel; i++) {
//				//celula.getRhs()[elemento[i]] += ef[i];
//				local_rhs[elemento[i]] += ef[i];
//
//				for (int j = 0; j < celula.nodpel; j++) {
//					//global_triplets.push_back(Triplet<double>(elemento[i], elemento[j], esm[i][j]));
//
//					global_triplets[elemIdx * celula.nodpel * celula.nodpel + i * celula.nodpel + j] 
//						= Triplet<double>(elemento[i], elemento[j], esm[i][j]);
//				}
//			}
//		}
//
//		#pragma omp critical(reduce_rhs)
//		for (int i = 0; i < celula.nNodes; i++) {
//			celula.getRhs()[i] += local_rhs[i];
//		}
//	}
//
//	celula.getMatriz().resize(celula.nNodes, celula.nNodes);
//	celula.getMatriz().setFromTriplets(global_triplets.begin(), global_triplets.end());
//
//	/* Resolución */
//	SimplicialLDLT<SparseMatrix<double>> cholesky(celula.getMatriz());
//	celula.setSolucion(cholesky.solve(celula.getRhs()));
//
//	//ConjugateGradient<SparseMatrix<double>> solver(celula.getMatriz());
//	//solver.setTolerance(1e-10);
//	//celula.setSolucion(solver.solveWithGuess(celula.getRhs(), celula.getSolucion()));
//	//assert(solver.info() == Success);
//
//	campo(celula);
//	corriente(celula);
//}

void Poisson::poisson(Celula &celula) {
	celula.getRhs().resize(celula.nNodes);
	celula.getRhs().fill(0);
	celula.getSolucion().resize(celula.nNodes);
	vector<Triplet<double>> triplets(celula.nElems * celula.nodpel * celula.nodpel);
		
	for (int elemIdx = 0; elemIdx < celula.nElems; elemIdx++) {
		Elemento elemento = celula.getElementos()[elemIdx];
		double sigma = elemento.sigma;
		double ef[MAXNPEL];

		Double2D pos[MAXNPEL];
		double esm[MAXNPEL][MAXNPEL];

		for (int i = 0; i < celula.nodpel; i++) {
			int j = elemento[i];
			pos[i].x = celula.getNodos()[j].x;
			pos[i].y = celula.getNodos()[j].y;
			ef[i] = 0.0;
		}

		Armado::armadoPoisson(pos, sigma, celula.nodpel, esm);

		/* Condiciones de contorno */
		for (int i = 0; i < celula.nodpel; i++) {
			Nodo nodo = celula.getNodos()[elemento[i]];
			double adiag = esm[i][i];

			if (nodo.esTierra) {
				for (int j = 0; j < celula.nodpel; j++) {
					esm[i][j] = 0.0;
					ef[j] -= esm[j][i] * TIERRA;
					esm[j][i] = 0.0;
				}
				esm[i][i] = adiag;
				ef[i] = adiag * TIERRA;
			}

			if (nodo.esPotencia) {
				for (int j = 0; j < celula.nodpel; j++) {
					esm[i][j] = 0.0;
					ef[j] -= esm[j][i] * celula.potencial;
					esm[j][i] = 0.0;
				}
				esm[i][i] = adiag;
				ef[i] = adiag * celula.potencial;
			}
		}

		/* Ensamblado */
		for (int i = 0; i < celula.nodpel; i++) {
			celula.getRhs()[elemento[i]] += ef[i];

			for (int j = 0; j < celula.nodpel; j++) {
				triplets.push_back(Triplet<double>(elemento[i], elemento[j], esm[i][j]));
			}
		}
	}

	celula.getMatriz().resize(celula.nNodes, celula.nNodes);
	celula.getMatriz().setFromTriplets(triplets.begin(), triplets.end());

	/* Resolución */
	SimplicialLDLT<SparseMatrix<double>> cholesky(celula.getMatriz());
	celula.setSolucion(cholesky.solve(celula.getRhs()));

	//ConjugateGradient<SparseMatrix<double>> solver(celula.getMatriz());
	//solver.setTolerance(1e-10);
	//celula.setSolucion(solver.solveWithGuess(celula.getRhs(), celula.getSolucion()));
	//assert(solver.info() == Success);

	campo(celula);
	corriente(celula);
}

void Poisson::campo(Celula &celula) {
	switch (celula.nodpel) {
	case 3:
		campo3(celula);
		break;
	case 4:
		campo4(celula);
		break;
	default: 
		assert(false);
	}
}

void Poisson::campo3(Celula &celula) {
	assert(celula.nodpel == 3);
	Double2D pos[3];
	double sol[3];
	celula.getGradElem().resize(celula.nElems);

	for (int iElem = 0; iElem < celula.nElems; iElem++) {
		Elemento elem = celula.getElementos()[iElem];
		double b[3], c[3];

		for (int i = 0; i < 3; i++) {
			int iNodo = elem[i];
			pos[i].x = celula.getNodos()[iNodo].x;
			pos[i].y = celula.getNodos()[iNodo].y;
			sol[i] = celula.getSolucion()[iNodo];
		}

		double det = Armado::determinante3(pos, b, c);

		celula.getGradElem()[iElem].x = (b[0] * sol[0] + b[1] * sol[1] + b[2] * sol[2]) / det;
		celula.getGradElem()[iElem].y = (c[0] * sol[0] + c[1] * sol[1] + c[2] * sol[2]) / det;
	}
}

void Poisson::campo4(Celula &celula) {
	assert(celula.nodpel == 4);
	const int NODPEL = 4;
	celula.getGradElem().resize(celula.nElems);

	for (int iElem = 0; iElem < celula.nElems; iElem++) {
		Elemento elem = celula.getElementos()[iElem];
		int kGauss = 0;
		double phi[2*NGAUSS][4];
		double dphi[NDIM][2*NGAUSS][4];
		double phidX[NDIM][2*NGAUSS][4];
		double sol[NODPEL];
		Double2D pos[4];
		Double2D eElem;
		eElem.x = 0;
		eElem.y = 0;

		for (int i = 0; i < NODPEL; i++) {
			int iNodo = elem[i];
			pos[i].x = celula.getNodos()[iNodo].x;
			pos[i].y = celula.getNodos()[iNodo].y;
			sol[i] = celula.getSolucion()[i];
		}

		for (int i = 0; i < NGAUSS; i++) for (int j = 0; j < NGAUSS; j++) {
			Armado::iteracion4(i, j, kGauss, pos, phi, dphi, phidX);

			for (int k = 0; k < NODPEL; k++) {
				eElem.x += dphi[0][kGauss][k] * sol[k];
				eElem.y += dphi[1][kGauss][k] * sol[k];
			}

			kGauss++;
		}

		celula.getGradElem()[iElem].x = -eElem.x * 0.25;
		celula.getGradElem()[iElem].y = -eElem.y * 0.25;
	}
}

void Poisson::corriente(Celula &celula) {
	celula.getCorrElem().resize(celula.nElems);

	for (int iElem = 0; iElem < celula.nElems; iElem++) {
		Elemento elem = celula.getElementos()[iElem];
		celula.getCorrElem()[iElem].x = -elem.sigma * celula.getGradElem()[iElem].x;
		celula.getCorrElem()[iElem].y = -elem.sigma * celula.getGradElem()[iElem].y;
	}
}
