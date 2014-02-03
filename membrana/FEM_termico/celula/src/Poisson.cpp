#include <iostream>
#include "Poisson.h"
#include "Armado.h"

using namespace std;

//TODO multiplicar corriente para obtener por m o cm
//TODO comparar resultados de campo para tri y quad con fortran

void Poisson::poisson(Celula &celula) {
	celula.getRhs().resize(celula.nNodes);
	celula.getRhs().fill(0.0);
	celula.getSolucion().resize(celula.nNodes);

	double error = 1.0;

	for (int contador = 0; error > EPSILON_POISSON && contador < N_COTA; contador++) {
		cout << "Armando matriz... " << flush;
		clock_t start = clock();
		vector< Triplet<double> > triplets;

		for (int elemIdx = 0; elemIdx < celula.nElems; elemIdx++) {

			Elemento elemento = celula.getElementos()[elemIdx];
			double sigma = celula.sigmas[elemento.material];
			double ef[celula.nodpel];

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
		int time = (clock() - start) / (CLOCKS_PER_SEC / 1000);
		cout << "OK\t\t" << time << "ms\n";

		/* Resolución */
		cout << "Resolviendo... " << flush;
		start = clock();

		SimplicialLDLT< SparseMatrix<double> > cholesky(celula.getMatriz());
		celula.setSolucion(cholesky.solve(celula.getRhs()));

		time = (clock() - start) / (CLOCKS_PER_SEC / 1000);
		cout << "OK\t\t" << time << "ms\n";

		//TODO siempre hace una sola iteración
		error = EPSILON_POISSON * .5;
	}

	cout << "Corriente y campo... " << flush;
	clock_t start = clock();

	campo(celula);
	corriente(celula);

	int time = (clock() - start) / (CLOCKS_PER_SEC / 1000);
	cout << "OK\t\t" << time << "ms\n";
}

void Poisson::campo(Celula &celula) {
	switch (celula.nodpel) {
	case 3:
		campo3(celula);
		break;
	case 4:
		campo4(celula);
		break;
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
		Double2D e[NODPEL];
		Double2D eElem;
		eElem.x = 0.;
		eElem.y = 0.;

		for (int i = 0; i < 3; i++) {
			int iNodo = elem[i];
			pos[i].x = celula.getNodos()[iNodo].x;
			pos[i].y = celula.getNodos()[iNodo].y;
			sol[i] = celula.getSolucion()[i];
			e[i].x = 0.;
			e[i].y = 0.;
		}

		for (int i = 0; i < NGAUSS; i++) for (int j = 0; j < NGAUSS; j++) {
			Armado::iteracion4(i, j, kGauss, pos, phi, dphi, phidX);

			for (int k = 0; k < NODPEL; k++) {
				e[k].x  += dphi[0][kGauss][k] * sol[k];
				e[k].y  += dphi[1][kGauss][k] * sol[k];
				eElem.x += dphi[0][kGauss][k] * sol[k];
				eElem.y += dphi[1][kGauss][k] * sol[k];
			}

			kGauss++;
		}

		celula.getGradElem()[iElem].x = -eElem.x * 0.25;
		celula.getGradElem()[iElem].y = -eElem.y * 0.25;
	}
}

//TODO multiplicar para que quede por metro?
void Poisson::corriente(Celula &celula) {
	celula.getCorrElem().resize(celula.nElems);

	for (int iElem = 0; iElem < celula.nElems; iElem++) {
		Elemento elem = celula.getElementos()[iElem];
		celula.getCorrElem()[iElem].x = -celula.sigmas[elem.material] * celula.getGradElem()[iElem].x;
		celula.getCorrElem()[iElem].y = -celula.sigmas[elem.material] * celula.getGradElem()[iElem].y;
	}
}
