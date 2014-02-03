#include <iostream>
#include <fstream>
#include <ctime>
#include <cassert>
#include "EntradaSalida.h"
#include "Celula.h"
#include "Armado.h"

//TODO ver cosa rara en armado4!!!
//TODO transporte
//TODO multiplicar corriente para obtener por m o cm
//TODO comparar resultados de campo para tri y quad con fortran
//TODO carpeta de salida
//TODO archivo input por parámetros
//TODO fort: input y qe = 0 -> ef
//TODO ignorar bien los comentarios
//TODO refactorizar a varios archivos/clases?

//TODO PERFO: ver como graba en windows
//TODO PERFO: compilar con NDEBUG
//TODO PERFO: cambiar flags de vectorización

//TODO FORT: anteriores en transporte
//TODO FORT: armado_t algo raro en esm con transporte

Celula::Celula() {
	potencial = 0;
	nNodes = nElems = nodpel = 0;

	EntradaSalida::leerInput(*this);
}

void Celula::poisson() {
	rhs.resize(nNodes);
	rhs.fill(0.0);
	solucion.resize(nNodes);

	double error = 1.0;

	for (int contador = 0; error > EPSILON_POISSON && contador < N_COTA; contador++) {
		cout << "Armando matriz... " << flush;
		clock_t start = clock();
		vector< Triplet<double> > triplets;

		for (uint elemIdx = 0; elemIdx < elementos.size(); elemIdx++) {

			Elemento elemento = elementos[elemIdx];
			double sigma = sigmas[elemento.material];
			double ef[nodpel];

			Double2D pos[MAXNPEL];
			double esm[MAXNPEL][MAXNPEL];

			for (int i = 0; i < nodpel; i++) {
				int j = elemento[i];
				pos[i].x = nodos[j].x;
				pos[i].y = nodos[j].y;
				ef[i] = 0.0;
			}

			Armado::armadoPoisson(pos, sigma, nodpel, esm);

			/* Condiciones de contorno */
			for (int i = 0; i < nodpel; i++) {
				Nodo nodo = nodos[elemento[i]];
				double adiag = esm[i][i];

				if (nodo.esTierra) {
					for (int j = 0; j < nodpel; j++) {
						esm[i][j] = 0.0;
						ef[j] -= esm[j][i] * TIERRA;
						esm[j][i] = 0.0;
					}
					esm[i][i] = adiag;
					ef[i] = adiag * TIERRA;
				}

				if (nodo.esPotencia) {
					for (int j = 0; j < nodpel; j++) {
						esm[i][j] = 0.0;
						ef[j] -= esm[j][i] * potencial;
						esm[j][i] = 0.0;
					}
					esm[i][i] = adiag;
					ef[i] = adiag * potencial;
				}
			}

			/* Ensamblado */
			for (int i = 0; i < nodpel; i++) {
				rhs[elemento[i]] += ef[i];

				for (int j = 0; j < nodpel; j++) {
					triplets.push_back(Triplet<double>(elemento[i], elemento[j], esm[i][j]));
				}
			}
		}

		matriz.resize(nNodes, nNodes);
		matriz.setFromTriplets(triplets.begin(), triplets.end());
		int time = (clock() - start) / (CLOCKS_PER_SEC / 1000);
		cout << "OK\t\t" << time << "ms\n";

		/* Resolución */
		cout << "Resolviendo... " << flush;
		start = clock();

		SimplicialLDLT< SparseMatrix<double> > cholesky(matriz);
		solucion = cholesky.solve(rhs);

		time = (clock() - start) / (CLOCKS_PER_SEC / 1000);
		cout << "OK\t\t" << time << "ms\n";

		//TODO siempre hace una sola iteración
		error = EPSILON_POISSON * .5;
	}

	cout << "Corriente y campo... " << flush;
	clock_t start = clock();

	campo();
	corriente();

	int time = (clock() - start) / (CLOCKS_PER_SEC / 1000);
	cout << "OK\t\t" << time << "ms\n";

	grabar();
}

void Celula::campo() {
	switch (nodpel) {
	case 3:
		campo3();
		break;
	case 4:
		campo4();
		break;
	}
}

void Celula::campo3() {
	assert(nodpel == 3);
	Double2D pos[3];
	double sol[3];
	gradElem.resize(nElems);

	for (int iElem = 0; iElem < nElems; iElem++) {
		Elemento elem = elementos[iElem];
		double b[3], c[3];

		for (int i = 0; i < 3; i++) {
			int iNodo = elem[i];
			pos[i].x = nodos[iNodo].x;
			pos[i].y = nodos[iNodo].y;
			sol[i] = solucion[iNodo];
		}

		double det = Armado::determinante3(pos, b, c);

		gradElem[iElem].x = (b[0] * sol[0] + b[1] * sol[1] + b[2] * sol[2]) / det;
		gradElem[iElem].y = (c[0] * sol[0] + c[1] * sol[1] + c[2] * sol[2]) / det;
	}
}

void Celula::campo4() {
	assert(nodpel == 4);
	gradElem.resize(nElems);

	for (int iElem = 0; iElem < nElems; iElem++) {
		Elemento elem = elementos[iElem];
		int kGauss = 0;
		double phi[2*NGAUSS][4];
		double dphi[NDIM][2*NGAUSS][4];
		double phidX[NDIM][2*NGAUSS][4];
		double sol[nodpel];
		Double2D pos[4];
		Double2D e[nodpel];
		Double2D eElem;
		eElem.x = 0.;
		eElem.y = 0.;

		for (int i = 0; i < 3; i++) {
			int iNodo = elem[i];
			pos[i].x = nodos[iNodo].x;
			pos[i].y = nodos[iNodo].y;
			sol[i] = solucion[i];
			e[i].x = 0.;
			e[i].y = 0.;
		}

		for (int i = 0; i < NGAUSS; i++) for (int j = 0; j < NGAUSS; j++) {
			Armado::iteracion4(i, j, kGauss, pos, phi, dphi, phidX);

			for (int k = 0; k < nodpel; k++) {
				e[k].x  += dphi[0][kGauss][k] * sol[k];
				e[k].y  += dphi[1][kGauss][k] * sol[k];
				eElem.x += dphi[0][kGauss][k] * sol[k];
				eElem.y += dphi[1][kGauss][k] * sol[k];
			}

			kGauss++;
		}

		gradElem[iElem].x = -eElem.x * 0.25;
		gradElem[iElem].y = -eElem.y * 0.25;
	}
}

//TODO multiplicar para que quede por metro?
void Celula::corriente() {
	corrElem.resize(nElems);

	for (int iElem = 0; iElem < nElems; iElem++) {
		Elemento elem = elementos[iElem];
		corrElem[iElem].x = -sigmas[elem.material] * gradElem[iElem].x;
		gradElem[iElem].y = -sigmas[elem.material] * gradElem[iElem].y;
	}
}

void Celula::grabar() {
	cout << "Grabando... " << flush;
	clock_t start = clock();

	/* Tensión */
	ofstream tension("tension.csv", ofstream::out);

	tension << "X, Y, V";

	for (int iNodo = 0; iNodo < nNodes; iNodo++) {
		Nodo nodo = nodos[iNodo];
		tension << "\n" << nodo.x << ", " << nodo.y << ", " << solucion[iNodo];
	}

	tension.close();

	/* Corriente y campo */
	if (nodpel == 3) {
		ofstream corriente("corriente.csv", ofstream::out);
		ofstream campo("campo.csv", ofstream::out);

		corriente << "X, Y, corriente";
		campo 	  << "X, Y, campo";

		for (int k = 0; k < nElems; k++) {
			double corr = sqrt(pow( gradElem[k].x, 2) + pow( gradElem[k].y, 2));
			double camp = sqrt(pow(corrElem[k].x, 2) + pow(corrElem[k].y, 2));
			double xMed = 0.0, yMed = 0.0;

			for (int j = 0; j < nodpel; j++) {
				int jNodo = elementos[k][j];
				xMed += nodos[jNodo].x;
				yMed += nodos[jNodo].y;
			}

			xMed /= nodpel;
			yMed /= nodpel;

			corriente << "\n" << xMed << ", " << yMed << ", " << corr;
			campo 	  << "\n" << xMed << ", " << yMed << ", " << camp;
		}

		corriente.close();
		campo.close();
	}

	int time = (clock() - start) / (CLOCKS_PER_SEC / 1000);
	cout << "OK\t\t\t" << time << "ms\n";
}

Elemento Celula::getElement(int i) {
	return elementos[i];
}

void Celula::chequearSimetria() {
	for (int i = 0; i < nNodes; i++) for (int j = i; j < nNodes; j++) {
		assert(abs(matriz.coeff(i, j) - matriz.coeff(j, i)) < 1e-9);
	}
}

void Celula::transporte() {
	const double T_CERO = 1.;

	for (int esp = 0; esp < NESPS; esp++) {
		concentraciones[esp].resize(nNodes);
//		anteriores[esp].resize(nNodes);
	}

	phAux[H_].resize(nNodes);
	phAux[OH].resize(nNodes);

	for (int iNode = 0; iNode < nNodes; iNode++) {
		for (int esp = 0; esp < NESPS; esp++) {
			concentraciones[esp][iNode] = CONCENTRACION_INICIAL[esp];
//			anteriores[esp][iNode] = CONCENTRACION_INICIAL[esp];
		}

//		phAux[H_][iNode] = -log10(concentraciones[H_][iNode] * 1e15 / 6.02e23);
//		phAux[OH][iNode] = -log10(concentraciones[OH][iNode] * 1e15 / 6.02e23);
	}

	masaDiag2D();
	cargas.resize(nNodes);

	for (double tt = 0.; tt < T_CERO; tt += DELTA_T) {
		poisson();

		for (int esp = 0; esp < NESPS; esp++) {
			concentracion(esp);
		}

		for (int kNodo = 0; kNodo < nNodes; kNodo++) {
			phAux[H_][kNodo] = -log10(concentraciones[H_][kNodo] * 1e15 / 6.02e23);
			phAux[OH][kNodo] = -log10(concentraciones[OH][kNodo] * 1e15 / 6.02e23);

			for (int esp = 0; esp < NESPS; esp++) {
				concentraciones[esp][kNodo] = RSA * concentraciones[esp][kNodo] + (1-RSA) * anteriores[esp][kNodo];
				if (concentraciones[esp][kNodo] < CONCENT_MINIMO) {
					concentraciones[esp][kNodo] = 0.;
				}
				anteriores[esp][kNodo] = concentraciones[esp][kNodo];
			}
		}

//		TODO escribir por salida
	}
}

void Celula::masaDiag2D() {
	assert(nodpel == 4);

	masas.resize(nNodes);
	for (int i = 0; i < nNodes; i++) masas[i] = 0.0;

	const int NLOCS = 2;
	const int INOGA[] = {0, 3, 1, 2};
	const double POSGL[] = {-1.0, 1.0};
	const double WEIGL[] = { 1.0, 1.0};

	for (int iElem = 0; iElem < nElems; iElem++) {
		Elemento elem = elementos[iElem];
		Nodo nodosElem[nodpel];
		for (int i = 0; i < nodpel; i++) nodosElem[i] = nodos[elem[i]];

		/*	subroutine armotodo(nope,x,y,deriv,weigc) */
		int iGauss = 0;
		double weigc[nodpel];
		double posgc[2][nodpel];
		double deriv[2][nodpel][nodpel];

		for (int ilocs = 0; ilocs < NLOCS; ilocs++) {
			for (int jlocs = 0; jlocs < NLOCS; jlocs++) {
				weigc[INOGA[iGauss]] = WEIGL[ilocs] * WEIGL[jlocs];
				posgc[0][INOGA[iGauss]] = POSGL[ilocs];
				posgc[1][INOGA[iGauss]] = POSGL[jlocs];
				iGauss++;
			}
		}

		for (int i = 0; i < nodpel; i++) {
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

		for (int iNode = 0; iNode < nodpel; iNode++) {
			double aJacob[2][2];
			for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++) aJacob[i][j] = 0.0;

			for (int jNode = 0; jNode < nodpel; jNode++) {
				aJacob[0][0] += nodosElem[jNode].x * deriv[0][jNode][iNode];
				aJacob[0][1] += nodosElem[jNode].x * deriv[1][jNode][iNode];
				aJacob[1][0] += nodosElem[jNode].y * deriv[0][jNode][iNode];
				aJacob[1][1] += nodosElem[jNode].y * deriv[1][jNode][iNode];
			}

			double gpDet = aJacob[0][0] * aJacob[1][1] - aJacob[0][1] * aJacob[1][0];
			double gpVol = weigc[iNode] * gpDet;
			masas[elem[iNode]] += gpVol;
		}

		for (int iPoint = 0; iPoint < nNodes; iPoint++)	assert(masas[iPoint] > TOLER_MASA);
	}
}

void Celula::carga() {
	const double CTE = 1e6 / 6.03e23;
	const double CTE_DILUCION = CONCENTRACION_INICIAL[H_] / CONCENTRACION_INICIAL[NA];

	for (int iNodo = 0; iNodo < nNodes; iNodo++) {
		cargas[iNodo] =	FARADAY / (EPSILON_TRANSPORTE * EPSILON_0) * CTE * (
			CARGA[H_] * concentraciones[H_][iNodo] +
			CARGA[OH] * concentraciones[OH][iNodo] +
			CARGA[NA] * concentraciones[NA][iNodo] * CTE_DILUCION +
			CARGA[CL] * concentraciones[CL][iNodo] * CTE_DILUCION
		);
	}
}

void Celula::concentracion(int esp) {
	double esm[nodpel][MAXNPEL];
	vector< Triplet<double> > triplets;

	for (uint kElem = 0; kElem < elementos.size(); kElem++) {
		Elemento elem = elementos[kElem];
		Double2D pos[nodpel];
		double sol[nodpel];
		double ef[nodpel];
		double ch_Med = 0.;
		double cohMed = 0.;

		for (int i = 0; i < nodpel; i++) {
			int iNodo = elem[i];
			Nodo nodo = nodos[iNodo];
			pos[i].x = nodo.x;
			pos[i].y = nodo.y;
			sol[i] = solucion[iNodo];

			if (esp == OH) {
				ch_Med += anteriores[H_][iNodo];
				cohMed += anteriores[OH][iNodo];
			}
		}

		if (esp == OH) {
			ch_Med /= nodpel;
			cohMed /= nodpel;
		}

		double sigma = sigmas[elem.material];

		double difElem = DIFUSION[esp];
		if (elem.material == MEMBRANA) difElem *= 1e-3;

		double mu = -difElem * CLAVE * CARGA[esp] * gradElem[kElem].y;
		double landa = 1.;
		double qe = 0.;

		if (esp == OH) {
			qe = KWB * CONCENT_H2O - KWF * ch_Med * cohMed;
		}

		Armado::armadoTransporte(pos, sigma, qe, landa, mu, sol, esm, ef);

		for (int i = 0; i < nodpel; i++) {
			int iNodo = elem[i];
			Nodo nodo = nodos[iNodo];

			if (nodo.esTierra) {
				double adiag = esm[i][i];
				for (int j = 0; j < nodpel; j++) {
					esm[i][j] = 0.;
					ef[j] -= esm[j][i] * CONCENTRACION_CATODO[esp];
					esm[j][i] = 0.;
				}
				esm[i][i] = adiag;
				ef[i] = adiag * CONCENTRACION_CATODO[esp];
			}

			if (nodo.esPotencia) {
				double adiag = esm[i][i];
				for (int j = 0; j < nodpel; j++) {
					esm[i][j] = 0.;
					ef[j] -= esm[j][i] * CONCENTRACION_ANODO[esp];
					esm[j][i] = 0.;
				}
				esm[i][i] = adiag;
				ef[iNodo] = adiag * CONCENTRACION_ANODO[esp];
			}
		}

		/* Ensamblado */
		for (int i = 0; i < nNodes; i++) rhs[i] = 0.;
		triplets.clear();

		for (int i = 0; i < nodpel; i++) {
			int iNodo = elem[i];
			rhs[iNodo] += ef[i];

			for (int j = 0; j < nodpel; j++) {
				int jNodo = elem[j];
				triplets.push_back(Triplet<double>(iNodo, jNodo, esm[i][j]));
			}
		}
	}

	matriz.resize(nNodes, nNodes);
	matriz.setFromTriplets(triplets.begin(), triplets.end());

	/* Resolución */
//	TODO usar LU
	BiCGSTAB< SparseMatrix<double> > solver;

//	TODO usar precondicionador
	solver.compute(matriz);

//	concentraciones[esp] = solver.solve(rhs);
	VectorXd guess = concentraciones[esp];
	concentraciones[esp] = solver.solveWithGuess(rhs, guess);
}

vector<Nodo>& Celula::getNodos() {
	return nodos;
}

vector<Elemento>& Celula::getElementos() {
	return elementos;
}

vector<Double2D>& Celula::getGradElem() {
	return gradElem;
}

vector<Double2D>& Celula::getCorrElem();

SparseMatrix<double>& getMatriz();

VectorXd& getRhs();
VectorXd& getSolucion();
