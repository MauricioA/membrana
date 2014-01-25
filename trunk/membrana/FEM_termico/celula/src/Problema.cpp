#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>

#include "Problema.h"

//TODO archivo input por parámetros
//TODO valgrind!
//TODO fort: input y qe = 0 -> ef

Problema::Problema() {
	string s, line, malla;
	vector<int> dirichV, dirichT;

	cout << "Leyendo archivos... ";
	ifstream input("input.in", ifstream::in);

	if (!input.is_open()) {
		cerr << "Error abriendo archivo input\n";
		exit(EXIT_FAILURE);
	}

	while (getline(input, line)) {
		istringstream iss(line);

		if (line.find("archivo1") != string::npos) {
			iss >> s >> malla;
		} else if (line.find("sigint") != string::npos) {
			iss >> s >> sigmas[INTERNO];
		} else if (line.find("sigext") != string::npos) {
			iss >> s >> sigmas[EXTERNO];
		} else if (line.find("sigmem") != string::npos) {
			iss >> s >> sigmas[MEMBRANA];
		} else if (line.find("Potencial") != string::npos) {
			iss >> s >> potencial;
		} else if (line.find("dirichV") != string::npos) {
			int nDirichV;
			iss >> s >> nDirichV;

			for (int i = 0; i < nDirichV; i++) {
				getline(input, line);
				istringstream iss2(line);
				int nodo;
				iss2 >> nodo;
				dirichV.push_back(nodo - 1);
			}
		} else if (line.find("dirichT") != string::npos) {
			int nDirichT;
			iss >> s >> nDirichT;

			for (int i = 0; i < nDirichT; i++) {
				getline(input, line);
				istringstream iss2(line);
				int nodo;
				iss2 >> nodo;
				dirichT.push_back(nodo - 1);
			}
		}
	}

	input.close();

	leerMalla(malla);

	for (uint i = 0; i < dirichV.size(); i++) {
		nodos[dirichV[i]].esPotencia = true;
	}

	for (uint i = 0; i < dirichT.size(); i++) {
		nodos[dirichT[i]].esTierra = true;
	}

	cout << "OK\n";
}

void Problema::leerMalla(string malla) {
	int n, nodo1, nodo2, nodo3;
	double x, y;
	int max = 256;
	char line[max];

	ifstream stream(malla.c_str(), ifstream::in);
	if (!stream.is_open()) {
		cerr << "Error abriendo malla\n";
		exit(EXIT_FAILURE);
	}

	/* *COORDINATES */
	stream.getline(line, max);

	/* Numero de nodos */
	stream >> nNodes;
	nodos.reserve(nNodes);

	/* Nodos */
	for (int i = 0; i < nNodes; i++) {
		stream >> n >> x >> y;
		Nodo nodo;
		nodo.x = x;
		nodo.y = y;
		nodo.esPotencia = false;
		nodo.esTierra = false;
		nodos.push_back(nodo);
	}

	/* *ELEMENT_GROUPS - 3 */
	stream.getline(line, max);
	stream.getline(line, max);
	stream.getline(line, max);

	/* Grupos */
	int elemsExt, elemsMemb, elemsInt;
	stream >> n >> elemsExt;
	stream.getline(line, max);
	stream >> n >> elemsMemb;
	stream.getline(line, max);
	stream >> n >> elemsInt;
	stream.getline(line, max);
	nElems = elemsExt + elemsMemb + elemsInt;
	elementos.reserve(nElems);

	/* *INCIDENCES */
	stream.getline(line, max);

	/* Elementos extrenos */
	for (int i = 0; i < elemsExt; i++) {
		stream >> n >> nodo1 >> nodo2 >> nodo3;
		int nodos[3] = {nodo1-1, nodo2-1, nodo3-1};
		elementos.push_back(Elemento(nodos, EXTERNO));
	}

	/* Elementos membrana */
	for (int i = 0; i < elemsMemb; i++) {
		stream >> n >> nodo1 >> nodo2 >> nodo3;
		int nodos[3] = {nodo1-1, nodo2-1, nodo3-1};
		elementos.push_back(Elemento(nodos, MEMBRANA));
	}

	/* Elementos internos */
	for (int i = 0; i < elemsInt; i++) {
		stream >> n >> nodo1 >> nodo2 >> nodo3;
		int nodos[3] = {nodo1-1, nodo2-1, nodo3-1};
		elementos.push_back(Elemento(nodos, INTERNO));
	}

	stream.close();
}

void Problema::poisson() {
	rhs.resize(nNodes);
	rhs.fill(0.0);
	solucion.resize(nNodes);
	matriz.resize(nNodes, nNodes);

	double error = 1.0;

	for (int contador = 0; error > EPSILON && contador < N_COTA; contador++) {
		vector< Triplet<double> > triplets;

		for (uint elemIdx = 0; elemIdx < elementos.size(); elemIdx++) {

			Elemento elemento = elementos[elemIdx];
			double sigma = sigmas[elemento.material];
			double ef[NODPEL];

			double x[NODPEL], y[NODPEL];
			double esm[3][3];

			for (int i = 0; i < NODPEL; i++) {
				int j = elemento[i];
				x[i] = nodos[j].x;
				y[i] = nodos[j].y;
				ef[i] = 0.0;
			}

			armado3(x, y, esm, sigma);

			/* Condiciones de contorno */
			for (int i = 0; i < NODPEL; i++) {
				Nodo nodo = nodos[elemento[i]];
				double adiag = esm[i][i];

				if (nodo.esTierra) {
					for (int j = 0; j < NODPEL; j++) {
						esm[i][j] = 0.0;
						ef[j] -= esm[j][i] * TIERRA;
						esm[j][i] = 0.0;
					}
					esm[i][i] = adiag;
					ef[i] = adiag * TIERRA;
				}

				if (nodo.esPotencia) {
					for (int j = 0; j < NODPEL; j++) {
						esm[i][j] = 0.0;
						ef[j] -= esm[j][i] * potencial;
						esm[j][i] = 0.0;
					}
					esm[i][i] = adiag;
					ef[i] = adiag * potencial;
				}
			}

			/* Ensamblado */
			for (int i = 0; i < NODPEL; i++) {
				rhs[elemento[i]] += ef[i];

				for (int j = 0; j < NODPEL; j++) {
					triplets.push_back(Triplet<double>(elemento[i], elemento[j], esm[i][j]));
				}
			}
		}

		matriz.setFromTriplets(triplets.begin(), triplets.end());

		/* Resolución */
		ConjugateGradient< SparseMatrix<double> > cg(matriz);

		cout << "Resolviendo... ";

		solucion = cg.solve(rhs);

		cout << "OK\terror: " << (double) cg.error() << " iters: " << cg.iterations() << "\n";

		//TODO siempre hace una sola iteración
		error = EPSILON * .5;
	}

	cout << "Corriente y campo... ";

	corriente();
	campo();

	cout << "OK\n";
	cout << "Grabando... ";

	grabar();

	cout << "OK\n";
}

void Problema::armado3(double x[], double y[], double esm[3][3], double sigma) {
	double b[3], c[3];

	double det = determinante3(x, y, b, c);

	double rMed = (x[0] + x[1] + x[2]) / 3;

	if (abs(det) < TOLER_AREA) {
		cerr << "Error area es cero\n";
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) {
		esm[i][j] = sigma * (b[i] * b[j] + c[i] * c[j]) * M_PI * rMed / det;
	}
}

void Problema::armado4(double x[], double y[], double esm[4][4], double sigma) {
	const int NGAUSS = 2, NDIM = 2;
	const double GAUSSPT[] = { (- 1 / M_SQRT3), (1 / M_SQRT3) };
	const double GAUSSWT[] = { 1.0, 1.0 };

	int kGauss = 0;
	for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) esm[i][j] = 0.0;

	double phi[2*NGAUSS][NODPEL];
	double dphi[NDIM][2*NGAUSS][NODPEL];
	double gxCod[NDIM][2*NGAUSS];
	double phidX[NDIM][2*NGAUSS][NODPEL];
	double cteI[2*NGAUSS];

	for (int j = 0; j < NGAUSS; j++) for (int i = 0; i < NGAUSS; i++) {
		double t = GAUSSPT[j];
		double s = GAUSSWT[i];

		double sm = 0.5 * (1.0 - s);
		double tm = 0.5 * (1.0 - t);
		double sq = 0.5 * (1.0 + s);
		double tp = 0.5 * (1.0 + t);

		int k = 0;
		phi[kGauss][k++] = sm * tm;
		phi[kGauss][k++] = sq * tm;
		phi[kGauss][k++] = sq * tp;
		phi[kGauss][k++] = sm * tp;

		k = 0;
		dphi[0][kGauss][k++] = -0.5 * tm;
		dphi[0][kGauss][k++] = -0.5 * tm;
		dphi[0][kGauss][k++] = -0.5 * tp;
		dphi[0][kGauss][k++] = -0.5 * tp;

		k = 0;
		dphi[0][kGauss][k++] = -0.5 * sm;
		dphi[0][kGauss][k++] = -0.5 * sq;
		dphi[0][kGauss][k++] = -0.5 * sq;
		dphi[0][kGauss][k++] = -0.5 * sm;

		for (int dim = 0; dim < NDIM; dim++) {
			gxCod[dim][kGauss] = 0.0;

			for (int i = 0; i < NODPEL; i++) {
				gxCod[dim][kGauss] += x[i] * phi[kGauss][i];
			}
		}

		double aJacob2D[2][2];
		for (k = 0; k < 2; k++) for (int l = 0; l < 2; l++) aJacob2D[k][l] = 0.0;

		for (k = 0; k < NODPEL; k++) {
			aJacob2D[0][0] += dphi[0][kGauss][k] * x[k];
			aJacob2D[0][1] += dphi[0][kGauss][k] * y[k];
			aJacob2D[1][0] += dphi[1][kGauss][k] * x[k];
			aJacob2D[1][1] += dphi[1][kGauss][k] * y[k];
		}

		double det =
			aJacob2D[0][0] * aJacob2D[1][1] -
			aJacob2D[0][1] * aJacob2D[1][0];

		double aJacobI2d[2][2] = {
			{
				aJacob2D[1][1] / det,
				aJacob2D[0][1] / det,
			}, {
				aJacob2D[1][0] / det,
				aJacob2D[0][0] / det,
			},
		};

		for (int k = 0; k < NODPEL; k++) {
			phidX[0][kGauss][k] =
				aJacobI2d[0][0] * dphi[0][kGauss][k] +
				aJacobI2d[0][1] * dphi[1][kGauss][k];

			phidX[1][kGauss][k] =
				aJacobI2d[1][0] * dphi[0][kGauss][k] +
				aJacobI2d[1][1] * dphi[1][kGauss][k];
		}

		cteI[kGauss] = det * GAUSSWT[i] * GAUSSWT[j] * 2 * M_PI * gxCod[0][kGauss];
		kGauss++;
	}

	for (kGauss = 0; kGauss < NGAUSS*NGAUSS; kGauss++) for (int i = 0; i < NODPEL; i++) for (int j = 0; j < NODPEL; j++) for (int dim = 0; dim < NDIM; dim++) {
		esm[i][j] += sigma * phidX[dim][kGauss][i] * phidX[dim][kGauss][j] * cteI[kGauss];
	}
}

void Problema::corriente() {
	double x[3], y[3], sol[3];
	corrElemX.resize(nElems);
	corrElemY.resize(nElems);

	for (int iElem = 0; iElem < nElems; iElem++) {
		Elemento elem = elementos[iElem];
		double b[3], c[3];

		for (int i = 0; i < NODPEL; i++) {
			int iNodo = elem[i];
			x[i] = nodos[iNodo].x;
			y[i] = nodos[iNodo].y;
			sol[i] = solucion[iNodo];
		}

		double det = determinante3(x, y, b, c);

		corrElemX[iElem] = (b[0] * sol[0] + b[1] * sol[1] + b[2] * sol[2]) / det;
		corrElemY[iElem] = (c[0] * sol[0] + c[1] * sol[1] + c[2] * sol[2]) / det;
	}
}

double Problema::determinante3(double x[], double y[], double b[], double c[]) {
	int i = 0;
	b[i++] = y[1] - y[2];
	b[i++] = y[2] - y[0];
	b[i++] = y[0] - y[1];

	i = 0;
	c[i++] = x[2] - x[1];
	c[i++] = x[0] - x[2];
	c[i++] = x[1] - x[0];

	return
		+ x[1]*y[2] + x[2]*y[0] + x[0]*y[1]
		- x[1]*y[0] - x[2]*y[1] - x[0]*y[2];
}

void Problema::campo() {
	campoElemX.resize(nElems);
	campoElemY.resize(nElems);

	for (int iElem = 0; iElem < nElems; iElem++) {
		Elemento elem = elementos[iElem];
		corrElemX[iElem] = -sigmas[elem.material] * corrElemX[iElem];
		corrElemY[iElem] = -sigmas[elem.material] * corrElemY[iElem];
	}
}

void Problema::grabar() {
	/* Tensión */
	ofstream tension("tension.csv", ofstream::out);

	tension << "X, Y, V";

	for (int iNodo = 0; iNodo < nNodes; iNodo++) {
		Nodo nodo = nodos[iNodo];
		tension << endl << nodo.x << ", " << nodo.y << ", " << solucion[iNodo];
	}

	tension.close();

	/* Corriente y campo */
	ofstream corriente("corriente.csv", ofstream::out);
	ofstream campo("campo.csv", ofstream::out);

	corriente << "X, Y, corriente";
	campo 	  << "X, Y, campo";

	for (int k = 0; k < nElems; k++) {
		double corr = sqrt(pow( corrElemX[k], 2) + pow( corrElemY[k], 2));
		double camp = sqrt(pow(campoElemX[k], 2) + pow(campoElemY[k], 2));
		double xMed = 0.0, yMed = 0.0;

		for (int j = 0; j < NODPEL; j++) {
			int jNodo = elementos[k][j];
			xMed += nodos[jNodo].x;
			yMed += nodos[jNodo].y;
		}

		xMed /= NODPEL;
		yMed /= NODPEL;

		corriente << "\n" << xMed << ", " << yMed << ", " << corr;
		campo 	  << "\n" << xMed << ", " << yMed << ", " << camp;
	}

	corriente.close();
	campo.close();
}

Elemento Problema::getElement(int i) {
	return elementos[i];
}

void Problema::chequearSimetria() {
	for (int i = 0; i < nNodes; i++) for (int j = i; j < nNodes; j++) {
		if (abs(matriz.coeff(i, j) - matriz.coeff(j, i)) > 1e-9) {
			cerr << "Error: matriz no es simétrica\n";
			exit(EXIT_FAILURE);
		}
	}

	cout << "chequeo completo\n";
}
