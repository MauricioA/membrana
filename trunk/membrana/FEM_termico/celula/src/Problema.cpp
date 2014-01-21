#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>

#include "Problema.h"

//TODO capacidad inicial de los vectores con reserve

Problema::Problema() {
	//TODO valores mockeados. Leer de input.in
	sigmas[INTERNO]  = 150e-9;
	sigmas[EXTERNO]  = 200e-9;
	sigmas[MEMBRANA] = 500e-15;
	potencial = 1.;

	leerMalla("celula.fem");

	nodos[   2].esPotencia = true;
	nodos[5005].esPotencia = true;
	nodos[5062].esPotencia = true;
	nodos[5076].esPotencia = true;

	nodos[   1].esPotencia = true;
	nodos[5004].esPotencia = true;
	nodos[5061].esPotencia = true;
	nodos[5070].esPotencia = true;
}

void Problema::leerMalla(const char* archivo) {
	int n, nodo1, nodo2, nodo3;
	double x, y;
	int max = 256;
	char line[max];
	ifstream stream(archivo, ifstream::in);

	/* *COORDINATES */
	stream.getline(line, max);

	/* Numero de nodos */
	stream >> nNodes;

	/* Nodos */
	for (int i = 0; i < nNodes; i++) {
		stream >> n >> x >> y;
		Nodo nodo;
		nodo.x = x;
		nodo.y = y;
		nodo.esPotencia = false;
		nodo.esTierra = false;
		nodo.potencial = 0.;
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

	nElems = elementos.size();

	stream.close();
}

void Problema::poisson() {
	control();

	double error = 1.;
	double ef[NODPEL];
	vector< Triplet<double> > triplets;

	for (int contador = 0; error > EPSILON && contador < N_COTA; contador++) {

		for (uint elemIdx = 0; elemIdx < elementos.size(); elemIdx++) {

			Elemento elemento = elementos[elemIdx];
			double sigma = sigmas[elemento.material];
			double qe = 0.;
			double x[NODPEL], y[NODPEL];
			double esm[3][3];

			for (int i = 0; i < NODPEL; i++) {
				int j = elemento[i];
				x[i] = nodos[j].x;
				y[i] = nodos[j].y;
			}

//			gradElemX[elemIdx] = gradElemY[elemIdx] = 0.;

			armado(x, y, ef, qe, esm, sigma);

			/* Condiciones de contorno */
			for (int i = 0; i < NODPEL; i++) {
				Nodo nodo = nodos[elemento[i]];
				double adiag = esm[i][i];

				if (nodo.esTierra) {
					for (int j = 0; j < NODPEL; j++) {
						esm[i][j] = 0.;
						ef[j] -= esm[j][i] * TIERRA;
						esm[j][i] = 0.;
					}
					esm[i][i] = adiag;
					ef[i] = adiag * TIERRA;
				}

				if (nodo.esPotencia) {
					for (int j = 0; j < NODPEL; j++) {
						esm[i][j] = 0.;
						ef[j] -= esm[j][i] * potencial;
						esm[j][i] = 0.;
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

		// chequearSimetria();

		/* Resolución */

		//TODO esto no sirve para nada
		//solucion = rhs;

		ConjugateGradient< SparseMatrix<double> > cg(matriz);

		cout << "resolviendo... ";

		solucion = cg.solve(rhs);

		cout << "done\n";

		double r = cg.error();
		cout << "error: " << r << endl;

		//TODO no tiene sentido, siempre hace una sola iteración
		error = EPSILON * .5;
	}

	VectorXd b1 = matriz * solucion;
	for (int i = 0; i < nNodes; i++) {
		if (abs(rhs.coeff(i) - b1.coeff(i)) > 1e-6) cerr << "ERROR!" << endl;
		cout << "(" << rhs.coeff(i) << ", " << b1.coeff(i) << ") " << endl;
	}

	corriente2D();

	grabar();
}

void Problema::control() {
	rhs.resize(nNodes);
	rhs.fill(0.0);
	solucion.resize(nNodes);
	solucion.fill(0.0);
	matriz.resize(nNodes, nNodes);
//	gradElemX.resize(nNodes);
//	gradElemY.resize(nNodes);
}

void Problema::armado(double x[], double y[], double ef[], double qe, double esm[3][3], double sigma) {
	double b[] = {
		y[1] - y[2],
		y[2] - y[0],
		y[0] - y[1],
	};

	double c[] = {
		x[2] - x[1],
		x[0] - x[2],
		x[1] - x[0],
	};

	double determinante =
		+ x[1]*y[2] + x[2]*y[0] + x[0]*y[1]
		- x[1]*y[0] - x[2]*y[1] - x[0]*y[2];

	double rMed = (x[0] + x[1] + x[2]) / 3.;

	if (abs(determinante) < TOLER_AREA) {
		cerr << "Error area es cero\n";
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < 3; i++) {
		ef[i] = 0.0;
		for (int j = 0; j < 3; j++) {
			double a = i == j ? 2. : 1.;
			ef[i] += (qe * determinante * M_PI / 12.) * (a * x[j]);
			esm[i][j] = (sigma * (b[i] * b[j] + c[i] * c[j])) * M_PI * rMed * determinante;
		}
	}
}

void Problema::corriente2D() {
	double x[3], y[3], sol[3];

	gradElemX.resize(nElems);
	gradElemY.resize(nElems);
	for (int i = 0; i < nElems; i++) {
		gradElemX[i] = 0.0;
		gradElemY[i] = 0.0;
	}

	gradX.resize(nNodes);
	gradY.resize(nNodes);
	for (int i = 0; i < nNodes; i++) {
		gradX[i] = 0.0;
		gradY[i] = 0.0;
	}

	for (int iElem = 0; iElem < nElems; iElem++) {
		Elemento elem = elementos[iElem];

		for (int i = 0; i < NODPEL; i++) {
			x[i] = nodos[elem[i]].x;
			y[i] = nodos[elem[i]].x;
			sol[i] = solucion.coeff(elem[i]);
		}

		double b[] = {
			y[1] - y[2],
			y[2] - y[0],
			y[0] - y[1],
		};

		double c[] = {
			x[2] - x[1],
			x[0] - x[2],
			x[1] - x[0],
		};

		double determinante =
			+ x[1]*y[2] + x[2]*y[0] + x[0]*y[1]
			- x[1]*y[0] - x[2]*y[1] - x[0]*y[2];

		gradElemX[iElem] = (b[0] * sol[0] + b[1] * sol[1] + b[2] * sol[2]) / determinante;
		gradElemY[iElem] = (c[0] * sol[0] + c[1] * sol[1] + c[2] * sol[2]) / determinante;

		gradElemX[iElem] = sigmas[elem.material] * gradElemX[iElem] * -1;
		gradElemY[iElem] = sigmas[elem.material] * gradElemY[iElem] * -1;

		for (int i = 0; i < NODPEL; i++) {
			int jNodo = elem[i];
			//TODO ESTO ESTÁ MAL SEGURO!!!!
			gradX[jNodo] += gradElemX[i] / NODPEL;
			gradY[jNodo] += gradElemY[i] / NODPEL;
		}
	}
}

void Problema::grabar() {
	/* Tensión */
	ofstream tension("tension.csv", ofstream::out);

	tension << "X, Y, V\n";

	for (int iNodo = 0; iNodo < nNodes; iNodo++) {
		Nodo nodo = nodos[iNodo];
		tension << nodo.x << ", " << nodo.y << ", " << solucion[iNodo] << "\n";
	}

	tension.close();
}

Elemento Problema::getElement(int i) {
	return elementos[i];
}

void Problema::chequearSimetria() {
	for (int i = 0; i < nNodes; i++) {
		for (int j = i; j < nNodes; j++) {
			if (abs(matriz.coeff(i, j) - matriz.coeff(j, i)) > 1e-4) {
				cout << "i = " << i << " j = " << j << " " << matriz.coeff(i, j) << " " << matriz.coeff(j, i) << endl;
			}
		}
	}

	cout << "chequeo completo\n";
}
