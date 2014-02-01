#include <iostream>
#include <fstream>
#include <ctime>
#include <cassert>

#include <Eigen/LU>

#include "Problema.h"

//TODO ver cosa rara en armado4!!!
//TODO transporte
//TODO corriente y campo para quad
//TODO multiplicar corriente para obtener por m o cm
//TODO carpeta de salida
//TODO archivo input por parámetros
//TODO fort: input y qe = 0 -> ef
//TODO ignorar bien los comentarios
//TODO refactorizar a varios archivos/clases?
//TODO PERFO: ver como graba en windows
//TODO PERFO: compilar con NDEBUG
//TODO PERFO: cambiar flags de vectorización
//TODO comparar resultados de campo para tri y quad con fortran

Problema::Problema() {
	cout << "Leyendo archivos... " << flush;
	clock_t start = clock();

	string s, line, malla;
	vector<int> dirichV, dirichT;
	ifstream input("input.in", ifstream::in);
	double alto;
	nodpel = 3;

	assert(input.is_open());

	while (getline(input, line)) {
		istringstream iss(line);

		if (line.find("malla") != string::npos) {
			iss >> s >> malla;
		} else if (line.find("nodpel") != string::npos) {
			iss >> s >> s;

			if (s == "tri") {
				nodpel = 3;
			} else if (s == "quad") {
				nodpel = 4;
			}
		} else if (line.find("sigint") != string::npos) {
			iss >> s >> sigmas[INTERNO];
		} else if (line.find("sigext") != string::npos) {
			iss >> s >> sigmas[EXTERNO];
		} else if (line.find("sigmem") != string::npos) {
			iss >> s >> sigmas[MEMBRANA];
		} else if (line.find("potencial") != string::npos) {
			iss >> s >> potencial;
		} else if (line.find("alto") != string::npos) {
			iss >> s >> alto;
		}
	}

	input.close();

	leerMalla(malla);

	for (int i = 0; i < nNodes; i++) {
		if (abs(nodos[i].y - alto) < EPSILON_POISSON) {
			nodos[i].esPotencia = true;
		} else if (abs(nodos[i].y) < EPSILON_POISSON) {
			nodos[i].esTierra = true;
		}
	}

	int time = (clock() - start) / (CLOCKS_PER_SEC / 1000);
	cout << "OK\t\t" << time << "ms\n";
}

void Problema::leerMalla(string malla) {
	int n;
	int nod[4];
	double x, y;
	string line;
	istringstream iss;

	ifstream stream(malla.c_str(), ifstream::in);
	assert(stream.is_open());

	/* Numero de nodos */
	dameLinea(stream, iss);
	iss >> nNodes;
	nodos.reserve(nNodes);

	/* Nodos */
	for (int i = 0; i < nNodes; i++) {
		dameLinea(stream, iss);
		iss >> n >> x >> y;
		Nodo nodo;
		nodo.x = x;
		nodo.y = y;
		nodo.esPotencia = false;
		nodo.esTierra = false;
		nodos.push_back(nodo);
	}

	/* nGrupos */
	int nGrupos;
	dameLinea(stream, iss);
	iss >> nGrupos;

	/* Grupos */
	int elemsExt, elemsMemb, elemsInt;
	dameLinea(stream, iss);
	iss >> n >> elemsExt;
	dameLinea(stream, iss);
	iss >> n >> elemsMemb;
	dameLinea(stream, iss);
	iss >> n >> elemsInt;
	nElems = elemsExt + elemsMemb + elemsInt;
	elementos.reserve(nElems);

	/* Elementos extrenos */
	for (int i = 0; i < elemsExt; i++) {
		dameLinea(stream, iss);
		iss >> n;
		for (int j = 0; j < nodpel; j++) iss >> nod[j];
		for (int j = 0; j < nodpel; j++) nod[j]--;
		elementos.push_back(Elemento(nod, nodpel, EXTERNO));
	}

	/* Elementos membrana */
	for (int i = 0; i < elemsMemb; i++) {
		dameLinea(stream, iss);
		iss >> n;
		for (int j = 0; j < nodpel; j++) iss >> nod[j];
		for (int j = 0; j < nodpel; j++) nod[j]--;
		elementos.push_back(Elemento(nod, nodpel, MEMBRANA));
	}

	/* Elementos internos */
	for (int i = 0; i < elemsInt; i++) {
		dameLinea(stream, iss);
		iss >> n;
		for (int j = 0; j < nodpel; j++) iss >> nod[j];
		for (int j = 0; j < nodpel; j++) nod[j]--;
		elementos.push_back(Elemento(nod, nodpel, INTERNO));
	}

	stream.close();
}

void Problema::dameLinea(ifstream& archivo, istringstream& iss) {
	string line;
	do {
		getline(archivo, line);
	} while(line[0] == '*');
	iss.clear();
	iss.str(line);
}

void Problema::poisson() {
	rhs.resize(nNodes);
	rhs.fill(0.0);
	solucion.resize(nNodes);
	matriz.resize(nNodes, nNodes);

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

			armado(pos, esm, sigma);

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

void Problema::armado(Double2D pos[], double esm[][MAXNPEL], double sigma) {
	switch (nodpel) {
	case 3:
		armado3(pos, esm, sigma);
		break;
	case 4:
		armado4(pos, esm, sigma, false, 0., 0., NULL, NULL);
		break;
	}
}

void Problema::armado3(Double2D pos[], double esm[][MAXNPEL], double sigma) {
	assert(nodpel == 3);

	const int NODPEL = 3;
	double b[3], c[3];

	double det = determinante3(pos, b, c);
	double rMed = (pos[0].x + pos[1].x + pos[2].x) / 3;

	assert(abs(det) > TOLER_AREA);

	for (int i = 0; i < NODPEL; i++) for (int j = 0; j < 3; j++) {
		esm[i][j] = sigma * (b[i] * b[j] + c[i] * c[j]) * M_PI * rMed / det;
	}
}

double Problema::determinante3(Double2D pos[], double b[], double c[]) {
	assert(nodpel == 3);

	int i = 0;
	b[i++] = pos[1].y - pos[2].y;
	b[i++] = pos[2].y - pos[0].y;
	b[i++] = pos[0].y - pos[1].y;

	i = 0;
	c[i++] = pos[2].x - pos[1].x;
	c[i++] = pos[0].x - pos[2].x;
	c[i++] = pos[1].x - pos[0].x;

	return
		+ pos[1].x*pos[2].y + pos[2].x*pos[0].y + pos[0].x*pos[1].y
		- pos[1].x*pos[0].y - pos[2].x*pos[1].y - pos[0].x*pos[2].y;
}

void Problema::armado4(Double2D pos[], double esm[][MAXNPEL], double sigma,
		bool transp, double landa, double mu, double est[][4], double mas[]) {
	assert(nodpel == 4);

	const int NODPEL = 4;
	double phi[2*NGAUSS][NODPEL];
	double dphi[NDIM][2*NGAUSS][NODPEL];
	double gxCod[NDIM][2*NGAUSS];
	double phidX[NDIM][2*NGAUSS][4];
	double cteI[2*NGAUSS];
	int kGauss = 0;

	for (int i = 0; i < NGAUSS; i++) for (int j = 0; j < NGAUSS; j++) {
		double det = iteracion4(phi, dphi, phidX, i, j, kGauss, pos);

		for (int dim = 0; dim < NDIM; dim++) {
			gxCod[dim][kGauss] = 0.0;

			for (int k = 0; k < NODPEL; k++) {
				gxCod[dim][kGauss] += pos[k].x * phi[kGauss][k];
			}
		}

		cteI[kGauss] = det * GAUSSWT[i] * GAUSSWT[j] * 2 * M_PI * gxCod[0][kGauss];
		kGauss++;
	}

	for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) esm[i][j] = 0.;

	for (int kGauss = 0; kGauss < NGAUSS*NGAUSS; kGauss++) for (int i = 0; i < nodpel; i++) {
		for (int j = 0; j < nodpel; j++) {
			for (int d = 0; d < NDIM; d++) {
				esm[i][j] += sigma * phidX[d][kGauss][i] * phidX[d][kGauss][j] * cteI[kGauss];
			}

			if (transp) {
//				TODO esto lo calcula de más porque despues lo sobrescribe en armadoTransporte()!!!
				esm[i][j] += mu * phidX[0][kGauss][i] * phi[kGauss][j] * cteI[kGauss];
			}
		}
		if (transp) {
			est[i][i] += mas[i] * landa * cteI[kGauss];
		}
	}
}

double Problema::iteracion4(double phi[2*NGAUSS][4], double dphi[NDIM][2*NGAUSS][4],
		double phidX[NDIM][2*NGAUSS][4], int i, int j, int kGauss, Double2D pos[4]) {
	assert(nodpel == 4);

	double t = GAUSSPT[i];
	double s = GAUSSPT[j];

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
	dphi[0][kGauss][k++] =  0.5 * tm;
	dphi[0][kGauss][k++] =  0.5 * tp;
	dphi[0][kGauss][k++] = -0.5 * tp;

	k = 0;
	dphi[1][kGauss][k++] = -0.5 * sm;
	dphi[1][kGauss][k++] = -0.5 * sq;
	dphi[1][kGauss][k++] =  0.5 * sq;
	dphi[1][kGauss][k++] =  0.5 * sm;

	double aJacob[2][2];
	for (int k = 0; k < 2; k++) for (int l = 0; l < 2; l++) aJacob[k][l] = 0.0;

	for (k = 0; k < nodpel; k++) {
		aJacob[0][0] += dphi[0][kGauss][k] * pos[k].x;
		aJacob[0][1] += dphi[0][kGauss][k] * pos[k].y;
		aJacob[1][0] += dphi[1][kGauss][k] * pos[k].x;
		aJacob[1][1] += dphi[1][kGauss][k] * pos[k].y;
	}

	double det =
		aJacob[0][0] * aJacob[1][1] -
		aJacob[0][1] * aJacob[1][0];

	double aJacobInv[2][2] = {
		{  aJacob[1][1] / det, -aJacob[0][1] / det, },
		{ -aJacob[1][0] / det,  aJacob[0][0] / det,	},
	};

	for (int k = 0; k < nodpel; k++) {
		phidX[0][kGauss][k] =
			aJacobInv[0][0] * dphi[0][kGauss][k] +
			aJacobInv[0][1] * dphi[1][kGauss][k];

		phidX[1][kGauss][k] =
			aJacobInv[1][0] * dphi[0][kGauss][k] +
			aJacobInv[1][1] * dphi[1][kGauss][k];
	}

	return det;
}

void Problema::campo() {
	switch (nodpel) {
	case 3:
		campo3();
		break;
	case 4:
		campo4();
		break;
	}
}

void Problema::campo3() {
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

		double det = determinante3(pos, b, c);

		gradElem[iElem].x = (b[0] * sol[0] + b[1] * sol[1] + b[2] * sol[2]) / det;
		gradElem[iElem].y = (c[0] * sol[0] + c[1] * sol[1] + c[2] * sol[2]) / det;
	}
}

void Problema::campo4() {
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
			iteracion4(phi, dphi, phidX, i, j, kGauss, pos);

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
void Problema::corriente() {
	corrElem.resize(nElems);

	for (int iElem = 0; iElem < nElems; iElem++) {
		Elemento elem = elementos[iElem];
		corrElem[iElem].x = -sigmas[elem.material] * gradElem[iElem].x;
		gradElem[iElem].y = -sigmas[elem.material] * gradElem[iElem].y;
	}
}

void Problema::grabar() {
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

Elemento Problema::getElement(int i) {
	return elementos[i];
}

void Problema::chequearSimetria() {
	for (int i = 0; i < nNodes; i++) for (int j = i; j < nNodes; j++) {
		assert(abs(matriz.coeff(i, j) - matriz.coeff(j, i)) < 1e-9);
	}
}

void Problema::transporte() {
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

		phAux[H_][iNode] = -log10(concentraciones[H_][iNode] * 1e15 / 6.02e23);
		phAux[OH][iNode] = -log10(concentraciones[OH][iNode] * 1e15 / 6.02e23);
	}

	masaDiag2D();
	cargas.resize(nNodes);

	for (double tt = 0.; tt < T_CERO; tt += DELTA_T) {
		poisson();


	}

}

void Problema::masaDiag2D() {
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

void Problema::carga() {
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

void Problema::concentracion(int esp) {
	double esm[nodpel][MAXNPEL];
	vector< Triplet<double> > triplets;

	for (uint kElem = 0; kElem < elementos.size(); kElem++) {
		Elemento elem = elementos[kElem];
		double yMed = 0.;
		Double2D pos[nodpel];
		double sol[nodpel];
		double ef[nodpel];

		for (int i = 0; i < nodpel; i++) {
			int iNodo = elem[i];
			Nodo nodo = nodos[iNodo];
			yMed += nodo.y;
			pos[i].x = nodo.x;
			pos[i].y = nodo.y;
			sol[i] = solucion[iNodo];
		}

		yMed /= nodpel;
		double sigma = sigmas[elem.material];

		double difElem = DIFUSION[esp];
		if (elem.material == MEMBRANA) difElem *= 1e-3;

		double mu = -difElem * CLAVE * CARGA[esp] * gradElem[kElem].y;

		double landa = 1.;

		armadoTransporte(pos, esm, sigma, landa, mu, sol, ef);

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

		for (int iNodo = 0; iNodo < nNodes; iNodo++) {
			concentraciones[esp][iNodo] = rhs[iNodo];
		}

		matriz.setFromTriplets(triplets.begin(), triplets.end());
		matriz.makeCompressed();

		/* Resolución */
		SparseLU< SparseMatrix<double, ColMajor>, COLAMDOrdering<Index> > solver;


//		SimplicialLDLT< SparseMatrix<double> > cholesky(matriz);
//		solucion = cholesky.solve(rhs);
//
//		time = (clock() - start) / (CLOCKS_PER_SEC / 1000);
//		cout << "OK\t\t" << time << "ms\n";
//
//		//TODO siempre hace una sola iteración
//		error = EPSILON_POISSON * .5;
	}

}

void Problema::armadoTransporte(Double2D pos[], double esm[][MAXNPEL], double sigma, double landa,
		double mu, double sol[], double ef[]) {
	assert(nodpel == 4);	//not implemented para nodpel == 3

	const int NODPEL = 4;
	const double th2 = 0.5;
	const double aCoef1 = DELTA_T * th2;
	const double aCoef2 = DELTA_T * (1. - th2);

	double est[NODPEL][NODPEL];
	double mas[NODPEL];

	for (int i = 0; i < NODPEL; i++) {
		mas[i] = 0.;
		for (int j = 0; j < NODPEL; j++) est[i][j] = 0.;
	}

	armado4(pos, esm, sigma, true, landa, mu, est, mas);

	for (int k = 1; k < NODPEL; k++) {
		double sum = 0.;
		for (int j = 0; j < NODPEL; j++) {
			sum += (est[k][j] - aCoef2 * esm[k][j]) * sol[j];
			esm[k][j] = est[k][j] + aCoef1 * esm[k][j];
		}
		ef[k] = (aCoef1 + aCoef2) * ef[k] + sum;
	}
}
