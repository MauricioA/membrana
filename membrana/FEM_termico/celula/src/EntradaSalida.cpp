#include <iostream>
#include <cassert>
#include <fstream>
#include "Celula.h"
#include "EntradaSalida.h"

using namespace std;

//TODO carpeta de salida
//TODO archivo input por parámetros
//TODO ignorar bien los comentarios

ofstream historial;
ofstream ph;
clock_t  EntradaSalida::start;
bool	 EntradaSalida::firstWrite = true;

void EntradaSalida::leerInput(Celula& celula) {
	printStart("Leyendo archivos...", true);

	string s, line, malla;
	vector<int> dirichV, dirichT;
	ifstream input("input.in", ifstream::in);
	double alto;

	assert(input.is_open());

	while (getline(input, line)) {
		istringstream iss(line);

		if (line.find("malla") != string::npos) {
			iss >> s >> malla;
		} else if (line.find("nodpel") != string::npos) {
			iss >> s >> s;

			if (s == "tri") {
				celula.nodpel = 3;
			} else if (s == "quad") {
				celula.nodpel = 4;
			} else {
				assert(false);
			}
		} else if (line.find("sigint") != string::npos) {
			iss >> s >> celula.sigmas[INTERNO];
		} else if (line.find("sigext") != string::npos) {
			iss >> s >> celula.sigmas[EXTERNO];
		} else if (line.find("sigmem") != string::npos) {
			iss >> s >> celula.sigmas[MEMBRANA];
		} else if (line.find("potencial") != string::npos) {
			iss >> s >> celula.potencial;
		} else if (line.find("alto") != string::npos) {
			iss >> s >> alto;
		}
	}

	input.close();

	leerMalla(celula, malla);

	for (int i = 0; i < celula.nNodes; i++) {
		if (abs(celula.getNodos()[i].y - alto) < EPSILON_POISSON) {
			celula.getNodos()[i].esPotencia = true;
		} else if (abs(celula.getNodos()[i].y) < EPSILON_POISSON) {
			celula.getNodos()[i].esTierra = true;
		}
	}

	printEnd();
}

void EntradaSalida::leerMalla(Celula& celula, string malla) {
	int n;
	int nod[4];
	double x, y;
	string line;
	istringstream iss;

	ifstream stream(malla.c_str(), ifstream::in);
	assert(stream.is_open());

	/* Numero de nodos */
	dameLinea(stream, iss);
	iss >> celula.nNodes;
	celula.getNodos().reserve(celula.nNodes);

	/* Nodos */
	for (int i = 0; i < celula.nNodes; i++) {
		dameLinea(stream, iss);
		iss >> n >> x >> y;
		Nodo nodo;
		nodo.x = x;
		nodo.y = y;
		nodo.esPotencia = false;
		nodo.esTierra = false;
		celula.getNodos().push_back(nodo);
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
	celula.nElems = elemsExt + elemsMemb + elemsInt;
	celula.getElementos().reserve(celula.nElems);

	/* Elementos extrenos */
	for (int i = 0; i < elemsExt; i++) {
		dameLinea(stream, iss);
		iss >> n;
		for (int j = 0; j < celula.nodpel; j++) iss >> nod[j];
		for (int j = 0; j < celula.nodpel; j++) nod[j]--;
		celula.getElementos().push_back(Elemento(nod, celula.nodpel, EXTERNO));
	}

	/* Elementos membrana */
	for (int i = 0; i < elemsMemb; i++) {
		dameLinea(stream, iss);
		iss >> n;
		for (int j = 0; j < celula.nodpel; j++) iss >> nod[j];
		for (int j = 0; j < celula.nodpel; j++) nod[j]--;
		celula.getElementos().push_back(Elemento(nod, celula.nodpel, MEMBRANA));
	}

	/* Elementos internos */
	for (int i = 0; i < elemsInt; i++) {
		dameLinea(stream, iss);
		iss >> n;
		for (int j = 0; j < celula.nodpel; j++) iss >> nod[j];
		for (int j = 0; j < celula.nodpel; j++) nod[j]--;
		celula.getElementos().push_back(Elemento(nod, celula.nodpel, INTERNO));
	}

	stream.close();
}

void EntradaSalida::dameLinea(ifstream& archivo, istringstream& iss) {
	string line;
	do {
		getline(archivo, line);
	} while(line[0] == '*');
	iss.clear();
	iss.str(line);
}

void EntradaSalida::grabarPoisson(Celula& celula) {
	printStart("Grabando...", true);

	/* Tensión */
	ofstream tension("tension.csv", ofstream::out);

	tension << "X, Y, V";

	for (int iNodo = 0; iNodo < celula.nNodes; iNodo++) {
		Nodo nodo = celula.getNodos()[iNodo];
		tension << "\n" << nodo.x << ", " << nodo.y << ", " << celula.getSolucion()[iNodo];
	}

	tension.close();

	/* Corriente y campo */
	if (celula.nodpel == 3) {
		ofstream corriente("corriente.csv", ofstream::out);
		ofstream campo("campo.csv", ofstream::out);

		corriente << "X, Y, corriente";
		campo 	  << "X, Y, campo";

		for (int k = 0; k < celula.nElems; k++) {
			double corr = sqrt(pow(celula.getGradElem()[k].x, 2) + pow(celula.getGradElem()[k].y, 2));
			double camp = sqrt(pow(celula.getCorrElem()[k].x, 2) + pow(celula.getCorrElem()[k].y, 2));
			double xMed = 0.0, yMed = 0.0;

			for (int j = 0; j < celula.nodpel; j++) {
				int jNodo = celula.getElementos()[k][j];
				xMed += celula.getNodos()[jNodo].x;
				yMed += celula.getNodos()[jNodo].y;
			}

			xMed /= celula.nodpel;
			yMed /= celula.nodpel;

			corriente << "\n" << xMed << ", " << yMed << ", " << corr;
			campo 	  << "\n" << xMed << ", " << yMed << ", " << camp;
		}

		corriente.close();
		campo.close();
	}

	printEnd(3);
}

void EntradaSalida::printStart(string message, bool verbose) {
	if (verbose) {
		start = clock();
		cout << message << " " << flush;
	}
}

void EntradaSalida::printEnd(int tabs, bool verbose) {
	if (verbose) {
		int time = (clock() - start) / (CLOCKS_PER_SEC / 1000);
		cout << "OK";
		for (int i = 0; i < tabs; i++) cout << "\t";
		cout << time << "ms" << endl;
	}
}

void EntradaSalida::grabarTransporte(Celula& cel, double time) {
	EntradaSalida::printStart("Grabando en disco...");

	_Ios_Openmode flags = firstWrite ? ios::out : ios::app;

	historial.open(RUTA_HISTORIAL.c_str(), flags);
	ph.open(RUTA_PH.c_str(), flags);

	assert(historial.is_open() && ph.is_open());
	ostringstream histSS, phSS;

	if (firstWrite) {
		histSS << cel.nNodes << "\n";
		phSS   << cel.nNodes << "\n";
	}
	histSS << "paso: " << time << "\n";
	phSS   << "paso: " << time << "\n";

	for (int jNodo = 0; jNodo < cel.nNodes; jNodo++) {
		Nodo nodo = cel.getNodos()[jNodo];

		histSS << jNodo+1 << "\t" << nodo.x << "\t" << nodo.y << "\t" << cel.getSolucion()[jNodo];
		phSS   << jNodo+1 << "\t" << nodo.x << "\t" << nodo.y << "\t" << cel.getSolucion()[jNodo];

		for (int esp = 0; esp < NESPS; esp++) {
			histSS << "\t" << cel.concentraciones[esp][jNodo];
		}

		phSS << "\t" << cel.phAux[H_][jNodo] << "\t" << cel.phAux[OH][jNodo];

		histSS << "\n";
		phSS   << "\n";
	}

	historial << histSS.str();
	ph << phSS.str();

	historial.close();
	ph.close();

	firstWrite = false;
	EntradaSalida::printEnd();
}
