#include <iostream>
#include <cassert>
#include <fstream>
#include "Celula.h"
#include "EntradaSalida.h"

using namespace std;

void EntradaSalida::leerInput(Celula &celula) {
	cout << "Leyendo archivos... " << flush;
	clock_t start = clock();

	string s, line, malla;
	vector<int> dirichV, dirichT;
	ifstream input("input.in", ifstream::in);
	double alto;
	celula.nodpel = 3;

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

	int time = (clock() - start) / (CLOCKS_PER_SEC / 1000);
	cout << "OK\t\t" << time << "ms\n";
}

void EntradaSalida::leerMalla(Celula &celula, string malla) {
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
