#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <cassert>
#include <fstream>
#include <cfloat>
#include <direct.h>
#include "EntradaSalida.h"
#include "Celula.h"

using namespace std;

EntradaSalida::EntradaSalida(Celula& celula) {
	_celula = &celula;

	leerInput();

	ftension	= fopen((getCelula().salida + "/tension.csv").c_str(), "w");
	fitv		= fopen((getCelula().salida + "/itv.csv").c_str(), "w");
	fcampo		= fopen((getCelula().salida + "/campo-corriente.csv").c_str(), "w");
	fporos		= fopen((getCelula().salida + "/poros.csv").c_str(), "w");
	ftransporte = fopen((getCelula().salida + "/transporte.csv").c_str(), "w");
	fph			= fopen((getCelula().salida + "/ph.csv").c_str(), "w");

	assert(ftension > 0);
	assert(fitv > 0);
	assert(fcampo > 0);
	assert(fporos > 0);
	assert(ftransporte > 0);
	assert(fph > 0);
}

EntradaSalida::~EntradaSalida() {
	fclose(ftension);
	fclose(fitv);
	fclose(fcampo);
	fclose(fporos);
	fclose(ftransporte);
	fclose(fph);
}

inline Celula& EntradaSalida::getCelula() {
	return *_celula;
}

void EntradaSalida::leerInput() {
	printStart("Leyendo archivos...", true);

	string s, line, malla;
	vector<int> dirichV, dirichT;
	ifstream input("input.in", ifstream::in);
	assert(input.is_open());

	while (getline(input, line)) {
		istringstream iss(line);

		if (line.find("malla") != string::npos) {
			iss >> s >> malla;
		} else if (line.find("nodpel") != string::npos) {
			iss >> s >> s;

			if (s == "tri") {
				getCelula().nodpel = 3;
			} else if (s == "quad") {
				getCelula().nodpel = 4;
			} else {
				assert(false);
			}
		} else if (line.find("salida") != string::npos) {
			iss >> s >> getCelula().salida;
		} else if (line.find("sigint") != string::npos) {
			iss >> s >> getCelula().sigmas[INTERNO];
		} else if (line.find("sigext") != string::npos) {
			iss >> s >> getCelula().sigmas[EXTERNO];
		} else if (line.find("sigmem") != string::npos) {
			iss >> s >> getCelula().sigmas[MEMBRANA];
		} else if (line.find("potencial") != string::npos) {
			iss >> s >> getCelula().potencial;
		} else if (line.find("radio") != string::npos) {
			iss >> s >> getCelula().radio;
		} else if (line.find("ancho") != string::npos) {
			iss >> s >> getCelula().ancho;
		} else if (line.find("threads") != string::npos) {
			iss >> s >> getCelula().threads;
		}
	}

	input.close();

	leerMalla(malla);

	getCelula().alto = -DBL_MAX;
	for (auto nodo : getCelula().getNodos()) {
		getCelula().alto = max(getCelula().alto, nodo.y);
	}

	for (int i = 0; i < getCelula().nNodes; i++) {
		if (abs(getCelula().getNodos()[i].y - getCelula().alto) < EPSILON_DIST) {
			getCelula().getNodos()[i].esPotencia = true;
		} else if (abs(getCelula().getNodos()[i].y) < EPSILON_DIST) {
			getCelula().getNodos()[i].esTierra = true;
		}
	}

	_mkdir(getCelula().salida.c_str());

	printEnd();
}

void EntradaSalida::leerMalla(string malla) {
	int n;
	int nod[4];
	double x, y;
	string line;
	istringstream iss;

	ifstream stream(malla.c_str(), ifstream::in);
	assert(stream.is_open());

	/* Numero de nodos */
	dameLinea(stream, iss);
	iss >> getCelula().nNodes;
	getCelula().getNodos().reserve(getCelula().nNodes);

	/* Nodos */
	for (int i = 0; i < getCelula().nNodes; i++) {
		dameLinea(stream, iss);
		iss >> n >> x >> y;
		Nodo nodo;
		nodo.x = x;
		nodo.y = y;
		nodo.esPotencia = false;
		nodo.esTierra = false;
		getCelula().getNodos().push_back(nodo);
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
	getCelula().nElems = elemsExt + elemsMemb + elemsInt;
	getCelula().getElementos().reserve(getCelula().nElems);

	/* Elementos extrenos */
	for (int i = 0; i < elemsExt; i++) {
		dameLinea(stream, iss);
		iss >> n;
		for (int j = 0; j < getCelula().nodpel; j++) iss >> nod[j];
		for (int j = 0; j < getCelula().nodpel; j++) nod[j]--;
		getCelula().getElementos().push_back(Elemento(
			nod, getCelula().nodpel, EXTERNO, getCelula().sigmas[EXTERNO]
		));
	}

	/* Elementos membrana */
	for (int i = 0; i < elemsMemb; i++) {
		dameLinea(stream, iss);
		iss >> n;
		for (int j = 0; j < getCelula().nodpel; j++) iss >> nod[j];
		for (int j = 0; j < getCelula().nodpel; j++) nod[j]--;
		getCelula().getElementos().push_back(Elemento(
			nod, getCelula().nodpel, MEMBRANA, getCelula().sigmas[MEMBRANA]
		));
	}

	/* Elementos internos */
	for (int i = 0; i < elemsInt; i++) {
		dameLinea(stream, iss);
		iss >> n;
		for (int j = 0; j < getCelula().nodpel; j++) iss >> nod[j];
		for (int j = 0; j < getCelula().nodpel; j++) nod[j]--;
		getCelula().getElementos().push_back(Elemento(
			nod, getCelula().nodpel, INTERNO, getCelula().sigmas[INTERNO]
		));
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

void EntradaSalida::grabarPoisson(bool verbose) {
	if (verbose) printStart("Grabando...", true);
	double time = getCelula().time;
	int nodpel = getCelula().nodpel;
	int nNodes = getCelula().nNodes;
	int nElems = getCelula().nElems;

	/* Tensión */
	for (int iNodo = 0; iNodo < nNodes; iNodo++) {
		Nodo nodo = getCelula().getNodos()[iNodo];
		fprintf(ftension, "%g,%.6g,%.6g,%.8g\n", time, nodo.x, nodo.y, getCelula().getSolucion()[iNodo]);
	}

	/* Campo y corriente */
	for (int k = 0; k < nElems; k++) {
		double corr = sqrt(pow(getCelula().getGradElem()[k].x, 2) + pow(getCelula().getGradElem()[k].y, 2));
		double camp = sqrt(pow(getCelula().getCorrElem()[k].x, 2) + pow(getCelula().getCorrElem()[k].y, 2));
		double xMed = 0, yMed = 0;

		for (int j = 0; j < nodpel; j++) {
			int jNodo = getCelula().getElementos()[k][j];
			xMed += getCelula().getNodos()[jNodo].x;
			yMed += getCelula().getNodos()[jNodo].y;
		}

		xMed /= nodpel;
		yMed /= nodpel;

		fprintf(fcampo, "%g,%.6g,%.6g,%.9g,%.9g\n", time, xMed, yMed, camp, corr);
	}

	if (verbose) printEnd(3);
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

void EntradaSalida::grabarTransporte(bool verbose) {
	if (verbose) EntradaSalida::printStart("Grabando en disco...");
	double time = getCelula().time;
	int nNodes = getCelula().nNodes;

	for (int jNodo = 0; jNodo < nNodes; jNodo++) {
		Nodo nodo = getCelula().getNodos()[jNodo];

		fprintf(ftransporte, "%g,%.6g,%.6g", time, nodo.x, nodo.y);
		fprintf(fph, "%g,%.6g,%.6g", time, nodo.x, nodo.y);

		for (int esp = 0; esp < NESPS; esp++) {
			fprintf(ftransporte, ",%.9e", getCelula().concentraciones[esp][jNodo]);
		}

		fprintf(fph, ",%.9e", getCelula().phAux[H_][jNodo]);
		fprintf(fph, ",%.9e", getCelula().phAux[OH][jNodo]);

		fprintf(ftransporte, "\n");
		fprintf(fph, "\n");
	}

	if (verbose) EntradaSalida::printEnd(); 
}

void EntradaSalida::grabarRadios(Poros& radios, bool verbose) {
	if (verbose) EntradaSalida::printStart("Grabando en disco...");
	double time = getCelula().time;

	for (auto& info : radios.getValores()) {
		for (auto& radio : info.porosGrandes) {
			fprintf(fporos, "%g,%.6g,%.6e\n", time, info.tita, radio.first);
		}

		for (int i = 0; i < info.porosChicos; i++) {
			fprintf(fporos, "%g,%.6g,%.6e\n", time, info.tita, info.radioChico);
		}
	}

	if (verbose) EntradaSalida::printEnd();
}

void EntradaSalida::grabarITV(Poros& poros) {
	double time = getCelula().time;
	vector<pair<double, double>> valores = poros.getITVs(time);

	for (auto& info : valores) {
		fprintf(fitv, "%g,%.9g,%.9g\n", time, info.first, info.second);
	}
}
