#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <cassert>
#include <fstream>
#include <cfloat>
#include <direct.h>
#include "Celula.h"
#include "EntradaSalida.h"

using namespace std;

ofstream historial;
ofstream ph;
FILE*    fPoros;
clock_t  EntradaSalida::start;
bool	 EntradaSalida::firstWriteTransporte = true;
bool	 EntradaSalida::firstWritePoros		 = true;
bool	 EntradaSalida::firstWriteITV		 = true;
bool	 EntradaSalida::firstWritePoisson	 = true;

void EntradaSalida::leerInput(Celula& celula) {
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
				celula.nodpel = 3;
			} else if (s == "quad") {
				celula.nodpel = 4;
			} else {
				assert(false);
			}
		} else if (line.find("salida") != string::npos) {
			iss >> s >> celula.salida;
		} else if (line.find("sigint") != string::npos) {
			iss >> s >> celula.sigmas[INTERNO];
		} else if (line.find("sigext") != string::npos) {
			iss >> s >> celula.sigmas[EXTERNO];
		} else if (line.find("sigmem") != string::npos) {
			iss >> s >> celula.sigmas[MEMBRANA];
		} else if (line.find("potencial") != string::npos) {
			iss >> s >> celula.potencial;
		} else if (line.find("radio") != string::npos) {
			iss >> s >> celula.radio;
		} else if (line.find("ancho") != string::npos) {
			iss >> s >> celula.ancho;
		} else if (line.find("threads") != string::npos) {
			iss >> s >> celula.threads;
		}
	}

	input.close();

	leerMalla(celula, malla);

	celula.alto = -DBL_MAX;
	for (auto nodo : celula.getNodos()) {
		celula.alto = max(celula.alto, nodo.y);
	}

	for (int i = 0; i < celula.nNodes; i++) {
		if (abs(celula.getNodos()[i].y - celula.alto) < EPSILON_DIST) {
			celula.getNodos()[i].esPotencia = true;
		} else if (abs(celula.getNodos()[i].y) < EPSILON_DIST) {
			celula.getNodos()[i].esTierra = true;
		}
	}

	_mkdir(celula.salida.c_str());

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
		celula.getElementos().push_back(Elemento(nod, celula.nodpel, EXTERNO, celula.sigmas[EXTERNO]));
	}

	/* Elementos membrana */
	for (int i = 0; i < elemsMemb; i++) {
		dameLinea(stream, iss);
		iss >> n;
		for (int j = 0; j < celula.nodpel; j++) iss >> nod[j];
		for (int j = 0; j < celula.nodpel; j++) nod[j]--;
		celula.getElementos().push_back(Elemento(nod, celula.nodpel, MEMBRANA, celula.sigmas[MEMBRANA]));
	}

	/* Elementos internos */
	for (int i = 0; i < elemsInt; i++) {
		dameLinea(stream, iss);
		iss >> n;
		for (int j = 0; j < celula.nodpel; j++) iss >> nod[j];
		for (int j = 0; j < celula.nodpel; j++) nod[j]--;
		celula.getElementos().push_back(Elemento(nod, celula.nodpel, INTERNO, celula.sigmas[INTERNO]));
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

void EntradaSalida::grabarPoisson(Celula& celula, double time, bool verbose) {
	if (verbose) printStart("Grabando...", true);

	/* Tensi�n */
	ios_base::open_mode flags = firstWritePoisson ? ios::out : ios::app;
	ofstream tension((celula.salida + "/tension.dat").c_str(), flags);

	tension << "paso: " << time << " " << celula.nNodes << "\n";

	for (int iNodo = 0; iNodo < celula.nNodes; iNodo++) {
		Nodo nodo = celula.getNodos()[iNodo];
		tension << nodo.x << ", " << nodo.y << ", " << celula.getSolucion()[iNodo] << "\n";
	}

	tension.close();

	/* Campo y corriente */
	ofstream campo_corriente((celula.salida + "/campo_corriente.dat").c_str(), flags);

	campo_corriente << "paso: " << time << " " << celula.nElems << "\n";

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

		campo_corriente << xMed << " " << yMed << " " << camp << " " << corr << "\n";
	}

	campo_corriente.close();
	firstWritePoisson = false;
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

void EntradaSalida::grabarTransporte(Celula& cel, double time, bool verbose) {
	if (verbose) EntradaSalida::printStart("Grabando en disco...");

	ios_base::open_mode flags = firstWriteTransporte ? ios::out : ios::app;
	historial.open((cel.salida + "/transporte.dat").c_str(), flags);
	ph.open((cel.salida + "/ph.dat").c_str(), flags);

	assert(historial.is_open());
	assert(ph.is_open());
	ostringstream histSS, phSS;

	if (firstWriteTransporte) {
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

	firstWriteTransporte = false;
	if (verbose) EntradaSalida::printEnd();
}

void EntradaSalida::grabarRadio(Celula& celula, Poros& radios, double time, bool verbose) {
	if (verbose) EntradaSalida::printStart("Grabando en disco...");
	
	char* flags = firstWritePoros ? "w" : "a";
	fPoros = fopen((celula.salida + "/poros.dat").c_str(), flags);
	assert(fPoros > 0);

	fprintf(fPoros, "paso %.9f %d\n", time, radios.getNPoros());

	for (auto& info : radios.getValores()) {
		for (auto& radio : info.porosGrandes) {
			fprintf(fPoros, "%.6f %.6e\n", info.tita, radio.first);
		}

		for (int i = 0; i < info.porosChicos; i++) {
			fprintf(fPoros, "%.6f %.6e\n", info.tita, info.radioChico);
		}
	}

	fclose(fPoros);
	firstWritePoros = false;
	if (verbose) EntradaSalida::printEnd();
}

void EntradaSalida::grabarITV(Celula& celula, Poros& poros, double time) {
	char* flags = firstWriteITV ? "w" : "a";
	FILE* file = fopen((celula.salida + "/itv.dat").c_str(), flags);
	assert(file > 0);

	vector<pair<double, double>> valores = poros.getITVs(time);

	fprintf(file, "paso %.9f %d\n", time, valores.size());

	for (auto& info : valores) {
		fprintf(file, "%.9f, %.9f\n", info.first, info.second);
	}

	fclose(file);
	firstWriteITV = false;
}
