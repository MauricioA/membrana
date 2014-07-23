#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <fstream>
#include <direct.h>
#include <cassert>
#include <cfloat>
#include "EntradaSalida.h"
#include "Celula.h"

using namespace std;

EntradaSalida::EntradaSalida(Celula& celula) {
	_celula = &celula;
	leerInput();

	/* Copiar input.in */
	ifstream src("input.in", ios::binary);
	ofstream dst(getCelula().salida + "/input.in", ios::binary);
	dst << src.rdbuf();

	fitv = fopen((getCelula().salida + "/itv.csv").c_str(), "w");
	assert(fitv > 0);
}

EntradaSalida::~EntradaSalida() {
	fclose(fitv);
}

inline Celula& EntradaSalida::getCelula() {
	return *_celula;
}

void EntradaSalida::leerInput() {
	printStart("Leyendo archivos... ", true);

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
	for (auto nodo : getCelula().nodos) {
		getCelula().alto = max(getCelula().alto, nodo.y);
	}

	for (int i = 0; i < getCelula().nNodes; i++) {
		if (abs(getCelula().nodos[i].y - getCelula().alto) < EPSILON_DIST) {
			getCelula().nodos[i].esPotencia = true;
		} else if (abs(getCelula().nodos[i].y) < EPSILON_DIST) {
			getCelula().nodos[i].esTierra = true;
		}
	}

	_mkdir(getCelula().salida.c_str());
	_mkdir((getCelula().salida + "/tension").c_str());
	_mkdir((getCelula().salida + "/campo").c_str());
	_mkdir((getCelula().salida + "/transporte").c_str());
	_mkdir((getCelula().salida + "/poros").c_str());

	printEnd(2, true);
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
	getCelula().nodos.reserve(getCelula().nNodes);

	/* Nodos */
	for (int i = 0; i < getCelula().nNodes; i++) {
		dameLinea(stream, iss);
		iss >> n >> x >> y;
		Nodo nodo;
		nodo.x = x;
		nodo.y = y;
		nodo.esPotencia = false;
		nodo.esTierra = false;
		getCelula().nodos.push_back(nodo);
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
	getCelula().elementos.reserve(getCelula().nElems);

	/* Elementos extrenos */
	for (int i = 0; i < elemsExt; i++) {
		dameLinea(stream, iss);
		iss >> n;
		for (int j = 0; j < getCelula().nodpel; j++) iss >> nod[j];
		for (int j = 0; j < getCelula().nodpel; j++) nod[j]--;
		getCelula().elementos.push_back(Elemento(
			nod, getCelula().nodpel, EXTERNO, getCelula().sigmas[EXTERNO]
		));
	}

	/* Elementos membrana */
	for (int i = 0; i < elemsMemb; i++) {
		dameLinea(stream, iss);
		iss >> n;
		for (int j = 0; j < getCelula().nodpel; j++) iss >> nod[j];
		for (int j = 0; j < getCelula().nodpel; j++) nod[j]--;
		getCelula().elementos.push_back(Elemento(
			nod, getCelula().nodpel, MEMBRANA, getCelula().sigmas[MEMBRANA]
		));
	}

	/* Elementos internos */
	for (int i = 0; i < elemsInt; i++) {
		dameLinea(stream, iss);
		iss >> n;
		for (int j = 0; j < getCelula().nodpel; j++) iss >> nod[j];
		for (int j = 0; j < getCelula().nodpel; j++) nod[j]--;
		getCelula().elementos.push_back(Elemento(
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
	printStart("Grabando poisson...", verbose);

	double time = getCelula().time;
	int nodpel = getCelula().nodpel;

	FILE* ftension = fopen((
		getCelula().salida + "/tension/tension-" + 
		to_string(nPoisson) + "-" + to_string(time) + ".csv"
	).c_str(), "w");

	assert(ftension > 0);
	fprintf(ftension, "         X,          Y,         tension\n");

	for (int iNodo = 0; iNodo < getCelula().nNodes; iNodo++) {
		Nodo nodo = getCelula().nodos[iNodo];
		fprintf(ftension, "%#10f, %#10f, %#15.9g\n", nodo.x, nodo.y, getCelula().solucion[iNodo]);
	}

	fclose(ftension);

	FILE* fcampo = fopen((
		getCelula().salida + "/campo/campo-" +
		to_string(nPoisson) + "-" + to_string(time) + ".csv"
	).c_str(), "w");
	
	assert(fcampo > 0);
	fprintf(ftension, "         X,          Y,           campo,       corriente\n");

	for (int k = 0; k < getCelula().nElems; k++) {
		double corr = sqrt(pow(getCelula().corrElem[k].x, 2) + pow(getCelula().gradElem[k].y, 2));
		double camp = sqrt(pow(getCelula().corrElem[k].x, 2) + pow(getCelula().corrElem[k].y, 2));
		double xMed = 0, yMed = 0;
	
		for (int j = 0; j < getCelula().nodpel; j++) {
			int jNodo = getCelula().elementos[k][j];
			xMed += getCelula().nodos[jNodo].x;
			yMed += getCelula().nodos[jNodo].y;
		}
	
		xMed /= nodpel;
		yMed /= nodpel;
		
		fprintf(fcampo, "%#10f, %#10f, %#15.9g, %#15.9g\n", xMed, yMed, camp, corr);
	}

	fclose(fcampo);

	nPoisson++;
	printEnd(2, verbose);
}

void EntradaSalida::printStart(string message, bool verbose) {
	if (verbose) {
		start = chrono::high_resolution_clock::now();
		cout << message << " " << flush;
	}
}

void EntradaSalida::printEnd(int tabs, bool verbose) {
	if (verbose) {
		auto end = chrono::high_resolution_clock::now();
		auto delta_ms = chrono::duration_cast<chrono::milliseconds>(end - start).count();
		cout << "OK";
		for (int i = 0; i < tabs; i++) cout << "\t";
		cout << delta_ms << "ms" << endl;
	}
}

void EntradaSalida::grabarTransporte(bool verbose) {
	printStart("Grabando transporte...", verbose);
	
	double time = getCelula().time;
	int nNodes = getCelula().nNodes;

	FILE* fpH = fopen((
		getCelula().salida + "/transporte/pH-" +
		to_string(nTransporte) + "-" + to_string(time) + ".csv"
	).c_str(), "w");
	
	FILE* fTrans = fopen((
		getCelula().salida + "/transporte/concent-" +
		to_string(nTransporte) + "-" + to_string(time) + ".csv"
	).c_str(), "w");
	
	assert(fpH > 0);
	assert(fTrans > 0);

	fprintf(fpH,    "         X,          Y,           pH_H+,          pH_OH-\n");
	fprintf(fTrans, "         X,          Y,             Na+,             Cl-\n");

	for (int jNodo = 0; jNodo < nNodes; jNodo++) {
		Nodo nodo = getCelula().nodos[jNodo];
		
		double pH_H  = -log10((getCelula().concentraciones[H_][jNodo] + 1e-18) / CONCENT);
		double ph_OH = -log10((getCelula().concentraciones[OH][jNodo] + 1e-18) / CONCENT);

		double molarNa = getCelula().concentraciones[NA][jNodo] / CONCENT;
		double molarCl = getCelula().concentraciones[CL][jNodo] / CONCENT;

		fprintf(fpH, "%#10f, %#10f, %#15.9g, %#15.9g\n", nodo.x, nodo.y, pH_H, ph_OH);
		fprintf(fTrans, "%#10f, %#10f, %#15.9g, %#15.9g\n", nodo.x, nodo.y, molarNa, molarCl);
	}

	fclose(fpH);
	fclose(fTrans);
	nTransporte++;
	printEnd(1, verbose);
}

void EntradaSalida::grabarPoros(Poros& radios, bool verbose) {
	printStart("Grabando poros...", verbose);
	double time = getCelula().time;

	FILE* fPoros = fopen((
		getCelula().salida + "/poros/poros-" +
		to_string(nPoros) + "-" + to_string(time) + ".csv"
	).c_str(), "w");

	assert(fPoros > 0);
	fprintf(fPoros, "      tita,           radio\n");

	for (auto& info : radios.getValores()) {
		for (auto& radio : info.porosGrandes) {
			fprintf(fPoros, "%#10f, %#15.9g\n", info.tita, radio.first);
		}

		for (int i = 0; i < info.porosChicos; i++) {
			fprintf(fPoros, "%#10f, %#15.9g\n", time, info.tita, info.radioChico);
		}
	}

	nPoros++;
	fclose(fPoros);
	printEnd(2, verbose);
}

void EntradaSalida::grabarITV(Poros& poros) {
	double time = getCelula().time;
	vector<pair<double, double>> valores = poros.getITVs(time);

	for (auto& info : valores) {
		fprintf(fitv, "%#10.6g, %#10f, %#10f\n", time, info.first, info.second);
	}
}
