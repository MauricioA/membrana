#ifdef _WIN32
#define _CRT_SECURE_NO_WARNINGS
#include <direct.h>
#endif

#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <cassert>
#include <cfloat>
#include "EntradaSalida.h"
#include "Celula.h"

using namespace std;

const double EPSILON_DIST = 1e-9;

EntradaSalida::EntradaSalida(Celula& celula) {
	#ifdef EIGEN_VECTORIZE
		cout << "Vectorization ON\n";
	#else
		cout << "Vectorization OFF\n";
	#endif

	_celula = &celula;
	leerInput();
	
	nPoisson = 0;
	nTransporte = 0;
	nPoros = 0;

	ifstream f(getCelula().salida + "/input.in");
	assert(!f.good() && "directorio de salida not empty!");
	f.close();

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
		} else if (line.find("problema") != string::npos) {
			iss >> s >> s;
			if (s == "potencial") {
				getCelula().soloPoisson = true;
				getCelula().calcularPoros = false;
				getCelula().calcularTransporte = false;
			} else if (s == "poros") {
				getCelula().soloPoisson = false;
				getCelula().calcularPoros = true;
				getCelula().calcularTransporte = false;
			} else if (s == "transporte") {
				getCelula().soloPoisson = false;
				getCelula().calcularPoros = false;
				getCelula().calcularTransporte = true;
			} else if (s == "acoplado") {
				getCelula().soloPoisson = false;
				getCelula().calcularPoros = true;
				getCelula().calcularTransporte = true;
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
		} else if (line.find("delta_t") != string::npos) {
			iss >> s >> getCelula().delta_t;
		} else if (line.find("pulsos") != string::npos) {
			iss >> s >> getCelula().pulsos;
		} else if (line.find("on_time") != string::npos) {
			iss >> s >> getCelula().times[ON];
		} else if (line.find("off_time") != string::npos) {
			iss >> s >> getCelula().times[OFF];
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
	
	#ifdef _WIN32
		_mkdir(getCelula().salida.c_str());
        _mkdir((getCelula().salida + "/tension").c_str());
        _mkdir((getCelula().salida + "/campo").c_str());
        _mkdir((getCelula().salida + "/transporte").c_str());
        _mkdir((getCelula().salida + "/poros").c_str());
	#else
		mkdir(getCelula().salida.c_str(), 0777);
        mkdir((getCelula().salida + "/tension").c_str(), 0777);
        mkdir((getCelula().salida + "/campo").c_str(), 0777);
        mkdir((getCelula().salida + "/transporte").c_str(), 0777);
        mkdir((getCelula().salida + "/poros").c_str(), 0777);
	#endif

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
	assert(nGrupos == 3);

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

	int nodpel = getCelula().nodpel;
	char buffer[512];

	sprintf(buffer, "%s/tension/tension.csv.%03d", getCelula().salida.c_str(), nPoisson);
	FILE* ftension = fopen(buffer, "w");

	assert(ftension > 0);
	fprintf(ftension, "         X,          Y,         tension\n");

	for (int iNodo = 0; iNodo < getCelula().nNodes; iNodo++) {
		Nodo nodo = getCelula().nodos[iNodo];
		fprintf(ftension, "%#10f, %#10f, %#15.9g\n", nodo.x, nodo.y, getCelula().potenciales[iNodo]);
	}

	fclose(ftension);

	sprintf(buffer, "%s/campo/campo.csv.%03d", getCelula().salida.c_str(), nPoisson);
	FILE* fcampo = fopen(buffer, "w");
	
	assert(fcampo > 0);
	fprintf(ftension, "         X,          Y,          campoX,          campoY,           campo,       corriente\n");

	for (int k = 0; k < getCelula().nElems; k++) {
		double corr = sqrt(pow(getCelula().corrElem[k].x, 2) + pow(getCelula().corrElem[k].y, 2));
		double camp = sqrt(pow(getCelula().gradElem[k].x, 2) + pow(getCelula().gradElem[k].y, 2));
		double xMed = 0, yMed = 0;
	
		for (int j = 0; j < getCelula().nodpel; j++) {
			int jNodo = getCelula().elementos[k][j];
			xMed += getCelula().nodos[jNodo].x;
			yMed += getCelula().nodos[jNodo].y;
		}
	
		xMed /= nodpel;
		yMed /= nodpel;
		
		fprintf(fcampo, "%#10f, %#10f, %#15.9g, %#15.9g, %#15.9g, %#15.9g\n", 
			xMed, yMed, getCelula().gradElem[k].x, getCelula().gradElem[k].y, camp, corr);
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
	
	int nNodes = getCelula().nNodes;
	char buffer[512];
	
	sprintf(buffer, "%s/transporte/concentracion.csv.%03d", getCelula().salida.c_str(), nTransporte);
	FILE* fConc = fopen(buffer, "w");
	
	sprintf(buffer, "%s/transporte/molar.csv.%03d", getCelula().salida.c_str(), nTransporte);
	FILE* fMolar = fopen(buffer, "w");

	sprintf(buffer, "%s/transporte/pH.csv.%03d", getCelula().salida.c_str(), nTransporte);
	FILE* fpH = fopen(buffer, "w");

	assert(fConc > 0 && fMolar > 0 && fpH > 0);

	fprintf(fConc,  "         X,          Y,              H+,             OH-,             Na+,             Cl-\n");
	fprintf(fMolar, "         X,          Y,              H+,             OH-,             Na+,             Cl-\n");
	fprintf(fpH,    "         X,          Y,           pH_H+,          pH_OH-\n");

	for (int jNodo = 0; jNodo < nNodes; jNodo++) {
		Nodo nodo = getCelula().nodos[jNodo];
		double conc[NESPS];
		double molar[NESPS];
		double ph[NESPS];

		for (int esp = 0; esp < NESPS; esp++) {
			conc[esp] = getCelula().concs[esp][jNodo];
			molar[esp] = (conc[esp] + 1e-18) / CONCENT;
			ph[esp] = -log10(molar[esp]);
		}

		fprintf(fConc,  "%#10f, %#10f, %#15.9g, %#15.9g, %#15.9g, %#15.9g\n", nodo.x, nodo.y, conc[H_], conc[OH], conc[NA], conc[CL]);
		fprintf(fMolar, "%#10f, %#10f, %#15.9g, %#15.9g, %#15.9g, %#15.9g\n", nodo.x, nodo.y, molar[H_], molar[OH], molar[NA], molar[CL]);
		fprintf(fpH, "%#10f, %#10f, %#15.9g, %#15.9g\n", nodo.x, nodo.y, ph[H_], ph[OH]);
	}

	fclose(fConc);
	fclose(fMolar);
	fclose(fpH);
	nTransporte++;
	printEnd(1, verbose);
}


void EntradaSalida::grabarPoros(Poros& radios, bool verbose) {
	printStart("Grabando poros...", verbose);
	double time = getCelula().time;
	char buffer[512];

	sprintf(buffer, "%s/poros/poros.csv.%03d", getCelula().salida.c_str(), nPoros);
	FILE* fPoros = fopen(buffer, "w");

	assert(fPoros > 0);
	fprintf(fPoros, "      tita,           radio\n");

	for (auto& info : radios.getValores()) {
		for (auto& radio : info.porosGrandes) {
			fprintf(fPoros, "%#10f, %#15.9g\n", info.tita, radio.first);
		}

		for (int i = 0; i < info.porosChicos; i++) {
			fprintf(fPoros, "%#10f, %#15.9g\n", info.tita, info.radioChico);
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

void EntradaSalida::grabarDebug(vector<Triplet<double>>& triplets, VectorXd& rhs, VectorXd& c_ant) {
	char buffer[512];

	sprintf(buffer, "%s/triplets.csv", getCelula().salida.c_str());
	FILE* fTrip = fopen(buffer, "w");
	for (auto& trip : triplets) {
		fprintf(fTrip, "%d, %d, %g\n", trip.col(), trip.row(), trip.value());
	}
	fclose(fTrip);

	sprintf(buffer, "%s/rhs.csv", getCelula().salida.c_str());
	FILE* fRhs = fopen(buffer, "w");
	for (int i = 0; i < rhs.size(); i++) {
		fprintf(fRhs, "%g\n", rhs[i]);
	}
	fclose(fRhs);

	sprintf(buffer, "%s/c_ant.csv", getCelula().salida.c_str());
	FILE* fCAnt = fopen(buffer, "w");
	for (int i = 0; i < c_ant.size(); i++) {
		fprintf(fRhs, "%g\n", c_ant[i]);
	}
	fclose(fCAnt);
}
