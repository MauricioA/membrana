#ifndef POROS_H_
#define POROS_H_

#include <unordered_map>
#include <utility>
#include "Celula.h"
#include "declaraciones.h"

using namespace std;
using namespace declaraciones;

class Poros {
public:
	Poros(Celula& celula);

	void iteracion(double deltaT);

	struct InfoAngulo {
		uint16_t nodosInternos[2];
		uint16_t nodosExternos[2];
		double densidad;
		double tita;
		double area;
		vector<pair<double, double>> porosGrandes;	// <radio, tiempo de creación>
		uint32_t porosChicos;
		double radioChico;
	};

	double getProporsionArea(int iElem);

	double getRadioMaximo();

	int getNPoros();

	int getNPorosChicos();

	int getNPorosGrandes();

	/* Máxima permeabilización */
	double getMaxPAD();

	vector<pair<double, double>> getITVs(double tiempo);

	/* Llamar al iniciar un pulso */
	void nuevoPulso();

	inline vector<InfoAngulo>& getValores() {
		return valores;
	}

private:
	Celula* _celula;

	/* <Elementos de la membrana, info*> */
	unordered_map<int, InfoAngulo*> mapa;

	inline Celula& getCelula();

	vector<InfoAngulo> valores;

	bool esNodoExterno(Nodo nodo);

	bool esNodoInterno(Nodo nodo);

	int getPorosEnTita(InfoAngulo& info);

	double getTita(Elemento nodo);

	double getTita(Nodo nodo);

	double getDensidadPromedio();

	double getUnRadio();

	double getITV(InfoAngulo& info);
	
	double inline actualizarRadio(double radio, double deltaT, double tensionEfectiva, double itv);

	double tau;
	double factorPulso;
	double comienzoPulso = 0;

	void actualizarSigmas();

};

#endif /* POROS_H_ */
