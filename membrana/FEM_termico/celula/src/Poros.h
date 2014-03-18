#ifndef POROS_H_
#define POROS_H_

#include <unordered_map>
#include "Celula.h"
#include "declaraciones.h"

using namespace std;
using namespace declaraciones;

class Poros {
public:
	Poros(Celula& celula);

	void iteracion(double deltaT, double tiempo);

	void loop();

	struct InfoAngulo {
		int nodosInternos[2];
		int nodosExternos[2];
		double densidad;
		double tita;
		double area;
		vector<double> poros;
	};

	double getProporsionArea(int iElem);

	double getRadioMaximo();

private:
	Celula* _celula;

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
	
	double tau;

	int getNPoros();

};

#endif /* POROS_H_ */
