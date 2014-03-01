#ifndef POROS_H_
#define POROS_H_

#include "Celula.h"
#include "declaraciones.h"

using namespace std;
using namespace declaraciones;

class Poros {
public:
	Poros(Celula& celula);

	void iteracion();

	void loop();

	struct InfoAngulo {
		int nodosInternos[2];
		int nodosExternos[2];
		double densidad;
		double tita;
		double constIntegral;	// 2 * pi * radio**2 * (cos(tita1) - cos(tita2))
		vector<double> poros;
	};

private:
	Celula* _celula;

	inline Celula& getCelula();

	vector<InfoAngulo> valores;

	bool esNodoExterno(Nodo nodo);

	bool esNodoInterno(Nodo nodo);

	int getPorosEnTita(InfoAngulo& info);

	double getTita(Nodo nodo);

	double getDensidadPromedio();

	double getRadioMaximo();

	double getUnRadio();

	int getNPoros();

};

#endif /* POROS_H_ */
