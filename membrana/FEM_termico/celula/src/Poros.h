#ifndef POROS_H_
#define POROS_H_

#include <map>
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
		double radios;
		double tita;
		double deltaTita;
		double sin;
	};

private:
	Celula* _celula;

	inline Celula& getCelula();

	vector<InfoAngulo> valores;

	bool esNodoExterno(Nodo nodo);

	bool esNodoInterno(Nodo nodo);

	double getTita(Nodo nodo);

	double getDensidadPromedio();

	double getNPoros();

};

#endif /* POROS_H_ */
