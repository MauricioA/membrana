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

//		const bool operator< (const infoAngulo other) {
//			return tita < other.tita;
//		}

//		inline bool operator< (const InfoAngulo& lhs, const InfoAngulo& rhs){ return lhs.tita < rhs.tita };
	};

private:
	Celula* _celula;

	inline Celula& getCelula();

//	map<int, infoAngulo> mapaMembrana;

	double densidadPromedio;

	vector<InfoAngulo> valores;

	bool esNodoExterno(Nodo nodo);

	bool esNodoInterno(Nodo nodo);

	double getTita(Nodo nodo);

	double getTita(Elemento elemento);

};

#endif /* POROS_H_ */
