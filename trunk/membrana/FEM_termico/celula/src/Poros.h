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

	struct ElementoMembrana {
		int NodosInternos[2];	//No necesariamente los nodos del elem!!
		int NodosExternos[2];
		double densidad;
		double radios;
	};

private:
	Celula* _celula;

	map<int, ElementoMembrana> mapaMembrana;

};

#endif /* POROS_H_ */
