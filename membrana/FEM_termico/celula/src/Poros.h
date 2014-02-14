#ifndef POROS_H_
#define POROS_H_

#include <map>
#include "declaraciones.h"

using namespace std;
using namespace declaraciones;

class Poros {
public:
	Poros();

	void iteracion();

	void loop();

private:
	map<int, ElementoMembrana> mapaMembrana;

};

#endif /* POROS_H_ */
