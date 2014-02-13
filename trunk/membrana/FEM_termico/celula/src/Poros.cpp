#include "Poros.h"

Poros& Poros::instance() {
	static Poros _instance;
	return _instance;
}

Poros::Poros() {
	s = "hola hola 1";
}
