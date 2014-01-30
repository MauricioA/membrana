#include "Problema.h"


#include <vector>
#include <iostream>

using namespace std;

void dameLinea(istringstream& iss) {
	string s = "hola1 hola2 hola3";
	iss.str(s);
}

int main() {

//	string s = "chau1 chau2 chau3";
//
//	istringstream iss(s);
//
//	string h;
//	iss >> h;
//	cout << h << endl;
//
//	dameLinea(iss);
//
//	iss >> h;
//	cout << h << endl;
//
//	return 0;

	int i = 5;

	Problema problema;

	problema.poisson();

	return EXIT_SUCCESS;
}
