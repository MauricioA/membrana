#include <stdlib.h>

//TODO borrar los que sobren
#include <iostream>
#include <vector>
#include <fstream>
#include "Problema.h"

using namespace std;

int main() {
	Problema problema;

	cout << problema.getElement(120)[0] << endl;
	cout << problema.getElement(120)[1] << endl;
	cout << problema.getElement(120)[2] << endl;

	/*problema.poisson();*/

	return EXIT_SUCCESS;
}



