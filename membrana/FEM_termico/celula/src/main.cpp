#include "Problema.h"
#include <iostream>





using namespace std;

class Test {
public:
	vector<int> vec;

	Test() {
		vec.push_back(1);
		vec.push_back(2);
		vec.push_back(3);
	}

	vector<int> &getVec() {
		return vec;
	}
};


int main() {
	BREAKPOINT

	Celula celula;

	cout << celula.nodpel << endl;
	celula.nodpel = 99;
	cout << celula.nodpel << endl;




//	Celula problema;
//
//	problema.transporte();
//
//	return EXIT_SUCCESS;
}
