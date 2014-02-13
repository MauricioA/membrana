#ifndef POROS_H_
#define POROS_H_

#include <string>

using namespace std;

class Poros {
public:
	static Poros& instance();

	string s;

private:
	Poros();
	Poros(Poros const&);
	void operator=(Poros const&);

};

#endif /* POROS_H_ */
