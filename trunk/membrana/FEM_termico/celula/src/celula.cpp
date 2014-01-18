#include <iostream>
#include <vector>
#include <cmath>
#include <stdlib.h>
//#include <Eigen/Dense>

//using Eigen::MatrixXd;
using namespace std;

const int 	 NODPEL 	= 3;
const double TOLER_AREA = 1e-6;

void poisson();
void armado(double x[], double y[], double ef[], double qe, double esm[3][3], double sigma);

int main() {

	poisson();

	return EXIT_SUCCESS;
}

class Elemento {
	private:
	int nodes[NODPEL];

	public:
	int& operator[](int i) {
		return nodes[i];
	}
};

void poisson() {
	//mockeados
	int nElems = 100;
	vector<int> materiales(nElems);
	for (int i = 0; i < nElems; i++) materiales[i] = 1;
	double sigmas[] = {1e-3, 2e-3, 3e-3};
	vector<Elemento> elementos(nElems);
	vector<double> coordX(nElems);
	vector<double> coordY(nElems);
	vector<double> solucionAnterior(nElems);
	vector<double> gradElem_x(nElems);
	vector<double> gradElem_y(nElems);

	const int nCota = 10;
	const double epsilon = 1e-3;

	double error = 1.;
	double ns[NODPEL];
	double ef[NODPEL];

	for (int contador = 0; error > epsilon && contador < nCota; contador++) {
		int an = 0, ad = 0, rhs = 0;

		for (int elem = 0; elem < nElems; elem++) {
			double sigma = sigmas[materiales[elem]];
			double qe = 0.;
			double x[NODPEL], y[NODPEL];
			double sol[NODPEL];
			double esm[3][3];

			for (int i = 0; i < NODPEL; i++) {
				int j = ns[i] = elementos[elem][i];
				x[i] = coordX[j];
				y[i] = coordY[j];
				sol[i] = solucionAnterior[j];
			}

			gradElem_x[elem] = gradElem_y[elem] = 0.;

			armado(x, y, ef, qe, esm, sigma);



		}
	}
}

void armado(double x[], double y[], double ef[], double qe, double esm[3][3], double sigma) {
	double b[] = {
		y[1] - y[2],
		y[2] - y[0],
		y[0] - y[1],
	};

	double c[] = {
		x[2] - x[1],
		x[0] - x[2],
		x[1] - x[0],
	};

	double determinante =
		+ x[1]*y[2] + x[2]*y[0] + x[0]*y[1]
		- x[1]*y[0] - x[2]*y[1] - x[0]*y[2];

	double rMed = (x[0] + x[1] + x[2]) / 3.;

	if (abs(determinante) < TOLER_AREA) {
		cerr << "Error area es cero\n";
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < 3; i++) {
		ef[i] = 0.;
		for (int j = 0; j < 3; j++) {
			double a = i == j ? 2. : 1.;
			ef[i] += (qe * determinante * M_PI / 12.) * (a * x[j]);
			esm[i][j] = (sigma * (b[i] * b[j] + c[i] * c[j])) * M_PI * rMed * determinante;
		}
	}
}

