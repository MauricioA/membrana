#include "Poisson.h"

void Poisson::Poisson(Celula &celula) {

	rhs.resize(nNodes);
	rhs.fill(0.0);
	solucion.resize(nNodes);

	double error = 1.0;

	for (int contador = 0; error > EPSILON_POISSON && contador < N_COTA; contador++) {
		cout << "Armando matriz... " << flush;
		clock_t start = clock();
		vector< Triplet<double> > triplets;

		for (uint elemIdx = 0; elemIdx < elementos.size(); elemIdx++) {

			Elemento elemento = elementos[elemIdx];
			double sigma = sigmas[elemento.material];
			double ef[nodpel];

			Double2D pos[MAXNPEL];
			double esm[MAXNPEL][MAXNPEL];

			for (int i = 0; i < nodpel; i++) {
				int j = elemento[i];
				pos[i].x = nodos[j].x;
				pos[i].y = nodos[j].y;
				ef[i] = 0.0;
			}

			Armado::armadoPoisson(pos, sigma, nodpel, esm);

			/* Condiciones de contorno */
			for (int i = 0; i < nodpel; i++) {
				Nodo nodo = nodos[elemento[i]];
				double adiag = esm[i][i];

				if (nodo.esTierra) {
					for (int j = 0; j < nodpel; j++) {
						esm[i][j] = 0.0;
						ef[j] -= esm[j][i] * TIERRA;
						esm[j][i] = 0.0;
					}
					esm[i][i] = adiag;
					ef[i] = adiag * TIERRA;
				}

				if (nodo.esPotencia) {
					for (int j = 0; j < nodpel; j++) {
						esm[i][j] = 0.0;
						ef[j] -= esm[j][i] * potencial;
						esm[j][i] = 0.0;
					}
					esm[i][i] = adiag;
					ef[i] = adiag * potencial;
				}
			}

			/* Ensamblado */
			for (int i = 0; i < nodpel; i++) {
				rhs[elemento[i]] += ef[i];

				for (int j = 0; j < nodpel; j++) {
					triplets.push_back(Triplet<double>(elemento[i], elemento[j], esm[i][j]));
				}
			}
		}

		matriz.resize(nNodes, nNodes);
		matriz.setFromTriplets(triplets.begin(), triplets.end());
		int time = (clock() - start) / (CLOCKS_PER_SEC / 1000);
		cout << "OK\t\t" << time << "ms\n";

		/* Resolución */
		cout << "Resolviendo... " << flush;
		start = clock();

		SimplicialLDLT< SparseMatrix<double> > cholesky(matriz);
		solucion = cholesky.solve(rhs);

		time = (clock() - start) / (CLOCKS_PER_SEC / 1000);
		cout << "OK\t\t" << time << "ms\n";

		//TODO siempre hace una sola iteración
		error = EPSILON_POISSON * .5;
	}

	cout << "Corriente y campo... " << flush;
	clock_t start = clock();

	campo();
	corriente();

	int time = (clock() - start) / (CLOCKS_PER_SEC / 1000);
	cout << "OK\t\t" << time << "ms\n";
}
