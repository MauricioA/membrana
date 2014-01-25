#TODO mejorar este makefile
#TODO ver el tema de los includes

celula: src/celula.cpp src/declaraciones.cpp src/Problema.h src/Problema.cpp
	g++ -s -O3 -msse3 -c -o celula.o src/celula.cpp -I eigen/
	g++ -s -O3 -msse3 -c -o declaraciones.o src/declaraciones.cpp -I eigen/
	g++ -s -O3 -msse3 -c -o problema.o src/Problema.cpp -I eigen/
	g++ -s -O3 -msse3 -o celula celula.o problema.o declaraciones.o 
