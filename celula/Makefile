#EN WINDOWS USAR VISUAL STUDIO
#Pueden faltar dependencias

COMP=g++
FLAGS=-O3 -std=c++0x -fopenmp -Wall -Werror -I eigen/
OBJS=bin/Armado.o bin/Celula.o bin/EntradaSalida.o bin/Poisson.o bin/Poros.o bin/Transporte.o

.PHONY: all clean

celula: src/main.cpp $(OBJS)
	$(COMP) $(FLAGS) -o $@ $< $(OBJS)

bin/%.o: src/%.cpp src/%.h
	$(COMP) $(FLAGS) -c -o $@ $<

clean:
	rm -f celula
	rm -f bin/*.o

all: clean celula
