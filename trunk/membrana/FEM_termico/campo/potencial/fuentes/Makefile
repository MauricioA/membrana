FC = gfortran
FLAGS = -O3 -s -ffree-line-length-none

celula: celula.o CONTROLA.o ENSAMBLO.o inicio.o lectura.o poisson.o campo.o finalice.o lbcg.o leemalla.o salidas.o def_constantes.o def_solver.o def_variables.o def_constantes.mod def_solver.mod def_variables.mod
	$(FC) $(FLAGS) -o celula celula.o CONTROLA.o ENSAMBLO.o inicio.o lectura.o poisson.o campo.o finalice.o lbcg.o leemalla.o salidas.o def_constantes.o def_solver.o def_variables.o
	cp celula ../celula

celula.o: celula.f90 def_variables.o
	$(FC) $(FLAGS) -c celula.f90

CONTROLA.o: CONTROLA.F90 def_solver.o def_constantes.o
	$(FC) $(FLAGS) -c CONTROLA.F90
	
ENSAMBLO.o: ENSAMBLO.F90
	$(FC) $(FLAGS) -c ENSAMBLO.F90
	
inicio.o: inicio.f90
	$(FC) $(FLAGS) -c inicio.f90
	
lectura.o: lectura.f90
	$(FC) $(FLAGS) -c lectura.f90
	
poisson.o: poisson.f90
	$(FC) $(FLAGS) -c poisson.f90
	
campo.o: campo.f90
	$(FC) $(FLAGS) -c campo.f90
	
finalice.o: finalice.f90
	$(FC) $(FLAGS) -c finalice.f90
	
lbcg.o: lbcg.f90
	$(FC) $(FLAGS) -c lbcg.f90
	
leemalla.o: leemalla.f90
	$(FC) $(FLAGS) -c leemalla.f90
	
salidas.o: salidas.f90
	$(FC) $(FLAGS) -c salidas.f90
	
def_constantes.o: def_constantes.f90
	$(FC) $(FLAGS) -c def_constantes.f90
	
def_solver.o: def_solver.f90
	$(FC) $(FLAGS) -c def_solver.f90
	
def_variables.o: def_variables.f90
	$(FC) $(FLAGS) -c def_variables.f90

cleanobjs:
	rm -f *.o *.mod
	
clean: cleanobjs
	rm -f celula ../celula
