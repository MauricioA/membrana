FC = gfortran
FLAGS = -O3 -s -ffree-line-length-none

celula: objs/celula.o objs/CONTROLA.o objs/ENSAMBLO.o objs/inicio.o objs/lectura.o objs/poisson.o objs/campo.o objs/finalice.o objs/lbcg.o objs/leemalla.o objs/salidas.o objs/def_constantes.o objs/def_solver.o objs/def_variables.o objs/def_constantes.mod objs/def_solver.mod objs/def_variables.mod
	$(FC) $(FLAGS) -o celula objs/celula.o objs/CONTROLA.o objs/ENSAMBLO.o objs/inicio.o objs/lectura.o objs/poisson.o objs/campo.o objs/finalice.o objs/lbcg.o objs/leemalla.o objs/salidas.o objs/def_constantes.o objs/def_solver.o objs/def_variables.o

objs/celula.o: fuentes/celula.f90 objs/def_variables.o
	$(FC) $(FLAGS) -c fuentes/celula.f90 -o objs/celula.o -J objs

objs/CONTROLA.o: fuentes/CONTROLA.F90 objs/def_solver.o objs/def_constantes.o
	$(FC) $(FLAGS) -c fuentes/CONTROLA.F90 -o objs/CONTROLA.o -J objs
	
objs/ENSAMBLO.o: fuentes/ENSAMBLO.F90
	$(FC) $(FLAGS) -c fuentes/ENSAMBLO.F90 -o objs/ENSAMBLO.o -J objs
	
objs/inicio.o: fuentes/inicio.f90
	$(FC) $(FLAGS) -c fuentes/inicio.f90 -o objs/inicio.o -J objs
	
objs/lectura.o: fuentes/lectura.f90
	$(FC) $(FLAGS) -c fuentes/lectura.f90 -o objs/lectura.o -J objs
	
objs/poisson.o: fuentes/poisson.f90
	$(FC) $(FLAGS) -c fuentes/poisson.f90 -o objs/poisson.o -J objs
	
objs/campo.o: fuentes/campo.f90
	$(FC) $(FLAGS) -c fuentes/campo.f90 -o objs/campo.o -J objs
	
objs/finalice.o: fuentes/finalice.f90
	$(FC) $(FLAGS) -c fuentes/finalice.f90 -o objs/finalice.o -J objs
	
objs/lbcg.o: fuentes/lbcg.f90
	$(FC) $(FLAGS) -c fuentes/lbcg.f90 -o objs/lbcg.o -J objs
	
objs/leemalla.o: fuentes/leemalla.f90
	$(FC) $(FLAGS) -c fuentes/leemalla.f90 -o objs/leemalla.o -J objs
	
objs/salidas.o: fuentes/salidas.f90
	$(FC) $(FLAGS) -c fuentes/salidas.f90 -o objs/salidas.o -J objs
	
objs/def_constantes.o: fuentes/def_constantes.f90
	$(FC) $(FLAGS) -c fuentes/def_constantes.f90 -o objs/def_constantes.o -J objs
	
objs/def_solver.o: fuentes/def_solver.f90
	$(FC) $(FLAGS) -c fuentes/def_solver.f90 -o objs/def_solver.o -J objs
	
objs/def_variables.o: fuentes/def_variables.f90
	$(FC) $(FLAGS) -c fuentes/def_variables.f90 -o objs/def_variables.o -J objs
	
clean: 
	rm -f objs/*.o objs/*.mod celula objs/celula
