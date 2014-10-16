FC = gfortran
FLAGS = -O3 -ffree-line-length-none
#FLAGS = -ffree-line-length-none -g

celula: objs/celula.o objs/controla.o objs/ensamblo.o objs/inicio.o objs/lectura.o objs/poisson.o objs/campo.o objs/finalice.o objs/lbcg.o objs/leemalla.o objs/salidas.o objs/def_constantes.o objs/def_solver.o objs/def_variables.o objs/def_constantes.mod objs/def_solver.mod objs/def_variables.mod objs/armomalla.o objs/calculo_carga.o objs/concentra_time.o objs/def_transpor.o objs/grad_concentra.o objs/transporte.o objs/densidad.o
	$(FC) $(FLAGS) -o celula objs/celula.o objs/controla.o objs/ensamblo.o objs/inicio.o objs/lectura.o objs/poisson.o objs/campo.o objs/finalice.o objs/lbcg.o objs/leemalla.o objs/salidas.o objs/def_constantes.o objs/def_solver.o objs/def_variables.o objs/armomalla.o objs/calculo_carga.o objs/concentra_time.o objs/def_transpor.o objs/grad_concentra.o objs/transporte.o objs/densidad.o

objs/celula.o: fuentes/celula.f90 objs/def_variables.o objs/def_transpor.o
	$(FC) $(FLAGS) -c fuentes/celula.f90 -o objs/celula.o -J objs

objs/controla.o: fuentes/controla.f90 objs/def_solver.o objs/def_constantes.o
	$(FC) $(FLAGS) -c fuentes/controla.f90 -o objs/controla.o -J objs

objs/ensamblo.o: fuentes/ensamblo.f90
	$(FC) $(FLAGS) -c fuentes/ensamblo.f90 -o objs/ensamblo.o -J objs

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

objs/armomalla.o: fuentes/armomalla.f90
	$(FC) $(FLAGS) -c fuentes/armomalla.f90 -o objs/armomalla.o -J objs 

objs/calculo_carga.o: fuentes/calculo_carga.f90
	$(FC) $(FLAGS) -c fuentes/calculo_carga.f90 -o objs/calculo_carga.o -J objs

objs/concentra_time.o: fuentes/concentra_time.f90
	$(FC) $(FLAGS) -c fuentes/concentra_time.f90 -o objs/concentra_time.o -J objs

#objs/def_deforma.o: fuentes/def_deforma.f90 
#	$(FC) $(FLAGS) -c fuentes/def_deforma.f90 -o objs/def_deforma.o -J objs

objs/def_transpor.o: fuentes/def_transpor.f90
	$(FC) $(FLAGS) -c fuentes/def_transpor.f90 -o objs/def_transpor.o -J objs
	
#objs/deforma.o: fuentes/deforma.f90 objs/def_variables.o
#	$(FC) $(FLAGS) -c fuentes/deforma.f90 -o objs/deforma.o -J objs

objs/densidad.o: fuentes/densidad.f90
	$(FC) $(FLAGS) -c fuentes/densidad.f90 -o objs/densidad.o -J objs

objs/grad_concentra.o: fuentes/grad_concentra.f90
	$(FC) $(FLAGS) -c fuentes/grad_concentra.f90 -o objs/grad_concentra.o -J objs

objs/transporte.o: fuentes/transporte.f90
	$(FC) $(FLAGS) -c fuentes/transporte.f90 -o objs/transporte.o -J objs
	
clean: 
	rm -f objs/*.o objs/*.mod celula objs/celula

run: celula
	./celula

new: clean celula
