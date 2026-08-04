all:
	make -C random
	make -C parameters
	make -C su3
	gfortran -I ../build/mod -c main.f90 -o ../build/obj/main.o
	gfortran ../build/obj/*.o -o ../build/2FLQCD
