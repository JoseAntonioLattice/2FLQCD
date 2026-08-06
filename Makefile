all:
	make -C src

run:
	echo "input/input_parameters.nml" | build/2FLQCD

test:
	make -C src/pbc
	make -C src/random
	make -C src/su3
	make -C src/starts
	gfortran -I build/mod -c tests/test.f90 -o build/obj/test.o
	gfortran build/obj/random.o build/obj/su3facts.o build/obj/test.o build/obj/pbc.o build/obj/starts.o -o build/test
	build/test
