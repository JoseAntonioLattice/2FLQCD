all:
	make -C src

run:
	echo "input/input_parameters.nml" | build/2FLQCD

test:
	make -C src/random
	make -C src/su3
	gfortran -I build/mod -c tests/test.f90 -o build/obj/test.o
	gfortran build/obj/random.o build/obj/su3facts.o build/obj/test.o -o build/test
	build/test
