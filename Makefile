all:
	make -C src

run:
	echo "input/input_parameters.nml" | build/2FLQCD

test:
	make -C src/parameters
	$(MAKE) -C src/pbc
	$(MAKE) -C src/random
	$(MAKE) -C src/su3
	$(MAKE) -C src/starts
	make -C src/gauge_links
	make -C src/HMC
	gfortran -I build/mod -c tests/test.f90 -o build/obj/test.o
	gfortran build/obj/parameters.o build/obj/random.o build/obj/su3facts.o build/obj/test.o build/obj/pbc.o build/obj/starts.o build/obj/gauge.o build/obj/hybridMC.o -o build/test
	build/test
