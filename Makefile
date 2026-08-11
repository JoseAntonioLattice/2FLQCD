FC := gfortran
OBJDIR := build/obj
MODDIR := build/mod

FLAGS += -I $(MODDIR) -J $(MODDIR)

SRCS := src/pbc/pbc.f90 src/parameters/parameters.f90 src/su3/su3facts.f90 src/arrays/arrays.f90 src/random/random.f90 \
src/starts/starts.f90 src/gauge_links/gauge.f90 src/HMC/hybridMC.f90 src/dynamics/dynamics.f90

OBJS := $(patsubst %.f90, $(OBJDIR)/%.o, $(notdir $(SRCS)))

EXE := build/2FLQCD 

SRC_DIRS := $(sort $(dir $(SRCS)))
VPATH := $(SRC_DIRS)


all: $(EXE)

$(EXE): $(OBJS) build/obj/main.o  
	$(FC) $^ -o $@

build/obj/main.o: src/main.f90
	$(FC) $(FLAGS) -c $< -o $@

$(OBJDIR)/%.o: %.f90
	$(FC) $(FLAGS) -c $< -o $@


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
