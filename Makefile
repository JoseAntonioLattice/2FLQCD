FC := gfortran
OBJDIR := build/obj
MODDIR := build/mod

FLAGS += -I $(MODDIR) -J $(MODDIR)

SRCS := src/pbc/pbc.f90 src/parameters/parameters.f90 src/su3/su3facts.f90 src/arrays/arrays.f90 src/random/random.f90 \
src/starts/starts.f90 src/gauge_links/gauge.f90 src/dirac/dirac.f90 /src/HMC/hybridMC.f90 src/lua/lua.f90 src/dynamics/dynamics.f90

OBJS := $(patsubst %.f90, $(OBJDIR)/%.o, $(notdir $(SRCS)))

EXE := build/2FLQCD
EXETEST := build/test

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

test: $(EXETEST)


$(EXETEST): $(OBJS) build/obj/test.o
	$(FC) $^ -o $@

build/obj/test.o: tests/test.f90
	$(FC) $(FLAGS) -c $< -o $@


test_run:
	build/test
