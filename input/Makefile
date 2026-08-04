src = parameters.f90
obj = $(patsubst %, ../../build/obj./%.o, $(src))

all : $(obj)

../../build/obj/%.o: %.f90
	gfortran -c -J ../../build/mod $< -o $@
