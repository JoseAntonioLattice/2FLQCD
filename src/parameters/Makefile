OBJDIR=../../build/obj
MODDIR=../../build/mod
FFLAGS = -I$(OBJDIR) -J$(MODDIR)

SRCS = parameters.f90
OBJS = $(patsubst %.f90,$(OBJDIR)/%.o,$(SRCS))

all: $(OBJS)

$(OBJDIR)/%.o: %.f90
	gfortran $(FFLAGS) -c $< -o $@
