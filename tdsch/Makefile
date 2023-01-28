
OBJS1 = numtype.o tdsch.o  # object files

PROG1 = tds   # code name

F90 = gfortran -fopenmp

F90FLAGS = -O3 -funroll-loops # -fexternal-blas # optimization 

LIBS = -llapack -lblas

LDFLAGS = $(LIBS)

all: $(PROG1) 

$(PROG1): $(OBJS1)
	$(F90) $(LDFLAGS) -o $@ $(OBJS1) 

clean:
	rm -f $(PROG1) *.{o,mod} fort.*

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<



