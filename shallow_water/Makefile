#LIB_TAG = d
LIB_TAG =
FC = gfortran
FCFLAGS = -g -fcheck=bounds -Wall
FCFLAGSPL = $(FCFLAGS) `pkg-config --cflags plplot$(LIB_TAG)-f95`
LDFLAGS = `pkg-config --libs plplot$(LIB_TAG)-f95`

EXE = shallow_water
all: $(EXE)

shallow_water.o: shallow_water_mo.o
shallow_water: shallow_water_mo.o

%.o: %.f90
	$(FC) -c $(FCFLAGSPL) $<

%.o: %.F90
	$(FC) -c $(FCFLAGSPL) $<

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

.SUFFIXES: .F90 .f90

clean:
	rm -f *.o *.mod

veryclean: clean
	rm -f *~ $(EXE)
