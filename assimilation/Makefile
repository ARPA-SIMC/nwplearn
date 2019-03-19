#LIB_TAG = d
LIB_TAG =
FC = gfortran
FCFLAGS = -g -fcheck=bounds -Wall `pkg-config --cflags plplot$(LIB_TAG)-f95`
LDFLAGS = -llapack `pkg-config --libs plplot$(LIB_TAG)-f95`

EXE = opt_int opt_int_multi
all: $(EXE)

missing_valuse.o: kinds.o
func_mo.o: kinds.o missing_values.o
coord_mo.o: kinds.o
obs_mo.o: coord_mo.o func_mo.o kinds.o
grid_mo.o: coord_mo.o kinds.o
an_opt_int_mo.o: grid_mo.o obs_mo.o coord_mo.o kinds.o
opt_int.o: func_mo.o coord_mo.o obs_mo.o grid_mo.o an_opt_int_mo.o
opt_int: func_mo.o coord_mo.o obs_mo.o grid_mo.o func_mo.o an_opt_int_mo.o kinds.o
opt_int1d.o: func_mo.o coord_mo.o obs_mo.o grid_mo.o func_mo.o an_opt_int_mo.o
opt_int1d: func_mo.o coord_mo.o obs_mo.o grid_mo.o func_mo.o an_opt_int_mo.o

opt_int_multi.o: func_mo.o coord_mo.o obs_mo.o grid_mo.o an_opt_int_mo.o
opt_int_multi: func_mo.o coord_mo.o obs_mo.o grid_mo.o func_mo.o an_opt_int_mo.o kinds.o

%.o: %.f90
	$(FC) -c $(FCFLAGS) $<

%.o: %.F90
	$(FC) -c $(FCFLAGS) $<

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

.SUFFIXES: .F90 .f90

clean:
	rm -f *.o *.mod

veryclean: clean
	rm -f *~ $(EXE)
