LDFLAGS = -lpgplot -lX11
FC = gfortran
#LIB_TAG = d
LIB_TAG =
EXE = oscill adv1d lorenz

all: $(EXE)

adv1d: adv1d.o
	$(FC) -o $@ $^ `pkg-config --libs plplot$(LIB_TAG)-f95`

adv1d.o: adv1d.f90
	$(FC) -c $(FCFLAGS) `pkg-config --cflags plplot$(LIB_TAG)-f95` $<

oscill: oscill.o
	$(FC) -o $@ $^ `pkg-config --libs plplot$(LIB_TAG)-f95`

oscill.o: oscill.f90
	$(FC) -c $(FCFLAGS) `pkg-config --cflags plplot$(LIB_TAG)-f95` $<

lorenz: lorenz.o
	$(FC) -o $@ $^ `pkg-config --libs plplot$(LIB_TAG)-f95`

lorenz.o: lorenz.f90
	$(FC) -c $(FCFLAGS) `pkg-config --cflags plplot$(LIB_TAG)-f95` $<

clean:
	rm -f *.o *.mod

veryclean: clean
	rm -f $(EXE) *.pdf *~

%.o: %.f90
	$(FC) -c $(FCFLAGS) $<
