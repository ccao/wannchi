include ../make.sys

objects = input.o flib.o  interpolate_bands.o calc_occ.o compute_chi.o output_chi.o

LIBS = ../modules/mod.a

$(LIBS):
	cd ../modules
	make mod.a
	cd ..

wannchi: wannchi.f90 $(objects) $(LIBS)
	$(F90) $(F90FLAGS) -fpp -o wannchi.x wannchi.f90 $(objects) $(LIBS) -I$(MKLINCLUDE) -I../modules -I. $(LAPACK95LIBS) $(MKLLIBS)

clean:
	rm -rf *.o *.mod *.x
