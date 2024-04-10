include ../make.sys

objects = input.o compute_chi.o green.o

LIBS = ../modules/mod.a

all: wannchi

wannchi: wannchi.f90 $(LIBS) $(objects) 
	$(F90) $(F90FLAGS) -fpp -o wannchi.x wannchi.f90 $(objects) $(LIBS) -I$(MKLINCLUDE) -I../modules -I. $(LAPACK95LIBS) $(MKLLIBS)

$(LIBS):
	cd ../modules
	make mod.a
	cd ..

clean:
	rm -rf *.o *.mod *.x *.i90
