include ../make.sys

objects = input.o flib.o  green.o output.o

LIBS = ../modules/mod.a

$(LIBS):
	cd ../modules
	make mod.a
	cd ..

wannband: wannband.f90 $(objects) $(LIBS)
	$(F90) $(F90FLAGS) -fpp -o wannband.x wannband.f90 $(objects) $(LIBS) -I$(MKLINCLUDE) -I../modules -I. $(LAPACK95LIBS) $(MKLLIBS)

wannchi: wannchi.f90 $(objects) $(LIBS)
	$(F90) $(F90FLAGS) -fpp -o wannchi.x wannchi.f90 $(objects) $(LIBS) -I$(MKLINCLUDE) -I../modules -I. $(LAPACK95LIBS) $(MKLLIBS)

clean:
	rm -rf *.o *.mod *.x
