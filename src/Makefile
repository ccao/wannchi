include ../make.sys

objects = input.o green.o output.o compute_chi.o

LIBS = ../modules/mod.a

all: wannband  wannchi

wannband: wannband.f90 $(LIBS) $(objects) 
	$(F90) $(F90FLAGS) -fpp -o wannband.x wannband.f90 $(objects) $(LIBS) -I$(MKLINCLUDE) -I../modules -I. $(LAPACK95LIBS) $(MKLLIBS)

wannchi: wannchi.f90 $(LIBS) $(objects) 
	$(F90) $(F90FLAGS) -fpp -o wannchi.x wannchi.f90 $(objects) $(LIBS) -I$(MKLINCLUDE) -I../modules -I. $(LAPACK95LIBS) $(MKLLIBS)

analyzechi: analyze_chi.f90 $(LIBS) $(objects) 
	$(F90) $(F90FLAGS) -fpp -o analyzechi.x analyze_chi.f90 $(objects) $(LIBS) -I$(MKLINCLUDE) -I../modules -I. $(LAPACK95LIBS) $(MKLLIBS)

$(LIBS):
	cd ../modules
	make mod.a
	cd ..

clean:
	rm -rf *.o *.mod *.x *.i90
