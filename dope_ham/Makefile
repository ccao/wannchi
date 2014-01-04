#F90=ifort
F90=ifort
#F90FLAGS=-D__MPI

MKLINCLUDE=$(MKLROOT)/include
MKLLIBS=-L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
LAPACK95LIB=$(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a


objects = constants.o wanndata.o orbitals.o dopedata.o translate_crystal.o read_ham.o write_ham.o mix_ham.o

$(objects): %.o:%.f90
	$(F90) $(F90FLAGS) -fpp -c -I$(MKLINCLUDE) $< -o $@

pandora: pandora.f90 $(objects)
	$(F90) $(F90FLAGS) -o pandora.x pandora.f90 $(objects) -I$(MKLINCLUDE) $(LAPACK95LIB) $(MKLLIBS)

hamdiff: hamdiff.f90 $(objects)
	$(F90) $(F90FLAGS) -o hamdiff.x hamdiff.f90 $(objects) -I$(MKLINCLUDE) $(LAPACK95LIB) $(MKLLIBS)

test: test.f90 $(objects)
	$(F90) $(F90FLAGS) -o test.x test.f90 $(objects) -I$(MKLINCLUDE) $(LAPACK95LIB) $(MKLLIBS)

clean:
	rm -rf *.o *.mod *.x
