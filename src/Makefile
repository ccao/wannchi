F90=mpif90

MKLINCLUDE=$(MKLROOT)/include
MKLLIBS=-L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
LAPACK95LIB=$(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a


objects = constants.o wanndata.o read_ham.o estate_mesh.o

$(objects): %.o:%.f90
	$(F90) -c $< -o $@

wannchi: wannchi.f90 $(objects)
	$(F90) -o wannchi.x wannchi.f90 $(objects) -I$(MKLINCLUDE) $(LAPACK95LIB) $(MKLLIBS)
