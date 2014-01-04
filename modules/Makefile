F90=mpif90
F90FLAGS=-D__MPI

MKLINCLUDE=$(MKLROOT)/include
MKLLIBS=-L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
LAPACK95LIBS=$(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a


objects = banddata.f90 chidata.f90 constants.f90 para.f90 wanndata.f90

$(objects): $.o:$.f90
	$(F90) $(F90FLAGS) -fpp -c -I$(MKLINCULDE) $< -o $@

clean:
	rm -rf *.o *.mod
