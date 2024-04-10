include ../make.sys

objects = constants.o linalgwrap.o symmetry.o para.o wanndata.o lattice.o intRPA.o pade.o

mod.a: $(objects)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean:
	rm -rf *.o *.mod *.a *.i90
