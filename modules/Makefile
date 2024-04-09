include ../make.sys

objects = constants.o linalgwrap.o symmetry.o para.o wanndata.o simp.o lattice.o intRPA.o

mod.a: $(objects)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean:
	rm -rf *.o *.mod *.a *.i90
