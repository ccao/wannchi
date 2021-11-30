include ../make.sys

objects = constants.o linalgwrap.o symmetry.o banddata.o chidata.o para.o wanndata.o impurity.o

mod.a: $(objects)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean:
	rm -rf *.o *.mod *.a *.i90
