include ../make.sys

objects = constants.o banddata.o chidata.o para.o wanndata.o impurity.o

mod.a: $(objects)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean:
	rm -rf *.o *.mod *.a
