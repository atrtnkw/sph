include ~/setting/setting.mk

ifeq ($(use_i9),yes)
FCC = gfortran
endif

ifeq ($(use_xc),yes)
FCC = ftn
endif

run: public_torch.o
	$(FCC) public_torch.o -o run

public_torch.o: public_torch.f90
	$(FCC) -c public_torch.f90

clean:
	rm -f *~

distclean:
	make clean
	rm -f *.o run
