#FCC = gfortran
FFLAGS += -O3

nse_main.o: nse_main.f
	$(FCC) $(FFLAGS) -c nse_main.f -o nse_main.o

nse.o: nse.f
	$(FCC) $(FFLAGS) -c nse.f -o nse.o

clean:
	rm -f *~

distclean:
	make -f Makefile.nse clean
	rm -f *.o

allclean:
	make -f Makefile.nse distclean
