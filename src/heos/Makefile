FFLAGS = -O3

eos.o: eos.f
	$(FCC) $(FFLAGS) -c eos.f -o eos.o

clean:
	rm -f *~
distclean:
	make clean
	rm -f *.o
allclean:
	make distclean
