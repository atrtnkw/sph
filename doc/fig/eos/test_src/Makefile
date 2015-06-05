#CCC = mpicxx
CCC = g++
CFLAGS += -std=c++0x -Wall
CFLAGS += -I. -I./eos
CFLAGS += -O3

FCC = gfortran
FFLAGS += -O3
FFLAGS += -Ieos

LIBS += -lstdc++ -lm

OBJ = main.o eos/eos.o
HDR = 
ETC = Makefile

run: $(OBJ) $(HDR) $(ETC)
	$(FCC) $(CFLAGS) -o run $(OBJ) $(LIBS)

main.o: main.cpp $(HDR) $(ETC)
	$(CCC) $(CFLAGS) -c main.cpp -o main.o

eos/eos.o: eos/eos.f
	make -C eos

clean:
	rm -f *~
	rm -f test.o*
	rm -f core
	make -C eos clean
distclean:
	make clean
	rm -f *.o
	make -C eos distclean
allclean:
	make distclean
	rm -f run EOS_3999.dat
	make -C eos allclean
