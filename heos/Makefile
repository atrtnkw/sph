#CCC = mpicxx
CCC = g++
CFLAGS += -std=c++0x -Wall
CFLAGS += -I. -I./eos
CFLAGS += -O3

FCC = gfortran
FFLAGS += -O3
FFLAGS += -Ieos

LIBS += -lstdc++ -lm

OBJ = main.o eos.o
HDR = 
ETC = Makefile

run: $(OBJ) $(HDR) $(ETC)
	$(FCC) $(CFLAGS) -o run $(OBJ) $(LIBS)

main.o: main.cpp $(HDR) $(ETC)
	$(CCC) $(CFLAGS) -c main.cpp -o main.o

eos.o: eos.f
	$(FCC) $(FFLAGS) -c eos.f -o eos.o

clean:
	rm -f *~
	rm -f test.o*
	rm -f core
distclean:
	make clean
	rm -f *.o
allclean:
	make distclean
	rm -f run EOS_3999.dat
