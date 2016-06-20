###### AT_PLATFORM ########
use_xc  = yes
###### AT_EOS #############
use_iso7  = no
use_ideal = no
###### AT_MPI_OPENMP ######
use_mpi = yes
use_omp = yes
###########################
###### AT_FLAGS ##############################
#CFLAGS += -DFOR_TUBE_TEST
CFLAGS += -DCOULOMB_CORRECTION
#CFLAGS += -DNBODYLIKE
##############################################
##############################################


ifeq ($(use_xc), yes)
FDPS_DIR = /home/tanikwat/fdps/src
SPH_DIR = /home/tanikwat/git-sph/std
CCC = CC
FCC = ftn
else
FDPS_DIR = ../../../../fdps-jenever0/src
SPH_DIR = ../../std
CCC = g++
FCC = gfortran
endif

PS_PATH  = -I${FDPS_DIR}

CFLAGS += -O3 -march=core-avx2 -ffast-math -funroll-loops #-Wall
CFLAGS += -I${SPH_DIR} -I${SPH_DIR}/heos
CLIBS  += -lstdc++ -lm
CLIBS  += -L./flash -lflash

OBJ  = main.o ${SPH_DIR}/heos/eos.o
ETC  = Makefile

run: $(OBJ) $(ETC)
	$(FCC) $(PS_PATH) $(CFLAGS) -o run $(OBJ) $(CLIBS)

main.o: main.cpp $(ETC)
	$(CCC) $(PS_PATH) $(CFLAGS) -c main.cpp $(CLIBS)

clean:
	rm -f *~
	rm -f *.o[0-9]*
	rm -f core
	rm -f *.s
distclean:
	make clean
	rm -f *.o
allclean:
	make distclean
	rm -f run
