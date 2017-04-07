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
CFLAGS += -DUSE_HELMHOLTZ
#CFLAGS += -DCOULOMB_CORRECTION
##############################################
##############################################

HOME_DIR = /home/tanikwat
FDPS_DIR = ${HOME_DIR}/fdps/src
SPH_DIR = ${HOME_DIR}/git-sph/std
FLASH_DIR = ${HOME_DIR}/work/git-sph/link/flash
CCC = CC
FCC = ftn

PS_PATH  = -I${FDPS_DIR}

CFLAGS += -O3 -march=core-avx2 -ffast-math -funroll-loops #-Wall
CFLAGS += -I${SPH_DIR} -I${SPH_DIR}/heos -I${SPH_DIR}/aprox13
CLIBS  += -lstdc++ -lm -L${FLASH_DIR} -lflash

OBJ  = main.o coldwd.o ${SPH_DIR}/heos/eos.o
ETC  = Makefile

run: $(OBJ) $(ETC)
	$(FCC) $(PS_PATH) $(CFLAGS) -o run $(OBJ) $(CLIBS)

main.o: main.cpp $(ETC)
	$(CCC) $(PS_PATH) $(CFLAGS) -c main.cpp $(CLIBS)

coldwd.o: coldwd.f $(ETC)
	$(FCC) $(PS_PATH) $(CFLAGS) -c coldwd.f

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
