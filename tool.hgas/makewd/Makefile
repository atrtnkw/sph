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

ifeq ($(use_xc), yes)
HOME_DIR = /home/tanikwat
FDPS_DIR = ${HOME_DIR}/fdps/src
SPH_DIR = ${HOME_DIR}/git-sph/std
FLASH_DIR = ${HOME_DIR}/work/git-sph/link/flash
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
CFLAGS += -I${SPH_DIR} -I${SPH_DIR}/heos -I${SPH_DIR}/aprox13
CLIBS  += -lstdc++ -lm -L${FLASH_DIR} -lflash

OBJ  = main.o ${SPH_DIR}/heos/eos.o #${SPH_DIR}/aprox13/aprox13_main.o \
	${SPH_DIR}/aprox13/aprox13_input.o ${SPH_DIR}/aprox13/aprox13.o \
	${SPH_DIR}/nse/nse.o ${SPH_DIR}/nse/nse_main.o
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
